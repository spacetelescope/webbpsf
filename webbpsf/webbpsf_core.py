"""
============
WebbPSF Core
============

An object-oriented modeling system for the JWST instruments.

Classes:
  * SpaceTelescopeInstrument
    * JWInstrument
      * MIRI
      * NIRCam
      * NIRSpec
      * NIRISS
      * FGS

WebbPSF makes use of python's ``logging`` facility for log messages, using
the logger name "webbpsf".

Code by Marshall Perrin <mperrin@stsci.edu>
"""
import os
import glob
from collections import namedtuple, OrderedDict
import numpy as np
import scipy.interpolate, scipy.ndimage

import astropy
import astropy.io.fits as fits
import astropy.io.ascii as ioascii
import astropy.units as units

import poppy

import pysiaf

from . import conf
from . import utils
from . import optics
from . import DATA_VERSION_MIN
from . import distortion
from . import gridded_library
from . import opds
from . import constants
import webbpsf.mast_wss

try:
    from .version import version
except ImportError:
    version = ''

try:
    _HAS_SYNPHOT = poppy.instrument._HAS_SYNPHOT
except AttributeError:
    _HAS_SYNPHOT = False
if _HAS_SYNPHOT:
    import synphot
import logging

_log = logging.getLogger('webbpsf')

Filter = namedtuple('Filter', ['name', 'filename', 'default_nlambda'])


class SpaceTelescopeInstrument(poppy.instrument.Instrument):
    """ A generic Space Telescope Instrument class.

    *Note*: Do not use this class directly - instead use one of the specific instrument subclasses!

    This class provides a simple interface for modeling PSF formation through the instrument,
    with configuration options and software interface loosely resembling the configuration of the instrument
    hardware mechanisms.

    This module currently only provides a modicum of error checking, and relies on the user
    being knowledgable enough to avoid trying to simulate some physically impossible or just plain silly
    configuration (such as trying to use a FQPM with the wrong filter).

    The instrument constructors do not take any arguments. Instead, create an instrument object and then
    configure the `filter` or other attributes as desired. The most commonly accessed parameters are
    available as object attributes: `filter`, `image_mask`, `pupil_mask`, `pupilopd`. More advanced
    configuration can be done by editing the :ref:`SpaceTelescopeInstrument.options` dictionary, either by
    passing options to ``__init__`` or by directly editing the dict afterwards.
    """
    telescope = "Generic Space Telescope"
    options = {}  # options dictionary
    """ A dictionary capable of storing other arbitrary options, for extensibility. The following are all optional, and
    may or may not be meaningful depending on which instrument is selected.

    This is a superset of the options provided in :py:attr:`poppy.Instrument.options`.

    Parameters
    ----------
    source_offset_r : float
        Radial offset of the target from the center, in arcseconds
    source_offset_theta : float
        Position angle for that offset, in degrees CCW.
    pupil_shift_x, pupil_shift_y : float
        Relative shift of the intermediate (coronagraphic) pupil in X and Y
        relative to the telescope entrance pupil, expressed as a decimal between -1.0-1.0
        Note that shifting an array too much will wrap around to the other side unphysically, but
        for reasonable values of shift this is a non-issue.  This option only has an effect for optical models that
        have something at an intermediate pupil plane between the telescope aperture and the detector.
    pupil_rotation : float
        Relative rotation of the intermediate (coronagraphic) pupil relative to
        the telescope entrance pupil, expressed in degrees counterclockwise.
        This option only has an effect for optical models that have something at
        an intermediate pupil plane between the telescope aperture and the detector.
    rebin : bool
        For output files, write an additional FITS extension including a version of the output array
        rebinned down to the actual detector pixel scale?
    jitter : string "gaussian" or None
        Type of jitter model to apply. Currently only convolution with a Gaussian kernel of specified
        width `jitter_sigma` is implemented. (default: None)
    jitter_sigma : float
        Width of the jitter kernel in arcseconds (default: 0.006 arcsec, 1 sigma per axis)
    parity : string "even" or "odd"
        You may wish to ensure that the output PSF grid has either an odd or even number of pixels.
        Setting this option will force that to be the case by increasing npix by one if necessary.
        Note that this applies to the number detector pixels, rather than the subsampled pixels if oversample > 1.
    force_coron : bool
        Set this to force full coronagraphic optical propagation when it might not otherwise take place
        (e.g. calculate the non-coronagraphic images via explicit propagation to all optical surfaces, FFTing
        to intermediate pupil and image planes whether or not they contain any actual optics, rather than
        taking the straight-to-MFT shortcut)
    no_sam : bool
        Set this to prevent the SemiAnalyticMethod coronagraph mode from being
        used when possible, and instead do the brute-force FFT calculations.
        This is usually not what you want to do, but is available for comparison tests.
        The SAM code will in general be much faster than the FFT method,
        particularly for high oversampling.

    """
    _detectors = {}
    """
    Dictionary mapping detector names to detector or wavefront information in some manner.
    The specific meaning of this mapping must be defined by subclasses as part of their
    implementation.

    (Subclasses must populate this at `__init__`.)
    """
    _detector = None
    """
    The name of the currently selected detector. Must be a key in _detectors, as validated by the
    `detector` property setter.

    (Subclasses must populate this at `__init__`.)
    """

    def _get_filters(self):
        filter_table = ioascii.read(os.path.join(self._WebbPSF_basepath, self.name, 'filters.tsv'))
        filter_info = {}
        filter_list = []  # preserve the order from the table

        for filter_row in filter_table:
            filter_filename = os.path.join(
                self._WebbPSF_basepath,
                self.name,
                'filters',
                '{}_throughput.fits'.format(filter_row['filter'])
            )
            filter_info[filter_row['filter']] = Filter(
                name=filter_row['filter'],
                filename=filter_filename,
                default_nlambda=filter_row['nlambda']
            )
            filter_list.append(filter_row['filter'])
        return filter_list, filter_info

    def _get_default_nlambda(self, filtername):
        """ Return the default # of wavelengths to be used for calculation by a given filter """
        return self._filters[filtername].default_nlambda

    def __init__(self, name="", pixelscale=0.064):
        self.name = name

        self._WebbPSF_basepath, self._data_version = utils.get_webbpsf_data_path(
            data_version_min=DATA_VERSION_MIN, return_version=True)

        self._datapath = os.path.join(self._WebbPSF_basepath, self.name)
        self._image_mask = None
        self._pupil_mask = None

        self.pupil = None
        "Filename *or* fits.HDUList for the pupil mask."
        self.pupilopd = None  # This can optionally be set to a tuple indicating (filename, slice in datacube)
        """Filename *or* fits.HDUList for pupil OPD.

        This can be either a full absolute filename, or a relative name in which case it is
        assumed to be within the instrument's `data/OPDs/` directory, or an actual
        fits.HDUList object corresponding to such a file. If the file contains a
        datacube, you may set this to a tuple (filename, slice) to select a
        given slice, or else the first slice will be used."""
        self.pupil_radius = None  # Set when loading FITS file in get_optical_system

        self.options = {}  # dict for storing other arbitrary options.

        # filter_list   available filter names in order by wavelength for public api
        # _filters      a dict of named tuples with name, filename, & default_nlambda
        #               with the filter name as the key
        self.filter_list, self._filters = self._get_filters()

        # choose a default filter, in case the user doesn't specify one
        self.filter = self.filter_list[0]

        self._rotation = None

        self._image_mask = None
        self.image_mask_list = []
        "List of available image_masks"

        self._pupil_mask = None
        self.pupil_mask_list = []
        "List of available pupil_masks"

        self.pixelscale = pixelscale
        "Detector pixel scale, in arcsec/pixel"
        self._spectra_cache = {}  # for caching synphot results.

        # n.b.STInstrument subclasses must set these
        self._detectors = {}
        self._detector = None
        self._detector_npixels = 2048

    @property
    def image_mask(self):
        """Currently selected image plane mask, or None for direct imaging"""
        return self._image_mask

    @image_mask.setter
    def image_mask(self, name):
        if name == "": name = None
        if name is not None:
            if name in self.image_mask_list:
                pass  # there's a perfect match, this is fine.
            else:
                name = name.upper()  # force to uppercase
                if name not in self.image_mask_list:  # if still not found, that's an error.
                    raise ValueError("Instrument %s doesn't have an image mask called '%s'." % (self.name, name))
        self._image_mask = name
        if hasattr(self, '_image_mask_apertures') and name in self._image_mask_apertures:
            self.set_position_from_aperture_name(self._image_mask_apertures[name])

    @property
    def pupil_mask(self):
        """Currently selected Lyot pupil mask, or None for direct imaging"""
        return self._pupil_mask

    @pupil_mask.setter
    def pupil_mask(self, name):
        if name == "":
            name = None
        if name is not None:
            if name in self.pupil_mask_list:
                pass  # there's a perfect match, this is fine.
            else:
                name = name.upper()  # force to uppercase
                if name not in self.pupil_mask_list:
                    raise ValueError("Instrument %s doesn't have a pupil mask called '%s'." % (self.name, name))

        self._pupil_mask = name

    def __str__(self):
        return "<{telescope}: {instrument_name}>".format(telescope=self.telescope, instrument_name=self.name)

    @property
    def detector(self):
        """Detector selected for simulated PSF

        Used in calculation of field-dependent aberrations. Must be
        selected from detectors in the `detector_list` attribute.
        """
        return self._detector

    @detector.setter
    def detector(self, value):
        if value.upper() not in self.detector_list:
            raise ValueError("Invalid detector. Valid detector names are: {}".format(', '.join(self.detector_list)))
        self._detector = value.upper()
        self._update_aperturename()  # automatically set an appropriate aperture name

    @property
    def detector_list(self):
        """Detectors on which the simulated PSF could lie"""
        return sorted(self._detectors.keys())

    @property
    def detector_position(self):
        """The pixel position in (X, Y) on the detector, relative to the currently-selected SIAF aperture subarray.
        By default the SIAF aperture will correspond to the full-frame detector, so (X,Y) will in that case be
        absolute (X,Y) pixels on the detector. But if you select a subarray aperture name from the SIAF, then
        the (X,Y) are interpreted as (X,Y) within that subarray.

        Please note, this is X,Y order - **not** a Pythonic y,x axes ordering.
        """
        return self._detector_position

    @detector_position.setter
    def detector_position(self, position):
        try:
            x, y = map(int, position)
        except ValueError:
            raise ValueError("Detector pixel coordinates must be pairs of nonnegative numbers, not {}".format(position))
        if x < 0 or y < 0:
            raise ValueError("Detector pixel coordinates must be nonnegative integers")
        if x > self._detector_npixels - 1 or y > self._detector_npixels - 1:
            raise ValueError("The maximum allowed detector pixel coordinate value is {}".format(
                self._detector_npixels - 1))

        self._detector_position = (int(position[0]), int(position[1]))

    @property
    def aperturename(self):
        """ SIAF aperture name for detector pixel to sky coords transformations"""
        return self._aperturename

    @aperturename.setter
    def aperturename(self, value):
        # Override in subclass to provide more specific functionality
        self._aperturename = value

    def _update_aperturename(self):
        """ Update SIAF aperture name after change in detector or other relevant properties
        """
        self.aperturename = self._detectors[self._detector]

    def _get_fits_header(self, result, options):
        """ populate FITS Header keywords """
        super(SpaceTelescopeInstrument, self)._get_fits_header(result, options)
        result[0].header['FILTER'] = (self.filter, 'Filter name')
        if self.image_mask is not None:
            result[0].header['CORONMSK'] = (self.image_mask, "Image plane mask")
        if self.pupil_mask is not None:
            result[0].header['PUPIL'] = (self.pupil_mask, "Pupil plane mask")

        result[0].header['VERSION'] = (version, "WebbPSF software version")
        result[0].header['DATAVERS'] = (self._data_version, "WebbPSF reference data files version")

        result[0].header['DET_NAME'] = (self.detector, "Name of detector on this instrument")

        # Correct detector pixel coordinates to allow for even arrays to be centered on half pixel boundary
        dpos = np.asarray(self.detector_position, dtype=float)
        oversamp = result[0].header['OVERSAMP']
        size = result[0].data.shape[0]

        if size / oversamp % 2 == 0: dpos += 0.5  # even arrays must be at a half pixel

        result[0].header['DET_X'] = (dpos[0], "Detector X pixel position of array center")
        result[0].header['DET_Y'] = (dpos[1], "Detector Y pixel position of array center")

        for key in self._extra_keywords:
            result[0].header[key] = self._extra_keywords[key]

    def _calc_psf_format_output(self, result, options):
        """ Apply desired formatting to output file:
                 - rebin to detector pixel scale if desired
                 - set up FITS extensions if desired
                 - output either the oversampled, rebinned, or both

            Modifies the 'result' HDUList object.
        """
        output_mode = options.get('output_mode', conf.default_output_mode)

        if output_mode == 'Mock JWST DMS Output':  # TODO:jlong: move out to JWInstrument
            # first rebin down to detector sampling
            # then call mockdms routines to embed in larger detector etc
            raise NotImplementedError('Not implemented yet')
        else:
            poppy.Instrument._calc_psf_format_output(self, result, options)

    def get_optical_system(self, fft_oversample=2, detector_oversample=None,
                            fov_arcsec=2, fov_pixels=None, options=None):
        """ Return an OpticalSystem instance corresponding to the instrument as currently configured.

        When creating such an OpticalSystem, you must specify the parameters needed to define the
        desired sampling, specifically the oversampling and field of view.


        Parameters
        ----------

        fft_oversample : int
            Oversampling factor for intermediate plane calculations. Default is 2
        detector_oversample: int, optional
            By default the detector oversampling is equal to the intermediate calculation oversampling.
            If you wish to use a different value for the detector, set this parameter.
            Note that if you just want images at detector pixel resolution you will achieve higher fidelity
            by still using some oversampling (i.e. *not* setting `oversample_detector=1`) and instead rebinning
            down the oversampled data.
        fov_pixels : float
            Field of view in pixels. Overrides fov_arcsec if both set.
        fov_arcsec : float
            Field of view, in arcseconds. Default is 2


        Returns
        -------
        osys : poppy.OpticalSystem
            an optical system instance representing the desired configuration.

        """

        _log.info("Creating optical system model:")

        self._extra_keywords = OrderedDict()  # Place to save info we later want to put
                                              # into the FITS header for each PSF.

        if options is None: options = self.options
        if detector_oversample is None: detector_oversample = fft_oversample

        _log.debug("Oversample: %d  %d " % (fft_oversample, detector_oversample))
        optsys = poppy.OpticalSystem(
            name='{telescope}+{instrument}'.format(telescope=self.telescope, instrument=self.name),
            oversample=fft_oversample
        )
        # For convenience offsets can be given in cartesian or radial coords
        if 'source_offset_x' in options or 'source_offset_y' in options:
            offx = options.get('source_offset_x', 0)
            offy = options.get('source_offset_y', 0)
            optsys.source_offset_r = np.sqrt(offx ** 2 + offy ** 2)
            optsys.source_offset_theta = np.rad2deg(np.arctan2(-offx, offy))
            _log.debug("Source offset from X,Y = ({}, {}) is (r,theta) = {},{}".format(
                offx, offy, optsys.source_offset_r, optsys.source_offset_theta))
        if 'source_offset_r' in options:
            optsys.source_offset_r = options['source_offset_r']
        if 'source_offset_theta' in options:
            optsys.source_offset_theta = options['source_offset_theta']

        # Telescope entrance pupil
        pupil_optic = self._get_telescope_pupil_and_aberrations()
        optsys.add_pupil(pupil_optic)

        pupil_rms_wfe_nm = np.sqrt(np.mean(pupil_optic.opd[pupil_optic.amplitude == 1] ** 2)) * 1e9
        self._extra_keywords['TEL_WFE'] = (pupil_rms_wfe_nm, '[nm] Telescope pupil RMS wavefront error')
        if hasattr(pupil_optic, 'header_keywords'):
            self._extra_keywords.update(pupil_optic.header_keywords())

        self.pupil_radius = pupil_optic.pupil_diam / 2.0

        # add coord transform from entrance pupil to exit pupil
        optsys.add_inversion(axis='y', name='OTE exit pupil', hide=True)

        # add rotation at this point, if present - needs to be after the
        # exit pupil inversion.
        # Sign convention: Here we are rotating the *wavefront* so the sign is opposite the _rotation attribute,
        # which gives the V3IdlYangle for the detector rotation.
        if self._rotation is not None:
            optsys.add_rotation(-self._rotation, hide=True)
            optsys.planes[-1].wavefront_display_hint = 'intensity'

        # Allow instrument subclass to add field-dependent aberrations
        aberration_optic = self._get_aberrations()
        if aberration_optic is not None:
            optsys.add_pupil(aberration_optic)

            try:
                # Calculate SI WFE over just the OTE entrance pupil aperture,
                # though with a flip in the Y axis to account for entrance vs. exit pupil conventions
                exit_pupil_mask = pupil_optic.amplitude[::-1] == 1
                inst_rms_wfe_nm = np.sqrt(np.mean(aberration_optic.opd[exit_pupil_mask] ** 2)) * 1e9
                self._extra_keywords['SI_WFE'] = (inst_rms_wfe_nm, '[nm] instrument pupil RMS wavefront error')
            except (TypeError, IndexError):
                # Currently the above does not work for Roman, but fixing this is deferred to future work
                pass

            if hasattr(aberration_optic, 'header_keywords'):
                self._extra_keywords.update(aberration_optic.header_keywords())

        # ---- Add defocus if requested
        if 'defocus_waves' in options:
            defocus_waves = options['defocus_waves']
            defocus_wavelength = float(options['defocus_wavelength']) if 'defocus_wavelength' in options else 2.0e-6
            _log.info(f"Adding defocus of {defocus_waves:.3f} waves at {defocus_wavelength*1e6:.3f} microns" )
            lens = poppy.ThinLens(
                name='Defocus',
                nwaves=defocus_waves,
                reference_wavelength=defocus_wavelength,
                radius=self.pupil_radius
            )
            optsys.add_pupil(optic=lens)
            self._extra_keywords['DEFOCUS'] = (defocus_waves, '# of waves of defocus added')
            self._extra_keywords['DEFOC_WL'] = (defocus_wavelength, 'Wavelength reference for defocus added')

        # ---- add coronagraph or spectrograph optics if requested,
        # and possibly flag to invoke semi-analytic coronagraphic propagation

        # first error check for null strings, which should be considered like None
        if self.image_mask == "": self.image_mask = None
        if self.pupil_mask == "": self.pupil_mask = None

        if (self.image_mask is not None or self.pupil_mask is not None or
                'WL' in self.filter or  # special case handling for NIRCam WLP4 filter that is also a lens
                ('force_coron' in options and options['force_coron'])):
            _log.debug("Adding coronagraph/spectrograph optics...")
            optsys, trySAM, SAM_box_size = self._addAdditionalOptics(optsys, oversample=fft_oversample)
        else:
            trySAM = False

        # --- add the detector element.
        if fov_pixels is None:
            if not np.isscalar(fov_arcsec): fov_arcsec = np.asarray(fov_arcsec)  # cast to ndarray if 2D
            fov_pixels = np.round(fov_arcsec / self.pixelscale)
            if 'parity' in options:
                if options['parity'].lower() == 'odd' and np.remainder(fov_pixels, 2) == 0: fov_pixels += 1
                if options['parity'].lower() == 'even' and np.remainder(fov_pixels, 2) == 1: fov_pixels += 1
        else:
            pass

        optsys.add_detector(self.pixelscale, fov_pixels=fov_pixels,
                            oversample=detector_oversample,
                            name=self.name + " detector")

        # ---  invoke semi-analytic coronagraphic propagation
        if trySAM and not ('no_sam' in self.options and self.options['no_sam']):
            # if this flag is set, try switching to SemiAnalyticCoronagraph mode.
            _log.info("Trying to invoke switch to Semi-Analytic Coronagraphy algorithm")
            try:
                SAM_optsys = poppy.SemiAnalyticCoronagraph(optsys,
                                                           oversample=fft_oversample,
                                                           occulter_box=SAM_box_size)
                _log.info("SAC OK")
                return SAM_optsys
            except ValueError as err:
                _log.warning(
                    "Could not switch to Semi-Analytic Coronagraphy mode; invalid set of optical planes? "
                    "Using default propagation instead.")
                _log.warning(str(err))
                # _log.warn("ERROR ({0}): {1}".format(errno, strerror))
                pass

        return optsys

    def _get_telescope_pupil_and_aberrations(self):
        """return OpticalElement modeling wavefront aberrations for the telescope.

        See also get_aberrations for the SI aberrations.
        """

        # ---- set pupil OPD
        if isinstance(self.pupilopd, str):  # simple filename
            opd_map = self.pupilopd if os.path.exists(self.pupilopd) else \
                      os.path.join(self._datapath, "OPD", self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], str):
            # tuple with filename and slice
            opd_map = (self.pupilopd[0] if os.path.exists(self.pupilopd[0])
                       else os.path.join(self._datapath, "OPD", self.pupilopd[0]),
                       self.pupilopd[1])
        elif isinstance(self.pupilopd, (fits.HDUList, poppy.OpticalElement)):
            opd_map = self.pupilopd  # not a path per se but this works correctly to pass it to poppy
        elif self.pupilopd is None:
            opd_map = None
        else:
            raise TypeError("Not sure what to do with a pupilopd of that type:" + str(type(self.pupilopd)))

        # ---- set pupil intensity
        if self.pupil is None:
            raise RuntimeError("The pupil shape must be specified in the "
                               "instrument class or by setting self.pupil")
        if isinstance(self.pupil, poppy.OpticalElement):
            # supply to POPPY as-is
            pupil_optic = self.pupil
        else:
            # wrap in an optic and supply to POPPY
            if isinstance(self.pupil, str):  # simple filename
                if os.path.exists(self.pupil):
                    pupil_transmission = self.pupil
                else:
                    pupil_transmission = os.path.join(
                        self._WebbPSF_basepath,
                        self.pupil
                    )
            elif isinstance(self.pupil, fits.HDUList):
                # POPPY can use self.pupil as-is
                pupil_transmission = self.pupil
            else:
                raise TypeError("Not sure what to do with a pupil of "
                                "that type: {}".format(type(self.pupil)))
            # ---- apply pupil intensity and OPD to the optical model
            pupil_optic = poppy.FITSOpticalElement(
                name='{} Entrance Pupil'.format(self.telescope),
                transmission=pupil_transmission,
                opd=opd_map,
                planetype=poppy.poppy_core.PlaneType.pupil
                # rotation=self._rotation
            )
        return pupil_optic

    def _addAdditionalOptics(self, optsys, oversample=2):
        """Add instrument-internal optics to an optical system, typically coronagraphic or
        spectrographic in nature. This method must be provided by derived instrument classes.

        Returns
        --------
        optsys : OpticalSystem
            modified to add coronagraph optics
        useSAM : bool
            flag that, after adding the Detector, the whole thing should be converted to
            a SemiAnalyticCoronagraph model
        SAM_box_size : float
            size of box that entirely encloses the image plane occulter, in arcsec.

        """
        raise NotImplementedError("needs to be subclassed.")

    def _get_synphot_bandpass(self, filtername):
        """ Return a synphot.spectrum.SpectralElement object for the given desired band.

        By subclassing this, you can define whatever custom bandpasses are appropriate for
        your instrument
        """

        # use our local throughput files and create a synphot
        # transmission object.
        try:
            filter_info = self._filters[filtername]
        except KeyError:
            msg = "Couldn't find filter '{}' for {} in PySynphot or local throughput files"
            raise RuntimeError(msg.format(filtername, self.name))

        # The existing FITS files all have wavelength in ANGSTROMS since that is
        # the pysynphot convention...
        filterfits = fits.open(filter_info.filename)
        waveunit = filterfits[1].header.get('WAVEUNIT')
        if waveunit is None:
            _log.warning('The supplied file, {}, does not have a WAVEUNIT keyword. Assuming it '
                         'is Angstroms.'.format(filter_info.filename))
            waveunit = 'angstrom'

        filterdata = filterfits[1].data
        try:
            band = synphot.SpectralElement(synphot.models.Empirical1D, points=filterdata.WAVELENGTH,
                                               lookup_table=filterdata.THROUGHPUT, keep_neg=False)

        except AttributeError:
            raise ValueError("The supplied file, %s, does not appear to be a FITS table "
                             "with WAVELENGTH and THROUGHPUT columns." % filter_info.filename)

        filterfits.close()
        return band

    def psf_grid(self, num_psfs=16, all_detectors=True, save=False,
                 outdir=None, outfile=None, overwrite=True, verbose=True,
                 use_detsampled_psf=False, single_psf_centered=True, **kwargs):
        """
        Create a PSF library in the form of a grid of PSFs across the detector
        based on the specified instrument, filter, and detector. The output
        GriddedPSFModel object will contain a 3D array with axes [i, y, x]
        where i is the PSF position on the detector grid and (y,x) is the 2D
        PSF.

        Parameters
        ----------
        num_psfs : int
            The total number of fiducial PSFs to be created and saved in the files.
            This number must be a square number. Default is 16.
            E.g. num_psfs = 16 will create a 4x4 grid of fiducial PSFs.
        all_detectors : bool
            If True, run all detectors for the instrument. If False, run for
            the detector set in the instance. Default is True
        save : bool
            True/False boolean if you want to save your file. Default is False.
        outdir : str
            If "save" keyword is set to True, your file will be saved in the
            specified directory. Default of None will save it in the current
            directory
        outfile : str
            If "save" keyword is set to True, your file will be saved as
            {outfile}_det.fits. Default of None will save it as
            instr_det_filt_fovp#_samp#_npsf#.fits
        overwrite : bool
            True/False boolean to overwrite the output file if it already exists.
            Default is True.
        verbose : bool
            True/False boolean to print status updates. Default is True.
        use_detsampled_psf : bool
            If True, the grid of PSFs returned will be detector sampled (made
            by binning down the oversampled PSF). If False, the PSFs will be
            oversampled by the factor defined by the
            oversample/detector_oversample/fft_oversample keywords. Default is False.
            This is rarely needed - if uncertain, leave this alone.
        single_psf_centered : bool
            If num_psfs is set to 1, this defines where that psf is located.
            If True it will be the center of the detector, if False it will
            be the location defined in the WebbPSF attribute detector_position
            (reminder - detector_position is (x,y)). Default is True
            This is also rarely needed.
        **kwargs
            Any extra arguments to pass the WebbPSF calc_psf() method call.

        Returns
        -------
        gridmodel : photutils GriddedPSFModel object or list of objects
            Returns a GriddedPSFModel object or a list of objects if more than one
            configuration is specified (1 per instrument, detector, and filter)
            User also has the option to save the grid as a fits.HDUlist object.

        Use
        ----
        nir = webbpsf.NIRCam()
        nir.filter = "F090W"
        list_of_grids = nir.psf_grid(all_detectors=True, num_psfs=4)

        wfi = webbpsf.WFI()
        wfi.filter = "Z087"
        wfi.detector = "SCA02"
        grid = wfi.psf_grid(all_detectors=False, oversample=5, fov_pixels=101)

        """

        # Keywords that could be set before the method call
        filt = self.filter

        if all_detectors is True:
            detectors = "all"
        else:
            detectors = self.detector

        if single_psf_centered is True:
            psf_location = (int((self._detector_npixels - 1) / 2), int((self._detector_npixels - 1) / 2))  # center pt
        else:
            psf_location = self.detector_position[::-1]  # (y,x)

        # add_distortion keyword is not implemented for WFI Class
        if self.name == "WFI" and "add_distortion" not in kwargs:
            kwargs["add_distortion"] = False
        elif self.name == "WFI" and kwargs["add_distortion"] == True:
            raise NotImplementedError("Geometric distortions are not implemented in WebbPSF for WFI Instrument. "
                                      "The add_distortion keyword must be set to False for this case.")

        # Call CreatePSFLibrary class
        inst = gridded_library.CreatePSFLibrary(instrument=self, filter_name=filt, detectors=detectors,
                                                num_psfs=num_psfs, psf_location=psf_location,
                                                use_detsampled_psf=use_detsampled_psf, save=save,
                                                outdir=outdir, filename=outfile, overwrite=overwrite,
                                                verbose=verbose, **kwargs)
        gridmodel = inst.create_grid()

        return gridmodel


#######  JWInstrument classes  #####

@utils.combine_docstrings
class JWInstrument(SpaceTelescopeInstrument):
    """ Superclass for all JWST instruments

    Notable attributes
    -------------------

    telescope : name of telescope
    pupilopd : filename or FITS file object

    include_si_wfe : boolean (default: True)
        Should SI internal WFE be included in models? Requires
        the presence of ``si_zernikes_isim_cv3.fits`` in the
        ``WEBBPSF_PATH``.
    """
    telescope = "JWST"
    pupilopd = None
    """Filename *or* fits.HDUList for JWST pupil OPD.

    This can be either a full absolute filename, or a relative name in which
    case it is assumed to be within the instrument's `data/OPDs/` directory,
    or an actual fits.HDUList object corresponding to such a file. If the file
    contains a datacube, you may set this to a tuple (filename, slice) to
    select a given slice, or else the first slice will be used."""

    def __init__(self, *args, **kwargs):
        super(JWInstrument, self).__init__(*args, **kwargs)

        self.siaf = pysiaf.Siaf(self.name)

        opd_path = os.path.join(self._datapath, 'OPD')
        self.opd_list = []
        for filename in glob.glob(os.path.join(opd_path, 'OPD*.fits*')):
            self.opd_list.append(os.path.basename(os.path.abspath(filename)))
        for filename in glob.glob(os.path.join(self._WebbPSF_basepath, 'JWST_OTE_OPD*.fits*')):
            self.opd_list.append(os.path.basename(os.path.abspath(filename)))

        if not len(self.opd_list) > 0:
            raise RuntimeError("No pupil OPD files found for {name} in {path}".format(name=self.name, path=opd_path))

        self.opd_list.sort()
        self.pupilopd = 'JWST_OTE_OPD_cycle1_example_2022-07-30.fits'    # Default is now an on-orbit measured example OPD

        self.pupil = os.path.abspath(os.path.join(
            self._WebbPSF_basepath,
            "jwst_pupil_RevW_npix1024.fits.gz"
        ))
        "Filename *or* fits.HDUList for JWST pupil mask. Usually there is no need to change this."

        self._aperturename = None
        self._detector = None

        # where is the source on the detector, in 'Science frame' pixels?
        self.detector_position = (1024, 1024)

        self.include_si_wfe = True
        self.include_ote_field_dependence = True  # Note, this will be implicitly ignored if pupilopd=None
        """Should calculations include the Science Instrument internal WFE?"""
        self.options['jitter'] = 'gaussian'
        self.options['jitter_sigma'] = constants.JWST_TYPICAL_LOS_JITTER_PER_AXIS

        # class name to use for SI internal WFE, which can be overridden in subclasses
        self._si_wfe_class = optics.WebbFieldDependentAberration

    def _get_default_fov(self):
        """ Return default FOV in arcseconds """
        return 5  # default for all NIR instruments

    def get_optical_system(self, fft_oversample=2, detector_oversample=None, fov_arcsec=2, fov_pixels=None, options=None):
        # invoke superclass version of this
        # then add a few display tweaks
        optsys = SpaceTelescopeInstrument.get_optical_system(self,
                                                              fft_oversample=fft_oversample,
                                                              detector_oversample=detector_oversample,
                                                              fov_arcsec=fov_arcsec, fov_pixels=fov_pixels,
                                                              options=options)
        # If the OTE model in the entrance pupil is a plain FITSOpticalElement, cast it to the linear model class
        if not isinstance(optsys.planes[0], opds.OTE_Linear_Model_WSS):
            lom_ote = opds.OTE_Linear_Model_WSS()
            # FIXME seems like some code is missing here...? But in practice this code path
            # never gets executed due to the _get_telescope_pupil_and_aberrations() function doing the right thing.
            lom_ote

        optsys.planes[0].display_annotate = utils.annotate_ote_pupil_coords
        return optsys

    def _get_aberrations(self):
        """ return OpticalElement modeling wavefront aberrations for a given instrument,
        including field dependence based on a lookup table of Zernike coefficients derived from
        ISIM cryovac test data.
        """
        if not self.include_si_wfe:
            return None

        optic = self._si_wfe_class(self)
        return optic

    def get_opd_file_full_path(self, opdfilename=None):
        """Return full path to the named OPD file.

        The OPD may be:
         - a local or absolute path,
         - or relative implicitly within an SI directory, e.g. $WEBBPSF_PATH/NIRCam/OPD
         - or relative implicitly within $WEBBPSF_PATH

        This function handles filling in the implicit path in the latter cases.
        """

        if opdfilename is None:
            opdfilename = self.pupilopd

        if os.path.exists(opdfilename):
            return opdfilename
        elif self.name in opdfilename:
            return os.path.join(self._datapath, "OPD", opdfilename)
        else:
            return os.path.join(self._WebbPSF_basepath, opdfilename)

    def _get_telescope_pupil_and_aberrations(self):
        """return OpticalElement modeling wavefront aberrations for the telescope.

        This is nearly identical to the version of this function in SpaceTelescopeInstrument, differing only at the
        very end. Here, we load the selected OPD file from disk into an instance of opds.OTE_Linear_Model_WSS if possible.
        It falls back to a plain FITSOpticalElement for nonstandard sizes of input pupil, since the linear model is not
        yet generalized to work on arbitrary sizes of pupil other than 1024 pixels.

        See also get_aberrations for the SI aberrations.
        """

        # ---- set pupil OPD



        opd_index = None  # default assumption: OPD file is not a datacube
        if isinstance(self.pupilopd, str):  # simple filename
            opd_map = self.get_opd_file_full_path(self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], str):
            # tuple with filename and slice, for a datacube
            opd_map = self.get_opd_file_full_path(self.pupilopd[0])
            opd_index = self.pupilopd[1]
        elif isinstance(self.pupilopd, (fits.HDUList, poppy.OpticalElement)):
            opd_map = self.pupilopd  # not a path per se but this works correctly to pass it to poppy
        elif self.pupilopd is None:
            opd_map = None
        else:
            raise TypeError("Not sure what to do with a pupilopd of that type:" + str(type(self.pupilopd)))

        # ---- set pupil intensity
        if self.pupil is None:
            raise RuntimeError("The pupil shape must be specified in the "
                               "instrument class or by setting self.pupil")
        if isinstance(self.pupil, poppy.OpticalElement):
            # supply to POPPY as-is
            pupil_optic = self.pupil
        else:
            # wrap in an optic and supply to POPPY
            if isinstance(self.pupil, str):  # simple filename
                if os.path.exists(self.pupil):
                    pupil_transmission = self.pupil
                else:
                    pupil_transmission = os.path.join(
                        self._WebbPSF_basepath,
                        self.pupil
                    )
                # Get npix from pupil_transmission
                npix = int(pupil_transmission.split('npix')[-1].split('.')[0])
            elif isinstance(self.pupil, fits.HDUList):
                # POPPY can use self.pupil as-is
                pupil_transmission = self.pupil
                # Get npix from the shape of the data
                npix = self.pupil[0].data.shape[0]
            else:
                raise TypeError("Not sure what to do with a pupil of "
                                "that type: {}".format(type(self.pupil)))

            # ---- apply pupil intensity and OPD to the optical model
            pupil_optic = opds.OTE_Linear_Model_WSS(
                name='{} Entrance Pupil'.format(self.telescope),
                transmission=pupil_transmission,
                opd=opd_map,
                opd_index=opd_index,
                v2v3=self._tel_coords(), npix=npix,
                include_nominal_field_dependence=self.include_ote_field_dependence
            )

        return pupil_optic


    @SpaceTelescopeInstrument.aperturename.setter
    def aperturename(self, value):
        """Set SIAF aperture name to new value, with validation
        """
        # Explicitly update detector reference coordinates to the default for the new selected aperture,
        # otherwise old coordinates can persist under certain circumstances

        try:
            ap = self.siaf[value]
        except KeyError:
            raise ValueError(f'Aperture name {value} not a valid SIAF aperture name for {self.name}')

        if self.detector not in value:
            raise ValueError(f'Aperture name {value} does not match currently selected detector {self.detector}. '
                             f'Change detector attribute first, then set desired aperture.')

        # Only update if new value is different
        if self._aperturename != value:
            self._aperturename = value
            # Update detector reference coordinates
            self.detector_position = (ap.XSciRef, ap.YSciRef)

            # Update DetectorGeometry class
            self._detector_geom_info = DetectorGeometry(self.siaf, self._aperturename)
            _log.info(f"{self.name} SIAF aperture name updated to {self._aperturename}")

    def _tel_coords(self):
        """ Convert from science frame coordinates to telescope frame coordinates using
        SIAF transformations. Returns (V2, V3) tuple, in arcminutes.

        Note that the astropy.units framework is used to return the result as a
        dimensional Quantity.
        """

        return self._detector_geom_info.pix2angle(self.detector_position[0], self.detector_position[1])

    def _xan_yan_coords(self):
        """ Convert from detector pixel coordinates to the XAN, YAN coordinate system
        which was used for much of ISIM optical testing. The origin of XAN, YAN is
        centered at the master chief ray, which passes through the ISIM focal plane
        between the NIRCam A3 and B4 detectors. The sign of YAN is flipped relative to V3.
        """
        coords = self._tel_coords()
        # XAN is the same as V2, therefore no change to first element
        # YAN is opposite direction as V3, and offset by 468 arcseconds
        coords[1] = -coords[1] - 468 * units.arcsec
        return coords

    def set_position_from_aperture_name(self, aperture_name):
        """ Set the simulated center point of the array based on a named SIAF aperture.
        This will adjust the detector and detector position attributes.
        """
        try:
            ap = self.siaf[aperture_name]

            # setting the detector must happen -before- we set the position
            detname = aperture_name.split('_')[0]
            self.detector = detname  # As a side effect this auto reloads SIAF info, see detector.setter

            self.aperturename = aperture_name

            if self.name != 'NIRSpec' and ap.AperType != 'SLIT':
                # Regular imaging apertures, so we can just look up the reference coords directly
                self.detector_position = (ap.XSciRef, ap.YSciRef)  # set this AFTER the SIAF reload
            else:
                # NIRSpec slit apertures need some separate handling, since they don't map directly to detector pixels
                ref_in_tel = ap.V2Ref, ap.V3Ref
                nrs_full_aperture = self.siaf[aperture_name.split('_')[0]+"_FULL"]
                ref_in_sci = nrs_full_aperture.tel_to_sci(*ref_in_tel)
                self.detector_position = ref_in_sci

            _log.debug("From {} set det. pos. to {} {}".format(aperture_name, detname, self.detector_position))

        except KeyError:
            raise ValueError("Not a valid aperture name for {}: {}".format(self.name, aperture_name))

    def _get_fits_header(self, result, options):
        """ populate FITS Header keywords """
        super(JWInstrument, self)._get_fits_header(result, options)

        # Add JWST-specific V2,V3 focal plane coordinate system.
        v2v3pos = self._tel_coords()
        result[0].header.insert("DET_Y", ('DET_V2', v2v3pos[0].value,
                                          "[arcmin] Det. pos. in telescope V2,V3 coord sys"), after=True)
        result[0].header.insert("DET_V2", ('DET_V3', v2v3pos[1].value,
                                           "[arcmin] Det. pos. in telescope V2,V3 coord sys"), after=True)
        result[0].header["APERNAME"] = (self._aperturename, "SIAF aperture name")

    def calc_psf(self, outfile=None, source=None, nlambda=None, monochromatic=None,
                 fov_arcsec=None, fov_pixels=None, oversample=None, detector_oversample=None, fft_oversample=None,
                 overwrite=True, display=False, save_intermediates=False, return_intermediates=False,
                 normalize='first', add_distortion=True, crop_psf=True):
        """
        Compute a PSF

        Parameters
        ----------
        add_distortion : bool
            If True, will add 2 new extensions to the PSF HDUlist object. The 2nd extension
            will be a distorted version of the over-sampled PSF and the 3rd extension will
            be a distorted version of the detector-sampled PSF.
        crop_psf : bool
            If True, when the PSF is rotated to match the detector's rotation in the focal
            plane, the PSF will be cropped so the shape of the distorted PSF will match it's
            undistorted counterpart. This will only be used for NIRCam, NIRISS, and FGS PSFs.

        """

        # Save new keywords to the options dictionary
        self.options['add_distortion'] = add_distortion
        self.options['crop_psf'] = crop_psf

        # UPDATE THE OPD V2V3 BASED ON DETECTOR POSITION, IN ORDER TO CALCULATE SM FIELD-DEPENDENT WFE.
        # SEE opds._apply_sm_field_dependence_model()
        #
        # v2v3 attribute exists only if using the linear model, so check first:
        if hasattr(self.pupil, 'v2v3'):
            if (self.pupil.v2v3 is None) or (not (self.pupil.v2v3 == self._tel_coords().to(units.arcsec)).all()):
                self.pupil.v2v3 = self._tel_coords().to(units.arcsec)
                self.pupil.update_opd()

        # Run poppy calc_psf
        psf = SpaceTelescopeInstrument.calc_psf(self, outfile=outfile, source=source, nlambda=nlambda,
                                                monochromatic=monochromatic, fov_arcsec=fov_arcsec,
                                                fov_pixels=fov_pixels, oversample=oversample,
                                                detector_oversample=detector_oversample, fft_oversample=fft_oversample,
                                                overwrite=overwrite, display=display,
                                                save_intermediates=save_intermediates,
                                                return_intermediates=return_intermediates, normalize=normalize)

        return psf

    def _calc_psf_format_output(self, result, options):
        """
        Add distortion to the created 1-extension PSF

        Apply desired formatting to output file:
                 - rebin to detector pixel scale if desired
                 - set up FITS extensions if desired
                 - output either the oversampled, rebinned, or both
        Which image(s) get output depends on the value of the options['output_mode']
        parameter. It may be set to 'Oversampled image' to output just the oversampled image,
        'Detector sampled image' to output just the image binned down onto detector pixels, or
        'Both as FITS extensions' to output the oversampled image as primary HDU and the
        rebinned image as the first image extension. For convenience, the option can be set
        to just 'oversampled', 'detector', or 'both'.

        Modifies the 'result' HDUList object.

        """
        # Pull values from options dictionary
        add_distortion = options.get('add_distortion', True)
        crop_psf = options.get('crop_psf', True)

        # Add distortion if set in calc_psf
        if add_distortion:
            _log.debug("Adding PSF distortion(s)")
            if self.image_mask == "LRS slit" and self.pupil_mask == "P750L":
                raise NotImplementedError("Distortion is not implemented yet for MIRI LRS mode.")

            # Set up new extensions to add distortion to:
            n_exts = len(result)
            for ext in np.arange(n_exts):
                hdu_new = fits.ImageHDU(result[ext].data, result[ext].header)  # these will be the PSFs that are edited
                result.append(hdu_new)
                ext_new = ext + n_exts
                result[ext_new].header["EXTNAME"] = result[ext].header["EXTNAME"][0:4] + "DIST"  # change extension name
                _log.debug("Appending new extension {} with EXTNAME = {}".format(ext_new, result[ext_new].header["EXTNAME"]))

            # Apply distortions based on the instrument
            if self.name in ["NIRCam", "NIRISS", "FGS"]:
                # Apply distortion effects: Rotation and optical distortion
                _log.debug("NIRCam/NIRISS/FGS: Adding rotation and optical distortion")
                psf_rotated = distortion.apply_rotation(result, crop=crop_psf)  # apply rotation
                psf_distorted = distortion.apply_distortion(psf_rotated)  # apply siaf distortion model
            elif self.name == "MIRI":
                # Apply distortion effects to MIRI psf: Distortion and MIRI Scattering
                _log.debug("MIRI: Adding optical distortion and Si:As detector internal scattering")
                psf_siaf = distortion.apply_distortion(result)  # apply siaf distortion
                psf_distorted = distortion.apply_miri_scattering(psf_siaf)  # apply scattering effect
            elif self.name == "NIRSpec":
                # Apply distortion effects to NIRSpec psf: Distortion only
                _log.debug("NIRSpec: Adding optical distortion")
                psf_distorted = distortion.apply_distortion(result)  # apply siaf distortion model

            # Edit the variable to match if input didn't request distortion
            # (cannot set result = psf_distorted due to return method)
            [result.append(fits.ImageHDU()) for i in np.arange(len(psf_distorted) - len(result))]
            for ext in np.arange(len(psf_distorted)): result[ext] = psf_distorted[ext]

        # Rewrite result variable based on output_mode set:
        SpaceTelescopeInstrument._calc_psf_format_output(self, result, options)

    def interpolate_was_opd(self, array, newdim):
        """ Interpolates an input 2D  array to any given size.

        Parameters
        ----------
        array: float
             input array to interpolate
        newdim: int
             new size of the 2D square array (newdim x newdim)

        Returns
        ---------
        newopd: new array interpolated to (newdim x newdim)

        """

        dim = array.shape[0]

        xmax, ymax = dim / 2, dim / 2
        x = np.arange(-xmax, xmax, 1)
        y = np.arange(-ymax, ymax, 1)
        X, Y = np.meshgrid(x, y)

        interp_spline = scipy.interpolate.RectBivariateSpline(y, x, array)

        dx, dy = float(dim) / float(newdim), float(dim) / float(newdim)

        x2 = np.arange(-xmax, xmax, dx)
        y2 = np.arange(-ymax, ymax, dy)
        X2, Y2 = np.meshgrid(x2, y2)
        newopd = interp_spline(y2, x2)
        newopd = np.reshape(newopd, (1, newdim, newdim))

        return newopd


    def _get_pupil_shift(self):
        """ Return a tuple of pupil shifts, for passing to OpticalElement constructors
        This is a minor utility function that gets used in most of the subclass optical
        system construction.

        For historical reasons, the pupil_shift_x and pupil_shift_y options are expressed
        in fractions of the pupil. The parameters to poppy should now be expressed in
        meters of shift. So the translation of that happens here.

        Returns
        -------
        shift_x, shift_y : floats or Nones
            Pupil shifts, expressed in meters.

        """
        if ('pupil_shift_x' in self.options and self.options['pupil_shift_x'] != 0) or \
                ('pupil_shift_y' in self.options and self.options['pupil_shift_y'] != 0):

            from .constants import JWST_CIRCUMSCRIBED_DIAMETER
            # missing values are treated as 0's
            shift_x = self.options.get('pupil_shift_x', 0)
            shift_y = self.options.get('pupil_shift_y', 0)
            # nones are likewise treated as 0's
            if shift_x is None: shift_x = 0
            if shift_y is None: shift_y = 0
            # Apply pupil scale
            shift_x *= JWST_CIRCUMSCRIBED_DIAMETER
            shift_y *= JWST_CIRCUMSCRIBED_DIAMETER
            _log.info("Setting Lyot pupil shift to ({}, {})".format(shift_x,shift_y))
        else:
            shift_x, shift_y = None, None
        return shift_x, shift_y


    def _apply_jitter(self,  result, local_options=None):
        """ Modify a PSF to account for the blurring effects of image jitter.
        Parameter arguments are taken from the options dictionary.

        This adds options to model JWST coarse point ("PCS=Coarse") under
        two sets of assumptions:
            "PCS=Coarse": 67 mas Gaussian jitter, as advised by Nelan & Maghami based on
                          detailed sims of observatory performance in coarse point mode.
            "PCS=Coarse_Like_ITM": Attempt to replicate same assumptions as in Ball's ITM tool.
                          This includes 200 mas sigma Gaussian jitter, plus a linear drift of
                          400 mas per exposure.

        Other types of jitter, in particular plain Gaussian jitter, are implemented by the
        superclass version of this function, in poppy.Instrument.

        Parameters
        -----------
        result : fits.HDUList
            HDU list containing a point spread function
        local_options : dict, optional
            Options dictionary. If not present, options will be taken from self.options.


        The image in the 'result' HDUlist will be modified by this function.
        """
        if local_options is None:
            local_options = self.options
        if 'jitter' not in local_options:
            result[0].header['JITRTYPE'] = ('None', 'Type of jitter applied')
            return

        _log.info("Calculating jitter using " + str(local_options['jitter']))

        def _linear_smear(smear_length, image):
            # Helper function, used below
            smear_length_pix = int(np.round(smear_length /  result[0].header['PIXELSCL']))
            if smear_length_pix % 2 ==0:
                smear_length_pix += 1   # Astropy convolution requires odd sized kernels only

            smear_model  = np.identity(smear_length_pix)
            _log.info("Jitter: Convolving with linear smear of {0:.3f} arcsec; {1:d} pixels".format(smear_length, smear_length_pix))
            kern = astropy.convolution.kernels.CustomKernel(smear_model)
            return astropy.convolution.convolve_fft(image, kern, allow_huge=True)

        if local_options['jitter'] is None:
            return
        elif local_options['jitter'].lower() == 'gaussian':
            # Regular version in poppy
            return super()._apply_jitter(result, local_options=local_options)
        elif local_options['jitter'].lower() == 'linear':
            # Drift by 0.12 arcsec (1 mas/second for 2 minutes)

            smear_length = 0.12 # arcsec

            out = _linear_smear(smear_length, result[0].data)
            result[0].header['JITRTYPE'] = ('Linear smear / drift', 'Type of jitter applied')
            result[0].header['JITSMEAR'] = (smear_length, 'Linear smear [arcsec]')

        elif local_options['jitter'].lower() == 'pcs=coarse':
            # JWST coarse point, current best estimate based on high fidelity monte carlo sims by Peiman Maghami

            cp_case = local_options.get('jitter_coarse_model_case', 2)      # Coarse pointing model case, 1 or 2
            exp_duration = local_options.get('exp_duration', 75)     # Duration in seconds
            exp_start_time = local_options.get('exp_start_time', 0)  # Start time in seconds

            offset, kernel = opds.get_coarse_blur_parameters(exp_start_time, exp_duration, result[0].header['PIXELSCL'], case=cp_case)

            kern = astropy.convolution.kernels.CustomKernel(kernel)
            out = astropy.convolution.convolve_fft(result[0].data, kern, allow_huge=True)

            result[0].header['JITRTYPE'] = ('PCS Coarse, high fidelity MC model results', 'Type of jitter applied')
            result[0].header['JITRCASE'] = (cp_case, 'PCS Coarse mode: Monte Carlo model case used')
            result[0].header['JITR_T0'] =  (exp_start_time, 'PCS Coarse mode: sim exposure start time [s]')
            result[0].header['JITRTEXP'] = (exp_duration, 'PCS Coarse mode: sim exposure duration [s]')
            result[0].header['JITRCPV2'] = (offset[0], "Coarse pointing offset in V2 [arcsec]")
            result[0].header['JITRCPV3'] = (offset[1], "Coarse pointing offset in V3 [arcsec]")

        elif local_options['jitter'].lower() == 'pcs=coarse_like_itm':
            # JWST coarse point, assumptions in ITM
            # Acton says:
            #  it is actually 0.4 for a boresight error, 0.4 smear, and 0.2 jitter. Boresight error is a random term for image placement, smear is mostly a linear uniform blur, and jitter is gaussian.

            # First we do the fast jitter part
            local_options['jitter_sigma'] = 0.2
            import scipy.ndimage

            sigma = local_options.get('jitter_sigma')

            # that will be in arcseconds, we need to convert to pixels:
            _log.info("Jitter: Convolving with Gaussian with sigma={0:.3f} arcsec".format(sigma))
            out = scipy.ndimage.gaussian_filter(result[0].data, sigma / result[0].header['PIXELSCL'])

            # Now we'll do the linear jitter part
            smear_length = 0.4 # arcsec
            out = _linear_smear(smear_length, out)

            result[0].header['JITRTYPE'] = ('PCS Coarse, like ITM', 'Type of jitter applied')
            result[0].header['JITRSIGM'] = (sigma, 'Gaussian sigma for jitter, per axis [arcsec]')
            result[0].header['JITSMEAR'] = (smear_length, 'Linear smear [arcsec]')

        elif local_options['jitter'].lower() == 'custom':
            # User-supplied arbitrary PSF convolution kernel

            if ('jitter_kernel' not in local_options) or (not local_options['jitter_kernel'].ndim==2):
                raise ValueError("You must supply an .options['jitter_kernel'] 2D array to use the custom jitter option")
            _log.info("Jitter: Convolving with user-supplied custom convolution kernel")
            kern = astropy.convolution.kernels.CustomKernel(local_options['jitter_kernel'])
            out = astropy.convolution.convolve_fft(result[0].data, kern, allow_huge=True)

            result[0].header['JITRTYPE'] = ('Custom jitter kernel', 'Type of jitter applied')

        else:
            raise ValueError('Unknown jitter option value: ' + local_options['jitter'])

        peak = result[0].data.max()
        newpeak = out.max()
        strehl = newpeak / peak  # not really the whole Strehl ratio, just the part due to jitter
        _log.info("        resulting image peak drops to {0:.3f} of its previous value".format(strehl))
        result[0].header['JITRSTRL'] = (strehl, 'Strehl reduction from jitter ')

        result[0].data = out

    def get_wfe(self, kind='si', wavelength=2e-6, plot=False):
        """Extract and return one component plane of the optical model for this instrument

        This is a utility function for convenience, making it easier to access and plot various OPD maps.
        It doesn't do anything unique which can't be done otherwise, and in particular this isn't used at all
        as part of the optical propagation calculations.

        Note, all WFE terms are returned in OTE entrance pupil orientation (i.e. as if you were in front
        of the OTE and looking at it), regardless of pupil flips and orientations in the optical propagation.

        Parameters
        ----------
        kind : string
            A type of WFE. Must be one of "SI", "OTE", "OTE_field_dep", or other values TBD.
            Case insensitive.
        plot : bool
            Make a quick plot of this WFE. Not very flexible or scriptable but useful for some interactive checks
        """
        osys = self.get_optical_system()
        wave = osys.input_wavefront(wavelength)
        ote = osys.planes[0]

        if kind.lower() =='total':
            # recursively get total OPD including SI plus OTE
            opd = self.get_wfe('ote') + self.get_wfe('si')
        elif kind.lower()=='si':
            aberration = self._get_aberrations()
            opd = aberration.get_opd(wave)
            if self.name.lower()=='nirspec':
                # For NIRSpec, the WFE is normally allocated to 1/3 before the MSA and 2/3 after the MSA.
                # The call to get_aberrations above just returns the foreoptics portion.
                # Multiply by 3x to get the total instrumental WFE.
                opd *= 3
            # Flip vertically to match OTE entrance pupil orientation
            opd = np.flipud(opd)
        elif kind.lower() == 'ote': # OTE *total* WFE including all terms
            opd = ote.get_opd(wave).copy()
            aperture = ote.get_transmission(wave)
            opd *= (aperture != 0)  # mask out to zero the global zernikes outside the aperture

        elif kind.lower() == 'ote_global': # OTE *global* WFE only, i.e. WFE common to all field points
            # This is done recursively, since that's a convenient way to code this up
            opd_ote_total = self.get_wfe('ote')
            opd_ote_fd = self.get_wfe('ote_field_dep')
            return opd_ote_total - opd_ote_fd
        elif kind.lower() == 'ote_field_dep':  # OTE field dependent variations
            wfe_ote_field_dep_nominal = ote._get_field_dependence_nominal_ote(ote.v2v3)
            wfe_ote_field_dep_mimf = ote._get_field_dependence_secondary_mirror(ote.v2v3)
            wfe_ote_field_dep = wfe_ote_field_dep_nominal + wfe_ote_field_dep_mimf
            aperture = ote.get_transmission(wave)
            opd = wfe_ote_field_dep * (aperture != 0)  # mask out to zero the global zernikes outside the aperture

        elif kind.lower() == 'ote_thermal_distortion':  # OTE temporal variations from backplane thermal distortion
            raise NotImplementedError(f"Not yet implemented: {kind}")
        else:
            raise NotImplementedError(f"Not a known kind of WFE: {kind}")

        if plot:
            import matplotlib, matplotlib.pyplot as plt
            plt.imshow(opd, vmin=-5e-7, vmax=5e-7, cmap=matplotlib.cm.RdBu_r, origin='lower')
            plt.title(kind+" WFE")
            mask = ote.get_transmission(wave) !=0
            plt.xlabel(f"RMS: {utils.rms(opd, mask)*1e9:.2f} nm")
            plt.colorbar(label='WFE [m]')

        return opd


    def visualize_wfe_budget(self, slew_delta_time=14*units.day, slew_case='EOL', ptt_only=False, verbose=True):
        """Display a visual WFE budget showing the various terms that sum into the overall WFE for a given instrument

        Compares a WebbPSF instrument instance with the JWST optical budget for that instrument

        Parameters
        ----------
        inst : webbpsf.JWInstrument
            A JWST instrument instance
        slew_delta_time : astropy.Quantity time
            Time duration for thermal slew model
        slew_case : basestring
            'BOL' or 'EOL' for beginning of life or end of life thermal slew model. EOL is about 3x higher amplitude
        ptt_only : bool
            When decomposing wavefront into controllable modes, use a PTT-only basis? The default is to use all
            controllable pose modes. (This is mostly a leftover debug option at this point, not likely useful in general)
        verbose : bool
            Be more verbose
        """
        import webbpsf.optical_budget
        webbpsf.optical_budget.visualize_wfe_budget(self,
                                                    slew_delta_time=slew_delta_time, slew_case=slew_case,
                                                    ptt_only=ptt_only, verbose=verbose)

    def load_wss_opd(self, filename, output_path = None, backout_si_wfe=True, verbose=True, plot=False, save_ote_wfe=False):
        """Load an OPD produced by the JWST WSS into this instrument instance, specified by filename

        This includes:
            - If necessary, downloading that OPD from MAST. Downloaded files are cached in $WEBBPSF_PATH/MAST_JWST_WSS_OPDs
            - calling `import_wss_opd` to load the OPD from the FITS file and
               perform some necessary format conversions
            - Subtract off the instrument WFE for the field point used in wavefront sensing, to get an
               OTE-only wavefront. WebbPSF will separately add back in the SI WFE for the appropriate
              field point, as usual.
            - Subtract off the modeled field dependence term in the OTE WFE for the sensing field point, to get
               an estimate of the OTE wavefront nominally at the master chief ray location (between the NIRCams).
               WebbPSF will automatically add back on top of this the OTE field dependent WFE for the appropriate
               field point. as usual.

        Parameters
        ----------
        filename : str
            Name of OPD file to load
        output_path : str
            Downloaded OPD are saved in this location. This option is convinient for STScI users using /grp/jwst/ote/webbpsf-data/.
            Default is $WEBBPSF_PATH/MAST_JWST_WSS_OPDs
        backout_si_wfe : bool
            Subtract model for science instrument WFE at the sensing field point? Generally this should be true
            which is the default.
        plot : bool
            Generate informative plots showing WFE, including the backout steps. Only works if backout_si_wfe is True.
        save_ote_wfe : bool
            Save OTE-only WFE model? This is not needed for calculations in WebbPSF, but can be used to export
            OTE WFE models for use with other software. The file will be saved in the WEBBPSF_DATA_PATH directory
            and a message will be printed on screen with the filename.
            Note that the exported OPD file will give the OTE estimated total WFE at the selected Instrument's field
            point, not the OTE global at master chief ray, since it is the OTE WFE at the selected field point
            which is most of use for some other tool.

        """

        # If the provided filename doesn't exist on the local disk, try retrieving it from MAST
        # Note, this will automatically use cached versions downloaded previously, if present
        if not os.path.exists(filename):
            filename = webbpsf.mast_wss.mast_retrieve_opd(filename, output_path = output_path, verbose=verbose)

        if verbose:
            print(f"Importing and format-converting OPD from {filename}")
        opdhdu = webbpsf.mast_wss.import_wss_opd(filename)

        # Mask out any pixels in the OPD array which are outside the OTE pupil.
        # This is mostly cosmetic, and helps mask out some edge effects from the extrapolation + interpolation in
        # resizing the OPDs
        ote_pupil_mask = utils.get_pupil_mask() != 0
        opdhdu[0].data *= ote_pupil_mask

        #opdhdu[0].header['RMS_OBS'] = (webbpsf.utils.rms(opdhdu[0].data, mask=ote_pupil_mask)*1e9,
        #                               "[nm] RMS Observatory WFE (i.e. OTE+SI) at sensing field pt")

        if plot:
            import matplotlib, matplotlib.pyplot as plt
            fig, axes = plt.subplots(figsize=(16, 9), ncols=3, nrows=2)
            vm = 2e-7
            plot_kwargs = {'vmin':-vm, 'vmax':vm, 'cmap':matplotlib.cm.RdBu_r, 'origin':'lower'}
            axes[0,0].imshow(opdhdu[0].data.copy() * ote_pupil_mask, **plot_kwargs)
            axes[0,0].set_title(f"OPD from\n{os.path.basename(filename)}")
            axes[0,0].set_xlabel(f"RMS: {utils.rms(opdhdu[0].data*1e9, ote_pupil_mask):.2f} nm rms")

        if backout_si_wfe:
            if verbose: print("Backing out SI WFE and OTE field dependence at the WF sensing field point")

            # Check which field point was used for sensing
            sensing_apername = opdhdu[0].header['APERNAME']

            # Create a temporary instance of an instrument, for the sensng instrument and field point,
            # in order to model and extract the SI WFE and OTE field dep WFE at the sensing field point.

            sensing_inst = instrument(sensing_apername[0:3])
            sensing_inst.pupil = self.pupil   # handle the case if the user has selected a different NPIX other than the default 1024
            if sensing_inst.name == 'NRC':
                sensing_inst.filter = 'F212N'
                # TODO: optionally check for the edge case in which the sensing was done in F187N
                # note that there is a slight focus offset between the two wavelengths, due to NIRCam's refractive design
            # Set to the sensing aperture, and retrieve the OPD there
            sensing_inst.set_position_from_aperture_name(sensing_apername)
            # special case: for the main sensing point FP1, we use the official WAS target phase map, rather than the
            # WebbPSF-internal SI WFE model.
            was_targ_file = os.path.join(utils.get_webbpsf_data_path(), 'NIRCam', 'OPD', 'wss_target_phase_fp1.fits')
            if sensing_apername == 'NRCA3_FP1' and os.path.exists(was_targ_file):
                sensing_fp_si_wfe = poppy.FITSOpticalElement(opd=was_targ_file).opd
            else:
                sensing_fp_si_wfe = sensing_inst.get_wfe('si')

            sensing_fp_ote_wfe = sensing_inst.get_wfe('ote_field_dep')


            sihdu = fits.ImageHDU(sensing_fp_si_wfe)
            sihdu.header['EXTNAME'] = 'SENSING_SI_WFE'
            sihdu.header['CONTENTS'] = 'Model of SI WFE at sensing field point'
            sihdu.header['BUNIT'] = 'meter'
            sihdu.header['APERNAME'] = sensing_apername
            sihdu.header.add_history("This model for SI WFE was subtracted from the measured total WFE")
            sihdu.header.add_history("to estimate the OTE-only portion of the WFE.")
            opdhdu.append(sihdu)

            otehdu = fits.ImageHDU(sensing_fp_ote_wfe)
            otehdu.header['EXTNAME'] = 'SENSING_OTE_FD_WFE'
            otehdu.header['CONTENTS'] = 'Model of OTE field dependent WFE at sensing field point'
            otehdu.header['BUNIT'] = 'meter'
            otehdu.header['APERNAME'] = sensing_apername
            otehdu.header.add_history("This model for OTE field dependence was subtracted from the measured total WFE")
            otehdu.header.add_history("to estimate the OTE global portion of the WFE, at the master chief ray")
            opdhdu.append(otehdu)

            # Subtract the SI WFE from the WSS OPD, to obtain an estimated OTE-only OPD
            opdhdu[0].data -= (sensing_fp_si_wfe + sensing_fp_ote_wfe) * ote_pupil_mask
            opdhdu[0].header['CONTENTS'] = "Estimated OTE WFE from Wavefront Sensing Measurements"
            opdhdu[0].header.add_history(f"Estimating SI WFE at sensing field point {sensing_apername}.")
            opdhdu[0].header.add_history('  See FITS extension SENSING_SI_WFE for the SI WFE model used.')
            opdhdu[0].header.add_history('  Subtracted SI WFE to estimate OTE-only global WFE.')
            opdhdu[0].header.add_history(f"Estimating OTE field dependence term at {sensing_apername}.")
            opdhdu[0].header.add_history(f"  Selected instrument field point is at V2,V3 = {sensing_inst._tel_coords()}.")
            opdhdu[0].header.add_history('  See FITS extension SENSING_OTE_FD_WFE for the WFE model used.')
            opdhdu[0].header.add_history('  Subtracted OTE field dependence to estimate OTE global WFE.')

            if plot or save_ote_wfe:
                # Either of these options will need the total OTE WFE.
                # Under normal circumstances webbpsf will compute this later automatically, but if needed we do it here too
                selected_fp_ote_wfe = self.get_wfe('ote_field_dep')
                total_ote_wfe_at_fp = opdhdu[0].data+(selected_fp_ote_wfe*ote_pupil_mask)

            if plot:
                axes[0,1].imshow(sensing_fp_si_wfe * ote_pupil_mask, **plot_kwargs)
                axes[0,1].set_title(f"SI OPD\nat {sensing_apername}")
                axes[0,1].set_xlabel(f"RMS: {utils.rms(sensing_fp_si_wfe * 1e9, ote_pupil_mask):.2f} nm rms")

                axes[0,2].imshow(opdhdu[0].data + sensing_fp_ote_wfe * ote_pupil_mask , **plot_kwargs)
                axes[0,2].set_title(f"OTE total OPD at sensing field point\ninferred from {os.path.basename(filename)}")
                axes[0,2].set_xlabel(f"RMS: {utils.rms(opdhdu[0].data*1e9, ote_pupil_mask):.2f} nm rms")

                axes[1,0].imshow(sensing_fp_ote_wfe * ote_pupil_mask, **plot_kwargs)
                axes[1,0].set_title(f"OTE field dependent OPD\nat {sensing_apername}")
                axes[1,0].set_xlabel(f"RMS: {utils.rms(sensing_fp_ote_wfe * 1e9, ote_pupil_mask):.2f} nm rms")

                axes[1,1].imshow(selected_fp_ote_wfe * ote_pupil_mask, **plot_kwargs)
                axes[1,1].set_title(f"OTE field dependent OPD\nat current field point in {self.name} {self.detector}")
                axes[1,1].set_xlabel(f"RMS: {utils.rms(selected_fp_ote_wfe * 1e9, ote_pupil_mask):.2f} nm rms")

                axes[1,2].imshow(total_ote_wfe_at_fp, **plot_kwargs)
                axes[1,2].set_title(f"Total OTE OPD at current FP in {self.name} {self.detector}\ninferred from {os.path.basename(filename)}")
                axes[1,2].set_xlabel(f"RMS: {utils.rms(total_ote_wfe_at_fp*1e9, ote_pupil_mask):.2f} nm rms")

                plt.tight_layout()

            if save_ote_wfe:
                # If requested, export the OPD for use in other external calculations.
                # We save out the total OTE WFE inferred at the selected instrument field point.
                outname = filename.replace(".fits", f"-ote-wfe-for-{self.name}-{self.detector}.fits")
                from copy import deepcopy
                opdhdu_at_si_fp = deepcopy(opdhdu)

                v2v3 = self._tel_coords()
                opdhdu_at_si_fp[0].header.add_history(f"Estimating OTE field dependence term in {self.name} {self.detector}.")
                opdhdu_at_si_fp[0].header.add_history(f"  Selected instrument field point is at V2,V3 = {v2v3}.")
                opdhdu_at_si_fp[0].header.add_history(f"Saving out total estimated OTE WFE (global+field dep) at that field point.")
                opdhdu_at_si_fp[0].header['INSTRUME'] = self.name
                opdhdu_at_si_fp[0].header['DETECTOR'] = self.detector
                opdhdu_at_si_fp[0].header['APERNAME'] = self.aperturename
                opdhdu_at_si_fp[0].header['V2'] = self.aperturename

                # Save files with output units of microns, for consistency with other OPD files
                opdhdu_at_si_fp[0].data = total_ote_wfe_at_fp * 1e6
                opdhdu_at_si_fp[0].header['BUNIT'] = 'micron'

                opdhdu_at_si_fp.writeto(outname, overwrite=True)
                print(f"*****\nSaving estimated OTE-only WFE to file:\n\t{outname}\n*****")

        self.pupilopd = opdhdu

    def load_wss_opd_by_date(self, date=None, choice='closest', verbose=True, plot=False, **kwargs):
        """Load an OPD produced by the JWST WSS into this instrument instance, specified by filename.

        This does a MAST query by date to identify the relevant OPD file, then calls load_wss_opd.

        Parameters
        ----------
        date: string
            Date time in UTC as ISO-format string, a la 2021-12-25T07:20:00
            Note, if date is left unspecified as None, the most recent
            available measurement will be retrieved.
        choice : string
            Method to choose which OPD file to use, e.g. 'before', 'after'

        Further keyword parameters may be passed via **kwargs to load_wss_opd


        """

        if date is None:
            date = astropy.time.Time.now().isot
        opd_fn = webbpsf.mast_wss.get_opd_at_time(date, verbose=verbose, choice=choice, **kwargs)
        self.load_wss_opd(opd_fn, verbose=verbose, plot=plot, **kwargs)




class MIRI(JWInstrument):
    """ A class modeling the optics of MIRI, the Mid-InfraRed Instrument.

    Relevant attributes include `filter`, `image_mask`, and `pupil_mask`.

    The pupil will auto-select appropriate values for the coronagraphic filters
    if the auto_pupil attribute is set True (which is the default).

    Special Options:

    The 'coron_shift_x' and 'coron_shift_y' options offset a coronagraphic mask in order to
    produce PSFs centered in the output image, rather than offsetting the PSF. This is useful
    for direct PSF convolutions. Values are in arcsec.
    ```
    miri.options['coron_shift_x'] = 3  # Shifts mask 3" to right; or source 3" to left.
    ```

    """

    def __init__(self):
        self.auto_pupil = True
        JWInstrument.__init__(self, "MIRI")
        self.pixelscale = 0.1108  # MIRI average of X and Y pixel scales. Source: SIAF PRDOPSSOC-031, 2021 April
        self._rotation = 4.834  # V3IdlYAngle, Source: SIAF PRDOPSSOC-031
                                # This is rotation counterclockwise; when summed with V3PA it will yield the Y axis PA on sky

        self.options['pupil_shift_x'] = -0.0069 # CV3 on-orbit estimate (RPT028027) + OTIS delta from predicted (037134)
        self.options['pupil_shift_y'] = -0.0027

        self.image_mask_list = ['FQPM1065', 'FQPM1140', 'FQPM1550', 'LYOT2300', 'LRS slit']
        self.pupil_mask_list = ['MASKFQPM', 'MASKLYOT', 'P750L']

        self._image_mask_apertures = {'FQPM1065': 'MIRIM_CORON1065',
                                      'FQPM1140': 'MIRIM_CORON1140',
                                      'FQPM1550': 'MIRIM_CORON1550',
                                      'LYOT2300': 'MIRIM_CORONLYOT'}
        self.auto_aperturename = True

        self.monochromatic = 8.0
        self._IFU_pixelscale = {
            'Ch1': (0.18, 0.19),
            'Ch2': (0.28, 0.19),
            'Ch3': (0.39, 0.24),
            'Ch4': (0.64, 0.27),
        }
        # The above tuples give the pixel resolution (perpendicular to the slice, along the slice).
        # The pixels are not square.

        self._detectors = {'MIRIM': 'MIRIM_FULL'}  # Mapping from user-facing detector names to SIAF entries.
        self.detector = self.detector_list[0]
        self._detector_npixels = 1024
        self.detector_position = (512, 512)

        self._si_wfe_class = optics.MIRIFieldDependentAberrationAndObscuration

    def _get_default_fov(self):
        """ Return default FOV in arcseconds """
        return 12

    @JWInstrument.filter.setter
    def filter(self, value):
        super(MIRI, self.__class__).filter.__set__(self, value)

        if self.auto_pupil:
            # set the pupil shape based on filter
            if self.filter.endswith('C'):
                # coronagraph masks
                if self.filter[1] == '1':
                    self.pupil_mask = 'MASKFQPM'
                else:
                    self.pupil_mask = 'MASKLYOT'
            else:
                # no mask, i.e. full pupil
                self.pupil_mask = None

    def _validate_config(self, **kwargs):
        """Validate instrument config for MIRI
        """
        if self.filter.startswith("MRS-IFU"):
            raise NotImplementedError("The MIRI MRS is not yet implemented.")
        return super(MIRI, self)._validate_config(**kwargs)

    def _addAdditionalOptics(self, optsys, oversample=2):
        """Add coronagraphic or spectrographic optics for MIRI.
        Semi-analytic coronagraphy algorithm used for the Lyot only.

        """

        # For MIRI coronagraphy, all the coronagraphic optics are rotated the same
        # angle as the instrument is, relative to the primary. So they see the unrotated
        # telescope pupil. Likewise the LRS grism is rotated but its pupil stop is not.
        #
        # We model this by just not rotating till after the coronagraph. Thus we need to
        # un-rotate the primary that was already created in get_optical_system.
        # This approach is required computationally so we can work in an unrotated frame
        # aligned with the FQPM axes.

        defaultpupil = optsys.planes.pop(2)  # throw away the rotation of the entrance pupil we just added

        if self.include_si_wfe:
            # temporarily remove the SI internal aberrations
            # from the system - will add back in after the
            # coronagraph planes.
            miri_aberrations = optsys.planes.pop(2)

        # Add image plane mask
        # For the MIRI FQPMs, we require the star to be centered not on the middle pixel, but
        # on the cross-hairs between four pixels. (Since that is where the FQPM itself is centered)
        # This is with respect to the intermediate calculation pixel scale, of course, not the
        # final detector pixel scale.
        if ((self.image_mask is not None and 'FQPM' in self.image_mask)
                or 'force_fqpm_shift' in self.options):
            optsys.add_pupil(poppy.FQPM_FFT_aligner())

        # Allow arbitrary offsets of the focal plane masks with respect to the pixel grid origin;
        # In most use cases it's better to offset the star away from the mask instead, using
        # options['source_offset_*'], but doing it this way instead is helpful when generating
        # the Pandeia ETC reference PSF library.
        offsets = {'shift_x': self.options.get('coron_shift_x', None),
                   'shift_y': self.options.get('coron_shift_y', None)}

        def make_fqpm_wrapper(name, wavelength):
            container = poppy.CompoundAnalyticOptic(name=name,
                                                    opticslist=[poppy.IdealFQPM(wavelength=wavelength,
                                                                                name=self.image_mask,
                                                                                **offsets),
                                                                poppy.SquareFieldStop(size=24,
                                                                                      rotation=self._rotation,
                                                                                      **offsets)])
            return container

        if self.image_mask == 'FQPM1065':
            optsys.add_image(make_fqpm_wrapper("MIRI FQPM 1065", 10.65e-6))
            trySAM = False
        elif self.image_mask == 'FQPM1140':
            optsys.add_image(make_fqpm_wrapper("MIRI FQPM 1140", 11.40e-6))
            trySAM = False
        elif self.image_mask == 'FQPM1550':
            optsys.add_image(make_fqpm_wrapper("MIRI FQPM 1550", 15.50e-6))
            trySAM = False
        elif self.image_mask == 'LYOT2300':
            # diameter is 4.25 (measured) 4.32 (spec) supposedly 6 lambda/D
            # optsys.add_image(function='CircularOcculter',radius =4.25/2, name=self.image_mask)
            # Add bar occulter: width = 0.722 arcsec (or perhaps 0.74, Dean says there is ambiguity)
            # optsys.add_image(function='BarOcculter', width=0.722, angle=(360-4.76))
            # position angle of strut mask is 355.5 degrees  (no = =360 -2.76 degrees
            # optsys.add_image(function='fieldstop',size=30)
            container = poppy.CompoundAnalyticOptic(name="MIRI Lyot Occulter",
                                            opticslist=[poppy.CircularOcculter(radius=4.25 / 2, name=self.image_mask, **offsets),
                                                        poppy.BarOcculter(width=0.722, height=31, **offsets),
                                                        poppy.SquareFieldStop(size=30, rotation=self._rotation, **offsets)])
            optsys.add_image(container)
            trySAM = False  # FIXME was True - see https://github.com/mperrin/poppy/issues/169
            SAM_box_size = [5, 20]
        elif self.image_mask == 'LRS slit':
            # one slit, 5.5 x 0.6 arcsec in height (nominal)
            #           4.7 x 0.51 arcsec (measured for flight model. See MIRI-TR-00001-CEA)
            #
            # Per Klaus Pontoppidan: The LRS slit is aligned with the detector x-axis, so that the
            # dispersion direction is along the y-axis.
            optsys.add_image(optic=poppy.RectangularFieldStop(width=4.7, height=0.51,
                                                              rotation=self._rotation, name=self.image_mask))
            trySAM = False
        else:
            optsys.add_image()
            trySAM = False

        if ((self.image_mask is not None and 'FQPM' in self.image_mask)
                or 'force_fqpm_shift' in self.options):
            optsys.add_pupil(poppy.FQPM_FFT_aligner(direction='backward'))

        # add pupil plane mask
        shift_x, shift_y = self._get_pupil_shift()
        rotation = self.options.get('pupil_rotation', None)

        if self.pupil_mask == 'MASKFQPM':
            optsys.add_pupil(transmission=self._datapath + "/optics/MIRI_FQPMLyotStop.fits.gz",
                             name=self.pupil_mask,
                             flip_y=True, shift_x=shift_x, shift_y=shift_y, rotation=rotation)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'MASKLYOT':
            optsys.add_pupil(transmission=self._datapath + "/optics/MIRI_LyotLyotStop.fits.gz",
                             name=self.pupil_mask,
                             flip_y=True, shift_x=shift_x, shift_y=shift_y, rotation=rotation)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'P750L':
            optsys.add_pupil(transmission=self._datapath + "/optics/MIRI_LRS_Pupil_Stop.fits.gz",
                             name=self.pupil_mask,
                             flip_y=True, shift_x=shift_x, shift_y=shift_y, rotation=rotation)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        else:  # all the MIRI filters have a tricontagon outline, even the non-coron ones.
            optsys.add_pupil(transmission=self._WebbPSF_basepath + "/tricontagon.fits.gz",
                             name='filter cold stop', shift_x=shift_x, shift_y=shift_y, rotation=rotation)
            # FIXME this is probably slightly oversized? Needs to have updated specifications here.

        if self.include_si_wfe:
            # now put back in the aberrations we grabbed above.
            optsys.add_pupil(miri_aberrations)

        optsys.add_rotation(-self._rotation, hide=True)
        optsys.planes[-1].wavefront_display_hint = 'intensity'

        return (optsys, trySAM, SAM_box_size if trySAM else None)


    def _update_aperturename(self):
        """Determine sensible SIAF aperture names for MIRI. Implements the auto_aperturename functionality.
        Called after detector is changed
        """

        str_debug = '_update_aperturename BEFORE - Det: {}, Ap: {}, ImMask: {}, PupMask: {}, DetPos: {}'.format(
            self._detector, self._aperturename, self.image_mask, self.pupil_mask, self.detector_position
        )
        _log.debug(str_debug)

        # Need to send correct aperture name for coronagraphic masks
        if (self._image_mask is not None):
            if 'LRS' in self._image_mask:
                apname = 'MIRIM_FULL'  # LRS slit uses full array readout
            else:
                apname = self._image_mask_apertures[self._image_mask]
        else:
            apname = 'MIRIM_FULL'

        # Call aperturename.setter to update ap ref coords and DetectorGeometry class
        self.aperturename = apname

        str_debug = '_update_aperturename AFTER  - Det: {}, Ap: {}, ImMask: {}, PupMask: {}, DetPos: {}'.format(
            self._detector, self._aperturename, self.image_mask, self.pupil_mask, self.detector_position
        )
        _log.debug(str_debug)


    def _get_fits_header(self, hdulist, options):
        """ Format MIRI-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        super(MIRI, self)._get_fits_header(hdulist, options)

        hdulist[0].header['GRATNG14'] = ('None', 'MRS Grating for channels 1 and 4')
        hdulist[0].header['GRATNG23'] = ('None', 'MRS Grating for channels 2 and 3')
        hdulist[0].header['FLATTYPE'] = ('?', 'Type of flat field to be used: all, one, principal')
        hdulist[0].header['CCCSTATE'] = ('open', 'Contamination Control Cover state: open, closed, locked')
        if self.image_mask is not None:
            hdulist[0].header['TACQNAME'] = ('None', 'Target acquisition file name')


class NIRCam(JWInstrument):
    """ A class modeling the optics of NIRCam.

    Relevant attributes include `filter`, `image_mask`, and `pupil_mask`.

    The NIRCam class is smart enough to automatically select the appropriate
    pixel scale for the short or long wavelength channel
    based on the selected detector (NRCA1 vs NRCA5, etc), and also on
    whether you request a short or long wavelength filter. The auto-selection
    based on filter name can be disabled, if necessary, by setting `.auto_channel = False`.
    Setting the detector name always toggles the channel regardless of `auto_channel`.

    Note, if you use the `monochromatic` option for calculating PSFs, that does not
    invoke the automatic channel selection. Make sure to set the correct channel *prior*
    to calculating any monochromatic PSFs.

    Similarly, SIAF aperture names are automatically chosen based on detector, filter,
    image mask, and pupil mask settings. The auto-selection can be disabled by
    setting `.auto_aperturename = False`. SIAF aperture information is mainly used for
    coordinate transformations between detector science pixels and telescope V2/V3.

    Special Options:
    The 'bar_offset' option allows specification of an offset position
    along one of the coronagraph bar occulters, in arcseconds.
    ```
    nc.image_mask = 'MASKLWB'
    nc.options['bar_offset'] = 3 # 3 arcseconds towards the right (narrow end on module A)
    ```

    Similarly, the 'coron_shift_x' and 'coron_shift_y' options will offset the mask in order
    to produce PSFs centered in the output image, rather than offsetting the PSF. This is useful
    for direct PSF convolutions of an image. Values are in arcsec. These options move the mask
    in the opposite sense as nc.options['bar_offset'].
    ```
    nc.options['coron_shift_x'] = 3  # Shifts mask 3" to right, equivalent to source 3" to left.
    ```

    The 'nd_squares' option allows toggling on and off the ND squares for TA in the simulation.
    Note that these of course aren't removable in the real instrument; this option exists solely for
    some simulation purposes.


    """
    SHORT_WAVELENGTH_MIN = 0.6 * 1e-6
    SHORT_WAVELENGTH_MAX = LONG_WAVELENGTH_MIN = 2.35 * 1e-6
    LONG_WAVELENGTH_MAX = 5.3 * 1e-6

    def __init__(self):
        self._pixelscale_short = 0.0311  # average over both X and Y for short-wavelen channels, SIAF PRDOPSSOC-031, 2021 April
        self._pixelscale_long = 0.0630  # average over both X and Y for long-wavelen channels, SIAF PRDOPSSOC-031, 2021 April
        self.pixelscale = self._pixelscale_short

        self.options['pupil_shift_x'] = 0  # Set to 0 since NIRCam FAM corrects for PM shear in flight
        self.options['pupil_shift_y'] = 0

        # need to set up a bunch of stuff here before calling superclass __init__
        # so the overridden filter setter will work successfully inside that.
        self.auto_channel = True
        self.auto_aperturename = True
        self._filter = 'F200W'
        self._detector = 'NRCA1'

        JWInstrument.__init__(self, "NIRCam")  # do this after setting the long & short scales.
        self._detector = 'NRCA1' # Must re-do this after superclass init since that sets it to None.
                                 # This is an annoying workaround to ensure all the auto-channel stuff is ok

        self.pixelscale = self._pixelscale_short  # need to redo 'cause the __init__ call will reset it to zero
        self._filter = 'F200W'  # likewise need to redo

        self.image_mask_list = ['MASKLWB', 'MASKSWB', 'MASK210R', 'MASK335R', 'MASK430R']
        self._image_mask_apertures = {'MASKLWB': 'NRCA5_MASKLWB',
                                      'MASKSWB': 'NRCA4_MASKSWB',
                                      'MASK210R': 'NRCA2_MASK210R',
                                      'MASK335R': 'NRCA5_MASK335R',
                                      'MASK430R': 'NRCA5_MASK430R'}

        self.pupil_mask_list = ['CIRCLYOT', 'WEDGELYOT', 'MASKRND', 'MASKSWB', 'MASKLWB',
                                # The last 3 of the above are synonyms for the first 2
                                'WEAK LENS +4', 'WEAK LENS +8', 'WEAK LENS -8', 'WEAK LENS +12 (=4+8)', 'WEAK LENS -4 (=4-8)',
                                'WLP4', 'WLM4', 'WLP8', 'WLM8', 'WLP12']

        self._detectors = dict()
        det_list = ['A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5']
        for name in det_list: self._detectors["NRC{0}".format(name)] = 'NRC{0}_FULL'.format(name)
        self.detector = self.detector_list[0]
        self._aperturename = '{}_FULL'.format(self._detector)  # SIAF aperture name

        self._si_wfe_class = optics.NIRCamFieldAndWavelengthDependentAberration

    def _update_aperturename(self):
        """Determine sensible SIAF aperture names for NIRCam. Implements the auto_aperturename functionality:
        when the detector is changed, the aperture updates to <det>_FULL, and coronagraph masks auto select the
        appropriate aperture. Other apertures can be selected using set_position_from_aperture_name

        Called after detector is changed; see detector.setter

        """

        str_debug = '_update_aperturename BEFORE - Det: {}, Ap: {}, ImMask: {}, PupMask: {}, DetPos: {}'.format(
            self._detector, self._aperturename, self.image_mask, self.pupil_mask, self.detector_position
        )
        _log.debug(str_debug)

        # Need to send correct aperture name for coronagraphic masks due to detector shift
        if (self._image_mask is not None):
            aps_modA = {'MASKLWB': 'NRCA5_FULL_MASKLWB',
                        'MASKSWB': 'NRCA4_FULL_MASKSWB',
                        'MASK210R': 'NRCA2_FULL_MASK210R',
                        'MASK335R': 'NRCA5_FULL_MASK335R',
                        'MASK430R': 'NRCA5_FULL_MASK430R'}
            # Choose coronagraphic subarray apertures for Module B
            aps_modB = {'MASKLWB': 'NRCB5_MASKLWB',
                        'MASKSWB': 'NRCB3_MASKSWB',
                        'MASK210R': 'NRCB1_MASK210R',
                        'MASK335R': 'NRCB5_MASK335R',
                        'MASK430R': 'NRCB5_MASK430R'}
            apname = aps_modA[self._image_mask] if self.module=='A' else aps_modB[self._image_mask]
            _log.debug(f"Inferred {apname} from coronagraph focal plane mask selected.")
        elif (self._pupil_mask is not None) and (('LYOT' in self._pupil_mask) or ('MASK' in self._pupil_mask)):
            # Want to use full frame apertures if only Lyot stops defined (no image mask)
            # Unfortunately, no full frame SIAF apertures are defined for Module B w/ Lyot
            # so we must select the subarray apertures as a special case.
            if 'long' in self.channel:
                if ('WEDGE' in self._pupil_mask) or ('LWB' in self._pupil_mask):
                    apname = 'NRCA5_FULL_WEDGE_BAR' if self.module=='A' else 'NRCB5_MASKLWB'
                else:
                    apname = 'NRCA5_FULL_WEDGE_RND' if self.module=='A' else 'NRCB5_MASK335R'
            else:
                if ('WEDGE' in self._pupil_mask) or ('SWB' in self._pupil_mask):
                    apname = 'NRCA4_FULL_WEDGE_BAR' if self.module=='A' else 'NRCB3_MASKSWB'
                else:
                    apname = 'NRCA2_FULL_WEDGE_RND' if self.module=='A' else 'NRCB1_MASK210R'
                    _log.debug(f"Inferred {apname} from coronagraph Lyot mask selected, and channel={self.channel}, module={self.module}")
        else:
            apname = self._detectors[self._detector]
            _log.debug(f"Inferred {apname} from selected detector.")

        # Call aperturename.setter to update ap ref coords and DetectorGeometry class
        self.aperturename = apname

        str_debug = '_update_aperturename AFTER  - Det: {}, Ap: {}, ImMask: {}, PupMask: {}, DetPos: {}'.format(
            self._detector, self._aperturename, self.image_mask, self.pupil_mask, self.detector_position
        )
        _log.debug(str_debug)

    @JWInstrument.aperturename.setter
    def aperturename(self, value):
        # Explicitly update detector reference coordinates,
        # otherwise old coordinates can persist under certain circumstances

        # Get NIRCam SIAF apertures
        try:
            ap = self.siaf[value]
        except KeyError:
            _log.warning(f'Aperture name {value} not a valid NIRCam pysiaf name')
            # Alternatives in case we are running an old pysiaf PRD
            if value=='NRCA5_FULL_WEDGE_BAR':
                newval = 'NRCA5_FULL_MASKLWB'
            elif value=='NRCA5_FULL_WEDGE_RND':
                newval = 'NRCA5_FULL_MASK335R'
            elif value=='NRCA4_FULL_WEDGE_BAR':
                newval = 'NRCA4_FULL_MASKSWB'
            elif value=='NRCA2_FULL_WEDGE_RND':
                newval = 'NRCA2_FULL_MASK210R'
            else:
                newval = None

            if newval is not None:
                # Set alternative aperture name as bandaid to continue
                value = newval
                _log.warning('Possibly running an old version of pysiaf missing some NIRCam apertures. Continuing with old aperture names.')
            else:
                return

        # Only update if new value is different
        if self._aperturename != value:
            self._aperturename = value
            # Update detector reference coordinates
            self.detector_position = (ap.XSciRef, ap.YSciRef)

            # Check if detector is correct
            new_det =  self._aperturename[0:5]
            if new_det != self._detector:
                new_channel = 'long' if new_det[-1] == '5' else 'short'
                self._switch_channel(new_channel)
                self._detector = new_det

            # Update DetectorGeometry class
            self._detector_geom_info = DetectorGeometry(self.siaf, self._aperturename)
            _log.info("NIRCam aperture name updated to {}".format(self._aperturename))

    @property
    def module(self):
        return self._detector[3]
        # note, you can't set module directly; it's inferred based on detector.

    @module.setter
    def module(self, value):
        raise RuntimeError("NIRCam module is not directly settable; set detector instead.")

    @property
    def channel(self):
        return 'long' if self.detector.endswith('5') else 'short'
        # note, you can't set channel directly; it's inferred based on detector.

    @channel.setter
    def channel(self, value):
        raise RuntimeError("NIRCam channel is not directly settable; set filter or detector instead.")

    @JWInstrument.detector.setter # override setter in this subclass, to implement auto channel switch
    def detector(self, value):
        """ Set detector, including reloading the relevant info from SIAF """
        if value.upper() not in self.detector_list:
            raise ValueError("Invalid detector. Valid detector names are: {}".format(', '.join(self.detector_list)))
        # set the channel based on the requested detector
        new_channel = 'long' if value[-1] == '5' else 'short'
        self._switch_channel(new_channel)
        self._detector = value.upper()
        self._update_aperturename()

    def _switch_channel(self,channel):
        """ Toggle to either SW or LW channel.
        This changes the detector name and the pixel scale,
        unless the user has set a custom/nonstandard pixel scale manually.
        """
        if self.channel == channel:
            return # nothing to do
        _log.debug("Automatically changing NIRCam channel SW/LW to "+channel)
        if channel=='long':
            # ensure long wave by switching to detector 5
            self._detector = self._detector[0:4] + '5'
            if self.pixelscale == self._pixelscale_short:
                self.pixelscale = self._pixelscale_long
                _log.info("NIRCam pixel scale switched to %f arcsec/pixel for the "
                          "long wave channel." % self.pixelscale)
        elif channel=='short':
            # only change detector if the detector was already LW;
            # don't override selection of a particular SW SCA otherwise
            if self._detector[-1] == '5':
                self._detector = self._detector[0:4] + '1'
            if self.pixelscale == self._pixelscale_long:
                self.pixelscale = self._pixelscale_short
                _log.info("NIRCam pixel scale switched to %f arcsec/pixel for the "
                          "short wave channel." % self.pixelscale)
        else:
            raise ValueError("Invalid NIRCam channel name: {}".format(channel))

    @JWInstrument.filter.setter
    def filter(self, value):
        super(NIRCam, self.__class__).filter.__set__(self, value)

        if self.auto_channel or self.auto_aperturename:
            # set the channel (via setting the detector) based on filter
            if self.filter=='WLP4':
                # special case, weak lens 4 is actually a filter too but isn't named like one
                wlnum = 212
            else:
                wlnum = int(self.filter[1:4])
            new_channel = 'long' if wlnum >= 250 else 'short'
            cur_channel = self.channel

            if self.auto_channel:
                self._switch_channel(new_channel)

            # Only change ap name if filter choice forces us to a different channel
            if self.auto_aperturename and (cur_channel != new_channel):
                self._update_aperturename()

    # Need to redefine image_mask.setter because _image_mask_apertures has limited aperture definitions
    @JWInstrument.image_mask.setter
    def image_mask(self, name):
        if name == "": name = None
        if name is not None:
            if name in self.image_mask_list:
                pass  # there's a perfect match, this is fine.
            else:
                name = name.upper()  # force to uppercase
                if name not in self.image_mask_list:  # if still not found, that's an error.
                    raise ValueError("Instrument %s doesn't have an image mask called '%s'." % (self.name, name))
        self._image_mask = name

        # Update aperture position, which updates detector and detector position
        self._update_aperturename()
        self.set_position_from_aperture_name(self._aperturename)

    @JWInstrument.pupil_mask.setter
    def pupil_mask(self, name):

        if name != self._pupil_mask:
            # only apply updates if the value is in fact new

            super(NIRCam, self.__class__).pupil_mask.__set__(self, name)
            _log.info(f"NIRCam pupil mask setter: aperturename {self._aperturename}")

            # infer a new aperture, since the coronagraph mask choice affects this
            self._update_aperturename()

            # Update aperture position, which updates detector and detector position
            self.set_position_from_aperture_name(self._aperturename)

    def _validate_config(self, **kwargs):
        """Validate instrument config for NIRCam

        For NIRCam, this automatically handles toggling between the short-wave and long-wave channels.
        I.e it selects a pixelscale based on the wavelengths requested
        """
        wavelengths = np.array(kwargs['wavelengths'])
        if np.min(wavelengths) < self.SHORT_WAVELENGTH_MIN:
            raise RuntimeError("The requested wavelengths are too short to be imaged with NIRCam")
        if np.max(wavelengths) > self.LONG_WAVELENGTH_MAX:
            raise RuntimeError("The requested wavelengths are too long to be imaged with NIRCam")
        if self.channel == 'short' and np.max(wavelengths) > self.SHORT_WAVELENGTH_MAX:
            raise RuntimeError("The requested wavelengths are too long for NIRCam short wave channel.")
        if self.channel == 'long' and np.min(wavelengths) < self.LONG_WAVELENGTH_MIN:
            raise RuntimeError("The requested wavelengths are too short for NIRCam long wave channel.")

        return super(NIRCam, self)._validate_config(**kwargs)

    def _addAdditionalOptics(self, optsys, oversample=2):
        """Add coronagraphic optics for NIRCam

        See Krist et al. 2007, 2009 SPIE

        Three circular occulters: HWHM = 6 lambda/D at 2.1, 3.35, 4.3
                                       = 0.4, 0.64, 0.8 arcsec (avg)
                                       assuming D_tel=6.5m exactly:
                                        = 0.3998, 0.6378, 0.8187 arcsec

        Two linear bar occulters: Wedges vary from HWHM = 2 lam/D to 6 lam/D at 2.1 and 4.6 micron
                    2.1e-6:    HWHM = 0.13327 to 0.3998
                    4.6e-6:    HWHM = 0.27290 to 0.8187
            The matching Lyot stop for the wedges are tuned for 4 lam/D.
            The linear ones have a fixed width at either side: maybe ~ 3-4 arcsec. Then a linear taper
            in between.


        Values of Sigma:
            For circular occulters, 0.3998 requires sigma = 5.253
                                    0.8187 requires sigma = 2.5652
                                    sigma = 2.10013932 / loc
                                    vs. Krist's statement sigma = 2.1001/hwhm

            For linear occulters, 0.3998 requires sigma = 4.5012
                                  0.13327 requires sigma = 13.5078

                        # This is NOT a linear relationship! It's a tricky inverse sin nonlinear thing.

        Empirical checks against John Krist's provided 430R and LWB files:
            430R should have sigma = 2.588496


        Since the Weak Lenses go in the pupil too, this function provides a convenient place to implement those as well.

        """

        # optsys.add_image(name='null for debugging NIRcam _addCoron') # for debugging
        from .optics import NIRCam_BandLimitedCoron

        nd_squares = self.options.get('nd_squares', True)

        SAM_box_size = None  # default

        # Allow arbitrary offsets of the focal plane masks with respect to the pixel grid origin;
        # In most use cases it's better to offset the star away from the mask instead, using
        # options['source_offset_*'], but doing it this way instead is helpful when generating
        # the Pandeia ETC reference PSF library.
        shifts = {'shift_x': self.options.get('coron_shift_x', None),
                  'shift_y': self.options.get('coron_shift_y', None)}

        if ((self.image_mask == 'MASK210R') or (self.image_mask == 'MASK335R') or
                (self.image_mask == 'MASK430R')):
            optsys.add_image(NIRCam_BandLimitedCoron(name=self.image_mask, module=self.module,
                                                     nd_squares=nd_squares, **shifts),
                             index=2)
            trySAM = False  # FIXME was True - see https://github.com/mperrin/poppy/issues/169
            SAM_box_size = 5.0
        elif ((self.image_mask == 'MASKSWB') or (self.image_mask == 'MASKLWB')):
            bar_offset = self.options.get('bar_offset', None)
            # If the bar offset is not provided, use the filter name to lookup the default
            # position. If an offset is provided and is a floating point value, use that
            # directly as the offset. Otherwise assume it's a filter name and try passing
            # that in to the auto offset. (that allows for selecting the narrow position, or
            # for simulating using a given filter at some other filter's position.)
            if bar_offset is None:
                auto_offset = self.filter
            else:
                try:
                    _ = float(bar_offset)
                    auto_offset = None
                except ValueError:
                    # If the "bar_offset" isn't a float, pass it to auto_offset instead
                    auto_offset = bar_offset
                    bar_offset = None

            optsys.add_image(NIRCam_BandLimitedCoron(name=self.image_mask, module=self.module,
                                                     nd_squares=nd_squares, bar_offset=bar_offset,
                                                     auto_offset=auto_offset, **shifts),
                             index=2)
            trySAM = False  # True FIXME
            SAM_box_size = [5, 20]
        elif ((self.pupil_mask is not None) and ('LENS' not in self.pupil_mask.upper())
                and ('WL' not in self.pupil_mask.upper() )):
            # no occulter selected but coronagraphic mode anyway. E.g. off-axis PSF
            # but don't add this image plane for weak lens calculations
            optsys.add_image(poppy.ScalarTransmission(name='No Image Mask Selected!'), index=2)
            trySAM = False
        else:
            trySAM = False

        # add pupil plane mask
        shift_x, shift_y = self._get_pupil_shift()
        rotation = self.options.get('pupil_rotation', None)

        # NIRCam as-built weak lenses, from WSS config file, PRDOPSFLT-027
        WLP4_diversity = 8.3443  # microns
        WLP8_diversity =  16.5932 # microns
        WLM8_diversity = -16.5593  # microns
        WL_wavelength = 2.12  # microns

        if self.pupil_mask == 'CIRCLYOT' or self.pupil_mask == 'MASKRND':
            optsys.add_pupil(transmission=self._datapath + "/optics/NIRCam_Lyot_Somb.fits.gz", name=self.pupil_mask,
                             flip_y=True, shift_x=shift_x, shift_y=shift_y, rotation=rotation, index=3)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'WEDGELYOT' or self.pupil_mask == 'MASKSWB' or self.pupil_mask == 'MASKLWB':
            optsys.add_pupil(transmission=self._datapath + "/optics/NIRCam_Lyot_Sinc.fits.gz", name=self.pupil_mask,
                             flip_y=True, shift_x=shift_x, shift_y=shift_y, rotation=rotation, index=3)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        # Note, for historical reasons there are multiple synonymous ways to specify the weak lenses
        # This includes versions that elide over the fact that WLP4 is in the filter wheel, plus
        # versions that take that into account explicitly.
        elif self.pupil_mask == 'WEAK LENS +4' or self.pupil_mask == 'WLP4' or (
                self.filter == 'WLP4' and self.pupil_mask is None) :
            optsys.add_pupil(optics.NIRCamFieldDependentWeakLens(name='WLP4', instrument=self,
                                                                 shift_x=shift_x, shift_y=shift_y, rotation=rotation,
                                                                ), index=3)
        elif self.pupil_mask == 'WEAK LENS +8' or (self.pupil_mask == 'WLP8' and self.filter != 'WLP4'):
            optsys.add_pupil(optics.NIRCamFieldDependentWeakLens(name='WLP8', instrument=self,
                                                                 shift_x=shift_x, shift_y=shift_y, rotation=rotation,
                                                                ), index=3)
        elif self.pupil_mask == 'WEAK LENS -8' or (self.pupil_mask == 'WLM8' and self.filter != 'WLP4'):
            optsys.add_pupil(optics.NIRCamFieldDependentWeakLens(name='WLM8', instrument=self,
                                                                 shift_x=shift_x, shift_y=shift_y, rotation=rotation,
                                                                ), index=3)
        elif self.pupil_mask == 'WEAK LENS +12 (=4+8)' or self.pupil_mask == 'WLP12' or (
                self.pupil_mask == 'WLP8' and self.filter == 'WLP4'):
            optsys.add_pupil(optics.NIRCamFieldDependentWeakLens(name='WLP12', instrument=self,
                                                                 shift_x=shift_x, shift_y=shift_y, rotation=rotation,
                                                                 ), index=3)
        elif self.pupil_mask == 'WEAK LENS -4 (=4-8)' or self.pupil_mask == 'WLM4' or (
                self.pupil_mask == 'WLM8' and self.filter == 'WLP4'):
            optsys.add_pupil(optics.NIRCamFieldDependentWeakLens(name='WLM4', instrument=self,
                                                                shift_x=shift_x, shift_y=shift_y, rotation=rotation,
                                                                 ), index=3)

        elif (self.pupil_mask is None and self.image_mask is not None):
            optsys.add_pupil(poppy.ScalarTransmission(name='No Lyot Mask Selected!'), index=3)
        else:
            optsys.add_pupil(transmission=self._WebbPSF_basepath + "/tricontagon_oversized_4pct.fits.gz",
                             name='filter stop', shift_x=shift_x, shift_y=shift_y, rotation=rotation)

        return (optsys, trySAM, SAM_box_size)

    def _get_fits_header(self, hdulist, options):
        """ Format NIRCam-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        super(NIRCam, self)._get_fits_header(hdulist, options)

        hdulist[0].header['MODULE'] = (self.module, 'NIRCam module: A or B')
        hdulist[0].header['CHANNEL'] = ('Short' if self.channel == 'short' else 'Long', 'NIRCam channel: long or short')
        # filter, pupil added by calc_psf header code
        hdulist[0].header['PILIN'] = ('False', 'Pupil imaging lens in optical path: T/F')

class NIRSpec(JWInstrument):
    """ A class modeling the optics of NIRSpec, in **imaging** mode.

    This is not a substitute for a spectrograph model, but rather a way of simulating a PSF as it
    would appear with NIRSpec in imaging mode (e.g. for target acquisition).  NIRSpec support is
    relatively simplistic compared to the other instruments at this point.

    Relevant attributes include `filter`. In addition to the actual filters, you may select 'IFU' to
    indicate use of the NIRSpec IFU, in which case use the `monochromatic` attribute to set the simulated wavelength.

    If a grating is selected in the pupil, then a rectangular pupil mask 8.41x7.91 m as projected onto the primary
    is added to the optical system. This is an estimate of the pupil stop imposed by the outer edge of the grating
    clear aperture, estimated based on optical modeling by Erin Elliot and Marshall Perrin.

    **Note: IFU to be implemented later**
    """

    def __init__(self):
        JWInstrument.__init__(self, "NIRSpec")
        self.pixelscale = 0.1043  # Average over both detectors.  SIAF PRDOPSSOC-031, 2021 April
        # Microshutters are 0.2x0.46 but we ignore that here.
        self._rotation = 138.4  # Average for both detectors in SIAF PRDOPSSOC-031
                                # This is rotation counterclockwise; when summed with V3PA it will yield the Y axis PA on sky
        self.filter_list.append("IFU")
        self._IFU_pixelscale = 0.1043  # same.
        self.monochromatic = 3.0
        self.filter = 'F110W'  # or is this called F115W to match NIRCam??

        self.options['pupil_shift_x'] = 0.0115  # CV3 on-orbit estimate (RPT028027) + OTIS delta from predicted (037134)
        self.options['pupil_shift_y'] = -0.0157

        # fixed slits
        self.image_mask_list = ['S200A1', 'S200A2', 'S400A1', 'S1600A1', 'S200B1',
                                'MSA all open', 'Single MSA open shutter',
                                'Three adjacent MSA open shutters', 'IFU']
        self.pupil_mask_list = ['NIRSpec grating']
        self.image_mask = 'MSA all open'
        self.pupil_mask = self.pupil_mask_list[-1]

        det_list = ['NRS1', 'NRS2']
        self._detectors = dict()
        for name in det_list: self._detectors[name] = '{0}_FULL'.format(name)
        self.detector = self.detector_list[0]
        self.detector_position = (1380, 1024)   # near S1600A1 square aperture / ISIM1 field point. see #348.
        self._si_wfe_class = optics.NIRSpecFieldDependentAberration  # note we end up adding 2 instances of this.

    def _validate_config(self, **kwargs):
        if self.filter.startswith("IFU"):
            raise NotImplementedError("The NIRSpec IFU is not yet implemented.")
        return super(NIRSpec, self)._validate_config(**kwargs)

    def _addAdditionalOptics(self, optsys, oversample=2):
        """ Add fixed slit optics for NIRSpec

        See Table 3-6 of NIRSpec Ops Concept Document, ESA-JWST-TN-0297 / JWST-OPS-003212

        """
        from .optics import NIRSpec_three_MSA_shutters, NIRSpec_MSA_open_grid
        trySAM = False  # semi-analytic method never applicable here.
        SAM_box_size = None

        if self.image_mask == 'S200A1' or self.image_mask == 'S200A2' or self.image_mask == 'S200B1':
            # three identical slits, 0.2 x 3.2 arcsec in length
            optsys.add_image(optic=poppy.RectangularFieldStop(width=0.2,
                                                              height=3.2, name=self.image_mask + " slit"))
        elif self.image_mask == 'S400A1':
            # one slit, 0.4 x 3.65 arcsec in height
            optsys.add_image(optic=poppy.RectangularFieldStop(width=0.4,
                                                              height=3.65, name=self.image_mask + " slit"))
        elif self.image_mask == 'S1600A1':
            # square aperture for exoplanet spectroscopy
            optsys.add_image(optic=poppy.RectangularFieldStop(width=1.6,
                                                              height=1.6, name=self.image_mask + " square aperture"))
        elif self.image_mask == 'IFU':
            # square aperture for the entrance to the slicer.
            # DOES NOT ACTUALLY MODEL THE SLICER OPTICS AT ALL!
            # Values talen from pre-flight SIAF, fall 2017
            optsys.add_image(optic=poppy.RectangularFieldStop(width=3.193,
                                                              height=3.097, name="IFU entrance"))
        elif self.image_mask == 'MSA all open':
            # all MSA shutters open
            optsys.add_image(optic=NIRSpec_MSA_open_grid(name=self.image_mask))
        elif self.image_mask == 'Single MSA open shutter':
            # one MSA open shutter aperture
            optsys.add_image(optic=poppy.RectangularFieldStop(width=0.2,
                                                              height=0.45, name=self.image_mask))
        elif self.image_mask == 'Three adjacent MSA open shutters':
            optsys.add_image(optic=NIRSpec_three_MSA_shutters(name=self.image_mask))

        if ((self.pupil_mask is not None) and ('grating' in self.pupil_mask.lower())):
            # NIRSpec pupil stop at the grating appears to be a rectangle.
            # see notes and ray trace from Erin Elliot in the webbpsf-data/NIRSpec/sources directory
            optsys.add_pupil(optic=poppy.RectangleAperture(height=8.41,
                                                           width=7.91, name='Pupil stop at grating wheel'))
            optsys.planes[-1].wavefront_display_hint = 'intensity'

        # Add here a second instance of the instrument WFE, representing the WFE in the
        # collimator and camera.
        if self.include_si_wfe:
            optsys.add_pupil(optic=self._si_wfe_class(self, where='spectrograph'))

        return (optsys, trySAM, SAM_box_size)

    def _get_fits_header(self, hdulist, options):
        """ Format NIRSpec-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        super(NIRSpec, self)._get_fits_header(hdulist, options)
        hdulist[0].header['GRATING'] = ('None', 'NIRSpec grating element name')
        hdulist[0].header['APERTURE'] = (str(self.image_mask), 'NIRSpec slit aperture name')


class NIRISS(JWInstrument):
    """ A class modeling the optics of the Near-IR Imager and Slit Spectrograph
        (formerly TFI)

    Relevant attributes include `image_mask`, and `pupil_mask`.

    **Imaging:**

    WebbPSF models the direct imaging and nonredundant aperture masking modes of NIRISS in the usual manner.

    Note that long wavelength filters (>2.5 microns) have a pupil which includes the pupil alignment reference.
    If auto_pupil is set, the pupil will be toggled between CLEAR and CLEARP automatically depending on filter.


    **Spectroscopy:**

    Added in version 0.3 is partial support for the single-object slitless spectroscopy ("SOSS") mode using the
    GR700XD cross-dispersed grating. Currently this includes the clipping of the pupil due to the undersized grating
    and its mounting hardware, and the cylindrical lens that partially defocuses the light in one direction.

    .. warning ::

        Prototype implementation - Not yet fully tested or verified.

    Note that WebbPSF does not model the spectral dispersion in any of NIRISS'
    slitless spectroscopy modes.  For wide-field slitless spectroscopy, this
    can best be simulated by using webbpsf output PSFs as input to the aXe
    spectroscopy code. Contact Van Dixon at STScI for further information.
    For SOSS mode, contact Loic Albert at Universite de Montreal.

    The other two slitless spectroscopy grisms use the regular pupil and do not require any special
    support in WebbPSF.

    """
    SHORT_WAVELENGTH_MIN = 0.6 * 1e-6
    # n.b., the SHORT/LONG distinction in NIRISS is not about
    # different detectors since it only has one of course,
    # rather it's about what's in each of the two wheels.
    SHORT_WAVELENGTH_MAX = LONG_WAVELENGTH_MIN = 2.35 * 1e-6
    LONG_WAVELENGTH_MAX = 5.3 * 1e-6

    def __init__(self, auto_pupil=True):
        self.auto_pupil = auto_pupil
        JWInstrument.__init__(self, "NIRISS")
        self.pixelscale = 0.0656  # Average of X and Y scales, SIAF PRDOPSSOC-031, 2021 April

        self.options['pupil_shift_x'] = 0.0243  # CV3 on-orbit estimate (RPT028027) + OTIS delta from predicted (037134)
        self.options['pupil_shift_y'] = -0.0141

        self.image_mask_list = ['CORON058', 'CORON075', 'CORON150', 'CORON200']  # available but unlikely to be used...
        self.pupil_mask_list = ['CLEARP', 'MASK_NRM', 'GR700XD']

        self._detectors = {'NIS': 'NIS_CEN'}
        self.detector = self.detector_list[0]

    def _addAdditionalOptics(self, optsys, oversample=2):
        """Add NRM or slitless spectroscopy optics for NIRISS.

            These are probably not going to be used much in practice for NIRISS, but they
            are present, so we might as well still provide the ability to simulate 'em.
        """

        from .optics import NIRISS_GR700XD_Grism, NIRISS_CLEARP
        if self.image_mask == 'CORON058':
            radius = 0.58 / 2
            optsys.add_image(function='CircularOcculter', radius=radius, name=self.image_mask)
            trySAM = True
        elif self.image_mask == 'CORON075':
            radius = 0.75 / 2
            optsys.add_image(function='CircularOcculter', radius=radius, name=self.image_mask)
            trySAM = True
        elif self.image_mask == 'CORON150':
            radius = 1.5 / 2
            optsys.add_image(function='CircularOcculter', radius=radius, name=self.image_mask)
            trySAM = True
        elif self.image_mask == 'CORON200':
            radius = 2.0 / 2
            optsys.add_image(function='CircularOcculter', radius=radius, name=self.image_mask)
            trySAM = True
        else:
            trySAM = False
            radius = 0.0  # irrelevant but variable needs to be initialized

        # add pupil plane mask
        shift_x, shift_y = self._get_pupil_shift()
        rotation = self.options.get('pupil_rotation', None)

        # Note - the syntax for specifying shifts is different between FITS files and
        # AnalyticOpticalElement instances. Annoying but historical.
        if self.pupil_mask == 'MASK_NRM':
            optsys.add_pupil(transmission=self._datapath + "/optics/MASK_NRM.fits.gz", name=self.pupil_mask,
                             flip_y=True, flip_x=True, shift_x=shift_x, shift_y=shift_y, rotation=rotation)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'CLEARP':
            optsys.add_pupil(optic=NIRISS_CLEARP(shift_x=shift_x, shift_y=shift_y, rotation=rotation))
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'GR700XD':
            optsys.add_pupil(optic=NIRISS_GR700XD_Grism(shift_x=shift_y, shift_y=shift_y, rotation=rotation))

        elif (self.pupil_mask is None and self.image_mask is not None):
            optsys.add_pupil(name='No Lyot Mask Selected!')

        return (optsys, trySAM, radius + 0.05)  # always attempt to cast this to a SemiAnalyticCoronagraph

    def _get_fits_header(self, hdulist, options):
        """ Format NIRISS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        super(NIRISS, self)._get_fits_header(hdulist, options)

        if self.image_mask is not None:
            hdulist[0].header['CORONPOS'] = (self.image_mask, 'NIRISS coronagraph spot location')
        hdulist[0].header['FOCUSPOS'] = (0, 'NIRISS focus mechanism not yet modeled.')

    @JWInstrument.filter.setter
    def filter(self, value):
        super(NIRISS, self.__class__).filter.__set__(self, value)
        # NIRISS pupils:
        # Short wave filters can be used with a full (clear) pupil
        # long filters have to be used with the CLEARP pupil that contains the
        # PAR reference.

        if self.auto_pupil:
            new_pupil_mask = self.pupil_mask  # default no change
            if self.filter == 'CLEAR':
                # The only science use case for the CLEAR filter position
                # is for GR700XD slitless spectroscopy, so we should set
                # the pupil mask appropriately
                new_pupil_mask = 'GR700XD'
            else:
                wlnum = int(self.filter[1:4])
                if wlnum >= 250:
                    # long wave - can't have clear pupil, it's NRM or GRISM or CLEARP
                    if self.pupil_mask is None:
                        new_pupil_mask = 'CLEARP'
                else:
                    # short wave filter - must have clear pupil
                    new_pupil_mask = None

            if new_pupil_mask != self.pupil_mask:
                _log.info("NIRISS pupil obscuration updated to {0} to match "
                          "the requested filter".format(new_pupil_mask))
                self.pupil_mask = new_pupil_mask

    def _validate_config(self, **kwargs):
        """Validate instrument config for NIRISS

        For NIRISS, this optionally adjusts the instrument pupil
        """
        wavelengths = np.array(kwargs['wavelengths'])
        if np.min(wavelengths) < self.SHORT_WAVELENGTH_MIN:
            raise RuntimeError("The requested wavelengths are too short to be imaged with NIRISS")
        if np.max(wavelengths) > self.LONG_WAVELENGTH_MAX:
            raise RuntimeError("The requested wavelengths are too long to be imaged with NIRISS")
        if (np.max(wavelengths) <= self.SHORT_WAVELENGTH_MAX and
                self.pupil == 'NRM'):
            raise RuntimeError('NRM pupil can only be used with long '
                               'wavelength filters (F277W and longer)')

        return super(NIRISS, self)._validate_config(**kwargs)


class FGS(JWInstrument):
    """ A class modeling the optics of the FGS.

    Not a lot to see here, folks: There are no selectable options, just a great big detector-wide bandpass
    and two detectors.

    The detectors are named as FGS1, FGS2 but may synonymously also be referred to as
    GUIDER1, GUIDER2 for compatibility with DMS convention
    """

    def __init__(self):
        JWInstrument.__init__(self, "FGS")
        self.pixelscale = 0.0691  # Average of X and Y scales for both detectors, SIAF PRDOPSSOC-031, 2021 April

        self.options['pupil_shift_x'] = 0.0041  # CV3 on-orbit estimate (RPT028027) + OTIS delta from predicted (037134)
        self.options['pupil_shift_y'] = -0.0023

        self._detectors = {'FGS1': 'FGS1_FULL', 'FGS2': 'FGS2_FULL'}
        self.detector = self.detector_list[0]

    def _addAdditionalOptics(self, optsys):
        raise NotImplementedError("No user-selectable optics in FGS.")

    def _get_fits_header(self, hdulist, options):
        """ Format FGS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        super(FGS, self)._get_fits_header(hdulist, options)
        hdulist[0].header['FOCUSPOS'] = (0, 'FGS focus mechanism not yet modeled.')

    @JWInstrument.detector.setter # override setter in this subclass
    def detector(self, value):
        # allow either FGS1 or GUIDER1 as synonyms
        if value.upper().startswith('GUIDER'):
            value = 'FGS'+value[-1]
        if value.upper() not in self.detector_list:
            raise ValueError("Invalid detector. Valid detector names are: {}".format(', '.join(self.detector_list)))
        self._detector = value.upper()
        self._update_aperturename()


###########################################################################
# Generic utility functions


def instrument(name):
    """This is just a convenience function, allowing one to access instrument objects based on a string.
    For instance,

    >>> t = instrument('NIRISS')

    Instruments can be referred to either as their full names or as the common three letter abbreviations,
    e.g. "NRC" for NIRCam

    Parameters
    ----------
    name : string
        Name of the instrument class to return. Case insensitive.

    """
    name = name.lower()
    if name == 'miri' or name == 'mir':
        return MIRI()
    elif name == 'nircam' or name == 'nrc':
        return NIRCam()
    elif name == 'nirspec' or name == 'nrs':
        return NIRSpec()
    elif name == 'niriss' or name == 'nis':
        return NIRISS()
    elif name == 'fgs':
        return FGS()
    else:
        raise ValueError("Incorrect instrument name " + name)


instrument.list = ['nircam', 'nirspec', 'niriss', 'miri']  # useful list for iteration


def calc_or_load_PSF(filename, inst, overwrite=False, **kwargs):
    """ Utility function for loading a precomputed PSF from disk, or else
    if that files does not exist, then compute it and save to disk.

    This is useful for writing scripts that cache results - i.e. calculate the
    PSF the first time through and save it, then just load the PSF on subsequent
    iterations.

    Parameters
    ------------
    filename : str
        Filename possibly including path
    inst : JWInstrument
        configured instance of a JWInstrument class
    **kwargs : dict
        Parameters to pass to calc_psf() of that instrument.

    Note that no validation is performed of the PSF loaded from disk to make sure it
    matches the desired properties.  This is just a quick-and-dirty unofficial/undocumented
    helper function.

    """
    if os.path.exists(filename) and not overwrite:
        _log.info("Already exists, no need to recalculate: "+filename)
        return fits.open(filename)
    else:
        return inst.calc_psf(outfile=filename, **kwargs)


#########################


class DetectorGeometry(object):
    """ Utility class for converting between detector coordinates
    in science frame pixels and field of view angular coordinates in arcminutes.


    This is an internal class used within webbpsf; most users will never need to
    interact directly with this class.

    Parameters
    ----------
    siaf : pysiaf.SIAF instance
        Instance of SIAF object for this instrument
    aperturename : string
        Name of SIAF aperture
    shortname : basestring
        Alternate short descriptive name for this aperture

    """

    def __init__(self, siaf, aperturename, shortname=None):
        self.name = aperturename
        if shortname is not None:
            self.name = shortname

        self.mysiaf = siaf
        self.aperture = self.mysiaf[aperturename]

    @property
    def shape(self):
        """ Return detector size in pixels """
        xdetsize = self.aperture.XDetSize
        ydetsize = self.aperture.YDetSize
        return (xdetsize, ydetsize)

    def validate_coords(self, x, y):
        """ Check if specified pixel coords are actually on the detector

        Parameters
        -----------
        x, y : floats
            coordinates in pixels
        """
        if x < 0:
            raise ValueError("Detector pixels X coordinate cannot be negative.")
        if y < 0:
            raise ValueError("Detector pixels Y coordinate cannot be negative.")
        if x > int(self.shape[0]) - 1:
            raise ValueError("Detector pixels X coordinate cannot be > {0}".format(int(self.shape[0]) - 1))
        if y > int(self.shape[1]) - 1:
            raise ValueError("Detector pixels Y coordinate cannot be > {0}".format(int(self.shape[1]) - 1))

    def pix2angle(self, xpix, ypix):
        """ Convert from science frame coordinates (in pixels) to telescope frame coordinates
        (in arcminutes) using SIAF transformations.

        See the pysiaf code for all the full details, or Lallo & Cox Tech Reports

        Parameters
        ------------
        xpix, ypix : floats
            X and Y pixel coordinates, 0 <= xpix, ypix < detector_size

        Returns
        --------
        V2, V3 : floats
            V2 and V3 coordinates, in arcMINUTES
            Note that the astropy.units framework is used to return the result as a
            dimensional Quantity.

        """

        tel_coords = np.asarray(self.aperture.sci_to_tel(xpix, ypix))
        tel_coords_arcmin = tel_coords / 60. * units.arcmin  # arcsec to arcmin

        return tel_coords_arcmin


#########################


def segname(val):
    """ Return WSS-compliant segment name for a variety of input formats

    For instance, one particular segment can be referred to as "B3", 11, "B3-11", etc.
    The WSS refers to this segment as "B3-11".  THis function will return the string
    "B3-11" for any of the above inputs, and similarly for any of the other segments.

    Parameters
    ------------
    val : string or int
        Something that can conceivably be the name or ID of a JWST PMSA.
    """

    try:
        intval = int(val)
        # Convert integer value to string name
        if intval < 1 or intval > 19: raise ValueError('Integer must be between 1 and 19')
        if intval < 7:
            return "A{0}-{0}".format(intval)
        elif intval == 19:
            return "SM-19"
        else:
            letter = 'B' if np.mod(intval, 2) == 1 else 'C'
            number = int(np.ceil((intval - 6) * 0.5))
            return "{0}{1}-{2}".format(letter, number, intval)
    except ValueError:
        # it had better be a letter string
        if val.startswith('SM'): return "SM-19"
        base = {'A': 0, 'B': 5, 'C': 6}
        try:
            offset = base[val[0]]
        except (KeyError, IndexError):
            raise ValueError("string must start with A, B, or C")
        try:
            num = int(val[1])
        except ValueError:
            raise ValueError("input string must have 2nd character as a number from 1-6")
        if num < 1 or num > 6:
            raise ValueError("input string must have 2nd character as a number from 1-6")
        if val[0] == 'A':
            return "{0}{1}-{1}".format(val[0], val[1])
        else:
            return "{0}{1}-{2}".format(val[0], val[1], offset + int(val[1]) * 2)


def one_segment_pupil(segmentname, npix=1024):
    """ Return a pupil image which corresponds to only a single
    segment of the telescope. This can be useful when simulating
    early stages of JWST alignment.


    Example
    -------
    nc = webbpsf.NIRCam()
    nc.pupil = webbpsf.one_segment_pupil('B1')

    """

    # get the master pupil file, which may or may not be gzipped
    segmap = os.path.join(utils.get_webbpsf_data_path(), f"JWpupil_segments_RevW_npix{npix}.fits.gz")
    if not os.path.exists(segmap):
        # try without .gz
        segmap = os.path.join(utils.get_webbpsf_data_path(), f"JWpupil_segments_RevW_npix{npix}.fits")

    newpupil = fits.open(segmap)
    if newpupil[0].header['VERSION'] < 2:
        raise RuntimeError(f"Expecting file version >= 2 for {segmap}")

    segment_official_name = segname(segmentname)
    num = int(segment_official_name.split('-')[1])

    newpupil[0].data = np.asarray(newpupil[0].data == num, dtype=int)

    newpupil[0].header['SEGMENT'] = segment_official_name
    return newpupil
