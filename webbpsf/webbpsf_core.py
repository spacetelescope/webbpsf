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
import functools
from abc import ABC, abstractmethod
import glob
import os
from collections import OrderedDict, namedtuple

import astropy
import astropy.io.ascii as ioascii
import astropy.io.fits as fits
import astropy.units as units
import numpy as np
import poppy
import pysiaf
import scipy.interpolate
import scipy.ndimage

import webbpsf.mast_wss
from webbpsf.utils import label_wavelength

from . import DATA_VERSION_MIN, constants, detectors, distortion, gridded_library, opds, optics, utils

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
    """A generic Space Telescope Instrument class.

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
    configuration can be done by editing the `options` dictionary attribute, either by
    passing options to ``__init__`` or by directly editing the dict afterwards.

    Attributes
    ----------
    telescope : str
        Name of selected telescope, JWST or Roman.
    filter : str
        Bandpass filter name
    image_mask : str
        Name of selected image plane mask, e.g. coronagraph mask or spectrograph slit
    pupil_mask : str
        Name of selected image plane mask, e.g. coronagraph mask or pupil stop
    pupilopd : str
        Filename for telescope pupil wavefront error Optical Path Difference data
    options : dict
        Dictionary for specifying additional specialized options, per each subclass and instance.
    """

    telescope = 'Generic Space Telescope'
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
        """Return the default # of wavelengths to be used for calculation by a given filter"""
        return self._filters[filtername].default_nlambda

    def __init__(self, name='', pixelscale=0.064):
        self.name = name

        self._WebbPSF_basepath, self._data_version = utils.get_webbpsf_data_path(
            data_version_min=DATA_VERSION_MIN, return_version=True
        )

        self._datapath = os.path.join(self._WebbPSF_basepath, self.name)
        self._image_mask = None
        self._pupil_mask = None

        self.pupil = None
        'Filename *or* fits.HDUList for the pupil mask.'
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
        'List of available image_masks'

        self._pupil_mask = None
        self.pupil_mask_list = []
        'List of available pupil_masks'

        self.pixelscale = pixelscale
        'Detector pixel scale, in arcsec/pixel'
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
        if name == '':
            name = None
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
        if name == '':
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
        return '<{telescope}: {instrument_name}>'.format(telescope=self.telescope, instrument_name=self.name)

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
            raise ValueError('Invalid detector. Valid detector names are: {}'.format(', '.join(self.detector_list)))
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
            raise ValueError('Detector pixel coordinates must be pairs of nonnegative numbers, not {}'.format(position))
        if x < 0 or y < 0:
            raise ValueError('Detector pixel coordinates must be nonnegative integers')
        if isinstance(self._detector_npixels, tuple):
            # A tuple has been provided for a non-square detector with different Y and X dimensions
            det_npix_y, det_npix_x = self._detector_npixels
        else:
            det_npix_y = det_npix_x = self._detector_npixels  # same dimensions in both X and Y

        if x > det_npix_x - 1 or y > det_npix_y - 1:
            raise ValueError(
                f'The maximum allowed detector pixel coordinate value is (X,Y) = ({det_npix_x-1}, {det_npix_y-1})'
            )

        self._detector_position = (int(position[0]), int(position[1]))

    @property
    def aperturename(self):
        """SIAF aperture name for detector pixel to sky coords transformations"""
        return self._aperturename

    @aperturename.setter
    def aperturename(self, value):
        # Override in subclass to provide more specific functionality
        self._aperturename = value

    def _update_aperturename(self):
        """Update SIAF aperture name after change in detector or other relevant properties"""
        self.aperturename = self._detectors[self._detector]

    def _get_fits_header(self, result, options):
        """populate FITS Header keywords"""
        super(SpaceTelescopeInstrument, self)._get_fits_header(result, options)
        result[0].header['FILTER'] = (self.filter, 'Filter name')
        if self.image_mask is not None:
            result[0].header['CORONMSK'] = (self.image_mask, 'Image plane mask')
        if self.pupil_mask is not None:
            result[0].header['PUPIL'] = (self.pupil_mask, 'Pupil plane mask')

        result[0].header['VERSION'] = (version, 'WebbPSF software version')
        result[0].header['DATAVERS'] = (self._data_version, 'WebbPSF reference data files version')

        result[0].header['DET_NAME'] = (self.detector, 'Name of detector on this instrument')

        # Correct detector pixel coordinates to allow for even arrays to be centered on half pixel boundary
        dpos = np.asarray(self.detector_position, dtype=float)
        oversamp = result[0].header['OVERSAMP']
        size = result[0].data.shape[0]

        if size / oversamp % 2 == 0:
            dpos += 0.5  # even arrays must be at a half pixel

        result[0].header['DET_X'] = (dpos[0], 'Detector X pixel position of array center')
        result[0].header['DET_Y'] = (dpos[1], 'Detector Y pixel position of array center')

        for key in self._extra_keywords:
            result[0].header[key] = self._extra_keywords[key]

    def get_optical_system(self, fft_oversample=2, detector_oversample=None, fov_arcsec=2, fov_pixels=None, options=None):
        """Return an OpticalSystem instance corresponding to the instrument as currently configured.

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

        _log.info('Creating optical system model:')

        self._extra_keywords = OrderedDict()  # Place to save info we later want to put into the FITS header for each PSF.

        if options is None:
            options = self.options
        if detector_oversample is None:
            detector_oversample = fft_oversample

        _log.debug('Oversample: %d  %d ' % (fft_oversample, detector_oversample))
        optsys = poppy.OpticalSystem(
            name='{telescope}+{instrument}'.format(telescope=self.telescope, instrument=self.name), oversample=fft_oversample
        )
        # For convenience offsets can be given in cartesian or radial coords
        if 'source_offset_x' in options or 'source_offset_y' in options:
            offx = options.get('source_offset_x', 0)
            offy = options.get('source_offset_y', 0)
            optsys.source_offset_r = np.sqrt(offx**2 + offy**2)
            optsys.source_offset_theta = np.rad2deg(np.arctan2(-offx, offy))
            _log.debug(
                'Source offset from X,Y = ({}, {}) is (r,theta) = {},{}'.format(
                    offx, offy, optsys.source_offset_r, optsys.source_offset_theta
                )
            )
        if 'source_offset_r' in options:
            optsys.source_offset_r = options['source_offset_r']
        if 'source_offset_theta' in options:
            optsys.source_offset_theta = options['source_offset_theta']

        # Telescope entrance pupil
        pupil_optic = self._get_telescope_pupil_and_aberrations()
        optsys.add_pupil(pupil_optic)

        pupil_rms_wfe_nm = np.sqrt(np.mean(pupil_optic.opd[pupil_optic.amplitude == 1] ** 2)) * 1e9
        self._extra_keywords['TEL_WFE'] = (float(pupil_rms_wfe_nm), '[nm] Telescope pupil RMS wavefront error')
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
                self._extra_keywords['SI_WFE'] = (float(inst_rms_wfe_nm), '[nm] instrument pupil RMS wavefront error')
            except (TypeError, IndexError):
                # Currently the above does not work for Roman, but fixing this is deferred to future work
                pass

            if hasattr(aberration_optic, 'header_keywords'):
                self._extra_keywords.update(aberration_optic.header_keywords())

        # ---- Add defocus if requested
        if 'defocus_waves' in options:
            defocus_waves = options['defocus_waves']
            defocus_wavelength = float(options['defocus_wavelength']) if 'defocus_wavelength' in options else 2.0e-6
            _log.info(f'Adding defocus of {defocus_waves:.3f} waves at {defocus_wavelength*1e6:.3f} microns')
            lens = poppy.ThinLens(
                name='Defocus', nwaves=defocus_waves, reference_wavelength=defocus_wavelength, radius=self.pupil_radius
            )
            optsys.add_pupil(optic=lens)
            self._extra_keywords['DEFOCUS'] = (defocus_waves, '# of waves of defocus added')
            self._extra_keywords['DEFOC_WL'] = (defocus_wavelength, 'Wavelength reference for defocus added')

        # ---- add coronagraph or spectrograph optics if requested,
        # and possibly flag to invoke semi-analytic coronagraphic propagation

        # first error check for null strings, which should be considered like None
        if self.image_mask == '':
            self.image_mask = None
        if self.pupil_mask == '':
            self.pupil_mask = None

        if (
            self.image_mask is not None
            or self.pupil_mask is not None
            or 'WL' in self.filter  # special case handling for NIRCam WLP4 filter that is also a lens
            or ('force_coron' in options and options['force_coron'])
        ):
            _log.debug('Adding coronagraph/spectrograph optics...')
            optsys, trySAM, SAM_box_size = self._addAdditionalOptics(optsys, oversample=fft_oversample)
        else:
            trySAM = False

        # --- add the detector element.
        if fov_pixels is None:
            if not np.isscalar(fov_arcsec):
                fov_arcsec = np.asarray(fov_arcsec)  # cast to ndarray if 2D
            fov_pixels = np.round(fov_arcsec / self.pixelscale)
            if 'parity' in options:
                if options['parity'].lower() == 'odd' and np.remainder(fov_pixels, 2) == 0:
                    fov_pixels += 1
                if options['parity'].lower() == 'even' and np.remainder(fov_pixels, 2) == 1:
                    fov_pixels += 1
        else:
            pass

        optsys.add_detector(
            self.pixelscale,
            fov_pixels=fov_pixels,
            oversample=detector_oversample,
            name=self.name + ' detector'
        )

        # ---  invoke semi-analytic coronagraphic propagation
        if trySAM and not ('no_sam' in self.options and self.options['no_sam']):
            # if this flag is set, try switching to SemiAnalyticCoronagraph mode.
            _log.info('Trying to invoke switch to Semi-Analytic Coronagraphy algorithm')
            try:
                SAM_optsys = poppy.SemiAnalyticCoronagraph(optsys, oversample=fft_oversample, occulter_box=SAM_box_size)
                _log.info('SAC OK')
                return SAM_optsys
            except ValueError as err:
                _log.warning(
                    'Could not switch to Semi-Analytic Coronagraphy mode; invalid set of optical planes? '
                    'Using default propagation instead.'
                )
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
            opd_map = self.pupilopd if os.path.exists(self.pupilopd) else os.path.join(self._datapath, 'OPD', self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], str):
            # tuple with filename and slice
            opd_map = (
                self.pupilopd[0]
                if os.path.exists(self.pupilopd[0])
                else os.path.join(self._datapath, 'OPD', self.pupilopd[0]),
                self.pupilopd[1],
            )
        elif isinstance(self.pupilopd, (fits.HDUList, poppy.OpticalElement)):
            opd_map = self.pupilopd  # not a path per se but this works correctly to pass it to poppy
        elif self.pupilopd is None:
            opd_map = None
        else:
            raise TypeError('Not sure what to do with a pupilopd of that type:' + str(type(self.pupilopd)))

        # ---- set pupil intensity
        if self.pupil is None:
            raise RuntimeError('The pupil shape must be specified in the ' 'instrument class or by setting self.pupil')
        if isinstance(self.pupil, poppy.OpticalElement):
            # supply to POPPY as-is
            pupil_optic = self.pupil
        else:
            # wrap in an optic and supply to POPPY
            if isinstance(self.pupil, str):  # simple filename
                if os.path.exists(self.pupil):
                    pupil_transmission = self.pupil
                else:
                    pupil_transmission = os.path.join(self._WebbPSF_basepath, self.pupil)
            elif isinstance(self.pupil, fits.HDUList):
                # POPPY can use self.pupil as-is
                pupil_transmission = self.pupil
            else:
                raise TypeError('Not sure what to do with a pupil of ' 'that type: {}'.format(type(self.pupil)))
            # ---- apply pupil intensity and OPD to the optical model
            pupil_optic = poppy.FITSOpticalElement(
                name='{} Entrance Pupil'.format(self.telescope),
                transmission=pupil_transmission,
                opd=opd_map,
                planetype=poppy.poppy_core.PlaneType.pupil,
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
        raise NotImplementedError('needs to be subclassed.')

    def _get_synphot_bandpass(self, filtername):
        """Return a synphot.spectrum.SpectralElement object for the given desired band.

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
            _log.warning(
                'The supplied file, {}, does not have a WAVEUNIT keyword. Assuming it ' 'is Angstroms.'.format(
                    filter_info.filename
                )
            )
            waveunit = 'angstrom'

        filterdata = filterfits[1].data
        try:
            band = synphot.SpectralElement(
                synphot.models.Empirical1D, points=filterdata.WAVELENGTH, lookup_table=filterdata.THROUGHPUT, keep_neg=False
            )

        except AttributeError:
            raise ValueError(
                'The supplied file, %s, does not appear to be a FITS table '
                'with WAVELENGTH and THROUGHPUT columns.' % filter_info.filename
            )

        filterfits.close()
        return band

    def psf_grid(
        self,
        num_psfs=16,
        all_detectors=True,
        save=False,
        outdir=None,
        outfile=None,
        overwrite=True,
        verbose=True,
        use_detsampled_psf=False,
        single_psf_centered=True,
        **kwargs,
    ):
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

        Examples
        --------
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
            detectors = 'all'
        else:
            detectors = self.detector

        if single_psf_centered is True:
            if isinstance(self._detector_npixels, tuple):
                # A tuple has been provided for a non-square detector with different Y and X dimensions
                det_npix_y, det_npix_x = self._detector_npixels
            else:
                det_npix_y = det_npix_x = self._detector_npixels  # same dimensions in both X and Y
            psf_location = (int(det_npix_x - 1) // 2, int(det_npix_y - 1) // 2)  # center pt
        else:
            psf_location = self.detector_position[::-1]  # (y,x)

        # Call CreatePSFLibrary class
        inst = gridded_library.CreatePSFLibrary(
            instrument=self,
            filter_name=filt,
            detectors=detectors,
            num_psfs=num_psfs,
            psf_location=psf_location,
            use_detsampled_psf=use_detsampled_psf,
            save=save,
            outdir=outdir,
            filename=outfile,
            overwrite=overwrite,
            verbose=verbose,
            **kwargs,
        )
        gridmodel = inst.create_grid()

        return gridmodel


#  JWInstrument classes  #####


@utils.combine_docstrings
class JWInstrument(SpaceTelescopeInstrument):
    """Superclass for all JWST instruments

    Attributes
    ----------
    telescope : str
        name of telescope
    pupilopd : file-like
        filename or FITS file object for the pupil Optical Path Difference
    include_si_wfe : boolean
        Should SI internal WFE be included in models? Requires
        the presence of ``si_zernikes_isim_cv3.fits`` in the
        ``WEBBPSF_PATH``. Default = True.
    """

    telescope = 'JWST'
    pupilopd = None
    """Filename *or* fits.HDUList for JWST pupil OPD.

    This can be either a full absolute filename, or a relative name in which
    case it is assumed to be within the instrument's `data/OPDs/` directory,
    or an actual fits.HDUList object corresponding to such a file. If the file
    contains a datacube, you may set this to a tuple (filename, slice) to
    select a given slice, or else the first slice will be used."""

    def __init__(self, *args, **kwargs):
        super(JWInstrument, self).__init__(*args, **kwargs)

        self.siaf = get_siaf_with_caching(self.name)

        opd_path = os.path.join(self._datapath, 'OPD')
        self.opd_list = []
        for filename in glob.glob(os.path.join(opd_path, 'OPD*.fits*')):
            self.opd_list.append(os.path.basename(os.path.abspath(filename)))
        for filename in glob.glob(os.path.join(self._WebbPSF_basepath, 'JWST_OTE_OPD*.fits*')):
            self.opd_list.append(os.path.basename(os.path.abspath(filename)))

        if not len(self.opd_list) > 0:
            raise RuntimeError('No pupil OPD files found for {name} in {path}'.format(name=self.name, path=opd_path))

        self.opd_list.sort()
        self.pupilopd = 'JWST_OTE_OPD_cycle1_example_2022-07-30.fits'  # Default is now an on-orbit measured example OPD

        self.pupil = os.path.abspath(os.path.join(self._WebbPSF_basepath, 'jwst_pupil_RevW_npix1024.fits.gz'))
        'Filename *or* fits.HDUList for JWST pupil mask. Usually there is no need to change this.'

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
        """Return default FOV in arcseconds"""
        return 5  # default for all NIR instruments

    def get_optical_system(self, fft_oversample=2, detector_oversample=None, fov_arcsec=2, fov_pixels=None, options=None):
        # invoke superclass version of this
        # then add a few display tweaks
        optsys = SpaceTelescopeInstrument.get_optical_system(
            self,
            fft_oversample=fft_oversample,
            detector_oversample=detector_oversample,
            fov_arcsec=fov_arcsec,
            fov_pixels=fov_pixels,
            options=options,
        )
        # If the OTE model in the entrance pupil is a plain FITSOpticalElement, cast it to the linear model class
        if not isinstance(optsys.planes[0], opds.OTE_Linear_Model_WSS):
            lom_ote = opds.OTE_Linear_Model_WSS()
            # FIXME seems like some code is missing here...? But in practice this code path
            # never gets executed due to the _get_telescope_pupil_and_aberrations() function doing the right thing.
            lom_ote

        optsys.planes[0].display_annotate = utils.annotate_ote_pupil_coords
        return optsys

    def _get_aberrations(self):
        """return OpticalElement modeling wavefront aberrations for a given instrument,
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
            return os.path.join(self._datapath, 'OPD', opdfilename)
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
            raise TypeError('Not sure what to do with a pupilopd of that type:' + str(type(self.pupilopd)))

        # ---- set pupil intensity
        if self.pupil is None:
            raise RuntimeError('The pupil shape must be specified in the ' 'instrument class or by setting self.pupil')
        if isinstance(self.pupil, poppy.OpticalElement):
            # supply to POPPY as-is
            pupil_optic = self.pupil
        else:
            # wrap in an optic and supply to POPPY
            if isinstance(self.pupil, str):  # simple filename
                if os.path.exists(self.pupil):
                    pupil_transmission = self.pupil
                else:
                    pupil_transmission = os.path.join(self._WebbPSF_basepath, self.pupil)
                # Get npix from pupil_transmission
                npix = int(pupil_transmission.split('npix')[-1].split('.')[0])
            elif isinstance(self.pupil, fits.HDUList):
                # POPPY can use self.pupil as-is
                pupil_transmission = self.pupil
                # Get npix from the shape of the data
                npix = self.pupil[0].data.shape[0]
            else:
                raise TypeError('Not sure what to do with a pupil of ' 'that type: {}'.format(type(self.pupil)))

            # ---- apply pupil intensity and OPD to the optical model
            pupil_optic = opds.OTE_Linear_Model_WSS(
                name='{} Entrance Pupil'.format(self.telescope),
                transmission=pupil_transmission,
                opd=opd_map,
                opd_index=opd_index,
                v2v3=self._tel_coords(),
                npix=npix,
                include_nominal_field_dependence=self.include_ote_field_dependence,
            )

        return pupil_optic

    @SpaceTelescopeInstrument.aperturename.setter
    def aperturename(self, value):
        """Set SIAF aperture name to new value, with validation.

        This also updates the pixelscale to the local value for that aperture, for a small precision enhancement.
        """
        # Explicitly update detector reference coordinates to the default for the new selected aperture,
        # otherwise old coordinates can persist under certain circumstances

        try:
            ap = self.siaf[value]
        except KeyError:
            raise ValueError(f'Aperture name {value} not a valid SIAF aperture name for {self.name}')

        # Only update if new value is different
        if self._aperturename != value:
            if ap.AperType == 'SLIT':
                # Special case for SLIT apertures (NIRSpec or MIRI). Note, this includes all IFU apertures for NIRSpec
                # apertures of type SLIT define V2,V3 position, but not pixel coordinates and pixelscale. So we
                # still have to use a full-detector aperturename for that subset of apertures
                # This code path also supports MIRI LRS
                detector_apername = self.detector + '_FULL'
                _log.info(f'Aperture {value} is of type SLIT; using {detector_apername} for detector geometry.')

                has_custom_pixelscale = self._aperturename and (
                    self.pixelscale != self._get_pixelscale_from_apername(detector_apername)
                )

                # Now apply changes:
                self._aperturename = value

                # Update DetectorGeometry class
                self._detector_geom_info = DetectorGeometry(self.siaf, self._aperturename)
                _log.info(
                    f'{self.name} SIAF aperture name updated to {self._aperturename} using geometry from {detector_apername}'
                )
                if not has_custom_pixelscale:
                    self.pixelscale = self._get_pixelscale_from_apername(detector_apername)
                    debug_message = (
                        f'Pixelscale updated to {self.pixelscale} '
                        f'based on average X+Y SciScale at SIAF aperture {detector_apername}'
                    )
                    _log.debug(debug_message)
            elif ap.AperType == 'COMPOUND' and self.name == 'MIRI':
                # For MIRI, many of the relevant IFU apertures are of COMPOUND type.
                has_custom_pixelscale = False  # custom scales not supported for MIRI IFU (yet?)
                # Unlike NIRSpec, there simply do not exist full-detector SIAF apertures for the MIRI IFU detectors
                info_message = (
                    f'Aperture {value} is of type COMPOUND for MIRI; '
                    'There do not exist corresponding SIAF apertures, so we ignore setting detector geometry.'
                )
                _log.info(info_message)

                # Now apply changes:
                self._aperturename = value
                # Update DetectorGeometry class
                self._detector_geom_info = DetectorGeometry(self.siaf, self._aperturename)
                if not has_custom_pixelscale:
                    self.pixelscale = self._get_pixelscale_from_apername(value)
                    _log.debug(f'Pixelscale updated to {self.pixelscale} based on IFU cubepars for {value}')

            else:
                if self.detector not in value:
                    error_message = (
                        f'Aperture name {value} does not match currently selected detector {self.detector}. '
                        f'Change detector attribute first, then set desired aperture.'
                    )
                    raise ValueError(error_message)

                # First, check some info from current settings, wich we will use below as part of auto pixelscale code
                # The point is to check if the pixel scale is set to a custom or default value,
                # and if it's custom then don't override that.
                # Note:
                #     check self._aperturename first to account for the edge case when
                #     this is called from __init__ before _aperturename is set
                #     and also check first that it's not a SLIT type aperture,
                #     for which the usual _get_pixelscale_from_apername won't work.
                #     and also check neither current nor requested aperture are of type
                #     SLIT since that doesn't have a pixelscale to get.
                has_custom_pixelscale = (
                    self._aperturename
                    and (self.siaf[self._aperturename].AperType != 'SLIT')
                    and (self.pixelscale != self._get_pixelscale_from_apername(self._aperturename))
                    and ap.AperType != 'SLIT'
                )

                # Now apply changes:
                self._aperturename = value
                # Update detector reference coordinates
                self.detector_position = (ap.XSciRef, ap.YSciRef)

                # Update DetectorGeometry class
                self._detector_geom_info = DetectorGeometry(self.siaf, self._aperturename)
                _log.info(f'{self.name} SIAF aperture name updated to {self._aperturename}')

                if not has_custom_pixelscale:
                    self.pixelscale = self._get_pixelscale_from_apername(self._aperturename)
                    _log.debug(
                        f'Pixelscale updated to {self.pixelscale}' +
                        f'based on average X+Y SciScale at SIAF aperture {self._aperturename}'
                    )

    def _tel_coords(self):
        """Convert from science frame coordinates to telescope frame coordinates using
        SIAF transformations. Returns (V2, V3) tuple, in arcminutes.

        Note that the astropy.units framework is used to return the result as a
        dimensional Quantity.
        """

        if self._detector_geom_info.aperture.AperType == 'SLIT':
            # These apertures don't map directly to particular detector position in the usual way
            # Return coords for center of the aperture reference location
            return (
                np.asarray((self._detector_geom_info.aperture.V2Ref, self._detector_geom_info.aperture.V3Ref))
                / 60
                * units.arcmin
            )
        elif self._detector_geom_info.aperture.AperType == 'COMPOUND':
            # handle MIRI MRS apertures, which don't have V2Ref,V3Ref defined, but this works:
            return np.asarray(self.siaf[self.aperturename].reference_point('tel')) / 60 * units.arcmin
        else:
            return self._detector_geom_info.pix2angle(self.detector_position[0], self.detector_position[1])

    def _xan_yan_coords(self):
        """Convert from detector pixel coordinates to the XAN, YAN coordinate system
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
        """Set the simulated center point of the array based on a named SIAF aperture.
        This will adjust the detector and detector position attributes.
        """
        try:
            ap = self.siaf[aperture_name]

            # setting the detector must happen -before- we set the position
            detname = aperture_name.split('_')[0]

            if detname.startswith('MIRIFU'):
                self._mode = 'IFU'
                if 'CHANNEL1' in aperture_name or 'CHANNEL2' in aperture_name:
                    self.detector = 'MIRIFUSHORT'
                else:
                    self.detector = 'MIRIFULONG'
            elif detname.startswith('MIRIM'):
                self._mode = 'imaging'
                self.detector = detname  # As a side effect this auto reloads SIAF info, see detector.setter

            elif detname != 'NRS':  # Many NIRSpec slit apertures are defined generally, not for a specific detector
                self.detector = detname  # As a side effect this auto reloads SIAF info, see detector.setter

            self.aperturename = aperture_name

            if (self.name == 'NIRSpec' or self.name == 'MIRI') and ap.AperType == 'SLIT':
                # NIRSpec slit apertures need some separate handling, since they don't map directly to detector pixels
                # In this case the detector position is not uniquely defined, but we ensure to get reasonable values by
                # using one of the full-detector NIRspec apertures
                # This code path also supports MIRI LRS SLIT
                _log.debug(f'Inferring detector position using V coords for SLIT aperture: {ap.V2Ref, ap.V3Ref}')
                ref_in_tel = ap.V2Ref, ap.V3Ref
                nrs_full_aperture = self.siaf[self.detector + '_FULL']
                ref_in_sci = nrs_full_aperture.tel_to_sci(*ref_in_tel)
                self.detector_position = ref_in_sci
            elif self.name == 'MIRI' and ap.AperType == 'COMPOUND':
                # MIRI IFU compound apertures need separate handling, since they don't map directoy to detector pixels
                # In this case the detector position is not uniquely defined, and there do not exist
                # in SIAF any full-detector MIRIFU apertures, so just set values to (512,512) as a placeholder.
                self.detector_position = [512, 512]
            else:
                # Regular imaging apertures, so we can just look up the reference coords directly
                self.detector_position = (ap.XSciRef, ap.YSciRef)  # set this AFTER the SIAF reload

            _log.debug('From {} set det. pos. to {} {}'.format(aperture_name, detname, self.detector_position))

        except KeyError:
            raise ValueError('Not a valid aperture name for {}: {}'.format(self.name, aperture_name))

    def _get_pixelscale_from_apername(self, apername):
        """Simple utility function to look up pixelscale from apername"""
        ap = self.siaf[apername]
        # Here we make the simplifying assumption of **square** pixels, which is true within 0.5%.
        # The slight departures from this are handled in the distortion model; see distortion.py
        return (ap.XSciScale + ap.YSciScale) / 2

    def _get_fits_header(self, result, options):
        """populate FITS Header keywords"""
        super(JWInstrument, self)._get_fits_header(result, options)

        # Add JWST-specific V2,V3 focal plane coordinate system.
        v2v3pos = self._tel_coords()
        result[0].header.insert(
            'DET_Y', ('DET_V2', v2v3pos[0].value, '[arcmin] Det. pos. in telescope V2,V3 coord sys'), after=True
        )
        result[0].header.insert(
            'DET_V2', ('DET_V3', v2v3pos[1].value, '[arcmin] Det. pos. in telescope V2,V3 coord sys'), after=True
        )
        result[0].header['APERNAME'] = (self._aperturename, 'SIAF aperture name')

    def calc_psf(
        self,
        outfile=None,
        source=None,
        nlambda=None,
        monochromatic=None,
        fov_arcsec=None,
        fov_pixels=None,
        oversample=None,
        detector_oversample=None,
        fft_oversample=None,
        overwrite=True,
        display=False,
        save_intermediates=False,
        return_intermediates=False,
        normalize='first',
        add_distortion=True,
        crop_psf=True,
    ):
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
        psf = SpaceTelescopeInstrument.calc_psf(
            self,
            outfile=outfile,
            source=source,
            nlambda=nlambda,
            monochromatic=monochromatic,
            fov_arcsec=fov_arcsec,
            fov_pixels=fov_pixels,
            oversample=oversample,
            detector_oversample=detector_oversample,
            fft_oversample=fft_oversample,
            overwrite=overwrite,
            display=display,
            save_intermediates=save_intermediates,
            return_intermediates=return_intermediates,
            normalize=normalize,
        )

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
        # you can turn on/off IPC corrections via the add_ipc option, default True.
        add_ipc = options.get('add_ipc', True)

        # Add distortion if set in calc_psf
        if add_distortion:
            _log.info('Adding PSF distortion(s) and detector effects')

            # Set up new extensions to add distortion to:
            n_exts = len(result)
            for ext in np.arange(n_exts):
                hdu_new = fits.ImageHDU(result[ext].data, result[ext].header)  # these will be the PSFs that are edited
                result.append(hdu_new)
                ext_new = ext + n_exts
                result[ext_new].header['EXTNAME'] = result[ext].header['EXTNAME'][0:4] + 'DIST'  # change extension name
                _log.debug('Appending new extension {} with EXTNAME = {}'.format(ext_new, result[ext_new].header['EXTNAME']))

            # Apply optical geometric distortions and detector systematic effects based on the instrument
            if self.name in ['NIRCam', 'NIRISS', 'FGS']:
                # Apply distortion effects: Rotation and optical distortion
                _log.debug('NIRCam/NIRISS/FGS: Adding rotation and optical distortion')
                psf_rotated = distortion.apply_rotation(result, crop=crop_psf)  # apply rotation
                psf_siaf_distorted = distortion.apply_distortion(psf_rotated)  # apply siaf distortion model
                psf_distorted = detectors.apply_detector_charge_diffusion(
                    psf_siaf_distorted, options
                )  # apply detector charge transfer model
            elif self.name == 'MIRI':
                # Apply distortion effects to MIRI psf: Distortion and MIRI Scattering
                _log.debug('MIRI: Adding optical distortion and Si:As detector internal scattering')
                if self.mode != 'IFU':
                    if self._detector_geom_info.aperture.AperType != 'SLIT':
                        psf_siaf = distortion.apply_distortion(result)  # apply siaf distortion
                    else:
                        # slit type aperture, specifically LRS SLIT, does not have distortion polynomials
                        # therefore omit apply_distortion if a SLIT aperture is selected.
                        psf_siaf = result
                    psf_siaf_rot = detectors.apply_miri_scattering(psf_siaf)  # apply scattering effect
                    psf_distorted = detectors.apply_detector_charge_diffusion(
                        psf_siaf_rot, options
                    )  # apply detector charge transfer model
                else:
                    # there is not yet any distortion calibration for the IFU, and
                    # we don't want to apply charge diffusion directly here
                    psf_distorted = detectors.apply_miri_ifu_broadening(result, options, slice_width=self._ifu_slice_width)
            elif self.name == 'NIRSpec':
                # Apply distortion effects to NIRSpec psf: Distortion only
                # (because applying detector effects would only make sense after simulating spectral dispersion)
                _log.debug('NIRSpec: Adding optical distortion')
                if self.mode != 'IFU':
                    psf_siaf = distortion.apply_distortion(result)  # apply siaf distortion model
                    psf_distorted = detectors.apply_detector_charge_diffusion(
                        psf_siaf, options
                    )  # apply detector charge transfer model

                else:
                    # there is not yet any distortion calibration for the IFU.
                    psf_distorted = detectors.apply_nirspec_ifu_broadening(result, options)

            # Edit the variable to match if input didn't request distortion
            # (cannot set result = psf_distorted due to return method)
            [result.append(fits.ImageHDU()) for i in np.arange(len(psf_distorted) - len(result))]
            for ext in np.arange(len(psf_distorted)):
                result[ext] = psf_distorted[ext]

        _log.info('Formatting output FITS extensions including for sampling.')
        # Rewrite result variable based on output_mode; this includes binning down to detector sampling.
        SpaceTelescopeInstrument._calc_psf_format_output(self, result, options)

        if add_ipc and add_distortion and ('DET_DIST' in result):
            result = detectors.apply_detector_ipc(result)  # apply detector IPC model (after binning to detector sampling)
        if add_ipc and add_distortion and ('OVERDIST' in result):
            result = detectors.apply_detector_ipc(result, extname='OVERDIST')  # apply detector IPC model to oversampled PSF

    def interpolate_was_opd(self, array, newdim):
        """Interpolates an input 2D  array to any given size.

        Parameters
        ----------
        array: float
             input array to interpolate
        newdim: int
             new size of the 2D square array (newdim x newdim)

        Returns
        -------
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
        """Return a tuple of pupil shifts, for passing to OpticalElement constructors
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
        if ('pupil_shift_x' in self.options and self.options['pupil_shift_x'] != 0) or (
            'pupil_shift_y' in self.options and self.options['pupil_shift_y'] != 0
        ):
            from .constants import JWST_CIRCUMSCRIBED_DIAMETER

            # missing values are treated as 0's
            shift_x = self.options.get('pupil_shift_x', 0)
            shift_y = self.options.get('pupil_shift_y', 0)
            # nones are likewise treated as 0's
            if shift_x is None:
                shift_x = 0
            if shift_y is None:
                shift_y = 0
            # Apply pupil scale
            shift_x *= JWST_CIRCUMSCRIBED_DIAMETER
            shift_y *= JWST_CIRCUMSCRIBED_DIAMETER
            _log.info('Setting Lyot pupil shift to ({}, {})'.format(shift_x, shift_y))
        else:
            shift_x, shift_y = None, None
        return shift_x, shift_y

    def _apply_jitter(self, result, local_options=None):
        """Modify a PSF to account for the blurring effects of image jitter.
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
        ----------
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

        _log.info('Calculating jitter using ' + str(local_options['jitter']))

        def _linear_smear(smear_length, image):
            # Helper function, used below
            smear_length_pix = int(np.round(smear_length / result[0].header['PIXELSCL']))
            if smear_length_pix % 2 == 0:
                smear_length_pix += 1  # Astropy convolution requires odd sized kernels only

            smear_model = np.identity(smear_length_pix)
            _log.info(
                'Jitter: Convolving with linear smear of {0:.3f} arcsec; {1:d} pixels'.format(smear_length, smear_length_pix)
            )
            kern = astropy.convolution.kernels.CustomKernel(smear_model)
            return astropy.convolution.convolve_fft(image, kern, allow_huge=True)

        if local_options['jitter'] is None:
            return
        elif local_options['jitter'].lower() == 'gaussian':
            # Regular version in poppy
            return super()._apply_jitter(result, local_options=local_options)
        elif local_options['jitter'].lower() == 'linear':
            # Drift by 0.12 arcsec (1 mas/second for 2 minutes)

            smear_length = 0.12  # arcsec

            out = _linear_smear(smear_length, result[0].data)
            result[0].header['JITRTYPE'] = ('Linear smear / drift', 'Type of jitter applied')
            result[0].header['JITSMEAR'] = (smear_length, 'Linear smear [arcsec]')

        elif local_options['jitter'].lower() == 'pcs=coarse':
            # JWST coarse point, current best estimate based on high fidelity monte carlo sims by Peiman Maghami

            cp_case = local_options.get('jitter_coarse_model_case', 2)  # Coarse pointing model case, 1 or 2
            exp_duration = local_options.get('exp_duration', 75)  # Duration in seconds
            exp_start_time = local_options.get('exp_start_time', 0)  # Start time in seconds

            offset, kernel = opds.get_coarse_blur_parameters(
                exp_start_time, exp_duration, result[0].header['PIXELSCL'], case=cp_case
            )

            kern = astropy.convolution.kernels.CustomKernel(kernel)
            out = astropy.convolution.convolve_fft(result[0].data, kern, allow_huge=True)

            result[0].header['JITRTYPE'] = ('PCS Coarse, high fidelity MC model results', 'Type of jitter applied')
            result[0].header['JITRCASE'] = (cp_case, 'PCS Coarse mode: Monte Carlo model case used')
            result[0].header['JITR_T0'] = (exp_start_time, 'PCS Coarse mode: sim exposure start time [s]')
            result[0].header['JITRTEXP'] = (exp_duration, 'PCS Coarse mode: sim exposure duration [s]')
            result[0].header['JITRCPV2'] = (offset[0], 'Coarse pointing offset in V2 [arcsec]')
            result[0].header['JITRCPV3'] = (offset[1], 'Coarse pointing offset in V3 [arcsec]')

        elif local_options['jitter'].lower() == 'pcs=coarse_like_itm':
            # JWST coarse point, assumptions in ITM
            # Acton says:
            #  it is actually 0.4 for a boresight error, 0.4 smear, and 0.2 jitter.
            #  Boresight error is a random term for image placement,
            #  smear is mostly a linear uniform blur, and jitter is gaussian.

            # First we do the fast jitter part
            local_options['jitter_sigma'] = 0.2
            import scipy.ndimage

            sigma = local_options.get('jitter_sigma')

            # that will be in arcseconds, we need to convert to pixels:
            _log.info('Jitter: Convolving with Gaussian with sigma={0:.3f} arcsec'.format(sigma))
            out = scipy.ndimage.gaussian_filter(result[0].data, sigma / result[0].header['PIXELSCL'])

            # Now we'll do the linear jitter part
            smear_length = 0.4  # arcsec
            out = _linear_smear(smear_length, out)

            result[0].header['JITRTYPE'] = ('PCS Coarse, like ITM', 'Type of jitter applied')
            result[0].header['JITRSIGM'] = (sigma, 'Gaussian sigma for jitter, per axis [arcsec]')
            result[0].header['JITSMEAR'] = (smear_length, 'Linear smear [arcsec]')

        elif local_options['jitter'].lower() == 'custom':
            # User-supplied arbitrary PSF convolution kernel

            if ('jitter_kernel' not in local_options) or (not local_options['jitter_kernel'].ndim == 2):
                raise ValueError("You must supply an .options['jitter_kernel'] 2D array to use the custom jitter option")
            _log.info('Jitter: Convolving with user-supplied custom convolution kernel')
            kern = astropy.convolution.kernels.CustomKernel(local_options['jitter_kernel'])
            out = astropy.convolution.convolve_fft(result[0].data, kern, allow_huge=True)

            result[0].header['JITRTYPE'] = ('Custom jitter kernel', 'Type of jitter applied')

        else:
            raise ValueError('Unknown jitter option value: ' + local_options['jitter'])

        peak = result[0].data.max()
        newpeak = out.max()
        strehl = newpeak / peak  # not really the whole Strehl ratio, just the part due to jitter
        _log.info('        resulting image peak drops to {0:.3f} of its previous value'.format(strehl))
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

        if kind.lower() == 'total':
            # recursively get total OPD including SI plus OTE
            opd = self.get_wfe('ote') + self.get_wfe('si')
        elif kind.lower() == 'si':
            aberration = self._get_aberrations()
            opd = aberration.get_opd(wave)
            if self.name.lower() == 'nirspec':
                # For NIRSpec, the WFE is normally allocated to 1/3 before the MSA and 2/3 after the MSA.
                # The call to get_aberrations above just returns the foreoptics portion.
                # Multiply by 3x to get the total instrumental WFE.
                opd *= 3
            # Flip vertically to match OTE entrance pupil orientation
            opd = np.flipud(opd)
        elif kind.lower() == 'ote':  # OTE *total* WFE including all terms
            opd = ote.get_opd(wave).copy()
            aperture = ote.get_transmission(wave)
            opd *= aperture != 0  # mask out to zero the global zernikes outside the aperture

        elif kind.lower() == 'ote_global':  # OTE *global* WFE only, i.e. WFE common to all field points
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
            raise NotImplementedError(f'Not yet implemented: {kind}')
        else:
            raise NotImplementedError(f'Not a known kind of WFE: {kind}')

        if plot:
            import matplotlib
            import matplotlib.pyplot as plt

            plt.imshow(opd, vmin=-5e-7, vmax=5e-7, cmap=matplotlib.cm.RdBu_r, origin='lower')
            plt.title(kind + ' WFE')
            mask = ote.get_transmission(wave) != 0
            plt.xlabel(f'RMS: {utils.rms(opd, mask)*1e9:.2f} nm')
            plt.colorbar(label='WFE [m]')

        return opd

    def visualize_wfe_budget(self, slew_delta_time=14 * units.day, slew_case='EOL', ptt_only=False, verbose=True):
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

        webbpsf.optical_budget.visualize_wfe_budget(
            self, slew_delta_time=slew_delta_time, slew_case=slew_case, ptt_only=ptt_only, verbose=verbose
        )

    def load_wss_opd(self, filename, output_path=None, backout_si_wfe=True, verbose=True, plot=False, save_ote_wfe=False):
        """Load an OPD produced by the JWST WSS into this instrument instance, specified by filename

        This includes:
            - If necessary, downloading that OPD from MAST. Downloaded files are cached in $WEBBPSF_PATH/MAST_JWST_WSS_OPDs
            - calling `import_wss_opd` to load the OPD from the FITS file and perform some necessary format conversions
            - Subtract off the instrument WFE for the field point used in wavefront sensing, to get an
                OTE-only wavefront. WebbPSF will separately add back in the SI WFE for the appropriate
                field point, as usual.
            - Subtract off the modeled field dependence term in the OTE WFE for the sensing field point, to get
                an estimate of the OTE wavefront nominally at the master chief ray location (between the NIRCams).
                WebbPSF will automatically add back on top of this the OTE field dependent WFE for the appropriate
                field point. as usual.
            - Scale the OPD to match the same size of the user provide pupil file

        Parameters
        ----------
        filename : str
            Name of OPD file to load

        output_path : str
            Downloaded OPD are saved in this location.
            This option is convinient for STScI users using /grp/jwst/ote/webbpsf-data/.
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
        # We use the size of the user supplied name of the JWST pupil in order to create the matching size OPD
        # The code assume the naming convention for the JWST pupil file: jwst_pupil_RevW_npix<size in pixels>.fits.gz
        npix_out = int(self.pupil[self.pupil.find('npix') + len('npix'):self.pupil.find('.fits')])

        if verbose and npix_out != 1024:
            print(
                  f'The size of the JWST pupil is different than nominal (1024px), {self.pupil}. '
                  f'The OPD will be scaled accordingly'
            )

        # If the provided filename doesn't exist on the local disk, try retrieving it from MAST
        # Note, this will automatically use cached versions downloaded previously, if present
        if not os.path.exists(filename):
            filename = webbpsf.mast_wss.mast_retrieve_opd(filename, output_path=output_path, verbose=verbose)

        if verbose:
            print(f'Importing and format-converting OPD from {filename}')
        opdhdu = webbpsf.mast_wss.import_wss_opd(filename, npix_out=npix_out)

        # Mask out any pixels in the OPD array which are outside the OTE pupil.
        # This is mostly cosmetic, and helps mask out some edge effects from the extrapolation + interpolation in
        # resizing the OPDs
        ote_pupil_mask = utils.get_pupil_mask(npix=npix_out) != 0
        opdhdu[0].data *= ote_pupil_mask

        # opdhdu[0].header['RMS_OBS'] = (webbpsf.utils.rms(opdhdu[0].data, mask=ote_pupil_mask)*1e9,
        #                               "[nm] RMS Observatory WFE (i.e. OTE+SI) at sensing field pt")

        if plot:
            import matplotlib
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(figsize=(16, 9), ncols=3, nrows=2)
            vm = 2e-7
            plot_kwargs = {'vmin': -vm, 'vmax': vm, 'cmap': matplotlib.cm.RdBu_r, 'origin': 'lower'}
            axes[0, 0].imshow(opdhdu[0].data.copy() * ote_pupil_mask, **plot_kwargs)
            axes[0, 0].set_title(f'OPD from\n{os.path.basename(filename)}')
            axes[0, 0].set_xlabel(f'RMS: {utils.rms(opdhdu[0].data*1e9, ote_pupil_mask):.2f} nm rms')

        if backout_si_wfe:
            # Check which field point was used for sensing
            sensing_apername = opdhdu[0].header['APERNAME']

            if verbose:
                print(f'Backing out SI WFE and OTE field dependence at the WF sensing field point ({sensing_apername})')

            # Create a temporary instance of an instrument, for the sensng instrument and field point,
            # in order to model and extract the SI WFE and OTE field dep WFE at the sensing field point.

            sensing_inst = instrument(sensing_apername[0:3])
            sensing_inst.pupil = (
                self.pupil
            )  # handle the case if the user has selected a different NPIX other than the default 1024
            sensing_inst.pupilopd = (
                opdhdu
            )  # handle the case if the user has selected a different NPIX other than the default 1024

            if sensing_inst.name == 'NRC':
                sensing_inst.filter = 'F212N'
                # TODO: optionally check for the edge case in which the sensing was done in F187N
                # note that there is a slight focus offset between the two wavelengths, due to NIRCam's refractive design
            # Set to the sensing aperture, and retrieve the OPD there
            sensing_inst.set_position_from_aperture_name(sensing_apername)
            # special case: for the main sensing points FP1 or FP6, we use the official WAS target phase map,
            # rather than the WebbPSF-internal SI WFE model.

            # Select correct target phase map based on sensing field point.
            # Note that the sensing maintenance program changed field point from NRC A3 to A1 around Dec 2024.
            if sensing_apername in ['NRCA3_FP1', 'NRCA1_FP6']:
                was_targ_file = utils.get_target_phase_map_filename(sensing_apername)
                sensing_fp_si_wfe = poppy.FITSOpticalElement(opd=was_targ_file).opd
            else:
                sensing_fp_si_wfe = sensing_inst.get_wfe('si')

            if npix_out != 1024:   # handle the case if the user has selected a different NPIX other than the default
                # the results from the zoom function preserve the STD between both phase maps and
                # the total sum between the phase maps is proportional to the zoom value
                sensing_fp_si_wfe = scipy.ndimage.zoom(sensing_fp_si_wfe, npix_out / 1024)

            sensing_fp_ote_wfe = sensing_inst.get_wfe('ote_field_dep')

            sihdu = fits.ImageHDU(sensing_fp_si_wfe)
            sihdu.header['EXTNAME'] = 'SENSING_SI_WFE'
            sihdu.header['CONTENTS'] = 'Model of SI WFE at sensing field point'
            sihdu.header['BUNIT'] = 'meter'
            sihdu.header['APERNAME'] = sensing_apername
            sihdu.header.add_history('This model for SI WFE was subtracted from the measured total WFE')
            sihdu.header.add_history('to estimate the OTE-only portion of the WFE.')
            opdhdu.append(sihdu)

            otehdu = fits.ImageHDU(sensing_fp_ote_wfe)
            otehdu.header['EXTNAME'] = 'SENSING_OTE_FD_WFE'
            otehdu.header['CONTENTS'] = 'Model of OTE field dependent WFE at sensing field point'
            otehdu.header['BUNIT'] = 'meter'
            otehdu.header['APERNAME'] = sensing_apername
            otehdu.header.add_history('This model for OTE field dependence was subtracted from the measured total WFE')
            otehdu.header.add_history('to estimate the OTE global portion of the WFE, at the master chief ray')
            opdhdu.append(otehdu)

            # Subtract the SI WFE from the WSS OPD, to obtain an estimated OTE-only OPD
            opdhdu[0].data -= (sensing_fp_si_wfe + sensing_fp_ote_wfe) * ote_pupil_mask
            opdhdu[0].header['CONTENTS'] = 'Estimated OTE WFE from Wavefront Sensing Measurements'
            opdhdu[0].header.add_history(f'Estimating SI WFE at sensing field point {sensing_apername}.')
            opdhdu[0].header.add_history('  See FITS extension SENSING_SI_WFE for the SI WFE model used.')
            opdhdu[0].header.add_history('  Subtracted SI WFE to estimate OTE-only global WFE.')
            opdhdu[0].header.add_history(f'Estimating OTE field dependence term at {sensing_apername}.')
            opdhdu[0].header.add_history(f'  Selected instrument field point is at V2,V3 = {sensing_inst._tel_coords()}.')
            opdhdu[0].header.add_history('  See FITS extension SENSING_OTE_FD_WFE for the WFE model used.')
            opdhdu[0].header.add_history('  Subtracted OTE field dependence to estimate OTE global WFE.')

            if plot or save_ote_wfe:
                # Either of these options will need the total OTE WFE.
                # Under normal circumstances webbpsf will compute this later automatically, but if needed we do it here too
                selected_fp_ote_wfe = sensing_inst.get_wfe('ote_field_dep')
                total_ote_wfe_at_fp = opdhdu[0].data + (selected_fp_ote_wfe * ote_pupil_mask)

            if plot:
                axes[0, 1].imshow(sensing_fp_si_wfe * ote_pupil_mask, **plot_kwargs)
                axes[0, 1].set_title(f'SI OPD\nat {sensing_apername}')
                axes[0, 1].set_xlabel(f'RMS: {utils.rms(sensing_fp_si_wfe * 1e9, ote_pupil_mask):.2f} nm rms')

                axes[0, 2].imshow(opdhdu[0].data + sensing_fp_ote_wfe * ote_pupil_mask, **plot_kwargs)
                axes[0, 2].set_title(f'OTE total OPD at sensing field point\ninferred from {os.path.basename(filename)}')
                axes[0, 2].set_xlabel(f'RMS: {utils.rms(opdhdu[0].data*1e9, ote_pupil_mask):.2f} nm rms')

                axes[1, 0].imshow(sensing_fp_ote_wfe * ote_pupil_mask, **plot_kwargs)
                axes[1, 0].set_title(f'OTE field dependent OPD\nat {sensing_apername}')
                axes[1, 0].set_xlabel(f'RMS: {utils.rms(sensing_fp_ote_wfe * 1e9, ote_pupil_mask):.2f} nm rms')

                axes[1, 1].imshow(selected_fp_ote_wfe * ote_pupil_mask, **plot_kwargs)
                axes[1, 1].set_title(f'OTE field dependent OPD\nat current field point in {self.name} {self.detector}')
                axes[1, 1].set_xlabel(f'RMS: {utils.rms(selected_fp_ote_wfe * 1e9, ote_pupil_mask):.2f} nm rms')

                axes[1, 2].imshow(total_ote_wfe_at_fp, **plot_kwargs)
                axes[1, 2].set_title(
                    f'Total OTE OPD at current FP in {self.name} {self.detector}\ninferred from {os.path.basename(filename)}'
                )
                axes[1, 2].set_xlabel(f'RMS: {utils.rms(total_ote_wfe_at_fp*1e9, ote_pupil_mask):.2f} nm rms')

                plt.tight_layout()

            if save_ote_wfe:
                # If requested, export the OPD for use in other external calculations.
                # We save out the total OTE WFE inferred at the selected instrument field point.
                outname = filename.replace('.fits', f'-ote-wfe-for-{self.name}-{self.detector}.fits')
                from copy import deepcopy

                opdhdu_at_si_fp = deepcopy(opdhdu)

                v2v3 = self._tel_coords()
                opdhdu_at_si_fp[0].header.add_history(
                    f'Estimating OTE field dependence term in {self.name} {self.detector}.'
                )
                opdhdu_at_si_fp[0].header.add_history(f'  Selected instrument field point is at V2,V3 = {v2v3}.')
                opdhdu_at_si_fp[0].header.add_history(
                    'Saving out total estimated OTE WFE (global+field dep) at that field point.'
                )
                opdhdu_at_si_fp[0].header['INSTRUME'] = self.name
                opdhdu_at_si_fp[0].header['DETECTOR'] = self.detector
                opdhdu_at_si_fp[0].header['APERNAME'] = self.aperturename
                opdhdu_at_si_fp[0].header['V2'] = self.aperturename

                # Save files with output units of microns, for consistency with other OPD files
                opdhdu_at_si_fp[0].data = total_ote_wfe_at_fp * 1e6
                opdhdu_at_si_fp[0].header['BUNIT'] = 'micron'

                opdhdu_at_si_fp.writeto(outname, overwrite=True)
                print(f'*****\nSaving estimated OTE-only WFE to file:\n\t{outname}\n*****')

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

    @poppy.utils.quantity_input(wavelengths=units.meter)
    def calc_datacube_fast(self, wavelengths, compare_methods=False, outfile=None,
                           add_distortion=True, *args, **kwargs):
        """Calculate a spectral datacube of PSFs: Simplified, much MUCH faster version.

        This is adapted from poppy.Instrument.calc_datacube, optimized and simplified
        for a substantial gain in speed at minimal reduction in accuracy for some use cases.

        ASSUMPTIONS:

        1) Assumes the wavefront error (OPD) and amplitude are independent of wavelength, such
            that we can do the expensive propagation from sky through the optics to the
            exit pupil of NIRSpec *only once*, save that, and reuse the same exit pupil wavefront
            many times changing only the wavelength for just the last DFT step to the detector.

        2) Assumes we do not need the binned-to-detector-resolution nor distorted versions;
            we just want the oversampled PSF datacube at many wavelengths as fast as possible.
            (If the binned output is also desired, it can be computed post facto.
            TODO: A future revision of this function may also add here an option for computing those
            derived versions as well.)

        Testing for NIRSpec IFU indicates this achieves ~150x speedup,
        and the differences in computed oversampled PSF are typically ~1/100th or less
        relative to the local PSF values in any given pixel.

        A consequence of the above assumption 1 is that this method is not well applicable
        for cases that have image plane masks, nor for NIRCam in general. It does seem to be
        reasonably applicable for NIRSpec IFU calculations within the current limited fidelity
        of webbpsf for that mode, IF we also neglect the image plane stop around the IFU FOV.

        Parameters
        ----------
        wavelengths : iterable of floats
            List or ndarray or tuple of floating point wavelengths in meters, such as
            you would supply in a call to calc_psf via the "monochromatic" option
        add_distortion : bool
            Same as for regular calc_psf.

        compare_methods : bool
            If true, compute the PSF **BOTH WAYS**, and return both for comparisons.
            This is of course much slower. Default is False. This is retained for
            test and debug usage for assessing cases in which this method is OK or not.


        Returns
        -------
        a PSF datacube, normally (with compare_methods=False)

        A list of two PSF datacubes and two exit wavefront objects, if compare_methods is True

        """

        nwavelengths = len(wavelengths)

        # Set up cube and initialize structure based on PSF at a representative wavelength
        _log.info('Starting fast/simplified multiwavelength data cube calculation.')
        ref_wave = np.mean(wavelengths)
        MIN_REF_WAVE = 2e-6 * units.meter  # This must not be too short, to avoid phase wrapping for the C3 bump
        if ref_wave < MIN_REF_WAVE:
            ref_wave = MIN_REF_WAVE
            log_message = (
                f'Performing initial propagation at minimum wavelength {MIN_REF_WAVE*1e6:.2f} microns; '
                'minimum set to avoid phase wrap of segment C3 surface.'
            )
            _log.info(log_message)
        else:
            _log.info(f'Performing initial propagation at average wavelength {ref_wave*1e6:.2f} microns.')

        psf, waves = self.calc_psf(*args, monochromatic=ref_wave, return_intermediates=True, **kwargs)
        from copy import deepcopy
        # Setup arrays to save data

        # Copy the first (oversampled) HDU only
        cubefast = astropy.io.fits.HDUList(deepcopy(psf[0]))
        try:
            # This is cosmetic only. Delete some exteraneous/redundant header keywords from poppy.
            # This function will below add a complete set of wavelength keywords
            del cubefast[0].header['WAVE0']
            del cubefast[0].header['WGHT0']
        except KeyError:
            pass

        ext = 0
        cubefast[ext].data = np.zeros((nwavelengths, psf[ext].data.shape[0], psf[ext].data.shape[1]))
        cubefast[ext].data[0] = psf[ext].data
        cubefast[ext].header[label_wavelength(nwavelengths, 0)] = wavelengths[0].to_value(units.meter)

        # Fast way. Assumes wavelength-independent phase and amplitude at the exit pupil!!
        if compare_methods:
            import time

            print('Running fast way')
            t0 = time.time()

        # Set up a simplified optical system just going from the exit pupil to the detector
        # Make the "entrance" pupil of this system replicate the exit pupl of the full calculation
        exitpupil = waves[-2]
        exit_opd = exitpupil.phase * exitpupil.wavelength.to_value(units.m) / (2 * np.pi)
        oversamp = psf[0].header['DET_SAMP']

        quickosys = poppy.OpticalSystem(
            npix=exitpupil.shape[0], pupil_diameter=exitpupil.shape[0] * units.pixel * exitpupil.pixelscale
        )
        quickosys.add_pupil(
            poppy.ArrayOpticalElement(opd=exit_opd, transmission=exitpupil.amplitude, pixelscale=exitpupil.pixelscale)
        )
        quickosys.add_detector(
            pixelscale=psf[0].header['PIXELSCL'] * oversamp,
            oversample=oversamp,
            fov_pixels=psf[0].header['NAXIS1'] // oversamp,
        )
        # Now do the propagations
        for i in range(0, nwavelengths):
            wl = wavelengths[i]
            psfw = quickosys.calc_psf(wavelength=wl, normalize='None')
            cubefast[0].data[i] = psfw[0].data
            cubefast[ext].header[label_wavelength(nwavelengths, i)] = wavelengths[i].to_value(units.meter)

        cubefast[0].header['NWAVES'] = nwavelengths

        # OPTIONAL
        # Also do the slower traditional way for comparison / debugging tests

        if compare_methods:
            psf2, waves2 = quickosys.calc_psf(wavelengths[0], return_intermediates=True)

            t1 = time.time()

            cube = deepcopy(psf)

            for ext in range(len(psf)):
                cube[ext].data = np.zeros((nwavelengths, psf[ext].data.shape[0], psf[ext].data.shape[1]))
                cube[ext].data[0] = psf[ext].data
                cube[ext].header[label_wavelength(nwavelengths, 0)] = wavelengths[0].to_value(units.meter)

            # iterate rest of wavelengths
            print('Running standard way')
            for i in range(0, nwavelengths):
                wl = wavelengths[i]
                psf = self.calc_psf(*args, monochromatic=wl, **kwargs)
                for ext in range(len(psf)):
                    cube[ext].data[i] = psf[ext].data
                    cube[ext].header[label_wavelength(nwavelengths, i)] = wl.to_value(units.meter)
                    cube[ext].header.add_history('--- Cube Plane {} ---'.format(i))
                    for h in psf[ext].header['HISTORY']:
                        cube[ext].header.add_history(h)
            t2 = time.time()
            cube[0].header['NWAVES'] = nwavelengths

            print(f'Fast way: {t1-t0:.3f} s')
            print(f'Standard way: {t2-t1:.3f} s')

            return cube, cubefast, waves, waves2  # return extra stuff for compariosns

        if outfile is not None:
            cubefast[0].header['FILENAME'] = (os.path.basename(outfile), 'Name of this file')
            cubefast.writeto(outfile, overwrite=True)
            _log.info('Saved result to ' + outfile)

        return cubefast


class JWInstrument_with_IFU(JWInstrument, ABC):
    """Subclass which adds some additional infrastructure for IFU sims"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # dict of modes and default aperture names
        self._modes_list = {'imaging': None, 'IFU': None}
        self._mode = 'imaging'

        self._IFU_bands_cubepars = {}  # placeholder, subclass should implement

    @property
    def mode(self):
        """Currently selected instrument major mode, imaging or IFU"""
        return self._mode

    @mode.setter
    def mode(self, value):
        if value not in self._modes_list:
            raise ValueError(f"'{value} is not an allowed mode for this instrument.")
        self._mode = value
        self.set_position_from_aperture_name(self._modes_list[value])

    @property
    @abstractmethod
    def band(self):
        # Subclass must implement this
        pass

    def get_IFU_wavelengths(self, nlambda=None):
        """Return an array of wavelengths spanning the currently selected IFU sub-band"""
        if self.mode != 'IFU':
            raise RuntimeError('This method only applies in IFU mode')
        spaxelsize, wavestep, minwave, maxwave = self._IFU_bands_cubepars[self.band]
        if nlambda:
            # Return the specified number of wavelengths, across that band
            return np.linspace(minwave, maxwave, nlambda) * units.micron
        else:
            # Return wavelength across that band using the same spectral sampling
            # as the instrument and pipeline
            return np.arange(minwave, maxwave, wavestep) * units.micron


class MIRI(JWInstrument_with_IFU):
    """A class modeling the optics of MIRI, the Mid-InfraRed Instrument.

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
        JWInstrument_with_IFU.__init__(self, 'MIRI')
        self.pixelscale = self._get_pixelscale_from_apername('MIRIM_FULL')
        self._rotation = 4.83544897  # V3IdlYAngle, Source: SIAF PRDOPSSOC-059
        # This is rotation counterclockwise; when summed with V3PA it will yield the Y axis PA on sky

        # Modes and default SIAF apertures for each
        self._modes_list = {
            'imaging': 'MIRIM_FULL',
            'IFU': 'MIRIFU_CHANNEL1A'
        }

        # Coordinate system note:
        # The pupil shifts get applied at the instrument pupil, which is an image of the OTE exit pupil
        # and is thus flipped in Y relative to the V frame entrance pupil. Therefore flip sign of pupil_shift_y
        self.options['pupil_shift_x'] = -0.0068  # In flight measurement. See Wright, Sabatke, Telfer 2022, Proc SPIE
        self.options['pupil_shift_y'] = -0.0110  # Sign intentionally flipped relative to that paper!! See note above.

        self.image_mask_list = ['FQPM1065', 'FQPM1140', 'FQPM1550', 'LYOT2300', 'LRS slit']
        self.pupil_mask_list = ['MASKFQPM', 'MASKLYOT', 'P750L']

        self._image_mask_apertures = {
            'FQPM1065': 'MIRIM_CORON1065',
            'FQPM1140': 'MIRIM_CORON1140',
            'FQPM1550': 'MIRIM_CORON1550',
            'LYOT2300': 'MIRIM_CORONLYOT',
        }
        self.auto_aperturename = True

        self.monochromatic = 8.0
        self._IFU_pixelscale = {
                  # slice width, pixel size.   Values from Argyriou et al. 2023 A&A 675
            'Ch1': (0.177, 0.196),
            'Ch2': (0.280, 0.196),
            'Ch3': (0.390, 0.245),
            'Ch4': (0.656, 0.273),
        }
        # The above tuples give the pixel resolution (first the 'alpha' direction, perpendicular to the slice,
        # then the 'beta' direction, along the slice).
        # The pixels are not square. See:
        # https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-medium-resolution-spectroscopy

        # Mappings between alternate names used for MRS subbands
        self._MRS_dichroic_to_subband = {'SHORT': 'A', 'MEDIUM': 'B', 'LONG': 'C'}
        self._MRS_subband_to_dichroic = {'A': 'SHORT', 'B': 'MEDIUM', 'C': 'LONG'}

        self._band = None
        #        self._MRS_bands = {"1A": [4.887326748103221, 5.753418963216559],   # Values provided by Polychronis Patapis
        #                           "1B": [5.644625711181792, 6.644794583147869],   # To-do: obtain from CRDS pipeline refs
        #                           "1C": [6.513777066360325, 7.669147994055998],
        #                           "2A": [7.494966046398437, 8.782517027772244],
        #                           "2B": [8.651469658142522, 10.168811217793243],
        #                           "2C": [9.995281242621394, 11.73039280033565],
        #                           "3A": [11.529088518317131, 13.491500288051483],
        #                           "3B": [13.272122736770127, 15.550153182343314],
        #                           "3C": [15.389530615108631, 18.04357852656418],
        #                           "4A": [17.686540162850203, 20.973301482912323],
        #                           "4B": [20.671069749545193, 24.476094964546686],
        #                           "4C": [24.19608171436692, 28.64871057821349]}
        self._IFU_bands_cubepars = {  # pipeline data cube parameters
            # Taken from ifucubepars_table in CRDS file 'jwst_miri_cubepar_0014.fits', current as of 2023 December
            # Each tuple gives pipeline spaxelsize, spectralstep, wave_min, wave_max
            '1A': (0.13, 0.0008, 4.90, 5.74),
            '1B': (0.13, 0.0008, 5.66, 6.63),
            '1C': (0.13, 0.0008, 6.53, 7.65),
            '2A': (0.17, 0.0013, 7.51, 8.77),
            '2B': (0.17, 0.0013, 8.67, 10.13),
            '2C': (0.17, 0.0013, 10.01, 11.70),
            '3A': (0.20, 0.0025, 11.55, 13.47),
            '3B': (0.20, 0.0025, 13.34, 15.57),
            '3C': (0.20, 0.0025, 15.41, 17.98),
            '4A': (0.35, 0.0060, 17.70, 20.95),
            '4B': (0.35, 0.0060, 20.69, 24.48),
            '4C': (0.35, 0.0060, 24.40, 28.70),
        }

        self._detectors = {
            'MIRIM': 'MIRIM_FULL',  # Mapping from user-facing detector names to SIAF entries.
            'MIRIFUSHORT': 'MIRIFU_CHANNEL1A',  # only applicable in IFU mode
            'MIRIFULONG': 'MIRIFU_CHANNEL3A',
        }  # ditto
        self.detector = 'MIRIM'
        self._detector_npixels = (1032, 1024)  # MIRI detector is not square
        self.detector_position = (512, 512)

        self._si_wfe_class = optics.MIRIFieldDependentAberrationAndObscuration

    def _get_default_fov(self):
        """Return default FOV in arcseconds"""
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
        """Validate instrument config for MIRI"""
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

        optsys.planes.pop(2)  # throw away the rotation of the entrance pupil we just added

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
        if (self.image_mask is not None and 'FQPM' in self.image_mask) or 'force_fqpm_shift' in self.options:
            optsys.add_pupil(poppy.FQPM_FFT_aligner())

        # Allow arbitrary offsets of the focal plane masks with respect to the pixel grid origin;
        # In most use cases it's better to offset the star away from the mask instead, using
        # options['source_offset_*'], but doing it this way instead is helpful when generating
        # the Pandeia ETC reference PSF library.
        offsets = {'shift_x': self.options.get('coron_shift_x', None), 'shift_y': self.options.get('coron_shift_y', None)}

        def make_fqpm_wrapper(name, wavelength):
            container = poppy.CompoundAnalyticOptic(
                name=name,
                opticslist=[
                    poppy.IdealFQPM(wavelength=wavelength, name=self.image_mask, **offsets),
                    poppy.SquareFieldStop(size=24, rotation=self._rotation, **offsets),
                ],
            )
            return container

        if self.image_mask == 'FQPM1065':
            optsys.add_image(make_fqpm_wrapper('MIRI FQPM 1065', 10.65e-6))
            trySAM = False
        elif self.image_mask == 'FQPM1140':
            optsys.add_image(make_fqpm_wrapper('MIRI FQPM 1140', 11.40e-6))
            trySAM = False
        elif self.image_mask == 'FQPM1550':
            optsys.add_image(make_fqpm_wrapper('MIRI FQPM 1550', 15.50e-6))
            trySAM = False
        elif self.image_mask == 'LYOT2300':
            # diameter is 4.25 (measured) 4.32 (spec) supposedly 6 lambda/D
            # optsys.add_image(function='CircularOcculter',radius =4.25/2, name=self.image_mask)
            # Add bar occulter: width = 0.722 arcsec (or perhaps 0.74, Dean says there is ambiguity)
            # optsys.add_image(function='BarOcculter', width=0.722, angle=(360-4.76))
            # position angle of strut mask is 355.5 degrees  (no = =360 -2.76 degrees
            # optsys.add_image(function='fieldstop',size=30)
            container = poppy.CompoundAnalyticOptic(
                name='MIRI Lyot Occulter',
                opticslist=[
                    poppy.CircularOcculter(radius=4.25 / 2, name=self.image_mask, **offsets),
                    poppy.BarOcculter(width=0.722, height=31, **offsets),
                    poppy.SquareFieldStop(size=30, rotation=self._rotation, **offsets),
                ],
            )
            optsys.add_image(container)
            trySAM = False  # FIXME was True - see https://github.com/mperrin/poppy/issues/169
            SAM_box_size = [5, 20]
        elif self.image_mask == 'LRS slit':
            # one slit, 5.5 x 0.6 arcsec in height (nominal)
            #           4.7 x 0.51 arcsec (measured for flight model. See MIRI-TR-00001-CEA)
            #
            # Per Klaus Pontoppidan: The LRS slit is aligned with the detector x-axis, so that the
            # dispersion direction is along the y-axis.
            # Slit width and height values derived from SIAF PRDOPSSOC-063, 2024 January
            # Undocumented options allow for offsetting the slit relative to the output pixel grid, to
            # more precisely match the actual instrument alignment
            lrs_slit = poppy.RectangularFieldStop(
                width=4.72345,
                height=0.51525,
                rotation=self._rotation,
                name=self.image_mask,
                shift_x=self.options.get('lrs_slit_offset_x', None),
                shift_y=self.options.get('lrs_slit_offset_y', None),
            )
            if self.options.get('lrs_use_mft', True):
                # Force the LRS slit to be rasterized onto a fine spatial sampling with gray subpixels
                # let's do a 3 arcsec box, sampled to 0.02 arcsec, with gray subpixels;
                # note poppy does not support non-square wavefront here
                lrs_pixscale = 0.02  # implicitly u.arcsec/u.pixel
                sampling = poppy.Wavefront(npix=int(5.5 / lrs_pixscale), pixelscale=lrs_pixscale)
                lrs_slit = poppy.fixed_sampling_optic(lrs_slit, sampling, oversample=8)

            optsys.add_image(optic=lrs_slit)
            trySAM = False
        else:
            optsys.add_image()
            trySAM = False

        if (self.image_mask is not None and 'FQPM' in self.image_mask) or 'force_fqpm_shift' in self.options:
            optsys.add_pupil(poppy.FQPM_FFT_aligner(direction='backward'))

        # add pupil plane mask
        shift_x, shift_y = self._get_pupil_shift()
        rotation = self.options.get('pupil_rotation', None)

        if self.options.get('coron_include_pre_lyot_plane', False) and self.pupil_mask.startswith('MASK'):
            optsys.add_pupil(poppy.ScalarTransmission(name='Pre Lyot Stop'))
            optsys.planes[3].wavefront_display_hint = 'intensity'

        if self.pupil_mask == 'MASKFQPM':
            optsys.add_pupil(
                transmission=self._datapath + '/optics/MIRI_FQPMLyotStop.fits.gz',
                name=self.pupil_mask,
                flip_y=True,
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
            )
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'MASKLYOT':
            optsys.add_pupil(
                transmission=self._datapath + '/optics/MIRI_LyotLyotStop.fits.gz',
                name=self.pupil_mask,
                flip_y=True,
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
            )
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'P750L' or self.image_mask == 'LRS slit':
            # This oversized pupil stop is present on all MIRI imaging filters, thus should
            # implicitly be included in all MIRI imager calculations, but in practice for
            # normal imaging modes, the system pupil stop is defined by the OTE primary, so this
            # stop has no effect. However for any light passing through the LRS slit, the spatial
            # filtering leads to diffractive spreading in the subsequen pupil which this should
            # be included for, in order to model slit losses correctly.
            optsys.add_pupil(
                transmission=self._datapath + '/optics/MIRI_LRS_Pupil_Stop.fits.gz',
                name=self.pupil_mask if self.pupil_mask else 'MIRI internal pupil stop',
                flip_y=True,
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
            )
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        else:  # all the MIRI filters have a tricontagon outline, even the non-coron ones.
            optsys.add_pupil(
                transmission=self._WebbPSF_basepath + '/tricontagon.fits.gz',
                name='filter cold stop',
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
            )
            # FIXME this is probably slightly oversized? Needs to have updated specifications here.

        optsys.add_rotation(-self._rotation, hide=True)
        optsys.planes[-1].wavefront_display_hint = 'intensity'

        if self.include_si_wfe:
            # now put back in the aberrations we grabbed above.
            # Note, the SI WFE models are in the detector coordinate frame, so this has to be added
            # *after* the rotation by ~5 degrees to that frame.
            optsys.add_pupil(miri_aberrations)

        # Special case for MIRI LRS slit spectroscopy. For this, we want to force the use of
        # MFT rather than FFT, for a small region, to ensure fine pixel sampling around the slit.
        # We can do this using poppy's MatrixFTCoronagraph class. No, the LRS is not a coronagraph;
        # but the desired handling of propagation steps and transforms is the same, so we can efficiently
        # reuse that existing poppy code path here.
        # The undocumented option for toggling this on/off mostly exists for testing as part of implementing
        # this enhancement, and could be removed from the code later.
        if self.image_mask == 'LRS slit' and self.options.get('lrs_use_mft', True):
            _log.info('Setting up special propagator for Matrix DFTs around MIRI LRS slit')

            # hard-coded values here are for a box encompassing the LRS slit, and sampled
            # sufficiently finely (8x Nyquist) to yield relatively precise and accurate results
            optsys = poppy.MatrixFTCoronagraph(optsys, occulter_box=[1, 3], oversample=8)
            trySAM = False

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
        if self._image_mask is not None:
            if 'LRS' in self._image_mask:
                apname = 'MIRIM_FULL'  # LRS slit uses full array readout
            else:
                apname = self._image_mask_apertures[self._image_mask]
        elif self.mode == 'imaging':
            apname = 'MIRIM_FULL'
        else:
            # IFU mode is complex, don't try to set a different apname here
            # (This gracefully works around an edge case in mode switching)
            apname = self.aperturename

        # Call aperturename.setter to update ap ref coords and DetectorGeometry class
        self.aperturename = apname

        str_debug = '_update_aperturename AFTER  - Det: {}, Ap: {}, ImMask: {}, PupMask: {}, DetPos: {}'.format(
            self._detector, self._aperturename, self.image_mask, self.pupil_mask, self.detector_position
        )
        _log.debug(str_debug)

    def _get_fits_header(self, hdulist, options):
        """Format MIRI-like FITS headers, based on JWST DMS SRD 1 FITS keyword info"""
        super(MIRI, self)._get_fits_header(hdulist, options)

        hdulist[0].header['GRATNG14'] = ('None', 'MRS Grating for channels 1 and 4')
        hdulist[0].header['GRATNG23'] = ('None', 'MRS Grating for channels 2 and 3')
        hdulist[0].header['FLATTYPE'] = ('?', 'Type of flat field to be used: all, one, principal')
        hdulist[0].header['CCCSTATE'] = ('open', 'Contamination Control Cover state: open, closed, locked')
        if self.image_mask is not None:
            hdulist[0].header['TACQNAME'] = ('None', 'Target acquisition file name')

    def _get_pixelscale_from_apername(self, apername):
        """Simple utility function to look up pixelscale from apername"""

        if 'MIRIFU' in apername:
            if apername.startswith('MIRIFU_CHANNEL'):
                band = apername[-2:]
                spaxelsize, _, _, _ = self._IFU_bands_cubepars[band]
                return spaxelsize
            else:
                raise RuntimeError(f'Not sure how to determine pixelscale for {apername}')
        else:
            return super()._get_pixelscale_from_apername(apername)

    def _get_aperture_rotation(self, apername):
        """Get the rotation angle of a given aperture, using values from SIAF.

        Returns ~ position angle counterclockwise from the V3 axis, in degrees
        (i.e. consistent with SIAF V3IdlYangle)

        Note, MIRIFU aperture geometry is extremely complex, and this oversimplifies.
        See https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-instrumentation/miri-mrs-field-and-coordinates
        Consistent with that reference, we compute the angle of the along-slice (alpha) direction relative to the
        horizontal (V2). Because the apertuers are skewed, this yields different values by several degrees
        than an equivalent calculation for beta.
        """
        if apername.startswith('MIRIM'):
            return self.siaf['MIRIM_FULL'].V3IdlYAngle
        elif apername.startswith('MIRIFU'):
            # These SLIT or COMPOUND apertures do not have a V3IdlYangle param defined in SIAF
            # But we can work out the angles from the aperture corner vertices which are defined.
            cx, cy = self.siaf[apername].corners('tel', rederive=False)
            # The aperture shapes are irregular quadrilaterals, not squares or rectangles
            # So, take the angles of both alpha-axis (~horizontal) sides, and average them to get an average rotation angle
            dx = cx[0] - cx[1]
            dy = cy[0] - cy[1]
            dx2 = cx[3] - cx[2]
            dy2 = cy[3] - cy[2]
            # take the average, and convert to degrees.
            # The signs and parity of the arctan2 are atypical here, to match the expected output convention of
            # rotation counterclockwise relative to the V2V3 frame.
            avg_V3IdlYangle = np.rad2deg((np.arctan2(dy, -dx) + np.arctan2(dy2, -dx2)) / 2)
            return avg_V3IdlYangle
        else:
            raise ValueError(f'Unexpected/invalid apername for MIRI: {apername}')

    @JWInstrument_with_IFU.aperturename.setter
    def aperturename(self, value):
        """Set aperturename, also update the rotation for MIRIM vs. IFU channel"""
        # apply the provided aperture name
        # Note, the syntax for calling a parent class property setter is... painful:
        super(MIRI, type(self)).aperturename.fset(self, value)
        # Update the rotation angle
        self._rotation = self._get_aperture_rotation(self.aperturename)

        # if it's an IFU aperture, we're now in IFU mode:
        self._mode = 'IFU' if value.startswith('MIRIFU') else 'imaging'
        # if in IFU mode, we probably also want to update the IFU band

        if self.mode == 'IFU':
            if self.band != value[-2:]:
                self.band = value[-2:]
            self._detector = 'MIRIFULONG' if self.band[0] in ['3', '4'] else 'MIRIFUSHORT'

    @property
    def band(self):
        """MRS IFU spectral band. E.g. '1A', '3B'. Only applicable in IFU mode."""
        if self.mode == 'IFU':
            return self._band
        else:
            return None

    @band.setter
    def band(self, value):
        if self.mode != 'IFU':
            if value is not None:
                raise RuntimeError("The 'band' property is only valid for IFU mode simulations.")
            return

        if value in self._IFU_bands_cubepars.keys():
            self._band = value
            self._ifu_slice_width = self._IFU_pixelscale[f"Ch{self._band[0]}"][0]
            self.aperturename = 'MIRIFU_CHANNEL' + value
            # setting aperturename will also auto update self._rotation
            # self._rotation = self.MRS_rotation[self._band]
            # update filter, image_mask and detector
            # self._filter = "D"+ self.subband_to_dichroic[self._band[1]]
            # self._image_mask = "MIRI-IFU_" + self._band[0]
            # self._update_detector()
        # if not (self.MRSbands[self.band][0] <= self._wavelength <= self.MRSbands[self.band][1]):
        #    self._wavelength = np.mean(self.MRSbands[self.band])
        else:
            raise ValueError(f'Not a valid MRS band: {value}')

    def _calc_psf_format_output(self, result, options):
        """Format output HDUList. In particular, add some extra metadata if in IFU mode"""
        super()._calc_psf_format_output(result, options)
        if self.mode == 'IFU':
            n_exts = len(result)
            for ext in np.arange(n_exts):
                result[ext].header['MODE'] = ('IFU', 'This is a MIRI MRS IFU mode simulation')
                result[ext].header['FILTER'] = ('MIRIFU_CHANNEL' + self.band, 'MIRI IFU sub-band simulated')
                result[ext].header['BAND'] = (self.band, 'MIRI IFU sub-band simulated')


class NIRCam(JWInstrument):
    """A class modeling the optics of NIRCam.

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
    nc.options['bar_offset'] = 3  # 3 arcseconds towards the right (narrow end on module A)
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
        # need to set up a bunch of stuff here before calling superclass __init__
        # so the overridden filter setter will not have errors when called from __init__
        self.auto_channel = False
        self.auto_aperturename = False
        JWInstrument.__init__(self, 'NIRCam')

        self._pixelscale_short = self._get_pixelscale_from_apername('NRCA1_FULL')
        self._pixelscale_long = self._get_pixelscale_from_apername('NRCA5_FULL')
        self.pixelscale = self._pixelscale_short

        self.options['pupil_shift_x'] = 0  # Set to 0 since NIRCam FAM corrects for PM shear in flight
        self.options['pupil_shift_y'] = 0

        # Enable the auto behaviours by default (after superclass __init__)
        self.auto_channel = True
        self.auto_aperturename = True
        self._filter = 'F200W'
        self._detector = 'NRCA1'

        self.image_mask_list = ['MASKLWB', 'MASKSWB', 'MASK210R', 'MASK335R', 'MASK430R']
        self._image_mask_apertures = {
            'MASKLWB': 'NRCA5_MASKLWB',
            'MASKSWB': 'NRCA4_MASKSWB',
            'MASK210R': 'NRCA2_MASK210R',
            'MASK335R': 'NRCA5_MASK335R',
            'MASK430R': 'NRCA5_MASK430R',
        }

        self.pupil_mask_list = [
            'CIRCLYOT',
            'WEDGELYOT',
            'MASKRND',
            'MASKSWB',
            'MASKLWB',
            # The last 3 of the above are synonyms for the first 2
            'WEAK LENS +4',
            'WEAK LENS +8',
            'WEAK LENS -8',
            'WEAK LENS +12 (=4+8)',
            'WEAK LENS -4 (=4-8)',
            'WLP4',
            'WLM4',
            'WLP8',
            'WLM8',
            'WLP12',
        ] + [f'DHS_{i+1:02d}' for i in range(10)]

        self._detectors = dict()
        det_list = ['A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5']
        for name in det_list:
            self._detectors['NRC{0}'.format(name)] = 'NRC{0}_FULL'.format(name)
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
        if self._image_mask is not None:
            aps_modA = {
                'MASKLWB': 'NRCA5_FULL_MASKLWB',
                'MASKSWB': 'NRCA4_FULL_MASKSWB',
                'MASK210R': 'NRCA2_FULL_MASK210R',
                'MASK335R': 'NRCA5_FULL_MASK335R',
                'MASK430R': 'NRCA5_FULL_MASK430R',
            }
            # Choose coronagraphic subarray apertures for Module B
            aps_modB = {
                'MASKLWB': 'NRCB5_MASKLWB',
                'MASKSWB': 'NRCB3_MASKSWB',
                'MASK210R': 'NRCB1_MASK210R',
                'MASK335R': 'NRCB5_MASK335R',
                'MASK430R': 'NRCB5_MASK430R',
            }
            apname = aps_modA[self._image_mask] if self.module == 'A' else aps_modB[self._image_mask]
            _log.debug(f'Inferred {apname} from coronagraph focal plane mask selected.')
        elif (self._pupil_mask is not None) and (('LYOT' in self._pupil_mask) or ('MASK' in self._pupil_mask)):
            # Want to use full frame apertures if only Lyot stops defined (no image mask)
            # Unfortunately, no full frame SIAF apertures are defined for Module B w/ Lyot
            # so we must select the subarray apertures as a special case.
            if 'long' in self.channel:
                if ('WEDGE' in self._pupil_mask) or ('LWB' in self._pupil_mask):
                    apname = 'NRCA5_FULL_WEDGE_BAR' if self.module == 'A' else 'NRCB5_MASKLWB'
                else:
                    apname = 'NRCA5_FULL_WEDGE_RND' if self.module == 'A' else 'NRCB5_MASK335R'
            else:
                if ('WEDGE' in self._pupil_mask) or ('SWB' in self._pupil_mask):
                    apname = 'NRCA4_FULL_WEDGE_BAR' if self.module == 'A' else 'NRCB3_MASKSWB'
                else:
                    apname = 'NRCA2_FULL_WEDGE_RND' if self.module == 'A' else 'NRCB1_MASK210R'
                    _log.debug(
                        f'Inferred {apname} from coronagraph Lyot mask selected,',
                        f'and channel={self.channel}, module={self.module}'
                    )
        else:
            apname = self._detectors[self._detector]
            _log.debug(f'Inferred {apname} from selected detector.')

        # Call aperturename.setter to update ap ref coords and DetectorGeometry class
        self.aperturename = apname

        str_debug = '_update_aperturename AFTER  - Det: {}, Ap: {}, ImMask: {}, PupMask: {}, DetPos: {}'.format(
            self._detector, self._aperturename, self.image_mask, self.pupil_mask, self.detector_position
        )
        _log.debug(str_debug)

    @JWInstrument.aperturename.setter
    def aperturename(self, value):
        """Set SIAF aperture name to new value, with validation.

        This also updates the pixelscale to the local value for that aperture, for a small precision enhancement.
        """
        # Explicitly update detector reference coordinates,
        # otherwise old coordinates can persist under certain circumstances

        # Get NIRCam SIAF apertures
        try:
            ap = self.siaf[value]
        except KeyError:
            _log.warning(f'Aperture name {value} not a valid NIRCam pysiaf name')
            # Alternatives in case we are running an old pysiaf PRD
            if value == 'NRCA5_FULL_WEDGE_BAR':
                newval = 'NRCA5_FULL_MASKLWB'
            elif value == 'NRCA5_FULL_WEDGE_RND':
                newval = 'NRCA5_FULL_MASK335R'
            elif value == 'NRCA4_FULL_WEDGE_BAR':
                newval = 'NRCA4_FULL_MASKSWB'
            elif value == 'NRCA2_FULL_WEDGE_RND':
                newval = 'NRCA2_FULL_MASK210R'
            else:
                newval = None

            if newval is not None:
                # Set alternative aperture name as bandaid to continue
                value = newval
                warning_message = (
                    'Possibly running an old version of pysiaf missing some NIRCam apertures. '
                    'Continuing with old aperture names.'
                )
                _log.warning(warning_message)
            else:
                return

        # Only update if new value is different
        if self._aperturename != value:
            # First, check some info from current settings, wich we will use below as part of auto pixelscale code
            # The point is to check if the pixel scale is set to a custom or default value,
            # and if it's custom then don't override that.
            # Note, check self._aperturename first to account for the edge case when
            # this is called from __init__ before _aperturename is set
            has_custom_pixelscale = self._aperturename and (
                self.pixelscale != self._get_pixelscale_from_apername(self._aperturename)
            )

            # Now apply changes:
            self._aperturename = value
            # Update detector reference coordinates
            self.detector_position = (ap.XSciRef, ap.YSciRef)

            # Check if detector is correct
            new_det = self._aperturename[0:5]
            if new_det != self._detector:
                new_channel = 'long' if new_det[-1] == '5' else 'short'
                self._switch_channel(new_channel)
                self._detector = new_det

            # Update DetectorGeometry class
            self._detector_geom_info = DetectorGeometry(self.siaf, self._aperturename)
            _log.info('NIRCam aperture name updated to {}'.format(self._aperturename))

            if not has_custom_pixelscale:
                self.pixelscale = self._get_pixelscale_from_apername(self._aperturename)
                debug_message = (
                    f'Pixelscale updated to {self.pixelscale} '
                    f'based on average X+Y SciScale at SIAF aperture {self._aperturename}'
                )
                _log.debug(debug_message)

    @property
    def module(self):
        return self._detector[3]
        # note, you can't set module directly; it's inferred based on detector.

    @module.setter
    def module(self, value):
        raise RuntimeError('NIRCam module is not directly settable; set detector instead.')

    @property
    def channel(self):
        return 'long' if self.detector.endswith('5') else 'short'
        # note, you can't set channel directly; it's inferred based on detector.

    @channel.setter
    def channel(self, value):
        raise RuntimeError('NIRCam channel is not directly settable; set filter or detector instead.')

    @JWInstrument.detector.setter  # override setter in this subclass, to implement auto channel switch
    def detector(self, value):
        """Set detector, including reloading the relevant info from SIAF"""
        if value.upper().endswith('LONG'):
            # treat NRCALONG and NRCBLONG as synonyms to NRCA5 and NRCB5
            value = value[:-4] + '5'
        if value.upper() not in self.detector_list:
            raise ValueError('Invalid detector. Valid detector names are: {}'.format(', '.join(self.detector_list)))
        # set the channel based on the requested detector
        new_channel = 'long' if value[-1] == '5' else 'short'
        self._switch_channel(new_channel)
        self._detector = value.upper()
        self._update_aperturename()

    def _switch_channel(self, channel):
        """Toggle to either SW or LW channel.
        This changes the detector name and the pixel scale,
        unless the user has set a custom/nonstandard pixel scale manually.
        """
        if self.channel == channel:
            return  # nothing to do
        _log.debug('Automatically changing NIRCam channel SW/LW to ' + channel)
        if channel == 'long':
            # ensure long wave by switching to detector 5
            self._detector = self._detector[0:4] + '5'
            if self.pixelscale == self._pixelscale_short:
                self.pixelscale = self._pixelscale_long
                _log.info('NIRCam pixel scale switched to %f arcsec/pixel for the ' 'long wave channel.' % self.pixelscale)
        elif channel == 'short':
            # only change detector if the detector was already LW;
            # don't override selection of a particular SW SCA otherwise
            if self._detector[-1] == '5':
                self._detector = self._detector[0:4] + '1'
            if self.pixelscale == self._pixelscale_long:
                self.pixelscale = self._pixelscale_short
                _log.info('NIRCam pixel scale switched to %f arcsec/pixel for the ' 'short wave channel.' % self.pixelscale)
        else:
            raise ValueError('Invalid NIRCam channel name: {}'.format(channel))

    @JWInstrument.filter.setter
    def filter(self, value):
        super(NIRCam, self.__class__).filter.__set__(self, value)

        if self.auto_channel or self.auto_aperturename:
            # set the channel (via setting the detector) based on filter
            if self.filter == 'WLP4':
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
        if name == '':
            name = None
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
            _log.info(f'NIRCam pupil mask setter: aperturename {self._aperturename}')

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
            raise RuntimeError('The requested wavelengths are too short to be imaged with NIRCam')
        if np.max(wavelengths) > self.LONG_WAVELENGTH_MAX:
            raise RuntimeError('The requested wavelengths are too long to be imaged with NIRCam')
        if self.channel == 'short' and np.max(wavelengths) > self.SHORT_WAVELENGTH_MAX:
            raise RuntimeError('The requested wavelengths are too long for NIRCam short wave channel.')
        if self.channel == 'long' and np.min(wavelengths) < self.LONG_WAVELENGTH_MIN:
            raise RuntimeError('The requested wavelengths are too short for NIRCam long wave channel.')

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
        shifts = {'shift_x': self.options.get('coron_shift_x', None), 'shift_y': self.options.get('coron_shift_y', None)}

        if (self.image_mask == 'MASK210R') or (self.image_mask == 'MASK335R') or (self.image_mask == 'MASK430R'):
            optsys.add_image(
                NIRCam_BandLimitedCoron(name=self.image_mask, module=self.module, nd_squares=nd_squares, **shifts), index=2
            )
            trySAM = False  # FIXME was True - see https://github.com/mperrin/poppy/issues/169
            SAM_box_size = 5.0
        elif (self.image_mask == 'MASKSWB') or (self.image_mask == 'MASKLWB'):
            bar_offset = self.options.get('bar_offset', None)
            # If the bar offset is not provided, use the SIAF aperture name, or else the filter name to lookup the default
            # position. If an offset is provided and is a floating point value, use that
            # directly as the offset. Otherwise assume it's a filter name and try passing
            # that in to the auto offset. (that allows for selecting the narrow position, or
            # for simulating using a given filter at some other filter's position.)
            # This code is somewhat convoluted, for historical reasons and back-compatibility
            if bar_offset is None:
                # Try to use the SIAF aperture name to determine the offset
                # This can help better automate simulations matching data, since match_data.py will
                # have copied the aperturename from the header, like NRCA5_MASKLWB_NARROW or similar
                if 'MASK' in self.aperturename:
                    apname_last_part = self.aperturename.split('_')[-1]
                    if apname_last_part == 'NARROW':
                        auto_offset = 'narrow'  # set to lower case for consistency with existing code in optics.py
                        _log.info(f'Set bar offset to {auto_offset} based on current aperture name {self.aperturename}')
                    elif apname_last_part.startswith('F'):
                        auto_offset = apname_last_part
                        _log.info(f'Set bar offset to {auto_offset} based on current aperture name {self.aperturename}')
                    else:
                        auto_offset = self.filter
                        _log.info(f'Set bar offset to {auto_offset} based on current filter {self.filter}')
                else:
                    auto_offset = self.filter
                    _log.info(f'Set bar offset to {auto_offset} based on current filter {self.filter}')

            else:
                try:
                    _ = float(bar_offset)
                    auto_offset = None
                except ValueError:
                    # If the "bar_offset" isn't a float, pass it to auto_offset instead
                    auto_offset = bar_offset
                    bar_offset = None

            optsys.add_image(
                NIRCam_BandLimitedCoron(
                    name=self.image_mask,
                    module=self.module,
                    nd_squares=nd_squares,
                    bar_offset=bar_offset,
                    auto_offset=auto_offset,
                    **shifts,
                ),
                index=2,
            )
            trySAM = False  # True FIXME
            SAM_box_size = [5, 20]
        elif (
            (self.pupil_mask is not None)
            and ('LENS' not in self.pupil_mask.upper())
            and ('WL' not in self.pupil_mask.upper())
            and ('DHS' not in self.pupil_mask.upper())
        ):
            # no occulter selected but coronagraphic mode anyway. E.g. off-axis PSF
            # but don't add this image plane for weak lens or DHS calculations
            optsys.add_image(poppy.ScalarTransmission(name='No Image Mask Selected!'), index=2)
            trySAM = False
        else:
            trySAM = False

        # add pupil plane mask
        shift_x, shift_y = self._get_pupil_shift()
        rotation = self.options.get('pupil_rotation', None)

        if self.pupil_mask == 'CIRCLYOT' or self.pupil_mask == 'MASKRND':
            optsys.add_pupil(
                transmission=self._datapath + '/optics/NIRCam_Lyot_Somb.fits.gz',
                name=self.pupil_mask,
                flip_y=True,
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
                index=3,
            )
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'WEDGELYOT' or self.pupil_mask == 'MASKSWB' or self.pupil_mask == 'MASKLWB':
            optsys.add_pupil(
                transmission=self._datapath + '/optics/NIRCam_Lyot_Sinc.fits.gz',
                name=self.pupil_mask,
                flip_y=True,
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
                index=3,
            )
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        # Note, for historical reasons there are multiple synonymous ways to specify the weak lenses
        # This includes versions that elide over the fact that WLP4 is in the filter wheel, plus
        # versions that take that into account explicitly.
        elif (
            self.pupil_mask == 'WEAK LENS +4'
            or self.pupil_mask == 'WLP4'
            or (self.filter == 'WLP4' and self.pupil_mask is None)
        ):
            optsys.add_pupil(
                optics.NIRCamFieldDependentWeakLens(
                    name='WLP4',
                    instrument=self,
                    shift_x=shift_x,
                    shift_y=shift_y,
                    rotation=rotation,
                ),
                index=3,
            )
        elif self.pupil_mask == 'WEAK LENS +8' or (self.pupil_mask == 'WLP8' and self.filter != 'WLP4'):
            optsys.add_pupil(
                optics.NIRCamFieldDependentWeakLens(
                    name='WLP8',
                    instrument=self,
                    shift_x=shift_x,
                    shift_y=shift_y,
                    rotation=rotation,
                ),
                index=3,
            )
        elif self.pupil_mask == 'WEAK LENS -8' or (self.pupil_mask == 'WLM8' and self.filter != 'WLP4'):
            optsys.add_pupil(
                optics.NIRCamFieldDependentWeakLens(
                    name='WLM8',
                    instrument=self,
                    shift_x=shift_x,
                    shift_y=shift_y,
                    rotation=rotation,
                ),
                index=3,
            )
        elif (
            self.pupil_mask == 'WEAK LENS +12 (=4+8)'
            or self.pupil_mask == 'WLP12'
            or (self.pupil_mask == 'WLP8' and self.filter == 'WLP4')
        ):
            optsys.add_pupil(
                optics.NIRCamFieldDependentWeakLens(
                    name='WLP12',
                    instrument=self,
                    shift_x=shift_x,
                    shift_y=shift_y,
                    rotation=rotation,
                ),
                index=3,
            )
        elif (
            self.pupil_mask == 'WEAK LENS -4 (=4-8)'
            or self.pupil_mask == 'WLM4'
            or (self.pupil_mask == 'WLM8' and self.filter == 'WLP4')
        ):
            optsys.add_pupil(
                optics.NIRCamFieldDependentWeakLens(
                    name='WLM4',
                    instrument=self,
                    shift_x=shift_x,
                    shift_y=shift_y,
                    rotation=rotation,
                ),
                index=3,
            )

        elif self.pupil_mask is None and self.image_mask is not None:
            optsys.add_pupil(poppy.ScalarTransmission(name='No Lyot Mask Selected!'), index=3)
        elif self.pupil_mask.startswith('DHS'):
            optsys.add_pupil(
                transmission=self._datapath + f'/optics/NIRCam_{self.pupil_mask}_npix1024.fits.gz',
                name=self.pupil_mask,
                flip_y=True,
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
                index=3,
            )
            optsys.planes[3].wavefront_display_hint = 'intensity'

        else:
            optsys.add_pupil(
                transmission=self._WebbPSF_basepath + '/tricontagon_oversized_4pct.fits.gz',
                name='filter stop',
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
            )

        if self.options.get('coron_include_pre_lyot_plane', False) and self.pupil_mask.startswith('MASK'):
            optsys.add_pupil(
                poppy.ScalarTransmission(name='Pre Lyot Stop'), index=3
            )  # this is before the above plane, but do the insertion here
            # because of all the hard-coded index=3 above

            optsys.planes[3].wavefront_display_hint = 'intensity'

        return (optsys, trySAM, SAM_box_size)

    def _get_fits_header(self, hdulist, options):
        """Format NIRCam-like FITS headers, based on JWST DMS SRD 1 FITS keyword info"""
        super(NIRCam, self)._get_fits_header(hdulist, options)

        hdulist[0].header['MODULE'] = (self.module, 'NIRCam module: A or B')
        hdulist[0].header['CHANNEL'] = ('Short' if self.channel == 'short' else 'Long', 'NIRCam channel: long or short')
        # filter, pupil added by calc_psf header code
        hdulist[0].header['PILIN'] = ('False', 'Pupil imaging lens in optical path: T/F')


class NIRSpec(JWInstrument_with_IFU):
    """A class modeling the optics of NIRSpec, in **imaging** mode.

    This is not a substitute for a spectrograph model, but rather a way of simulating a PSF as it
    would appear with NIRSpec in imaging mode (e.g. for target acquisition).  NIRSpec support is
    relatively simplistic compared to the other instruments at this point.

    Relevant attributes include `filter`. In addition to the actual filters, you may select 'IFU' to
    indicate use of the NIRSpec IFU, in which case use the `monochromatic` attribute to set the simulated wavelength.

    If a grating is selected in the pupil, then a rectangular pupil mask 8.41x7.91 m as projected onto the primary
    is added to the optical system. This is an estimate of the pupil stop imposed by the outer edge of the grating
    clear aperture, estimated based on optical modeling by Erin Elliot and Marshall Perrin.

    Notes on IFU support:
        Additional features for modeling NRS IFU PSFs are enabled by setting the .mode attribute to 'IFU'.

        The pipeline-output data products, assuming the 'ifualign' frame is used in the cube build step, which
        is rotated relative to the typical 'sci' output frame used in all other webbpsf sim outputs.
        For convenience, for IFU-mode simulations an extra rotation is included in the PSF calculation
        such that the output product orientation matches the IFUalign s3d cube orientation. This happens
        automatically and transparently to the user (and source offset parameters, e.g. options['source_offset_x']
        will automatically be interpreted as X and Y position in that output frame, including the effects of the
        rotation). If the rotation to ifualign frame is for some reason not desired, it can be disabled by
        setting nrs.options['ifualign_rotation'] = False

    """

    def __init__(self):
        JWInstrument_with_IFU.__init__(self, 'NIRSpec')
        self.pixelscale = 0.10435  # Average over both detectors.  SIAF PRDOPSSOC-059, 2022 Dec
        # Microshutters are 0.2x0.46 but we ignore that here.
        self._rotation = 138.5  # Average for both detectors in SIAF PRDOPSSOC-059
        # This is rotation counterclockwise; when summed with V3PA it will yield the Y axis PA on sky

        # Modes and default SIAF apertures for each
        self._modes_list = {
            'imaging': 'NRS1_FULL',
            'IFU': 'NRS_FULL_IFU'
        }

        self._IFU_pixelscale = 0.1043  # same.
        self.monochromatic = 3.0
        self.filter = 'F110W'  # or is this called F115W to match NIRCam??

        self.options['pupil_shift_x'] = 0.0115  # CV3 on-orbit estimate (RPT028027) + OTIS delta from predicted (037134)
        self.options['pupil_shift_y'] = -0.0157

        # fixed slits
        self.image_mask_list = [
            'S200A1',
            'S200A2',
            'S400A1',
            'S1600A1',
            'S200B1',
            'MSA all open',
            'Single MSA open shutter',
            'Three adjacent MSA open shutters',
            'IFU',
        ]
        self.pupil_mask_list = ['NIRSpec grating']
        self.image_mask = 'MSA all open'
        self.pupil_mask = self.pupil_mask_list[-1]

        self.disperser_list = ['PRISM', 'G140M', 'G140H', 'G235M', 'G235H', 'G395M', 'G395H']
        self._disperser = None
        self._IFU_bands_cubepars = {
            'PRISM/CLEAR': (0.10, 0.0050, 0.60, 5.30),
            'G140M/F070LP': (0.10, 0.0006, 0.70, 1.27),
            'G140M/F100LP': (0.10, 0.0006, 0.97, 1.89),
            'G140H/F070LP': (0.10, 0.0002, 0.70, 1.27),
            'G140H/F100LP': (0.10, 0.0002, 0.97, 1.89),
            'G235M/F170LP': (0.10, 0.0011, 1.66, 3.17),
            'G235H/F170LP': (0.10, 0.0004, 1.66, 3.17),
            'G395M/F290LP': (0.10, 0.0018, 2.87, 5.27),
            'G395H/F290LP': (0.10, 0.0007, 2.87, 5.27),
        }

        det_list = ['NRS1', 'NRS2']
        self._detectors = dict()
        for name in det_list:
            self._detectors[name] = '{0}_FULL'.format(name)
        self.detector = self.detector_list[0]
        self.detector_position = (1380, 1024)  # near S1600A1 square aperture / ISIM1 field point. see #348.
        self._si_wfe_class = optics.NIRSpecFieldDependentAberration  # note we end up adding 2 instances of this.

    def _validate_config(self, **kwargs):
        return super(NIRSpec, self)._validate_config(**kwargs)

    def _addAdditionalOptics(self, optsys, oversample=2):
        """Add fixed slit optics for NIRSpec

        See Table 3-6 of NIRSpec Ops Concept Document, ESA-JWST-TN-0297 / JWST-OPS-003212

        """
        from .optics import NIRSpec_MSA_open_grid, NIRSpec_three_MSA_shutters

        trySAM = False  # semi-analytic method never applicable here.
        SAM_box_size = None

        if self.image_mask == 'S200A1' or self.image_mask == 'S200A2' or self.image_mask == 'S200B1':
            # three identical slits, 0.2 x 3.2 arcsec in length
            optsys.add_image(optic=poppy.RectangularFieldStop(width=0.2, height=3.2, name=self.image_mask + ' slit'))
        elif self.image_mask == 'S400A1':
            # one slit, 0.4 x 3.65 arcsec in height
            optsys.add_image(optic=poppy.RectangularFieldStop(width=0.4, height=3.65, name=self.image_mask + ' slit'))
        elif self.image_mask == 'S1600A1':
            # square aperture for exoplanet spectroscopy
            optsys.add_image(
                optic=poppy.RectangularFieldStop(width=1.6, height=1.6, name=self.image_mask + ' square aperture')
            )
        elif self.image_mask == 'IFU':
            # square aperture for the entrance to the slicer.
            # DOES NOT ACTUALLY MODEL THE SLICER OPTICS AT ALL!
            # Values talen from pre-flight SIAF, fall 2017
            optsys.add_image(optic=poppy.RectangularFieldStop(width=3.193, height=3.097, name='IFU entrance'))
        elif self.image_mask == 'MSA all open':
            # all MSA shutters open
            optsys.add_image(optic=NIRSpec_MSA_open_grid(name=self.image_mask))
        elif self.image_mask == 'Single MSA open shutter':
            # one MSA open shutter aperture
            optsys.add_image(optic=poppy.RectangularFieldStop(width=0.2, height=0.45, name=self.image_mask))
        elif self.image_mask == 'Three adjacent MSA open shutters':
            optsys.add_image(optic=NIRSpec_three_MSA_shutters(name=self.image_mask))

        if (self.pupil_mask is not None) and ('grating' in self.pupil_mask.lower()):
            # NIRSpec pupil stop at the grating appears to be a rectangle.
            # see notes and ray trace from Erin Elliot in the webbpsf-data/NIRSpec/sources directory
            optsys.add_pupil(optic=poppy.RectangleAperture(height=8.41, width=7.91, name='Pupil stop at grating wheel'))
            optsys.planes[-1].wavefront_display_hint = 'intensity'

        # Add here a second instance of the instrument WFE, representing the WFE in the
        # collimator and camera.
        if self.include_si_wfe:
            optsys.add_pupil(optic=self._si_wfe_class(self, where='spectrograph'))

        if self.mode == 'IFU' and self.options.get('ifualign_rotation', True):
            optsys.add_rotation(
                90, hide=True
            )  # Rotate by 90 degrees clockwise to match the IFUalign output convention, with slices horizontal.
            optsys.planes[-1].wavefront_display_hint = 'intensity'

        return (optsys, trySAM, SAM_box_size)

    def _get_fits_header(self, hdulist, options):
        """Format NIRSpec-like FITS headers, based on JWST DMS SRD 1 FITS keyword info"""
        super(NIRSpec, self)._get_fits_header(hdulist, options)
        hdulist[0].header['GRATING'] = ('None', 'NIRSpec grating element name')
        hdulist[0].header['APERTURE'] = (str(self.image_mask), 'NIRSpec slit aperture name')

    @JWInstrument.aperturename.setter
    def aperturename(self, value):
        """Set SIAF aperture name to new value, with validation.

        This also updates the pixelscale to the local value for that aperture, for a small precision enhancement.

        Similar to superclass function, but handles the more complex situation with NIRSpec apertures and detectors
        """
        # Explicitly update detector reference coordinates to the default for the new selected aperture,
        # otherwise old coordinates can persist under certain circumstances

        try:
            ap = self.siaf[value]
        except KeyError:
            raise ValueError(f'Aperture name {value} not a valid SIAF aperture name for {self.name}')

        # NIRSpec apertures can either be per detector (i.e. "NRS1_FULL")
        # or for the focal plane but not per detector (i.e. "NRS_FULL_IFU")

        if value[0:4] in ['NRS1', 'NRS2']:
            # this is a regular per-detector aperture, so just call the regular code in the superclass
            JWInstrument.aperturename.fset(self, value)
        else:
            # apertures that start with NRS define V2,V3 position, but not pixel coordinates and pixelscale. So we
            # still have to use a full-detector aperturename for that.
            detector_apername = self.detector + '_FULL'

            # Only update if new value is different
            if self._aperturename != value:
                # First, check some info from current settings, which we will use below as part of auto pixelscale code
                # The point is to check if the pixel scale is set to a custom or default value,
                # and if it's custom then don't override that.
                # Note, check self._aperturename first to account for the edge case when this is
                # called from __init__ before _aperturename is set
                has_custom_pixelscale = self._aperturename and (
                    self.pixelscale != self._get_pixelscale_from_apername(detector_apername)
                )

                # Now apply changes:
                self._aperturename = value
                # Update detector reference coordinates
                # self.detector_position = (ap.XSciRef, ap.YSciRef)

                # Update DetectorGeometry class
                self._detector_geom_info = DetectorGeometry(self.siaf, self._aperturename)
                _log.info(f'{self.name} SIAF aperture name updated to {self._aperturename}')

                if not has_custom_pixelscale:
                    self.pixelscale = self._get_pixelscale_from_apername(detector_apername)
                    debug_message = (
                        f'Pixelscale updated to {self.pixelscale} '
                        f'based on average X+Y SciScale at SIAF aperture {self._aperturename}'
                    )
                    _log.debug(debug_message)

                if 'IFU' in self.aperturename:
                    self._mode = 'IFU'
                    if self._disperser is None:
                        self.disperser = 'PRISM'  # Set some default spectral mode
                        self.filter = 'CLEAR'
                    if self.image_mask not in ['IFU', None]:
                        info_message = (
                            'The currently-selected image mask (slit) is not compatible with IFU mode. '
                            'Setting image_mask=None'
                        )
                        _log.info(info_message)
                        self.image_mask = None
                else:
                    self._mode = 'imaging'  # More to implement here later!
        # Update the rotation angle
        # This works the same for both regular and IFU modes
        self._rotation = self._get_aperture_rotation(self.aperturename)

    def _tel_coords(self):
        """Convert from science frame coordinates to telescope frame coordinates using
        SIAF transformations. Returns (V2, V3) tuple, in arcminutes.

        Note that the astropy.units framework is used to return the result as a
        dimensional Quantity.

        Some extra steps for NIRSpec to handle the more complicated/flexible mapping between detector and sky coordinates
        """

        if self.aperturename.startswith('NRS_'):
            # These apertures don't map directly to particular detector position in the usual way
            # Return coords for center of the aperture reference location
            return (
                np.asarray((self._detector_geom_info.aperture.V2Ref, self._detector_geom_info.aperture.V3Ref))
                / 60
                * units.arcmin
            )
        else:
            return super()._tel_coords()

    def _get_pixelscale_from_apername(self, apername):
        """Simple utility function to look up pixelscale from apername"""
        if 'IFU' in apername:
            return super()._get_pixelscale_from_apername('NRS1_FULL')
        else:
            return super()._get_pixelscale_from_apername(apername)

    def _get_aperture_rotation(self, apername):
        """Get the rotation angle of a given aperture, using values from SIAF.

        Returns ~ position angle counterclockwise from the V3 axis, in degrees
        (i.e. SIAF V3IdlYangle)

        For NIRSpec this is simple, since even the SLIT type apertures have
        V3IdlYAngle values defined.  And we don't have the complexity of
        COMPOUND type apertures that MIRI has to deal with.

        """
        return self.siaf[apername].V3IdlYAngle

    @property
    def disperser(self):
        """NIRSpec spectral dispersing element (grating or prism).
        Only applies for IFU mode sims, currently; used to help set the
        wavelength range to simulate
        """
        if self.mode == 'IFU':
            return self._disperser
        else:
            return None

    @disperser.setter
    def disperser(self, value):
        if (value is None) or (value in self.disperser_list):
            self._disperser = value
        else:
            raise RuntimeError(f'Not a valid NIRSpec disperser name: {value}')

    def _calc_psf_format_output(self, result, options):
        """Format output HDUList. In particular, add some extra metadata if in IFU mode"""
        super()._calc_psf_format_output(result, options)
        if self.mode == 'IFU':
            n_exts = len(result)
            for ext in np.arange(n_exts):
                result[ext].header['MODE'] = ('IFU', 'This is a NIRSpec IFU mode simulation')
                result[ext].header['GRATING'] = (self.disperser, 'Name of the grating (or prism) element simulated.')

    @property
    def band(self):
        if self.mode != 'IFU':
            return None

        return self.disperser + '/' + self.filter

    @band.setter
    def band(self, value):
        raise RuntimeError('This is a read-only property. Set grating and/or filter attributes instead.')


class NIRISS(JWInstrument):
    """A class modeling the optics of the Near-IR Imager and Slit Spectrograph
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
        JWInstrument.__init__(self, 'NIRISS')
        self.pixelscale = 0.065657  # Average of X and Y scales, SIAF PRDOPSSOC-059, 2022 Dec

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

        from .optics import NIRISS_CLEARP, NIRISS_GR700XD_Grism

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
            optsys.add_pupil(
                transmission=self._datapath + '/optics/MASK_NRM.fits.gz',
                name=self.pupil_mask,
                flip_y=True,
                flip_x=True,
                shift_x=shift_x,
                shift_y=shift_y,
                rotation=rotation,
            )
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'CLEARP':
            optsys.add_pupil(optic=NIRISS_CLEARP(shift_x=shift_x, shift_y=shift_y, rotation=rotation))
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'GR700XD':
            optsys.add_pupil(optic=NIRISS_GR700XD_Grism(shift_x=shift_y, shift_y=shift_y, rotation=rotation))

        elif self.pupil_mask is None and self.image_mask is not None:
            optsys.add_pupil(name='No Lyot Mask Selected!')

        return (optsys, trySAM, radius + 0.05)  # always attempt to cast this to a SemiAnalyticCoronagraph

    def _get_fits_header(self, hdulist, options):
        """Format NIRISS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info"""
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
                _log.info('NIRISS pupil obscuration updated to {0} to match ' 'the requested filter'.format(new_pupil_mask))
                self.pupil_mask = new_pupil_mask

    def _validate_config(self, **kwargs):
        """Validate instrument config for NIRISS

        For NIRISS, this optionally adjusts the instrument pupil
        """
        wavelengths = np.array(kwargs['wavelengths'])
        if np.min(wavelengths) < self.SHORT_WAVELENGTH_MIN:
            raise RuntimeError('The requested wavelengths are too short to be imaged with NIRISS')
        if np.max(wavelengths) > self.LONG_WAVELENGTH_MAX:
            raise RuntimeError('The requested wavelengths are too long to be imaged with NIRISS')
        if np.max(wavelengths) <= self.SHORT_WAVELENGTH_MAX and self.pupil == 'NRM':
            raise RuntimeError('NRM pupil can only be used with long ' 'wavelength filters (F277W and longer)')

        return super(NIRISS, self)._validate_config(**kwargs)


class FGS(JWInstrument):
    """A class modeling the optics of the FGS.

    Not a lot to see here, folks: There are no selectable options, just a great big detector-wide bandpass
    and two detectors.

    The detectors are named as FGS1, FGS2 but may synonymously also be referred to as
    GUIDER1, GUIDER2 for compatibility with DMS convention
    """

    def __init__(self):
        JWInstrument.__init__(self, 'FGS')
        self.pixelscale = 0.068991  # Average of X and Y scales for both detectors, SIAF PRDOPSSOC-059, 2022 Dec

        self.options['pupil_shift_x'] = 0.0041  # CV3 on-orbit estimate (RPT028027) + OTIS delta from predicted (037134)
        self.options['pupil_shift_y'] = -0.0023

        self._detectors = {'FGS1': 'FGS1_FULL', 'FGS2': 'FGS2_FULL'}
        self.detector = self.detector_list[0]

    def _addAdditionalOptics(self, optsys):
        raise NotImplementedError('No user-selectable optics in FGS.')

    def _get_fits_header(self, hdulist, options):
        """Format FGS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info"""
        super(FGS, self)._get_fits_header(hdulist, options)
        hdulist[0].header['FOCUSPOS'] = (0, 'FGS focus mechanism not yet modeled.')

    @JWInstrument.detector.setter  # override setter in this subclass
    def detector(self, value):
        # allow either FGS1 or GUIDER1 as synonyms
        if value.upper().startswith('GUIDER'):
            value = 'FGS' + value[-1]
        if value.upper() not in self.detector_list:
            raise ValueError('Invalid detector. Valid detector names are: {}'.format(', '.join(self.detector_list)))
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
        raise ValueError('Incorrect instrument name ' + name)


instrument.list = ['nircam', 'nirspec', 'niriss', 'miri']  # useful list for iteration


def calc_or_load_PSF(filename, inst, overwrite=False, **kwargs):
    """Utility function for loading a precomputed PSF from disk, or else
    if that files does not exist, then compute it and save to disk.

    This is useful for writing scripts that cache results - i.e. calculate the
    PSF the first time through and save it, then just load the PSF on subsequent
    iterations.

    Parameters
    ----------
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
        _log.info('Already exists, no need to recalculate: ' + filename)
        return fits.open(filename)
    else:
        return inst.calc_psf(outfile=filename, **kwargs)


#########################


@functools.lru_cache
def get_siaf_with_caching(instrname):
    """Parsing and loading the SIAF information is particularly time consuming,
    (can be >0.1 s per call, so multiple invokations can be a large overhead)
    Therefore avoid unnecessarily reloading it by caching results.
    This is a small speed optimization."""
    return pysiaf.Siaf(instrname)


class DetectorGeometry(object):
    """Utility class for converting between detector coordinates
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
        """Return detector size in pixels"""
        xdetsize = self.aperture.XDetSize
        ydetsize = self.aperture.YDetSize
        return (xdetsize, ydetsize)

    def validate_coords(self, x, y):
        """Check if specified pixel coords are actually on the detector

        Parameters
        ----------
        x, y : floats
            coordinates in pixels
        """
        if x < 0:
            raise ValueError('Detector pixels X coordinate cannot be negative.')
        if y < 0:
            raise ValueError('Detector pixels Y coordinate cannot be negative.')
        if x > int(self.shape[0]) - 1:
            raise ValueError('Detector pixels X coordinate cannot be > {0}'.format(int(self.shape[0]) - 1))
        if y > int(self.shape[1]) - 1:
            raise ValueError('Detector pixels Y coordinate cannot be > {0}'.format(int(self.shape[1]) - 1))

    def pix2angle(self, xpix, ypix):
        """Convert from science frame coordinates (in pixels) to telescope frame coordinates
        (in arcminutes) using SIAF transformations.

        See the pysiaf code for all the full details, or Lallo & Cox Tech Reports

        Parameters
        ----------
        xpix, ypix : floats
            X and Y pixel coordinates, 0 <= xpix, ypix < detector_size

        Returns
        -------
        V2, V3 : floats
            V2 and V3 coordinates, in arcMINUTES
            Note that the astropy.units framework is used to return the result as a
            dimensional Quantity.

        """

        tel_coords = np.asarray(self.aperture.sci_to_tel(xpix, ypix))
        tel_coords_arcmin = tel_coords / 60.0 * units.arcmin  # arcsec to arcmin

        return tel_coords_arcmin


#########################


def segname(val):
    """Return WSS-compliant segment name for a variety of input formats

    For instance, one particular segment can be referred to as "B3", 11, "B3-11", etc.
    The WSS refers to this segment as "B3-11".  THis function will return the string
    "B3-11" for any of the above inputs, and similarly for any of the other segments.

    Parameters
    ----------
    val : string or int
        Something that can conceivably be the name or ID of a JWST PMSA.
    """

    try:
        intval = int(val)
        # Convert integer value to string name
        if intval < 1 or intval > 19:
            raise ValueError('Integer must be between 1 and 19')
        if intval < 7:
            return 'A{0}-{0}'.format(intval)
        elif intval == 19:
            return 'SM-19'
        else:
            letter = 'B' if np.mod(intval, 2) == 1 else 'C'
            number = int(np.ceil((intval - 6) * 0.5))
            return '{0}{1}-{2}'.format(letter, number, intval)
    except ValueError:
        # it had better be a letter string
        if val.startswith('SM'):
            return 'SM-19'
        base = {'A': 0, 'B': 5, 'C': 6}
        try:
            offset = base[val[0]]
        except (KeyError, IndexError):
            raise ValueError('string must start with A, B, or C')
        try:
            num = int(val[1])
        except ValueError:
            raise ValueError('input string must have 2nd character as a number from 1-6')
        if num < 1 or num > 6:
            raise ValueError('input string must have 2nd character as a number from 1-6')
        if val[0] == 'A':
            return '{0}{1}-{1}'.format(val[0], val[1])
        else:
            return '{0}{1}-{2}'.format(val[0], val[1], offset + int(val[1]) * 2)


def one_segment_pupil(segmentname, npix=1024):
    """Return a pupil image which corresponds to only a single
    segment of the telescope. This can be useful when simulating
    early stages of JWST alignment.


    Example
    -------
    nc = webbpsf.NIRCam()
    nc.pupil = webbpsf.one_segment_pupil('B1')

    """

    # get the master pupil file, which may or may not be gzipped
    segmap = os.path.join(utils.get_webbpsf_data_path(), f'JWpupil_segments_RevW_npix{npix}.fits.gz')
    if not os.path.exists(segmap):
        # try without .gz
        segmap = os.path.join(utils.get_webbpsf_data_path(), f'JWpupil_segments_RevW_npix{npix}.fits')

    newpupil = fits.open(segmap)
    if newpupil[0].header['VERSION'] < 2:
        raise RuntimeError(f'Expecting file version >= 2 for {segmap}')

    segment_official_name = segname(segmentname)
    num = int(segment_official_name.split('-')[1])

    newpupil[0].data = np.asarray(newpupil[0].data == num, dtype=int)

    newpupil[0].header['SEGMENT'] = segment_official_name
    return newpupil
