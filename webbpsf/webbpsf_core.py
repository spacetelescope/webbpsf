from __future__ import division, print_function, absolute_import, unicode_literals

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
import time
import six
import copy
from collections import namedtuple, OrderedDict
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate, scipy.ndimage
import matplotlib

import astropy.io.fits as fits
import astropy.io.ascii as ioascii
import astropy.units as units

import poppy

import pysiaf

from . import conf
from . import utils
from . import version
from . import optics
from . import DATA_VERSION_MIN
from . import distortion

try:
    import pysynphot

    _HAS_PYSYNPHOT = True
except ImportError:
    _HAS_PYSYNPHOT = False

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
        the telescope entrace pupil, expressed in degrees counterclockwise.
        This option only has an effect for optical models that have something at 
        an intermediate pupil plane between the telescope aperture and the detector.
    rebin : bool
        For output files, write an additional FITS extension including a version of the output array
        rebinned down to the actual detector pixel scale?
    jitter : string "gaussian" or None
        Type of jitter model to apply. Currently only convolution with a Gaussian kernel of specified
        width `jitter_sigma` is implemented. (default: None)
    jitter_sigma : float
        Width of the jitter kernel in arcseconds (default: 0.007 arcsec)
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
        self.pupil_radius = None  # Set when loading FITS file in _get_optical_system

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
        self._spectra_cache = {}  # for caching pysynphot results.

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
        if name is "": name = None
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
        if name is "":
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

    @property
    def detector_list(self):
        """Detectors on which the simulated PSF could lie"""
        return sorted(self._detectors.keys())

    @property
    def detector_position(self):
        """The pixel position in (X, Y) on the detector"""
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

    def _get_fits_header(self, result, options):
        """ populate FITS Header keywords """
        super(SpaceTelescopeInstrument, self)._get_fits_header(result, options)
        result[0].header['FILTER'] = (self.filter, 'Filter name')
        if self.image_mask is not None:
            result[0].header['CORONMSK'] = (self.image_mask, "Image plane mask")
        if self.pupil_mask is not None:
            result[0].header['PUPIL'] = (self.pupil_mask, "Pupil plane mask")

        result[0].header['VERSION'] = (version.version, "WebbPSF software version")
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

    def _get_optical_system(self, fft_oversample=2, detector_oversample=None,
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

        # ---- set pupil OPD
        if isinstance(self.pupilopd, six.string_types):  # simple filename
            opd_map = self.pupilopd if os.path.exists(self.pupilopd) else \
                      os.path.join(self._datapath, "OPD", self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], six.string_types):
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
            pupil_optic = optsys.add_pupil(self.pupil)
        else:
            # wrap in an optic and supply to POPPY
            if isinstance(self.pupil, six.string_types):  # simple filename
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
            pupil_optic = optsys.add_pupil(
                name='{} Entrance Pupil'.format(self.telescope),
                transmission=pupil_transmission,
                opd=opd_map,
                # rotation=self._rotation
            )
        pupil_rms_wfe_nm = np.sqrt(np.mean(pupil_optic.opd[pupil_optic.amplitude == 1] ** 2)) * 1e9
        self._extra_keywords['TEL_WFE'] = (pupil_rms_wfe_nm, '[nm] Telescope pupil RMS wavefront error')
        self.pupil_radius = pupil_optic.pupil_diam / 2.0

        # add coord transform from entrance pupil to exit pupil
        optsys.add_inversion(axis='y', name='OTE exit pupil', hide=True)

        # add rotation at this point, if present - needs to be after the
        # exit pupil inversion.
        if self._rotation is not None:
            optsys.add_rotation(self._rotation, hide=True)
            optsys.planes[-1].wavefront_display_hint = 'intensity'

        # Allow instrument subclass to add field-dependent aberrations
        aberration_optic = self._get_aberrations()
        if aberration_optic is not None:
            optsys.add_pupil(aberration_optic)

            try:
                inst_rms_wfe_nm = np.sqrt(np.mean(aberration_optic.opd[aberration_optic.amplitude == 1] ** 2)) * 1e9
                self._extra_keywords['SI_WFE'] = (inst_rms_wfe_nm, '[nm] instrument pupil RMS wavefront error')
            except TypeError:
                pass

            if hasattr(aberration_optic, 'header_keywords'):
                self._extra_keywords.update(aberration_optic.header_keywords())

        # ---- Add defocus if requested
        if 'defocus_waves' in options:
            defocus_waves = options['defocus_waves']
            defocus_wavelength = float(options['defocus_wavelength']) if 'defocus_wavelength' in options else 2.0e-6
            _log.info("Adding defocus of %d waves at %.2f microns" % (defocus_waves, defocus_wavelength * 1e6))
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

        if self.image_mask is not None or self.pupil_mask is not None or \
                ('force_coron' in options and options['force_coron']):
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
        """ Return a pysynphot.ObsBandpass object for the given desired band.

        By subclassing this, you can define whatever custom bandpasses are appropriate for
        your instrument
        """

        # Excise never-in-practice-used code path with ObsBandpass
        # see https://github.com/mperrin/webbpsf/issues/51
        #  Leaving this code here for now, just commented out, in case we ever decide to
        #  implement HST modes a la effectively porting TinyTim to Python...
        #
        # obsmode = '{instrument},im,{filter}'.format(instrument=self.name, filter=filtername)
        # try:
        #    band = pysynphot.ObsBandpass(obsmode.lower())
        #    return band
        # except (ValueError, TypeError) as e:
        #    _log.debug("Couldn't find filter '{}' in PySynphot, falling back to "
        #               "local throughput files".format(filtername))
        #    _log.debug("Underlying PySynphot exception was: {}".format(e))

        # the requested band is not yet supported in synphot/CDBS. (those files are still a
        # work in progress...). Therefore, use our local throughput files and create a synphot
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
            band = pysynphot.spectrum.ArraySpectralElement(
                throughput=filterdata.THROUGHPUT, wave=filterdata.WAVELENGTH,
                waveunits=waveunit, name=filtername
            )
        except AttributeError:
            raise ValueError("The supplied file, %s, does not appear to be a FITS table "
                             "with WAVELENGTH and THROUGHPUT columns." % filter_info.filename)

        filterfits.close()
        return band


#######  JWInstrument classes  #####


class JWInstrument(SpaceTelescopeInstrument):
    """ Superclass for all JWST instruments

    Notable attributes:

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

        opd_path = os.path.join(self._datapath, 'OPD')
        self.opd_list = []
        for filename in glob.glob(os.path.join(opd_path, 'OPD*.fits*')):
            self.opd_list.append(os.path.basename(os.path.abspath(filename)))

        if not len(self.opd_list) > 0:
            raise RuntimeError("No pupil OPD files found for {name} in {path}".format(name=self.name, path=opd_path))

        self.opd_list.sort()
        self.pupilopd = self.opd_list[-1]

        self.pupil = os.path.abspath(os.path.join(
            self._WebbPSF_basepath,
            "jwst_pupil_RevW_npix1024.fits.gz"
        ))
        "Filename *or* fits.HDUList for JWST pupil mask. Usually there is no need to change this."

        self._detector = None

        # where is the source on the detector, in 'Science frame' pixels?
        self.detector_position = (1024, 1024)

        self.include_si_wfe = True
        """Should calculations include the Science Instrument internal WFE?"""
        self.options['jitter'] = 'gaussian'
        self.options['jitter_sigma'] = 0.007

        # class name to use for SI internal WFE, which can be overridden in subclasses
        self._si_wfe_class = optics.WebbFieldDependentAberration

    def _get_default_fov(self):
        """ Return default FOV in arcseconds """
        return 5  # default for all NIR instruments

    def _get_optical_system(self, fft_oversample=2, detector_oversample=None, fov_arcsec=2, fov_pixels=None, options=None):
        # invoke superclass version of this
        # then add a few display tweaks
        optsys = SpaceTelescopeInstrument._get_optical_system(self,
                                                              fft_oversample=fft_oversample,
                                                              detector_oversample=detector_oversample,
                                                              fov_arcsec=fov_arcsec, fov_pixels=fov_pixels,
                                                              options=options)
        optsys.planes[0].display_annotate = utils.annotate_ote_entrance_coords
        return optsys

    def _get_aberrations(self):
        """ Compute field-dependent aberration for a given instrument
        based on a lookup table of Zernike coefficients derived from
        ISIM cryovac test data.

        This is a very preliminary version!
        """
        if not self.include_si_wfe:
            return None

        optic = self._si_wfe_class(self)

        return optic

    @SpaceTelescopeInstrument.detector.setter  # override setter in this subclass
    def detector(self, value):
        if value.upper() not in self.detector_list:
            raise ValueError("Invalid detector. Valid detector names are: {}".format(', '.join(self.detector_list)))
        self._detector = value.upper()
        self._detector_geom_info = DetectorGeometry(self.name, self._detectors[self._detector])

    def _tel_coords(self):
        """ Convert from detector pixel coordinates to SIAF aperture coordinates,

        Returns (V2, V3) tuple, in arcminutes.
        Note that the astropy.units framework is used to return the result as a
        dimensional Quantity. """

        return self._detector_geom_info.pix2angle(self.detector_position[0], self.detector_position[1])

    def _xan_yan_coords(self):
        """ Convert from detector pixel coordinates to the XAN, YAN coordinate syste
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
        siaf = pysiaf.Siaf(self.name)
        try:
            ap = siaf[aperture_name]

            self.detector_position = (ap.XDetRef, ap.YDetRef)
            detname = aperture_name.split('_')[0]
            self.detector = detname # As a side effect this auto reloads SIAF info, see detector.setter
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
        result[0].header["APERNAME"] = (self._detectors[self._detector], "SIAF aperture name")

    def calc_psf(self, outfile=None, source=None, nlambda=None, monochromatic=None,
                 fov_arcsec=None, fov_pixels=None, oversample=None, detector_oversample=None, fft_oversample=None,
                 overwrite=True, display=False, save_intermediates=False, return_intermediates=False,
                 normalize='first', add_distortion=True, crop_psf=True):
        """ Compute a PSF and apply distortion effects.

        Parameters
        ----------
       add_distortion : bool
            If True, will add 2 new extensions to the PSF HDUlist object. The 2nd extension
            will be a distorted version of the over-sampled PSF and the 3rd extension will
            be a distorted version of bthe detector-sampled PSF.
        crop_psf : bool
            If True, when the PSF is rotated to match the detector's rotation in the focal
            plane, the PSF will be cropped so the shape of the distorted PSF will match it's
            undistorted counterpart. This will only be used for NIRCam, NIRISS, and FGS PSFs.


        """

        # Run poppy calc_psf
        psf = SpaceTelescopeInstrument.calc_psf(self, outfile=outfile, source=source, nlambda=nlambda,
                                                monochromatic=monochromatic, fov_arcsec=fov_arcsec,
                                                fov_pixels=fov_pixels, oversample=oversample,
                                                detector_oversample=detector_oversample, fft_oversample=fft_oversample,
                                                overwrite=overwrite, display=display,
                                                save_intermediates=save_intermediates,
                                                return_intermediates=return_intermediates, normalize=normalize)

        if return_intermediates:
            psf, intermediates = psf

        # If chosen to add distortion
        if add_distortion:
            if self.image_mask == "LRS slit" and self.pupil_mask == "P750L LRS grating":
                raise NotImplementedError("Distortion is not implemented yet for MIRI LRS mode.")

            # Set up new extensions to add distortion to
            for ext in [0, 1]:
                hdu_new = fits.ImageHDU(psf[ext].data, psf[ext].header)  # these will be the PSFs that are edited
                psf.append(hdu_new)
                ext_new = ext + 2
                psf[ext_new].header["EXTNAME"] = psf[ext].header["EXTNAME"][0:4] + "DIST"  # change extension name

            # Apply distortions based on the instrument
            if self.name in ["NIRCam", "NIRISS", "FGS"]:
                # Apply distortion effects: Rotation and Detector Distortion
                psf_rotated = distortion.apply_rotation(psf, crop=crop_psf)  # apply rotation
                psf_distorted = distortion.apply_distortion(psf_rotated)  # apply siaf distortion
            elif self.name == "MIRI":
                # Apply distortion effects to MIRI psf: SIAF and MIRI Scattering
                psf_siaf = distortion.apply_distortion(psf)  # apply siaf distortion
                psf_distorted = distortion.apply_miri_scattering(psf_siaf)  # apply miri scattering detector_effect
            elif self.name == "NIRSpec":
                # Apply distortion effects to NIRSpec psf: SIAF only
                psf_distorted = distortion.apply_distortion(psf)  # apply siaf distortion

            # If keyword set, save out the PSF object again, with all it's extensions
            if outfile is not None:
                psf_distorted.writeto(outfile, overwrite=True, output_verify='ignore')
                # already created in the 1st calc_psf
                _log.info("Re-saved result with distortion to " + outfile)

            psf = psf_distorted

        if return_intermediates:
            return psf, intermediates

        else:
            return psf

    calcPSF = calc_psf

    # Allow users to see poppy calc_psf docstring too
    ind0 = calc_psf.__doc__.index("add_distortion")  # pull the new parameters
    ind1 = SpaceTelescopeInstrument.calc_psf.__doc__.index("Returns")  # pull where the parameters list ends
    calc_psf.__doc__ = SpaceTelescopeInstrument.calc_psf.__doc__[0:ind1] + calc_psf.__doc__[ind0:] + \
                       SpaceTelescopeInstrument.calc_psf.__doc__[ind1:]

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

    def load_was_opd(self, inputWasOpd, size=1024, save=False, filename='new_was_opd.fits'):
        """ Load and interpolate an OPD from the WAS.

        Ingests a WAS OPD and interpolates it to the proper size for WebbPSF.


        Parameters
        ----------
        HDUlist_or_filename : string
            Either a fits.HDUList object or a filename of a FITS file on disk
        size: int, optional
            Desired size of the output OPD. Default is 1024.
        save: bool, optional
            Save the interpolated OPD if True. Default is False.
        filename : string, optional
            Filename of the output OPD, if 'save' is True. Default is 'new_was_opd.fits'.

        Returns
        ---------
        HDUlist : string
           fits.HDUList object of the interpolated OPD


        """

        wasopd = fits.open(inputWasOpd)
        arrayOPD = wasopd[1].data
        dim = arrayOPD.shape[0]
        hdr = wasopd[0].header
        print("Converting {:s} from {:d}x{:d} to 1024x1024".format(inputWasOpd, dim, dim))

        hdr["BUNIT"] = 'micron'
        newopd = -1. * self.interpolate_was_opd(arrayOPD, size)  # negative sign by convention

        if save:
            outhdu = fits.HDUList()
            outhdu.append(fits.ImageHDU(newopd, header=hdr))
            outhdu.writeto(filename, clobber=True)
            outhdu.close()

        return fits.HDUList(fits.ImageHDU(newopd, header=hdr))


class MIRI(JWInstrument):
    """ A class modeling the optics of MIRI, the Mid-InfraRed Instrument.

    Relevant attributes include `filter`, `image_mask`, and `pupil_mask`.

    The pupil will auto-select appropriate values for the coronagraphic filters
    if the auto_pupil attribute is set True (which is the default).

    """

    def __init__(self):
        self.auto_pupil = True
        JWInstrument.__init__(self, "MIRI")
        self.pixelscale = 0.1110  # Source: SIAF PRDDEVSOC-D-012, 2016 April
        self._rotation = 4.4497  # Source: SIAF PRDOPSSOC-H-014

        self.options['pupil_shift_x'] = -0.0069 # CV3 on-orbit estimate (RPT028027) + OTIS delta from predicted (037134)
        self.options['pupil_shift_y'] = -0.0027

        self.image_mask_list = ['FQPM1065', 'FQPM1140', 'FQPM1550', 'LYOT2300', 'LRS slit']
        self.pupil_mask_list = ['MASKFQPM', 'MASKLYOT', 'P750L LRS grating']

        self._image_mask_apertures = {'FQPM1065': 'MIRIM_CORON1065',
                                      'FQPM1140': 'MIRIM_CORON1140',
                                      'FQPM1550': 'MIRIM_CORON1550',
                                      'LYOT2300': 'MIRIM_CORONLYOT'}
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
        # un-rotate the primary that was already created in _get_optical_system.
        # This approach is required computationally so we can work in an unrotated frame
        # aligned with the FQPM axes.

        defaultpupil = optsys.planes.pop(2)  # throw away the rotation of the entrance pupil we just added

        if self.include_si_wfe:
            # temporarily remove the SI internal aberrations
            # from the system - will add back in after the
            # coronagraph planes.
            miri_aberrations = optsys.planes.pop(2)

        # _log.debug('Amplitude:'+str(defaultpupil.amplitude_file))
        # _log.debug('OPD:'+str(defaultpupil.opd_file))
        # opd = defaultpupil.opd_file
        # if hasattr(defaultpupil,'opd_slice'):
        #    opd = (defaultpupil.opd_file, defaultpupil.opd_slice) # rebuild tuple if needed to slice
        # optsys.add_pupil(name='JWST Entrance Pupil',
        #        transmission=defaultpupil.amplitude_file, opd=opd, rotation=None,
        #        index=0)

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
        offsets = {'shift_x': self.options.get('coron_offset_x', None),
                   'shift_y': self.options.get('coron_offset_y', None)}

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
            container = poppy.CompoundAnalyticOptic(name="MIRI FQPM 1140",
                                            opticslist=[poppy.IdealFQPM(wavelength=11.40e-6, name=self.image_mask),
                                                        poppy.SquareFieldStop(size=24, rotation=self._rotation)])
            optsys.add_image(container)
            trySAM = False
        elif self.image_mask == 'FQPM1550':
            container = poppy.CompoundAnalyticOptic(name="MIRI FQPM 1550",
                                            opticslist=[poppy.IdealFQPM(wavelength=15.50e-6, name=self.image_mask),
                                                        poppy.SquareFieldStop(size=24, rotation=self._rotation)])
            optsys.add_image(container)
            trySAM = False
        elif self.image_mask == 'LYOT2300':
            # diameter is 4.25 (measured) 4.32 (spec) supposedly 6 lambda/D
            # optsys.add_image(function='CircularOcculter',radius =4.25/2, name=self.image_mask)
            # Add bar occulter: width = 0.722 arcsec (or perhaps 0.74, Dean says there is ambiguity)
            # optsys.add_image(function='BarOcculter', width=0.722, angle=(360-4.76))
            # position angle of strut mask is 355.5 degrees  (no = =360 -2.76 degrees
            # optsys.add_image(function='fieldstop',size=30)
            container = poppy.CompoundAnalyticOptic(name="MIRI Lyot Occulter",
                                            opticslist=[poppy.CircularOcculter(radius=4.25 / 2, name=self.image_mask),
                                                        poppy.BarOcculter(width=0.722),
                                                        poppy.SquareFieldStop(size=30, rotation=self._rotation)])
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
        if ('pupil_shift_x' in self.options and self.options['pupil_shift_x'] != 0) or \
                ('pupil_shift_y' in self.options and self.options['pupil_shift_y'] != 0):

            shift = (self.options.get('pupil_shift_x', 0),
                     self.options.get('pupil_shift_y', 0))
            _log.info("Setting Lyot pupil shift to %s" % (str(shift)))
        else:
            shift = None
        rotation = self.options.get('pupil_rotation', None)

        if self.pupil_mask == 'MASKFQPM':
            optsys.add_pupil(transmission=self._datapath + "/optics/MIRI_FQPMLyotStop.fits.gz",
                             name=self.pupil_mask,
                             flip_y=True, shift=shift, rotation=rotation)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'MASKLYOT':
            optsys.add_pupil(transmission=self._datapath + "/optics/MIRI_LyotLyotStop.fits.gz",
                             name=self.pupil_mask,
                             flip_y=True, shift=shift, rotation=rotation)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'P750L LRS grating' or self.pupil_mask == 'P750L':
            optsys.add_pupil(transmission=self._datapath + "/optics/MIRI_LRS_Pupil_Stop.fits.gz",
                             name=self.pupil_mask,
                             flip_y=True, shift=shift, rotation=rotation)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        else:  # all the MIRI filters have a tricontagon outline, even the non-coron ones.
            optsys.add_pupil(transmission=self._WebbPSF_basepath + "/tricontagon.fits.gz",
                             name='filter cold stop', shift=shift, rotation=rotation)
            # FIXME this is probably slightly oversized? Needs to have updated specifications here.

        if self.include_si_wfe:
            # now put back in the aberrations we grabbed above.
            optsys.add_pupil(miri_aberrations)

        optsys.add_rotation(self._rotation, hide=True)
        optsys.planes[-1].wavefront_display_hint = 'intensity'

        return (optsys, trySAM, SAM_box_size if trySAM else None)

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

    Special Options:
    The 'bar_offset' option allows specification of an offset position
    along one of the coronagraph bar occulters, in arcseconds.
    ```
    nc.image_mask = 'MASKLWB'
    nc.options['bar_offset'] = 3 # 3 arcseconds towards the right (narrow end on module A)
    ```

    The 'nd_squares' option allows toggling on and off the ND squares for TA in the simulation.
    Note that these of course aren't removable in the real instrument; this option exists solely for
    some simulation purposes.


    """
    SHORT_WAVELENGTH_MIN = 0.6 * 1e-6
    SHORT_WAVELENGTH_MAX = LONG_WAVELENGTH_MIN = 2.35 * 1e-6
    LONG_WAVELENGTH_MAX = 5.3 * 1e-6

    def __init__(self):
        self._module = 'A'  # NIRCam A or B?
        self._pixelscale_short = 0.0311  # for short-wavelen channels, SIAF PRDDEVSOC-D-012, 2016 April
        self._pixelscale_long = 0.0630  # for long-wavelen channels,  SIAF PRDDEVSOC-D-012, 2016 April
        self.pixelscale = self._pixelscale_short

        self.options['pupil_shift_x'] = 0  # Set to 0 since NIRCam FAM corrects for PM shear in flight
        self.options['pupil_shift_y'] = 0

        # need to set up a bunch of stuff here before calling superclass __init__
        # so the overridden filter setter will work successfully inside that.
        self.auto_channel = True
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

        self._si_wfe_class = optics.NIRCamFieldAndWavelengthDependentAberration

    @property
    def module(self):
        return self._detector[3]

    @property
    def channel(self):
        return 'long' if self.detector.endswith('5') else 'short'
        # note, you can't set channel directly; it's inferred based on detector.

    @channel.setter
    def channel(self, value):
        raise RuntimeError("NIRCam channel is not directly settable; set filter or detector instead.")

    @JWInstrument.detector.setter # override setter in this subclass
    def detector(self, value):
        """ Set detector, including reloading the relevant info from SIAF """
        if value.upper() not in self.detector_list:
            raise ValueError("Invalid detector. Valid detector names are: {}".format(', '.join(self.detector_list)))
        # set the channel based on the requested detector
        new_channel = 'long' if value[-1] == '5' else 'short'
        self._switch_channel(new_channel)
        self._detector = value.upper()

        self._detector_geom_info = DetectorGeometry(self.name, self._detectors[self._detector])

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

        if self.auto_channel:
            # set the channel (via setting the detector) based on filter
            wlnum = int(self.filter[1:4])
            new_channel = 'long' if wlnum >= 250 else 'short'
            self._switch_channel(new_channel)

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
        elif ((self.pupil_mask is not None) and ('LENS' not in self.pupil_mask.upper())):
            # no occulter selected but coronagraphic mode anyway. E.g. off-axis PSF
            # but don't add this image plane for weak lens calculations
            optsys.add_image(poppy.ScalarTransmission(name='No Image Mask Selected!'), index=2)
            trySAM = False
        else:
            trySAM = False

        # add pupil plane mask
        if ('pupil_shift_x' in self.options and self.options['pupil_shift_x'] != 0) or \
                ('pupil_shift_y' in self.options and self.options['pupil_shift_y'] != 0):
            shift = (self.options.get('pupil_shift_x', 0),
                     self.options.get('pupil_shift_y', 0))
        else:
            shift = None
        rotation = self.options.get('pupil_rotation', None)

        # NIRCam as-built weak lenses, from WSS config file
        WLP4_diversity = 8.27309  # microns
        WLP8_diversity = 16.4554  # microns
        WLM8_diversity = -16.4143  # microns
        WL_wavelength = 2.12  # microns

        if self.pupil_mask == 'CIRCLYOT' or self.pupil_mask == 'MASKRND':
            optsys.add_pupil(transmission=self._datapath + "/optics/NIRCam_Lyot_Somb.fits.gz", name=self.pupil_mask,
                             flip_y=True, shift=shift, rotation=rotation, index=3)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'WEDGELYOT' or self.pupil_mask == 'MASKSWB' or self.pupil_mask == 'MASKLWB':
            optsys.add_pupil(transmission=self._datapath + "/optics/NIRCam_Lyot_Sinc.fits.gz", name=self.pupil_mask,
                             flip_y=True, shift=shift, rotation=rotation, index=3)
            optsys.planes[-1].wavefront_display_hint = 'intensity'
        elif self.pupil_mask == 'WEAK LENS +4' or self.pupil_mask == 'WLP4':
            optsys.add_pupil(poppy.ThinLens(
                name='Weak Lens +4',
                nwaves=WLP4_diversity / WL_wavelength,
                reference_wavelength=WL_wavelength * 1e-6,  # convert microns to meters
                radius=self.pupil_radius,
                shift=shift, rotation=rotation,
            ), index=3)
        elif self.pupil_mask == 'WEAK LENS +8' or self.pupil_mask == 'WLP8':
            optsys.add_pupil(poppy.ThinLens(
                name='Weak Lens +8',
                nwaves=WLP8_diversity / WL_wavelength,
                reference_wavelength=WL_wavelength * 1e-6,
                radius=self.pupil_radius,
                shift=shift, rotation=rotation,
            ), index=3)
        elif self.pupil_mask == 'WEAK LENS -8' or self.pupil_mask == 'WLM8':
            optsys.add_pupil(poppy.ThinLens(
                name='Weak Lens -8',
                nwaves=WLM8_diversity / WL_wavelength,
                reference_wavelength=WL_wavelength * 1e-6,
                radius=self.pupil_radius,
                shift=shift, rotation=rotation,
            ), index=3)
        elif self.pupil_mask == 'WEAK LENS +12 (=4+8)' or self.pupil_mask == 'WLP12':
            stack = poppy.CompoundAnalyticOptic(name='Weak Lens Pair +12', opticslist=[
                poppy.ThinLens(
                    name='Weak Lens +4',
                    nwaves=WLP4_diversity / WL_wavelength,
                    reference_wavelength=WL_wavelength * 1e-6,
                    radius=self.pupil_radius,
                    shift=shift, rotation=rotation,
                ),
                poppy.ThinLens(
                    name='Weak Lens +8',
                    nwaves=WLP8_diversity / WL_wavelength,
                    reference_wavelength=WL_wavelength * 1e-6,
                    radius=self.pupil_radius,
                    shift=shift, rotation=rotation,
                )]
                                                )
            optsys.add_pupil(stack, index=3)
        elif self.pupil_mask == 'WEAK LENS -4 (=4-8)' or self.pupil_mask == 'WLM4':
            stack = poppy.CompoundAnalyticOptic(name='Weak Lens Pair -4', opticslist=[
                poppy.ThinLens(
                    name='Weak Lens +4',
                    nwaves=WLP4_diversity / WL_wavelength,
                    reference_wavelength=WL_wavelength * 1e-6,
                    radius=self.pupil_radius,
                    shift=shift, rotation=rotation,
                ),
                poppy.ThinLens(
                    name='Weak Lens -8',
                    nwaves=WLM8_diversity / WL_wavelength,
                    reference_wavelength=WL_wavelength * 1e-6,
                    radius=self.pupil_radius,
                    shift=shift, rotation=rotation,
                )]
                                                )
            optsys.add_pupil(stack, index=3)


        elif (self.pupil_mask is None and self.image_mask is not None):
            optsys.add_pupil(poppy.ScalarTransmission(name='No Lyot Mask Selected!'), index=3)
        else:
            optsys.add_pupil(transmission=self._WebbPSF_basepath + "/tricontagon_oversized_4pct.fits.gz",
                             name='filter stop', shift=shift, rotation=rotation)

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
        self.pixelscale = 0.1043  # Average over both detectors.  SIAF PRDDEVSOC-D-012, 2016 April
        # Microshutters are 0.2x0.46 but we ignore that here.
        self._rotation = 138.4  # Average for both detectors in SIAF PRDOPSSOC-H-014
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
        self.pixelscale = 0.0656  # SIAF PRDDEVSOC-D-012, 2016 April

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
        # if ('pupil_shift_x' in self.options and self.options['pupil_shift_x'] != 0) or \
        # ('pupil_shift_y' in self.options and self.options['pupil_shift_y'] != 0):
        # shift_x = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        # else:
        # shift = None
        shift_x = self.options.get('pupil_shift_x', None)
        shift_y = self.options.get('pupil_shift_y', None)
        rotation = self.options.get('pupil_rotation', None)

        # Note - the syntax for specifying shifts is different between FITS files and
        # AnalyticOpticalElement instances. Annoying but historical.
        if self.pupil_mask == 'MASK_NRM':
            optsys.add_pupil(transmission=self._datapath + "/optics/MASK_NRM.fits.gz", name=self.pupil_mask,
                             flip_y=True, shift=(shift_y, shift_y), rotation=rotation)
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
    """

    def __init__(self):
        JWInstrument.__init__(self, "FGS")
        self.pixelscale = 0.0691  # SIAF PRDDEVSOC-D-012, 2016 April

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


###########################################################################
# Generic utility functions


def Instrument(name):
    """This is just a convenience function, allowing one to access instrument objects based on a string.
    For instance,

    >>> t = Instrument('NIRISS')



    Parameters
    ----------
    name : string
        Name of the instrument class to return. Case insensitive.

    """
    name = name.lower()
    if name == 'miri':
        return MIRI()
    elif name == 'nircam':
        return NIRCam()
    elif name == 'nirspec':
        return NIRSpec()
    elif name == 'niriss':
        return NIRISS()
    elif name == 'fgs':
        return FGS()
    else:
        raise ValueError("Incorrect instrument name " + name)


Instrument.list = ['nircam', 'nirspec', 'niriss', 'miri']  # useful list for iteration


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
    """

    def __init__(self, instrname, aperturename, shortname=None):
        self.instrname = instrname
        self.name = aperturename
        if shortname is not None:
            self.name = shortname

        self.mysiaf = pysiaf.Siaf(self.instrname)
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
        """ Convert  from detector coordinates to telescope frame coordinates using SIAF transformations
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
        Something that can concievably be the name or ID of a JWST PMSA.
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


def one_segment_pupil(segmentname):
    """ Return a pupil image which corresponds to only a single
    segment of the telescope. This can be useful when simulating
    early stages of JWST alignment.


    Example
    -------
    nc = webbpsf.NIRCam()
    nc.pupil = webbpsf.one_segment_pupil('B1')

    """

    # get the master pupil file

    segmap = os.path.join(utils.get_webbpsf_data_path(), "JWpupil_segments.fits")

    newpupil = fits.open(segmap)
    if newpupil[0].header['VERSION'] < 2:
        raise RuntimeError("Expecting file version >= 2 for JWpupil_segments.fits")

    segment_official_name = segname(segmentname)
    num = int(segment_official_name.split('-')[1])

    newpupil[0].data = np.asarray(newpupil[0].data == num, dtype=int)

    newpupil[0].header['SEGMENT'] = segment_official_name
    return newpupil
