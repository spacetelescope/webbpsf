#!/usr/bin/env python
"""
=======
WebbPSF
=======

An object-oriented modeling system for the JWST instruments.

Full documentation at http://www.stsci.edu/~mperrin/software/webbpsf/

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
import types
import glob
import time
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate, scipy.ndimage
import matplotlib

import astropy.io.fits as fits
import astropy.io.ascii as ioascii


import poppy

#from . import config
from . import conf
from . import utils


try: 
    import pysynphot
    _HAS_PYSYNPHOT = True
except:
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
    configuration can be done by editing the :ref:`SpaceTelescopeInstrument.options` dictionary, either by passing options to __init__ or by directly editing the dict afterwards.
    """
    telescope = "Generic Space Telescope"
    options = {} # options dictionary
    """ A dictionary capable of storing other arbitrary options, for extensibility. The following are all optional, and
    may or may not be meaningful depending on which instrument is selected.

    Parameters
    ----------
    source_offset_r : float
        Radial offset of the target from the center, in arcseconds
    source_offset_theta : float
        Position angle for that offset
    pupil_shift_x, pupil_shift_y : float
        Relative shift of the intermediate (coronagraphic) pupil in X and Y relative to the telescope entrace pupil, expressed as a decimal between 0.0-1.0
        Note that shifting an array too much will wrap around to the other side unphysically, but
        for reasonable values of shift this is a non-issue.  This option only has an effect for optical models that
        have something at an intermediate pupil plane between the telescope aperture and the detector. 
    pupil_rotation : float
        Relative rotation of the intermediate (coronagraphic) pupil relative to the telescope entrace pupil, expressed in degrees counterclockwise. 
        This option only has an effect for optical models that have something at an intermediate pupil plane between the telescope aperture and the detector.
    rebin : bool
        For output files, write an additional FITS extension including a version of the output array 
        rebinned down to the actual detector pixel scale?
    jitter : string
        Type of jitter model to apply. Currently not implemented
    parity : string "even" or "odd"
        You may wish to ensure that the output PSF grid has either an odd or even number of pixels.
        Setting this option will force that to be the case by increasing npix by one if necessary.
        Note that this applies to the number detector pixels, rather than the subsampled pixels if oversample>1. 
    force_coron : bool
        Set this to force full coronagraphic optical propagation when it might not otherwise take place
        (e.g. calculate the non-coronagraphic images via explicit propagation to all optical surfaces, FFTing 
        to intermediate pupil and image planes whether or not they contain any actual optics, rather than
        taking the straight-to-MFT shortcut)
    no_sam : bool
        Set this to prevent the SemiAnalyticMethod coronagraph mode from being used when possible, and instead do
        the brute-force FFT calculations. This is usually not what you want to do, but is available for comparison tests.
        The SAM code will in general be much faster than the FFT method, particularly for high oversampling.

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

    def __init__(self, name="", pixelscale = 0.064):
        self.name = name

        self._WebbPSF_basepath = utils.get_webbpsf_data_path()

        self._datapath = self._WebbPSF_basepath + os.sep + self.name + os.sep
        self._image_mask = None
        self._pupil_mask = None

        self.pupil = None
        "Filename *or* fits.HDUList for JWST pupil mask. Usually there is no need to change this."
        #TODO:jlong: is it enough to move this to JWInstr?
        self.pupilopd = None   # This can optionally be set to a tuple indicating (filename, slice in datacube)
        """Filename *or* fits.HDUList for JWST pupil OPD.


        This can be either a full absolute filename, or a relative name in which case it is
        assumed to be within the instrument's `data/OPDs/` directory, or an actual fits.HDUList object corresponding to such a file.
        If the file contains a datacube, you may set this to a tuple (filename, slice) to select a given slice, or else
        the first slice will be used."""
        self.pupil_radius = None  # Set when loading FITS file in _getOpticalSystem

        self.options = {}  # dict for storing other arbitrary options.

        # filter_list   available filter names in order by wavelength for public api
        # _filters      a dict of named tuples with name, filename, & default_nlambda
        #               with the filter name as the key
        self.filter_list, self._filters = self._get_filters()

        # choose a default filter, in case the user doesn't specify one
        self.filter = self.filter_list[0]

        self._rotation = None

        self.opd_list = [os.path.basename(os.path.abspath(f)) for f in glob.glob(self._datapath+os.sep+'OPD/OPD*.fits')]
        self.opd_list.sort()
        if len(self.opd_list) > 0:
            self.pupilopd = self.opd_list[-1]

        self._image_mask=None
        self.image_mask_list=[]
        "List of available image_masks"

        self._pupil_mask=None
        self.pupil_mask_list=[]
        "List of available pupil_masks"

        self.pixelscale = pixelscale
        "Detector pixel scale, in arcsec/pixel"
        self._spectra_cache = {}  # for caching pysynphot results.


        self.detector_list = ['Default']
        self._detector = None

        # where is the source on the detector, in 'Science frame' pixels?
        self.detector_coordinates = (0, 0)

    def _validate_config(self):
        """Validate the configuration in some sense; details depend on instrument.

        This is called automatically when assembling an OpticalSystem object prior to any
        calculation.
        """
        raise NotImplementedError("Subclasses must implement _validate_config")

    @property
    def image_mask(self):
        """Currently selected image plane mask, or None for direct imaging"""
        return self._image_mask

    @image_mask.setter
    def image_mask(self, name):
        if name is "": name = None
        if name is not None:
            if name in self.image_mask_list:
                pass # there's a perfect match, this is fine.
            else:
                name = name.upper() # force to uppercase
                if name not in self.image_mask_list: # if still not found, that's an error.
                    raise ValueError("Instrument %s doesn't have an image mask called '%s'." % (self.name, name))
        self._image_mask = name

    @property
    def pupil_mask(self):
        """Currently selected Lyot pupil mask, or None for direct imaging"""
        return self._pupil_mask

    @pupil_mask.setter
    def pupil_mask(self,name):
        if name is "":
            name = None
        if name is not None:
            if name in self.pupil_mask_list:
                pass  # there's a perfect match, this is fine.
            else:
                name = name.upper() # force to uppercase
                if name not in self.pupil_mask_list:
                    raise ValueError("Instrument %s doesn't have a pupil mask called '%s'." % (self.name, name))

        self._pupil_mask = name
 
    @property
    def detector(self):
        """Currently selected detector name (for instruments with multiple detectors)"""
        return self._detector.name

    @detector.setter
    def detector(self, detname):
        if detname is not None:
            detname = detname.upper()  # force to uppercase
            try:
                siaf_aperture_name = self._detector2siaf[detname]
            except:
                raise ValueError("Unknown name: {0} is not a valid known name for a detector "
                                 "for instrument {1}".format(detname, self.name))
            self._detector = DetectorGeometry(self.name, siaf_aperture_name, shortname=detname)

    def __str__(self):
        return "<{instrument_name}>".format(instrument_name=self.name)

    #----- actual optical calculations follow here -----
    def calcPSF(self, outfile=None, source=None, filter=None, nlambda=None, monochromatic=None,
                fov_arcsec=None, fov_pixels=None,  oversample=None, detector_oversample=None,
                fft_oversample=None, calc_oversample=None, rebin=True, clobber=True, display=False,
                return_intermediates=False, **kwargs):
        """ Compute a PSF.

        The result can either be written to disk (set outfile="filename") or else will be returned as
        an astropy.io.fits HDUList object.


        Output sampling may be specified in one of two ways: 

        1) Set `oversample=<number>`. This will use that oversampling factor beyond detector pixels
           for output images, and beyond Nyquist sampling for any FFTs to prior optical planes. 
        2) set `detector_oversample=<number>` and `fft_oversample=<other_number>`. This syntax lets
           you specify distinct oversampling factors for intermediate and final planes. This is generally
           only relevant in the case of coronagraphic calculations.

        By default, both oversampling factors are set equal to 4. This default can be changed in your
        webbpsf configuration file.

        Notes
        -----
        More advanced PSF computation options (pupil shifts, source positions, jitter, ...)
        may be set by configuring the `.options` dictionary attribute of this class.

        Parameters
        ----------
        filter : string, optional
            Filter name. Setting this is just a shortcut for setting the object's filter first, then
            calling calcPSF afterwards.
        source : pysynphot.SourceSpectrum or dict or tuple
            specification of source input spectrum. Default is a 5700 K sunlike star.
        nlambda : int
            How many wavelengths to model for broadband? 
            The default depends on how wide the filter is, as set by a lookup table in the webbpsf data distribution.
        monochromatic : float, optional
            Setting this to a wavelength value (in meters) will compute a monochromatic PSF at that 
            wavelength, overriding filter and nlambda parameters.
        fov_arcsec : float
            field of view in arcsec. Default=5
        fov_pixels : int
            field of view in pixels. This is an alternative to fov_arcsec.
        outfile : string
            Filename to write. If None, then result is returned as an HDUList
        oversample, detector_oversample, fft_oversample : int
            How much to oversample. Default=4. By default the same factor is used for final output 
            pixels and intermediate optical planes, but you may optionally use different factors 
            if so desired.
        rebin : bool, optional
            If set, the output file will contain a FITS image extension containing the PSF rebinned
            onto the actual detector pixel scale. Thus, setting oversample=<N> and rebin=True is
            the proper way to obtain high-fidelity PSFs computed on the detector scale. Default is True.
        clobber : bool
            overwrite output FITS file if it already exists?
        display : bool
            Whether to display the PSF when done or not.
        save_intermediates, return_intermediates : bool
            Options for saving to disk or returning to the calling function the intermediate optical planes during the propagation. 
            This is useful if you want to e.g. examine the intensity in the Lyot plane for a coronagraphic propagation. These have no
            effect for simple direct imaging calculations.


        For additional arguments, see the documentation for poppy.OpticalSystem.calcPSF()


        Returns
        -------
        outfits : fits.HDUList
            The output PSF is returned as a fits.HDUlist object.
            If `outfile` is set to a valid filename, the output is also written to that file.


        """

        if calc_oversample is not None: 
            raise DeprecationWarning("The calc_oversample parameter is deprecated and will be removed in webbpsf 0.4. User fft_oversample instead.")
            fft_oversample = calc_oversample # back compatibility hook for deprecated arg name.

        _log.info("Setting up PSF calculation for "+self.name)

        # first make sure that webbpsf's configuration is used to override any of the
        # same configuration options in poppy. This is admittedly perhaps overbuilt to have identical
        # settings in both packages, but the intent is to shield typical users of webbpsf
        # from having to think about the existence of the underlying library. They can 
        # just deal with one set of settings.
        #config._apply_settings_to_poppy()

        if filter is not None:
            self.filter = filter

        local_options = self.options.copy()  # all local state should be stored in a dict, for
                                      # ease of handing off to the various subroutines of
                                      # calcPSF. Don't just modify the global self.options
                                      # structure since that would pollute it with temporary
                                      # state as well as persistent state.
        local_options['monochromatic'] = monochromatic


    
        #----- choose # of wavelengths intelligently. Do this first before generating the source spectrum weighting.
        if nlambda is None or nlambda==0:
            # Automatically determine number of appropriate wavelengths.
            # Make selection based on filter configuration file
            try:
                nlambda = self._filter_nlambda_default[self.filter]
                _log.debug("Automatically selecting # of wavelengths: %d" % nlambda)
            except:
                nlambda=10
                _log.warn("Filter %s not found in lookup table for default number of wavelengths to use.. setting default nlambda=%d" % (self.filter, nlambda))
        local_options['nlambda'] = nlambda



        #----- calculate field of view depending on supplied parameters
        if fov_arcsec is None and fov_pixels is None:  #pick decent defaults.
            if self.name =='MIRI': fov_arcsec=12.
            else: fov_arcsec=5.
            fov_spec = 'arcsec = %f' % fov_arcsec
        elif fov_pixels is not None:

            if np.isscalar(fov_pixels): 
                fov_spec = 'pixels = %d' % fov_pixels
            else:
                fov_spec = 'pixels = (%d, %d)' % (fov_pixels[0], fov_pixels[1])
        elif fov_arcsec is not None:
            if np.isscalar(fov_arcsec): 
                fov_spec = 'arcsec = %f' % fov_arcsec
            else:
                fov_spec = 'arcsec = (%.3f, %.3f)' % (fov_arcsec[0], fov_arcsec[1])

        _log.debug('FOV set to '+fov_spec)

        #---- Implement the semi-convoluted logic for the oversampling options. See docstring above
        if oversample is not None and detector_oversample is not None and fft_oversample is not None:
            # all options set, contradictorily -> complain!
            raise ValueError("You cannot specify simultaneously the oversample= option with the detector_oversample and fft_oversample options. Pick one or the other!")
        elif oversample is None and detector_oversample is None and fft_oversample is None:
            # nothing set -> set oversample = 4
            oversample = conf.default_oversampling
        if detector_oversample is None: detector_oversample = oversample
        if fft_oversample is None: fft_oversample = oversample
        local_options['detector_oversample']=detector_oversample
        local_options['fft_oversample']=fft_oversample


        _log.info("PSF calc using fov_%s, oversample = %d, nlambda = %d" % (fov_spec, detector_oversample, nlambda) )

        #----- compute weights for each wavelength based on source spectrum
        wavelens, weights = self._getWeights(source=source, nlambda=nlambda, monochromatic=monochromatic)


        #---- now at last, actually do the PSF calc:
        #  instantiate an optical system using the current parameters
        self.optsys = self._getOpticalSystem(fov_arcsec=fov_arcsec, fov_pixels=fov_pixels,
            fft_oversample=fft_oversample, detector_oversample=detector_oversample, options=local_options)
        # and use it to compute the PSF (the real work happens here, in code in poppy.py)
        #result = self.optsys.calcPSF(source, display_intermediates=display, save_intermediates=save_intermediates, display=display)
        #if _USE_MULTIPROC and monochromatic is None :
            #result = self.optsys.calcPSFmultiproc(source, nprocesses=_MULTIPROC_NPROCESS) # no fancy display args for multiproc.
        #else:
        result = self.optsys.calcPSF(wavelens, weights, display_intermediates=display, display=display, return_intermediates=return_intermediates, **kwargs)

        if return_intermediates: # this implies we got handed back a tuple, so split it apart
            result, intermediates = result


        self._getFITSHeader(result, local_options)

        self._calcPSF_format_output(result, local_options)


        if display:
            f = plt.gcf()
            #p.text( 0.1, 0.95, "%s, filter= %s" % (self.name, self.filter), transform=f.transFigure, size='xx-large')

            if monochromatic is None:
                plt.suptitle( "%s, filter= %s" % (self.name, self.filter), size='xx-large')
                plt.text( 0.99, 0.04, "Calculation with %d wavelengths (%g - %g um)" % (nlambda, wavelens[0]*1e6, wavelens[-1]*1e6), transform=f.transFigure, horizontalalignment='right')
            else:
                plt.suptitle( "{self.name},  $\lambda$ = {wavelen} um".format(self=self, wavelen = monochromatic*1e6), size='xx-large')
 
        if outfile is not None:
            result[0].header["FILENAME"] = (os.path.basename (outfile), "Name of this file")
            result.writeto(outfile, clobber=clobber)
            _log.info("Saved result to "+outfile)

        if return_intermediates:
            return result, intermediates
        else:
            return result

    def _getFITSHeader(self, result, options):
        """ populate FITS Header keywords """
        poppy.Instrument._getFITSHeader(self,result, options)
        result[0].header['FILTER'] = (self.filter, 'Filter name')
        if self.image_mask is not None:
            result[0].header['CORONMSK'] = ( self.image_mask, "Image plane mask")
        if self.pupil_mask is not None:
            result[0].header['PUPIL'] = ( self.pupil_mask, "Pupil plane mask")


    def _calcPSF_format_output(self, result, options):
        """ Apply desired formatting to output file:
                 - rebin to detector pixel scale if desired
                 - set up FITS extensions if desired
                 - output either the oversampled, rebinned, or both

            Modifies the 'result' HDUList object.
        """
        output_mode = options.get('output_mode',conf.default_output_mode)

        if output_mode == 'Mock JWST DMS Output': #TODO:jlong: move out to JWInstrument
            # first rebin down to detector sampling
            # then call mockdms routines to embed in larger detector etc
            raise NotImplementedError('Not implemented yet')
        else:
            poppy.Instrument._calcPSF_format_output(self, result, options)


    def _getOpticalSystem(self,fft_oversample=2, detector_oversample = None, fov_arcsec=2, fov_pixels=None, options=None):
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

        self._validate_config()

        _log.info("Creating optical system model:")

        if options is None: options = self.options 
        if detector_oversample is None: detector_oversample = fft_oversample

        _log.debug("Oversample: %d  %d " % (fft_oversample, detector_oversample))
        optsys = poppy.OpticalSystem(
            name='{telescope}+{instrument}'.format(telescope=self.telescope, instrument=self.name),
            oversample=fft_oversample
        )
        if 'source_offset_r' in options.keys():
            optsys.source_offset_r = options['source_offset_r']
        if 'source_offset_theta' in options.keys():
            optsys.source_offset_theta = options['source_offset_theta']


        #---- set pupil OPD
        if isinstance(self.pupilopd, str):  # simple filename
            opd_map = self.pupilopd if os.path.exists( self.pupilopd) else os.path.join(self._datapath, "OPD",self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], basestring): # tuple with filename and slice
            opd_map =  (self.pupilopd[0] if os.path.exists( self.pupilopd[0]) else os.path.join(self._datapath, "OPD",self.pupilopd[0]), self.pupilopd[1])
        elif isinstance(self.pupilopd, (fits.HDUList, poppy.OpticalElement)):
            opd_map = self.pupilopd # not a path per se but this works correctly to pass it to poppy
        elif self.pupilopd is None: 
            opd_map = None
        else:
            raise TypeError("Not sure what to do with a pupilopd of that type:"+str(type(self.pupilopd)))

        #---- set pupil intensity
        if self.pupil is None:
            raise RuntimeError("The pupil shape must be specified in the "
                               "instrument class or by setting self.pupil")
        if isinstance(self.pupil, poppy.OpticalElement):
            # supply to POPPY as-is
            pupil_optic = optsys.addPupil(self.pupil)
        else:
            # wrap in an optic and supply to POPPY
            if isinstance(self.pupil, str): # simple filename
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
            #---- apply pupil intensity and OPD to the optical model
            pupil_optic = optsys.addPupil(
                name='{} Pupil'.format(self.telescope),
                transmission=pupil_transmission,
                opd=opd_map,
                opdunits='micron',
                rotation=self._rotation
            )
        self.pupil_radius = pupil_optic.pupil_diam / 2.0


        #---- Add defocus if requested
        if 'defocus_waves' in options.keys(): 
           defocus_waves = options['defocus_waves'] 
           defocus_wavelength = float(options['defocus_wavelength']) if 'defocus_wavelength' in options.keys() else 2.0e-6
           _log.info("Adding defocus of %d waves at %.2f microns" % (defocus_waves, defocus_wavelength *1e6))
           lens = poppy.ThinLens(
               name='Defocus',
               nwaves=defocus_waves,
               reference_wavelength=defocus_wavelength,
               radius=jwpupil.pupil_diam
           )
           optsys.addPupil(optic=lens)


        #---- add coronagraph or spectrograph optics if requested, and possibly flag to invoke semi-analytic coronagraphic propagation

        # first error check for null strings, which should be considered like None
        if self.image_mask == "": self.image_mask = None
        if self.pupil_mask == "": self.pupil_mask = None


        if self.image_mask is not None or self.pupil_mask is not None or ('force_coron' in options.keys() and options['force_coron']):
            _log.debug("Adding coronagraph/spectrograph optics...")
            optsys, trySAM, SAM_box_size = self._addAdditionalOptics(optsys, oversample=fft_oversample)
        else: trySAM = False

        #--- add the detector element. 
        if fov_pixels is None:
            if not np.isscalar(fov_arcsec): fov_arcsec = np.asarray(fov_arcsec) # cast to ndarray if 2D
            fov_pixels = np.round(fov_arcsec/self.pixelscale)
            if 'parity' in options.keys():
                if options['parity'].lower() == 'odd'  and np.remainder(fov_pixels,2)==0: fov_pixels +=1
                if options['parity'].lower() == 'even' and np.remainder(fov_pixels,2)==1: fov_pixels +=1
        else:
            pass

        optsys.addDetector(self.pixelscale, fov_pixels = fov_pixels, oversample = detector_oversample, name=self.name+" detector")

        #---  invoke semi-analytic coronagraphic propagation
        if trySAM and not ('no_sam' in self.options.keys() and self.options['no_sam']): # if this flag is set, try switching to SemiAnalyticCoronagraph mode. 
            _log.info("Trying to invoke switch to Semi-Analytic Coronagraphy algorithm")
            try: 
                SAM_optsys = poppy.SemiAnalyticCoronagraph(optsys, oversample=fft_oversample, occulter_box = SAM_box_size )
                _log.info("SAC OK")
                return SAM_optsys
            except ValueError as err: 
                _log.warn("Could not switch to Semi-Analytic Coronagraphy mode; invalid set of optical planes? Using default propagation instead.")
                _log.warn(str(err))
                #_log.warn("ERROR ({0}): {1}".format(errno, strerror))
                pass
 


        return optsys

    def _addAdditionalOptics(self,optsys, oversample=2):
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

    def _getSynphotBandpass(self, filtername):
        """ Return a pysynphot.ObsBandpass object for the given desired band. 

        By subclassing this, you can define whatever custom bandpasses are appropriate for
        your instrument
        """
        obsmode = '{instrument},im,{filter}'.format(instrument=self.name, filter=filtername)
        try:
            band = pysynphot.ObsBandpass(obsmode.lower())
            return band
        except:
            #TODO:jlong: what exceptions can this raise?
            _log.warn("Couldn't find filter '{}' in PySynphot, falling back to "
                      "local throughput files".format(filtername))

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
            _log.warn('The supplied file, {}, does not have a WAVEUNIT keyword. Assuming it '
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

        return band

#######  JWInstrument classes  #####


class JWInstrument(SpaceTelescopeInstrument):
    telescope = "JWST"

    def __init__(self, *args, **kwargs):
        super(JWInstrument, self).__init__(*args, **kwargs)

        self.pupil = os.path.abspath(self._datapath+"../pupil_RevV.fits")
        "Filename *or* fits.HDUList for JWST pupil mask. Usually there is no need to change this."
        self.pupilopd = None   # This can optionally be set to a tuple indicating (filename, slice in datacube)
        """Filename *or* fits.HDUList for JWST pupil OPD.

        This can be either a full absolute filename, or a relative name in which case it is
        assumed to be within the instrument's `data/OPDs/` directory, or an actual fits.HDUList object corresponding to such a file.
        If the file contains a datacube, you may set this to a tuple (filename, slice) to select a given slice, or else
        the first slice will be used."""

class MIRI(JWInstrument):
    """ A class modeling the optics of MIRI, the Mid-InfraRed Instrument.
    
    Relevant attributes include `filter`, `image_mask`, and `pupil_mask`.

    In addition to the actual filters, you may select 'MRS-IFU Ch1' to
    indicate use of the MIRI IFU in Channel 1, and so forth. In this case, the `monochromatic` attribute controls the simulated wavelength.
    Note that the pixel scale varies with channel, which is why they are implemented separately. 
    **Note: IFU to be implemented later**


    """
    def __init__(self):
        JWInstrument.__init__(self, "MIRI")
        self.pixelscale = 0.11
        self._rotation = 4.561 # Source: MIRI OBA DD, page 3-16

        self.image_mask_list = ['FQPM1065', 'FQPM1140', 'FQPM1550', 'LYOT2300', 'LRS slit']
        self.pupil_mask_list = ['MASKFQPM', 'MASKLYOT', 'P750L LRS grating']

        for i in range(4):
            self.filter_list.append('MRS-IFU Ch%d'% (i+1) )
        self.monochromatic = 8.0
        self._IFU_pixelscale = {
            'Ch1': (0.18, 0.19),
            'Ch2': (0.28, 0.19),
            'Ch3': (0.39, 0.24),
            'Ch4': (0.64, 0.27),
        }
        # The above tuples give the pixel resolution (perpendicular to the slice, along the slice).
        # The pixels are not square.

        #self._default_aperture='MIRIM_center' # reference into SIAF for ITM simulation V/O coords
        self.detector_list = ['MIRIM']
        self._detector2siaf = {'MIRIM': 'MIRIM_FULL_ILLCNTR'}
        self.detector = self.detector_list[0]

    def _validate_config(self):
        """Validate instrument config for MIRI
        """
        if self.filter.startswith("MRS-IFU"):
            raise NotImplementedError("The MIRI MRS is not yet implemented.")
        return
        #TODO:jlong: is this still useful?
        # if self.image_mask is not None or self.pupil_mask is not None:
        #     if self.filter == 'F1065C':
        #         assert self.image_mask == 'FQPM1065', 'Invalid configuration'
        #         assert self.pupil_mask == 'MASKFQPM', 'Invalid configuration'
        #     elif self.filter == 'F1140C':
        #         assert self.image_mask == 'FQPM1140', 'Invalid configuration'
        #         assert self.pupil_mask == 'MASKFQPM', 'Invalid configuration'
        #     elif self.filter == 'F1550C':
        #         assert self.image_mask == 'FQPM1550', 'Invalid configuration'
        #         assert self.pupil_mask == 'MASKFQPM', 'Invalid configuration'
        #     elif self.filter == 'F2300C':
        #         assert self.image_mask == 'LYOT2300', 'Invalid configuration'
        #         assert self.pupil_mask == 'MASKLYOT', 'Invalid configuration'
        #     else:
        #         #raise ValueError("Invalid configuration selected!")
        #         _log.warn("*"*80)
        #         _log.warn("WARNING: you appear to have selected an invalid/nonphysical "
        #                   "configuration of that instrument!")
        #         _log.warn("")
        #         _log.warn("I'm going to continue trying the calculation, but YOU are responsible "
        #                   "for interpreting")
        #         _log.warn("any results in a meaningful fashion or discarding them..")
        #         _log.warn("*"*80)


    def _addAdditionalOptics(self,optsys, oversample=2):
        """Add coronagraphic or spectrographic optics for MIRI.
        Semi-analytic coronagraphy algorithm used for the Lyot only.

        """

        # For MIRI coronagraphy, all the coronagraphic optics are rotated the same
        # angle as the instrument is, relative to the primary. So they see the unrotated
        # telescope pupil. Likewise the LRS grism is rotated but its pupil stop is not.
        #
        # We model this by just not rotating till after the coronagraph. Thus we need to
        # un-rotate the primary that was already created in _getOpticalSystem.


        defaultpupil = optsys.planes.pop() # throw away the rotated pupil we just previously added
        _log.debug('Amplitude:'+str(defaultpupil.amplitude_file))
        _log.debug('OPD:'+str(defaultpupil.opd_file))
        opd = defaultpupil.opd_file
        if hasattr(defaultpupil,'opd_slice'): opd = (defaultpupil.opd_file, defaultpupil.opd_slice) # rebuild tuple if needed to slice
        optsys.addPupil(name='JWST Pupil', transmission=defaultpupil.amplitude_file, opd=opd, opdunits='micron', rotation=None)
        #optsys.addPupil('Circle', radius=6.5/2)


        # Add image plane mask
        # For the MIRI FQPMs, we require the star to be centered not on the middle pixel, but
        # on the cross-hairs between four pixels. (Since that is where the FQPM itself is centered)
        # This is with respect to the intermediate calculation pixel scale, of course, not the
        # final detector pixel scale. 
        if (self.image_mask is not None and 'FQPM' in self.image_mask) or 'force_fqpm_shift' in self.options.keys() : optsys.addPupil( poppy.FQPM_FFT_aligner() )

        if self.image_mask == 'FQPM1065':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1065",
                opticslist = [  poppy.IdealFQPM(wavelength=10.65e-6, name=self.image_mask),
                                poppy.SquareFieldStop(size=24, angle=-self._rotation)])
            optsys.addImage(container)
            trySAM = False
        elif self.image_mask == 'FQPM1140':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1140",
                opticslist = [  poppy.IdealFQPM(wavelength=11.40e-6, name=self.image_mask),
                                poppy.SquareFieldStop(size=24, angle=-self._rotation)])
            optsys.addImage(container)
            trySAM = False
        elif self.image_mask == 'FQPM1550':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1550",
                opticslist = [  poppy.IdealFQPM(wavelength=15.50e-6, name=self.image_mask),
                                poppy.SquareFieldStop(size=24, angle=-self._rotation)])
            optsys.addImage(container)
            trySAM = False
        elif self.image_mask =='LYOT2300':
            #diameter is 4.25 (measured) 4.32 (spec) supposedly 6 lambda/D
            #optsys.addImage(function='CircularOcculter',radius =4.25/2, name=self.image_mask) 
            # Add bar occulter: width = 0.722 arcsec (or perhaps 0.74, Dean says there is ambiguity)
            #optsys.addImage(function='BarOcculter', width=0.722, angle=(360-4.76))
            # position angle of strut mask is 355.5 degrees  (no = =360 -2.76 degrees
            #optsys.addImage(function='fieldstop',size=30)
            container = poppy.CompoundAnalyticOptic(name = "MIRI Lyot Occulter",
                opticslist = [poppy.CircularOcculter(radius =4.25/2, name=self.image_mask),
                              poppy.BarOcculter(width=0.722), 
                              poppy.SquareFieldStop(size=30, angle=-self._rotation)] )
            optsys.addImage(container)
            trySAM = True
            SAM_box_size = [5,20]
        elif self.image_mask == 'LRS slit':
            # one slit, 5.5 x 0.6 arcsec in height (nominal)
            #           4.7 x 0.51 arcsec (measured for flight model. See MIRI-TR-00001-CEA)
            # 
            # Per Klaus Pontoppidan: The LRS slit is aligned with the detector x-axis, so that the dispersion direction is along the y-axis. 
            optsys.addImage(optic=poppy.RectangularFieldStop(width=5.5, height=0.6, angle=self._rotation, name= self.image_mask))
            trySAM = False
        else:
            optsys.addImage()
            trySAM = False

        if (self.image_mask is not None and 'FQPM' in self.image_mask)  or 'force_fqpm_shift' in self.options.keys() : optsys.addPupil( poppy.FQPM_FFT_aligner(direction='backward'))

        # add pupil plane mask
        if ('pupil_shift_x' in self.options.keys() and self.options['pupil_shift_x'] != 0) or \
           ('pupil_shift_y' in self.options.keys() and self.options['pupil_shift_y'] != 0):

            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
            _log.info("Setting Lyot pupil shift to %s" % (str(shift)))
        else: 
            shift = None
            #_log.info('no pupil shift!')


        #optsys.addPupil('Circle', radius=6.5/2)

        if self.pupil_mask == 'MASKFQPM':
            optsys.addPupil(transmission=self._datapath+"/optics/MIRI_FQPMLyotStop.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKLYOT':
            optsys.addPupil(transmission=self._datapath+"/optics/MIRI_LyotLyotStop.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'P750L LRS grating' or self.pupil_mask == 'P750L':
            optsys.addPupil(transmission=self._datapath+"/optics/MIRI_LRS_Pupil_Stop.fits.gz", name=self.pupil_mask, shift=shift)
        else: # all the MIRI filters have a tricontagon outline, even the non-coron ones.
            optsys.addPupil(transmission=self._WebbPSF_basepath+"/tricontagon.fits", name = 'filter cold stop', shift=shift)
            # FIXME this is probably slightly oversized? Needs to have updated specifications here.

        optsys.addRotation(self._rotation)

        return (optsys, trySAM, SAM_box_size if trySAM else None)

    def _getFITSHeader(self, hdulist, options):
        """ Format MIRI-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)

        hdulist[0].header['GRATNG14'] = ('None', 'MRS Grating for channels 1 and 4')
        hdulist[0].header['GRATNG23'] = ('None', 'MRS Grating for channels 2 and 3')
        hdulist[0].header['FLATTYPE'] = ('?', 'Type of flat field to be used: all, one, principal')
        hdulist[0].header['CCCSTATE'] = ('open', 'Contamination Control Cover state: open, closed, locked')
        if self.image_mask is not None:
            hdulist[0].header['TACQNAME'] = ('None', 'Target acquisition file name')


class NIRCam(JWInstrument):
    """ A class modeling the optics of NIRCam. 
    
    Relevant attributes include `filter`, `image_mask`, and `pupil_mask`.

    The NIRCam class is smart enough to automatically select the appropriate pixel scale for the short or long wavelength channel 
    based on whether you request a short or long wavelength filter.
 
    """
    def __init__(self):
        self.module='A'          # NIRCam A or B?
        self.pixelscale = 0.0317 # for short-wavelen channels
        self._pixelscale_short = 0.0317 # for short-wavelen channels
        self._pixelscale_long = 0.0648 # for short-wavelen channels
        JWInstrument.__init__(self, "NIRCam") # do this after setting the long & short scales.
        self.pixelscale = 0.0317 # need to redo 'cause the __init__ call will reset it to zero.

        #self.image_mask_list = ['BLC2100','BLC3350','BLC4300','WEDGESW','WEDGELW']
        self.image_mask_list = ['MASKLWB','MASKSWB','MASK210R','MASK335R','MASK430R']

        self.pupil_mask_list = ['CIRCLYOT','WEDGELYOT', 'WEAK LENS +4', 'WEAK LENS +8', 'WEAK LENS -8', 'WEAK LENS +12 (=4+8)','WEAK LENS -4 (=4-8)']

        self.filter = 'F200W' # default
        self._default_aperture='NIRCam A1 center' # reference into SIAF for ITM simulation V/O coords

        self.detector_list = ['A1','A2','A3','A4','A5', 'B1','B2','B3','B4','B5']
        self._detector2siaf = dict()
        for name in self.detector_list: self._detector2siaf[name] = 'NRC{0}_FULL_CNTR'.format(name)
        self.detector=self.detector_list[0]



    def _validate_config(self):
        """Validate instrument config for NIRCam

        For NIRCam, this checks whenever you change a filter and updates the pixelscale appropriately"""
        filtwave = float(self.filter[1:4])/100
        newscale = self._pixelscale_short if filtwave < 2.4 else self._pixelscale_long
        # update the pixel scale if it has changed *and*
        # only if the user has not already set the pixel scale to some custom value
        if newscale != self.pixelscale and (self.pixelscale == self._pixelscale_short or self.pixelscale==self._pixelscale_long):
            self.pixelscale = newscale
            _log.info("NIRCam pixel scale updated to %f arcsec/pixel to match channel for the selected filter." % self.pixelscale)


    def _addAdditionalOptics(self,optsys, oversample=2):
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

        #optsys.addImage(name='null for debugging NIRcam _addCoron') # for debugging

        if self.image_mask == 'MASK210R':
            optsys.addImage( poppy.BandLimitedCoron( kind='nircamcircular', sigma=5.253 , name=self.image_mask))
            trySAM = True
            SAM_box_size = 5.0
        elif self.image_mask == 'MASK335R':
            optsys.addImage( poppy.BandLimitedCoron(kind='nircamcircular', sigma=3.2927866 , name=self.image_mask))
            trySAM = True
            SAM_box_size = 5.0
        elif self.image_mask == 'MASK430R':
            optsys.addImage( poppy.BandLimitedCoron(kind='nircamcircular', sigma=2.588496*0.99993495 , name=self.image_mask))
            trySAM = True
            SAM_box_size = 5.0
        elif self.image_mask == 'MASKSWB':
            optsys.addImage( poppy.BandLimitedCoron(kind='nircamwedge', wavelength=2.1e-6, name=self.image_mask))
            trySAM = False #True FIXME
            SAM_box_size = [5,20]
        elif self.image_mask == 'MASKLWB':
            optsys.addImage( poppy.BandLimitedCoron(kind='nircamwedge', wavelength=4.6e-6, name=self.image_mask))
            trySAM = False #True FIXME
            SAM_box_size = [5,20]
        else:
            # no occulter selected but coronagraphic mode anyway.
            trySAM = False
            SAM_box_size = 1.0 # irrelevant but variable still needs to be set.
 
        # add pupil plane mask
        if ('pupil_shift_x' in self.options.keys() and self.options['pupil_shift_x'] != 0) or \
           ('pupil_shift_y' in self.options.keys() and self.options['pupil_shift_y'] != 0):
            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        else: shift = None


        #NIRCam as-built weak lenses, from WSS config file
        WLP4_diversity =  8.27398 # microns
        WLP8_diversity = 16.4554  # microns
        WLM8_diversity =-16.4143  # microns
        WL_wavelength =   2.12    # microns 

        #optsys.addPupil( name='null for debugging NIRcam _addCoron') # debugging
        if self.pupil_mask == 'CIRCLYOT':
            optsys.addPupil(transmission=self._datapath+"/optics/NIRCam_Lyot_Somb.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'WEDGELYOT':
            optsys.addPupil(transmission=self._datapath+"/optics/NIRCam_Lyot_Sinc.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'WEAK LENS +4':
            optsys.addPupil(poppy.ThinLens(
                name='Weak Lens +4',
                nwaves=WLP4_diversity / WL_wavelength,
                reference_wavelength=WL_wavelength*1e-6, #convert microns to meters
                radius=self.pupil_radius
            ))
        elif self.pupil_mask == 'WEAK LENS +8':
            optsys.addPupil(poppy.ThinLens(
                name='Weak Lens +8',
                nwaves=WLP8_diversity / WL_wavelength,
                reference_wavelength=WL_wavelength*1e-6,
                radius=self.pupil_radius
            ))
        elif self.pupil_mask == 'WEAK LENS -8':
            optsys.addPupil(poppy.ThinLens(
                name='Weak Lens -8',
                nwaves=WLM8_diversity / WL_wavelength,
                reference_wavelength=WL_wavelength*1e-6,
                radius=self.pupil_radius
            ))
        elif self.pupil_mask == 'WEAK LENS +12 (=4+8)':
            stack = poppy.CompoundAnalyticOptic(name='Weak Lens Stack +12', opticslist=[
                poppy.ThinLens(
                    name='Weak Lens +4',
                    nwaves=WLP4_diversity / WL_wavelength,
                    reference_wavelength=WL_wavelength*1e-6,
                    radius=self.pupil_radius
                ),
                poppy.ThinLens(
                    name='Weak Lens +8',
                    nwaves=WLP8_diversity / WL_wavelength,
                    reference_wavelength=WL_wavelength*1e-6,
                    radius=self.pupil_radius
                )]
            )
            optsys.addPupil(stack)
        elif self.pupil_mask == 'WEAK LENS -4 (=4-8)':
            stack = poppy.CompoundAnalyticOptic(name='Weak Lens Stack -4', opticslist=[
                poppy.ThinLens(
                    name='Weak Lens +4',
                    nwaves=WLP4_diversity / WL_wavelength,
                    reference_wavelength=WL_wavelength*1e-6,
                    radius=self.pupil_radius
                ),
                poppy.ThinLens(
                    name='Weak Lens -8',
                    nwaves=WLM8_diversity / WL_wavelength,
                    reference_wavelength=WL_wavelength*1e-6,
                    radius=self.pupil_radius
                )]
            )
            optsys.addPupil(stack)


        elif (self.pupil_mask is None and self.image_mask is not None):
            optsys.addPupil(name='No Lyot Mask Selected!')

        return (optsys, trySAM, SAM_box_size)

    def _getFITSHeader(self, hdulist, options):
        """ Format NIRCam-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self,hdulist, options)

        hdulist[0].header['MODULE'] = (self.module, 'NIRCam module: A or B')
        hdulist[0].header['CHANNEL'] = ( 'Short' if self.pixelscale == self._pixelscale_short else 'Long', 'NIRCam channel: long or short')
        # filter, pupil added by calcPSF header code
        hdulist[0].header['PILIN'] = ( 'False', 'Pupil imaging lens in optical path: T/F')


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
        self.pixelscale = 0.100 #  100 mas pixels. Microshutters are 0.2x0.46 but we ignore that here. 
        self._rotation = None
        self._rotation = 90+41.5  # based on SIAF docs by M. Lallo & WFS FOV doc by S. Knight
        self.filter_list.append("IFU")
        self._IFU_pixelscale = 0.1 # same.
        self.monochromatic= 3.0
        self.filter = 'F110W' # or is this called F115W to match NIRCam??

        # fixed slits
        self.image_mask = None
        self.image_mask_list = ['S200A1','S200A2','S400A1','S1600A1','S200B1', 'MSA all open', 'Single MSA open shutter', 'Three adjacent MSA open shutters']
        self.pupil_mask_list = ['NIRSpec grating']
        self.pupil_mask = self.pupil_mask_list[-1]

        #self._default_aperture='NIRSpec A center' # reference into SIAF for ITM simulation V/O coords
        self.detector_list = ['1','2']
        self._detector2siaf = dict()
        for name in self.detector_list: self._detector2siaf[name] = 'NRS{0}_FULL_CNTR'.format(name)
        self.detector=self.detector_list[0]


    def _validate_config(self):
        #if (not self.image_mask is None) or (not self.pupil_mask is None):
        #    raise ValueError('NIRSpec does not have image or pupil masks!')
        #    self.image_mask = None
        #    self.pupil_mask = None
        if self.filter.startswith("IFU"): raise NotImplementedError("The NIRSpec IFU is not yet implemented.")

    def _addAdditionalOptics(self,optsys, oversample=2):
        """ Add fixed slit optics for NIRSpec

        See Table 3-6 of NIRSpec Ops Concept Document, ESA-JWST-TN-0297 / JWST-OPS-003212

        """
        from .optics import NIRSpec_three_MSA_shutters, NIRSpec_MSA_open_grid
        trySAM = False # semi-analytic method never applicable here. 
        SAM_box_size = None

        if self.image_mask == 'S200A1' or self.image_mask == 'S200A2' or self.image_mask == 'S200B1':
            # three identical slits, 0.2 x 3.2 arcsec in length
            optsys.addImage(optic=poppy.RectangularFieldStop(width=0.2, height=3.2, name= self.image_mask + " slit"))
        elif self.image_mask == 'S400A1':
            # one slit, 0.4 x 3.65 arcsec in height
            optsys.addImage(optic=poppy.RectangularFieldStop(width=0.4, height=3.65, name= self.image_mask + " slit"))
        elif self.image_mask == 'S1600A1':
            # square aperture for exoplanet spectroscopy
            optsys.addImage(optic=poppy.RectangularFieldStop(width=1.6, height=1.6, name= self.image_mask + " square aperture"))
        elif self.image_mask == 'MSA all open':
            # all MSA shutters open 
            optsys.addImage(optic=NIRSpec_MSA_open_grid(name= self.image_mask))
        elif self.image_mask == 'Single MSA open shutter':
            # one MSA open shutter aperture 
            optsys.addImage(optic=poppy.RectangularFieldStop(width=0.2, height=0.45, name= self.image_mask))
        elif self.image_mask == 'Three adjacent MSA open shutters':
            optsys.addImage(optic=NIRSpec_three_MSA_shutters(name=self.image_mask))

 


        if ((self.pupil_mask is not None) and  ('grating' in self.pupil_mask.lower())):
            # NIRSpec pupil stop at the grating appears to be a rectangle.
            # see notes and ray trace from Erin Elliot in the webbpsf-data/NIRSpec/sources directory
            optsys.addPupil(optic=poppy.RectangleAperture(height=8.41, width=7.91,  name='Pupil stop at grating wheel'))

        #if (self.pupil_mask is None and self.image_mask is not None):
            # if we don't have a specific pupil stop, just assume for now we're
            # stopped down to a JWST like geometry
            # FIXME this is not really right - should be updated for the NIRSpec grating wheels
            #optsys.addPupil(optic=optsys[0], name='No Pupil stop provided')
            #optsys.addPupil(optic=poppy.SquareAperture(size=3.5,  name='Pupil stop at grating wheel'))



        return (optsys, trySAM, SAM_box_size)



    def _getFITSHeader(self, hdulist, options):
        """ Format NIRSpec-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)
        hdulist[0].header['GRATING'] = ( 'None', 'NIRSpec grating element name')
        hdulist[0].header['APERTURE'] = ( str(self.image_mask), 'NIRSpec slit aperture name')


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
    def __init__(self, auto_pupil=True):
        JWInstrument.__init__(self, "NIRISS", pixelscale=0.064)

        self.image_mask_list = ['CORON058', 'CORON075','CORON150','CORON200'] # available but unlikely to be used...
        self.pupil_mask_list = ['CLEAR', 'CLEARP', 'MASK_NRM','GR700XD']

        self._detector2siaf = {'NIRISS':'NIS_FULL_CNTR'}
        self.detector_list = ['NIRISS']
        self.detector=self.detector_list[0]


        self._default_aperture='NIRISS center' # reference into SIAF for ITM simulation V/O coords

        self.auto_pupil = auto_pupil

    def _validate_config(self):
        """Validate instrument config for NIRISS
        """
        pass

    def _addAdditionalOptics(self,optsys, oversample=2):
        """Add NRM or slitless spectroscopy optics for NIRISS. 

            These are probably not going to be used much in practice for NIRISS, but they
            are present, so we might as well still provide the ability to simulate 'em. 
        """

        from .optics import NIRISS_GR700XD_Grism, NIRISS_CLEARP
        if self.image_mask == 'CORON058':
            radius = 0.58/2
            optsys.addImage(function='CircularOcculter', radius=radius, name=self.image_mask)
            trySAM = True
        elif self.image_mask == 'CORON075':
            radius=0.75/2
            optsys.addImage(function='CircularOcculter', radius=radius, name=self.image_mask)
            trySAM = True
        elif self.image_mask == 'CORON150':
            radius=1.5/2
            optsys.addImage(function='CircularOcculter', radius=radius, name=self.image_mask)
            trySAM = True
        elif self.image_mask == 'CORON200':
            radius=2.0/2
            optsys.addImage(function='CircularOcculter', radius=radius, name=self.image_mask)
            trySAM = True
        else:
            trySAM = False
            radius = 0.0 # irrelevant but variable needs to be initialized

        # add pupil plane mask
        if ('pupil_shift_x' in self.options.keys() and self.options['pupil_shift_x'] != 0) or \
           ('pupil_shift_y' in self.options.keys() and self.options['pupil_shift_y'] != 0):
            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        else: shift = None

        #if self.pupil_mask == 'MASKC21N':
            #optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC21np.fits", name=self.pupil_mask, shift=shift)
        #elif self.pupil_mask == 'MASKC66N':
            #optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC66np.fits", name=self.pupil_mask, shift=shift)
        #elif self.pupil_mask == 'MASKC71N':
            #optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC71np.fits", name=self.pupil_mask, shift=shift)
        if self.pupil_mask == 'MASK_NRM':
            optsys.addPupil(transmission=self._datapath+"/optics/MASK_NRM.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'CLEAR':
            optsys.addPupil(transmission=self._datapath+"/optics/MASKCLEAR.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'CLEARP':
            #fn = os.path.join(self._datapath,"optics/MASK_CLEARP.fits.gz")
            #print fn
            #optsys.addPupil(transmission=fn, name=self.pupil_mask, shift=shift)

            # CLEARP pupil info from:
            #   MODIFIED CALIBRATION OPTIC HOLDER - NIRISS
            #   DRAWING NO 196847  REV 0  COMDEV
            #   Design file name 196847Rev0.pdf sent by Loic Albert
            # Properties:
            #  39 mm outer diam, corresponds to the circumscribing pupil of JWST
            #  2.0 mm vane width
            #  6.0 mm radius for central obstruction
            # Note the circumscribing pupil of JWST is 6603.464 mm in diameter
            #  (Ball SER on geometric optics model: BALL-JWST-SYST-05-003)

            optsys.addPupil(optic = NIRISS_CLEARP())
        elif self.pupil_mask == 'GR700XD':
            optsys.addPupil(optic = NIRISS_GR700XD_Grism(shift=shift))
 
        elif (self.pupil_mask  is None and self.image_mask is not None):
            optsys.addPupil(name='No Lyot Mask Selected!')

        return (optsys, trySAM, radius+0.05) # always attempt to cast this to a SemiAnalyticCoronagraph

    def _getFITSHeader(self, hdulist, options):
        """ Format NIRISS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)

        if self.image_mask is not None:
            hdulist[0].header['CORONPOS'] = ( self.image_mask, 'NIRISS coronagraph spot location')
        hdulist[0].header['FOCUSPOS'] = (0,'NIRISS focus mechanism not yet modeled.')


class FGS(JWInstrument):
    """ A class modeling the optics of the FGS.
    
    Not a lot to see here, folks: There are no selectable options, just a great big detector-wide bandpass.
    """
    def __init__(self):
        JWInstrument.__init__(self, "FGS")
        self.pixelscale = 0.069 # for FGS

        self.detector_list = ['1','2']
        self._detector2siaf = dict()
        for name in self.detector_list: self._detector2siaf[name] = 'FGS{0}_FULL_CNTR'.format(name)
        self.detector=self.detector_list[0]


    def _validate_config(self):
        """Validate instrument config for FGS 
        """
        # No user configurable options.
        pass
 
    def _addAdditionalOptics(self,optsys):
        raise NotImplementedError("No user-selectable optics in FGS.")

    def _getFITSHeader(self, hdulist, options):
        """ Format FGS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)
        hdulist[0].header['FOCUSPOS'] = (0,'FGS focus mechanism not yet modeled.')

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
    if name == 'miri': return MIRI()
    if name == 'nircam': return NIRCam()
    if name == 'nirspec': return NIRSpec()
    if name == 'niriss': return NIRISS()
    if name == 'tfi': return TFI()
    if name == 'fgs': return FGS()
    else: raise ValueError("Incorrect instrument name "+name)
Instrument.list = ['nircam', 'nirspec', 'niriss', 'miri'] # useful list for iteration


def calc_or_load_PSF(filename, inst, clobber=False, **kwargs):
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
        Parameters to pass to calcPSF() of that instrument. 

    Note that no validation is performed of the PSF loaded from disk to make sure it
    matches the desired properties.  This is just a quick-and-dirty unofficial/undocumented
    helper function.

    """
    if os.path.exists(filename) and not clobber:
        return fits.open(filename)
    else: 
        return inst.calcPSF(outfile = filename, **kwargs)




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
        if shortname is not None: self.name=shortname
        from jwxml import SIAF

        self.mysiaf = SIAF(instr=self.instrname, basepath=os.path.join( utils.get_webbpsf_data_path(), self.instrname) )
        self.aperture = self.mysiaf[aperturename]


    @property
    def shape(self):
        """ Return detector size in pixels """
        xdetsize = self.aperture.XDetSize
        ydetsize = self.aperture.YDetSize
        return (xdetsize,ydetsize)

    def validate_coords(self, x, y):
        """ Check if specified pixel coords are actually on the detector

        Parameters
        -----------
        x, y : floats
            coordinates in pixels
        """
        if x < 0: raise ValueError("Detector pixels X coordinate cannot be negative.")
        if y < 0: raise ValueError("Detector pixels Y coordinate cannot be negative.")
        if x > int(self.shape[0])-1: raise ValueError("Detector pixels X coordinate cannot be > {0}".format( int(self.shape[0])-1 ))
        if y > int(self.shape[1])-1: raise ValueError("Detector pixels Y coordinate cannot be > {0}".format( int(self.shape[1])-1 ))


    def pix2angle(self, xpix, ypix):
        """ Convert  from detector coordinates to telescope frame coordinates using SIAF transformations 
        See the SIAF code in jwxml for all the full details, or Lallo & Cox Tech Reports 
        
        Parameters
        ------------
        xpix, ypix : floats
            X and Y pixel coordinates, 0 <= xpix, ypix < detector_size

        Returns
        --------
        V2, V3 : floats
            V2 and V3 coordinates, in arcMINUTES

        """
        tel_coords = np.asarray( self.aperture.Sci2Tel(xpix, ypix) )
        tel_coords_arcmin = tel_coords / 60. # arcsec to arcmin



#########################3



