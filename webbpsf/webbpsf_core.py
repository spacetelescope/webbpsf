#!/usr/bin/env python
"""
=======
WebbPSF
=======

An object-oriented modeling system for the JWST instruments.

Full documentation at http://www.stsci.edu/~mperrin/software/webbpsf/

Classes:
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
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate, scipy.ndimage
import matplotlib

import astropy.io.fits as fits
import astropy.io.ascii as ioascii
from astropy.config import ConfigurationItem, get_config_dir, save_config


import poppy

from . import conf


try: 
    import pysynphot
    _HAS_PYSYNPHOT = True
except:
    _HAS_PYSYNPHOT = False
 

import logging
_log = logging.getLogger('webbpsf')



class JWInstrument(poppy.instrument.Instrument):
    """ A generic JWST Instrument class.

    *Note*: Do not use this class directly - instead use one of the specific instrument subclasses!

    This class provides a simple interface for modeling PSF formation through the JWST instruments, 
    with configuration options and software interface loosely resembling the configuration of the instrument 
    hardware mechanisms.   
    
    This module currently only provides a modicum of error checking, and relies on the user
    being knowledgable enough to avoid trying to simulate some physically impossible or just plain silly
    configuration (such as trying to use a FQPM with the wrong filter).

    The instrument constructors do not take any arguments. Instead, create an instrument object and then
    configure the `filter` or other attributes as desired. The most commonly accessed parameters are 
    available as object attributes: `filter`, `image_mask`, `pupil_mask`, `pupilopd`. More advanced
    configuration can be done by editing the :ref:`JWInstrument.options` dictionary, either by passing options to __init__ or by directly editing the dict afterwards.
    """

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

    def __init__(self, name="", pixelscale = 0.064):
        self.name=name
        self.pixelscale = pixelscale

        self._WebbPSF_basepath = conf.get_webbpsf_data_path()

        self._datapath = self._WebbPSF_basepath + os.sep + self.name + os.sep
        self._filter = None
        self._image_mask = None
        self._pupil_mask = None
        self.pupil = os.path.abspath(self._datapath+"../pupil_RevV.fits")
        "Filename *or* fits.HDUList for JWST pupil mask. Usually there is no need to change this."
        self.pupilopd = None   # This can optionally be set to a tuple indicating (filename, slice in datacube)
        """Filename *or* fits.HDUList for JWST pupil OPD. 
        
        This can be either a full absolute filename, or a relative name in which case it is
        assumed to be within the instrument's `data/OPDs/` directory, or an actual fits.HDUList object corresponding to such a file.
        If the file contains a datacube, you may set this to a tuple (filename, slice) to select a given slice, or else
        the first slice will be used."""


        self.options = {} # dict for storing other arbitrary options. 

        #create private instance variables. These will be
        # wrapped just below to create properties with validation.
        self._filter=None

        filter_table = ioascii.read(self._WebbPSF_basepath + os.sep+ 'filters.txt')
        wmatch = np.where(filter_table['instrument'] == self.name)
        self.filter_list = filter_table['filter'][wmatch].tolist()
        "List of available filters"
        self._filter_nlambda_default = dict(zip(filter_table['filter'][wmatch], filter_table['nlambda'][wmatch]))

        #self._filter_files= [os.path.abspath(f) for f in glob.glob(self._datapath+os.sep+'filters/*_thru.fits')]
        #self.filter_list=[os.path.basename(f).split("_")[0] for f in self._filter_files]
        if len(self.filter_list) ==0: 
            #self.filter_list=[''] # don't crash for FGS which lacks filters in the usual sense
            raise ValueError("No filters available!")

        def sort_filters(filtname):
            try:
                if name =='MIRI': return int(filtname[1:-1]) # MIRI filters have variable length number parts
                else: return int(filtname[1:4]) # the rest do not, but have variable numbers of trailing characters
            except:
                return filtname
        self.filter_list.sort(key=sort_filters)
        self._filter_files = [self._datapath+os.sep+'filters/'+f+"_throughput.fits" for f in self.filter_list]

        self.filter = self.filter_list[0]
        self._rotation = None


        #self.opd_list = [os.path.basename(os.path.abspath(f)) for f in glob.glob(self._datapath+os.sep+'OPD/*.fits')]
        self.opd_list = [os.path.basename(os.path.abspath(f)) for f in glob.glob(self._datapath+os.sep+'OPD/OPD*.fits')]
        if len(self.opd_list) > 0:
            self.pupilopd = self.opd_list[-1]
            #self.pupilopd = self.opd_list[len(self.opd_list)/2]
        #self.opd_list.insert(0,"Zero OPD (Perfect)")

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

        self.detector_coordinates = (0,0) # where is the source on the detector, in 'Science frame' pixels?

    def _validate_config(self):
        pass

    # create properties with error checking
#    @property
#    def filter(self):
#        'Currently selected filter name (e.g. "F200W")'
#        return self._filter
#    @filter.setter
#    def filter(self, value):
#        value = value.upper() # force to uppercase
#        if value not in self.filter_list:
#            raise ValueError("Instrument %s doesn't have a filter called %s." % (self.name, value))
#        self._filter = value
#        self._validate_config()

    @property
    def image_mask(self):
        'Currently selected image plane mask, or None for direct imaging'
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
                    raise ValueError("Instrument %s doesn't have an image mask called %s." % (self.name, name))
        self._image_mask = name

    @property
    def pupil_mask(self):
        'Currently selected Lyot pupil mask, or None for direct imaging'
        return self._pupil_mask
    @pupil_mask.setter
    def pupil_mask(self,name):
        if name is "": name = None
        if name is not None:
            if name in self.pupil_mask_list:
                pass # there's a perfect match, this is fine.
            else:
                name = name.upper() # force to uppercase
                if name not in self.pupil_mask_list:
                    raise ValueError("Instrument %s doesn't have an pupil mask called %s." % (self.name, name))

        self._pupil_mask = name

 
    @property
    def detector(self):
        'Currently selected detector name (for instruments with multiple detectors)'
        return self._detector.name
    @detector.setter
    def detector(self,detname):
        #if name is "": name = None
        if detname is not None:
            detname = detname.upper() # force to uppercase
            try:
                siaf_aperture_name = self._detector2siaf[detname]
            except:
                raise ValueError("Unknown name: {0} is not a valid known name for a detector for instrument {1}".format(detname, self.name))
            self._detector = DetectorGeometry(self.name, siaf_aperture_name, shortname=detname)

    def __str__(self):
        return "JWInstrument name="+self.name

    #----- actual optical calculations follow here -----
    def calcPSF(self, outfile=None, source=None, filter=None,  nlambda=None, monochromatic=None ,
            fov_arcsec=None, fov_pixels=None,  oversample=None, detector_oversample=None, fft_oversample=None, calc_oversample=None, rebin=True,
            clobber=True, display=False, return_intermediates=False, **kwargs):
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

        if calc_oversample is not None: fft_oversample = calc_oversample # back compatibility hook for deprecated arg name.

        _log.info("Setting up PSF calculation for "+self.name)

        # first make sure that webbpsf's configuration is used to override any of the
        # same configuration options in poppy. This is admittedly perhaps overbuilt to have identical
        # settings in both packages, but the intent is to shield typical users of webbpsf
        # from having to think about the existence of the underlying library. They can 
        # just deal with one set of settings.
        conf._apply_settings_to_poppy()

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
                _log.warn("unrecognized filter %s. setting default nlambda=%d" % (self.filter, nlambda))
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
            oversample = conf.default_oversampling()
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
        output_mode = options.get('output_mode',conf.default_output_mode())

        if output_mode == 'Mock JWST DMS Output':
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
        optsys = poppy.OpticalSystem(name='JWST+'+self.name, oversample=fft_oversample)
        if 'source_offset_r' in options.keys(): optsys.source_offset_r = options['source_offset_r']
        if 'source_offset_theta' in options.keys(): optsys.source_offset_theta = options['source_offset_theta']


        #---- set pupil OPD
        if isinstance(self.pupilopd, str):  # simple filename
            full_opd_path = self.pupilopd if os.path.exists( self.pupilopd) else os.path.join(self._datapath, "OPD",self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], basestring): # tuple with filename and slice
            full_opd_path =  (self.pupilopd[0] if os.path.exists( self.pupilopd[0]) else os.path.join(self._datapath, "OPD",self.pupilopd[0]), self.pupilopd[1])
        elif isinstance(self.pupilopd, fits.HDUList) or isinstance(self.pupilopd, poppy.OpticalElement): # OPD supplied as FITS object
            full_opd_path = self.pupilopd # not a path per se but this works correctly to pass it to poppy
        elif self.pupilopd is None: 
            full_opd_path = None
        else:
            raise TypeError("Not sure what to do with a pupilopd of that type:"+str(type(self.pupilopd)))

        #---- set pupil intensity
        if isinstance(self.pupil, str): # simple filename
            full_pupil_path = self.pupil if os.path.exists( self.pupil) else os.path.join(self._WebbPSF_basepath,self.pupil)
        elif isinstance(self.pupil, fits.HDUList): # pupil supplied as FITS object
            full_pupil_path = self.pupil
        else: 
            raise TypeError("Not sure what to do with a pupil of that type:"+str(type(self.pupil)))


        #---- apply pupil intensity and OPD to the optical model
        optsys.addPupil(name='JWST Pupil', transmission=full_pupil_path, opd=full_opd_path, opdunits='micron', rotation=self._rotation)

        #---- Add defocus if requested
        if 'defocus_waves' in options.keys(): 
           defocus_waves = options['defocus_waves'] 
           defocus_wavelength = float(options['defocus_wavelength']) if 'defocus_wavelength' in options.keys() else 2.0e-6
           _log.info("Adding defocus of %d waves at %.2f microns" % (defocus_waves, defocus_wavelength *1e6))
           lens = poppy.ThinLens(name='Defocus', nwaves=defocus_waves, reference_wavelength=defocus_wavelength)
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
        """Add instrument-internal optics to an optical system, typically coronagraphic or spectrographic in nature. 
        This method must be provided by derived instrument classes. 

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

        By subclassing this, you can define whatever custom bandpasses are appropriate for your instrument

        """
        try:
            band = pysynphot.ObsBandpass( ('%s,im,%s'%(self.name, filtername)).lower())
        except:
            _log.warn("Filter %s not supported in available pysynphot/CDBS. Falling back to local filter transmission files" % filtername)
            _log.warn("These may be less accurate.")

            # the requested band is not yet supported in synphot/CDBS. (those files are still a
            # work in progress...). Therefore, use our local throughput files and create a synphot
            # transmission object.
            wf = np.where(np.asarray(self.filter_list)== filtername)[0]

            if len(wf) != 1:
                _log.error("Could not find a match for filter name = %s in the filters for %s." % (filtername, self.name))
            else:
                wf = wf[0]
            # The existing FITS files all have wavelength in ANGSTROMS since that is the pysynphot convention...
            filterfits = fits.open(self._filter_files[wf])
            filterdata = filterfits[1].data 
            try:
                f1 = filterdata.WAVELENGTH
                d2 = filterdata.THROUGHPUT
            except:
                raise ValueError("The supplied file, %s, does not appear to be a FITS table with WAVELENGTH and THROUGHPUT columns." % self._filter_files[wf] )
            try:
                if filterfits[1].header['WAVEUNIT'] != 'Angstrom': raise ValueError("The supplied file, %s, does not have WAVEUNIT = Angstrom as expected." % self._filter_files[wf] )
            except:
                _log.warn('The supplied file, %s, does not have a WAVEUNIT keyword. Assuming it is Angstroms.' %  self._filter_files[wf])
 
            band = pysynphot.spectrum.ArraySpectralElement(throughput=filterdata.THROUGHPUT,
                                wave=filterdata.WAVELENGTH, waveunits='angstrom',name=filtername)
        return band



#######  JWInstrument classes  #####

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
        self.monochromatic= 8.0
        self._IFU_pixelscale = {'Ch1':(0.18, 0.19), 'Ch2':(0.28, 0.19), 'Ch3': (0.39, 0.24), 'Ch4': (0.64, 0.27) }
            # The above tuples give the pixel resolution (perpendicular to the slice, along the slice). 
            # The pixels are not square.

        #self._default_aperture='MIRIM_center' # reference into SIAF for ITM simulation V/O coords
        self.detector_list = ['MIRIM']
        self._detector2siaf = {'MIRIM':'MIRIM_FULL_ILLCNTR'}
        self.detector=self.detector_list[0]

    def _validate_config(self):
        #_log.debug("MIRI validating:    %s, %s, %s " % (self.filter, self.image_mask, self.pupil_mask))
        if self.filter.startswith("MRS-IFU"): raise NotImplementedError("The MIRI MRS is not yet implemented.")

        #_log.warn("MIRI config validation disabled for now - TBD rewrite ")
        return

        if self.image_mask is not None or self.pupil_mask is not None:
            if self.filter == 'F1065C':
                assert self.image_mask == 'FQPM1065', 'Invalid configuration'
                assert self.pupil_mask == 'MASKFQPM', 'Invalid configuration'
            elif self.filter == 'F1140C':
                assert self.image_mask == 'FQPM1140', 'Invalid configuration'
                assert self.pupil_mask == 'MASKFQPM', 'Invalid configuration'
            elif self.filter == 'F1550C':
                assert self.image_mask == 'FQPM1550', 'Invalid configuration'
                assert self.pupil_mask == 'MASKFQPM', 'Invalid configuration'
            elif self.filter == 'F2300C':
                assert self.image_mask == 'LYOT2300', 'Invalid configuration'
                assert self.pupil_mask == 'MASKLYOT', 'Invalid configuration'
            else:
                #raise ValueError("Invalid configuration selected!")
                _log.warn("*"*80)
                _log.warn("WARNING: you appear to have selected an invalid/nonphysical configuration of that instrument!")
                _log.warn("")
                _log.warn("I'm going to continue trying the calculation, but YOU are responsible for interpreting")
                _log.warn("any results in a meaningful fashion or discarding them..")
                _log.warn("*"*80)


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
        if (self.image_mask is not None and 'FQPM' in self.image_mask) or 'force_fqpm_shift' in self.options.keys() : optsys.addPupil("FQPM_FFT_aligner")

        if self.image_mask == 'FQPM1065':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1065",
                opticslist = [  poppy.IdealFQPM(wavelength=10.65e-6, name=self.image_mask),
                                poppy.IdealFieldStop(size=24, angle=-self._rotation)])
            optsys.addImage(container)
            trySAM = False
        elif self.image_mask == 'FQPM1140':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1140",
                opticslist = [  poppy.IdealFQPM(wavelength=11.40e-6, name=self.image_mask),
                                poppy.IdealFieldStop(size=24, angle=-self._rotation)])
            optsys.addImage(container)
            trySAM = False
        elif self.image_mask == 'FQPM1550':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1550",
                opticslist = [  poppy.IdealFQPM(wavelength=15.50e-6, name=self.image_mask),
                                poppy.IdealFieldStop(size=24, angle=-self._rotation)])
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
                opticslist = [poppy.IdealCircularOcculter(radius =4.25/2, name=self.image_mask),
                              poppy.IdealBarOcculter(width=0.722), 
                              poppy.IdealFieldStop(size=30, angle=-self._rotation)] )
            optsys.addImage(container)
            trySAM = True
            SAM_box_size = [5,20]
        elif self.image_mask == 'LRS slit':
            # one slit, 5.5 x 0.6 arcsec in height (nominal)
            #           4.7 x 0.51 arcsec (measured for flight model. See MIRI-TR-00001-CEA)
            # 
            # Per Klaus Pontoppidan: The LRS slit is aligned with the detector x-axis, so that the dispersion direction is along the y-axis. 
            optsys.addImage(optic=poppy.IdealRectangularFieldStop(width=5.5, height=0.6, angle=self._rotation, name= self.image_mask))
            trySAM = False
        else:
            optsys.addImage()
            trySAM = False

        if (self.image_mask is not None and 'FQPM' in self.image_mask)  or 'force_fqpm_shift' in self.options.keys() : optsys.addPupil("FQPM_FFT_aligner", direction='backward')

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
        """For NIRCam, this checks whenever you change a filter and updates the pixelscale appropriately"""
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
            optsys.addImage(function='BandLimitedCoron', kind='nircamcircular', sigma=5.253 , name=self.image_mask)
            trySAM = True
            SAM_box_size = 5.0
        elif self.image_mask == 'MASK335R':
            optsys.addImage(function='BandLimitedCoron', kind='nircamcircular', sigma=3.2927866 , name=self.image_mask)
            trySAM = True
            SAM_box_size = 5.0
        elif self.image_mask == 'MASK430R':
            optsys.addImage(function='BandLimitedCoron', kind='nircamcircular', sigma=2.588496*0.99993495 , name=self.image_mask)
            trySAM = True
            SAM_box_size = 5.0
        elif self.image_mask == 'MASKSWB':
            optsys.addImage(function='BandLimitedCoron', kind='nircamwedge', wavelength=2.1e-6, name=self.image_mask)
            trySAM = False #True FIXME
            SAM_box_size = [5,20]
        elif self.image_mask == 'MASKLWB':
            optsys.addImage(function='BandLimitedCoron', kind='nircamwedge', wavelength=4.6e-6, name=self.image_mask)
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


        #optsys.addPupil( name='null for debugging NIRcam _addCoron') # debugging
        if self.pupil_mask == 'CIRCLYOT':
            optsys.addPupil(transmission=self._datapath+"/optics/NIRCam_Lyot_Somb.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'WEDGELYOT':
            optsys.addPupil(transmission=self._datapath+"/optics/NIRCam_Lyot_Sinc.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'WEAK LENS +4':
            optsys.addPupil(poppy.ThinLens(name='Weak Lens +4', nwaves=4, reference_wavelength=2e-6))
        elif self.pupil_mask == 'WEAK LENS +8':
            optsys.addPupil(poppy.ThinLens(name='Weak Lens +8', nwaves=8, reference_wavelength=2e-6))
        elif self.pupil_mask == 'WEAK LENS -8':
            optsys.addPupil(poppy.ThinLens(name='Weak Lens -8', nwaves=-8, reference_wavelength=2e-6))
        elif self.pupil_mask == 'WEAK LENS +12 (=4+8)':
            stack = poppy.CompoundAnalyticOptic(name='Weak Lens Stack +12', opticslist=[
                poppy.ThinLens(name='Weak Lens +4', nwaves=4, reference_wavelength=2e-6),
                poppy.ThinLens(name='Weak Lens +8', nwaves=8, reference_wavelength=2e-6)])
            optsys.addPupil(stack)
        elif self.pupil_mask == 'WEAK LENS -4 (=4-8)':
            stack = poppy.CompoundAnalyticOptic(name='Weak Lens Stack -4', opticslist=[
                poppy.ThinLens(name='Weak Lens +4', nwaves=4, reference_wavelength=2e-6),
                poppy.ThinLens(name='Weak Lens -8', nwaves=-8, reference_wavelength=2e-6)])
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
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('NIRSpec does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
    def _addAdditionalOptics(self,optsys, oversample=2):
        """ Add fixed slit optics for NIRSpec

        See Table 3-6 of NIRSpec Ops Concept Document, ESA-JWST-TN-0297 / JWST-OPS-003212

        """
        trySAM = False # semi-analytic method never applicable here. 
        SAM_box_size = None
        if self.image_mask == 'S200A1' or self.image_mask == 'S200A2' or self.image_mask == 'S200B1':
            # three identical slits, 0.2 x 3.2 arcsec in length
            optsys.addImage(optic=poppy.IdealRectangularFieldStop(width=0.2, height=3.2, name= self.image_mask + " slit"))
        elif self.image_mask == 'S400A1':
            # one slit, 0.4 x 3.65 arcsec in height
            optsys.addImage(optic=poppy.IdealRectangularFieldStop(width=0.4, height=3.65, name= self.image_mask + " slit"))
        elif self.image_mask == 'S1600A1':
            # square aperture for exoplanet spectroscopy
            optsys.addImage(optic=poppy.IdealRectangularFieldStop(width=1.6, height=1.6, name= self.image_mask + " square aperture"))
        elif self.image_mask == 'MSA all open':
            # all MSA shutters open 
            optsys.addImage(optic=NIRSpec_MSA_open_grid(name= self.image_mask))
        elif self.image_mask == 'Single MSA open shutter':
            # one MSA open shutter aperture 
            optsys.addImage(optic=poppy.IdealRectangularFieldStop(width=0.2, height=0.45, name= self.image_mask))
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

    def _validate_config(self):
        if self.filter.startswith("IFU"): raise NotImplementedError("The NIRSpec IFU is not yet implemented.")


    def _getFITSHeader(self, hdulist, options):
        """ Format NIRSpec-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)
        hdulist[0].header['GRATING'] = ( 'None', 'NIRSpec grating element name')
        hdulist[0].header['APERTURE'] = ( str(self.image_mask), 'NIRSpec slit aperture name')


class NIRISS(JWInstrument):
    """ A class modeling the optics of the Near-IR Imager and Slit Spectrograph
        (formerly TFI)
    
    Relevant attributes include `image_mask`, and `pupil_mask`.

    WebbPSF models the direct imaging and nonredundant aperture masking modes of NIRISS in the usual manner. 

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

    """
    def __init__(self):
        JWInstrument.__init__(self, "NIRISS", pixelscale=0.064)

        self.image_mask_list = ['CORON058', 'CORON075','CORON150','CORON200'] # available but unlikely to be used...
        self.pupil_mask_list = ['MASK_NRM','CLEAR', 'GR700XD']

        self._detector2siaf = {'NIRISS':'NIS_FULL_CNTR'}
        self.detector_list = ['NIRISS']
        self.detector=self.detector_list[0]


        self._default_aperture='NIRISS center' # reference into SIAF for ITM simulation V/O coords


    def _validate_config(self):
        pass

    def _addAdditionalOptics(self,optsys, oversample=2):
        """Add NRM or slitless spectroscopy optics for NIRISS. 

            These are probably not going to be used much in practice for NIRISS, but they
            are present, so we might as well still provide the ability to simulate 'em. 
        """
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
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASK_NRM.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'CLEAR':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKCLEAR.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'GR700XD':
            #optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKSOSS.fits.gz", name=self.pupil_mask, shift=shift)
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
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('FGS does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
        #TODO only one possible filter fot the FGS, too. 
    def _addAdditionalOptics(self,optsys):
        raise NotImplementedError("No Coronagraph in FGS!")

    def _getFITSHeader(self, hdulist, options):
        """ Format FGS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)
        hdulist[0].header['FOCUSPOS'] = (0,'FGS focus mechanism not yet modeled.')


#######  Custom Optics used in JWInstrument classes  #####


class NIRSpec_three_MSA_shutters(poppy.AnalyticOpticalElement):
    """ Three NIRSpec MSA shutters, adjacent vertically."""

    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the field stop.

        The area of an open shutter is 0.2 x 0.45, while the shutter pitch is 0.26x0.51
        The walls separating adjaced shutters are 0.06 arcsec wide.
        """

        msa_width = 0.2
        msa_height = 0.45
        msa_wall = 0.06

        if not isinstance(wave, poppy.Wavefront):
            raise ValueError("IdealFieldStop getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._IMAGE)

        y, x= wave.coordinates()
        #xnew =  x*np.cos(np.deg2rad(self.angle)) + y*np.sin(np.deg2rad(self.angle))
        #ynew = -x*np.sin(np.deg2rad(self.angle)) + y*np.cos(np.deg2rad(self.angle))
        #x,y = xnew, ynew


        self.transmission = np.zeros(wave.shape)
        # get the innermost shutter than spans the Y axis
        w_inside_1 = np.where( (abs(y) < (msa_height/2))  & (abs(x) < (msa_width/2)))
        self.transmission[w_inside_1] = 1
        # get the adjacent shutters one above and one below.
        w_inside_2 = np.where( (abs(y) > (msa_height/2)+msa_wall) & (abs(y) < msa_height*1.5+msa_wall)  & (abs(x) < (msa_width/2)))
        self.transmission[w_inside_2] = 1

        return self.transmission


class NIRSpec_MSA_open_grid(poppy.AnalyticOpticalElement):
    """ An infinite repeating region of the NIRSpec MSA grid"""

    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the field stop.

        The area of an open shutter is 0.2 x 0.45, while the shutter pitch is 0.26x0.51
        The walls separating adjaced shutters are 0.06 arcsec wide.
        """

        msa_width = 0.2
        msa_height = 0.45
        msa_wall = 0.06
        msa_x_pitch = 0.26
        msa_y_pitch = 0.51

        if not isinstance(wave, poppy.Wavefront):
            raise ValueError("IdealFieldStop getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._IMAGE)

        y, x= wave.coordinates()
        #xnew =  x*np.cos(np.deg2rad(self.angle)) + y*np.sin(np.deg2rad(self.angle))
        #ynew = -x*np.sin(np.deg2rad(self.angle)) + y*np.cos(np.deg2rad(self.angle))
        #x,y = xnew, ynew

        mask_vert_walls  = np.abs(np.mod(np.abs(x), msa_x_pitch) - (msa_x_pitch/2)) < msa_wall/2
        mask_horz_walls  = np.abs(np.mod(np.abs(y), msa_y_pitch) - (msa_y_pitch/2)) < msa_wall/2



        self.transmission = np.ones(wave.shape)
        self.transmission[mask_vert_walls] = 0
        self.transmission[mask_horz_walls] = 0

        return self.transmission


class NIRISS_GR700XD_Grism(poppy.FITSOpticalElement):
    """ Custom optic class to model the NIRISS SOSS grim GR700XD

    This includes both the pupil mask file and the cylindrical lens

    Based on inputs from Loic Albert, Anand Sivaramakrishnan, and Andre Martel
    In particular see FGS_TFI_UdM_035_RevD for details of the NIRISS GR700XD
    measurement, and JWST-STScI-003338 for detector orientation and layout.

    GRISM DESIGN:

    The grism (and cylinder) are per design rotated by 2 degrees so as to be able
    to sample an emission line across different pixel position along the spatial
    direction (kind of resampling the line and not be limited by intra pixel
    response).  

    From Loic Albert's NIRISS technical report:

        * surface sag for the cylinder: 3.994 micron peak
        * limited to 3.968 microns for the 26 mm FOV mask

    From Loic Albert's email to Marshall 2013-07-18:

            I do have an update concerning the geometry of the GR700XD pupil
            mask. It turns out that they clocked the grism by about 2.25 degrees wrt the
            OTE system of coordinates. However, the field mask did not follow and is still
            aligned along the OTE s.o.c. That was a mistake that fortunately does have much
            impact.

            Comdev is in the process of modelling a new mask for the
            Spare grism. Remember that we will swap the current FLight GR700XD for
            its Spare which offers much improved performances. The mask size will
            be a little different (rectangular) and this time will be clocked 2.25
            degrees along with the grism.

            The sign of the rotation of the grism will have to be
            devised by trying the 2 possibilities an looking at the resulting tilt
            of the monochromatic PSF and the position of that PSF on the detector.
            Attached is a simulation of what we expect based on my own PSF
            generator.

            The cylinder lens has a well characterized power (actually radius of curvature). The values are:
                current Flight: 22.85 meters
                Spare: 22.39 meters 

            Prism physical size: pupil is 26 mm on a side for the current prism, will be 28 mm for the spare

    From Loic Albert's email to Marshall 2013-09-19:

        The latest news on this front are:

        1 - The current Flight mask is attached. It is 26x26 mm. The mask and grism are
            *not* aligned along the same coordinates. That was a mistake. I'll forward you
            a message from Michael M., our optics expert at CSA. 

        2 - The Spare mask (likely the one which will fly) is not built yet. The mask
            will be aligned along the grism coordinate and both will be clocked 2.2 deg wrt
            the OTE.

        3 - A ghost analysis showed that the current grism clocking will suffer from
            large ghosts. So we are studying how to clock the Spare grism in its cell to
            minimize ghosts. Likely a 90 degrees rotation will be applied to baseline of
            point 2.

        From Michael.Maszkiewicz@asc-csa.gc.ca:

            As far as I understand now, we have two rotations in the as-built
            GR700. One rotation is for the prism-grism combo by 2 deg CCW, looking along
            the local +z axis, and the second rotation is for  the mask by 3.05 deg  but
            CW. As a result there is total 5.05 deg rotation between grism and its mask.
            See my annotations to your drawing attached.
   
     ORIENTATION:

        See Figure 2 of JWST-STScI-003338
        In "DMS" coordinates, as projected looking outwards onto the sky,
        The GR700XD grating trace is near the extreme right edge of the detector
        with long wavelengths closest to (2048,2048) and short wavelengths nearest (2048,0)
        (The raw detector coordinates are very different from this due to a 180 degree rotation)

        **PLEASE NOTE** that the DMS when processing spectral data performs an additional transformation:
            For spectral data, the science X-axis is aligned with the detector
            dispersion direction and the science frame Y-axis is at a right angle
            to the X-axis in a right-handed coordinate system (Swade 2003)

        We choose here to ignore that complication; WebbPSF simulates the 2D sky projected
        image in "Sci" coordinates in the terminology for SIAF from Lallo et al. 
        In this coordinate system, the dispersion from the cylinder lens is aligned 
        almost along V2 and the longer wavelengths are oriented toward +V3. 
        

   

    Parameters
    ----------
    which : string
        'flight' or 'spare'. Properties are hard coded. 
    """
    #
    #    transmission : string filename
    #        file for the pupil transmission function
    #    cylinder_sag_mm : float
    #        physical thickness of the cylindrical lens, in millimeters
    #    rotation_angle : float
    #        degrees clockwise for the orientation of the cylinder's dispersing axis. Default
    #        of 92.25 should be consistent with current NIRISS flight and spare, except for
    #        sign ambiguity.
    #    rotate_mask : bool
    #        should the field mask be rotated along with the cylinder? False for first gen flight
    #        prism, true for expected spare replacement.

    def __init__(self, name='GR700XD', which='flight',
            #cylinder_radius=22.85,  cylinder_sag_mm=4.0, rotation_angle=92.25, rotate_mask=False, transmission=None, 
            shift=None):
        # Initialize the base optical element with the pupil transmission and zero OPD

        

        if which=='spare':
            raise NotImplementedError("Rotated field mask for spare grism not yet implemented!")
        else:
            transmission=os.path.join( conf.get_webbpsf_data_path(), "NIRISS/optics/MASKGR700XD.fits.gz")

        self.shift=shift
        poppy.FITSOpticalElement.__init__(self, name=name, transmission=transmission, planetype=poppy.poppy_core._PUPIL, shift=shift)

        # UPDATED NUMBERS 2013-07:
        # See Document FGS_TFI_UdM_035_RevD

        if which =='flight':
            # 3.994 microns P-V over 27.02 mm measured (Loic's email)
            # This is **surface sag**, corresponding to P-V of 6.311 waves at lambda=632.8 nm.
            # should correspond to 3.698 microns over 26 mm clear aperture. 
            self.cylinder_radius = 22.85 # radius of curvature
            self.prism_size = 0.02702 # 27.02 millimeters for the physical prism
            self.prism_clear_aperture = 0.0260 # 26 mm clear aperture for the prism + mount
            self.cylinder_rotation_angle = 2.25

            # pupil magnification computed from 22 mm clear aperture reported = 
            # 857-169 pixels = 699 pixels in the 2D array which has scale =.00645604
            # = 4.44175 meters projected on the primary

            # therefore the magnification is 0.1708 meters projected on the primary / mm in the NIRISS pupil
            self.pupil_demagnification =  170.8367 # meters on the primary / meters in the NIRISS pupil
        else:
            # 5.8 microns P-V over 32.15 mm (Loic's email)
            # should correspond to 4.38 microns over 28 mm clear aperture
            self.cylinder_radius = 22.39 # radius of curvature
            self.prism_size = 0.03215 # millimeters for the physical prism
            self.prism_clear_aperture = 0.0280 # clear aperture for the prism + mount
            self.cylinder_rotation_angle = 2.25

        # initial population of the OPD array for display etc.
        self.makeCylinder( 2.0e-6) 
    def makeCylinder(self, wave):
        if isinstance(wave, poppy.Wavefront):
            wavelength=wave.wavelength
        else:
            wavelength=wave

        # compute indices in pixels, relative to center of plane, with rotation
        # units of these are meters
        y, x = np.indices(self.opd.shape, dtype=float)
        y-= (self.opd.shape[0]-1)/2.
        x-= (self.opd.shape[1]-1)/2.  
 
        ang = np.deg2rad(self.cylinder_rotation_angle )
        x = np.cos(ang)*x - np.sin(ang)*y
        y = np.sin(ang)*x + np.cos(ang)*y


        

        # From IDL code by David Lafreniere:
        #  ;the cylindrical defocus
        #x=(dindgen(pupdim)-pupdim/2)#replicate(1,pupdim)
        #y0=(rpuppix^2+sag[s]^2)/(2*sag[s])
        #wfe1=y0-sqrt(y0^2-x^2)
        #if sag[s] lt 1.e-5 then wfe1=0.d0

        # Here I will just translate that to Python exactly, making use of the
        # variables here:

        # rpuppix = radius of pupil in pixels
        #rpuppix = self.amplitude_header['DIAM'] / self.amplitude_header['PUPLSCAL'] / 2
        # Calculate the radius of curvature of the cylinder, bsaed on 
        # the chord length and height 

        # In this case we're assuming the cylinder is precisely as wide as the projected
        # telescope pupil. This doesn't seem guaranteed:
        #  * actual chord length across cylinder: 27.02 mm. 
        #  * projected primary scale at NIRISS = ?



        # Compute the overall sag of the cylinder lens at its outer edge. This is not actually used, it's
        # just for cross-check of the values
        # the sag will depend on half the pupil size since that's the offset from center to edge
        sag0 = np.sqrt(self.cylinder_radius**2 - (self.prism_size/2)**2) - self.cylinder_radius
        _log.debug(" Computed GR700XD cylinder sag: {0:.3g} meters".format(sag0))

        # now compute the spatially dependent sag of the cylinder, as projected onto the primary
        sag = np.sqrt(self.cylinder_radius**2 - (x*self.amplitude_header['PUPLSCAL']/self.pupil_demagnification)**2) - self.cylinder_radius


        # what we really want to do is take the physical properties of the as-built optic, and interpolate into that
        # to compute the OPD after remapping based on the pupil scale (and distortion?)
        #y0=(rpuppix**2+self.cylinder_sag**2)/(2*self.cylinder_sag)
        
        #wfe1=y0-np.sqrt(y0**2-x**2)

        _log.debug(" Cylinder P-V: {0:.4g} meters physical sag across full array".format(sag.max()-sag.min()) )


        # remove piston offset 
        wnz = np.where(self.amplitude != 0)
        sag -= sag[wnz].min()   # normalize to 0 at the minimum
        #sag -= sag[wnz].mean()    # normalize around the mean
        sag[self.amplitude == 0] = 0 # no OPD in opaque regions (makes no difference in propagation but improves display)
        _log.debug(" Cylinder P-V: {0:.4g} meters physical sag across clear aperture".format(sag[wnz].max()-sag[wnz].min()) )

        # scale for ZnSe index of refraction, 
        self.opd = sag *  (self.ZnSe_index(wavelength) -1)
        _log.debug(" Cylinder P-V: {0:.4g} meters optical sag at {1:.3g} microns across clear aperture".format(self.opd[wnz].max()-self.opd[wnz].min(), wavelength*1e6) )

        #stop()

    def ZnSe_index(self, wavelength):
        """ Return cryogenic index of refraction of ZnSe at an arbitrary wavelength
        """
        # From ZnSe_index.txt provided by Loic Albert
        #from Michael M. Nov 9 2012 in excel table],
        #ZnSe-40K index,  ],
        # ZnSe
        ZnSe_data =np.asarray([[500,  2.7013],
                                [540,  2.6508],
                                [600,  2.599],
                                [644,  2.56937],
                                [688,  2.54709],
                                [732,  2.52977],
                                [776,  2.51596],
                                [820,  2.50472],
                                [864,  2.49542],
                                [900,  2.4876],
                                [908,  2.48763],
                                [952,  2.48103],
                                [996,  2.47537],
                                [1040,  2.47048],
                                [1084,  2.46622],
                                [1128,  2.46249],
                                [1172,  2.4592],
                                [1216,  2.45628],
                                [1260,  2.45368],
                                [1304,  2.45134],
                                [1348,  2.44924],
                                [1392,  2.44734],
                                [1436,  2.44561],
                                [1480,  2.44405],
                                [1524,  2.44261],
                                [1568,  2.4413],
                                [1612,  2.44009],
                                [1656,  2.43897],
                                [1700,  2.43794],
                                [1744,  2.43699],
                                [1788,  2.4361],
                                [1832,  2.43527],
                                [1876,  2.4345],
                                [1920,  2.43378],
                                [1964,  2.4331],
                                [2008,  2.43247],
                                [2052,  2.43187],
                                [2096,  2.4313],
                                [2140,  2.43077],
                                [2184,  2.43026],
                                [2228,  2.42978],
                                [2272,  2.42933],
                                [2316,  2.4289],
                                [2360,  2.42848],
                                [2404,  2.42809],
                                [2448,  2.42771],
                                [2492,  2.42735],
                                [2536,  2.42701],
                                [2580,  2.42667],
                                [2624,  2.42635],
                                [2668,  2.42604],
                                [2712,  2.42575],
                                [2756,  2.42546],
                                [2800,  2.42518],
                                [2844,  2.42491],
                                [2888,  2.42465],
                                [2932,  2.4244],
                                [2976,  2.42416],
                                [3020,  2.42392]] )

        interpol_znse = scipy.interpolate.interp1d( ZnSe_data[:,0]*1e-9, ZnSe_data[:,1] )
        return interpol_znse(wavelength)

    def getPhasor(self, wave):
        """ Scale the cylindrical lens OPD appropriately for the current wavelength
            Then call the regular getphasor method of the parent class 

        """
        self.makeCylinder(wave)
        return poppy.FITSOpticalElement.getPhasor(self, wave)

    def display(self, opd_vmax=6e-6, *args, **kwargs):
        "Same as regular display for any other optical element, except opd_vmax default changed"
        poppy.FITSOpticalElement.display(self,*args, opd_vmax=opd_vmax, **kwargs)



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
Instrument.list = ['nircam', 'nirspec', 'tfi', 'miri'] # useful list for iteration


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


def MakePSF(self, instrument=None, pupil_file=None, phase_file=None, output=None,
                  diameter=None, oversample=4, type=np.float64,
                  filter=((1.,),(1.,)),
                  output_size=512, pixel_size=None, verbose=False):
    """This is a wrapper function to provide back-compatibility with the
    interface of the original JWPSF. New code should make use of the
    object-oriented interface instead. 
    """
    instr = Instrument(name)
    raise NotImplementedError("finish me!")

    if diameter is not None or float(diameter) != 6.5:
        raise NotImplementedError("Changing diameters is not supported by this version of JWPSF. Nor by JWST itself, for that matter!")

    if pixel_size is not None:
        instr.pixelscale = pixel_size
    # set pupil files
    instr.pupil = pupil_file
    instr.pupilopd = phase_file

    # deal with filter? 
    # pick the mean wavelen and use the closest filter to that?

    return instr.calcPSF(oversample=oversample, fov_pixels=output_size)



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

        self.mysiaf = SIAF(instr=self.instrname, basepath=os.path.join( conf.get_webbpsf_data_path(), self.instrname) )
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


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')

#
#    nc = NIRCam()
#    nc.filter = 'F460M'
#    nc.image_mask = 'MASK430R'
#    nc.pupil_mask = 'CIRCLYOT'
#    #nc.calcPSF('test_nircam.fits', mono=False)
#
    miri=MIRI()
    miri.image_mask = 'LYOT2300'
    miri.pupil_mask = 'MASKLYOT'
    miri.filter='F2300C'
    plt.clf()
    miri.display()
#
#    #miri.display()
#    nircam=NIRCam()
#    tfi = TFI()
#    tfi.image_mask = "CORON058"
#    tfi.pupil_mask = 'MASKC66N'
#    nirspec = NIRSpec()
