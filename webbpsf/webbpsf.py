#!/usr/bin/env python
"""

=======
WebbPSF
========

An object-oriented modeling system for the JWST instruments.

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
import atpy
import pyfits

import poppy
#from poppy.utils import * # functions for measuring the PSF's properties, etc.
#from fwcentroid import fwcentroid

#from _version import __version__


try: 
    import pysynphot
    _HAS_PYSYNPHOT = True
except:
    _HAS_PYSYNPHOT = False
 

import logging
_log = logging.getLogger('webbpsf')
_log.setLevel(logging.DEBUG)
_log.setLevel(logging.INFO)



try:
    __IPYTHON__
    from IPython.Debugger import Tracer; stop = Tracer()
except:
    def stop(): 
        pass




class JWInstrument(poppy.Instrument):
    """ A generic JWST Instrument class.

    *Note*: Do not use this class directly - instead use one of the :ref:`specific_instrument` subclasses!

    This class provides a simple interface for modeling PSF formation through the JWST instruments, 
    with configuration options and software interface loosely resembling the configuration of the instrument 
    mechanisms.   
    
    This module currently only provides a modicum of error checking, and relies on the user
    being knowledgable enough to avoid trying to simulate some physically impossible or just plain silly
    configuration (such as trying to use a FQPM with the wrong filter).

    The instrument constructors do not take any arguments. Instead, create an instrument object and then
    configure the `filter` or other attributes as desired. The most commonly accessed parameters are 
    available as object attributes: `filter`, `image_mask`, `pupil_mask`, `pupilopd`. More advanced
    configuration can be done by editing the :ref:`JWInstrument.options` dictionary, either by passing options to __init__ or by directly editing the dict afterwards.
    """

    def __init__(self, name=""):
        self.name=name

        self._WebbPSF_basepath = os.getenv('WEBBPSF_PATH', default= os.path.dirname(os.path.dirname(os.path.abspath(poppy.__file__))) +os.sep+"data" )

        self._datapath = self._WebbPSF_basepath + os.sep + self.name + os.sep
        self._filter = None
        self._image_mask = None
        self._pupil_mask = None
        self.pupil = os.path.abspath(self._datapath+"../pupil_RevV.fits")
        "Filename *or* pyfits.HDUList for JWST pupil mask. Usually there is no need to change this."
        self.pupilopd = None   # This can optionally be set to a tuple indicating (filename, slice in datacube)
        """Filename *or* pyfits.HDUList for JWST pupil OPD. 
        
        This can be either a full absolute filename, or a relative name in which case it is
        assumed to be within the instrument's `data/OPDs/` directory, or an actual pyfits.HDUList object corresponding to such a file.
        If the file contains a datacube, you may set this to a tuple (filename, slice) to select a given slice, or else
        the first slice will be used."""

        self.options = {} # dict for storing other arbitrary options. 
        """ A dictionary capable of storing other arbitrary options, for extensibility. The following are all optional, and
        may or may not be meaningful depending on which instrument is selected.

        Parameters
        ----------
        source_offset_r : float
            Radial offset of the target from the center, in arcseconds
        source_offset_theta : float
            Position angle for that offset
        pupil_shift_x, pupil_shift_y : float
            Relative shift of a coronagraphic pupil in X and Y, expressed as a decimal between 0.0-1.0
            Note that shifting an array too much will wrap around to the other side unphysically, but
            for reasonable values of shift this is a non-issue.
        rebin : bool
            For output files, write an additional FITS extension including a version of the output array 
            rebinned down to the actual detector pixel scale?
        jitter : string
            Type of jitter model to apply. Currently not implemented
        parity : string "even" or "odd"
            You may wish to ensure that the output PSF grid has either an odd or even number of pixels.
            Setting this option will force that to be the case by increasing npix by one if necessary.
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

        #create private instance variables. These will be
        # wrapped just below to create properties with validation.
        self._filter=None

        filter_table = atpy.Table(self._WebbPSF_basepath + os.sep+ 'filters.txt',type='ascii',delimiter='\t')
        wmatch = np.where(filter_table.instrument == self.name)
        self.filter_list = filter_table.filter[wmatch].tolist()
        "List of available filters"
        self._filter_nlambda_default = dict(zip(filter_table.filter[wmatch], filter_table.nlambda[wmatch]))

        #self._filter_files= [os.path.abspath(f) for f in glob.glob(self._datapath+os.sep+'filters/*_thru.fits')]
        #self.filter_list=[os.path.basename(f).split("_")[0] for f in self._filter_files]
        if len(self.filter_list) ==0: 
            self.filter_list=[''] # don't crash for FGS which lacks filters in the usual sense
            stop()

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

        self.pixelscale = 0.0
        "Detector pixel scale, in arcsec/pixel"
        self._spectra_cache = {}  # for caching pysynphot results.

    def _validate_config(self):
        pass

    # create properties with error checking
    @property
    def filter(self):
        'Currently selected filter name (e.g. "F200W")'
        return self._filter
    @filter.setter
    def filter(self, value):
        value = value.upper() # force to uppercase
        if value not in self.filter_list:
            raise ValueError("Instrument %s doesn't have a filter called %s." % (self.name, value))
        self._filter = value
        self._validate_config()

    @property
    def image_mask(self):
        'Currently selected coronagraphic image plane mask, or None for direct imaging'
        return self._image_mask
    @image_mask.setter
    def image_mask(self, name):
        if name is "": name = None
        if name is not None:
            name = name.upper() # force to uppercase
            if name not in self.image_mask_list:
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
            name = name.upper() # force to uppercase
            if name not in self.pupil_mask_list:
                raise ValueError("Instrument %s doesn't have an pupil mask called %s." % (self.name, name))

        self._pupil_mask = name

    def __str__(self):
        return "JWInstrument name="+self.name

    #----- actual optical calculations follow here -----
    def calcPSF(self, outfile=None, source=None, filter=None,  nlambda=None, monochromatic=None ,
            fov_arcsec=None, fov_pixels=None,  oversample=None, detector_oversample=None, fft_oversample=None, rebin=True,
            clobber=True, display=False, save_intermediates=False, return_intermediates=False):
        """ Compute a PSF.
        The result can either be written to disk (set outfile="filename") or else will be returned as
        a pyfits HDUlist object.


        Output sampling may be specified in one of two ways: 

        1) Set `oversample=<number>`. This will use that oversampling factor beyond detector pixels
           for output images, and beyond Nyquist sampling for any FFTs to prior optical planes. 
        2) set `detector_oversample=<number>` and `fft_oversample=<other_number>`. This syntax lets
           you specify distinct oversampling factors for intermediate and final planes. 

        By default, both oversampling factors are set equal to 2.

        Notes
        -----
        More advanced PSF computation options (pupil shifts, source positions, jitter, ...)
        may be set by configuring the `.options` dictionary attribute of this class.

        Parameters
        ----------
        filter : string, optional
            Filter name. Setting this is just a shortcut for setting the object's filter first, then
            calling calcPSF afterwards.
        source : pysynphot.SourceSpectrum or dict
            specification of source input spectrum. Default is a 5700 K sunlike star.
        nlambda : int
            How many wavelengths to model for broadband? 
            The default depends on how wide the filter is: (5,3,1) for types (W,M,N) respectively
        monochromatic : float, optional
            Setting this to a wavelength value (in meters) will compute a monochromatic PSF at that 
            wavelength, overriding filter and nlambda settings.
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
            This is useful if you want to e.g. examine the intensity in the Lyot plane for a coronagraphic propagation.

        Returns
        -------
        outfits : pyfits.HDUList
            The output PSF is returned as a pyfits.HDUlist object.
            If `outfile` is set to a valid filename, the output is also written to that file.


        """
        if filter is not None:
            self.filter = filter

        local_options = self.options  # all local state should be stored in a dict, for
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
                #if self.name=='TFI':    # filter names are irrelevant for TFI.
                    #nlambda=5
                #else:
                    #filt_width = self.filter[-1]
                    #lookup_table = {'NIRCam': {'2': 10, 'W':20,'M':3,'N':1}, 
                                    #'NIRSpec':{'W':5,'M':3,'N':1}, 
                                    #'MIRI':{'W':5,'M':3,'N':1}, 
                                    #'FGS':{'W':5,'M':3,'N':1}}

                    #nlambda = lookup_table[self.name][filt_width]
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
            fov_spec = 'pixels = %d' % fov_pixels
        elif fov_arcsec is not None:
            fov_spec = 'arcsec = %f' % fov_arcsec

        _log.debug('FOV set to '+fov_spec)

        #---- Implement the semi-convoluted logic for the oversampling options. See docstring above
        if oversample is not None and detector_oversample is not None and fft_oversample is not None:
            # all options set, contradictorily -> complain!
            raise ValueError("You cannot specify simultaneously the oversample= option with the detector_oversample and fft_oversample options. Pick one or the other!")
        elif oversample is None and detector_oversample is None and fft_oversample is None:
            # nothing set -> set oversample = 4
            oversample = 4
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
            fft_oversample=fft_oversample, detector_oversample=detector_oversample)
        # and use it to compute the PSF (the real work happens here, in code in poppy.py)
        #result = self.optsys.calcPSF(source, display_intermediates=display, save_intermediates=save_intermediates, display=display)
        #if _USE_MULTIPROC and monochromatic is None :
            #result = self.optsys.calcPSFmultiproc(source, nprocesses=_MULTIPROC_NPROCESS) # no fancy display args for multiproc.
        #else:
        result = self.optsys.calcPSF(wavelens, weights, display_intermediates=display, display=display, save_intermediates=save_intermediates, return_intermediates=return_intermediates)

        if return_intermediates: # this implies we got handed back a tuple, so split it apart
            result, intermediates = result


        self._getFITSHeader(result, local_options)

        self._calcPSF_format_output(result, local_options)


        if display:
            f = plt.gcf()
            #p.text( 0.1, 0.95, "%s, filter= %s" % (self.name, self.filter), transform=f.transFigure, size='xx-large')
            plt.suptitle( "%s, filter= %s" % (self.name, self.filter), size='xx-large')
            plt.text( 0.99, 0.04, "Calculation with %d wavelengths (%g - %g um)" % (nlambda, wavelens[0]*1e6, wavelens[-1]*1e6), transform=f.transFigure, horizontalalignment='right')

        if outfile is not None:
            result[0].header.update ("FILENAME", os.path.basename (outfile),
                           comment="Name of this file")
            result.writeto(outfile, clobber=clobber)
            _log.info("Saved result to "+outfile)

        if return_intermediates:
            return result, intermediates
        else:
            return result

    def _getFITSHeader(self, result, options):
        """ populate FITS Header keywords """
        poppy.Instrument._getFITSHeader(self,result, options)
        result[0].header.update('FILTER', self.filter, 'Filter name')
        if self.image_mask is not None:
            result[0].header.update('CORONMSK', self.image_mask)
        if self.pupil_mask is not None:
            result[0].header.update('PUPIL', self.pupil_mask)


    def _calcPSF_format_output(self, result, options):
        """ Apply desired formatting to output file:
                 - rebin to detector pixel scale if desired
                 - set up FITS extensions if desired
                 - output either the oversampled, rebinned, or both

            Modifies the 'result' HDUList object.
        """
        output_mode = options.get('output_mode','Both as FITS extensions')

        if output_mode == 'Mock JWST DMS Output':
            # first rebin down to detector sampling
            # then call mockdms routines to embed in larger detector etc
            raise NotImplementedError('Not implemented yet')
        else:
            poppy.Instrument._calcPSF_format_output(self, result, options)


    def _getOpticalSystem(self,fft_oversample=2, detector_oversample = None, fov_arcsec=2, fov_pixels=None):
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


        if detector_oversample is None: detector_oversample = fft_oversample

        _log.debug("Oversample: %d  %d " % (fft_oversample, detector_oversample))
        optsys = poppy.OpticalSystem(name='JWST+'+self.name, oversample=fft_oversample)
        if 'source_offset_r' in self.options.keys(): optsys.source_offset_r = self.options['source_offset_r']
        if 'source_offset_theta' in self.options.keys(): optsys.source_offset_theta = self.options['source_offset_theta']


        #---- set pupil OPD
        if isinstance(self.pupilopd, str):  # simple filename
            full_opd_path = self.pupilopd if os.path.exists( self.pupilopd) else os.path.join(self._datapath, "OPD",self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], basestring): # tuple with filename and slice
            full_opd_path =  (self.pupilopd[0] if os.path.exists( self.pupilopd[0]) else os.path.join(self._datapath, "OPD",self.pupilopd[0]), self.pupilopd[1])
        elif isinstance(self.pupilopd, pyfits.HDUList): # OPD supplied as FITS object
            full_opd_path = self.pupilopd # not a path per se but this works correctly to pass it to poppy
        elif self.pupilopd is None: 
            full_opd_path = None
        else:
            raise TypeError("Not sure what to do with a pupilopd of that type:"+str(type(self.pupilopd)))

        #---- set pupil intensity
        if isinstance(self.pupil, str): # simple filename
            full_pupil_path = self.pupil if os.path.exists( self.pupil) else os.path.join(self._WebbPSF_basepath,self.pupil)
        elif isinstance(self.pupil, pyfits.HDUList): # pupil supplied as FITS object
            full_pupil_path = self.pupil
        else: 
            raise TypeError("Not sure what to do with a pupil of that type:"+str(type(self.pupil)))


        #---- apply pupil intensity and OPD to the optical model
        optsys.addPupil(name='JWST Pupil', transmission=full_pupil_path, opd=full_opd_path, opdunits='micron', rotation=self._rotation)

        #---- add coronagraphy if requested, and flag to invoke semi-analytic coronagraphic propagation
        if self.image_mask is not None or self.pupil_mask is not None or ('force_coron' in self.options.keys() and self.options['force_coron']):
            optsys, trySAM, SAM_box_size = self._addCoronagraphOptics(optsys, oversample=fft_oversample)
        else: trySAM = False

        #--- add the detector element. 
        if fov_pixels is None:
            fov_pixels = np.round(fov_arcsec/self.pixelscale)
            if 'parity' in self.options.keys():
                if self.options['parity'].lower() == 'odd'  and np.remainder(fov_pixels,2)==0: fov_pixels +=1
                if self.options['parity'].lower() == 'even' and np.remainder(fov_pixels,2)==1: fov_pixels +=1

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

    def _addCoronagraphOptics(self,optsys, oversample=2):
        """Add coronagraphic optics to an optical system. 
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
            # the requested band is not yet supported in synphot/CDBS. (those files are still a
            # work in progress...). Therefore, use our local throughput files and create a synphot
            # transmission object.
            wf = np.where(filtername == np.asarray(self.filter_list))[0]
            # The existing FITS files all have wavelength in ANGSTROMS since that is the pysynphot convention...
            filterdata = atpy.Table(self._filter_files[wf], type='fits')

            _log.warn("Filter %s not supported in pysynphot/CDBS yet. Falling back to local filter transmission files" % filtername)
            _log.warn("These may be less accurate.")
            band = pysynphot.spectrum.ArraySpectralElement(throughput=filterdata.THROUGHPUT,
                                wave=filterdata.WAVELENGTH, waveunits='angstrom',name=filtername)
        return band



###########

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

        self.image_mask_list = ['FQPM1065', 'FQPM1140', 'FQPM1550', 'LYOT2300']
        self.pupil_mask_list = ['MASKFQPM', 'MASKLYOT']

        for i in range(4):
            self.filter_list.append('MRS-IFU Ch%d'% (i+1) )
        self.monochromatic= 8.0
        self._IFU_pixelscale = {'Ch1':(0.18, 0.19), 'Ch2':(0.28, 0.19), 'Ch3': (0.39, 0.24), 'Ch4': (0.64, 0.27) }
            # The above tuples give the pixel resolution (perpendicular to the slice, along the slice). 
            # The pixels are not square.

        self._default_aperture='MIRIM_center' # reference into SIAF for ITM simulation V/O coords


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


    def _addCoronagraphOptics(self,optsys, oversample=2):
        """Add coronagraphic optics for MIRI.
        Semi-analytic coronagraphy algorithm used for the Lyot only.

        """

        # For MIRI coronagraphy, all the coronagraphic optics are rotated the same
        # angle as the instrument is, relative to the primary. So they see the unrotated
        # telescope pupil.
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
        if (self.image_mask is not None and 'FQPM' in self.image_mask) or 'force_fqpm_shift' in self.options.keys() : optsys.addPupil("FQPM FFT aligner")

        if self.image_mask == 'FQPM1065':
            #optsys.addImage() # null debugging image plane FIXME
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1065",
                opticslist = [  poppy.IdealFQPM(wavelength=10.65e-6, name=self.image_mask),
                                poppy.IdealFieldStop(size=24, angle=-4.56)])
            optsys.addImage(container)
            trySAM = False
            SAM_box_size = 1.0 # irrelevant but variable still needs to be set.
        elif self.image_mask == 'FQPM1140':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1140",
                opticslist = [  poppy.IdealFQPM(wavelength=11.40e-6, name=self.image_mask),
                                poppy.IdealFieldStop(size=24, angle=-4.56)])
            optsys.addImage(container)
            trySAM = False
            SAM_box_size = 1.0 # irrelevant but variable still needs to be set.
        elif self.image_mask == 'FQPM1550':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1550",
                opticslist = [  poppy.IdealFQPM(wavelength=15.50e-6, name=self.image_mask),
                                poppy.IdealFieldStop(size=24, angle=-4.56)])
            optsys.addImage(container)
            trySAM = False
            SAM_box_size = 1.0 # irrelevant but variable still needs to be set.
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
                              poppy.IdealFieldStop(size=30, angle=-4.56)] )
            optsys.addImage(container)
            trySAM = True
            SAM_box_size = [5,20]
        else:
            optsys.addImage()
            trySAM = False
            SAM_box_size= 1.0 # irrelevant but variable still needs to be set.

        if (self.image_mask is not None and 'FQPM' in self.image_mask)  or 'force_fqpm_shift' in self.options.keys() : optsys.addPupil("FQPM FFT aligner", direction='backward')

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
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MIRI_FQPMLyotStop.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKLYOT':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MIRI_LyotLyotStop.fits.gz", name=self.pupil_mask, shift=shift)
        else: # all the MIRI filters have a tricontagon outline, even the non-coron ones.
            optsys.addPupil(transmission=self._WebbPSF_basepath+"/tricontagon.fits", name = 'filter cold stop', shift=shift)

        optsys.addRotation(self._rotation)

        return (optsys, trySAM, SAM_box_size)

    def _getFITSHeader(self, hdulist, options):
        """ Format MIRI-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)

        hdulist[0].header.update('GRATNG14','None', 'MRS Grating for channels 1 and 4')
        hdulist[0].header.update('GRATNG23','None', 'MRS Grating for channels 2 and 3')
        hdulist[0].header.update('FLATTYPE','?', 'Type of flat field to be used: all, one, principal')
        hdulist[0].header.update('CCCSTATE','open', 'Contamination Control Cover state: open, closed, locked')
        if self.image_mask is not None:
            hdulist[0].header.update('TACQNAME','None', 'Target acquisition file name')


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



    def _validate_config(self):
        """For NIRCam, this checks whenever you change a filter and updates the pixelscale appropriately"""
        filtwave = float(self.filter[1:4])/100
        newscale = self._pixelscale_short if filtwave < 2.4 else self._pixelscale_long
        # update the pixel scale if it has changed *and*
        # only if the user has not already set the pixel scale to some custom value
        if newscale != self.pixelscale and (self.pixelscale == self._pixelscale_short or self.pixelscale==self._pixelscale_long):
            self.pixelscale = newscale
            _log.info("NIRCam pixel scale updated to %f arcsec/pixel to match channel for the selected filter." % self.pixelscale)


    def _addCoronagraphOptics(self,optsys, oversample=2):
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
            optsys.addPupil(transmission=self._datapath+"/coronagraph/NIRCam_Lyot_Somb.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'WEDGELYOT':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/NIRCam_Lyot_Sinc.fits", name=self.pupil_mask, shift=shift)
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
        elif self.pupil_mask == 'WEAK LENS -4 (=4+8)':
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

        hdulist[0].header.update('MODULE',self.module, 'NIRCam module: A or B')
        hdulist[0].header.update('CHANNEL', 'Short' if self.pixelscale == self._pixelscale_short else 'Long', 'NIRCam channel: long or short')
        # filter, pupil added by calcPSF header code
        hdulist[0].header.update('PILIN', 'False', 'Pupil imaging lens in optical path: T/F')



class NIRSpec(JWInstrument):
    """ A class modeling the optics of NIRSpec, in **imaging** mode. 

    This is not a substitute for a spectrograph model, but rather a way of simulating a PSF as it
    would appear with NIRSpec in imaging mode (e.g. for target acquisition).  NIRSpec support is 
    relatively simplistic compared to the other instruments at this point.
    
    Relevant attributes include `filter`. In addition to the actual filters, you may select 'IFU' to
    indicate use of the NIRSpec IFU, in which case use the `monochromatic` attribute to set the simulated wavelength.
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


        self._default_aperture='NIRSpec A center' # reference into SIAF for ITM simulation V/O coords

    def _validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('NIRSpec does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
    def _addCoronagraphOptics(self,optsys, oversample=2):
        raise NotImplementedError("No Coronagraph in NIRSpec!")


    def _validate_config(self):
        if self.filter.startswith("IFU"): raise NotImplementedError("The NIRSpec IFU is not yet implemented.")


    def _getFITSHeader(self, hdulist, options):
        """ Format NIRSpec-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)
        hdulist[0].header.update('GRATING', 'None', 'NIRSpec grating element name')
        hdulist[0].header.update('APERTURE', 'None', 'NIRSpec slit aperture name')


class NIRISS(JWInstrument):
    """ A class modeling the optics of the Near-IR Imager and Slit Spectrograph
        (formerly nTFI)
    
    Relevant attributes include `image_mask`, and `pupil_mask`.

    """
    def __init__(self):
        JWInstrument.__init__(self, "NIRISS")
        self.pixelscale = 0.064 

        self.image_mask_list = ['CORON058', 'CORON075','CORON150','CORON200'] # available but unlikely to be used...
        self.pupil_mask_list = ['MASK_NRM','CLEAR']
        self._default_aperture='NIRISS center' # reference into SIAF for ITM simulation V/O coords


    def _validate_config(self):
        pass

    def _addCoronagraphOptics(self,optsys, oversample=2):
        """Add coronagraphic optics for TFI/NIRISS. 

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
        elif (self.pupil_mask  is None and self.image_mask is not None):
            optsys.addPupil(name='No Lyot Mask Selected!')

        return (optsys, trySAM, radius+0.05) # always attempt to cast this to a SemiAnalyticCoronagraph

    def _getFITSHeader(self, hdulist, options):
        """ Format NIRISS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)

        if self.image_mask is not None:
            hdulist[0].header.update('CORONPOS', self.image_mask, 'NIRISS coronagraph spot location')
        hdulist[0].header.update('FOCUSPOS',0,'NIRISS focus mechanism not yet modeled.')


class TFI(JWInstrument):
    """ A class modeling the optics of the Tunable Filter Imager

    ** This class is preserved here for archival/historical purposes in this version of WebbPSF. It is now
    deprecated in favor of NIRISS. **
    
    Relevant attributes include `image_mask`, and `pupil_mask`.

    Because of its tunable etalon, wavelength selection for TFI is handled a bit differently than
    for the other SIs.  The `filter` attribute, while present, is not used. Instead, there is an
    `etalon_wavelength` attribute, which is the wavelength in microns that the etalon is tuned to.
    Acceptable values are between 1.5 - 2.7 and 3.0 - 5.0 microns. The effective resolution for the TFI
    at any given resolution is obtained from a lookup table and used to calculate the PSF across the
    resulting bandpass.

    You may also use the `monochromatic=` option to `calcPSF()` to calculate a PSF at a single wavelength.
    """
    def __init__(self):
        JWInstrument.__init__(self, "TFI")
        self.pixelscale = 0.064 # for TFI

        self.filter_list = [""]
        self.filter=""

        self.image_mask_list = ['CORON058', 'CORON075','CORON150','CORON200']
        self.pupil_mask_list = ['MASKC21N','MASKC66N','MASKC71N','MASK_NRM','CLEAR']
        self.etalon_wavelength = 2.0
        """ Tunable filter etalon wavelength setting """
        self.resolution_table = atpy.Table(self._datapath+os.sep+"filters/TFI_resolution.txt", type='ascii',names=('wavelength','resolution'))



    def _validate_config(self):
        pass

    def _addCoronagraphOptics(self,optsys, oversample=2):
        """Add coronagraphic optics for TFI
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

        if self.pupil_mask == 'MASKC21N':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC21N.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKC66N':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC66N.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKC71N':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC71N.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASK_NRM':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASK_NRM.fits.gz", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'CLEAR':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKCLEAR.fits.gz", name=self.pupil_mask, shift=shift)
        elif (self.pupil_mask  is None and self.image_mask is not None):
            optsys.addPupil(name='No Lyot Mask Selected!')

        return (optsys, trySAM, radius+0.05) # always attempt to cast this to a SemiAnalyticCoronagraph

    def _getSynphotBandpass(self, filtername):
        """ Return a pysynphot.ObsBandpass object for the given desired band.
        This uses a lookup table to predict the properties of the TFI tunable filter etalon.
        """

        # filter name parameter is ignored, instead this uses etalon_wavelength
        
        if (self.etalon_wavelength < 1.5 or self.etalon_wavelength > 5.0 or
            (self.etalon_wavelength > 2.7 and self.etalon_wavelength < 3.0)):
            raise ValueError("Invalid value for etalon wavelength: %f. Please set a value in 1.5-2.7 or 3.0-5.0 microns." % self.etalon_wavelength)


        match_index = np.abs(self.resolution_table.wavelength - self.etalon_wavelength).argmin()
        resolution = self.resolution_table.resolution[match_index]
        _log.info("Etalon wavelength %.3f has resolution %.2f" % (self.etalon_wavelength, resolution))
        wavelen = np.linspace(1.0, 5.0, 1000)

        fwhm = self.etalon_wavelength/resolution
        sigma = fwhm / 2.35482

        transmission = np.exp(- (wavelen - self.etalon_wavelength)**2/ (2*sigma**2))

        #plt.plot(wavelen, transmission)
        band = pysynphot.ArrayBandpass(wave=wavelen*1e4, throughput=transmission, waveunits='angstrom',name='TFI-etalon-%.3d' % self.etalon_wavelength)

        self.filter = '%.3f um' % self.etalon_wavelength # update this. used for display and FITS header info

        return band


    def _getSpecCacheKey(self, source, nlambda):
        """ return key for the cache of precomputed spectral weightings.
        This is a separate function so the TFI subclass can override it.
        """
        return ("%.3f" %self.etalon_wavelength, source.name, nlambda)


    def filter(self, value): # we just store a string here for the wavelength... don't worry about validation.
        self._filter = value

class FGS(JWInstrument):
    """ A class modeling the optics of the FGS.
    
    Not a lot to see here, folks: There are no selectable options, just a great big detector-wide bandpass.
    """
    def __init__(self):
        JWInstrument.__init__(self, "FGS")
        self.pixelscale = 0.069 # for FGS
        self._default_aperture='FGS1 center' # reference into SIAF for ITM simulation V/O coords

    def _validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('FGS does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
        #TODO only one possible filter fot the FGS, too. 
    def _addCoronagraphOptics(self,optsys):
        raise NotImplementedError("No Coronagraph in FGS!")

    def _getFITSHeader(self, hdulist, options):
        """ Format FGS-like FITS headers, based on JWST DMS SRD 1 FITS keyword info """
        JWInstrument._getFITSHeader(self, hdulist, options)
        hdulist[0].header.update('FOCUSPOS',0,'FGS focus mechanism not yet modeled.')



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
        return pyfits.open(filename)
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











def makeFakeFilter(filename, lcenter, dlam, instrument='NIRCam', name='filter', source='Fake top hat filter', clobber=False):
    """ arguments in microns, but file written in angstroms """

    lstart = lcenter - dlam/2
    lstop = lcenter + dlam/2

    nlambda = 40

    _log.info("Filter from %f - %f " % (lstart, lstop))
    wavelength = np.linspace( lstart-dlam*0.1, lstop+dlam*0.1, nlambda)
    _log.debug(wavelength)
    transmission = np.zeros_like(wavelength)
    transmission[np.where( (wavelength > lstart) & (wavelength < lstop) )] = 1.0

    t = atpy.Table()
    t.add_column('WAVELENGTH', wavelength*1e4, unit='angstrom')
    t.add_column('THROUGHPUT', transmission)

    t.add_keyword('TELESCOP','JWST')
    t.add_keyword('INSTRUME',instrument)
    t.add_keyword('FILTER',name)

    t.add_keyword('SOURCE', source)
    t.add_comment("This is a fake filter profile, represented as a top-hat function.")

    t.add_keyword("LAMBDA0",lcenter)
    t.add_keyword("DELTALAM",dlam)

    t.write(filename, overwrite=clobber)
    _log.info("Created fake filter profile in "+filename)


    return t
def _makeMIRIfilters():
    makeFakeFilter('F560W_throughput.fits', 5.6, 1.2, clobber=True)
    makeFakeFilter('F770W_throughput.fits', 7.7, 2.2, clobber=True)
    makeFakeFilter('F1000W_throughput.fits', 10, 2, clobber=True)
    makeFakeFilter('F1130W_throughput.fits', 11.3, 0.7, clobber=True)
    makeFakeFilter('F1280W_throughput.fits', 12.8, 2.4, clobber=True)
    makeFakeFilter('F1500W_throughput.fits', 15, 3, clobber=True)
    makeFakeFilter('F1800W_throughput.fits', 18, 3, clobber=True)
    makeFakeFilter('F2100W_throughput.fits', 21, 5, clobber=True)
    makeFakeFilter('F2550W_throughput.fits', 25.5, 4, clobber=True)
    makeFakeFilter('F1065C_throughput.fits',10.65, 0.53,clobber=True)
    makeFakeFilter('F1140C_throughput.fits',11.40, 0.57,clobber=True)
    makeFakeFilter('F1550C_throughput.fits',15.50, 0.78,clobber=True)
    makeFakeFilter('F2300C_throughput.fits',23.00, 4.60,clobber=True)
    makeFakeFilter('FND_throughput.fits', 11.5, 7.0,clobber=True)



def _makeNIRCamFilters():
    "Create nircam filters based on http://ircamera.as.arizona.edu/nircam/features.html "
    makeFakeFilter('F070W_thru.fits', 0.7000, 0.1750, clobber=True)
    makeFakeFilter('F090W_thru.fits', 0.9000, 0.2250, clobber=True)
    makeFakeFilter('F115W_thru.fits', 1.1500, 0.2875, clobber=True)
    makeFakeFilter('F150W_thru.fits', 1.5000, 0.3750, clobber=True)
    makeFakeFilter('F150W2_thru.fits', 1.5000, 1.0000, clobber=True)
    makeFakeFilter('F200W_thru.fits', 2.0000, 0.5000, clobber=True)
    makeFakeFilter('F277W_thru.fits', 2.7700, 0.6925, clobber=True)
    makeFakeFilter('F322W2_thru.fits', 3.2200, 1.6100, clobber=True)
    makeFakeFilter('F356W_thru.fits', 3.5600, 0.8900, clobber=True)
    makeFakeFilter('F444W_thru.fits', 4.4400, 1.1100, clobber=True)
    makeFakeFilter('F140M_thru.fits', 1.4000, 0.1400, clobber=True)
    makeFakeFilter('F162M_thru.fits', 1.6200, 0.1510, clobber=True)
    makeFakeFilter('F182M_thru.fits', 1.8200, 0.2210, clobber=True)
    makeFakeFilter('F210M_thru.fits', 2.1000, 0.2100, clobber=True)
    makeFakeFilter('F250M_thru.fits', 2.5000, 0.1667, clobber=True)
    makeFakeFilter('F300M_thru.fits', 3.0000, 0.3000, clobber=True)
    makeFakeFilter('F335M_thru.fits', 3.3500, 0.3350, clobber=True)
    makeFakeFilter('F360M_thru.fits', 3.6000, 0.3600, clobber=True)
    makeFakeFilter('F410M_thru.fits', 4.1000, 0.4100, clobber=True)
    makeFakeFilter('F430M_thru.fits', 4.3000, 0.2000, clobber=True)
    makeFakeFilter('F460M_thru.fits', 4.6000, 0.2000, clobber=True)
    makeFakeFilter('F480M_thru.fits', 4.8000, 0.4000, clobber=True)
    makeFakeFilter('F164N_thru.fits', 1.6440, 0.0164, clobber=True)
    makeFakeFilter('F187N_thru.fits', 1.8756, 0.0188, clobber=True)
    makeFakeFilter('F212N_thru.fits', 2.1218, 0.0212, clobber=True)
    makeFakeFilter('F225N_thru.fits', 2.2477, 0.0225, clobber=True)
    makeFakeFilter('F323N_thru.fits', 3.2350, 0.0324, clobber=True)
    makeFakeFilter('F405N_thru.fits', 4.0523, 0.0405, clobber=True)
    makeFakeFilter('F418N_thru.fits', 4.1813, 0.0418, clobber=True)
    makeFakeFilter('F466N_thru.fits', 4.6560, 0.0466, clobber=True)
    makeFakeFilter('F470N_thru.fits', 4.7050, 0.0471, clobber=True)

def _makeNIRSpecFilters():
    makeFakeFilter('F110W_throughput.fits', 1.4,1.2, clobber=True, source='NIRSpec Acq Filter docs',name='F110W',instrument='NIRSpec')
    makeFakeFilter('F140X_throughput.fits', 1.1,0.2, clobber=True, source='NIRSpec Acq Filter docs',name='F110W',instrument='NIRSpec')

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
