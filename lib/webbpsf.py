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
    * TFI
    * FGS


WebbPSF makes use of python's ``logging`` facility for log messages, using
the logger name "webbpsf".



Code by Marshall Perrin <mperrin@stsci.edu>

"""
import os
import types
import glob
import time
import numpy as N
import scipy.interpolate, scipy.ndimage
import pylab as P
import matplotlib
import atpy
import pyfits
from matplotlib.colors import LogNorm  # for log scaling of images, with automatic colorbar support

import poppy
from fwcentroid import fwcentroid

__version__ = poppy.__version__


try: 
    import pysynphot
    _HAS_PYSYNPHOT = True
except:
    _HAS_PYSYNPHOT = False
 

import logging
_log = logging.getLogger('webbpsf')
_log.setLevel(logging.DEBUG)
_log.setLevel(logging.INFO)
#_log.addHandler(logging.NullHandler())



try:
    __IPYTHON__
    from IPython.Debugger import Tracer; stop = Tracer()
except:
    def stop(): 
        pass




class JWInstrument(object):
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

    def __init__(self, name=None):
        self.name=name

        self._WebbPSF_basepath = os.getenv('WEBBPSF_PATH', default= os.path.dirname(os.path.dirname(os.path.abspath(poppy.__file__))) +os.sep+"data" )

        self._datapath = self._WebbPSF_basepath + os.sep + self.name + os.sep
        self._filter = None
        self._image_mask = None
        self._pupil_mask = None
        self.pupil = os.path.abspath(self._datapath+"../pupil_RevV.fits")
        "Filename for JWST pupil mask. Usually there is no need to change this."
        self.pupilopd = None   # This can optionally be set to a tuple indicating (filename, slice in datacube)
        """Filename for JWST pupil OPD. This can be either a full absolute filename, or a relative name in which case it is
        assumed to be within the instrument's `data/OPDs/` directory.
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
            Type of jitter model to apply.
        parity : string "even" or "odd"
            You may wish to ensure that the output PSF grid has either an odd or even number of pixels.
            Setting this option will force that to be the case by increasing npix by one if necessary.
        force_coron : bool
            Set this to force full coronagraphic optical propagation when it might not otherwise take place
            (e.g. calculate the non-coronagraphic images in the same way as coronagraphy is done, rather than
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
        wmatch = N.where(filter_table.instrument == self.name)
        self.filter_list = filter_table.filter[wmatch].tolist()
        "List of available filters"
        self._filter_nlambda_default = dict(zip(filter_table.filter[wmatch], filter_table.nlambda[wmatch]))

        #self._filter_files= [os.path.abspath(f) for f in glob.glob(self._datapath+os.sep+'filters/*_thru.fits')]
        #self.filter_list=[os.path.basename(f).split("_")[0] for f in self._filter_files]
        if len(self.filter_list) ==0: self.filter_list=[''] # don't crash for TFI or FGS which lack filters in the usual sense

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
            fov_arcsec=None, fov_pixels=None,  oversample=None, detector_oversample=None, calc_oversample=None, rebin=True,
            clobber=True, display=False, save_intermediates=False):
        """ Compute a PSF.
        The result can either be written to disk (set outfile="filename") or else will be returned as
        a pyfits HDUlist object.


        Output sampling may be specified in one of two ways: 

        1) Set `oversample=<number>`. This will use that oversampling factor beyond detector pixels
           for output images, and beyond Nyquist sampling for any FFTs to prior optical planes. 
        2) set `detector_oversample=<number>` and `calc_oversample=<other_number>`. This syntax lets
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
        oversample, detector_oversample, calc_oversample : int
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

        Returns
        -------
        outfits : pyfits.HDUList
            The output PSF is returned as a pyfits.HDUlist object.
            If `outfile` is set to a valid filename, the output is also written to that file.


        """
        if filter is not None:
            self.filter = filter

    
        #----- choose # of wavelengths intelligently. Do this first before generating the source spectrum weighting.
        if nlambda is None or nlambda==0:
            # Automatically determine number of appropriate wavelengths.
            # Make selection based on filter configuration file
            try:
                if self.name=='TFI':    # filter names are irrelevant for TFI.
                    nlambda=5
                else:
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



        #----- calculate field of view depending on supplied parameters
        if fov_arcsec is None and fov_pixels is None:  #pick decent defaults.
            if self.name =='MIRI': fov_arcsec=12.
            else: fov_arcsec=5.
            fov_spec = 'arcsec = %f' % fov_arcsec
        elif fov_pixels is not None:
            fov_spec = 'pixels = %d' % fov_pixels
        elif fov_arcsec is not None:
            fov_spec = 'arcsec = %f' % fov_arcsec


        #---- Implement the semi-convoluted logic for the oversampling options. See docstring above
        if oversample is not None and detector_oversample is not None and calc_oversample is not None:
            # all options set, contradictorily -> complain!
            raise ValueError("You cannot specify simultaneously the oversample= option with the detector_oversample and calc_oversample options. Pick one or the other!")
        elif oversample is None and detector_oversample is None and calc_oversample is None:
            # nothing set -> set oversample = 4
            oversample = 4
        if detector_oversample is None: detector_oversample = oversample
        if calc_oversample is None: calc_oversample = oversample

        _log.info("PSF calc using fov_%s, oversample = %d, nlambda = %d" % (fov_spec, detector_oversample, nlambda) )

        #----- compute weights for each wavelength based on source spectrum
        wavelens, weights = self._getWeights(source=source, nlambda=nlambda, monochromatic=monochromatic)


        #---- now at last, actually do the PSF calc:
        #  instantiate an optical system using the current parameters
        self.optsys = self._getOpticalSystem(fov_arcsec=fov_arcsec, fov_pixels=fov_pixels,
            calc_oversample=calc_oversample, detector_oversample=detector_oversample)
        # and use it to compute the PSF (the real work happens here, in code in poppy.py)
        #result = self.optsys.calcPSF(source, display_intermediates=display, save_intermediates=save_intermediates, display=display)
        #if _USE_MULTIPROC and monochromatic is None :
            #result = self.optsys.calcPSFmultiproc(source, nprocesses=_MULTIPROC_NPROCESS) # no fancy display args for multiproc.
        #else:
        result = self.optsys.calcPSF(wavelens, weights, display_intermediates=display, save_intermediates=save_intermediates, display=display)


        #---  update FITS header, display, and output.
        result[0].header.update('PUPIL', os.path.basename(self.pupil))
        if self.pupilopd is None:
            result[0].header.update('PUPILOPD', "NONE - perfect telescope! ")
        else:
            if isinstance(self.pupilopd, basestring):
                result[0].header.update('PUPILOPD', os.path.basename(self.pupilopd))
            else:
                result[0].header.update('PUPILOPD', "%s slice %d" % (os.path.basename(self.pupilopd[0]), self.pupilopd[1]))
        result[0].header.update('INSTRUME', self.name)
        result[0].header.update('FILTER', self.filter)
        if self.image_mask is not None:
            result[0].header.update('CORON', self.image_mask)
        if self.pupil_mask is not None:
            result[0].header.update('LYOTMASK', self.pupil_mask)
        result[0].header.update('EXTNAME', 'OVERSAMP')
        result[0].header.add_history('Created by JWPSF v4 ')
        result[0].header.update('OVERSAMP', calc_oversample, 'Oversampling factor for FFTs in computation')
        result[0].header.update('DET_SAMP', detector_oversample, 'Oversampling factor for MFT to detector plane')
        (year, month, day, hour, minute, second, weekday, DOY, DST) =  time.gmtime()
        result[0].header.update ("DATE", "%4d-%02d-%02dT%02d:%02d:%02d" % (year, month, day, hour, minute, second), "Date of calculation")
        result[0].header.update ("AUTHOR", "%s@%s" % (os.getenv('USER'), os.getenv('HOST')), "username@host for calculation")


        # Should we downsample? 
        if rebin and detector_oversample > 1:
            _log.info(" Downsampling to detector pixel scale.")
            rebinned_result = result[0].copy()
            rebinned_result.data = rebin_array(rebinned_result.data, rc=(detector_oversample, detector_oversample))
            rebinned_result.header.update('OVERSAMP', 1, 'These data are rebinned to detector pixels')
            rebinned_result.header.update('CALCSAMP', detector_oversample, 'This much oversampling used in calculation')
            rebinned_result.header.update('EXTNAME', 'DET_SAMP')
            rebinned_result.header['PIXELSCL'] *= detector_oversample
            result.append(rebinned_result)


        if display:
            f = P.gcf()
            #p.text( 0.1, 0.95, "%s, filter= %s" % (self.name, self.filter), transform=f.transFigure, size='xx-large')
            P.suptitle( "%s, filter= %s" % (self.name, self.filter), size='xx-large')
            P.text( 0.99, 0.04, "Calculation with %d wavelengths (%g - %g um)" % (nlambda, wavelens[0]*1e6, wavelens[-1]*1e6), transform=f.transFigure, horizontalalignment='right')

        if outfile is not None:
            result[0].header.update ("FILENAME", os.path.basename (outfile),
                           comment="Name of this file")
            result.writeto(outfile, clobber=clobber)
            _log.info("Saved result to "+outfile)
        return result

    def _getSpecCacheKey(self, source, nlambda):
        """ return key for the cache of precomputed spectral weightings.
        This is a separate function so the TFI subclass can override it.
        """
        return (self.filter, source.name, nlambda)


    def _getWeights(self, source=None, nlambda=5, monochromatic=None, verbose=False):
        if monochromatic is not None:
            _log.info(" monochromatic calculation requested.")
            return (N.asarray([monochromatic]),  N.asarray([1]) )

        elif _HAS_PYSYNPHOT and (isinstance(source, pysynphot.spectrum.SourceSpectrum)  or source is None):
            """ Given a pysynphot.SourceSpectrum object, perform synthetic photometry for
            nlambda bins spanning the wavelength range of interest.

            Because this calculation is kind of slow, cache results for reuse in the frequent
            case where one is computing many PSFs for the same spectral source.
            """
            if source is None:
                try:
                    source = pysynphot.Icat('ck04models',5700,0.0,2.0)
                except:
                    _log.error("Could not load Castelli & Kurucz stellar model from disk; falling back to 5700 K blackbody")
                    source = pysynphot.Blackbody(5700)

            try:
                key = self._getSpecCacheKey(source, nlambda)
                if key in self._spectra_cache.keys(): 
                    return self._spectra_cache[keys]
            except:
                pass  # in case sourcespectrum lacks a name element - just do the below calc.

            _log.info("Computing wavelength weights using synthetic photometry for %s..." % self.filter)
            band = self._getSynphotBandpass()
            # choose reasonable min and max wavelengths
            w_above10 = N.where(band.throughput > 0.10*band.throughput.max())

            minwave = band.wave[w_above10].min()
            maxwave = band.wave[w_above10].max()
            _log.debug("Min, max wavelengths = %f, %f" % (minwave/1e4, maxwave/1e4))
            # special case: ignore red leak for MIRI F560W, which has a negligible effect in practice
            if self.filter == 'F560W':
                _log.debug("Special case: setting max wavelength to 6.38 um to ignore red leak")
                maxwave = 63800.0
            elif self.filter == 'F1280W':
                _log.debug("Special case: setting max wavelength to 14.32 um to ignore red leak")
                maxwave = 143200.0
 
            wave_bin_edges =  N.linspace(minwave,maxwave,nlambda+1)
            wavesteps = (wave_bin_edges[:-1] +  wave_bin_edges[1:])/2
            deltawave = wave_bin_edges[1]-wave_bin_edges[0]
            effstims = []

            #t0= time.time()
            for wave in wavesteps:
                if verbose: _log.info("using band centered at %f with width %f" % (wave,deltawave))
                box = pysynphot.Box(wave, deltawave) * band
                if box.throughput.max() == 0:  # watch out for pathological cases with no overlap (happens with MIRI FND at high nlambda)
                    result = 0.0
                else:
                    result = pysynphot.Observation(source, box).effstim('counts')
                effstims.append(result)
            #t1 = time.time()
            #print "  that took %f seconds for %d wavelengths" % (t1-t0, nlambda)

            effstims = N.array(effstims)
            effstims /= effstims.sum()
            wave_m =  band.waveunits.Convert(wavesteps,'m') # convert to meters

            #newsource = {'wavelengths': wave_m, 'weights':effstims}
            newsource = (wave_m, effstims)
            if verbose: print newsource
            self._spectra_cache[ self._getSpecCacheKey(source,nlambda)] = newsource
            return newsource

        else:  #Fallback simple code for if we don't have pysynphot.
            _log.warning("Pysynphot unavailable (or invalid source supplied)!   Assuming flat # of counts versus wavelength.")
            # compute a source spectrum weighted by the desired filter curves.
            # TBD this will eventually use pysynphot, so don't write anything fancy for now!
            wf = N.where(self.filter == N.asarray(self.filter_list))[0]
            # The existing FITS files all have wavelength in ANGSTROMS since that is the pysynphot convention...
            filterdata = atpy.Table(self._filter_files[wf], type='fits')
            _log.warn("CAUTION: Just interpolating rather than integrating filter profile, over %d steps" % nlambda)
            wtrans = N.where(filterdata.THROUGHPUT > 0.4)
            if self.filter == 'FND':  # special case MIRI's ND filter since it is < 0.1% everywhere...
                wtrans = N.where(  ( filterdata.THROUGHPUT > 0.0005)  & (filterdata.WAVELENGTH > 7e-6*1e10) & (filterdata.WAVELENGTH < 26e-6*1e10 ))
            lrange = filterdata.WAVELENGTH[wtrans] *1e-10  # convert from Angstroms to Meters
            lambd = N.linspace(N.min(lrange), N.max(lrange), nlambda)
            filter_fn = scipy.interpolate.interp1d(filterdata.WAVELENGTH*1e-10, filterdata.THROUGHPUT,kind='cubic', bounds_error=False)
            weights = filter_fn(lambd)
            return (lambd,weights)
            #source = {'wavelengths': lambd, 'weights': weights}



    def _getSynphotBandpass(self):
        """ Return a pysynphot.ObsBandpass object for the given desired band
        This is split out as its own function so that the TFI subclass can override
        it with an etalon-aware version.
        """
        try:
            band = pysynphot.ObsBandpass( ('%s,im,%s'%(self.name, self.filter)).lower())
        except:
            # the requested band is not yet supported in synphot/CDBS. (those files are still a
            # work in progress...). Therefore, use our local throughput files and create a synphot
            # transmission object.
            wf = N.where(self.filter == N.asarray(self.filter_list))[0]
            # The existing FITS files all have wavelength in ANGSTROMS since that is the pysynphot convention...
            filterdata = atpy.Table(self._filter_files[wf], type='fits')

            _log.warn("Filter %s not supported in pysynphot/CDBS yet. Falling back to local filter transmission files" % self.filter)
            _log.warn("These may be less accurate.")
            band = pysynphot.spectrum.ArraySpectralElement(throughput=filterdata.THROUGHPUT,
                                wave=filterdata.WAVELENGTH, waveunits='angstrom',name=self.filter)
        return band



    def _getOpticalSystem(self,calc_oversample=2, detector_oversample = None, fov_arcsec=2, fov_pixels=None):
        """ Return an OpticalSystem instance corresponding to the instrument as currently configured.

        When creating such an OpticalSystem, you must specify the parameters needed to define the 
        desired sampling, specifically the oversampling and field of view. 


        Parameters
        ----------

        calc_oversample : int
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


        if detector_oversample is None: detector_oversample = calc_oversample

        _log.debug("Oversample: %d  %d " % (calc_oversample, detector_oversample))
        optsys = poppy.OpticalSystem(name='JWST+'+self.name, oversample=calc_oversample)
        if 'source_offset_r' in self.options.keys(): optsys.source_offset_r = self.options['source_offset_r']
        if 'source_offset_theta' in self.options.keys(): optsys.source_offset_theta = self.options['source_offset_theta']

        if isinstance(self.pupilopd, str):
            full_opd_path = self.pupilopd if os.path.exists( self.pupilopd) else os.path.join(self._datapath, "OPD",self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__'):
            full_opd_path =  (self.pupilopd[0] if os.path.exists( self.pupilopd[0]) else os.path.join(self._datapath, "OPD",self.pupilopd[0]), self.pupilopd[1])
        elif self.pupilopd is None: full_opd_path = None

        full_pupil_path = self.pupil if os.path.exists( self.pupil) else os.path.join(self._WebbPSF_basepath,self.pupil)
        optsys.addPupil(name='JWST Pupil', transmission=full_pupil_path, opd=full_opd_path, opdunits='micron', rotation=self._rotation)

        if self.image_mask is not None or self.pupil_mask is not None or ('force_coron' in self.options.keys() and self.options['force_coron']):
            optsys, trySAM, SAM_box_size = self._addCoronagraphOptics(optsys, oversample=calc_oversample)
        else: trySAM = False

        if fov_pixels is None:
            fov_pixels = N.round(fov_arcsec/self.pixelscale)
            if 'parity' in self.options.keys():
                if self.options['parity'].lower() == 'odd'  and N.remainder(fov_pixels,2)==0: fov_pixels +=1
                if self.options['parity'].lower() == 'even' and N.remainder(fov_pixels,2)==1: fov_pixels +=1

        optsys.addDetector(self.pixelscale, fov_pixels = fov_pixels, oversample = detector_oversample, name=self.name+" detector")

        if trySAM and not ('no_sam' in self.options.keys() and self.options['no_sam']): # if this flag is set, try switching to SemiAnalyticCoronagraph mode. 
            try: 
                SAM_optsys = poppy.SemiAnalyticCoronagraph(optsys, oversample=calc_oversample, occulter_box = SAM_box_size )
                return SAM_optsys
            except: 
                pass
 


        return optsys


    def display(self):
        """Display the currently configured optical system on screen """
        optsys = self._getOpticalSystem()
        optsys.display()

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

        self.apertures =  [ {'name': 'Imager', 'size': (768,1024), 'avail_filt': [f for f in self.filter_list if 'C' in f]},
                {'name': 'Cor-1065', 'size': (256,256), 'avail_filt': ['F1065C']},
                {'name': 'Cor-1140', 'size': (256,256), 'avail_filt': ['F1140C']},
                {'name': 'Cor-1150', 'size': (256,256), 'avail_filt': ['F1550C']},
                {'name': 'Cor-2300', 'size': (256,256), 'avail_filt': ['F2300C']},
                {'name': 'MRS', 'size': (100,100), 'avail_filt': 'IFU'}]


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
        """Add coronagraphic optics for MIRI 
        """

        # For MIRI coronagraphy, all the coronagraphic optics are rotated the same
        # angle as the instrument is, relative to the primary. So they see the unrotated
        # telescope pupil.
        # We model this by just not rotating till after the coronagraph. Thus we need to
        # un-rotate the primary that was already created in _getOpticalSystem.

        defaultpupil = optsys.planes.pop() # throw away the rotated pupil we just previously added
        optsys.addPupil(name='JWST Pupil', transmission=defaultpupil.amplitude_file, opd=defaultpupil.opd_file, opdunits='micron', rotation=None)
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

            #optsys.addImage(poppy.IdealFQPM(wavelength=10.65e-6, name=self.image_mask) )
        elif self.image_mask == 'FQPM1140':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1140",
                opticslist = [  poppy.IdealFQPM(wavelength=11.40e-6, name=self.image_mask),
                                poppy.IdealFieldStop(size=24, angle=-4.56)])
            optsys.addImage(container)
        elif self.image_mask == 'FQPM1550':
            container = poppy.CompoundAnalyticOptic(name = "MIRI FQPM 1550",
                opticslist = [  poppy.IdealFQPM(wavelength=15.50e-6, name=self.image_mask),
                                poppy.IdealFieldStop(size=24, angle=-4.56)])
            optsys.addImage(container)
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
        else: 
            optsys.addImage()

        if (self.image_mask is not None and 'FQPM' in self.image_mask)  or 'force_fqpm_shift' in self.options.keys() : optsys.addPupil("FQPM FFT aligner", direction='backward')

        # add pupil plane mask
        if ('pupil_shift_x' in self.options.keys() and self.options['pupil_shift_x'] != 0) or \
           ('pupil_shift_y' in self.options.keys() and self.options['pupil_shift_y'] != 0):
            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        else: shift = None

        #optsys.addPupil('Circle', radius=6.5/2)

        if self.pupil_mask == 'MASKFQPM':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MIRI_FQPMLyotStop.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKLYOT':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MIRI_LyotLyotStop.fits", name=self.pupil_mask, shift=shift)
        else: # all the MIRI filters have a tricontagon outline, even the non-coron ones.
            optsys.addPupil(transmission=self._WebbPSF_basepath+"/tricontagon.fits", name = 'filter cold stop', shift=shift)

        optsys.addRotation(self._rotation)

        return (optsys, False, 0) # don't ever try SemiAnalyticCoronagraph for MIRI? 
            # TBD maybe for Lyot mode, later?


class NIRCam(JWInstrument):
    """ A class modeling the optics of NIRCam. 
    
    Relevant attributes include `filter`, `image_mask`, and `pupil_mask`.

    The NIRCam class is smart enough to automatically select the appropriate pixel scale for the short or long wavelength channel 
    based on whether you request a short or long wavelength filter.
 
    """
    def __init__(self):
        self.pixelscale = 0.0317 # for short-wavelen channels
        self._pixelscale_short = 0.0317 # for short-wavelen channels
        self._pixelscale_long = 0.0648 # for short-wavelen channels
        JWInstrument.__init__(self, "NIRCam") # do this after setting the long & short scales.
        self.pixelscale = 0.0317 # need to redo 'cause the __init__ call will reset it to zero.

        #self.image_mask_list = ['BLC2100','BLC3350','BLC4300','WEDGESW','WEDGELW']
        self.image_mask_list = ['MASKLWB','MASKSWB','MASK210R','MASK335R','MASK430R']

        self.pupil_mask_list = ['CIRCLYOT','WEDGELYOT']

        self.filter = 'F200W' # default

        self.apertures = [
            {'name': 'Imager-SW A', 'size': (2048,2048), 'avail_filt': self.filter_list}, 
            {'name': 'Imager-SW B', 'size': (2048,2048), 'avail_filt': self.filter_list}, 
            {'name': 'Imager-LW A', 'size': (2048,2048), 'avail_filt': self.filter_list}, 
            {'name': 'Imager-LW B', 'size': (2048,2048), 'avail_filt': self.filter_list}, 
            {'name': 'Coron-BLC 2.1', 'size': (256,256), 'avail_filt': self.filter_list}, 
            {'name': 'Coron-BLC 2.1', 'size': (256,256), 'avail_filt': self.filter_list}]


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
            trySAM = True
            SAM_box_size = [5,20]
        elif self.image_mask == 'MASKLWB':
            optsys.addImage(function='BandLimitedCoron', kind='nircamwedge', wavelength=4.6e-6, name=self.image_mask)
            trySAM = True
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
        elif (self.pupil_mask is None and self.image_mask is not None):
            optsys.addPupil(name='No Lyot Mask Selected!')

        return (optsys, trySAM, SAM_box_size) # don't ever try SemiAnalyticCoronagraph for NIRCam? 


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

    def _validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('NIRSpec does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
    def _addCoronagraphOptics(self,optsys, oversample=2):
        raise NotImplementedError("No Coronagraph in NIRSpec!")


    def _validate_config(self):
        if self.filter.startswith("IFU"): raise NotImplementedError("The NIRSpec IFU is not yet implemented.")


class TFI(JWInstrument):
    """ A class modeling the optics of the Tunable Filter Imager
    
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
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC21N.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKC66N':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC66N.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKC71N':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC71N.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASK_NRM':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASK_NRM.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'CLEAR':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKCLEAR.fits", name=self.pupil_mask, shift=shift)
        elif (self.pupil_mask  is None and self.image_mask is not None):
            optsys.addPupil(name='No Lyot Mask Selected!')

        return (optsys, trySAM, radius+0.05) # always attempt to cast this to a SemiAnalyticCoronagraph

    def _getSynphotBandpass(self):
        """ Return a pysynphot.ObsBandpass object for the given desired band.
        This uses a lookup table to predict the properties of the TFI tunable filter etalon.
        """

        
        if (self.etalon_wavelength < 1.5 or self.etalon_wavelength > 5.0 or
            (self.etalon_wavelength > 2.7 and self.etalon_wavelength < 3.0)):
            raise ValueError("Invalid value for etalon wavelength: %f. Please set a value in 1.5-2.7 or 3.0-5.0 microns." % self.etalon_wavelength)


        match_index = N.abs(self.resolution_table.wavelength - self.etalon_wavelength).argmin()
        resolution = self.resolution_table.resolution[match_index]
        _log.info("Etalon wavelength %.3f has resolution %.2f" % (self.etalon_wavelength, resolution))
        wavelen = N.linspace(1.0, 5.0, 1000)

        fwhm = self.etalon_wavelength/resolution
        sigma = fwhm / 2.35482

        transmission = N.exp(- (wavelen - self.etalon_wavelength)**2/ (2*sigma**2))

        #P.plot(wavelen, transmission)
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

    def _validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('FGS does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
        #TODO only one possible filter fot the FGS, too. 
    def _addCoronagraphOptics(self,optsys):
        raise NotImplementedError("No Coronagraph in FGS!")


###########################################################################
# Generic utility functions

def Instrument(name):
    """This is just a convenience function, allowing one to access instrument objects based on a string.
    For instance, 

    >>> t = Instrument('TFI')


    
    Parameters
    ----------
    name : string
        Name of the instrument class to return. Case insensitive.
    
    """
    name = name.lower()
    if name == 'miri': return MIRI()
    if name == 'nircam': return NIRCam()
    if name == 'nirspec': return NIRSpec()
    if name == 'tfi': return TFI()
    if name == 'fgs': return FGS()
    else: raise ValueError("Incorrect instrument name "+name)
Instrument.list = ['nircam', 'nirspec', 'tfi', 'miri'] # useful list for iteration


def calc_or_load_psf(filename, inst, clobber=False, **kwargs):
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
                  diameter=None, oversample=4, type=N.float64,
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








###########################################################################
#
#    Display and PSF evaluation functions follow..
#
def display_PSF(HDUlist_or_filename=None, ext=0,
    vmin=1e-8,vmax=1e-1, scale='log', cmap = matplotlib.cm.jet, 
        title=None, imagecrop=None, adjust_for_oversampling=False, normalize='None', crosshairs=False, markcentroid=False, colorbar=True, colorbar_orientation='vertical',
        pixelscale='PIXELSCL'):
    """Display nicely a PSF from a given HDUlist or filename 

    This is extensively configurable. In addition to making an attractive display, for
    interactive usage this function provides a live display of the pixel value at a
    given (x,y) as you mouse around the image. 
    
    Parameters
    ----------
    HDUlist_or_filename : pyfits.HDUlist or string
        FITS file containing image to display.
    ext : int
        FITS extension. default = 0
    vmin, vmax : float
        min and max for image display scaling
    scale : str
        'linear' or 'log', default is log
    cmap : matplotlib.cm.Colormap instance
        Colormap to use. Default is matplotlib.cm.jet
    title : string, optional
    imagecrop : float
        size of region to display (default is whole image)
    adjust_for_oversampling : bool
        rescale to conserve surface brightness for oversampled PSFs? 
        (making this True conserves surface brightness but not total flux)
        default is False, to conserve total flux.
    markcentroid : bool
        Draw a crosshairs at the image centroid location?
        Centroiding is computed with the JWST-standard moving box algorithm.
    colorbar : bool
        Draw a colorbar?
    colorbar_orientation : str
        either 'horizontal' or 'vertical'; default is vertical.
    pixelscale : str or float
        if str, interpreted as the FITS keyword name for the pixel scale in arcsec/pixels.
        if float, used as the pixelscale directly.



    """
    if isinstance(HDUlist_or_filename, str):
        HDUlist = pyfits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, pyfits.HDUList):
        HDUlist = HDUlist_or_filename
    else: raise ValueError("input must be a filename or HDUlist")

    if adjust_for_oversampling:

        try:
            scalefactor = HDUlist[ext].header['OVERSAMP']**2
        except:
            _log.error("Could not determine oversampling scale factor; therefore NOT rescaling fluxes.")
            scalefactor=1
        im = HDUlist[ext].data *scalefactor
    else: im = HDUlist[ext].data

    if normalize.lower() == 'peak':
        _log.debug("Displaying image normalized to peak = 1")
        im /= im.max()
    elif normalize.lower() =='total':
        _log.debug("Displaying image normalized to PSF total = 1")
        im /= im.sum()
 

    if scale == 'linear':
        norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    else: 
        norm=LogNorm(vmin=vmin, vmax=vmax)

    if type(pixelscale) is str:
        halffov = HDUlist[ext].header[pixelscale]*HDUlist[ext].data.shape[0]/2
    else:
        try: 
            pixelscale = float(pixelscale)
        except:
            _log.warning("Provided pixelscale is neither float nor str; cannot use it. Using default=1 instead.")
            pixelscale = 1.0
        halffov = pixelscale*HDUlist[ext].data.shape[0]/2
    unit="arcsec"
    extent = [-halffov, halffov, -halffov, halffov]


    ax = poppy.imshow_with_mouseover( im   ,extent=extent,cmap=cmap, norm=norm)
    if imagecrop is not None:
        halffov = min( (imagecrop/2, halffov))
    ax.set_xbound(-halffov, halffov)
    ax.set_ybound(-halffov, halffov)
    if crosshairs: 
        ax.axhline(0,ls=":", color='k')
        ax.axvline(0,ls=":", color='k')


    if title is None:
        try:
            fspec = "%s, %s" % (HDUlist[ext].header['INSTRUME'], HDUlist[ext].header['FILTER'])
        except: 
            fspec= str(HDUlist_or_filename)
        title="PSF sim for "+fspec
    P.title(title)

    if colorbar:
        cb = P.colorbar(ax.images[0], orientation=colorbar_orientation)
        if scale.lower() =='log':
            ticks = N.logspace(N.log10(vmin), N.log10(vmax), N.log10(vmax/vmin)+1)
            if colorbar_orientation=='horizontal' and vmax==1e-1 and vmin==1e-8: ticks = [1e-8, 1e-6, 1e-4,  1e-2, 1e-1] # looks better
            cb.set_ticks(ticks)
            cb.set_ticklabels(ticks)
        if normalize =='peak':
            cb.set_label('Intensity relative to peak pixel')
        else: 
            cb.set_label('Fractional intensity per pixel')

    if markcentroid:
        _log.info("measuring centroid to mark on plot...")
        ceny, cenx = measure_centroid(HDUlist, ext=ext, units='arcsec', relativeto='center', boxsize=20, threshhold=0.1)
        ax.plot(cenx, ceny, 'k+', markersize=15, markeredgewidth=1)
        _log.info("centroid: (%f, %f) " % (cenx, ceny))
        P.draw()


def display_PSF_difference(HDUlist_or_filename1=None, HDUlist_or_filename2=None, ext1=0, ext2=0, vmax=1e-4, title=None, imagecrop=None, adjust_for_oversampling=False, normalize=False, crosshairs=False, colorbar=True, colorbar_orientation='vertical', print_=False):
    """Display nicely the difference of two PSFs from given files 
    
    Parameters
    ----------
    HDUlist_or_filename1,2 : pyfits.HDUlist or string
        FITS files containing image to difference
    ext1, ext2 : int
        FITS extension. default = 0
    vmax : float
        for the  scaling
    title : string, optional
    imagecrop : float
        size of region to display (default is whole image)
    normalize : bool
        Display (difference image)/(mean image) instead of just the difference image
    adjust_for_oversampling : bool
        rescale to conserve surface brightness for oversampled PSFs? 
        (making this True conserves surface brightness but not total flux)
        default is False, to conserve total flux.
    """
    if isinstance(HDUlist_or_filename1, str):
        HDUlist1 = pyfits.open(HDUlist_or_filename1)
    elif isinstance(HDUlist_or_filename1, pyfits.HDUList):
        HDUlist1 = HDUlist_or_filename1
    else: raise ValueError("input must be a filename or HDUlist")
    if isinstance(HDUlist_or_filename2, str):
        HDUlist2 = pyfits.open(HDUlist_or_filename2)
    elif isinstance(HDUlist_or_filename2, pyfits.HDUList):
        HDUlist2 = HDUlist_or_filename2
    else: raise ValueError("input must be a filename or HDUlist")


    if adjust_for_oversampling:
        scalefactor = HDUlist1[ext1].header['OVERSAMP']**2
        im1 = HDUlist1[ext1].data *scalefactor
        scalefactor = HDUlist2[ext2].header['OVERSAMP']**2
        im2 = HDUlist1[ext2].data *scalefactor
    else: 
        im1 = HDUlist1[ext1].data
        im2 = HDUlist2[ext2].data

    diff_im = im1-im2

    if normalize:
        avg_im = (im1+im2)/2
        diff_im /= avg_im
        cbtitle = 'Image difference / average  (per pixel)' #Relative intensity difference per pixel'
    else:
        cbtitle = 'Intensity difference per pixel'


    if print_:
        rms_diff = N.sqrt((diff_im**2).mean())
        print "RMS of difference image: %f" % rms_diff

    norm=matplotlib.colors.Normalize(vmin=-vmax, vmax=vmax)
    cmap = matplotlib.cm.gray
    halffov = HDUlist1[ext1].header['PIXELSCL']*HDUlist1[ext1].data.shape[0]/2
    unit="arcsec"
    extent = [-halffov, halffov, -halffov, halffov]


    ax = poppy.imshow_with_mouseover( diff_im   ,extent=extent,cmap=cmap, norm=norm)
    if imagecrop is not None:
        halffov = min( (imagecrop/2, halffov))
    ax.set_xbound(-halffov, halffov)
    ax.set_ybound(-halffov, halffov)
    if crosshairs: 
        ax.axhline(0,ls=":", color='k')
        ax.axvline(0,ls=":", color='k')


    if title is None:
        try:
            fspec= str(HDUlist_or_filename1) +"-"+str(HDUlist_or_filename2)
            #fspec = "Difference Image " # "%s, %s" % (HDUlist[ext].header['INSTRUME'], HDUlist[ext].header['FILTER'])
        except: 
            fspec= str(HDUlist_or_filename1) +"-"+str(HDUlist_or_filename2)
        title="Difference of "+fspec
    P.title(title)

    if colorbar:
        cb = P.colorbar(ax.images[0], orientation=colorbar_orientation)
        #ticks = N.logspace(N.log10(vmin), N.log10(vmax), N.log10(vmax/vmin)+1)
        #if vmin == 1e-8 and vmax==1e-1: 
            #ticks = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
        ticks = [-vmax, -0.5*vmax, 0, 0.5*vmax, vmax]
        cb.set_ticks(ticks)
        cb.set_ticklabels(ticks)
        #stop()
        cb.set_label(cbtitle)



def display_profiles(HDUlist_or_filename=None,ext=0, overplot=False ):
    if isinstance(HDUlist_or_filename, str):
        HDUlist = pyfits.open(filename,ext=ext)
    elif isinstance(HDUlist_or_filename, pyfits.HDUList):
        HDUlist = HDUlist_or_filename
    else: raise ValueError("input must be a filename or HDUlist")


    radius, profile, EE = radial_profile(HDUlist, EE=True)

    if not overplot:
        P.clf()
        P.title("PSF sim for %s, %s" % (HDUlist[ext].header['INSTRUME'], HDUlist[ext].header['FILTER']))
        P.xlabel("Radius [arcsec]")
        P.ylabel("PSF radial profile")
    P.subplot(2,1,1)
    P.semilogy(radius, profile)

    fwhm = 2*radius[N.where(profile < profile[0]*0.5)[0][0]]
    P.text(fwhm, profile[0]*0.5, 'FWHM = %.3f"' % fwhm)

    P.subplot(2,1,2)
    #P.semilogy(radius, EE, nonposy='clip')
    P.plot(radius, EE, color='r') #, nonposy='clip')
    if not overplot:
        P.xlabel("Radius [arcsec]")
        P.ylabel("Encircled Energy")

    for level in [0.5, 0.8, 0.95]:
        EElev = radius[N.where(EE > level)[0][0]]
        yoffset = 0 if level < 0.9 else -0.05 
        P.text(EElev+0.1, level+yoffset, 'EE=%2d%% at r=%.3f"' % (level*100, EElev))


def radial_profile(HDUlist_or_filename=None, ext=0, EE=False, center=None, stddev=False, binsize=None, maxradius=None):
    """ Compute a radial profile of the image. 

    This computes a discrete radial profile evaluated on the provided binsize. For a version
    interpolated onto a continuous curve, see measure_radial().

    Code taken pretty much directly from pydatatut.pdf

    Parameters
    ----------
    HDUlist_or_filename : string
        what it sounds like.
    ext : int
        Extension in FITS file
    EE : bool
        Also return encircled energy (EE) curve in addition to radial profile?
    center : tuple of floats
        Coordinates (x,y) of PSF center. Default is image center. 
    binsize : float
        size of step for profile. Default is pixel size.
    stddev : bool
        Compute standard deviation in each radial bin, not average?


    Returns
    --------
    results : tuple
        Tuple containing (radius, profile) or (radius, profile, EE) depending on what is requested.
        The radius gives the center radius of each bin, while the EE is given inside the whole bin
        so you should use (radius+binsize/2) for the radius of the EE curve if you want to be
        as precise as possible.
    """
    if isinstance(HDUlist_or_filename, str):
        HDUlist = pyfits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, pyfits.HDUList):
        HDUlist = HDUlist_or_filename
    else: raise ValueError("input must be a filename or HDUlist")

    image = HDUlist[ext].data
    pixelscale = HDUlist[ext].header['PIXELSCL']


    if maxradius is not None:
        raise NotImplemented("add max radius")
        stop()


    if binsize is None:
        binsize=pixelscale

    y,x = N.indices(image.shape)
    if center is None:
        # get exact center of image
        #center = (image.shape[1]/2, image.shape[0]/2)
        center = tuple( (a-1)/2.0 for a in image.shape[::-1])

    r = N.sqrt( (x-center[0])**2 + (y-center[1])**2) *pixelscale / binsize # radius in bin size steps
    ind = N.argsort(r.flat)


    sr = r.flat[ind]
    sim = image.flat[ind]
    ri = sr.astype(int)
    deltar = ri[1:]-ri[:-1] # assume all radii represented (more work if not)
    rind = N.where(deltar)[0]
    nr = rind[1:] - rind[:-1] # number in radius bin
    csim = N.cumsum(sim, dtype=float) # cumulative sum to figure out sums for each bin
    tbin = csim[rind[1:]] - csim[rind[:-1]] # sum for image values in radius bins
    radialprofile=tbin/nr

    #pre-pend the initial element that the above code misses.
    radialprofile2 = N.empty(len(radialprofile)+1)
    if rind[0] != 0:
        radialprofile2[0] =  csim[rind[0]] / (rind[0]+1)  # if there are multiple elements in the center bin, average them
    else:
        radialprofile2[0] = csim[0]                       # otherwise if there's just one then just take it. 
    radialprofile2[1:] = radialprofile
    rr = N.arange(len(radialprofile2))*binsize + binsize*0.5  # these should be centered in the bins, so add a half.

    if stddev:
        stddevs = N.zeros_like(radialprofile2)
        r_pix = r * binsize
        for i, radius in enumerate(rr):
            if i == 0: wg = N.where(r < radius+ binsize/2)
            else: 
                wg = N.where( (r_pix >= (radius-binsize/2)) &  (r_pix < (radius+binsize/2)))
                #print radius-binsize/2, radius+binsize/2, len(wg[0])
                #wg = N.where( (r >= rr[i-1]) &  (r <rr[i] )))
            stddevs[i] = image[wg].std()
        #stop()
        return (rr, stddevs)

    if not EE:
        return (rr, radialprofile2)
    else:
        #weighted_profile = radialprofile2*2*N.pi*(rr/rr[1])
        #EE = N.cumsum(weighted_profile)
        EE = csim[rind]
        return (rr, radialprofile2, EE) 


def measure_EE(HDUlist_or_filename=None, ext=0, center=None, binsize=None):
    """ Returns a function object which when called returns the Encircled Energy inside a given radius.



    Parameters
    ----------
    HDUlist_or_filename : string
        what it sounds like.
    ext : int
        Extension in FITS file
    center : tuple of floats
        Coordinates (x,y) of PSF center. Default is image center. 
    binsize: 
        size of step for profile. Default is pixel size.

    Returns
    --------
    encircled_energy: function
        A function which will return the encircled energy interpolated to any desired radius.


    Examples
    --------
    >>> EE = measure_EE("someimage.fits")
    >>> print "The EE at 0.5 arcsec is ", EE(0.5)

    """

    rr, radialprofile2, EE = radial_profile(HDUlist_or_filename, ext, EE=True, center=center, binsize=binsize)

    # append the zero at the center
    rr_EE = rr + (rr[1]-rr[0])/1  # add half a binsize to this, because the EE is measured inside the
                                  # outer edge of each annulus. 
    rr0 = N.concatenate( ([0], rr_EE)) 
    EE0 = N.concatenate( ([0], EE))


    EE_fn = scipy.interpolate.interp1d(rr0, EE0,kind='cubic', bounds_error=False)

    return EE_fn
    

def measure_radial(HDUlist_or_filename=None, ext=0, center=None, binsize=None):
    """ Returns a function object which when called returns the mean value at a given radius.

    Parameters
    ----------
    HDUlist_or_filename : string
        what it sounds like.
    ext : int
        Extension in FITS file
    center : tuple of floats
        Coordinates (x,y) of PSF center. Default is image center. 
    binsize: 
        size of step for profile. Default is pixel size.

    Returns
    --------
    radial_profile: function
        A function which will return the mean PSF value at any desired radius.


    Examples
    --------
    >>> rp = measure_radial("someimage.fits")
    >>> radius = N.linspace(0, 5.0, 100)
    >>> plot(radius, rp(radius), label="PSF")

    """

    rr, radialprofile, EE = radial_profile(HDUlist_or_filename, ext, EE=True, center=center, binsize=binsize)

    radial_fn = scipy.interpolate.interp1d(rr, radialprofile,kind='cubic', bounds_error=False)

    return radial_fn
    

def measure_fwhm(HDUlist_or_filename=None, ext=0, center=None, level=0.5):
    """ Measure FWHM* by interpolation of the radial profile 
    (* or full width at some other fraction of max.)

    Parameters
    ----------
    HDUlist_or_filename, ext : string, int
        Same as above
    center : tuple
        center to compute around.  Default is image center.
    level : float
        Fraction of max to compute for; default is 0.5 for Half Max. 
        You can also measure widths at other levels e.g. FW at 10% max
        by setting level=0.1


    """

    rr, radialprofile, EE = radial_profile(HDUlist_or_filename, ext, EE=True, center=center)
    rpmax = radialprofile.max()

    wlower = N.where(radialprofile < rpmax *level)
    wmin = N.min(wlower[0])
    # go just a bit beyond the half way mark
    winterp = N.arange(0, wmin+2, dtype=int)[::-1]

    if len(winterp) < 6: kind='linear'
    else: kind = 'cubic'

    interp_hw = scipy.interpolate.interp1d(radialprofile[winterp], rr[winterp], kind=kind)
    return 2*interp_hw(rpmax*level)
 

def measure_sharpness(HDUlist_or_filename=None, ext=0):
    """ Compute image sharpness, the sum of pixel squares.

    See Makidon et al. JWST-STScI-001157 for a discussion of this image metric
    and its relationship to noise equivalent pixels.

    Parameters
    ----------
    HDUlist_or_filename, ext : string, int
        Same as above
 
    """
    if isinstance(HDUlist_or_filename, str):
        HDUlist = pyfits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, pyfits.HDUList):
        HDUlist = HDUlist_or_filename
    else: raise ValueError("input must be a filename or HDUlist")


    # TODO or check that the oversampling factor is 1
    try:
        detpixels = HDUlist['DET_SAMP']
    except:
        raise ValueError("You can only measure sharpness for an image with an extension giving the rebinned actual detector pixel values.""")

    sharpness =  (detpixels.data**2).sum()
    return sharpness

def measure_centroid(HDUlist_or_filename=None, ext=0, slice=0, boxsize=50, print_=False, units='pixels', relativeto='origin', **kwargs):
    """ Measure the center of an image via center-of-mass

    Parameters
    ----------
    HDUlist_or_filename, ext : string, int
        Same as above
    boxsize : int
        Half box size for centroid

    relativeto : string
        either 'origin' for relative to pixel (0,0) or 'center' for relative to image center. Default is 'origin'
    units : string
        either 'pixels' for position in pixels or 'arcsec' for arcseconds. 
        Relative to the relativeto parameter point in either case.
 

    Returns
    -------
    CoM : array_like
        [Y, X] coordinates of center of mass.

    """
    if isinstance(HDUlist_or_filename, str):
        HDUlist = pyfits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, pyfits.HDUList):
        HDUlist = HDUlist_or_filename
    else: raise ValueError("input must be a filename or HDUlist")

    image = HDUlist[ext].data
    
    if image.ndim >=3:  # handle datacubes gracefully
        image = image[slice,:,:]


    if 0: 
        y, x= N.indices(image.shape)
        wpeak = N.where(image == image.max())
        cy, cx = y[wpeak][0], x[wpeak][0]
        print "Peak pixel: (%d, %d)" % (cx, cy)


        cutout = image[cy-boxsize:cy+boxsize+1, cx-boxsize:cx+boxsize+1]
        cent_of_mass_cutout = N.asarray(scipy.ndimage.center_of_mass(cutout))
        cent_of_mass =  cent_of_mass_cutout + N.array([cy-boxsize, cx-boxsize])
    else:
        cent_of_mass = fwcentroid(image, halfwidth=boxsize, **kwargs)

    if print_: print("Center of mass: (%.4f, %.4f)" % (cent_of_mass[1], cent_of_mass[0]))

    if relativeto == 'center':
        imcen = N.array([ (image.shape[0]-1)/2., (image.shape[1]-1)/2. ])
        cent_of_mass  = tuple( N.array(cent_of_mass) -  imcen)


    if units == 'arcsec':
        pixelscale = HDUlist[ext].header['PIXELSCL']
        cent_of_mass = tuple( N.array(cent_of_mass) *pixelscale)

    return cent_of_mass


def measure_strehl(HDUlist_or_filename=None, ext=0, center=None, display=True, print_=True):
    """ Estimate the Strehl ratio for a PSF.
    
    This requires computing a simulated PSF with the same
    properties as the one under analysis.

    Note that this calculation will not be very accurate unless both PSFs are well sampled,
    preferably several times better than Nyquist. See 
    `Roberts et al. 2004 SPIE 5490 <http://adsabs.harvard.edu/abs/2004SPIE.5490..504R>`_
    for a discussion of the various possible pitfalls when calculating Strehl ratios. 

    Parameters
    ----------
    HDUlist_or_filename, ext : string, int
        Same as above

    center : tuple
        center to compute around.  Default is image center. If the center is on the
        crosshairs between four pixels, then the mean of those four pixels is used.
        Otherwise, if the center is in a single pixel, then that pixel is used. 

    print_, display : bool
        control whether to print the results or display plots on screen. 


    Returns
    ---------
    strehl : float
        Strehl ratio as a floating point number between 0.0 - 1.0
  
    """
    if isinstance(HDUlist_or_filename, str):
        HDUlist = pyfits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, pyfits.HDUList):
        HDUlist = HDUlist_or_filename
    else: raise ValueError("input must be a filename or HDUlist")

    image = HDUlist[ext].data
    header = HDUlist[ext].header
 
    if center is None:
        # get exact center of image
        #center = (image.shape[1]/2, image.shape[0]/2)
        center = tuple( (a-1)/2.0 for a in image.shape[::-1])



    # Compute a comparison image
    _log.info("Now computing image with zero OPD for comparison...")
    inst = Instrument(header['INSTRUME'])
    inst.filter = header['FILTER']
    inst.pupilopd = None # perfect image
    inst.pixelscale = header['PIXELSCL'] * header['OVERSAMP'] # same pixel scale pre-oversampling
    comparison_psf = inst.calcPSF(fov_arcsec = header['FOV'], oversample=header['OVERSAMP'], nlambda=header['NWAVES'])
    comparison_image = comparison_psf[0].data

    
    if (int(center[1]) == center[1]) and (int(center[0]) == center[0]):
        # individual pixel
        meas_peak =           image[center[1], center[0]]
        ref_peak = comparison_image[center[1], center[0]]
    else:
        # average across a group of 4
        bot = [N.floor(f) for f in center]
        top = [N.ceil(f)+1 for f in center]
        meas_peak =           image[bot[1]:top[1], bot[0]:top[0]].mean()
        ref_peak = comparison_image[bot[1]:top[1], bot[0]:top[0]].mean()
    strehl = (meas_peak/ref_peak)

    if display:
        P.clf()
        P.subplot(121)
        display_PSF(HDUlist, title="Observed PSF")
        P.subplot(122)
        display_PSF(comparison_psf, title="Perfect PSF")
        P.gcf().suptitle("Strehl ratio = %.3f" % strehl) 


    if print_:

        print "Measured peak:  %.3g" % meas_peak
        print "Reference peak: %.3g" % ref_peak
        print "  Strehl ratio = %.3f " % strehl

    return strehl


def measure_anisotropy(HDUlist_or_filename=None, ext=0, slice=0, boxsize=50):
    pass


 

#########################3

def rebin_array(a = None, rc=(2,2), verbose=False):
	"""  
	Perform simple-minded flux-conserving binning... clip trailing
	size mismatch: eg a 10x3 array binned by 3 results in a 3x1 array

    Parameters
    ----------
    a : array_like
        input array
    rc : two-element tuple 
        (nrows, ncolumns) desired for rebinned array
    verbose : bool
        print additional status text?


	anand@stsci.edu

	"""

	r, c = rc

	R = a.shape[0]
	C = a.shape[1]

	nr = int(R / r)
	nc = int(C / c)

	b = a[0:nr, 0:nc].copy()
	b = b * 0

	for ri in range(0, nr):
		Rlo = ri * r
		if verbose:
			print "row loop"
		for ci in range(0, nc):
			Clo = ci * c
			b[ri, ci] = N.add.reduce(a[Rlo:Rlo+r, Clo:Clo+c].copy().flat)
			if verbose:
				print "    [%d:%d, %d:%d]" % (Rlo,Rlo+r, Clo,Clo+c),
				print "%4.0f"  %   N.add.reduce(a[Rlo:Rlo+r, Clo:Clo+c].copy().flat)
	return b


def makeFakeFilter(filename, lcenter, dlam, instrument='NIRCam', name='filter', source='Fake top hat filter', clobber=False):
    """ arguments in microns, but file written in angstroms """

    lstart = lcenter - dlam/2
    lstop = lcenter + dlam/2

    nlambda = 40

    _log.info("Filter from %f - %f " % (lstart, lstop))
    wavelength = N.linspace( lstart-dlam*0.1, lstop+dlam*0.1, nlambda)
    _log.debug(wavelength)
    transmission = N.zeros_like(wavelength)
    transmission[N.where( (wavelength > lstart) & (wavelength < lstop) )] = 1.0

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


    nc = NIRCam()
    nc.filter = 'F460M'
    nc.image_mask = 'MASK430R'
    nc.pupil_mask = 'CIRCLYOT'
    #nc.calcPSF('test_nircam.fits', mono=False)

    miri=MIRI()
    miri.image_mask = 'FQPM1065'
    miri.pupil_mask = 'MASKFQPM'
    miri.filter='F1065C'

    #miri.display()
    nircam=NIRCam()
    tfi = TFI()
    tfi.image_mask = "CORON058"
    tfi.pupil_mask = 'MASKC66N'
    nirspec = NIRSpec()
