#!/usr/bin/env python
"""

=====
JWOPT
=====

An object-oriented modeling system for the JWST instruments.

Classes:
  * JWInstrument
    * MIRI
    * NIRCam
    * NIRSpec
    * TFI
    * FGS

  * kurucz_stars  (a wrapper for the Kurucz spectra library)



    Code by Marshall Perrin <mperrin@stsci.edu>

"""

import numpy as N
import poppy
from poppy import *
import os
import glob
import atpy
import types
import utils
import time


try:
    __IPYTHON__
    from IPython.Debugger import Tracer; stop = Tracer()
except:
    pass




class JWInstrument(object):
    """ A generic JWST Instrument class.
    
    *Note*: Do not use this class directly - instead use one of the instrument subclasses!
    
    This class provides a simple interface for modeling PSF formation through the JWST instruments, 
    with configuration options and software interface loosely resembling the configuration of the instrument 
    mechanisms.   For instance, one can create an instance of MIRI and configure it for coronagraphic observations, 
    and then compute a PSF:

    >>> miri = MIRI()
    >>> miri.filter = 'F1065C'
    >>> miri.image_mask = 'FQPM1065'
    >>> miri.pupil_mask = 'MASKFQPM'
    >>> miri.calcPSF('outfile.fits')

    Note that the interface currently only provides a modicum of error checking, and relies on the user
    being knowledgable enough to avoid trying to simulate some physically impossible or just plain silly
    configuration (such as trying to use a FQPM with the wrong filter).


    The instrument constructors do not take any arguments. Instead, create an instrument object and then
    configure the filter or other attributes as desired. 


    """

    def __init__(self, name=None):
        self.name=name

        self._JWPSF_basepath = os.path.dirname(os.path.dirname(os.path.abspath(poppy.__file__))) +os.sep+"data"

        self._datapath = self._JWPSF_basepath + os.sep + self.name + os.sep
        self.pupil = os.path.abspath(self._datapath+"../pupil_RevV.fits")
        "Filename for JWST pupil mask. Usually there is no need to change this."
        self.pupilopd = None   # This can optionally be set to a tuple indicating (filename, slice in datacube)
        "Filename for JWST pupil OPD. If the file contains a datacube, set this to a tuple (filename, slice) "

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
        downsample : bool
            For output files, write an additional FITS extension including a version of the output array 
            rebinned down to the actual detector pixel scale?
        jitter : string
            Type of jitter model to apply.
        """

        #create private instance variables. These will be
        # wrapped just below to create properties with validation.
        self._filter=None
        self._filter_files= [os.path.abspath(f) for f in glob.glob(self._datapath+os.sep+'filters/*_thru.fits')]
        self.filter_list=[os.path.basename(f).split("_")[0] for f in self._filter_files]
        "List of available filters"

        def sort_filters(filtname):
            try:
                return int(filtname[1:4])
            except:
                return filtname
        self.filter_list.sort(key=sort_filters)
        self._filter_files = [self._datapath+os.sep+'filters/'+f+"_thru.fits" for f in self.filter_list]

        self.filter = self.filter_list[0]


        #self.opd_list = [os.path.basename(os.path.abspath(f)) for f in glob.glob(self._datapath+os.sep+'OPD/*.fits')]
        self.opd_list = [os.path.basename(os.path.abspath(f)) for f in glob.glob(self._datapath+os.sep+'OPD/OPD*.fits')]
        #self.opd_list.insert(0,"Zero OPD (Perfect)")
        
        self._image_mask=None
        self.image_mask_list=[]
        "List of available image_masks"

        self._pupil_mask=None
        self.pupil_mask_list=[]
        "List of available pupil_masks"

        self.pixelscale = 0.0
        "Detector pixel scale, in arcsec/pixel"

    def _validate_config(self):
        pass

    # create properties with error checking
    @property
    def filter(self):
        'Currently selected filter name (e.g. "F200W")'
        return self._filter
    @filter.setter
    def filter(self, value):
        if value not in self.filter_list:
            raise ValueError("Instrument %s doesn't have a filter called %s." % (self.name, value))
        self._filter = value
        #self._validate_config()
    @property
    def image_mask(self):
        'Currently selected coronagraphic image plane mask, or None for direct imaging'
        return self._image_mask
    @image_mask.setter
    def image_mask(self, name):
        if name is "": name = None
        if name is not None:
            if name not in self.image_mask_list:
                raise ValueError("Instrument %s doesn't have an image mask called %s." % (self.name, value))
        self._image_mask = name
    @property
    def pupil_mask(self):
        'Currently selected Lyot pupil mask, or None for direct imaging'
        return self._pupil_mask
    @pupil_mask.setter
    def pupil_mask(self,name):
        if name is "": name = None
        if name is not None:
            if name not in self.pupil_mask_list:
                raise ValueError("Instrument %s doesn't have an pupil mask called %s." % (self.name, value))

        self._pupil_mask = name
            #self._validate_config()

    def __str__(self):
        return "JWInstrument name="+self.name

    #----- actual optical calculations follow here -----
    def calcPSF(self, filter=None, outfile=None, oversample=None, detector_oversample=None, calc_oversample=None, fov_arcsec=5., display=True, clobber=True, nlambda=5 ):
        """ Compute a PSF.
        The result can either be written to disk (set outfile="filename") or else will be returned as
        a pyfits HDUlist object.


        Output sampling may be specified in one of two ways: 
        1) Set `oversample=<number>`. This will use that oversampling factor beyond detector pixels
           for output images, and beyond Nyquist sampling for any FFTs to prior optical planes. 
        2) set `detector_oversample=<number>` and `calc_oversample=<other_number>`. This syntax lets
           you specify distinct oversampling factors for intermediate and final planes. 

        By default, both oversampling factors are set equal to 2.

        Parameters
        ----------
        outfile : string
            Filename to write. If None, then result is returned as an HDUList
        oversample : int
            How much to oversample. Default=2.
        fov_arcsec: float
            field of view in arcsec. Default=5
        clobber: bool
            overwrite output FITS file if it already exists?
        nlambda : int
            How many wavelengths to model for broadband? Default=5
        display : bool
            Whether to display the PSF when done or not.
        filter : string, optional
            Filter name. Setting this is just a shortcut for setting the object's filter first, then
            calling calcPSF afterwards.

        Returns
        -------
        outfits : pyfits.HDUList
            If `outfile` is set, the output is written to that file. Otherwise it is
            returned as a pyfits.HDUList object.


        """
        if filter is not None:
            self.filter = filter


        #if outfile is None: 
            #outfile = "PSF_%s_%s.fits" % (self.name, self.filter)
            #raise ValueError("You must specify an output file name.")


        # Implement the semi-convoluted logic for the oversampling options. See docstring above
        if oversample is not None and detector_oversample is not None and calc_oversample is not None:
            # all set -> complain!
            raise ValueError("You cannot specify simultaneously the oversample= option with the detector_oversample and calc_oversample options. Pick one or the other!")
        elif oversample is None and detector_oversample is None and calc_oversample is None:
            # nothing set -> set oversample = 2
            oversample = 2
        if detector_oversample is None: detector_oversample = oversample
        if calc_oversample is None: calc_oversample = oversample

        # instantiate an optical system using the current parameters
        self.optsys = self.getOpticalSystem(fov_arcsec=fov_arcsec, calc_oversample=calc_oversample, detector_oversample=detector_oversample)


        # compute a source spectrum weighted by the desired filter curves.
        # TBD this will eventually use pysynphot, so don't write anything fancy for now!
        wf = N.where(self.filter_list == self.filter)
        wf = N.where(self.filter == N.asarray(self.filter_list))[0]
        filterdata = atpy.Table(self._filter_files[wf])

        print "CAUTION: Really basic top-hat function for filter profile, with %d steps" % nlambda
        wtrans = N.where(filterdata.THROUGHPUT > 0.5)
        lrange = filterdata.WAVELENGTH[wtrans] *1e-10
        lambd = N.linspace(N.min(lrange), N.max(lrange), nlambda)
        weights = N.ones(nlambda)
        source = {'wavelengths': lambd, 'weights': weights}
        result = self.optsys.calcPSF(source, display_intermediates=True, save_intermediates=False, display=display)

        f = p.gcf()
        #p.text( 0.1, 0.95, "%s, filter= %s" % (self.name, self.filter), transform=f.transFigure, size='xx-large')
        p.suptitle( "%s, filter= %s" % (self.name, self.filter), size='xx-large')
        p.text( 0.7, 0.95, "Calculation with %d wavelengths (%g - %g um)" % (nlambda, lambd[0]*1e6, lambd[-1]*1e6), transform=f.transFigure)



        # TODO update FITS header here
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

        (year, month, day, hour, minute, second, weekday, DOY, DST) =  time.gmtime()
        result[0].header.update ("DATE", "%4d-%02d-%02dT%02d:%02d:%02d" % (year, month, day, hour, minute, second), "Date of calculation")
        result[0].header.update ("AUTHOR", "%s@%s" % (os.getenv('USER'), os.getenv('HOST')), "username@host for calculation")


        # Should we downsample? 
        if 'downsample' in self.options.keys() and self.options['downsample'] == True:
            print "** Downsampling to detector pixel scale."
            downsampled_result = result[0].copy()
            downsampled_result.data = utils.rebin(downsampled_result.data, rc=(detector_oversample, detector_oversample))
            downsampled_result.header.update('OVERSAMP', 1, 'These data are rebinned to detector pixels')
            downsampled_result.header.update('CALCSAMP', detector_oversample, 'This much oversampling used in calculation')
            downsampled_result.header.update('EXTNAME', 'DET_SAMP')
            downsampled_result.header['PIXELSCL'] *= detector_oversample
            result.append(downsampled_result)




        if outfile is not None:
            result[0].header.update ("FILENAME", os.path.basename (outfile),
                           comment="Name of this file")
            result.writeto(outfile, clobber=clobber)
            print "Saved result to "+outfile
        else:
            return result



    def getOpticalSystem(self,calc_oversample=2, detector_oversample = None, fov_arcsec=2):
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
        fov_arcsec : float
            Field of view, in arcseconds. Default is 2


        Returns
        -------
        osys : poppy.OpticalSystem 
            an optical system instance representing the desired configuration.

        """

        self._validate_config()

        print calc_oversample, detector_oversample
        optsys = OpticalSystem(name='JWST+'+self.name, oversample=calc_oversample)
        optsys.addPupil(name='JWST Pupil', transmission= self.pupil, opd=self.pupilopd, opdunits='micron')

        if self.image_mask is not None:
            optsys = self.addCoronagraphOptics(optsys)

        optsys.addDetector(self.pixelscale, fov_npix = fov_arcsec/self.pixelscale, oversample = detector_oversample, name=self.name+" detector")
        return optsys


    def display(self):
        """Display the currently configured optical system on screen """
        optsys = self.getOpticalSystem()
        optsys.display()

    def addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics to an optical system. 
        This method must be provided by derived instrument classes. 
        """
        raise NotImplementedError("needs to be subclassed.")


    def getFilter(self,filtername):
        """ Given a filter name, load the actual response curve and return it.  (depreciated??)"""
        if filtername not in self.filter_list:
            raise ValueError("Unknown/incorrect filter name for %s: %s" % (self.name, filtername))

        wm = N.where(N.asarray(self.filter_list) == filtername)
        filtfile = self._filter_files[wm[0]]
        print "Loading filter %s from %s" % (filtername, filtfile)
        t = atpy.Table(filtfile)

        t.WAVELENGTH /= 1e4 # convert from angstroms to microns
        return t


class MIRI(JWInstrument):
    """ A class modeling the optics of MIRI, the Mid-InfraRed Instrument"""
    def __init__(self):
        JWInstrument.__init__(self, "MIRI")
        self.pixelscale = 0.11

        self.pupilopd = None
        self.image_mask_list = ['FQPM1065', 'FQPM1140', 'FQPM1550', 'LYOT2300']
        self.pupil_mask_list = ['MASKFQPM', 'MASKLYOT']

        self.apertures =  [ {'name': 'Imager', 'size': (768,1024), 'avail_filt': [f for f in self.filter_list if 'C' in f]},
                {'name': 'Cor-1065', 'size': (256,256), 'avail_filt': ['F1065C']},
                {'name': 'Cor-1140', 'size': (256,256), 'avail_filt': ['F1140C']},
                {'name': 'Cor-1150', 'size': (256,256), 'avail_filt': ['F1550C']},
                {'name': 'Cor-2300', 'size': (256,256), 'avail_filt': ['F2300C']},
                {'name': 'MRS', 'size': (100,100), 'avail_filt': 'IFU'}]


    def _validate_config(self):
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
            elif self.filter == 'F1550C':
                assert self.image_mask == 'FQPM1550', 'Invalid configuration'
                assert self.pupil_mask == 'MASKFQPM', 'Invalid configuration'
            elif self.filter == 'F2300C':
                assert self.image_mask == 'LYOT2300', 'Invalid configuration'
                assert self.pupil_mask == 'MASKLYOT', 'Invalid configuration'
            else:
                raise ValueError("Invalid configuration selected!")


    def addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics for MIRI 
        """

        # Add image plane mask
        # For the MIRI FQPMs, we require the star to be centered not on the middle pixel, but
        # on the cross-hairs between four pixels. (Since that is where the FQPM itself is centered)
        # This is with respect to the intermediate calculation pixel scale, of course, not the
        # final detector pixel scale. 
        if self.image_mask == 'FQPM1065':
            optsys.addImage(function='FQPM',wavelength=10.65e-6, name=self.image_mask)
            optsys.addImage(function='fieldstop',size=24)
            optsys.source_position = N.array([0.5, 0.5]) * intermediate_pixel_scale
        elif self.image_mask == 'FQPM1140':
            optsys.addImage(function='FQPM',wavelength=11.40e-6, name=self.image_mask)
            optsys.addImage(function='fieldstop',size=24)
        elif self.image_mask == 'FQPM1550':
            optsys.addImage(function='FQPM',wavelength=15.50e-6, name=self.image_mask)
            optsys.addImage(function='fieldstop',size=24)
        elif self.image_mask == 'FQPM1550':
            optsys.addImage(function='CircularOcculter',radius =1., name=self.image_mask) 
            optsys.addImage(function='fieldstop',size=30)


        # add pupil plane mask
        if ('pupil_shift_x' in self.options.keys() and self.options['pupil_shift_x'] != 0) or \
           ('pupil_shift_y' in self.options.keys() and self.options['pupil_shift_y'] != 0):
            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        else: shift = None

        if self.pupil_mask == 'MASKFQPM':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MIRI_FQPMLyotStop.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKLYOT':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MIRI_LyotLyotStop.fits", name=self.pupil_mask, shift=shift)


        return optsys



class NIRCam(JWInstrument):
    """ A class modeling the optics of NIRCam"""
    def __init__(self):
        JWInstrument.__init__(self, "NIRCam")
        self.pixelscale = 0.0317 # for short-wavelen channels

        self.image_mask_list = ['BLC2100','BLC3350','BLC4300','WEDGESW','WEDGELW']
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
        pass

    def addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics for NIRCam
        """

        if self.image_mask is 'BLC2100':
            optsys.addImage(function='BandLimitedCoron', kind='circular', sigma=1, name=self.image_mask)
        if self.image_mask is 'BLC3350':
            optsys.addImage(function='BandLimitedCoron', kind='circular', sigma=1, name=self.image_mask)
        if self.image_mask is 'BLC4300':
            optsys.addImage(function='BandLimitedCoron', kind='circular', sigma=1, name=self.image_mask)
        elif self.image_mask is 'WEDGESW':
            optsys.addImage(function='BandLimitedCoron', kind='linear', sigma=1, name=self.image_mask)
        elif self.image_mask is 'WEDGELW':
            optsys.addImage(function='BandLimitedCoron', kind='linear', sigma=1, name=self.image_mask)

        # add pupil plane mask
        if ('pupil_shift_x' in self.options.keys() and self.options['pupil_shift_x'] != 0) or \
           ('pupil_shift_y' in self.options.keys() and self.options['pupil_shift_y'] != 0):
            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        else: shift = None


        if self.pupil_mask == 'MASKFQPM':
            optsys.addPupil(self._datapath+"/coronagraph/MIRI_FQPMLyotStop.fits", name=self.pupil_mask, shift=shift)
        elif self.pupil_mask == 'MASKLYOT':
            optsys.addPupil(self._datapath+"/coronagraph/MIRI_LyotLyotStop.fits", name=self.pupil_mask, shift=shift)

        return optsys


class NIRSpec(JWInstrument):
    """ A class modeling the optics of NIRSpec"""
    def __init__(self):
        JWInstrument.__init__(self, "NIRSpec")
        self.pixelscale = 0.0317 # for NIRCAM short-wavelen channels

    def _validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('NIRSpec does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
    def addCoronagraphOptics(self,optsys):
        raise NotImplementedError("No Coronagraph in NIRSpec!")


class TFI(JWInstrument):
    """ A class modeling the optics of the Tunable Filter Imager """
    def __init__(self):
        JWInstrument.__init__(self, "TFI")
        self.pixelscale = 0.064 # for TFI

        self.image_mask_list = ['CORON058', 'CORON075','CORON150','CORON200']
        self.pupil_mask_list = ['MASKC21N','MASKC66N','MASKC71N','MASK_NRM','CLEAR']

    def _validate_config(self):
        pass

    def addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics for TFI
        """
        if self.image_mask is 'CORON058':
            optsys.addImage(function='CircularOcculter', radius=0.58/2, name=self.image_mask)
        if self.image_mask is 'CORON075':
            optsys.addImage(function='CircularOcculter', radius=0.75/2, name=self.image_mask)
        if self.image_mask is 'CORON150':
            optsys.addImage(function='CircularOcculter', radius=1.5/2, name=self.image_mask)
        if self.image_mask is 'CORON200':
            optsys.addImage(function='CircularOcculter', radius=2.0/2, name=self.image_mask)

        # add pupil plane mask
        if ('pupil_shift_x' in self.options.keys() and self.options['pupil_shift_x'] != 0) or \
           ('pupil_shift_y' in self.options.keys() and self.options['pupil_shift_y'] != 0):
            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        else: shift = None

        if self.pupil_mask == 'MASKC21N':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC21N.fits", name=self.pupil_mask, shift=shift)
        if self.pupil_mask == 'MASKC66N':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC66N.fits", name=self.pupil_mask, shift=shift)
        if self.pupil_mask == 'MASKC71N':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKC71N.fits", name=self.pupil_mask, shift=shift)
        if self.pupil_mask == 'CLEAR':
            optsys.addPupil(transmission=self._datapath+"/coronagraph/MASKCLEAR.fits", name=self.pupil_mask, shift=shift)


        return optsys


class FGS(JWInstrument):
    """ A class modeling the optics of the FGS """
    def __init__(self):
        JWInstrument.__init__(self, "FGS")
        self.pixelscale = 0.069 # for FGS

    def _validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('FGS does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
        #TODO only one possible filter fot the FGS, too. 
    def addCoronagraphOptics(self,optsys):
        raise NotImplementedError("No Coronagraph in FGS!")


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




#########################3

class kurucz_stars(object):
    "A simple access object for a library of stellar spectra from the Kurucz models"
    def __init__(self):
        # The keys are spectral types; the values are tuples of the
        # table name and column name.
        self.Kurucz_filenames = {
            "O3V":   ("kp00_50000.fits", "g50"),
            "O5V":   ("kp00_45000.fits", "g50"),
            "O6V":   ("kp00_40000.fits", "g45"),
            "O8V":   ("kp00_35000.fits", "g40"),
            "O5I":   ("kp00_40000.fits", "g45"),
            "O6I":   ("kp00_40000.fits", "g45"),
            "O8I":   ("kp00_34000.fits", "g40"),
            "B0V":   ("kp00_30000.fits", "g40"),
            "B3V":   ("kp00_19000.fits", "g40"),
            "B5V":   ("kp00_15000.fits", "g40"),
            "B8V":   ("kp00_12000.fits", "g40"),
            "B0III": ("kp00_29000.fits", "g35"),
            "B5III": ("kp00_15000.fits", "g35"),
            "B0I":   ("kp00_26000.fits", "g30"),
            "B5I":   ("kp00_14000.fits", "g25"),
            "A0V":   ("kp00_9500.fits", "g40"),
            "A5V":   ("kp00_8250.fits", "g45"),
            "A0I":   ("kp00_9750.fits", "g20"),
            "A5I":   ("kp00_8500.fits", "g20"),
            "F0V":   ("kp00_7250.fits", "g45"),
            "F5V":   ("kp00_6500.fits", "g45"),
            "F0I":   ("kp00_7750.fits", "g20"),
            "F5I":   ("kp00_7000.fits", "g15"),
            "G0V":   ("kp00_6000.fits", "g45"),
            "G5V":   ("kp00_5750.fits", "g45"),
            "G0III": ("kp00_5750.fits", "g30"),
            "G5III": ("kp00_5250.fits", "g25"),
            "G0I":   ("kp00_5500.fits", "g15"),
            "G5I":   ("kp00_4750.fits", "g10"),
            "K0V":   ("kp00_5250.fits", "g45"),
            "K5V":   ("kp00_4250.fits", "g45"),
            "K0III": ("kp00_4750.fits", "g20"),
            "K5III": ("kp00_4000.fits", "g15"),
            "K0I":   ("kp00_4500.fits", "g10"),
            "K5I":   ("kp00_3750.fits", "g05"),
            "M0V":   ("kp00_3750.fits", "g45"),
            "M2V":   ("kp00_3500.fits", "g45"),
            "M5V":   ("kp00_3500.fits", "g50"),
            "M0III": ("kp00_3750.fits", "g15"),
            "M0I":   ("kp00_3750.fits", "g00"),
            "M2I":   ("kp00_3500.fits", "g00")}

        self.sptype_list = self.Kurucz_filenames.keys()

        def sort_sptype(typestr):
            letter = typestr[0]
            lettervals = {'O':0, 'B': 10, 'A': 20,'F': 30, 'G':40, 'K': 50, 'M':60}
            value = lettervals[letter]*1.0
            value += int(typestr[1])
            if "III" in typestr: value += .3
            elif "I" in typestr: value += .1
            elif "V" in typestr: value += .5
            return value
        
        self.sptype_list.sort(key=sort_sptype)

        self._JWPSF_basepath = os.path.dirname(os.path.dirname(os.path.abspath(poppy.__file__))) +os.sep+"data"




    def wavelengthUnits (self, hdu, column):
        """Interpret the units string in the table header.

        The function value will be the multiplicative factor needed
        to convert the wavelengths to microns.  If no units are
        specified (or can't be interpreted) for the wavelength column,
        the units will be assumed to be Angstroms.
        """

        ANGSTROMStoMICRONS = 0.0001
        NANOMETERStoMICRONS = 0.001
        METERStoMICRONS = 1.e6

        coldefs = hdu.get_coldefs()
        if isinstance (column, types.IntType):
            column_units = coldefs.units[column]
        else:
            column = column.lower()
            found = False
            for i in range (len (coldefs.names)):
                column_name = coldefs.names[i].lower()
                if column_name == column:
                    column_units = coldefs.units[i]
                    found = True
                    break
            if not found:
                print "warning:  can't find %s column" % column
                return ANGSTROMStoMICRONS

        if column_units is None:
            units = "angstrom"          # default
        else:
            units = column_units.lower()
        if units == "a" or units == "angstrom" or units == "angstroms":
            factor = ANGSTROMStoMICRONS
        elif units == "nm" or units == "nanometer" or units == "nanometers":
            factor = NANOMETERStoMICRONS
        elif units == "micron" or units == "microns":
            factor = 1.
        elif units == "m" or units == "meter" or units == "meters":
            factor = METERStoMICRONS
        else:
            print " wavelength units '%s' not given; " \
                  "Angstroms assumed" % column_units
            factor = ANGSTROMStoMICRONS

        return factor


    def specFromSpectralType (self, spectral_type):
        """Get spectrum specified by spectral type."""

        startype = spectral_type.upper()
        (fname, gcol) = self.Kurucz_filenames[startype]
        fullname = os.path.join (self._JWPSF_basepath, 'k93models', fname[0:4], fname)
        self.spectrum_file = fullname
        fd = pyfits.open (fullname)
        hdu = fd[1]
        data = hdu.data
        fd.close()
        factor = self.wavelengthUnits (hdu, "WAVELENGTH")
        wave = data.field('WAVELENGTH') * factor
        flux = data.field(gcol)
        self.spectrum = (wave, flux)
        self.spectral_type = startype

        spectrum = N.rec.fromarrays([wave,flux], names='wavelength_um, flux')

        return spectrum


#########################3

def test():
    jwst = JWST_OTE("path/to/some/OPDs")
    miri = jwst.MIRI
    gstar = pysynphot('G2 star')

    psf2 = miri.psf('imaging', center=(512,512), filter='F1500W', oversample=4, spectrum=gstar)

    corPSF = miri.psf('lyot', filter='F2550W', decenter=0.01, oversample=4)


def makeFakeFilter(filename, lcenter, dlam, clobber=False):
    """ arguments in microns, but file written in angstroms """

    lstart = lcenter - dlam/2
    lstop = lcenter + dlam/2

    nlambda = 40

    print "Filter from %f - %f " % (lstart, lstop)
    wavelength = N.linspace( lstart-dlam*0.1, lstop+dlam*0.1, nlambda)
    print wavelength
    transmission = N.zeros_like(wavelength)
    transmission[N.where( (wavelength > lstart) & (wavelength < lstop) )] = 1.0

    t = atpy.Table()
    t.add_column('WAVELENGTH', wavelength*1e4, unit='angstrom')
    t.add_column('THROUGHPUT', transmission)

    t.add_comment("This is a fake filter profile, represented as a top-hat function.")
    t.add_keyword("LAMBDA0",lcenter)
    t.add_keyword("DELTALAM",dlam)

    t.write(filename, overwrite=clobber)
    print("Created fake filter profile in "+filename)


    return t
def makeMIRIfilters():
    makeFakeFilter('F1065C_thru.fits',10.65, 0.53,clobber=True)
    makeFakeFilter('F1140C_thru.fits',11.40, 0.57,clobber=True)
    makeFakeFilter('F1550C_thru.fits',15.50, 0.78,clobber=True)
    makeFakeFilter('F2300C_thru.fits',23.00, 4.60,clobber=True)

    makeFakeFilter('FGS_thru.fits', 2.8, 4.40,clobber=True)
def makeNIRCamFilters():
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
#########################3


if __name__ == "__main__":

    if 0: 
        m = MIRI()
        m.filter = 'F1000W'
        m.calcPSF('test1.fits', clobber=True)
    #nc = NIRCam()
    #nc.filter = 'F200W'
    #nc.calcPSF('test_nircam.fits', mono=False)

    miri=MIRI()
    miri.filter='F1065C'
    miri.image_mask = 'FQPM1065'
    miri.pupil_mask = 'MASKFQPM'
    #nircam=NIRCam()
    tfi = TFI()
    nirspec = NIRSpec()

    if 0:
        miri.options = {'pupil_shift_y': 0.2, 'pupil_shift_x': 0}
        miri.calcPSF(nlambda=1)


    tfi.image_mask = 'CORON075'
    tfi.pupil_mask = 'MASKC21N'
    tfi.calcPSF(nlambda=1)
