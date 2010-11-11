#!/usr/bin/env python
import numpy as N
from poppy import *



class JWInstrument():
    def __init__(self, name=name):
        self.name=name

        #create private instance variables. These will be
        # wrapped just below to create properties with validation.
        self.__filter=None
        self.__filter_list=[]

        self.__image_mask=None
        self.__image_mask_list=[]

        self.__pupil_mask=None
        self.__pupil_mask_list=[]

        self.pixelscale = 0.0

    def __validate_config(self):
        pass

    # create 'filter' property with error checking
    @property
    def filter(self):
        doc='Currently selected filter'
        def fget(self):
            return self.__filter
        def fset(self, name):
            assert name in self.__filter_list
            self.__filter = name
            self.__validate_config()
        def fdel(self):
            del self.__filter
    # create 'image_mask' property with error checking
    @property
    def image_mask(self):
        doc='Currently selected image_mask'
        def fget(self):
            return self.__image_mask
        def fset(self, name):
            if name is not None:
                assert name in self.__image_mask_list
            self.__image_mask = name
            self.__validate_config()
        def fdel(self):
            del self.__image_mask
    # create 'pupil_mask' property with error checking
    @property
    def pupil_mask(self):
        doc='Currently selected pupil_mask'
        def fget(self):
            return self.__pupil_mask
        def fset(self, name):
            if name is not None:
                assert name in self.__pupil_mask_list
            self.__pupil_mask = name
            self.__validate_config()
        def fdel(self):
            del self.__pupil_mask
 
    def __str__(self):
        return "JWInstrument name="+self.name


    #----- actual optical calculations follow here -----
    def calcPSF(self, filter=None,oversample=2, fov_arcsec=5.):

        self.optsys = self.__getOpticalSystem()

        # load filter profile
        # load stellar/target profile
        # compute output


    def __getOpticalSystem(self,oversample=2):
        if filter is not None:
            self.filter=filter

        optsys = OpticalSystem('JWST+'+self.name, oversample=oversample)
        optsys.addPupil('JWST Pupil', amplitude = self.pupil, opd=self.pupilopd)

        if self.image_mask is not None:
            optsys = self.__addCoronagraphOptics(optsys)

        optsys.addDetector(self.pixelscale, fov_npix = fov_arcsec/self.pixelscale)


    def __addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics to an optical system. 
        This should be replaced by derived instrument classes. 
        """
        if self.image_mask is not None:
            optsys.addImage(self.image_mask)
        if self.pupil_mask is not None:
            optsys.addPupil(self.pupil_mask)
        return optsys


    



class MIRI(JWInstrument):

    """ An instrument observation for simulating a PSF can be specified by
     - aperture
     - pixel location
     - filter

     In some cases, specifying an aperture is sufficient to imply the pixel location and filter too. 

     """
    def __init__(self):

        self.JWInstrument.__init__(self, "MIRI")
        self.pixelscale = 0.11

        self.datapath = "./data/MIRI/"
        self.pupil = "./data/pupil.fits"
        self.pupilopd = self.datapath+"OPD/MIRI_OPDisim1.fits"

    def __validate_config(self):
        pass
        if self.filter == 'F1065C':
            self.image_mask = 'FQPM'
            self.pupil_mask = 'MIRI_FQPMLyotStop.fits'
        elif self.filter == 'F1550C':
            self.image_mask = 'FQPM'
            self.pupil_mask = 'MIRI_FQPMLyotStop.fits'
        else:
            self.image_mask = None
            self.pupil_mask = None

        # TODO if coronagraphic mode, check the filter+mask combo is legit.

    def __addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics for MIRI 
        """
        if self.image_mask is 'FQPM':
            optsys.addImage(self.image_mask)
            optsys.addPupil(self.pupil_mask)
        elif self.image_mask is 'Lyot':
            optsys.addImage(self.image_mask)
            optsys.addPupil(self.pupil_mask)

        return optsys



class NIRCam(JWInstrument):
    def __init__(self):

        self.JWInstrument.__init__(self, "NIRCam")

        self.pixelscale = 0.0317 # for short-wavelen channels

    def __validate_config(self):
        pass
        if self.filter == 'F1065C':
            self.image_mask = 'FQPM'
            self.pupil_mask = 'MIRI_FQPMLyotStop.fits'
        elif self.filter == 'F1550C':
            self.image_mask = 'FQPM'
            self.pupil_mask = 'MIRI_FQPMLyotStop.fits'
        else:
            self.image_mask = None
            self.pupil_mask = None


    def __addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics for MIRI 
        """
        if self.image_mask is 'FQPM':
            optsys.addImage(self.image_mask)
            optsys.addPupil(self.pupil_mask)
        elif self.image_mask is 'Lyot':
            optsys.addImage(self.image_mask)
            optsys.addPupil(self.pupil_mask)

        return optsys


class NIRSpec(JWInstrument):
    def __init__(self):

        self.JWInstrument.__init__(self, "NIRSpec")

        self.pixelscale = 0.0317 # for NIRCAM short-wavelen channels

    def __validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('NIRSpec does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
    def __addCoronagraphOptics(self,optsys):
        return optsys # do nothing!


class TFI(JWInstrument):
    def __init__(self):

        self.JWInstrument.__init__(self, "TFI")

        self.pixelscale = 0.064 # for TFI

    def __validate_config(self):
        pass


    def __addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics for MIRI 
        """
        if self.image_mask is 'FQPM':
            optsys.addImage(self.image_mask)
            optsys.addPupil(self.pupil_mask)
        elif self.image_mask is 'Lyot':
            optsys.addImage(self.image_mask)
            optsys.addPupil(self.pupil_mask)

        return optsys


class FGS(JWInstrument):
    def __init__(self):

        self.JWInstrument.__init__(self, "TFI")
        self.pixelscale = 0.069 # for FGS

    def __validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('FGS does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
        #TODO only one possible filter fot the FGS, too. 
    def __addCoronagraphOptics(self,optsys):
        return optsys # do nothing!




def test():
    jwst = JWST_OTE("path/to/some/OPDs")
    miri = jwst.MIRI
    gstar = pysynphot('G2 star')

    psf2 = miri.psf('imaging', center=(512,512), filter='F1500W', oversample=4, spectrum=gstar)

    corPSF = miri.psf('lyot', filter='F2550W', decenter=0.01, oversample=4)


if __name__ == "__main__":
        



