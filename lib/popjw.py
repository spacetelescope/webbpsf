#!/usr/bin/env python
import numpy as N
import poppy
from poppy import *
import os
import glob
import atpy



class JWInstrument(object):
    def __init__(self, name=None):
        self.name=name

        self.JWPSF_basepath = os.path.dirname(os.path.dirname(os.path.abspath(poppy.__file__))) +os.sep+"data"

        self.datapath = self.JWPSF_basepath + os.sep + self.name + os.sep
        self.pupil = os.path.abspath(self.datapath+"../pupil.fits")
        self.pupilopd = None

        #create private instance variables. These will be
        # wrapped just below to create properties with validation.
        self.__filter=None
        self.__filter_files= [os.path.abspath(f) for f in glob.glob(self.datapath+os.sep+'filters/*.fits')]
        self.filter_list=[os.path.basename(f).split("_")[0] for f in self.__filter_files]
        self.filter = self.filter_list[0]

        self.__image_mask=None
        self.image_mask_list=[]

        self.__pupil_mask=None
        self.pupil_mask_list=[]

        self.pixelscale = 0.0

    def __validate_config(self):
        pass

    # create properties with error checking
    @property
    def filter(self):
        'Currently selected filter'
        return self.__filter
    @filter.setter
    def filter(self, value):
        if value not in self.filter_list:
            raise ValueError("Instrument %s doesn't have a filter called %s." % (self.name, value))
        self.__filter = value
        #self.__validate_config()
    @property
    def image_mask(self):
        'Currently selected image_mask'
        return self.__image_mask
    @image_mask.setter
    def image_mask(self, name):
        if name is not None:
            if name not in self.image_mask_list:
                raise ValueError("Instrument %s doesn't have an image mask called %s." % (self.name, value))
        self.__image_mask = name
    @property
    def pupil_mask(self):
        'Currently selected pupil_mask'
        return self.__pupil_mask
    @pupil_mask.setter
    def pupil_mask(self,name):
        if name is not None:
            if name not in self.pupil_mask_list:
                raise ValueError("Instrument %s doesn't have an pupil mask called %s." % (self.name, value))

        self.__pupil_mask = name
            #self.__validate_config()

    def __str__(self):
        return "JWInstrument name="+self.name

    #----- actual optical calculations follow here -----
    def calcPSF(self, filter=None, outfile=None,oversample=2, fov_arcsec=5., clobber=True, mono=False, nlambda=5 ):
        """ Compute a PSF
        """
        
        if filter is not None:
            self.filter = filter

        if outfile is None: 
            outfile = "PSF_%s_%s.fits" % (self.name, self.filter)
            #raise ValueError("You must specify an output file name.")

        self.optsys = self.getOpticalSystem(fov_arcsec=fov_arcsec, oversample=oversample)

        # for now, just calc monochromatic, centered on the filter transmission curve

        wf = N.where(self.filter_list == self.filter)
        wf = N.where(self.filter == N.asarray(self.filter_list))[0]
        filterdata = atpy.Table(self.__filter_files[wf])

        if mono:
            centerwave = (filterdata.THROUGHPUT*filterdata.WAVELENGTH).sum() / filterdata.THROUGHPUT.sum() / 1e10  # convert from angstroms to meters
            result = self.optsys.propagate(centerwave, display_intermediates=True, save_intermediates=True)
        else:
            print "CAUTION: Really basic top-hat function for filter profile, with %d steps" % nlambda
            wtrans = N.where(filterdata.THROUGHPUT > 0.5)
            lrange = filterdata.WAVELENGTH[wtrans] *1e-10
            lambd = N.linspace(N.min(lrange), N.max(lrange), nlambda)
            weights = N.ones(nlambda)
            source = {'wavelengths': lambd, 'weights': weights}
            result = self.optsys.calcPSF(source, display_intermediates=True, save_intermediates=True)

            f = p.gcf()
            p.text( 0.1, 0.95, "%s, filter= %s" % (self.name, self.filter), transform=f.transFigure, size='xx-large')

        result.writeto(outfile, clobber=clobber)


        # load filter profile
        # load stellar/target profile
        # compute output


    def getOpticalSystem(self,oversample=2, fov_arcsec=2):
        self.__validate_config()

        optsys = OpticalSystem(name='JWST+'+self.name, oversample=oversample)
        optsys.addPupil(name='JWST Pupil', transmission= self.pupil, opd=self.pupilopd)

        if self.image_mask is not None:
            optsys = self.addCoronagraphOptics(optsys)

        optsys.addDetector(self.pixelscale, fov_npix = fov_arcsec/self.pixelscale, name=self.name+" detector")
        return optsys


    def addCoronagraphOptics(self,optsys):
        """Add coronagraphic optics to an optical system. 
        This should be replaced by derived instrument classes. 
        """
        raise NotImplementedError("needs to be subclassed.")



class MIRI(JWInstrument):
    def __init__(self):
        JWInstrument.__init__(self, "MIRI")
        self.pixelscale = 0.11

        self.pupilopd = self.datapath+"OPD/MIRI_OPDisim1.fits"
        self.pupilopd = None
        self.image_mask_list = ['FQPM1065', 'FQPM1140', 'FQPM1550', 'LYOT2300']
        self.pupil_mask_list = ['MASKFQPM', 'MASKLYOT']

    def __validate_config(self):
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

        if self.image_mask == 'FQPM1065':
            optsys.addImage(function='FQPM',wavelength=10.65e-6, name=self.image_mask)
            optsys.addImage(function='fieldstop',size=24)
        elif self.image_mask == 'FQPM1140':
            optsys.addImage(function='FQPM',wavelength=11.40e-6, name=self.image_mask)
            optsys.addImage(function='fieldstop',size=24)
        elif self.image_mask == 'FQPM1550':
            optsys.addImage(function='FQPM',wavelength=15.50e-6, name=self.image_mask)
            optsys.addImage(function='fieldstop',size=24)
        elif self.image_mask == 'FQPM1550':
            optsys.addImage(function='CircularOcculter',radius =1., name=self.image_mask) 
            optsys.addImage(function='fieldstop',size=30)

        if self.pupil_mask == 'MASKFQPM':
            optsys.addPupil(transmission=self.datapath+"/coronagraph/MIRI_FQPMLyotStop.fits", name=self.pupil_mask)
        elif self.pupil_mask == 'MASKLYOT':
            optsys.addPupil(transmission=self.datapath+"/coronagraph/MIRI_LyotLyotStop.fits", name=self.pupil_mask)


        return optsys



class NIRCam(JWInstrument):
    def __init__(self):
        JWInstrument.__init__(self, "NIRCam")
        self.pixelscale = 0.0317 # for short-wavelen channels

        self.image_mask_list = ['BLC2100','BLC3350','BLC4300','WEDGESW','WEDGELW']
        self.pupil_mask_list = ['CIRCLYOT','WEDGELYOT']

    def __validate_config(self):
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

        if self.pupil_mask == 'MASKFQPM':
            optsys.addPupil(self.datapath+"/coronagraph/MIRI_FQPMLyotStop.fits", name=self.pupil_mask)
        elif self.pupil_mask == 'MASKLYOT':
            optsys.addPupil(self.datapath+"/coronagraph/MIRI_LyotLyotStop.fits", name=self.pupil_mask)

        return optsys


class NIRSpec(JWInstrument):
    def __init__(self):
        JWInstrument.__init__(self, "NIRSpec")
        self.pixelscale = 0.0317 # for NIRCAM short-wavelen channels

    def __validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('NIRSpec does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
    def addCoronagraphOptics(self,optsys):
        raise NotImplementedError("No Coronagraph in NIRSpec!")


class TFI(JWInstrument):
    def __init__(self):
        JWInstrument.__init__(self, "TFI")
        self.pixelscale = 0.064 # for TFI

        self.image_mask_list = ['CORON058', 'CORON075','CORON150','CORON200']
        self.pupil_mask_list = ['MASKC21N','MASKC66N','MASKC71N']

    def __validate_config(self):
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

        if self.pupil_mask == 'MASKC21N':
            optsys.addPupil(self.datapath+"/coronagraph/MIRI_FQPMLyotStop.fits", name=self.pupil_mask)
        if self.pupil_mask == 'MASKC66N':
            optsys.addPupil(self.datapath+"/coronagraph/MIRI_FQPMLyotStop.fits", name=self.pupil_mask)
        if self.pupil_mask == 'MASKC71N':
            optsys.addPupil(self.datapath+"/coronagraph/MIRI_FQPMLyotStop.fits", name=self.pupil_mask)


        return optsys


class FGS(JWInstrument):
    def __init__(self):
        JWInstrument.__init__(self, "FGS")
        self.pixelscale = 0.069 # for FGS

    def __validate_config(self):
        if (not self.image_mask is None) or (not self.pupil_mask is None):
            raise ValueError('FGS does not have image or pupil masks!')
            self.image_mask = None
            self.pupil_mask = None
        #TODO only one possible filter fot the FGS, too. 
    def addCoronagraphOptics(self,optsys):
        raise NotImplementedError("No Coronagraph in FGS!")




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

    print "Filter from %f - %f " % (lstart, lstop)
    wavelength = N.linspace( lstart-dlam*0.1, lstop+dlam*0.1, 20)
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


if __name__ == "__main__":

    p.clf()
    if 0: 
        m = MIRI()
        m.filter = 'F1000W'
        m.calcPSF('test1.fits', clobber=True)
    nc = NIRCam()
    nc.filter = 'F200W'
    #nc.calcPSF('test_nircam.fits', mono=False)

    miri=MIRI()
    miri.filter='F1065C'
    miri.image_mask = 'FQPM1065'
    miri.pupil_mask = 'MASKFQPM'
    nircam=NIRCam()
    tfi = TFI()
    nirspec = NIRSpec()



