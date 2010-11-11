#!/usr/bin/env  python

import utils
import numarray as NIX
from NFFT2 import  *
from SimpleFits import *
from MakeDisk import *
# from MakeJinc2 import *

INIT = 0x01
PROP = 0x02
AMPL = 0x03
PHASE = 0x04 # not used yet....


class OpticalTrain:

    def __init__(self, name=None,  pupil=None, pupilFile=None, oversample=None,  \
                       phase=None, phaseFile=None, wavelength=None,\
                       fftwshape=None, FFTWstuff=None,  
                       plantype=FFTW_ESTIMATE, calctype=Float64, verbose=None, \
                       savetype=Float64):
        """ 
        Stay w/even square arrays
        anand@stsci.edu  September 2003 
        """

        self.name = name
        self.calctype = calctype
        self.oversample = oversample
        self.mycomplex = type(1.0 + 1j * 1.0)
        self.t = savetype
    
        if pupilFile:
            self.pupil = SimpleFitsRead(fn=pupilFile,  \
                                report=0xffff,  verbose=verbose).data
            self.pupilFile = pupilFile
        else:
            self.pupil = pupil

        if phaseFile:
            self.phase = SimpleFitsRead(fn=phaseFile,  \
                                report=0xffff,  verbose=verbose).data
            self.phaseFile = phaseFile
        else:
            self.phase = phase

        # set up oversampling, and fftw shape, create FFTW plan
        # To use existing FFTWstuff be sure array sizes match!!!!
        # caveat usor...
        if FFTWstuff is None:
            if oversample is None:
                oversample = 4
            self.oversample = int(oversample)
            self.fftwshape = (self.pupil.shape[0] * oversample, \
                              self.pupil.shape[1] * oversample)
            self.FFTWstuff = NFFT2(self.fftwshape, plantype=plantype)
        else:
            self.fftwshape = fftwshape
            self.FFTWstuff = FFTWstuff
            print """
           WARNING FROM %s initialization:
           I will use your user-supplied FFTWstuff. 
           It had better be correct because I am not checking for array size
           matches etc etc due to different input pupils/phases/whatnot...

            Remember to give me a non-None fftwshape too, friend....

           """  % self.name


        (R, C) = self.fftwshape
        (r, c) = self.pupil.shape
        if r%2 == 1:
            evenr = 1
            print "caution: odd-sized r: offsetting pupil in oversized array by one pixel"
        else:
            evenr = 0
        if c%2 == 1:
            evenc = 1
            print "caution: odd-sized c: offsetting pupil in oversized array by one pixel"
        else:
            evenc = 0

        # default is to save the whole image plane
        self.savebox = R
        # default is to save the active pupil plane
        self.psavebox = r

        # create large pupil array and embed pupil in middle
        fftwpupil = NIX.zeros(self.fftwshape, calctype)
        #print "R/2 - r/2 : R/2 + r/2, C/2 - c/2 : C/2 + c/2"
        #print "%9d:%9d, %9d:%9d" % (R/2 - r/2, R/2 + r/2, C/2 - c/2, C/2 + c/2)
        #print "    %9d         %9d" % (R/2 - r/2 - ( R/2 + r/2), C/2 - c/2 - (C/2 + c/2))
        fftwpupil[R/2 - r/2 : R/2 + r/2 + evenr, C/2 - c/2 : C/2 + c/2 + evenc] \
            = self.pupil.copy()

        # create large phase array and embed phase in middle
        fftwphase = zeros(self.fftwshape, calctype)
        fftwphase[R/2 - r/2 : R/2 + r/2 + evenr, C/2 - c/2 : C/2 + c/2 + evenc] \
            = self.phase.copy()
        
        # create large first plane field strength
        E = fftwpupil * NIX.exp(1j * fftwphase)

        # create list of actions in train...
        self.planes = [E,]
        self.oplist = ['entrance pupil',]
        self.splist = [PUPIL,]
        self.planedict = {'entrance pupil':0}


    def amplitudeMask(self, mask=None, operation=None):
        """ 
        take last field E to new field F with correct flopness
        using the (physical space flopness) mask array
        and append to optical train's planes
        """
        if mask.shape != self.fftwshape:
            print "Optical Train FATAL ERROR --- "
            print "input array in %s wrong shape: "  %  operation,
            print self.fftwshape,
            print " != ",
            print mask.shape
        # The field strength as we find it...
        E = self.planes[-1]

        # decide if we're working in pupil or image space
        if self.splist[-1] == SPACE:
            spacetype = PUPIL
            fmask = mask
        else:
            fmask = utils.flop(mask)
            spacetype = IMAGE
            
        # store the mask as per [un]flopped state as appropriate
        self.planes.append(fmask)
        self.oplist.append(operation)
        self.splist.append(spacetype)
        self.planedict[operation] =  len(self.splist) - 1

        # multiply by the mask with correct flopness
        # if fmask.type() is Complex64:
        #   print "The Lyot Stop is of type %s" % fmask.type()
        #   F = E * sqrt(fmask.real*fmask.real + fmask.imag*fmask.imag)
        #   print "E is of type %s" % E.type()
        #   print "F is of type %s" % F.type()
        #   SimpleFitsWrite(fn='LyotStop.fits', data=F, clobber='y')
        # else:
        F = E * fmask
        # now save the field strength and appropriate bookkeeping data
        self.planes.append(F)
        self.oplist.append('after ' + operation)
        self.splist.append(spacetype)
        self.planedict['after ' + operation] =  len(self.splist) - 1


    def toConjugate(self, operation=None):
        """ 
        take last field E to new field F with correct transform direction
        and append to optical train's planes
        """

        if self.splist[-1] == PUPIL:
            direction = FFTW_FORWARD
        else:
            direction = FFTW_BACKWARD

        F = self.FFTWstuff.transform(self.planes[-1], direction)
        self.planes.append(F)
        self.oplist.append(operation)

        if direction == FFTW_FORWARD:
            self.splist.append(IMAGE)
        else:
            self.splist.append(PUPIL)
        self.planedict[operation] =  len(self.splist) - 1



    def setsavetype(self, savetype=None):
        print "%s saving files as " % self.name, self.t
        self.t = savetype

    def setsavebox_pixels(self, saveboxpixels=None):
        # sets size of image plane (diam or fullsize) to be saved to files
        self.savebox = int(saveboxpixels)


    def setpsavebox(self, psavebox=None):
        # sets size of pupil plane (diam or full size) to be saved to files
        self.psavebox = int(psavebox)


    def displayproperties(self):
        print ""
        print """Optical train "%s" properties are:""" % self.name
        print "\tPupil shape = ",
        print self.pupil.shape
        print "\tOversample = %d" % self.oversample
        print "\tFFTW shape  = ",
        print self.fftwshape

        if hasattr(self, 'pupilFile'):
            print "\tPupil file was = %s" % self.pupilFile
        else:
            print "\tPupil was initialized w/ array"
        print "\tIncoming power in pupil = %.4g" % ((self.pupil*self.pupil).sum())

        if hasattr(self, 'phaseFile'):
            print "\tPhase file was = %s" % self.phaseFile
        else:
            print "\tPhase was initialized w/ array"

        if hasattr(self, 'FFTWstuff'):
            print "\tFFTWstuff is set up"

        if hasattr(self, 'savebox'):
            print "\tSaving +/- %d reselts image box" % self.savebox

        print ""


    def saveIntensity(self, planename=None, fn=None, clobber='n', kwdict=None, \
                      kwlist=None, multfactor=None, logThresh=None):

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            E = self.planes[planenum]
            if self.splist[planenum] == PUPIL:
                tmp = utils.centralSection(E, npoints=self.psavebox)
            else:
                tmp = utils.centralSection(utils.flop(E), npoints=self.savebox)

            pow = tmp.real * tmp.real  + tmp.imag * tmp.imag
            if multfactor:
                pow = pow*multfactor
            if logThresh:
                pow = NIX.log10(pow + logThresh)
            SimpleFitsWrite(fn=fn, data=pow.astype(self.t), clobber=clobber, kwdict=kwdict, \
                            kwlist=kwlist)
        else:
            print "plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)


    def saveComplexAmplitude(self, planename=None, fn=None, clobber='n', \
                             kwdict=None, kwlist=None):

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            E = self.planes[planenum]
            if self.splist[planenum] == PUPIL:
                tmp = utils.centralSection(E, npoints=self.psavebox)
            else:
                tmp = utils.centralSection(utils.flop(E), npoints=self.savebox)

            camp = NIX.zeros((2, tmp.shape[0], tmp.shape[1]), type=Float64)
            camp[0,:,:] = tmp.real.copy()
            camp[1,:,:] = tmp.imag.copy()
            SimpleFitsWrite(fn=fn, data=camp.astype(self.t), clobber=clobber, kwdict=kwdict, \
                            kwlist=kwlist)
        else:
            print "plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)

    def saveComponent(self, planename=None, fn=None, clobber='n', kwdict=None, \
                      kwlist=None, component="real"):

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            if component == "real":
                E = self.planes[planenum].real.copy()
            if component == "imag":
                E = self.planes[planenum].imag.copy()
            if self.splist[planenum] == PUPIL:
                tmp = utils.centralSection(E, npoints=self.psavebox)
            else:
                tmp = utils.centralSection(utils.flop(E), npoints=self.savebox)

            SimpleFitsWrite(fn=fn, data=tmp.astype(self.t), clobber=clobber, kwdict=kwdict, \
                            kwlist=kwlist)
        else:
            print "OT: saveComponent: plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)


    def getComplexAmplitude(self, planename=None):
        # returns the complex amplitude array as is...

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            E = self.planes[planenum].copy()
            return E
        else:
            print "plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)
            return None


    def getIntensity(self, planename=None):
        # returns the intensity of a complex amplitude  array as is...

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            E = self.planes[planenum].copy()
            I = E.real*E.real + E.imag*E.imag
            if self.splist[planenum] == PUPIL:
                tmp = utils.centralSection(I, npoints=self.psavebox)
            else:
                tmp = utils.centralSection(utils.flop(I), npoints=self.savebox)
            return tmp
        else:
            print "plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)
            return None


    def planeinfo(self, planename=None):
        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            s1 = "Power in %s  %.3g" % \
                    (planename, utils.power2d(self.planes[planenum]))
            a = self.planes[planenum]
            if type(a[0,0]) == self.mycomplex:
                intensity = a.real*a.real + a.imag*a.imag
            else:
                intensity = a * a
            s2 = "Maximum intensity %.3g" %  intensity.max()
        else:
            s1 =  "OT.totalPower: plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)
            s2 = " "
        rstr = """
   %s
   %s"""  %  (s1, s2)
        return rstr


    def totalPower(self, planename=None):

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            return utils.power2d(self.planes[planenum])
        else:
            print "OT.totalPower: plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)
            return None


    def maxIntensity(self, planename=None):

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            a = self.planes[planenum]
            intensity = a.real*a.real + a.imag*a.imag
            return intensity.max()
        else:
            print "OT.maxIntensity: plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)
            return None


    def totalArea(self, planename=None):

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            return add.reduce(abs(self.planes[planenum].flat))
        else:
            print "OT.totalArea: plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)
            return None




    def saveTransmissionMask(self, planename=None, fn=None, clobber='n', \
                             kwdict=None, kwlist=None):
    # flop the in-mem mask if it is imagesque

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            M = self.planes[planenum]
            if self.splist[planenum] == PUPIL:
                tmp = utils.centralSection(M,  npoints=self.psavebox)
            else:
                tmp = utils.centralSection(utils.flop(M), npoints=self.savebox)

            if tmp.type() is Complex64:
                print "Writing real and complex parts of transmission Mask"
                SimpleFitsWrite(fn=fn[:-4]+"real.fits", data=tmp.real.astype(self.t), \
                    clobber=clobber, kwdict=kwdict, kwlist=kwlist)
                SimpleFitsWrite(fn=fn[:-4]+"imag.fits", data=tmp.imag.astype(self.t), \
                    clobber=clobber, kwdict=kwdict, kwlist=kwlist)
            else:
                SimpleFitsWrite(fn=fn, data=tmp.astype(self.t), clobber=clobber, kwdict=kwdict, kwlist=kwlist)
        else:
            print "plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)


    def deletePlane(self, planename=None):

        if self.planedict.has_key(planename):
            planenum = self.planedict[planename]
            self.planes[planenum] = None
        else:
            print "plane '%s' not found in optical train '%s'"  %  \
                   (planename, self.name)
            print "Cannot delete plane '%s' "


    def displaytrain(self):

        print """\ntrain "%s" is currently...""" %  self.name
        #for op in self.oplist:
        #   print "\t%s" % op
        for i in range(len(self.oplist)):
            if self.splist[i] == PUPIL:
                print "\tPupil: ",
            if self.splist[i] == IMAGE:
                print "\tImage: ",
            print "%s" % self.oplist[i]


if __name__ == '__main__':

    print """

    Please use the Python module 'trainDriverMW.py' to run this code!

    """

