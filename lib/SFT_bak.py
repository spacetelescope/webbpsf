#! /usr/bin/env  python 
import numpy as np
import utils
import SimpleFits as SF


# FFTSTYLE
def SFT1(pupil, u, nimg):
    """

    Compute an "FFTSTYLE" matrix fourier transform. 

    Parameters
    ----------
    pupil
        pupil array (n by n)
    u
        size of focal plane array, in units of lam/D
    nimg
        number of pixels per side side of focal plane array
    """

    npup = pupil.shape[0]

    dp = u / float(nimg)
    dq = u / float(nimg)

    dx = 1.0/float(npup)
    dy = 1.0/float(npup)

    xs = (np.arange(npup) - (npup/2)) * dx
    ys = (np.arange(npup) - (npup/2)) * dy
    #print "xs ", xs.shape,

    ps = (np.arange(nimg) - nimg/2) * dp
    qs = (np.arange(nimg) - nimg/2) * dq

    #print "ps", ps.shape
    xp = np.outer(xs, ps)
    yq = np.outer(ys, qs)
    #print "xp = out(xs,ps)", xp.shape, "yq = out(ys,qs)", yq.shape


    expxp = np.exp(-2.0 * np.pi * 1j * xp)
    expyq = np.exp(-2.0 * np.pi * 1j * yq)

    expxp = expxp.T.copy()
    t1 = np.dot(expxp, pupil)

    t2 = np.dot(t1, expyq)

    return t2 * dx * dy


# SYMMETRIC
def SFT2(pupil, u, nimg):
    """
    pupil array (n by n)
    u foc plane size lam/D
    nimg num pix on a side in focalplane
    """

    npup = pupil.shape[0]

    dp = u / float(nimg)
    dq = u / float(nimg)

    dx = 1.0/float(npup)
    dy = 1.0/float(npup)

    xs = (np.arange(npup) - float(npup)/2.0 + 0.5) * dx
    #print "xs ", xs.shape,
    #xs.shape = (1,npup)
    #xs.transpose()

    ys = (np.arange(npup) - float(npup)/2.0 + 0.5) * dy
    #ys.shape = (1, npup)
    #ys.transpose()

    ps = (np.arange(nimg) - float(nimg)/2.0 + 0.5) * dp
    #print "ps", ps.shape
    #ps.shape = (1,nimg)
    qs = (np.arange(nimg) - float(nimg)/2.0 + 0.5) * dq
    #qs.shape = (1,nimg)

    xp = np.outer(xs, ps)
    yq = np.outer(ys, qs)

    expxp = np.exp(-2.0 * np.pi * 1j * xp)
    expyq = np.exp(-2.0 * np.pi * 1j * yq)

    expxp = expxp.T.copy()
    t1 = np.dot(expxp, pupil)

    t2 = np.dot(t1, expyq)

    return t2 * dx * dy



class SlowFourierTransform:
    """Implements a discrete matrix Fourier transform for optical 
    propagation, following the algorithms discussed in 
    Soummer et al. 2007 JOSA 15 24

    Parameters
    ----------
    choice
        Either 'SYMMETRIC' or 'FFTSTYLE'. Sets whether the DFT result is
        centered at pixel n/2+1 (FFTSTYLE) or on the crosshairs between
        the central pixels (SYMMETRIC). Default is FFTSTYLE. 


    Example
    -------
    sft = SlowFourierTransform()
    sft.perform(pupilArray, focalplane_size, focalplane_npix)


    History
    -------
    Code by Sivaramakrishnan based on Soummer et al.
    2010-01 Documentation updated by Perrin

    """

    def __init__(self, choice="FFTSTYLE", verbose=False):

        self.verbose=verbose

        self.choix = ("FFTSTYLE", "SYMMETRIC")
        self.correctoffset = {self.choix[0]: 0.5, self.choix[1]: 0.0}
        if choice not in self.choix:
            print "Error: choice must be one of ", selff.choix
        self.choice = choice

        if self.choice == self.choix[0]:
            if self.verbose:
                print "Announcement  - This instance of SlowFourierTransform uses SFT1"
                print "This is a FFTSTYLE sft set-up"
            self.perform = SFT1
        if self.choice == self.choix[1]:
            if self.verbose:
                print "Announcement  - This instance of SlowFourierTransform uses SFT2"
                print "This is a SYMMETRIC sft set-up"
            self.perform = SFT2

    def offset(self):
            return self.correctoffset[self.choice]


    def performFITS(hdulist, focalplane_size, focalplane_npix):
        """ Perform an MFT, and return the result as a pyfits.HDUlist """
        newHDUlist = hdulist.copy()
        newim = self.perform(hdulist[0].data, focalplane_size, focalplane_npix)

        newHDUlist[0].data = newim
        #TODO fits header keyword updates

        return newHDUlist



def test_SFT(choice='FFTSTYLE', outdir='testdata', outname='SFT1'):
    import os

    def complexinfo(a, str=None):

        if str:
            print 
            print "\t", str
        re = a.real.copy()
        im = a.imag.copy()
        print "\t%.2e  %.2g  =  re.sum im.sum" % (re.sum(), im.sum())
        print "\t%.2e  %.2g  =  abs(re).sum abs(im).sum" % (abs(re).sum(), abs(im).sum())


    npupil = 156
    pctr = int(npupil/2)
    nimg = 256
    nimg = 256
    u = 16.0
    s = (npupil,npupil)


    # FFT style
    sft1 = SlowFourierTransform(choice=choice)

    ctr = (float(npupil)/2.0 + sft1.offset(), float(npupil)/2.0 + sft1.offset())
    #print ctr
    pupil = utils.makedisk(s=s, c=ctr, r=float(npupil)/2.0001, t=np.float64, grey=0)

    pupil /= np.sqrt(pupil.sum())

    SF.SimpleFitsWrite(fn=outdir+os.sep+outname+"pupil.fits", data=pupil.astype(np.float32), clobber='y')

    a = sft1.perform(pupil, u, nimg)

    pre = (abs(pupil)**2).sum() 
    post = (abs(a)**2).sum() 
    ratio = post / pre
    calcr = 1./(u**2 *nimg**2)     # multiply post by this to make them equal
    #print "Pre-FFT  total: "+str( pre)
    #print "Post-FFT total: "+str( post )
    print "Ratio:          "+str( ratio)
    print "Calc ratio  :   "+str( calcr)
    print "uncorrected:    "+str( ratio/calcr)


    complexinfo(a, str="sft1 asf")
    #print 
    asf = a.real.copy()
    SF.SimpleFitsWrite(fn=outdir+os.sep+outname+"asf.fits", data=asf.astype(np.float32), clobber='y')
    cpsf = a * a.conjugate()
    psf = cpsf.real.copy()
    SF.SimpleFitsWrite(fn=outdir+os.sep+outname+"psf.fits", data=psf.astype(np.float32), clobber='y')





if __name__ == "__main__":

    test_SFT('FFTSTYLE', outname='SFT1')
    #test_SFT('SYMMETRIC', outname='SFT2')
