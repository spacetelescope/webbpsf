#! /usr/bin/env  python 
"""
    "Slow Fourier Transform"

    Matrix-based Fourier transforms for computing PSFs. 
    See Soummer et al. 2007 JOSA


    Example
    -------
    sf = SFT.SlowFourierTransform()
    result = sf.perform(pupilArray, focalplane_size, focalplane_npix)




    History
    -------
    Code originally by A. Sivaramakrishnan
    2010-11-05 Revised normalizations for flux conservation consistent
        with Soummer et al. 2007. Updated documentation.  -- M. Perrin

"""


import numpy as np
import utils
import pyfits
#import SimpleFits as SF


# FFTSTYLE
def SFT1(pupil, nlamD, npix):
    """

    Compute an "FFTSTYLE" matrix fourier transform. 

    Parameters
    ----------
    pupil
        pupil array (n by n)
    nlamD
        size of focal plane array, in units of lam/D
        (corresponds to 'm' in Soummer et al. 2007 4.2)
    npix
        number of pixels per side side of focal plane array
        (corresponds to 'N_B' in Soummer et al. 2007 4.2)


    """

    npup = pupil.shape[0]

    du = nlamD / float(npix)
    dv = nlamD / float(npix)

    dx = 1.0/float(npup)
    dy = 1.0/float(npup)

    Xs = (np.arange(npup) - (npup/2)) * dx
    Ys = (np.arange(npup) - (npup/2)) * dy

    Us = (np.arange(npix) - npix/2) * du
    Vs = (np.arange(npix) - npix/2) * dv

    XU = np.outer(Xs, Us)
    YV = np.outer(Ys, Vs)


    expXU = np.exp(-2.0 * np.pi * 1j * XU)
    expYV = np.exp(-2.0 * np.pi * 1j * YV)

    expXU = expXU.T.copy()
    t1 = np.dot(expXU, pupil)

    t2 = np.dot(t1, expYV)

    #return  nlamD/(npup*npix) *   t2 * dx * dy
    return  float(nlamD)/(npup*npix) *   t2 


# SYMMETRIC
def SFT2(pupil, nlamD, npix):
    """
    Compute a "SYMMETRIC" matrix fourier transform. 

    Parameters
    ----------
    pupil
        pupil array (n by n)
    nlamD
        size of focal plane array, in units of lam/D
        (corresponds to 'm' in Soummer et al. 2007 4.2)
    npix
        number of pixels per side side of focal plane array
        (corresponds to 'N_B' in Soummer et al. 2007 4.2)


    """



    npup = pupil.shape[0]

    du = nlamD / float(npix)
    dv = nlamD / float(npix)

    dx = 1.0/float(npup)
    dy = 1.0/float(npup)

    Xs = (np.arange(npup) - float(npup)/2.0 + 0.5) * dx
    #print "xs ", xs.shape,
    #xs.shape = (1,npup)
    #xs.transpose()

    Ys = (np.arange(npup) - float(npup)/2.0 + 0.5) * dy
    #ys.shape = (1, npup)
    #ys.transpose()

    Us = (np.arange(npix) - float(npix)/2.0 + 0.5) * du
    #print "ps", ps.shape
    #ps.shape = (1,npix)
    Vs = (np.arange(npix) - float(npix)/2.0 + 0.5) * dv
    #qs.shape = (1,npix)

    XU = np.outer(Xs, Us)
    YV = np.outer(Ys, Vs)

    expXU = np.exp(-2.0 * np.pi * 1j * XU)
    expYV = np.exp(-2.0 * np.pi * 1j * YV)

    expXU = expXU.T.copy()
    t1 = np.dot(expXU, pupil)

    t2 = np.dot(t1, expYV)

    #return t2 * dx * dy
    return  float(nlamD)/(npup*npix) *   t2 



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



def test_SFT(choice='FFTSTYLE', outdir='.', outname='SFT1'):
    import os

    print "Testing SFT, style = "+choice

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
    npix = 1024
    u = 100    # of lam/D
    s = (npupil,npupil)


    # FFT style
    sft1 = SlowFourierTransform(choice=choice)

    ctr = (float(npupil)/2.0 + sft1.offset(), float(npupil)/2.0 + sft1.offset())
    #print ctr
    pupil = utils.makedisk(s=s, c=ctr, r=float(npupil)/2.0001, t=np.float64, grey=0)

    pupil /= np.sqrt(pupil.sum())

    pyfits.PrimaryHDU(pupil.astype(np.float32)).writeto(outdir+os.sep+outname+"pupil.fits", clobber=True)

    a = sft1.perform(pupil, u, npix)

    pre = (abs(pupil)**2).sum() 
    post = (abs(a)**2).sum() 
    ratio = post / pre
    calcr = 1./(u**2 *npix**2)     # multiply post by this to make them equal
    print "Pre-FFT  total: "+str( pre)
    print "Post-FFT total: "+str( post )
    print "Ratio:          "+str( ratio)
    #print "Calc ratio  :   "+str( calcr)
    #print "uncorrected:    "+str( ratio/calcr)


    complexinfo(a, str="sft1 asf")
    #print 
    asf = a.real.copy()
    #SF.SimpleFitsWrite(fn=outdir+os.sep+outname+"asf.fits", data=asf.astype(np.float32), clobber='y')
    pyfits.PrimaryHDU(asf.astype(np.float32)).writeto(outdir+os.sep+outname+"asf.fits", clobber=True)
    cpsf = a * a.conjugate()
    psf = cpsf.real.copy()
    #SF.SimpleFitsWrite(fn=outdir+os.sep+outname+"psf.fits", data=psf.astype(np.float32), clobber='y')
    pyfits.PrimaryHDU(psf.astype(np.float32)).writeto(outdir+os.sep+outname+"psf.fits", clobber=True)





if __name__ == "__main__":

    test_SFT('FFTSTYLE', outname='SFT1')
    test_SFT('SYMMETRIC', outname='SFT2')
