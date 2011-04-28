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
    Ys = (np.arange(npup) - float(npup)/2.0 + 0.5) * dy

    Us = (np.arange(npix) - float(npix)/2.0 + 0.5) * du
    Vs = (np.arange(npix) - float(npix)/2.0 + 0.5) * dv

    XU = np.outer(Xs, Us)
    YV = np.outer(Ys, Vs)


    expXU = np.exp(-2.0 * np.pi * 1j * XU)
    expYV = np.exp(-2.0 * np.pi * 1j * YV)
    expXU = expXU.T.copy()

    t1 = np.dot(expXU, pupil)
    t2 = np.dot(t1, expYV)

    #return t2 * dx * dy
    return  float(nlamD)/(npup*npix) *   t2 


# ADJUSTIBLE
def SFT3(pupil, nlamD, npix, offset=(0.0,0.0)):
    """
    Compute an adjustible-center matrix fourier transform. 

    For an output array with ODD size n,
    the PSF center will be at the center of pixel (n-1)/2
    
    For an output array with EVEN size n, 
    the PSF center will be in the corner between pixel (n/2-1,n/2-1) and (n/2,n/2)

    Those coordinates all assume IDL or Python style pixel coordinates running from
    (0,0) up to (n-1, n-1). 

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
    offset
        an offset in pixels relative to the above

    """



    npup = pupil.shape[0]

    du = nlamD / float(npix)
    dv = nlamD / float(npix)

    dx = 1.0/float(npup)
    dy = 1.0/float(npup)

    Xs = (np.arange(npup) - float(npup)/2.0 - offset[1] + 0.5) * dx
    Ys = (np.arange(npup) - float(npup)/2.0 - offset[0] + 0.5) * dy

    Us = (np.arange(npix) - float(npix)/2.0 - offset[1] + 0.5) * du
    Vs = (np.arange(npix) - float(npix)/2.0 - offset[0] + 0.5) * dv

    XU = np.outer(Xs, Us)
    YV = np.outer(Ys, Vs)


    expXU = np.exp(-2.0 * np.pi * 1j * XU)
    expYV = np.exp(-2.0 * np.pi * 1j * YV)
    expXU = expXU.T.copy()

    t1 = np.dot(expXU, pupil)
    t2 = np.dot(t1, expYV)

    #return t2 * dx * dy
    return  float(nlamD)/(npup*npix) *   t2 



# ROTATABLE 
def SFT4(pupil, nlamD, npix, offset=(0.0,0.0), angle=0.0):
    """
    Compute an adjustible-center, rotatable matrix fourier transform. 

    For an output array with ODD size n,
    the PSF center will be at the center of pixel (n-1)/2
    
    For an output array with EVEN size n, 
    the PSF center will be in the corner between pixel (n/2-1,n/2-1) and (n/2,n/2)

    Those coordinates all assume IDL or Python style pixel coordinates running from
    (0,0) up to (n-1, n-1). 

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
    offset
        an offset in pixels relative to the above

    """

    rotation = 45.0
    cosr = np.cos( np.radians(rotation))
    sinr = np.sin( np.radians(rotation))


    npup = pupil.shape[0]

    du = nlamD / float(npix)
    dv = nlamD / float(npix)

    dx = 1.0/float(npup)
    dy = 1.0/float(npup)

    # Xs and Ys are, unsurprisingly, the X and Y coordinates for the pupil array.
    # I think? Though then why does it make it work properly to shift them here like this?
    # Actually it should make NO difference. A shift in the pupil plane just induces a phase tilt
    # in the image plane, which we don't care about for PSFs since we just measure total intensity.
      #Yep, confirmed numerically with some tests. Applying shifts here is unnecessary!
        #Xs = (np.arange(npup) - float(npup)/2.0 - offset[1] + 0.5) * dx
        #Ys = (np.arange(npup) - float(npup)/2.0 - offset[0] + 0.5) * dy

    Xs = (np.arange(npup) - float(npup)/2.0 ) * dx
    Ys = (np.arange(npup) - float(npup)/2.0 ) * dy



    # OK, a 2D FFT can be computed as a the result of two separate 1D transforms...


    # Aaaargh this is not going to work.
    #
    Us = (np.arange(npix) - float(npix)/2.0 - offset[1] + 0.5) * du
    Vs = (np.arange(npix) - float(npix)/2.0 - offset[0] + 0.5) * dv
    Us.shape = (1,npix)
    Vs.shape = (npix,1)

    UsR =  cosr*Us + sinr*Vs
    VsR = -sinr*Us + cosr*Vs


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
    choice : string
        Either 'SYMMETRIC', 'FFTSTYLE', or 'ADJUSTIBLE'. 
        Sets whether the DFT result is centered at pixel n/2+1 (FFTSTYLE) 
        or on the crosshairs between the central pixels (SYMMETRIC),
        or exactly centered in the array no matter what (ADJUSTIBLE). Default is FFTSTYLE. 


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

        self.choices = ("FFTSTYLE", "SYMMETRIC", "ADJUSTIBLE")
        self.correctoffset = {self.choices[0]: 0.5, self.choices[1]: 0.0}
        if choice not in self.choices:
            raise ValueError("Error: choice must be one of " % self.choices)
        self.choice = choice

        fns = {'FFTSTYLE':SFT1, "SYMMETRIC":SFT2, "ADJUSTIBLE":SFT3}

        if self.verbose:
            #print choice
            #print "Announcement  - This instance of SlowFourierTransform uses SFT2"
            print "This instance of SFT is a(n) %s  set-up calling %s " % (choice, fns[choice])
        self.perform = fns[choice]

    def offset(self):
            return self.correctoffset[self.choice]


    def performFITS(hdulist, focalplane_size, focalplane_npix):
        """ Perform an MFT, and return the result as a pyfits.HDUlist """
        newHDUlist = hdulist.copy()
        newim = self.perform(hdulist[0].data, focalplane_size, focalplane_npix)

        newHDUlist[0].data = newim
        #TODO fits header keyword updates

        return newHDUlist


#---------------------------------------------------------------------
#  Test functions 

def euclid2(s, c=None):

	if c is None:
		c = (0.5*float(s[0]),  0.5*float(s[1]))

	y, x = N.indices(s)
	r2 = (x - c[0])**2 + (y - c[1])**2

	return r2

def makedisk(s=None, c=None, r=None, inside=1.0, outside=0.0, grey=None, t=None):
	
	# fft style or sft asymmetric style - center = nx/2, ny/2
	# see ellipseDriver.py for details on symm...

	disk = N.where(euclid2(s, c=c) <= r*r, inside, outside)
	return disk



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
    pupil = makedisk(s=s, c=ctr, r=float(npupil)/2.0001, t=np.float64, grey=0)

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



def test_SFT_center( npix=100, outdir='.', outname='SFT1'):
    choice='ADJUSTIBLE'
    import os

    npupil = 156
    pctr = int(npupil/2)
    npix = 1024
    u = 100    # of lam/D
    s = (npupil,npupil)


    # FFT style
    sft1 = SlowFourierTransform(choice=choice)

    ctr = (float(npupil)/2.0 + sft1.offset(), float(npupil)/2.0 + sft1.offset())
    #print ctr
    pupil = makedisk(s=s, c=ctr, r=float(npupil)/2.0001, t=np.float64, grey=0)

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
