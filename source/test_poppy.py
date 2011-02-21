from jwopt import *
import jwopt
import poppy
import numpy as N
import pyfits
import unittest
import matplotlib.pyplot as P

import optics

import logging
_log = logging.getLogger('tester')
_log.addHandler(logging.NullHandler())


""" Testing code for Poppy & JWST PSF models

    This is a suite of test cases for the various Poppy optical propagation tools, 
    right now with emphasis on verifying coronagraphy & FQPMs. 
    
    This is written using Python's unittest framework. This is my first attempt using
    unittest, and it seems somewhat convoluted in terms of setup, particularly since there
    is no easy way to iterate over test cases with variable arguments like wavelength.
    So there's some work-arounds using global variables here.


    Main routines to look at:  test_run() and test_multiwave() 

    test_run is a convenience function to let you run all the numbered test suites in order,
    or to pick just one at a time. 

    MDP 2011-02-10

"""


def basic_test():
    jwst = JWST_OTE("path/to/some/OPDs")
    miri = jwst.MIRI
    gstar = pysynphot('G2 star')

    psf2 = miri.psf('imaging', center=(512,512), filter='F1500W', oversample=4, spectrum=gstar)

    corPSF = miri.psf('lyot', filter='F2550W', decenter=0.01, oversample=4)


def check_wavefront(filename, slice=0, test='nearzero', comment=""):
    """ A helper routine to verify certain properties of a wavefront FITS file, 
    as requested by some test routine. """
    imstack = pyfits.getdata(filename)
    im = imstack[slice,:,:]

    try:

        if test=='nearzero':
            assert(  N.all(N.abs(im) < N.finfo(im.dtype).eps*10))
            _log.info("Slice %d of %s %s is all essentially zero" % (slice, filename, comment))
            return True
        elif test == 'is_real':
            #assumes output type = 'all'
            cplx_im = imstack[1,:,:] * N.exp(1j*imstack[2,:,:])
            assert(  N.all( cplx_im.imag < N.finfo(im.dtype).eps*10))
            _log.info("File %s %s is essentially all real " % (filename, comment))
            return True

    except:
        _log.error("Test %s failed for %s " % (test, filename))
        return False




####################################################33
#
#  Test Cases
#


_TEST_WAVELENGTH = 2e-6
_TEST_OVERSAMP = 2

class TestPoppy(unittest.TestCase):
    """ Base test case that allows setting wavelength and oversampling for the tests, 
    and cleans up temp files afterwards """
    # Note: you CANNOT override the __init__ method for a TestCase - this behaves in some
    # semi-horrible way due to complexities of the guts of how unittest works. 
    #
    def setUp(self):
        self.wavelength = _TEST_WAVELENGTH
        self.oversample = _TEST_OVERSAMP
        self.pixelscale = 0.010
        self.test_wavelengths = [1e-6, 2e-6, 3e-6]
    def tearDown(self):
        """ Clean up after tests """
        if 0: 
            os.delete('wavefront_plane*fits') # delete intermediate wavefront files
            os.delete('test*fits')

    def iter_wavelengths(self, function, *args, **kwargs):
        """ Iterate over the set of wavelengths to test, repeatedly calling a given function 
        This is implemented as a generator function, to allow the calling routine to examine intermediate
        data products before the next step of the iteration """

        for wl in self.test_wavelengths:
            _log.info("---- Testing wavelength = %e " % wl)
            self.wavelength=wl
            yield function(*args, **kwargs)

class Test_FOV_size(TestPoppy):
    def test_0_fov_size(self):
        """ Test the PSF normalization """

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil(function='Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_pixels=100, oversample=1)
        
        psf = osys.calcPSF(wavelen=self.wavelength, normalize=norm)

        self.assertEqual(psf[0].data.shape[0], 100)


class Test1(TestPoppy):
    def do_test_1_normalization(self):
        """ Test the PSF normalization """

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil(function='Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=10.0)
        
        for norm in ['first', 'last']:
            P.clf()
            psf = osys.calcPSF(wavelen=self.wavelength, normalize=norm)
            jwopt.display_psf(psf)
            tot = psf[0].data.sum()
            _log.info("Using normalization method=%s, the PSF total is\t%f" % (norm, tot))
            if norm =='last':
                self.assertAlmostEqual(abs(tot), 1 ) # the PSF's total on the computed array should be 1, or very close to it.
            else: 
                self.assertAlmostEqual(abs(tot), 1, delta=0.01 )  # this should be very roughly 1.
    def test_1_normalization_multiwave(self):
        results = [res for res in self.iter_wavelengths(self.do_test_1_normalization)]

class Test2(TestPoppy):
    def test_2_lyotreal(self):
        """ Test  no FQPM, no field mask. Verify proper behavior in Lyot plane"""

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addImage()  # perfect image plane
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)
        
        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = False
        poppy._USE_FFTW3 = True
        psf = osys.calcPSF(wavelen=self.wavelength, save_intermediates=True, display_intermediates=True)
        psf.writeto('test2_psf.fits', clobber=True)

        # after the Lyot plane, the wavefront should be all real. 
        self.assertTrue(check_wavefront('wavefront_plane_002.fits', test='is_real', comment='(Lyot Plane)'))
    def test_multiwave(self):
        self.iter_wavelengths(self.test_2_lyotreal)

class Test3(TestPoppy):
    def test_3_fqpm_tilt(self):
        """ Test FQPM tilting (no FQPM yet), no field mask. Verify proper behavior in Lyot plane"""

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addPupil('FQPM_FFT_aligner')
        osys.addImage()  # perfect image plane
        osys.addPupil('FQPM_FFT_aligner', direction='backward')
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)
            #TODO testing of odd and even focal plane sizes?
        
        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = False
        psf = osys.calcPSF(wavelen=self.wavelength, save_intermediates=True, display_intermediates=True)
        psf.writeto('test3a_psf.fits', clobber=True)

        # after the Lyot plane, the wavefront should be all real. 
        check_wavefront('wavefront_plane_004.fits', test='is_real', comment='(Lyot Plane)')

        cen = jwopt.measure_center('wavefront_plane_002.fits', boxsize=50)
        head = pyfits.getheader('wavefront_plane_002.fits')
        desired_pos = (head['NAXIS1']-1)/2.0
        self.assertAlmostEqual( cen[0], desired_pos, delta=0.025) #within 1/50th of a pixel of desired pos?
        self.assertAlmostEqual( cen[1], desired_pos, delta=0.025) #within 1/50th of a pixel of desired pos?
                # This is likely dominated by uncertainties in the simple center measuring algorithm...

        _log.info("FQPM FFT half-pixel tilting is working properly in intermediate image plane")

        cen2 = jwopt.measure_center('wavefront_plane_005.fits', boxsize=50)
        head2 = pyfits.getheader('wavefront_plane_005.fits')
        desired_pos2 = (head2['NAXIS1']-1)/2.0
        self.assertAlmostEqual( cen2[0], desired_pos2, delta=0.05) #within 1/20th of a pixel of desired pos?
                                    


    def do_test_3_ideal_fqpm(self):
        """ Test  ideal FQPM, no field mask. Verify proper behavior in Lyot plane"""

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addPupil('FQPM_FFT_aligner')
        osys.addImage('FQPM', wavelength=self.wavelength)  # perfect FQPM for this wavelength
        osys.addPupil('FQPM_FFT_aligner', direction='backward')
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)
        
        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = False
        #poppy._USE_FFTW3 = True
        #logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')
        psf = osys.calcPSF(wavelen=self.wavelength, save_intermediates=True, display_intermediates=True)
        psf.writeto('test3_psf.fits', clobber=True)

        # after the Lyot plane, the wavefront should be all real. 
        check_wavefront('wavefront_plane_004.fits', test='is_real', comment='(Lyot Plane)')

    def test_3_multiwave_multiwave(self):
        """ Verify that the fluxes at the lyot planes are small and independent of wavelength"""
        lyot_fluxes = []
        for result in self.iter_wavelengths(self.do_test_3_ideal_fqpm):
            # The Lyot plane flux should be indep of wavelength
            im = pyfits.getdata('wavefront_plane_004.fits')
            lyot_fluxes.append(im[0,:,:].sum())
            self.assertLess(lyot_fluxes[-1], 0.005) 
        _log.info("Lyot plane fluxes : "+str(lyot_fluxes))

        for i in range(len(lyot_fluxes)-1):
            self.assertAlmostEqual(lyot_fluxes[i], lyot_fluxes[i+1])
        _log.info("Lyot plane is independent of wavelength. ")

class Test4(TestPoppy):
    def do_test_4_source_offset_0(self):
        do_4_source_offset(angle=.00)
    #def do_test_4_source_offset_45(self):
        #do_4_source_offset(angle=45.0)

    def do_4_source_offset(self, angle=0):
        #oversample=2, verbose=True, wavelen=2e-6, angle=0):
        """ Perfect circular case  no FQPM no field mask, off-axis PSF location
        
        Test point source shifting. no FQPM, no field mask. Verify proper behavior in Lyot plane"""

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addImage()
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)
        
        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = False

        osys.source_offset_theta = angle
        for i in range(15):
            osys.source_offset_r = i * 0.1
            psf = osys.calcPSF(wavelen=self.wavelength, display_intermediates=True)
            #psf.writeto('test3_psf.fits', clobber=True)
            # TODO check position

        # after the Lyot plane, the wavefront should be all real. 
        #check_wavefront('wavefront_plane_004.fits', test='is_real')

    def test_multiwave_offsets(self):
        self.iter_wavelengths(self.do_4_source_offsets, angle=0.0)
        self.iter_wavelengths(self.do_4_source_offsets, angle=45.0)

class Test5(TestPoppy):
    def test_5(self):
        """ Perfect circular case  no FQPM with fieldMask
        
        Test  ideal FQPM, with field mask. Verify proper behavior in Lyot plane"""

        poppy._IMAGECROP = 30 # plot images w/out zooming in on just the center.

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addPupil('FQPM_FFT_aligner')
        osys.addImage() 
        osys.addImage('fieldstop', size=20.0) # do we need to worry about the half-pixel offset here?
        osys.addPupil('FQPM_FFT_aligner', direction='backward')
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)
        
        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = False
        #poppy._USE_FFTW3 = True
        #logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')
        psf = osys.calcPSF(wavelen=self.wavelength, display_intermediates=True)

class Test6(TestPoppy):
    def test_6(self): #oversample=2, verbose=True, wavelen=2e-6):
        """ Perfect circular case  with FQPM with fieldMask
        
        Test  ideal FQPM, with field mask. Verify proper behavior in Lyot plane"""

        poppy._IMAGECROP = 30 # plot images w/out zooming in on just the center.

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addPupil('FQPM_FFT_aligner')
        osys.addImage('FQPM', wavelength=wavelen)  # perfect FQPM for this wavelength
        osys.addImage('fieldstop', size=20.0)  
        osys.addPupil('FQPM_FFT_aligner', direction='backward')
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)
        
        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = False
        #poppy._USE_FFTW3 = True
        #logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')
        psf = osys.calcPSF(wavelen=self.wavelength, display_intermediates=True)

def test_7(oversample=2, verbose=True, wavelen=10.65e-6, radius=0.0):
    """ Perfect circular case  with FQPM with fieldMask off-axis
    
    Test  ideal FQPM, with field mask. Verify proper behavior in Lyot plane"""

    poppy._IMAGECROP = 5 # plot images w/out zooming in on just the center.
    oversample = 2
    pixelscale = 0.1

    osys = poppy.OpticalSystem("test", oversample=oversample)

    osys.addPupil('Circle', radius=6.5/2)
    osys.addPupil('FQPM_FFT_aligner')
    osys.addImage('FQPM', wavelength=wavelen)  # perfect FQPM for this wavelength
    osys.addImage('fieldstop', size=20.0)  
    osys.addPupil('FQPM_FFT_aligner', direction='backward')
    osys.addPupil('Circle', radius=6.5/2)
    osys.addDetector(pixelscale=pixelscale, fov_arcsec=3.0)


    osys.source_offset_r = radius
    osys.source_offset_theta = 45.
    
    P.clf()
    poppy._FLUXCHECK=True
    poppy._USE_FFTW3 = False
    #poppy._USE_FFTW3 = True
    #logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')
    psf = osys.calcPSF(wavelen=wavelen, display_intermediates=True, save_intermediates=True)

    return osys

def test_8(oversample=2, verbose=True, wavelen=2e-6, angle=0):
    """ Test point source shifting, given variable pupil array padding. """

    for pad_factor in [1.0, 1.1, 1.5, 2.0]: 
        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2, pad_factor = pad_factor)
        osys.addImage()
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)
        
        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = False


        osys.source_offset_theta = angle
        osys.source_offset_r =  0.5
        psf = osys.calcPSF(wavelen=self.wavelength, display_intermediates=True)
            #psf.writeto('test3_psf.fits', clobber=True)
            # TODO check position

        # after the Lyot plane, the wavefront should be all real. 
        #check_wavefront('wavefront_plane_004.fits', test='is_real')

def test_blc2(oversample=2, verbose=True, wavelen=2e-6, angle=0, kind='circular', sigma=1.0, loc = 0.3998):
    import scipy
    x = N.linspace(-5, 5, 401)
    sigmar = sigma*x
    if kind == 'circular':
        trans = (1-  (2*scipy.special.jn(1,sigmar)/sigmar)**2)**2
    else: 
        trans = (1-  (N.sin(sigmar)/sigmar)**2)**2
    P.clf()
    P.plot(x, trans)
    P.axhline(0.5, ls='--', color='k')


    P.axvline(loc, ls='--', color='k')
    #P.gca().set_xbound(loc*0.98, loc*1.02)
    wg = N.where(sigmar > 0.01)
    intfn = scipy.interpolate.interp1d(x[wg], trans[wg])
    print "Value at %.4f :\t%.4f" % (loc, intfn(loc))
    #stop()

def width_blc(desired_width, approx=None, plot=False):
    """ The calculation of sigma parameters for the wedge BLC function is not straightforward.

    This function numerically solves the relevant equation to determine the sigma required to 
    acheive a given HWHM.

    It uses recursion to iterate to a higher precision level.
    """

    loc = desired_width

    if approx is None:
        sigma = N.linspace(0, 20, 5000)
    else: 
        sigma = N.linspace(approx*0.9, approx*1.1, 100000.)
    lhs = loc* N.sqrt(1 - N.sqrt(0.5))
    rhs = N.sin(sigma * loc) / sigma
    diff = N.abs(lhs - rhs)
    wmin = N.where(diff == N.nanmin(diff))
    sig_ans = sigma[wmin][0]

    if approx: 
        return sig_ans
    else:
        # use recursion
        sig_ans = width_blc(loc, sig_ans)

    if plot:
        check =  (1-  (N.sin(sig_ans * loc)/sig_ans/loc)**2)**2
        #P.plot(sigma, lhs)
        P.clf()
        P.plot(sigma, rhs)
        P.axhline(lhs)

        print "sigma = %f implies HWHM = %f" % (sig_ans, loc)
        print " check: 0.5 == %f" % (check)
    return sig_ans


def calc_blc_wedge(deg=4, wavelen=2.1e-6):
    """ This function determines the desired sigma coefficients required to 
    achieve a wedge from 2 to 6 lam/D.

    It returns the coefficients of a polynomial fit that maps from
    nlambda/D to sigma. 

    """
    import scipy
    r = N.linspace(2, 6, 161)
    difflim = wavelen / 6.5 * 180.*60*60/N.pi 
    sigs = [width_blc(difflim * ri) for ri in r]

    pcs = scipy.polyfit(r, sigs, deg)
    p = scipy.poly1d(pcs)
    P.plot(r, sigs, 'b')
    P.plot(r, p(r), "r--")
    diffs = (sigs - p(r))
    print "Poly fit:" +repr(pcs)
    print "  fit rms: "+str(diffs.std())



def test_blc(oversample=2, verbose=True, wavelength=2.1e-6, angle=0, kind='nircamcircular'):

    if wavelength == 2.1e-6: sigma=5.253
    elif wavelength == 3.35e-6: sigma=3.2927866
    else: sigma=2.5652

    blc = poppy.BandLimitedCoron('myBLC', kind, sigma=sigma, wavelength=wavelength)
    #wf = poppy.Wavefront( wavelength = wavelen, pixelscale = 0.010, npix=512)
    #blcphase = blc.getPhasor(wf)


    P.clf()
    #P.subplot(121)
    blc.display()

    #stop()


def test_nc_corons(oversample=2, verbose=True, wavelen=2e-6, angle=0):
    """ Test point source shifting, given variable pupil array padding. """

    oversample=2
    pixelscale = 0.010

    osys = poppy.OpticalSystem("test", oversample=oversample)
    osys.addPupil('Circle', radius=6.5/2)
    osys.addImage()
    osys.addImage('BandLimitedCoron', 'circular',  sigma=1.0)
    osys.addPupil('Circle', radius=6.5/2)
    osys.addDetector(pixelscale=pixelscale, fov_arcsec=3.0)
    
    P.clf()
    poppy._FLUXCHECK=True
    poppy._USE_FFTW3 = False


    osys.source_offset_theta = angle
    osys.source_offset_r =  0.0
    psf = osys.calcPSF(wavelen=wavelen, display_intermediates=True)
            #psf.writeto('test3_psf.fits', clobber=True)
            # TODO check position

        # after the Lyot plane, the wavefront should be all real. 
        #check_wavefront('wavefront_plane_004.fits', test='is_real')




class TestPupils(TestPoppy):
    """ Test circular, square, hexagonal pupils and their PSFs """
    def test_pupils(self):

        def image_cut_1d(image, angle=0):
            """ Make a quick 1D cut through an image starting at the center """
            #y, x = N.indices(image)
            #y-= (image.shape[0]-1)/2
            #x-= (image.shape[1]-1)/2

            t = N.arange(image.shape[0])
            cx = N.cos(angle*N.pi/180)*t +  (image.shape[0]-1)/2
            cy = N.sin(angle*N.pi/180)*t +  (image.shape[1]-1)/2
            cx = N.asarray(N.round(cx), dtype=int)
            cy = N.asarray(N.round(cy), dtype=int)

            wg = N.where( (cx >=0) & (cy >=0) & (cx < image.shape[1]) & (cy < image.shape[0]))
            return image[cy[wg],cx[wg]]


        pupils = ['Circle', 'Hexagon', 'Square']
        angles = [[0], [0, 30], [0, 45]]
        effective_diams = [[2], [2,N.sqrt(3)], [2,2*N.sqrt(2)]]

        P.clf()
        cuts = []
        for i in range(3):
            P.subplot(2,3, i+1)
            osys = poppy.OpticalSystem("test", oversample=self.oversample)
            osys.addPupil(pupils[i])
            osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)

            psf = osys.calcPSF(wavelen=self.wavelength)

            jwopt.display_psf(psf)
            P.title(pupils[i])

            P.subplot(2,3, i+4)
            for ang, diam in zip(angles[i], effective_diams[i]):
                cut = image_cut_1d(psf[0].data, ang)
                r = N.arange(cut.size) * 0.010
                cut /= cut.max() # normalize to peak=1
                P.semilogy(r, cut, label='$\\theta = %d^o$' % ang )
                if i == 0:
                    radius, airyfn = optics.airy_1d(diam, self.wavelength, pixelscale=self.pixelscale)
                    P.plot(radius, airyfn, "k--", label='analytic')
            P.gca().set_xbound(0,3)
            P.gca().set_ybound(1e-13,1.5)

            
            # TODO - overplot perfect analytical PSFs.
            P.legend(loc='upper right', frameon=False)

        self.assertTrue(True)
        #FIXME - compare cuts to airy function etc.



def test_run(index=None, wavelength=2e-6):
    """ This function provides a simple interface for running all available tests, or just one """
    #tests = [Test1]
    global _TEST_WAVELENGTH
    _TEST_WAVELENGTH = wavelength
    tests = [TestPupils, TestPoppy, Test1, Test2, Test3, Test4, Test5]

    if index is not None:
        if not hasattr(index, '__iter__') : index=[index]
        tests = [tests[i] for i in index]

    suite = unittest.TestSuite()
    for t in tests:
        suite.addTest( unittest.TestLoader().loadTestsFromTestCase( t) )
        #suite.addTest( t()  )
    unittest.TextTestRunner(verbosity=2).run(suite)



def test_multiwave(*args):
    global _TEST_WAVELENGTH
    for w in [1e-6, 2e-6, 5e-6]:
        print("*"*70)
        print("  Running tests with wavelength = %e m" % w)
        test_run(wavelength=w, *args)



if __name__== "__main__":
    logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')

    #reload poppy
    #reload jwopt
    #nc = jwopt.NIRCam()
    #nc.pupilopd=None
    #osys = nc.getOpticalSystem()
    pass



