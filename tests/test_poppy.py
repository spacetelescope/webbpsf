from webbpsf import *
import webbpsf
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


def check_wavefront(filename_or_hdulist, slice=0, ext=0, test='nearzero', comment=""):
    """ A helper routine to verify certain properties of a wavefront FITS file, 
    as requested by some test routine. """
    if isinstance(filename_or_hdulist, str):
        hdulist = pyfits.open(filename_or_hdulist)
        filename = filename_or_hdulist
    elif isinstance(filename_or_hdulist, pyfits.HDUList):
        hdulist = filename_or_hdulist
        filename = 'input HDUlist'
    imstack = hdulist[ext].data
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
        """ Test the PSF field of view size"""

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil(function='Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_pixels=100, oversample=1)
        
        psf = osys.calcPSF(wavelength=self.wavelength, normalize=norm)

        self.assertEqual(psf[0].data.shape[0], 100)



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

            psf = osys.calcPSF(wavelength=self.wavelength)

            webbpsf.display_psf(psf)
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



class Test1(TestPoppy):
    def do_test_1_normalization(self):
        """ Test the PSF normalization """

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil(function='Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=10.0)
        
        for norm in ['first', 'last']:
            P.clf()
            psf = osys.calcPSF(wavelength=self.wavelength, normalize=norm)
            webbpsf.display_psf(psf)
            tot = psf[0].data.sum()
            _log.info("Using normalization method=%s, the PSF total is\t%f" % (norm, tot))
            if norm =='last':
                self.assertAlmostEqual(abs(tot), 1 ) # the PSF's total on the computed array should be 1, or very close to it.
            else: 
                self.assertAlmostEqual(abs(tot), 1, delta=0.01 )  # this should be very roughly 1.
    def test_1_normalization_multiwave(self):
        results = [res for res in self.iter_wavelengths(self.do_test_1_normalization)]

class Test2(TestPoppy):
    """ Is the wave in the Lyot plane essentially all real? i.e. negligible imaginary part """
    def test_2_lyotreal_numpyfft(self):
        poppy._USE_FFTW3 = False
        self.do_test_2_lyotreal()
    def test_2_lyotreal_fftw(self):
        poppy._USE_FFTW3 = True
        self.do_test_2_lyotreal()
    def do_test_2_lyotreal(self):
        """ Test  no FQPM, no field mask. Verify proper behavior in Lyot plane"""
        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addImage()  # perfect image plane
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)
        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = True
        psf = osys.calcPSF(wavelength=self.wavelength, save_intermediates=True, display_intermediates=True)
        psf.writeto('test2_psf.fits', clobber=True)

        # after the Lyot plane, the wavefront should be all real. 
        self.assertTrue(check_wavefront('wavefront_plane_002.fits', test='is_real', comment='(Lyot Plane)'))
    def test_multiwave(self):
        self.iter_wavelengths(self.test_2_lyotreal_fftw)

class Test3(TestPoppy):
    """ First, verify the FQPM tilt behavior works as desired. 
        Then, test an ideal FQPM  """

    def test_3_fqpm_tilt_numpyfft(self):
        poppy._USE_FFTW3 = False
        self.do_test_3_fqpm_tilt()
    def test_3_fqpm_tilt_fftw(self):
        poppy._USE_FFTW3 = True
        self.do_test_3_fqpm_tilt()
    def test_3_ideal_fqpm_numpyfft(self):
        poppy._USE_FFTW3 = False
        self.do_test_3_ideal_fqpm()
    def test_3_ideal_fqpm_fftw(self):
        poppy._USE_FFTW3 = True
        self.do_test_3_ideal_fqpm()
         
    def do_test_3_fqpm_tilt(self):
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
        psf = osys.calcPSF(wavelength=self.wavelength, save_intermediates=True, display_intermediates=True)
        psf.writeto('test3a_psf.fits', clobber=True)

        # after the Lyot plane, the wavefront should be all real. 
        check_wavefront('wavefront_plane_004.fits', test='is_real', comment='(Lyot Plane)')

        cen = webbpsf.measure_centroid('wavefront_plane_002.fits', boxsize=50)
        head = pyfits.getheader('wavefront_plane_002.fits')
        desired_pos = (head['NAXIS1']-1)/2.0
        self.assertAlmostEqual( cen[0], desired_pos, delta=0.025) #within 1/50th of a pixel of desired pos?
        self.assertAlmostEqual( cen[1], desired_pos, delta=0.025) #within 1/50th of a pixel of desired pos?
                # This is likely dominated by uncertainties in the simple center measuring algorithm...

        _log.info("FQPM FFT half-pixel tilting is working properly in intermediate image plane")

        cen2 = webbpsf.measure_centroid('wavefront_plane_005.fits', boxsize=50)
        head2 = pyfits.getheader('wavefront_plane_005.fits')
        desired_pos2 = (head2['NAXIS1']-1)/2.0
        self.assertAlmostEqual( cen2[0], desired_pos2, delta=0.05) #within 1/20th of a pixel of desired pos?
                                    
        _log.info("FQPM FFT half-pixel tilting is working properly in final image plane")


    def do_test_3_ideal_fqpm(self):
        """ Test  ideal FQPM, no field mask. Verify proper behavior in Lyot plane"""


        #self.wavelength = 8e-6 # for ease of seeing details on screen make it bigger
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
        psf, int_wfs = osys.calcPSF(wavelength=self.wavelength, save_intermediates=False, display_intermediates=True, return_intermediates=True)
        #psf.writeto('test3_psf.fits', clobber=True)
        lyot_wf = int_wfs[-2]
        lyot_wf.writeto("wavefront_plane_004.fits", what='all', clobber=True) # need to save this for the multiwave comparison in test_3_multiwave()

        # after the Lyot plane, the wavefront should be all real. 
        self.assertTrue(check_wavefront(lyot_wf.asFITS(what='all'), test='is_real', comment='(Lyot Plane)'))
        self.assertLess(psf[0].data.sum(), 0.002) 
        _log.info("post-FQPM flux is appropriately low.")

    def test_3_multiwave(self):
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
    """ Verify ability to shift point sources.
    The Lyot plane should still be all real. """
    def test_offsets_numpyfft(self):
        poppy._USE_FFTW3 = False
        self.do_4_source_offsets()
    def test_offsets_fftw(self):
        poppy._USE_FFTW3 = True
        self.do_4_source_offsets()

    def do_4_source_offsets(self, angle=0):
        #oversample=2, verbose=True, wavelength=2e-6, angle=0):
        """ Perfect circular case  no FQPM no field mask, off-axis PSF location

        Test point source shifting. no FQPM, no field mask. Verify proper behavior in Lyot plane"""

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addImage()
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=3.0)

        P.clf()
        poppy._FLUXCHECK=True

        osys.source_offset_theta = angle
        for i in range(15):
            osys.source_offset_r = i * 0.1
            psf = osys.calcPSF(wavelength=self.wavelength, display_intermediates=True)

            pos = webbpsf.fwcentroid(psf[0].data, halfwidth=10, threshhold=1e-2)
            # pos is the pixel coords in y,x in pixel units.
            cenx = (psf[0].data.shape[0]-1)/2.0
            offset = N.sqrt( (pos[0]-cenx)**2 + (pos[1]-cenx)**2)* self.pixelscale/self.oversample
            _log.info("Desired offset is %f, measured is %f " % (osys.source_offset_r, offset))
            self.assertAlmostEqual(osys.source_offset_r, offset, 3)


    def test_multiwave_offsets(self):
        self.iter_wavelengths(self.do_4_source_offsets, angle=0.0)
        self.iter_wavelengths(self.do_4_source_offsets, angle=45.0)

class Test5(TestPoppy):
    """ Test the Field Mask works properly """
    def test_5(self):
        """ Perfect circular case  no FQPM with fieldMask
        
        Test  ideal FQPM, with field mask. Verify proper behavior in Lyot plane"""

        poppy._IMAGECROP = 30 # plot images w/out zooming in on just the center.

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addPupil('FQPM_FFT_aligner')
        osys.addImage() 
        osys.addImage('fieldstop', size=10.0)  # do we need to worry about the half-pixel offset here  A: no, now handed in wavefront.coordinates
        osys.addPupil('FQPM_FFT_aligner', direction='backward')
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=12)
        
        P.clf()
        poppy._FLUXCHECK=True
        #poppy._USE_FFTW3 = False
        #poppy._USE_FFTW3 = True
        #logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')
        psf = osys.calcPSF(wavelength=self.wavelength, display_intermediates=True)

        # TODO need to do some kind of evaluation here!

class Test6(TestPoppy):
    """ Test FQPM plus field mask together """
    #def test

    def test_6(self): #oversample=2, verbose=True, wavelength=2e-6):
        """ Perfect circular case  with FQPM with fieldMask
        Test  ideal FQPM, with field mask. Verify proper behavior in Lyot plane"""

        poppy._IMAGECROP = 5 # plot images w/out zooming in on just the center.

        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=6.5/2)
        osys.addPupil('FQPM_FFT_aligner')
        osys.addImage('FQPM', wavelength=self.wavelength)  # perfect FQPM for this wavelength
        osys.addImage('fieldstop', size=6.0)  
        osys.addPupil('FQPM_FFT_aligner', direction='backward')
        osys.addPupil('Circle', radius=6.5/2)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=10.0)

        P.clf()
        poppy._FLUXCHECK=True
        poppy._USE_FFTW3 = False
        #poppy._USE_FFTW3 = True
        #logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')
        psf = osys.calcPSF(wavelength=self.wavelength, display_intermediates=True)
        self.assertLess(psf[0].data.sum(), 0.002) 
        _log.info("post-FQPM flux is appropriately low.")


class Test7(TestPoppy):
    """ Test the FQPM with field mask, off axis 
    
    When you are off axis the FQPM should not appreciably mess with the PSF.

    """
    def test_7(self): #oversample=2, verbose=True, wavelength=2e-6):
 
    #def test_7(oversample=2, verbose=True, wavelength=10.65e-6, radius=0.0):
        """ Perfect circular case  with FQPM with fieldMask off-axis
        
        Test  ideal FQPM, with field mask. Verify proper behavior in Lyot plane"""

        poppy._IMAGECROP = 5 # plot images w/out zooming in on just the center.
        #oversample = 2
        #pixelscale = 0.1

        fov = 6

        osys1 = poppy.OpticalSystem(" no FQPM, offset", oversample=self.oversample)
        osys1.addPupil('Circle', radius=6.5/2)
        osys1.addPupil('FQPM_FFT_aligner')
        osys1.addImage('fieldstop', size=20.0)  
        osys1.addPupil('FQPM_FFT_aligner', direction='backward')
        osys1.addPupil('Circle', radius=6.5/2)
        osys1.addDetector(pixelscale=self.pixelscale, fov_arcsec=fov)

        osys2 = poppy.OpticalSystem("FQPM offset", oversample=self.oversample)
        osys2.addPupil('Circle', radius=6.5/2)
        osys2.addPupil('FQPM_FFT_aligner')
        osys2.addImage('FQPM', wavelength=self.wavelength)  # perfect FQPM for this wavelength
        osys2.addImage('fieldstop', size=20.0)  
        osys2.addPupil('FQPM_FFT_aligner', direction='backward')
        osys2.addPupil('Circle', radius=6.5/2)
        osys2.addDetector(pixelscale=self.pixelscale, fov_arcsec=fov)

        myoffset = 3.0
        osys1.source_offset_r = myoffset
        osys1.source_offset_theta = 45.
        osys2.source_offset_r = myoffset
        osys2.source_offset_theta = 45.
 
        poppy._FLUXCHECK=True
        P.figure(1)
        P.clf()
        psf1 = osys1.calcPSF(wavelength=self.wavelength, display_intermediates=True, save_intermediates=False)
        P.figure(2)
        P.clf()
        psf2 = osys2.calcPSF(wavelength=self.wavelength, display_intermediates=True, save_intermediates=False)

        P.figure(3)
        P.subplot(211)
        webbpsf.display_psf(psf1, title=osys1.name)
        P.subplot(212)
        webbpsf.display_psf(psf2, title=osys2.name)

        pos1 = webbpsf.fwcentroid(psf1[0].data, halfwidth=10, threshhold=1e-2)
        pos2 = webbpsf.fwcentroid(psf2[0].data, halfwidth=10, threshhold=1e-2)

        rel_offset = N.sqrt(((N.array(pos1) - N.array(pos2))**2).sum())
        self.assertTrue(rel_offset < 1e-3 ) 
        _log.info("Source position does not appear to be affected by FQPMs for far off-axis sources")



class Test8(TestPoppy):
    "Verify that extra padding around the aperture makes no difference "
    def test_padding_numpyfft(self):
        poppy._USE_FFTW3=False
        self.do_test_8()
    def test_padding_fftw(self):
        poppy._USE_FFTW3=True
        self.do_test_8()

    def do_test_8(self):
        """ Test point source shifting, given variable pupil array padding. """

        angle = 36.
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
            psf = osys.calcPSF(wavelength=self.wavelength, display_intermediates=True)
            #psf.writeto('test3_psf.fits', clobber=True)
            # TODO check position
            pos = webbpsf.fwcentroid(psf[0].data, halfwidth=10, threshhold=1e-2)
            # pos is the pixel coords in y,x in pixel units.
            cenx = (psf[0].data.shape[0]-1)/2.0
            offset = N.sqrt( (pos[0]-cenx)**2 + (pos[1]-cenx)**2)* self.pixelscale/self.oversample
            _log.info("Desired offset is %f, measured is %f " % (osys.source_offset_r, offset))
            self.assertAlmostEqual(osys.source_offset_r, offset, 3)
 
        # after the Lyot plane, the wavefront should be all real. 
        #check_wavefront('wavefront_plane_004.fits', test='is_real')

class Test9(TestPoppy):
    "Test BLC corons. Verify reasonable on- and off- axis behavior. "
    def test_9_circ(self):
        self.do_test_9(kind='circular')
    def test_9_linear(self):
        self.do_test_9(kind='linear')
    def test_9_circ_offset(self):
        self.do_test_9(kind='circular', offset=True)
    def test_9_linear_offset(self):
        self.do_test_9(kind='linear', offset=True)
  
    def do_test_9(self, kind='circular', offset=False):
        _log.info("Testing BLC kind = "+kind)
        
        radius = 6.5/2
        lyot_radius = 6.5/2.5
        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=radius)
        osys.addImage()
        osys.addImage('BandLimitedCoron', kind=kind, sigma=5.0)
        osys.addPupil()
        osys.addPupil('Circle', radius=lyot_radius)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=5.0)
        poppy._FLUXCHECK=True


        if offset: 
            osys.source_offset_r =  2.0
        else: 
            osys.source_offset_r =  0.0
        poppy._FLUXCHECK= True
        P.clf()
        psf, int_wfs = osys.calcPSF(wavelength=self.wavelength, display_intermediates=True, return_intermediates=True)

        if offset:
            # the flux should be relatively high
            # with the vast majority of the loss due just to the undersized Lyot stop
            self.assertGreater(psf[0].data.sum(), (lyot_radius/radius)**2 *0.95 ) 
            _log.info("For offset source, post-BLC flux is appropriately high")
        else:
            # after the Lyot plane, the wavefront should be all real. 
            lyot_wf = int_wfs[-2]
            lyot_wf_fits = lyot_wf.asFITS(what='all') # need to save this for the multiwave comparison in test_3_multiwave()
            self.assertTrue(check_wavefront(lyot_wf_fits, test='is_real', comment='(Lyot Plane)'))

            # and the flux should be low.
            self.assertLess(psf[0].data.sum(), 1e-4) 
            _log.info("post-BLC flux is appropriately low.")

 

class Test10(TestPoppy):
    "Test multiwavelength multiprocessor propagation "

    def test_10_multiproc_numpyfft(self):
        poppy._USE_FFTW3 = False
        # multiprocessor and FFTW not sure if they play nicely all the time?
        self.do_test_10_multiproc()


    #def test_10_multiproc_fftw3(self):
        #poppy._USE_FFTW3 = True
        ## multiprocessor and FFTW not sure if they play nicely all the time?
        #self.do_test_10_multiproc()
 
    def do_test_10_multiproc(self):

        # for the sake of a reasonably realistic usage case, we use the BLC setup from Test 9
        kind = 'circular'
        radius = 6.5/2
        lyot_radius = 6.5/2.5
        osys = poppy.OpticalSystem("test", oversample=self.oversample)
        osys.addPupil('Circle', radius=radius)
        osys.addImage()
        osys.addImage('BandLimitedCoron', kind=kind, sigma=5.0)
        osys.addPupil()
        osys.addPupil('Circle', radius=lyot_radius)
        osys.addDetector(pixelscale=self.pixelscale, fov_arcsec=5.0)
        osys.source_offset_r =  1.5 # make the PSF easy to see...
 

        nlam= 6
        source = {'weights': [0.1]*nlam, 'wavelengths': N.linspace(2.0e-6, 3.0e-6, nlam)}

        _log.info("Calculating multiprocess PSF")
        times = []
        times.append(time.time())
        psf2 = osys.calcPSFmultiproc(source)
        times.append(time.time())
        tmulti =  times[-1]-times[-2]
        _log.info(" Time for multiprocessor: %f s " % (tmulti))

        _log.info("Calculating single process PSF")
        times.append(time.time())
        psf1 = osys.calcPSF(source['wavelengths'], source['weights'])
        times.append(time.time())
        tsing =  times[-1]-times[-2]
        _log.info(" Time for single processor: %f s " % (tsing))


        _log.info(" Speedup factor: %f " % (tsing/tmulti))
        P.clf()
        ax = P.subplot(1,3,1)
        webbpsf.display_psf(psf1)
        ax = P.subplot(1,3,2)
        webbpsf.display_psf(psf2)
        ax = P.subplot(1,3,3)

        poppy.imshow_with_mouseover(psf1[0].data - psf2[0].data, ax=ax)

        self.assertTrue(  (psf1[0].data == psf2[0].data).all()  )
        _log.info("Exact same result achieved both ways.")



        

#################################################################################

def test_blc2(oversample=2, verbose=True, wavelength=2e-6, angle=0, kind='circular', sigma=1.0, loc = 0.3998):
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


def calc_blc_wedge(deg=4, wavelength=2.1e-6):
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


def test_blc_corons(oversample=2, verbose=True, wavelength=2e-6, angle=0, kind='circular'):
    """ Test point source shifting, given variable pupil array padding. """

    oversample=2
    pixelscale = 0.010
    wavelength = 4.6e-6
    poppy._IMAGECROP = 5

    osys = poppy.OpticalSystem("test", oversample=oversample)
    osys.addPupil('Circle', radius=6.5/2)
    osys.addImage()
    osys.addImage('BandLimitedCoron', kind=kind,  sigma=5.0)
    osys.addPupil()
    osys.addPupil('Circle', radius=6.5/2.5)
    osys.addDetector(pixelscale=pixelscale, fov_arcsec=3.0)
    
    P.clf()
    poppy._FLUXCHECK=True
    poppy._USE_FFTW3 = False


    osys.source_offset_theta = angle
    osys.source_offset_r =  0.0
    psf = osys.calcPSF(wavelength=wavelength, display_intermediates=True)
            #psf.writeto('test3_psf.fits', clobber=True)
            # TODO check position

        # after the Lyot plane, the wavefront should be all real. 
        #check_wavefront('wavefront_plane_004.fits', test='is_real')


#################################################################################


def test_run(index=None, wavelength=2e-6):
    """ This function provides a simple interface for running all available tests, or just one """
    #tests = [Test1]
    global _TEST_WAVELENGTH
    _TEST_WAVELENGTH = wavelength
    tests = [TestPupils, Test1, Test2, Test3, Test4, Test5, Test6, Test7, Test8, Test9, Test10]

    if index is not None:
        if not hasattr(index, '__iter__') : index = [index]
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
    #reload webbpsf
    #nc = webbpsf.NIRCam()
    #nc.pupilopd=None
    #osys = nc.getOpticalSystem()
    pass



