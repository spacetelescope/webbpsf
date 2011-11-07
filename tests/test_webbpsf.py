import sys, os
if "../lib" not in sys.path: sys.path.append(os.path.abspath('../lib'))
import webbpsf
import webbpsf as jw
import poppy
import unittest
import numpy as np
import matplotlib.pyplot as plt
try:
    __IPYTHON__
    from IPython.Debugger import Tracer; stop = Tracer()
except:
    pass


import logging
_log = logging.getLogger('tester-jw')
_log.addHandler(logging.NullHandler())

poppy._log.setLevel(logging.INFO)


def make_fig_instrument_comparison():
    """ Make the instrument comparison figure for the webbPSF web page.
    Argh I already wrote this once, now I need to write it again 'cause I can't find it...
    """

    # make figure size = array([ 14.325,   2.975]) in inches

    
    plt.clf()
    f, axarr = plt.subplots(1, 6, sharex=True, num=1, sharey=True)

    nc = webbpsf.NIRCam()
    nc.filter='F210M'
    webbpsf.calc_or_load_PSF('NIRCam_demo_F212.fits', nc, oversample=4, rebin=True, fov_arcsec=5)
    webbpsf.display_PSF('NIRCam_demo_F212.fits', ext=1, ax= axarr[0], colorbar=False, title="NIRCam F210M")

    nc = webbpsf.NIRCam()
    nc.filter='F444W'
    webbpsf.calc_or_load_PSF('NIRCam_demo_F444.fits', nc, oversample=4, rebin=True, fov_arcsec=5)
    webbpsf.display_PSF('NIRCam_demo_F444.fits', ext=1, ax= axarr[1], colorbar=False, title="NIRCam F444W")

    ns = webbpsf.NIRSpec()
    ns.filter='F110W'
    webbpsf.calc_or_load_PSF('NIRSpec_demo.fits', ns, oversample=4, rebin=True, fov_arcsec=5)
    webbpsf.display_PSF('NIRSpec_demo.fits', ext=1, ax= axarr[2], colorbar=False, title="NIRspec F110W")

    ni = webbpsf.NIRISS()
    ni.filter='F380M'
    webbpsf.calc_or_load_PSF('NIRISS_demo.fits', ni, oversample=4, rebin=True, fov_arcsec=5)
    webbpsf.display_PSF('NIRISS_demo.fits', ext=1, ax= axarr[3], colorbar=False, title="NIRISS F380M")


    mi = webbpsf.MIRI()
    mi.filter='F1000W'
    webbpsf.calc_or_load_PSF('MIRI_demo.fits', mi, oversample=4, rebin=True, fov_arcsec=5)
    webbpsf.display_PSF('MIRI_demo.fits', ext=1, ax= axarr[4], colorbar=False, title="MIRI F1000W")

    fg = webbpsf.FGS()
    webbpsf.calc_or_load_PSF('FGS_demo.fits', fg, oversample=4, rebin=True, fov_arcsec=5)
    webbpsf.display_PSF('FGS_demo.fits', ext=1, ax= axarr[5], colorbar=False, title="FGS")


    plt.subplots_adjust(hspace=0.01)
    for ax in axarr:
        ax.set_ylim(-2.5, 2.5)
        ax.set_xlim(-2.5, 2.5)


    plt.savefig('fig_instrument_comparison.pdf')
    #plt.tight_layout()


class Test_Image_Size(unittest.TestCase):
    def generic_test(self, iname):

        _log.info("Testing image output sizes for %s " % iname)
        inst = webbpsf.Instrument(iname)
        pxscale = inst.pixelscale
        fov_arcsec = 5.0

        PSF = inst.calcPSF(nlambda=1, fov_pixels = 100, oversample=1)
        self.assertEqual(PSF[0].data.shape[0], 100)

        PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
        fov_pix = int(np.round(fov_arcsec / pxscale))
        self.assertEqual(PSF[0].data.shape[0], fov_pix)

        inst.options['parity'] = 'odd'
        PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
        self.assertTrue( np.remainder(PSF[0].data.shape[0],2) == 1)

        inst.options['parity'] = 'even'
        PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
        self.assertTrue( np.remainder(PSF[0].data.shape[0],2) == 0)

        # odd array, even oversampling = even
        inst.options['parity'] = 'odd'
        PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=2)
        self.assertTrue( np.remainder(PSF[0].data.shape[0],2) == 0)

        # odd array, odd oversampling = odd
        inst.options['parity'] = 'odd'
        PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=3)
        self.assertTrue( np.remainder(PSF[0].data.shape[0],2) == 1)



        
    test_nircam = lambda self : self.generic_test('NIRCam')
    test_miri= lambda self : self.generic_test('MIRI')
    test_nirspec= lambda self : self.generic_test('NIRSpec')
    test_tfi= lambda self : self.generic_test('TFI')
    test_fgs= lambda self : self.generic_test('FGS')


class Test_Source_Offset(unittest.TestCase):
    def do_test_source_offset(self, iname, theta=0.0):

        nc = webbpsf.Instrument(iname)
        nc.pupilopd=None

        nsteps = 3
        oversample = 2

        shift_req = []
        psfs = []

        for i in range(nsteps+1):
            nc.options['source_offset_r'] = i*0.1
            nc.options['source_offset_theta'] = theta
            nc.options['source_offset_r'] = i*nc.pixelscale*5
            shift_req.append(nc.options['source_offset_r'])
            psfs.append(  nc.calcPSF(nlambda=1, oversample=oversample) )

        webbpsf.display_PSF(psfs[0])

        cent0 = np.asarray(webbpsf.measure_centroid(psfs[0]))
        center_pix = (psfs[0][0].data.shape[0]-1)/2.0
        self.assertAlmostEqual(cent0[0], center_pix, 3)
        self.assertAlmostEqual(cent0[1], center_pix, 3)
        _log.info("Center of unshifted image: (%d, %d)" % tuple(cent0))


        for i in range(1, nsteps+1):
            webbpsf.display_PSF(psfs[i])
            cent = webbpsf.measure_centroid(psfs[i])
            rx = shift_req[i] * (-np.sin(theta*np.pi/180))
            ry = shift_req[i] * (np.cos(theta*np.pi/180))
            _log.info("   Shift_requested:\t(%10.3f, %10.3f)" % (rx, ry))
            shift = (cent-cent0) * (nc.pixelscale/oversample)
            _log.info("   Shift_achieved: \t(%10.3f, %10.3f)" % (shift[1], shift[0]))
            self.assertAlmostEqual(rx, shift[1], 3)
            self.assertAlmostEqual(ry, shift[0], 3)

    test_nircam_00 = lambda self : self.do_test_source_offset('NIRCam', theta=0.0)
    test_nircam_45 = lambda self : self.do_test_source_offset('NIRCam', theta=45.0)

    test_miri_00 = lambda self : self.do_test_source_offset('MIRI', theta=0.0)
    test_miri_45 = lambda self : self.do_test_source_offset('MIRI', theta=45.0)


class Test_MIRI_FQPM(unittest.TestCase):
    def test_fqpm(self, theta=0.0, clobber=True):
        poppy._FLUXCHECK=True
        miri = webbpsf.MIRI()
        miri.pupilopd = None
        miri.filter='F1065C'
        miri.image_mask = 'FQPM1065'
        miri.pupil_mask = 'MASKFQPM'
        
        nlam = 20
        oversample=2


        for offset in np.linspace(0.0, 1.0, 100):
            miri.options['source_offset_theta'] = 0.0
            miri.options['source_offset_r'] = offset

            if not os.path.exists('test_miri_fqpm_t0_r%.2f.fits' % offset) or clobber:
                psf = miri.calcPSF(oversample=oversample, nlambda=nlam, save_intermediates=False, display=True)#, monochromatic=10.65e-6)
                psf.writeto('test_miri_fqpm_t0_r%.2f.fits' % offset, clobber=clobber)
            if not os.path.exists('test_miri_fqpm_t45_r%.2f.fits' % offset) or clobber:
                miri.options['source_offset_theta'] = 45#np.pi/4
                psf = miri.calcPSF(oversample=oversample, nlambda=nlam, save_intermediates=False, display=True)#, monochromatic=10.65e-6)
                psf.writeto('test_miri_fqpm_t45_r%.2f.fits' % offset, clobber=clobber)
 
class Test_nircam_coron(unittest.TestCase):
    " Test NIRCam coronagraph by computing a whole bunch of models "

    def test_blc_circ(self):
        self.do_test_blc(kind='circular')
    def test_blc_wedge(self):
        self.do_test_blc(kind='linear')


    def do_test_blc(self, clobber=False, kind='circular'):
        poppy._FLUXCHECK=True
        nc = webbpsf.NIRCam()
        nc.pupilopd = None
        nc.filter='F210M'
        if kind =='circular':
            nc.image_mask = 'MASK210R'
            nc.pupil_mask = 'CIRCLYOT'
            fn = 'm210r'
        else:
            nc.image_mask = 'MASKSWB'
            nc.pupil_mask = 'WEDGELYOT'
            fn ='mswb'
 
        nlam = 1 #20
        oversample=2


        #for offset in [0]:
        for offset in np.linspace(0.0, 0.5, 20):
            for angle in [0, 45]:
                nc.options['source_offset_theta'] = angle
                nc.options['source_offset_r'] = offset

                fnout = 'test_nircam_%s_t%d_r%.2f.fits' % (fn, angle, offset)
                if not os.path.exists(fnout) or clobber:
                    psf = nc.calcPSF(oversample=oversample, nlambda=nlam, save_intermediates=False, display=True)#, monochromatic=10.65e-6)
                    psf.writeto(fnout, clobber=clobber)
 
        _log.info("Lots of test files output as test_nircam_*.fits")


def test_run(index=None, wavelength=2e-6):
    """ This function provides a simple interface for running all available tests, or just one """
    #tests = [TestPupils, TestPoppy, Test1, Test2, Test3, Test4, Test5]
    tests = [Test_nircam_coron, Test_MIRI_FQPM, Test_Source_Offset, Test_Image_Size]

    if index is not None:
        if not hasattr(index, '__iter__') : index=[index]
        tests = [tests[i] for i in index]

    suite = unittest.TestSuite()
    for t in tests:
        suite.addTest( unittest.TestLoader().loadTestsFromTestCase( t) )
    unittest.TextTestRunner(verbosity=2).run(suite)


if __name__== "__main__":
    logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')

    unittest.main()
