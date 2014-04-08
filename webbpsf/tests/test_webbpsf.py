import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits


import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())


from .. import webbpsf_core
import poppy

#poppy._log.setLevel(logging.INFO)


def generic_output_test(iname):

    _log.info("Testing image output sizes for %s " % iname)
    inst = webbpsf_core.Instrument(iname)
    pxscale = inst.pixelscale
    fov_arcsec = 5.0

    PSF = inst.calcPSF(nlambda=1, fov_pixels = 100, oversample=1)
    assert(PSF[0].data.shape[0] == 100)

    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    fov_pix = int(np.round(fov_arcsec / pxscale))
    assert(PSF[0].data.shape[0] == fov_pix)

    inst.options['parity'] = 'odd'
    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    assert( np.remainder(PSF[0].data.shape[0],2) == 1)

    inst.options['parity'] = 'even'
    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    assert( np.remainder(PSF[0].data.shape[0],2) == 0)

    # odd array, even oversampling = even
    inst.options['parity'] = 'odd'
    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=2)
    assert( np.remainder(PSF[0].data.shape[0],2) == 0)

    # odd array, odd oversampling = odd
    inst.options['parity'] = 'odd'
    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=3)
    assert( np.remainder(PSF[0].data.shape[0],2) == 1)

test_nircam = lambda : generic_output_test('NIRCam')
test_miri= lambda : generic_output_test('MIRI')
test_nirspec= lambda : generic_output_test('NIRSpec')
test_niriss= lambda : generic_output_test('NIRISS')
test_fgs= lambda : generic_output_test('FGS')


def do_test_source_offset(iname, distance=0.5,  nsteps=1, theta=0.0, display=False):
    """ Test source offsets  
    """

    nc = webbpsf_core.Instrument(iname)
    nc.pupilopd=None

    oversample = 2

    # Calculations
    shift_req = []
    psfs = []

    # unshifted PSF
    #psfs.append(  nc.calcPSF(nlambda=1, oversample=oversample) )
    #shift_req.append(0)

    steps = np.linspace(0, distance, nsteps+1)
    for i, value in enumerate(steps):
        nc.options['source_offset_r'] =  steps[i]
        nc.options['source_offset_theta'] = theta
        #nc.options['source_offset_r'] = i*nc.pixelscale*5
        shift_req.append(nc.options['source_offset_r'])
        psfs.append(  nc.calcPSF(nlambda=1, oversample=oversample) )


    # Control case: an unshifted image
    cent0 = np.asarray(poppy.measure_centroid(psfs[0]))
    center_pix = (psfs[0][0].data.shape[0]-1)/2.0
    assert( abs(cent0[0] == center_pix) < 1e-3 )
    assert( abs(cent0[1] == center_pix) < 1e-3 )
    _log.info("Center of unshifted image: (%d, %d)" % tuple(cent0))

    if display:
        poppy.display_PSF(psfs[0])

    # Compare to shifted case(s)
    for i in range(1, nsteps+1):

        if display:
            poppy.display_PSF(psfs[i])

        cent = poppy.measure_centroid(psfs[i])
        rx = shift_req[i] * (-np.sin(theta*np.pi/180))
        ry = shift_req[i] * (np.cos(theta*np.pi/180))
        _log.info("   Shift_requested:\t(%10.3f, %10.3f)" % (rx, ry))
        shift = (cent-cent0) * (nc.pixelscale/oversample)
        _log.info("   Shift_achieved: \t(%10.3f, %10.3f)" % (shift[1], shift[0]))
        assert( abs(rx -  shift[1]) < 1e-3 )
        assert( abs(ry -  shift[0]) < 1e-3 )

test_nircam_00 = lambda : do_test_source_offset('NIRCam', theta=0.0)
test_nircam_45 = lambda : do_test_source_offset('NIRCam', theta=45.0)
test_miri_00 = lambda : do_test_source_offset('MIRI', theta=0.0)
test_miri_45 = lambda : do_test_source_offset('MIRI', theta=45.0)
test_fgs_00 = lambda : do_test_source_offset('FGS', theta=0.0)
test_fgs_45 = lambda : do_test_source_offset('FGS', theta=45.0)


#------------------    MIRI Tests    ----------------------------

def test_miri_fqpm(theta=0.0, nsteps=3, nlambda=1, clobber=True):
    #poppy._FLUXCHECK=True
    miri = webbpsf_core.MIRI()
    miri.pupilopd = None
    miri.filter='F1065C'
    miri.image_mask = 'FQPM1065'
    miri.pupil_mask = 'MASKFQPM'
    
    oversample=2


    for offset in np.linspace(0.0, 1.0, nsteps):
        miri.options['source_offset_theta'] = 0.0
        miri.options['source_offset_r'] = offset

        if not os.path.exists('test_miri_fqpm_t0_r%.2f.fits' % offset) or clobber:
            psf = miri.calcPSF(oversample=oversample, nlambda=nlambda, save_intermediates=False, display=True)#, monochromatic=10.65e-6)
            psf.writeto('test_miri_fqpm_t0_r%.2f.fits' % offset, clobber=clobber)
        if not os.path.exists('test_miri_fqpm_t45_r%.2f.fits' % offset) or clobber:
            miri.options['source_offset_theta'] = 45#np.pi/4
            psf = miri.calcPSF(oversample=oversample, nlambda=nlambda, save_intermediates=False, display=True)#, monochromatic=10.65e-6)
            psf.writeto('test_miri_fqpm_t45_r%.2f.fits' % offset, clobber=clobber)

    #FIXME - add some assertion tests here. 

#------------------    NIRCam Tests    ----------------------------


test_nircam_blc_circ_45 =  lambda : do_test_nircam_blc(kind='circular', angle=45)
test_nircam_blc_circ_0 =   lambda : do_test_nircam_blc(kind='circular', angle=0)
test_nircam_blc_wedge_0 =  lambda : do_test_nircam_blc(kind='linear', angle=0)
test_nircam_blc_wedge_45 = lambda : do_test_nircam_blc(kind='linear', angle=45)


def do_test_nircam_blc(clobber=False, kind='circular', angle=0, save=False, display=False):
    """ Test NIRCam BLC coronagraphs

    Compute BLC PSFs on axis and offset and check the values against the expectation.
    Note that the 'correct' values are just prior calculations with WebbPSF; the purpose of
    this routine is just to check for basic functionaltiy of the code and consistency with
    prior results. See the validate_* tests instead for validation against independent
    models of JWST coronagraphy performance - that is NOT what we're trying to do here.
    """

    nc = webbpsf_core.NIRCam()
    nc.pupilopd = None
    nc.filter='F210M'
    offsets = [0, 0.25, 0.50]
    if kind =='circular':
        nc.image_mask = 'MASK210R'
        nc.pupil_mask = 'CIRCLYOT'
        fn = 'm210r'
        expected_total_fluxes=[1.35e-5, 0.0237, 0.1367]  # Based on a prior calculation with WebbPSF
    else:
        nc.image_mask = 'MASKSWB'
        nc.pupil_mask = 'WEDGELYOT'
        fn ='mswb'
        if angle==0:
            expected_total_fluxes=[0.0012, 0.0606, 0.1396]  # Based on a prior calculation with WebbPEF
        else:
            expected_total_fluxes=[0.0012, 0.0219, 0.1146]  # Based on a prior calculation 

    nlam = 3 #20
    oversample=2



    #for offset in [0]:
    for offset, exp_flux in zip(offsets, expected_total_fluxes): #np.linspace(0.0, 0.5, nsteps):
        nc.options['source_offset_theta'] = angle
        nc.options['source_offset_r'] = offset

        fnout = 'test_nircam_%s_t%d_r%.2f.fits' % (fn, angle, offset)

        # We can save the outputs; this is not recommended or useful for general testing but is
        # helpful when/if debugging this test routine itself.
        if not os.path.exists(fnout) or clobber:
            psf = nc.calcPSF(oversample=oversample, nlambda=nlam, save_intermediates=False, display=display)#, monochromatic=10.65e-6)
            if save:
                psf.writeto(fnout, clobber=clobber)
        else:
            psf = fits.open(fnout)
        totflux = psf[0].data.sum()

        assert( abs(totflux - exp_flux) < 1e-4 )
        _log.info("File {0} has the expected total flux based on prior reference calculation: {1}".format(fnout, totflux))

    #_log.info("Lots of test files output as test_nircam_*.fits")
#--------------------------------------------------------------------------------

