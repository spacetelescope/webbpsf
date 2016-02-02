import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())

from .. import webbpsf_core
import poppy
from .test_errorhandling import _exception_message_starts_with

import pytest


#------------------    NIRCam Tests    ----------------------------
from .test_webbpsf import generic_output_test, do_test_source_offset
test_nircam = lambda : generic_output_test('NIRCam')
test_nircam_source_offset_00 = lambda : do_test_source_offset('NIRCam', theta=0.0, monochromatic=2e-6)
test_nircam_source_offset_45 = lambda : do_test_source_offset('NIRCam', theta=45.0, monochromatic=2e-6)

test_nircam_blc_circ_45 =  lambda : do_test_nircam_blc(kind='circular', angle=45)
test_nircam_blc_circ_0 =   lambda : do_test_nircam_blc(kind='circular', angle=0)

@pytest.mark.xfail
def test_nircam_blc_wedge_0():
    return do_test_nircam_blc(kind='linear', angle=0)

@pytest.mark.xfail
def test_nircam_blc_wedge_45():
    return do_test_nircam_blc(kind='linear', angle=45)


# The test setup for this one is not quite right yet
#  See https://github.com/mperrin/webbpsf/issues/30
#  and https://github.com/mperrin/poppy/issues/29
@pytest.mark.xfail
def test_nircam_SAMC(oversample=4):

    _log.info("Comparing semi-analytic and direct FFT calculations for NIRCam coronagraphy")
    nc = webbpsf_core.NIRCam()

    nc.pupilopd=None
    nc.filter='F212N'
    nc.image_mask='MASK210R'
    nc.pupil_mask='CIRCLYOT'
    nc.options['output_mode'] = 'Detector'


    psf_sam = nc.calcPSF(oversample=oversample, nlambda=1) # should use semi-analytic coronagraph method by default.
    nc.options['no_sam']=True
    psf_fft  = nc.calcPSF(oversample=oversample, nlambda=1)


    maxdiff = np.abs(psf_fft[0].data - psf_sam[0].data).max()
    _log.info("Max difference between results: {0} cts/pixel".format( maxdiff))
    assert( maxdiff < 1e-7)


    # The pixel by pixel difference should be small:
    maxdiff = np.abs(psf_fft[0].data - psf_sam[0].data).max()
    #print "Max difference between results: ", maxdiff
    assert( maxdiff < 1e-7)


    # and the overall flux difference should be small also:
    if oversample<=4:
        thresh = 1e-4
    elif oversample==6:
        thresh=5e-5
    elif oversample>=8:
        thresh = 4e-6
    else:
        raise NotImplementedError("Don't know what threshold to use for oversample="+str(oversample))

    # What is a reasonable degree of agreement between the SAM and classical FFT based approaches?
    # Here are the results for total intensity that I get in the two cases, for varying degrees of
    # oversampling:

    # Oversampling          sum(SAM)        sum(FFT)    Ratio
    #    2                  3.213e-4        4.383e-3    13.6
    #    4                  1.494e-5        3.693e-4    24.7
    #    6                  2.856e-5        1.117e-4     3.92
    #    8                  1.494e-5        5.741e-5     3.84
    #   16                  1.655e-5        2.266e-5     1.37
    #   20                  1.596e-5        1.965e-5     1.23


    assert False

    assert np.abs(psf_sam[0].data.sum() - psf_fft[0].data.sum()) < thresh




def do_test_nircam_blc(clobber=False, kind='circular', angle=0, save=False, display=False, outputdir=None):
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
            expected_total_fluxes=[2.09e-6, .0415, 0.1442]  # Based on a prior calculation with WebbPSF
            #expected_total_fluxes=[0.0012, 0.0606, 0.1396]  # Based on a prior calculation with WebbPSF
                                 #2e-6
        elif angle==45:
            expected_total_fluxes=[2.09e-6, 0.0219, 0.1173]  # Based on a prior calculation
            #expected_total_fluxes=[0.0012, 0.0219, 0.1146]  # Based on a prior calculation
        else:
            raise ValueError("Don't know how to check fluxes for angle={0}".format(angle))

    # If you change either of the following, the expected flux values will need to be updated:
    nlam = 1
    oversample=4

    if outputdir is None:
        import tempfile
        outputdir = tempfile.gettempdir()



    #for offset in [0]:
    for offset, exp_flux in zip(offsets, expected_total_fluxes): #np.linspace(0.0, 0.5, nsteps):
        nc.options['source_offset_theta'] = angle
        nc.options['source_offset_r'] = offset

        fnout = os.path.join(outputdir,'test_nircam_%s_t%d_r%.2f.fits' % (fn, angle, offset))

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

@pytest.mark.xfail
def test_nircam_get_detector():
    nc=webbpsf_core.NIRCam()

    detname = nc.detector
    assert detname=='A1'



def test_nircam_auto_pixelscale():
    nc = webbpsf_core.NIRCam()

    nc.filter='F200W'
    wavelengths, _ = nc._getWeights()
    nc._validateConfig(wavelengths=wavelengths)
    assert nc.pixelscale == nc._pixelscale_short

    # auto switch to long
    nc.filter='F444W'
    wavelengths, _ = nc._getWeights()
    nc._validateConfig(wavelengths=wavelengths)
    assert nc.pixelscale == nc._pixelscale_long

    # and it can switch back to short:
    nc.filter='F200W'
    wavelengths, _ = nc._getWeights()
    nc._validateConfig(wavelengths=wavelengths)
    assert nc.pixelscale == nc._pixelscale_short

    nc.pixelscale = 0.0123  # user is allowed to set something custom
    nc.filter='F444W'
    wavelengths, _ = nc._getWeights()
    nc._validateConfig(wavelengths=wavelengths)
    assert nc.pixelscale == 0.0123  # and that persists & overrides the default switching.

    # restore a standard pixel scale (long since we're testing short next)
    nc.pixelscale = nc._pixelscale_long

    # wavelengths fit on shortwave channel -> short channel pixel scale
    nc._validateConfig(wavelengths=np.linspace(nc.SHORT_WAVELENGTH_MIN, nc.SHORT_WAVELENGTH_MAX, 3))
    assert nc.pixelscale == nc._pixelscale_short

    # wavelengths fit on long channel -> long pixel scale
    nc._validateConfig(wavelengths=np.linspace(nc.LONG_WAVELENGTH_MIN, nc.LONG_WAVELENGTH_MAX, 3))
    assert nc.pixelscale == nc._pixelscale_long

    # wavelengths don't fit on one channel -> exception
    with pytest.raises(RuntimeError) as excinfo:
        nc._validateConfig(wavelengths=np.linspace(nc.SHORT_WAVELENGTH_MIN - 1e-6, nc.SHORT_WAVELENGTH_MIN, 3))
    assert _exception_message_starts_with(excinfo,"The requested wavelengths are too short to be imaged with NIRCam")

    with pytest.raises(RuntimeError) as excinfo:
        nc._validateConfig(wavelengths=np.linspace(nc.LONG_WAVELENGTH_MAX, nc.LONG_WAVELENGTH_MAX + 1e-6, 3))
    assert _exception_message_starts_with(excinfo,"The requested wavelengths are too long to be imaged with NIRCam")

    with pytest.raises(RuntimeError) as excinfo:
        nc._validateConfig(wavelengths=np.linspace(nc.SHORT_WAVELENGTH_MAX - 1e-6, nc.LONG_WAVELENGTH_MIN + 1e-6, 3))
    assert _exception_message_starts_with(excinfo,"Wavelengths requested don't fit entirely on either NIRCam short"
                                     " or long wavelength channels")

