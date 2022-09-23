import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import copy

import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())

import webbpsf
from .. import webbpsf_core
from .test_errorhandling import _exception_message_starts_with

import pytest


#------------------    NIRCam Tests    ----------------------------
from .test_webbpsf import generic_output_test, do_test_source_offset, do_test_set_position_from_siaf
test_nircam = lambda : generic_output_test('NIRCam')
test_nircam_source_offset_00 = lambda : do_test_source_offset('NIRCam', theta=0.0, monochromatic=2e-6)
test_nircam_source_offset_45 = lambda : do_test_source_offset('NIRCam', theta=45.0, monochromatic=2e-6)

test_nircam_set_siaf = lambda : do_test_set_position_from_siaf('NIRCam',
        ['NRCA5_SUB160', 'NRCA3_DHSPIL_SUB96','NRCA5_MASKLWB_F300M', 'NRCA2_TAMASK210R'])

test_nircam_blc_circ_45 =  lambda : do_test_nircam_blc(kind='circular', angle=45)
test_nircam_blc_circ_0 =   lambda : do_test_nircam_blc(kind='circular', angle=0)


def test_nircam_blc_wedge_0(**kwargs):
    return do_test_nircam_blc(kind='linear', angle=0, **kwargs)

def test_nircam_blc_wedge_45(**kwargs):
    return do_test_nircam_blc(kind='linear', angle=-45, **kwargs)


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


    psf_sam = nc.calc_psf(oversample=oversample, nlambda=1) # should use semi-analytic coronagraph method by default.
    nc.options['no_sam']=True
    psf_fft  = nc.calc_psf(oversample=oversample, nlambda=1)


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
    this routine is just to check for basic functionality of the code and consistency with
    prior results. See the validate_* tests instead for validation against independent
    models of JWST coronagraphy performance - that is NOT what we're trying to do here.

    """

    nc = webbpsf_core.NIRCam()
    nc.pupilopd = None

    nc,ote = webbpsf.enable_adjustable_ote(nc)
    ote._include_nominal_field_dep = False  # disable OTE field dependence model for this test
                                            # for consistency with expected values prepared before that model existed

    nc.filter='F210M'
    offsets = [0, 0.25, 0.50]
    if kind =='circular':
        nc.image_mask = 'MASK210R'
        nc.pupil_mask = 'CIRCLYOT'
        fn = 'm210r'
        expected_total_fluxes=[1.84e-5, 0.0240, 0.1376]  # Based on a prior calculation with WebbPSF
        # values updated slightly for Rev W aperture results
        # Updated 2019-05-02 for coron WFE - changes from [1.35e-5, 0.0240, 0.1376] to [1.84e-5, 0.0240, 0.1376]
    else:
        nc.image_mask = 'MASKSWB'
        nc.pupil_mask = 'WEDGELYOT'
        nc.options['bar_offset'] = 0  # For consistency with how this test was developed
                                      # FIXME update the expected fluxes for the offset positions
                                      # which are now the default.
        fn ='mswb'
        if angle==0:
            expected_total_fluxes=[3.71e-6, .0628, 0.1449]  # Based on a prior calculation with WebbPSF
            # Updated 2019-05-02 for coron WFE - changes from [2.09e-6, .0415, 0.1442] to [3.71e-6, .0628, 0.1449]
        elif angle==45 or angle==-45:
            expected_total_fluxes=[3.71e-6, 0.0221, 0.1192]  # Based on a prior calculation
            # Updated 2016-09-29 for Rev W results - slight change from 0.1171 to 0.1176
            # Updated 2018-02-20 for recoded MASKSWB - changes from 0.0219 to 0.0220; 0.1176 to 0.1192 ??
            # Updated 2019-05-02 for coron WFE - changes from [2.09e-6, 0.0220, 0.1192] to [3.71e-6, 0.0221, 0.1192]
        else:
            raise ValueError("Don't know how to check fluxes for angle={0}".format(angle))

    # If you change either of the following, the expected flux values will need to be updated:
    nlam = 1
    oversample=4

    if outputdir is None:
        import tempfile
        outputdir = tempfile.gettempdir()

    if display:
        nc.display()
        plt.figure()


    #for offset in [0]:
    for offset, exp_flux in zip(offsets, expected_total_fluxes): #np.linspace(0.0, 0.5, nsteps):
        nc.options['source_offset_theta'] = angle
        nc.options['source_offset_r'] = offset

        fnout = os.path.join(outputdir,'test_nircam_%s_t%d_r%.2f.fits' % (fn, angle, offset))

        # We can save the outputs; this is not recommended or useful for general testing but is
        # helpful when/if debugging this test routine itself.
        if not os.path.exists(fnout) or clobber:
            psf = nc.calc_psf(oversample=oversample, nlambda=nlam, save_intermediates=False, display=display)#, monochromatic=10.65e-6)
            if save:
                plt.savefig(fnout+".pdf")
                psf.writeto(fnout, clobber=clobber)
        else:
            psf = fits.open(fnout)
        totflux = psf[0].data.sum()

        #print("Offset: {}    Expected Flux: {}  Calc Flux: {}".format(offset,exp_flux,totflux))

        # FIXME tolerance temporarily increased to 1% in final flux, to allow for using
        # regular propagation rather than semi-analytic. See poppy issue #169
        assert abs(totflux - exp_flux) < 1e-4, f"Total flux {totflux} is out of tolerance relative to expectations {exp_flux}, for offset={offset}, angle={angle}"
        #assert( abs(totflux - exp_flux) < 1e-2 )
        _log.info("File {0} has the expected total flux based on prior reference calculation: {1}".format(fnout, totflux))

    #_log.info("Lots of test files output as test_nircam_*.fits")

def test_nircam_get_detector():
    nc=webbpsf_core.NIRCam()

    detname = nc.detector
    assert detname=='NRCA1'



def test_nircam_auto_pixelscale():
    nc = webbpsf_core.NIRCam()

    nc.filter='F200W'
    assert nc.pixelscale == nc._pixelscale_short
    assert nc.channel == 'short'

    # auto switch to long
    nc.filter='F444W'
    assert nc.pixelscale == nc._pixelscale_long
    assert nc.channel == 'long'

    # and it can switch back to short:
    nc.filter='F200W'
    assert nc.pixelscale == nc._pixelscale_short
    assert nc.channel == 'short'

    nc.pixelscale = 0.0123  # user is allowed to set something custom
    nc.filter='F444W'
    assert nc.pixelscale == 0.0123  # and that persists & overrides the default switching.


    # back to standard scale
    nc.pixelscale = nc._pixelscale_long
    # switch short again
    nc.filter='F212N'
    assert nc.pixelscale == nc._pixelscale_short
    assert nc.channel == 'short'

    # And test we can switch based on detector names too
    nc.detector ='NRCA5'
    assert nc.pixelscale == nc._pixelscale_long
    assert nc.channel == 'long'

    nc.detector ='NRCB1'
    assert nc.pixelscale == nc._pixelscale_short
    assert nc.channel == 'short'

    nc.detector ='NRCA3'
    assert nc.pixelscale == nc._pixelscale_short
    assert nc.channel == 'short'


    nc.auto_channel = False
    # now we can switch filters and nothing else should change:
    nc.filter='F480M'
    assert nc.pixelscale == nc._pixelscale_short
    assert nc.channel == 'short'

    # but changing the detector explicitly always updates pixelscale, regardless
    # of auto_channel being False

    nc.detector = 'NRCA5'
    assert nc.pixelscale == nc._pixelscale_long
    assert nc.channel == 'long'


def test_validate_nircam_wavelengths():
    nc = webbpsf_core.NIRCam()

    # wavelengths fit on shortwave channel -> no exception
    nc.filter='F200W'
    nc._validate_config(wavelengths=np.linspace(nc.SHORT_WAVELENGTH_MIN, nc.SHORT_WAVELENGTH_MAX, 3))
    assert nc.pixelscale == nc._pixelscale_short

    # short wave is selected but user tries a long wave calculation
    with pytest.raises(RuntimeError) as excinfo:
        nc._validate_config(wavelengths=np.linspace(nc.LONG_WAVELENGTH_MIN + 1e-6, nc.LONG_WAVELENGTH_MAX, 3))
    assert _exception_message_starts_with(excinfo,"The requested wavelengths are too long for NIRCam short wave channel")

    # wavelengths fit on long channel -> no exception
    nc.filter='F444W'
    nc._validate_config(wavelengths=np.linspace(nc.LONG_WAVELENGTH_MIN, nc.LONG_WAVELENGTH_MAX, 3))
    assert nc.pixelscale == nc._pixelscale_long

    # long wave is selected but user tries a short  wave calculation
    with pytest.raises(RuntimeError) as excinfo:
        nc._validate_config(wavelengths=np.linspace(nc.SHORT_WAVELENGTH_MIN + 1e-8, nc.SHORT_WAVELENGTH_MAX, 3))
    assert _exception_message_starts_with(excinfo,"The requested wavelengths are too short for NIRCam long wave channel")


    # too short or long for NIRCAM at all
    with pytest.raises(RuntimeError) as excinfo:
        nc._validate_config(wavelengths=np.linspace(nc.SHORT_WAVELENGTH_MIN - 1e-6, nc.SHORT_WAVELENGTH_MIN, 3))
    assert _exception_message_starts_with(excinfo,"The requested wavelengths are too short to be imaged with NIRCam")

    with pytest.raises(RuntimeError) as excinfo:
        nc._validate_config(wavelengths=np.linspace(nc.LONG_WAVELENGTH_MAX, nc.LONG_WAVELENGTH_MAX + 1e-6, 3))
    assert _exception_message_starts_with(excinfo,"The requested wavelengths are too long to be imaged with NIRCam")

def test_nircam_coron_unocculted(plot=False):
    """ NIRCam with lyot mask but not an occulter
    See https://github.com/mperrin/webbpsf/issues/157
    """

    nc = webbpsf_core.NIRCam()
    nc.pupilopd = None
    nc.filter='F212N'
    nc.pupil_mask='WEDGELYOT'
    nc.image_mask=None

    if plot:
        nc.display()

    psf = nc.calc_psf(monochromatic=2.12e-6)
    return(psf)

def test_defocus(fov_arcsec=1, display=False):
    """Test we can apply a defocus to a PSF
    via either a weak lens, or via the options dict,
    and we get consistent results either way.

    Note this is now an *inexact* comparison, because the weak lenses now include non-ideal effects, in particular field dependent astigmatism

    Test for #59 among other things
    """
    nrc = webbpsf_core.NIRCam()
    nrc.set_position_from_aperture_name('NRCA3_FP1')
    nrc.pupilopd=None
    nrc.include_si_wfe=False

    # Calculate defocus with a weak lens
    nrc.pupil_mask = 'WEAK LENS +4'
    psf = nrc.calc_psf(nlambda=1, fov_arcsec=fov_arcsec, oversample=1, display=False, add_distortion=False)

    # Calculate equivalent via the options structure
    nrc.pupil_mask = None
    nrc.options['defocus_waves']=3.9024 # as measured
    nrc.options['defocus_wavelength']=2.12e-6
    psf_2 = nrc.calc_psf(nlambda=1, fov_arcsec=fov_arcsec, oversample=1, display=False, add_distortion=False)

    assert np.allclose(psf[0].data, psf_2[0].data, atol=1e-4), "Defocused PSFs calculated two ways don't agree as precisely as expected"

    if display:
        import webbpsf
        plt.figure()
        webbpsf.display_psf(psf)
        plt.figure()
        webbpsf.display_psf(psf_2)


def test_ways_to_specify_weak_lenses():
    """ There are multiple ways to specify combinations of weak lenses. Test they work as expected.

    """


    testcases = (
        # FILTER  PUPIL   EXPECTED_DEFOCUS
        # Test methods directly specifying a single element
        ('F212N', 'WLM8', 'WLM8'),
        ('F200W', 'WLP8', 'WLP8'),
        ('F187N', 'WLP8', 'WLP8'),
        # Note WLP4 can be specified as filter or pupil element or both
        ('WLP4', 'WLP4', 'WLP4'),
        (None, 'WLP4', 'WLP4'),
        ('WLP4', None, 'WLP4'),
        # Test methods directly specifying a pair of elements stacked together
        ('WLP4', 'WLM8', 'WLM4'),
        ('WLP4', 'WLP8', 'WLP12'),
        # Test methods using virtual pupil elements WLM4 and WLP12
        ('WLP4', 'WLM4', 'WLM4'),
        ('WLP4', 'WLP12', 'WLP12'),
        ('F212N', 'WLM4', 'WLM4'),
        ('F212N', 'WLP12', 'WLP12'),
    )

    nrc = webbpsf_core.NIRCam()
    nrc.pupilopd = None # irrelevant for this test and slows it down
    nrc.include_si_wfe = False # irrelevant for this test and slows it down
    for filt, pup, expected in testcases:
        nrc.pupil_mask = pup
        if filt is not None: nrc.filter = filt

        assert expected in [p.name for p in nrc.get_optical_system().planes], "Optical system did not contain expected plane {} for {}, {}".format(expected, filt, pup)


def test_nircam_coron_wfe_offset(fov_pix=15, oversample=2, fit_gaussian=True):
    """
    Test offset of LW coronagraphic PSF w.r.t. wavelength due to optical wedge dispersion.
    Option to fit a Gaussian to PSF core in order to better determine peak position.
    Difference from 2.5 to 3.3 um should be ~0.015mm.
    Difference from 3.3 to 5.0 um should be ~0.030mm.
    """

    # Disable Gaussian fit if astropy not installed
    if fit_gaussian:
        try:
            from astropy.modeling import models, fitting
        except ImportError:
            fit_gaussian = False

    # Ensure oversample to >1 no Gaussian fitting
    if fit_gaussian == False:
        oversample = 2 if oversample<2 else oversample
        rtol = 0.2
    else:
        rtol = 0.1

    # Set up an off-axis coronagraphic PSF
    inst = webbpsf_core.NIRCam()
    inst.filter = 'F335M'
    inst.pupil_mask = 'CIRCLYOT'
    inst.image_mask = None
    inst.include_si_wfe = True
    inst.options['jitter'] = None

    # size of an oversampled pixel in mm (detector pixels are 18um)
    mm_per_pix = 18e-3/oversample

    # Investigate the differences between three wavelengths
    warr = np.array([2.5,3.3,5.0])

    # Find PSF position for each wavelength
    yloc = []
    for w in warr:
        hdul = inst.calc_psf(monochromatic=w*1e-6, oversample=oversample, add_distortion=False, fov_pixels=fov_pix)

        # Vertical image cross section of oversampled PSF
        im = hdul[0].data
        sh = im.shape
        xvals = mm_per_pix * (np.arange(sh[0]) - sh[0]/2)
        yvals = im[:,int(sh[1]/2)]

        # Fit 1D Gaussian to vertical cross section of PSF
        if fit_gaussian:
            # Create Gaussian model fit of PSF core to determine y offset
            g_init = models.Gaussian1D(amplitude=yvals.max(), mean=0, stddev=0.01)
            fit_g = fitting.LevMarLSQFitter()
            g = fit_g(g_init, xvals, yvals)
            yloc.append(g.mean.value)
        else:
            # Just use PSF max location
            yloc.append(xvals[yvals==yvals.max()][0])
    yloc = np.array(yloc)

    # Difference from 2.5 to 3.3 um should be ~0.015mm
    diff_25_33 = np.abs(yloc[0] - yloc[1])
    assert np.allclose( diff_25_33, 0.016, rtol=rtol), "PSF shift between {:.2f} and {:.2f} um of {:.3f} mm does not match expected value (~0.016 mm).".format(warr[1], warr[0], diff_25_33)
    # Difference from 3.3 to 5.0 um should be ~0.030mm
    diff_50_33 = np.abs(yloc[2] - yloc[1])
    assert np.allclose( diff_50_33, 0.032, rtol=rtol), "PSF shift between {:.2f} and {:.2f} um of {:.3f} mm does not match expected value (~0.032 mm).".format(warr[1], warr[2], diff_50_33)

def test_nircam_auto_aperturename():
    """
    Test that correct apertures are chosen depending on channel, module, detector, mode, etc.
    """
    import pysiaf 

    nc = webbpsf_core.NIRCam()

    nc.filter='F200W'
    assert nc.aperturename == 'NRCA1_FULL'

    # auto switch to long
    nc.filter = 'F444W'
    assert nc.aperturename == 'NRCA5_FULL'

    # and it can switch back to short
    nc.filter='F200W'
    assert nc.aperturename == 'NRCA1_FULL'

    # user is allowed to set custom something
    nc.aperturename = 'NRCA3_SUB400P'
    assert nc.aperturename == 'NRCA3_SUB400P'
    assert nc.detector == 'NRCA3'
    # switching to a different SW filter should not change apname
    nc.filter = 'F212N'
    assert nc.aperturename == 'NRCA3_SUB400P'
    assert nc.detector == 'NRCA3'

    # changing detector will update apname
    nc.detector = 'NRCB5'
    assert nc.aperturename == 'NRCB5_FULL'

    # and switch it back to A1
    nc.detector = 'NRCA1'
    assert nc.aperturename == 'NRCA1_FULL'

    # check Lyot stop modes but no coronagraphic occulter
    nc.pupil_mask = 'CIRCLYOT'
    assert (nc.aperturename == 'NRCA2_FULL_WEDGE_RND') or (nc.aperturename == 'NRCA2_FULL_MASK210R')
    nc.pupil_mask = 'MASKRND'
    assert (nc.aperturename == 'NRCA2_FULL_WEDGE_RND') or (nc.aperturename == 'NRCA2_FULL_MASK210R')
    # if we switch to LW we should get an aperture on A5
    nc.detector='NRCA5'
    assert (nc.aperturename == 'NRCA5_FULL_WEDGE_RND')

    # Add in coronagraphic occulter
    nc.image_mask = 'MASK210R'
    assert nc.aperturename == 'NRCA2_FULL_MASK210R'

    # Change to LW A
    nc.filter = 'F444W'
    nc.image_mask = 'MASK430R'
    assert nc.aperturename == 'NRCA5_FULL_MASK430R'
    nc.image_mask = None
    assert (nc.aperturename == 'NRCA5_FULL_WEDGE_RND') or (nc.aperturename == 'NRCA5_FULL_MASK335R')
    nc.pupil_mask = 'WEDGELYOT'
    assert (nc.aperturename == 'NRCA5_FULL_WEDGE_BAR') or (nc.aperturename == 'NRCA5_FULL_MASKLWB')

    # Switch back to LW imaging
    nc.pupil_mask = None
    assert nc.aperturename == 'NRCA5_FULL'
    # Switch back to SW imaging
    nc.filter='F200W'
    assert nc.aperturename == 'NRCA1_FULL'

    # Turn off auto_aperturename
    nc.auto_aperturename = False
    # now we can switch filters and apname should not change
    nc.filter='F480M'
    assert nc.aperturename == 'NRCA1_FULL'

    # but changing the detector explicitly always updates apname
    nc.detector = 'NRCA5'
    assert nc.aperturename == 'NRCA5_FULL'

    # Test the ability to set arbitrary apertures by name, using a representative subset
    names = ['NRCA1_FULL', 'NRCA3_FP1', 'NRCA2_FP4MIMF', 'NRCB2_FP3MIMF', 'NRCB5_FULL', 'NRCB5_TAMASK430R']

    for apname in names:
        nc.set_position_from_aperture_name(apname)
        # CHeck aperture name and detector name
        assert nc.aperturename == apname, "Aperture name did not match"
        assert nc.detector == apname.split('_')[0], "Detector name did not match"

        # check the detector positions match (to integer pixel precision)
        assert nc.detector_position[0] == int(nc.siaf[
                                                   apname].XSciRef), f"Aperture XSci did not match for {apname}: {nc.detector_position}, {nc.siaf[apname].XSciRef} "
        assert nc.detector_position[1] == int(nc.siaf[
                                                   apname].YSciRef), f"Aperture YSci did not match for {apname}: {nc.detector_position}, {nc.siaf[apname].YSciRef} "

    # Test that switching any detector sets to the FULL aperture for that.
    for det in nc.detector_list:
        nc.detector = det
        assert nc.aperturename == f'{det}_FULL'


def test_coron_shift(offset_npix_x=4, offset_npix_y=-3, plot=False):
    """Test that the two different ways of specifying an offset coronagraphic PSF yield consistent results.
    We can compute an offset PSF either by offsetting the source ('source_offset_{x/y}') or
    by offsetting the coronagraph mask in the opposite direction ('coron_offset_{x/y}').

    Test this by computing an offset PSF in both ways, and checking that the output results are consistent.
    """

    # Make two NIRCam instances, identical instrument config
    nrc1 = webbpsf.NIRCam()
    nrc1.image_mask = 'MASK335R'
    nrc1.pupil_mask = 'MASKRND'
    nrc1.filter='F335M'
    nrc2 = copy.deepcopy(nrc1)

    # Set one up to use source_offset and the other to use coron_shift
    offset_x = nrc1.pixelscale*offset_npix_x
    nrc1.options['source_offset_x'] = offset_x
    nrc2.options['coron_shift_x'] = -offset_x  # Note opposite sign convention!

    offset_y = nrc1.pixelscale*offset_npix_y
    nrc1.options['source_offset_y'] = offset_y
    nrc2.options['coron_shift_y'] = -offset_y  # Note opposite sign convention!


    # Compute PSFs
    psf1, waves1 = nrc1.calc_psf(monochromatic=3.3e-6, fov_pixels=101, add_distortion=False, return_intermediates=True)
    psf2, waves2 = nrc2.calc_psf(monochromatic=3.3e-6, fov_pixels=101, add_distortion=False, return_intermediates=True)

    # Register the one using coron_shift to align with the one that used source_offset
    shifted_psf2 = copy.deepcopy(psf2)
    shifted_psf2[1].data = np.roll(shifted_psf2[1].data, offset_npix_y, axis=0)
    shifted_psf2[1].data = np.roll(shifted_psf2[1].data, offset_npix_x, axis=1)

    # Cut out the central region of the image
    #  (we want to discard the edge regions where the np.roll will wrap around)
    cutout_1 = psf1[1].data[25:75, 25:75]
    cutout_2 = shifted_psf2[1].data[25:75, 25:75]

    if plot:
        fig, axes = plt.subplots(figsize=(16,9), ncols=2)
        webbpsf.display_psf(psf1, ax=axes[0], vmax=1e-4, ext=1, title='Using source_offset', colorbar_orientation='horizontal')
        webbpsf.display_psf(shifted_psf2, ax=axes[1], vmax=1e-4, ext=1, title='Using coron_shift', colorbar_orientation='horizontal')
        for ax in axes:
            ax.axvline(offset_x, ls='--', color='blue')
            ax.axhline(offset_y, ls='--', color='blue')
            ax.axvline(0, ls='-', color='cyan')
            ax.axhline(0, ls='-', color='cyan')

    # Check the two ways of computing this yield consistent results
    assert np.isclose(cutout_1.sum(), cutout_2.sum()), "PSF cutout sums should be consistent"

    assert np.allclose(cutout_1, cutout_2), "PSF cutouts should be consistent"
