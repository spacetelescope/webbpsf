import logging
import os

import astropy.units as u
import numpy as np
import pysiaf

import webbpsf
from .. import webbpsf_core
from .test_webbpsf import do_test_set_position_from_siaf, do_test_source_offset, generic_output_test

_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())


# ------------------    MIRI Tests    ----------------------------


def test_miri():
    return generic_output_test('MIRI')


def test_miri_source_offset_00():
    return do_test_source_offset('MIRI', theta=0.0, monochromatic=8e-6)


def test_miri_source_offset_45():
    return do_test_source_offset('MIRI', theta=45.0, monochromatic=8e-6)


def test_miri_set_siaf():
    return do_test_set_position_from_siaf(
        'MIRI',
        [
            'MIRIM_SUB128',
            'MIRIM_FP1MIMF',
            'MIRIM_BRIGHTSKY',
            'MIRIM_TASLITLESSPRISM',
        ],
    )


def do_test_miri_fqpm(
    nlambda=1, clobber=True, angle=0.0, offset=0.0, oversample=2, outputdir=None, display=False, save=False
):
    miri = webbpsf_core.MIRI()
    miri.pupilopd = None
    miri.filter = 'F1065C'
    miri.image_mask = 'FQPM1065'
    miri.pupil_mask = 'MASKFQPM'

    # for offset in np.linspace(0.0, 1.0, nsteps):
    # miri.options['source_offset_theta'] = 0.0
    miri.options['source_offset_r'] = offset

    # for angle in [0,45]:
    miri.options['source_offset_theta'] = angle
    psf = miri.calc_psf(oversample=oversample, nlambda=nlambda, save_intermediates=False, display=display)

    if save:
        if outputdir is None:
            import tempfile

            outputdir = tempfile.gettempdir()

        fn = os.path.join(outputdir, 'test_miri_fqpm_t{0}_r{1:.2f}.fits'.format(angle, offset))
        psf.writeto(fn, clobber=clobber)

    # FIXME - add some assertion tests here.


def test_miri_fqpm_centered(*args, **kwargs):
    do_test_miri_fqpm(angle=0.0, offset=0.0)


def test_miri_fqpm_offset_00(*args, **kwargs):
    do_test_miri_fqpm(angle=0.0, offset=1.0)


def test_miri_fqpm_offset_45(*args, **kwargs):
    do_test_miri_fqpm(angle=45.0, offset=1.0)


def test_miri_aperturename():
    """Test aperture name functionality"""
    miri = webbpsf_core.MIRI()
    assert miri.aperturename == miri._detectors[miri.detector], 'Default SIAF aperture is not as expected'

    ref_tel_coords = miri._tel_coords()

    miri.aperturename = 'MIRIM_SUB256'
    assert miri.detector_position == (128, 128), (
        "Changing to a subarray aperture didn't change the " 'reference pixel coords as expected'
    )
    assert np.any(
        miri._tel_coords() != ref_tel_coords
    ), "Changing to a subarray aperture didn't change the V2V3 coords as expected."


def test_miri_slit_apertures():
    """Test that we can use slit and aperture names that don't map to a specific detector
    Verify that the V2 and V3 coordinates are reported as expected.
    """
    miri = webbpsf_core.MIRI()

    apname = 'MIRIM_SLIT'  # this is the only slit aperture on the MIRI imager
    miri.set_position_from_aperture_name(apname)

    ap = pysiaf.Siaf('MIRI')[apname]

    assert np.isclose(miri._tel_coords()[0].to_value(u.arcsec), ap.V2Ref)
    assert np.isclose(miri._tel_coords()[1].to_value(u.arcsec), ap.V3Ref)

    # Test we can switch back from SLIT to regular type apertures without any error
    miri.set_position_from_aperture_name('MIRIM_SLIT')
    miri.set_position_from_aperture_name('MIRIM_FULL')


def test_lrs_psf():
    """Test that we can at least run an LRS mode PSF calculation

    This is not yet a very good test; does not test the correctness of the output at all.
    """
    miri = webbpsf_core.MIRI()
    miri.set_position_from_aperture_name('MIRIM_SLIT')
    miri.image_mask = 'LRS slit'
    miri.pupil_mask = 'P750L'
    psf = miri.calc_psf(nlambda=1)

    assert psf[0].header['APERNAME'] == miri.aperturename == 'MIRIM_SLIT'


def test_miri_nonsquare_detector():
    """Test that we can handle the slightly different
    dimenssions in X and Y of the MIRI detector"""
    miri = webbpsf_core.MIRI()
    miri.detector_position = (1023, 1031)  # recall this is X, Y order
    assert miri.detector_position == (1023, 1031)


def test_mode_switch():
    """Test switching between imaging and IFU modes, and switching IFU bands
    Also checks this works to switch aperturenane, and conversely setting aperturename switches mode if needed.
    Also checks that this automatically changes the rotation and pixelscale properties, as expected.
    """
    miri = webbpsf_core.MIRI()
    imager_rotation = miri._rotation
    imager_pixelscale = miri.pixelscale

    # Explicitly switch mode to IFU
    miri.mode = 'IFU'
    assert 'IFU' in miri.aperturename
    assert miri.detector == 'MIRIFUSHORT'
    assert miri.aperturename.startswith('MIRIFU_CH')
    assert miri._rotation != imager_rotation
    assert miri.pixelscale > imager_pixelscale
    # Explicitly switch back to imaging
    miri.mode = 'imaging'
    assert 'IFU' not in miri.aperturename
    assert miri.detector == 'MIRIM'
    assert miri.aperturename.startswith('MIRIM_')
    assert miri._rotation == imager_rotation
    assert miri.pixelscale == imager_pixelscale

    # Implicitly switch to IFU
    miri.set_position_from_aperture_name('MIRIFU_CHANNEL3B')
    assert 'IFU' in miri.aperturename
    assert miri.detector == 'MIRIFULONG'
    assert miri.aperturename == 'MIRIFU_CHANNEL3B'
    assert miri._rotation != imager_rotation
    assert miri.pixelscale > imager_pixelscale

    # implicitly switch to imaging
    # LRS is an odd case, SLIT aper type but operates like in imaging mode
    miri.set_position_from_aperture_name('MIRIM_SLIT')
    assert 'IFU' not in miri.aperturename
    assert miri.detector == 'MIRIM'
    assert miri.aperturename.startswith('MIRIM_')
    assert miri._rotation == imager_rotation
    assert miri.pixelscale == imager_pixelscale

    # And back to IFU again:
    miri.mode = 'IFU'
    assert 'IFU' in miri.aperturename
    assert miri.detector == 'MIRIFUSHORT'
    assert miri.aperturename.startswith('MIRIFU_CH')
    assert miri._rotation != imager_rotation
    assert miri.pixelscale > imager_pixelscale

    # band switching should toggle detector and aper name
    miri.band = '4C'
    assert miri.detector == 'MIRIFULONG'
    assert miri.aperturename == 'MIRIFU_CHANNEL4C'
    assert miri.pixelscale > 3 * imager_pixelscale

    miri.band = '2A'
    assert miri.detector == 'MIRIFUSHORT'
    assert miri.aperturename == 'MIRIFU_CHANNEL2A'
    assert imager_pixelscale < miri.pixelscale < 2 * imager_pixelscale

    # Test also we can switch to LRS, and then back to a regular imaging aperture
    # this tests another edge case for aperture and pixel scale switches
    miri.set_position_from_aperture_name('MIRIM_SLIT')
    assert 'IFU' not in miri.aperturename
    assert miri.detector == 'MIRIM'
    assert miri.aperturename == 'MIRIM_SLIT'
    assert miri.pixelscale == imager_pixelscale
    miri.set_position_from_aperture_name('MIRIM_FULL')
    assert 'IFU' not in miri.aperturename
    assert miri.detector == 'MIRIM'
    assert miri.aperturename == 'MIRIM_FULL'
    assert miri.pixelscale == imager_pixelscale


def test_IFU_wavelengths():
    """Test computing the wqvelength sampling for a sim IFU cube"""
    miri = webbpsf_core.MIRI()
    # check mode swith to IFU
    miri.mode = 'IFU'
    miri.band = '2A'
    waves = miri.get_IFU_wavelengths()
    assert isinstance(waves, u.Quantity)

    assert len(waves) > 900  # there are lots of wavelengths in MRScubes
    # and test we can specify a reduced wavelength sampling:
    for n in (10, 100):
        assert len(miri.get_IFU_wavelengths(n)) == n


def test_miri_ifu_broadening():
    """ Basic functional test for the code that adjusts PSF outputs to better match empirical IFU PSFs
    """

    miri = webbpsf_core.MIRI()
    miri.mode = 'IFU'
    psf = miri.calc_psf(monochromatic=6.8e-6, fov_pixels=20)

    fwhm_oversamp = webbpsf.measure_fwhm(psf, ext='OVERSAMP')
    fwhm_overdist = webbpsf.measure_fwhm(psf, ext='OVERDIST')
    assert fwhm_overdist > fwhm_oversamp, "IFU broadening model should increase the FWHM for the distorted extensions"

    fwhm_detsamp = webbpsf.measure_fwhm(psf, ext='DET_SAMP')
    fwhm_detdist = webbpsf.measure_fwhm(psf, ext='DET_DIST')
    assert fwhm_overdist > fwhm_oversamp, "IFU broadening model should increase the FWHM for the distorted extensions"

    # test flux conservation. THis is close but not exact, due to the optical distortion part which is distinct from but
    # happens at same point in the calculation as the IFU broadening effect.
    assert np.isclose(psf['DET_SAMP'].data.sum(), psf['DET_DIST'].data.sum()), "IFU broadening should not change total flux much"


    # Now test that we can also optionally turn off that effect
    miri.options['ifu_broadening'] = None
    psf_nb = miri.calc_psf(monochromatic=6.8e-6, fov_pixels=20)

    fwhm_oversamp = webbpsf.measure_fwhm(psf_nb, ext='OVERSAMP')
    fwhm_overdist = webbpsf.measure_fwhm(psf_nb, ext='OVERDIST')
    assert fwhm_overdist == fwhm_oversamp, "IFU broadening model should be disabled for this test case"

    # test flux conservation with and without the IFU broadening. This should be a more exact match
    assert np.isclose(psf_nb['DET_DIST'].data.sum(), psf['DET_DIST'].data.sum() ), "IFU broadening should not change total flux much"


