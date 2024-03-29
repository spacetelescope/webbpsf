import logging
import os
import astropy.units as u
import pysiaf
import numpy as np

_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())

from .. import webbpsf_core

# ------------------    MIRI Tests    ----------------------------

from .test_webbpsf import generic_output_test, do_test_source_offset, do_test_set_position_from_siaf

test_miri = lambda: generic_output_test('MIRI')
test_miri_source_offset_00 = lambda: do_test_source_offset('MIRI', theta=0.0, monochromatic=8e-6)
test_miri_source_offset_45 = lambda: do_test_source_offset('MIRI', theta=45.0, monochromatic=8e-6)

test_miri_set_siaf = lambda: do_test_set_position_from_siaf('MIRI',
                                                            ['MIRIM_SUB128', 'MIRIM_FP1MIMF', 'MIRIM_BRIGHTSKY',
                                                             'MIRIM_TASLITLESSPRISM', ])


def do_test_miri_fqpm(nlambda=1, clobber=True, angle=0.0, offset=0.0, oversample=2, outputdir=None, display=False,
                      save=False):
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
    """ Test aperture name functionality """
    miri = webbpsf_core.MIRI()
    assert miri.aperturename == miri._detectors[miri.detector], "Default SIAF aperture is not as expected"

    ref_tel_coords = miri._tel_coords()

    miri.aperturename = 'MIRIM_SUB256'
    assert miri.detector_position == (128, 128), "Changing to a subarray aperture didn't change the " \
                                                 "reference pixel coords as expected"
    assert np.any( miri._tel_coords() != ref_tel_coords), "Changing to a subarray aperture didn't change the V2V3 coords as expected."


def test_miri_slit_apertures():
    """Test that we can use slit and aperture names that don't map to a specific detector
    Verify that the V2 and V3 coordinates are reported as expected.
    """
    miri = webbpsf_core.MIRI()

    apname = "MIRIM_SLIT"  # this is the only slit aperture on the MIRI imager
    miri.set_position_from_aperture_name(apname)

    ap = pysiaf.Siaf('MIRI')[apname]

    assert np.isclose(miri._tel_coords()[0].to_value(u.arcsec), ap.V2Ref)
    assert np.isclose(miri._tel_coords()[1].to_value(u.arcsec), ap.V3Ref)

def test_miri_nonsquare_detector():
    """ Test that we can handle the slightly different
    dimenssions in X and Y of the MIRI detector"""
    miri = webbpsf_core.MIRI()
    miri.detector_position = (1023, 1031)  # recall this is X, Y order
    assert miri.detector_position == (1023, 1031)
