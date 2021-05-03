import os
import numpy as np
import pytest
from webbpsf import roman, measure_fwhm
from numpy import allclose


GRISM_FILTER = roman.GRISM_FILTER
PRISM_FILTER = roman.PRISM_FILTER
MASKED_FLAG = "FULL_MASK"
UNMASKED_FLAG = "RIM_MASK"
AUTO_FLAG = "AUTO"

def test_WFI_psf():
    """
    Just test that instantiating WFI works and can compute a PSF without raising
    any exceptions
    """
    wi = roman.WFI()
    wi.calc_psf(fov_pixels=4)


def test_WFI_filters():
    wi = roman.WFI()
    filter_list = wi.filter_list
    for filter in filter_list:
        wi = roman.WFI()
        wi.filter = filter
        wi.calc_psf(fov_pixels=4, oversample=1, nlambda=3)

def test_aberration_detector_position_setter():
    detector = roman.FieldDependentAberration(4096, 4096)

    with pytest.raises(ValueError) as excinfo:
        detector.field_position = (-1, 1)
    assert 'pixel_x' in str(excinfo.value), 'Failed to raise exception for small out-of-bounds ' \
                                            'x pixel position'
    with pytest.raises(ValueError) as excinfo:
        detector.field_position = (4096+1, 1)
    assert 'pixel_x' in str(excinfo.value), 'Failed to raise exception for large out-of-bounds ' \
                                            'x pixel position'
    with pytest.raises(ValueError) as excinfo:
        detector.field_position = (1, -1)
    assert 'pixel_y' in str(excinfo.value), 'Failed to raise exception for small out-of-bounds ' \
                                            'y pixel position'
    with pytest.raises(ValueError) as excinfo:
        detector.field_position = (1, 4096+1)
    assert 'pixel_y' in str(excinfo.value), 'Failed to raise exception for large out-of-bounds ' \
                                            'y pixel position'

    valid_pos = (1.0, 1.0)
    detector.field_position = valid_pos
    assert detector._field_position == valid_pos, 'Setting field position through setter did not ' \
                                                  'update private `_field_position` value'


def test_WFI_fwhm():
    """
    Test that computed PSFs are physically realistic, at least relatively.
    Loose test...
    """
    wfi = roman.WFI()

    wfi.pupilopd = None
    wfi.options['jitter'] = None

    wfi.filter = 'F062'
    fwhm_f062 = measure_fwhm(wfi.calc_psf(oversample= 6))

    wfi.filter = 'F184'
    fwhm_f184 = measure_fwhm(wfi.calc_psf(oversample= 6))

    assert (4.0 > fwhm_f184/fwhm_f062 > 2.0)


def test_WFI_pupil_controller():
    wfi = roman.WFI()

    for detector in wfi.detector_list:
        wfi.detector = detector

        detector_cropped = detector[:3] + str(int((detector[3:])))  # example "SCA01" -> "SCA1"

        unmasked_pupil_path = os.path.join(wfi._pupil_controller._pupil_basepath,
                                           '{}_rim_mask.fits.gz'.format(detector_cropped))

        masked_pupil_path = os.path.join(wfi._pupil_controller._pupil_basepath,
                                         '{}_full_mask.fits.gz'.format(detector_cropped))

        assert os.path.isfile(unmasked_pupil_path), "Pupil file missing {}".format(unmasked_pupil_path)
        assert os.path.isfile(masked_pupil_path), "Pupil file missing {}".format(masked_pupil_path)

        # Test detector change was successful
        assert wfi.detector == detector, "WFI detector was not set correctly"
        assert wfi._unmasked_pupil_path == unmasked_pupil_path, "unmasked_pupil_path was not set correctly"
        assert wfi._masked_pupil_path == masked_pupil_path, "masked_pupil_path was not set correctly"
        assert wfi.pupil in [unmasked_pupil_path, masked_pupil_path], "pupil was not set correctly"

        # Test mask overriding
        wfi.pupil_mask = MASKED_FLAG
        assert wfi.pupil == masked_pupil_path, "pupil was not set correctly"
        assert wfi._pupil_controller.auto_pupil is False, "auto_pupil is active after user override"
        assert wfi._pupil_controller._pupil_mask == wfi.pupil_mask, "pupil mask was not set correctly"

        wfi.pupil_mask = UNMASKED_FLAG
        assert wfi.pupil == unmasked_pupil_path, "pupil was not set correctly"
        assert wfi._pupil_controller.auto_pupil is False, "auto_pupil is active after user override"
        assert wfi._pupil_controller._pupil_mask == wfi.pupil_mask, "pupil mask was not set correctly"

        # Outdated mask overriding backward comparability test:
        wfi.pupil_mask = "COLD_PUPIL"
        assert wfi.pupil == masked_pupil_path, "pupil was not set correctly"
        assert wfi._pupil_controller.auto_pupil is False, "auto_pupil is active after user override"
        assert wfi._pupil_controller._pupil_mask == wfi.pupil_mask, "pupil mask was not set correctly"

        wfi.pupil_mask = "UNMASKED"
        assert wfi.pupil == unmasked_pupil_path, "pupil was not set correctly"
        assert wfi._pupil_controller.auto_pupil is False, "auto_pupil is active after user override"
        assert wfi._pupil_controller._pupil_mask == wfi.pupil_mask, "pupil mask was not set correctly"

        wfi.pupil_mask = AUTO_FLAG
        assert wfi._pupil_controller.auto_pupil is True, "auto_pupil is inactive after mask is set to AUTO"
        assert wfi._pupil_controller._pupil_mask == wfi.pupil_mask, "pupil mask was not set correctly"

        # Test filters
        for filter in wfi.filter_list:
            wfi.filter = filter
            if filter in wfi._pupil_controller._masked_filters:
                assert wfi.pupil == masked_pupil_path, \
                    "Pupil did not set to correct value according to filter {}".format(filter)
            else:
                assert wfi.pupil == unmasked_pupil_path, \
                    "Pupil did not set to correct value according to filter {}".format(filter)

    # Test calculating a single PSF
    wfi = roman.WFI()
    wfi.detector = detector
    valid_pos = (4000, 1000)
    wfi.detector_position = valid_pos
    wfi.pupil_mask = "COLD_PUPIL"
    assert wfi.pupil == masked_pupil_path, "Pupil did not set to correct value according to override"
    wfi.calc_psf(fov_pixels=4)
    assert wfi.pupil == masked_pupil_path, "Pupil did not set to correct value according to override"


def test_WFI_detector_position_setter():
    wfi = roman.WFI()
    wfi.detector = 'SCA01'
    valid_pos = (4000, 1000)
    wfi.detector_position = valid_pos
    assert wfi._detectors[wfi._detector].field_position == valid_pos, (
        "Setting field position through Instrument.detector_position did not update field_position "
        "for the detector's aberration optic"
    )
    assert wfi.detector_position == valid_pos, "`detector_position` getter doesn't reflect " \
                                               "assignment to setter"

def test_WFI_includes_aberrations():
    wfi = roman.WFI()
    wfi.detector = 'SCA01'
    osys = wfi.get_optical_system()
    assert isinstance(osys[2], roman.FieldDependentAberration), (
        "Third plane of Roman WFI optical system should be the "
        "field dependent aberration virtual optic"
    )

def test_WFI_chooses_pupil_masks():
    wfi = roman.WFI()

    def autopupil():
        """Helper to trigger pupil selection in testing"""
        wavelengths, _ = wfi._get_weights()
        wfi._validate_config(wavelengths=wavelengths)
    wfi.filter = 'F087'
    autopupil()
    assert wfi.pupil == wfi._unmasked_pupil_path, "WFI did not select unmasked pupil for F087"
    wfi.filter = 'F184'
    autopupil()
    assert wfi.pupil == wfi._masked_pupil_path, "WFI did not select masked pupil for F158"
    wfi.filter = 'F087'
    autopupil()
    assert wfi.pupil == wfi._unmasked_pupil_path, "WFI did not re-select unmasked pupil for F087"

    def _test_filter_pupil(filter_name, expected_pupil):
        wfi.filter = 'F087'
        autopupil()
        wfi.filter = filter_name
        autopupil()
        assert wfi.pupil == expected_pupil, "Expected pupil {} " \
                                            "for filter {}".format(filter_name, expected_pupil)

    _test_filter_pupil('F106', wfi._unmasked_pupil_path)
    _test_filter_pupil('F129', wfi._unmasked_pupil_path)
    _test_filter_pupil('F062', wfi._unmasked_pupil_path)
    _test_filter_pupil('F158', wfi._unmasked_pupil_path)
    _test_filter_pupil('F146', wfi._unmasked_pupil_path)
    _test_filter_pupil(PRISM_FILTER, wfi._unmasked_pupil_path)

    _test_filter_pupil('F184', wfi._masked_pupil_path)
    _test_filter_pupil(GRISM_FILTER, wfi._masked_pupil_path)

def test_swapping_modes(wfi=None):

    if wfi is None:
        wfi = roman.WFI()

    tests = [
        # [filter, mode, pupil_file]
        ['F062', 'imaging',  wfi._unmasked_pupil_path],
        ['F184', 'imaging', wfi._masked_pupil_path],
        [PRISM_FILTER, 'prism', wfi._unmasked_pupil_path],
        [GRISM_FILTER, 'grism', wfi._masked_pupil_path],
    ]

    for test_filter, test_mode, test_pupil in tests:
        wfi.filter = test_filter
        assert wfi.filter == test_filter
        assert wfi.mode == test_mode
        assert wfi._current_aberrations_file == wfi._aberrations_files[test_mode]
        assert wfi.pupil == test_pupil

def test_custom_aberrations():

    wfi = roman.WFI()

    # Use grism aberrations_file for testing
    test_aberrations_file = wfi._aberrations_files['grism']

    # Test override
    # -------------
    wfi.override_aberrations(test_aberrations_file)

    for filter in wfi.filter_list:
        wfi.filter = filter
        assert wfi._current_aberrations_file == test_aberrations_file, "Filter change caused override to fail"

    # Test Release Override
    # ---------------------
    wfi.reset_override_aberrations()
    assert wfi._aberrations_files['custom'] is None, "Custom aberrations file not deleted on override release."
    test_swapping_modes(wfi)

def test_WFI_limits_interpolation_range():
    wfi = roman.WFI()
    det = wfi._detectors['SCA01']
    det.get_aberration_terms(1.29e-6)
    det.field_position = (0, 0)
    det.get_aberration_terms(1.29e-6)

    with pytest.raises(ValueError) as excinfo:
        det.field_position = (500000, 0)
    assert 'Requested pixel_x position' in str(excinfo.value), (
        "FieldDependentAberration did not error on out-of-bounds field point"
    )

    with pytest.raises(ValueError) as excinfo:
        det.field_position = (-1, 0)
    assert 'Requested pixel_x position' in str(excinfo.value), (
        "FieldDependentAberration did not error on out-of-bounds field point"
    )

    with pytest.raises(ValueError) as excinfo:
        det.field_position = (0, 500000)
    assert 'Requested pixel_y position' in str(excinfo.value), (
        "FieldDependentAberration did not error on out-of-bounds field point"
    )

    with pytest.raises(ValueError) as excinfo:
        det.field_position = (0, -1)
    assert 'Requested pixel_y position' in str(excinfo.value), (
        "FieldDependentAberration did not error on out-of-bounds field point"
    )

    det.field_position = (2048, 2048)

    # Test the get_aberration_terms function uses approximated wavelength when
    # called with an out-of-bound wavelength.
    assert allclose(det.get_aberration_terms(2.0e-6), det.get_aberration_terms(2.5e-6)), (
        "Aberration outside wavelength range did not return closest value."
    )

    assert allclose(det.get_aberration_terms(0.48e-6), det.get_aberration_terms(0.40e-6)), (
        "Aberration outside wavelength range did not return closest value."
    )

    # Test border pixels that are outside of the ref data
    # As of cycle 8 and 9, (4, 4) is the first pixel so we
    # check if (0, 0) is approximated to (4, 4) via nearest point
    # approximation:

    det.field_position = (0, 0)
    coefficients_outlier = det.get_aberration_terms(1e-6)

    det.field_position = (4, 4)
    coefficients_data = det.get_aberration_terms(1e-6)

    assert np.allclose(coefficients_outlier, coefficients_data), "nearest point extrapolation " \
                                                                 "failed for outlier field point"

def test_CGI_detector_position():
    """ Test existence of the CGI detector position etc, and that you can't set it."""
    cgi = roman.CGI()

    valid_pos = (512,512)
    assert cgi.detector_position == valid_pos, "CGI detector position isn't as expected"

    with pytest.raises(RuntimeError) as excinfo:
        cgi.detector_position = valid_pos
    assert 'not adjustable' in str(excinfo.value), ("Failed to raise exception for"\
                                                        "trying to change CGI detector position.")

def test_CGI_psf(display=False):
    """
    Just test that instantiating CGI works and can compute a PSF without raising
    any exceptions
    """
    char_spc = roman.CGI()
    char_spc.mode = 'CHARSPC_F660'

    #print('Reading instrument data from {:s}'.format(charspc._WebbPSF_basepath)
    #print('Filter list: {:}'.format(charspc.filter_list))

    monopsf = char_spc.calc_psf(nlambda=1, display=False)
    if display:
        roman.poppy.display_psf(monopsf)
