from __future__ import division, print_function, absolute_import, unicode_literals
import pytest
from webbpsf import wfirst

def test_WFI_psf():
    """
    Just test that instantiating WFI works and can compute a PSF without raising
    any exceptions
    """
    wi = wfirst.WFI()
    wi.calc_psf(fov_pixels=4)

def test_detector_position_setter():
    detector = wfirst.FieldDependentAberration(4096, 4096)

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

def test_WFI_detector_position_setter():
    wfi = wfirst.WFI()
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
    wfi = wfirst.WFI()
    wfi.detector = 'SCA01'
    osys = wfi._getOpticalSystem()
    assert isinstance(osys[2], wfirst.FieldDependentAberration), (
        "Third plane of WFIRST WFI optical system should be the "
        "field dependent aberration virtual optic"
    )

def test_WFI_chooses_pupil_masks():
    wfi = wfirst.WFI()

    def autopupil():
        """Helper to trigger pupil selection in testing"""
        wavelengths, _ = wfi._getWeights()
        wfi._validateConfig(wavelengths=wavelengths)
    wfi.filter = 'Z087'
    autopupil()
    assert wfi.pupil == wfi._unmasked_pupil_path, "WFI did not select unmasked pupil for Z087"
    wfi.filter = 'H158'
    autopupil()
    assert wfi.pupil == wfi._masked_pupil_path, "WFI did not select masked pupil for H158"
    wfi.filter = 'Z087'
    autopupil()
    assert wfi.pupil == wfi._unmasked_pupil_path, "WFI did not re-select unmasked pupil for Z087"

    def _test_filter_pupil(filter_name, expected_pupil):
        wfi.filter = 'Z087'
        autopupil()
        wfi.filter = filter_name
        autopupil()
        assert wfi.pupil == expected_pupil, "Expected pupil {} " \
                                            "for filter {}".format(filter_name, expected_pupil)

    _test_filter_pupil('Y106', wfi._unmasked_pupil_path)
    _test_filter_pupil('J129', wfi._unmasked_pupil_path)
    _test_filter_pupil('H158', wfi._masked_pupil_path)
    _test_filter_pupil('F184', wfi._masked_pupil_path)
    _test_filter_pupil('W149', wfi._masked_pupil_path)

def test_WFI_limits_interpolation_range():
    wfi = wfirst.WFI()
    det = wfi._detectors['SCA01']
    det.get_aberration_terms(1.29e-6)
    det.field_position = (0, 0)
    with pytest.raises(RuntimeError) as excinfo:
        det.get_aberration_terms(1.29e-6)
    assert 'out-of-bounds field point' in str(excinfo.value), (
        "FieldDependentAberration did not error on out-of-bounds field point"
    )
    with pytest.raises(RuntimeError) as excinfo:
        det.get_aberration_terms(1.29e-6)
    assert 'out-of-bounds field point' in str(excinfo.value), (
        "FieldDependentAberration did not error on out-of-bounds field point"
    )
    det.field_position = (2048, 2048)
    with pytest.raises(RuntimeError) as excinfo:
        det.get_aberration_terms(5e-6)
    assert 'wavelength outside the range' in str(excinfo.value), (
        "FieldDependentAberration did not error on out-of-bounds wavelength"
    )
