import os
import numpy as np
import pytest
from webbpsf import roman, measure_fwhm
from astropy.table import Table
from numpy import allclose


GRISM_FILTERS = roman.GRISM_FILTERS
PRISM_FILTERS = roman.PRISM_FILTERS

def detector_substr(detector):
    """
    change detector string to match file format
    (e.g., "SCA01" -> "SCA_1")
    """
    return f"{detector[:3]}_{str(int((detector[3:])))}"

def pupil_path(wfi, mask=None):
    """
    dynamically generate current pupil path for a given WFI instance
    """
    mask = (wfi._pupil_controller._get_filter_mask(wfi.filter) if mask is None
            else mask)
    detector = detector_substr(wfi.detector)

    base = wfi._pupil_controller._pupil_basepath
    file = wfi._pupil_controller.pupil_file_formatters[mask]

    return os.path.join(base, file).format(detector)

def test_WFI_psf():
    """
    Test that instantiating WFI works and can compute a PSF without
    raising any exceptions
    """
    wfi = roman.WFI()
    wfi.calc_psf(fov_pixels=4)

def test_WFI_filters():
    wfi = roman.WFI()
    filter_list = wfi.filter_list

    for filter in filter_list:
        wfi.filter = filter
        wfi.calc_psf(fov_pixels=4, oversample=1, nlambda=3)

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

        assert os.path.isfile(pupil_path(wfi)), f"Pupil file missing: {pupil_path(wfi)}"

        # Test detector change was successful
        assert wfi.detector == detector, "WFI detector was not set correctly"
        assert wfi.pupil == pupil_path(wfi), "pupil path was not set correctly"

        # Test pupil mask lock/unlock
        for mask in wfi.pupil_mask_list:
            # test lock
            wfi.lock_pupil_mask(mask)

            assert wfi.pupil == pupil_path(wfi, mask), "Pupil path was not set correctly"

            # introduce differing filter to modify
            wfi.filter = "PRISM" if mask != "PRISM" else "F062"

            assert wfi._pupil_controller._pupil_mask == wfi.pupil_mask, "Pupil mask was not set correctly"

            # test unlock
            wfi.unlock_pupil_mask()

            assert wfi.pupil == pupil_path(wfi), f"Pupil mask unlock failed"

        assert wfi._pupil_controller._auto_pupil, "Pupil is locked and should not be"
        assert wfi._pupil_controller._auto_pupil_mask, "Pupil mask is locked and should not be"

        # Test pupil lock/unlock
        with pytest.raises(FileNotFoundError) as err:
            assert wfi.lock_pupil("file_that_does_not_exist.fits"), "FileNotFoundError was not raised"

        this_file = __file__
        wfi.lock_pupil(this_file)
        assert wfi.pupil == this_file, "Pupil did not lock to proper file."

        wfi.unlock_pupil()
        assert wfi.pupil == pupil_path(wfi), f"Pupil unlock failed."

        assert wfi._pupil_controller._auto_pupil, "Pupil is locked and should  not be"
        assert wfi._pupil_controller._auto_pupil_mask, "Pupil mask is locked and should not be"

        # Test effect of changing the filter on pupil path
        for filter in wfi.filter_list:
            wfi.filter = filter

            assert wfi.pupil == pupil_path(wfi), f"Pupil was not set to correct value for filter {filter}"

    # Test persistence of pupil and pupil mask locks through a PSF calculation
    wfi2 = roman.WFI()
    wfi2.detector = detector
    valid_pos = (4000, 1000)
    wfi2.detector_position = valid_pos

    wfi2.filter = "F129"
    wfi2.lock_pupil_mask("GRISM")
    wfi2.filter = "F129"
    assert wfi2.pupil == pupil_path(wfi2, "GRISM"), "Pupil path was not set correctly"
    wfi2.calc_psf(monochromatic=1.3e-6, fov_pixels=4)

    assert wfi.pupil_mask == "GRISM", "Pupil mask changed during PSF calculation"
    assert wfi2.pupil == pupil_path(wfi2, "GRISM"), "Pupil path changed during PSF calculation"

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

def test_swapping_modes(wfi=None):

    if wfi is None:
        wfi = roman.WFI()

    # change detector string to match file format (e.g., "SCA01" -> "SCA_1")
    detector_substr = lambda det: f"{det[:3]}_{str(int((det[3:])))}"

    # dynamically generate current pupil path for a given WFI instance
    pupil_path = (
        lambda self, mask=None: os.path.join(
            self._pupil_controller._pupil_basepath,
            self._pupil_controller.pupil_file_formatters[self._pupil_controller._get_filter_mask(self.filter) if mask is None else mask]
        ).format(detector_substr(self.detector))
    )

    tests = [
        # [filter, mode, pupil_file]
        ['F146', 'imaging', pupil_path],
        ['F213', 'imaging', pupil_path],
        [PRISM_FILTERS[0], 'prism', pupil_path],
        [GRISM_FILTERS[0], 'grism', pupil_path],
    ]

    for test_filter, test_mode, test_pupil in tests:
        wfi.filter = test_filter

        fail_str = (f"failed on {test_filter}, {test_mode}, "
                    f"{test_pupil(wfi).split('/')[-1]}")

        assert wfi.filter == test_filter, fail_str
        assert wfi.mode == test_mode, fail_str
        assert wfi._current_aberration_file == wfi._aberration_files[test_mode], fail_str
        assert wfi.pupil == test_pupil(wfi), fail_str

def test_custom_aberrations():

    wfi = roman.WFI()

    # Use grism aberration_file for testing
    test_aberration_file = wfi._aberration_files['grism']

    # Test override
    # -------------
    wfi.lock_aberrations(test_aberration_file)

    for filter in wfi.filter_list:
        wfi.filter = filter
        assert wfi._current_aberration_file == test_aberration_file, "Filter change caused override to fail"

    # Test Release Override
    # ---------------------
    wfi.unlock_aberrations()
    assert wfi._aberration_files['custom'] is None, "Custom aberration file not deleted on override release."
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

    # Get min and max valid wavelengths from aberration file
    zern = Table.read(wfi._aberration_files[wfi.mode], format='ascii.csv')
    min_wv = zern['wavelength'][0] * 1e-6 # convert from micron to meter
    max_wv = zern['wavelength'][-1] * 1e-6

    # Test that get_aberration_terms() uses an approximated wavelength when
    # called with an out-of-bounds wavelength.
    too_lo_wv = min_wv * .9; too_hi_wv = max_wv / .9
    valid_wv = np.mean([min_wv, max_wv])

    assert allclose(det.get_aberration_terms(min_wv),
                    det.get_aberration_terms(too_lo_wv)), (
        "Aberration below wavelength range did not return closest value."
    )

    assert allclose(det.get_aberration_terms(max_wv),
                    det.get_aberration_terms(too_hi_wv)), (
        "Aberration above wavelength range did not return closest value."
    )

    # Test border pixels outside the ref data. In Cycle 9, (0, 37) is the first
    # pixel, so we check if (0, 0) is approximated to it as the nearest point.
    det.field_position = (0, 0)
    coefficients_outlier = det.get_aberration_terms(valid_wv)

    det.field_position = (0, 37)
    coefficients_data = det.get_aberration_terms(valid_wv)

    assert np.allclose(coefficients_outlier, coefficients_data), "nearest point extrapolation " \
                                                                 "failed for outlier field point"

def test_coronagraph_detector_position():
    """ Test existence of the Coronagraph detector position etc, and that you can't set it."""
    cor = roman.RomanCoronagraph()

    valid_pos = (512,512)
    assert cor.detector_position == valid_pos, "Coronagraph detector position isn't as expected"

    with pytest.raises(RuntimeError) as excinfo:
        cor.detector_position = valid_pos
    assert 'not adjustable' in str(excinfo.value), ("Failed to raise exception for"\
                                                        "trying to change Coronagraph detector position.")

def test_coronagraph_psf(display=False):
    """
    Just test that instantiating RomanCoronagraph works and can compute a PSF
    without raising any exceptions
    """
    char_spc = roman.RomanCoronagraph()
    char_spc.mode = 'CHARSPC_F660'

    #print('Reading instrument data from {:s}'.format(charspc._WebbPSF_basepath)
    #print('Filter list: {:}'.format(charspc.filter_list))

    monopsf = char_spc.calc_psf(nlambda=1, display=False)
    if display:
        roman.poppy.display_psf(monopsf)
