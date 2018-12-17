"""
Tests for opds.py
"""
from astropy.io import fits
import astropy.units as u
import numpy as np
import pytest
import webbpsf

def test_enable_adjustable_ote():
    """ Some basic tests of the OTE LOM"""
    nc = webbpsf.NIRCam()
    nc, ote = webbpsf.enable_adjustable_ote(nc)

    # did this produce an OTE object?
    assert isinstance(ote, webbpsf.opds.OTE_Linear_Model_WSS), "Didn't get an OTE object back"

    # can we compute the rms?
    rms = ote.rms()

    # and can we move a mirror?

    ote.move_seg_local('B1', piston=10, clocking=200)

    assert ote.segment_state[6, 2] == 10, "Couldn't piston"
    assert ote.segment_state[6, 3] == 200, "Couldn't clock"

    # did that misalignment make things much worse?

    assert ote.rms() > rms*10, "Huge piston offset didn't make the WFE much worse"



GLOBAL_FOCUS = [-0.018251043541410904e-08]
GLOBAL_FOCUS_SCALED = [-0.018251043541410904e-08 * 0.5]
COEFFS_A1 = np.array([-3.52633363e-09, -2.90050902e-09, 1.25432196e-09, -7.43319098e-12,
                      -5.82462948e-11, -1.27115922e-10, -1.91541104e-12, 3.64760396e-11,
                      4.97176630e-13])
thermal_model_parameters = ([1 * u.day, 'SM', GLOBAL_FOCUS, None],
                            [1 * u.day, 'SM', GLOBAL_FOCUS_SCALED, .5],
                            [1440, 'SM', GLOBAL_FOCUS, None],
                            [1 * u.day, 'A1', COEFFS_A1, None],
                            [0.0 * u.day, 'SM', [0.0], None],
                            [0.0 * u.day, 'A1', np.zeros(9), None],
                            [1.0 * u.day, 'D1', [0.0], None])
@pytest.mark.parametrize('time, seg, coeff_truth, scaling', thermal_model_parameters)
def test_thermal_slew_opd(time, seg, coeff_truth, scaling):
    """ Test that the OTE Thermal model is outputting the correct values """
    delta_time = time
    # Create the thermal model
    coeffs = webbpsf.opds.thermal_slew_opd(delta_time, segid=seg, start_angle=-5.,
                                           end_angle=45., scaling=scaling)
    #thermal_model = webbpsf.opds.OteThermalModel(delta_time)
    # Pull out coefficients
    #coeffs = thermal_model.get_coeffs(seg)
    if isinstance (coeffs, float):
        coeffs = [coeffs]
    # Assert the coefficents
        for coeff, truth in zip(coeffs, coeff_truth):
                assert np.round(coeff, decimals=4) == np.round(truth, decimals=4)


def test_thermal_slew_update_opd():
    ''' Test that running webbpsf.opds.OTE_Linear_Model_WSS.thermal_slew() will
        give the expected output'''
    otelm = webbpsf.opds.OTE_Linear_Model_WSS()
    otelm.thermal_slew(delta_time=1.0*u.day)
    assert np.max(otelm.opd) == 4.1333852193439e-08

def test_update_opd():
    ''' The start of what should be many tests of this function'''
    otelm = webbpsf.opds.OTE_Linear_Model_WSS()
    otelm.update_opd()
    assert np.max(otelm.opd) == 0.0
