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


# The following "truth" values" are based off the global focus and A1 Hexike coeffs
#   that will be returned for some time after a maximum slew using the time
#   constants and amplitudes established by the model in Fall 2018 (last updates
#   to otelm/thermal_OPD_fitting_parameters_9H_um.fits)
# Random scaling factor used below
SCALING_FACTOR = 0.5
# Coefficients for SM based on 1 day after maximum slew, first with no scaling
#   factor, and second with a scaling factor as specified above, all predicted with
#   above file. The below truth values are in units of METERS
GLOBAL_FOCUS = [-1.8251043541410904e-08]
GLOBAL_FOCUS2 = [GLOBAL_FOCUS[0] * SCALING_FACTOR]
# Coeffcients for A1 based on 1 day after maximum slew with no scaling predicted
#   using above file
COEFFS_A1 = np.array([-3.52633363e-09, -2.90050902e-09, 1.25432196e-09, -7.43319098e-12,
                      -5.82462948e-11, -1.27115922e-10, -1.91541104e-12, 3.64760396e-11,
                      4.97176630e-13])
# Coeffcients for A4 based on 5 days after maximum slew with no scaling,
#   start_angle=5. and end_angle=15., predicted using above file
COEFFS_A4 = np.array([4.08243373e-10, 1.89137803e-10, 1.24425015e-10, 4.63693785e-13,
                      -3.38635705e-11, 7.27484711e-12, -1.13484927e-13, 1.20634228e-12,
                      5.51330010e-14])
# Default slew angles
START_ANGLE = -5.
END_ANGLE = 45.

# Parameters to test the thermal model
tm_parameters = ([1 * u.day, 'SM', None, START_ANGLE, END_ANGLE, GLOBAL_FOCUS],
                 [1 * u.day, 'SM', SCALING_FACTOR, START_ANGLE, END_ANGLE, GLOBAL_FOCUS2],
                 [24, 'SM', None, START_ANGLE, END_ANGLE, GLOBAL_FOCUS],
                 [0.0 * u.day, 'SM', None, START_ANGLE, END_ANGLE, [0.0]],
                 [1 * u.day, 'A1', None, START_ANGLE, END_ANGLE, COEFFS_A1],
                 [0.0 * u.day, 'A1', None, START_ANGLE, END_ANGLE, np.zeros(9)],
                 [1.0 * u.day, 'D1', None, START_ANGLE, END_ANGLE, [0.0]],
                 [5 * u.hour, 'A4', None, 5., 15., COEFFS_A4])
@pytest.mark.parametrize('time, seg, scaling, start_angle, end_angle, coeff_truth', tm_parameters)
def test_get_thermal_slew_coeffs(time, seg, scaling, start_angle, end_angle,
                                 coeff_truth):
    """ Test that the OTE Thermal model is outputting the correct values
    These tests will go through the following (in order as listed in
    thermal_model_parameters):

     1. Test for SM with defaults
     2. Test for SM with a scaling factor
     3. Test for SM if no units specified for delta_time
     4. Test for SM with no delta_time
     5. Test for PM segment with defaults
     6. Test for PM segment with no delta_time
     7. Test for PM segment that is not in list of segnames
     8. Test for PM segment with start and end angles
    """
    delta_time = time
    # Create the thermal model
    otelm = webbpsf.opds.OTE_Linear_Model_WSS()
    otelm.thermal_slew(delta_time, start_angle, end_angle, scaling)
    coeffs = otelm._get_thermal_slew_coeffs(segid=seg)
    # Pull out coefficients
    if isinstance(coeffs, float):
        coeffs = [coeffs]
    # Assert the coefficents
    for coeff, truth in zip(coeffs, coeff_truth):
        #assert np.round(coeff, decimals=4) == np.round(truth, decimals=4)
        coeff /= 1e-9 # Convert to nm so we are not dealing with such small numbers
        truth /= 1e-9 # Convert to nm so we are not dealing with such small numbers
        assert np.isclose(coeff, truth), "Coeffs do not match expected value after day slew."


def test_thermal_slew_update_opd():
    ''' Test that running webbpsf.opds.OTE_Linear_Model_WSS.thermal_slew() will
        give the expected output'''
    otelm = webbpsf.opds.OTE_Linear_Model_WSS()
    otelm.thermal_slew(delta_time=1.0*u.day)
    max_truth = 4.13338e-08 / 1e-9 # Convert the max truth to units of nm
    assert np.isclose(np.max(otelm.opd)/1e-9, max_truth), "OPD max does not match expected value after 1 day slew."

def test_update_opd():
    ''' The start of what should be many tests of this function'''
    otelm = webbpsf.opds.OTE_Linear_Model_WSS()
    otelm.update_opd()
    assert np.max(otelm.opd) == 0.0
