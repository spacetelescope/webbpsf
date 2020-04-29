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
    otelm.thermal_slew(delta_time, start_angle, end_angle, scaling, case='EOL')
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
    otelm.thermal_slew(delta_time=1.0*u.day, case='EOL')
    max_truth = 4.13338e-08 / 1e-9 # Convert the max truth to units of nm
    assert np.isclose(np.max(otelm.opd)/1e-9, max_truth), "OPD max does not match expected value after 1 day slew."


def test_thermal_slew_reproducibility():
    """ If you call the thermal slew model multiple times, the OPD values should depend
    only on the LAST set of function call parameters. Not on the full time history.

    See issue #338
    """
    ote = webbpsf.opds.OTE_Linear_Model_WSS()

    ote.thermal_slew(12*u.hour, start_angle=-5, end_angle=45, case='EOL')
    opd1 = ote.opd.copy()

    ote.thermal_slew(24*u.hour, start_angle=-5, end_angle=45, case='EOL')
    opd2 = ote.opd.copy()

    ote.thermal_slew(12*u.hour, start_angle=-5, end_angle=45, case='EOL')
    opd3 = ote.opd.copy()

    assert np.allclose(opd1, opd2)==False, "OPDs expected to differ didn't"
    assert np.allclose(opd1, opd3), "OPDs expected to match didn't"


def test_update_opd():
    ''' The start of what should be many tests of this function'''

    # Test the very basics
    ote = webbpsf.opds.OTE_Linear_Model_WSS()
    ote.update_opd()
    assert np.max(ote.opd) == 0.0

    # can we add a deterministic frill drift?
    requested_wfe = 5
    ote.apply_frill_drift(requested_wfe)
    assert np.allclose(ote.rms(), requested_wfe, rtol=0.1), "Frill WFE amplitude not as expected"
    ote.apply_frill_drift(0.0)

    # can we add a deterministic IEC drift?
    requested_wfe = 15
    ote.apply_iec_drift(requested_wfe)
    assert np.allclose(ote.rms(), requested_wfe, rtol=0.1), "IEC WFE amplitude not as expected"

    # Todo test random drifts


def test_move_sur(plot=False):
    """ Test we can move mirrors using Segment Update Requests
    """
    import webbpsf
    import os
    import glob
    surdir = os.path.join(webbpsf.__path__[0], 'tests', 'surs')
    surs = glob.glob(surdir+'/*sur.xml')

    nrc = webbpsf.NIRCam()
    nrc.filter='F212N'
    nrc, ote = webbpsf.enable_adjustable_ote(nrc)
    ote.zero(zero_original=True)

    for s in surs:
        print("Testing "+s)
        ote.reset()
        ote.move_sur(s)
        # the coarse phasing SUR is a no-op after 3 groups; all others have some effect
        if 'coarse_phasing' not in s:
            assert not np.allclose(ote.segment_state, 0), "Expected some segments to be moved"

        ote.move_sur(s, reverse=True)
        assert np.allclose(ote.segment_state, 0), "Reversing moves didn't bring us back to zero"

    # Test moving one at a time. This test relies on specifics of what's in the image stacking SUR.
    s = glob.glob(surdir+'/example_image_stacking*sur.xml')[0]
    print("Testing moving one group at a time with "+s)
    ote.reset()
    import jwxml
    sur = jwxml.SUR(s)

    ngroups = len(sur.groups)
    oldstate = ote.segment_state.copy()

    for igrp in range(1,ngroups+1):
        print("Group {} should move segment {}".format(igrp, 2*igrp+6))
        ote.move_sur(s, group=igrp)

        movedsegs = np.abs((ote.segment_state - oldstate).sum(axis=1))
        assert (movedsegs!=0).sum()==1, "Only expected one segment to move"
        whichmoved = np.argmax(movedsegs)+1
        print ("Moved segment", whichmoved)
        assert whichmoved == 2*igrp+6, "An unexpected segment moved"
        oldstate = ote.segment_state.copy()
        if plot:
            psf = nrc.calc_psf(fov_pixels=256, add_distortion=False)
            plt.figure()
            ote.display_opd(title="After Group {}".format(igrp))
            plt.figure()
            webbpsf.display_psf(psf, ext=1, title="After Group {}".format(igrp))


