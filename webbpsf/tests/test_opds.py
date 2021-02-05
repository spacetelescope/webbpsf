"""
Tests for opds.py
"""
from astropy.io import fits
import astropy.units as u
import numpy as np
import pysiaf
import pytest
import webbpsf
import matplotlib.pyplot as plt

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
# Coeffcients for A4 based on 5 hours after maximum slew with no scaling,
#   start_angle=5. and end_angle=15., predicted using above file
# Updated on 9/18/2020
COEFFS_A4 = np.array([ 3.89238932e-10,  1.80333109e-10,  1.18632814e-10,  4.42108030e-13,
                      -3.22871622e-11,  6.93619028e-12, -1.08202005e-13,  1.15018494e-12,
                      5.25664635e-14])
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


def test_thermal_slew_partial_angle():
    """ total slew shoudl give same total amplitude if broken into smaller jumps"""

    otelm = webbpsf.opds.OTE_Linear_Model_WSS()

    start_angle = -5
    mid_angle = 20
    end_angle = 45

    delta_time = 1 * u.hr

    # One large slew
    otelm.thermal_slew(delta_time, start_angle, end_angle, case='EOL')
    cf_full = np.array([otelm._get_thermal_slew_coeffs(segid=seg) for seg in otelm.segnames[0:18]])

    # Small slew 1
    otelm.thermal_slew(delta_time, start_angle, mid_angle, case='EOL')
    cf_all1 = np.array([otelm._get_thermal_slew_coeffs(segid=seg) for seg in otelm.segnames[0:18]])

    # Small slew 2
    otelm.thermal_slew(delta_time, mid_angle, end_angle, case='EOL')
    cf_all2 = np.array([otelm._get_thermal_slew_coeffs(segid=seg) for seg in otelm.segnames[0:18]])
    cf_tot = cf_all1 + cf_all2

    # Multiply by 1E9 so we're not dealing with small numbers
    assert np.allclose(1e9*cf_full, 1e9*cf_tot), "should get same total coefficients for one big slew or if broken into two parts"



def test_thermal_slew_update_opd():
    ''' Test that running webbpsf.opds.OTE_Linear_Model_WSS.thermal_slew() will
        give the expected output

        '''
    otelm = webbpsf.opds.OTE_Linear_Model_WSS()
    otelm.thermal_slew(delta_time=1.0*u.day, case='EOL')

    # the exact value expected is affected by which version of the linear model is used.
    if otelm._segment_masks_version < 3:
        # rev V pupil segment mask file, labeled as VERSION=2 in jwpupil_segments.fits
        expected_max = 41.3338  # nanometers, expected value for peak.
                                # value derived by kjbrooks based on thermal model coefficients
        expected_rms = 11.13    # nm
                                # value derived by mperrin based on evaluation of opd map in this case
    else:
        # rev W pupil segment mask file, labeled as VERSION=3 in jwpupil_segments.fits
        # Values here are by mperrin based on evaluation of the same exact linear model code as above
        # changing only the data file $WEBBPSF_DATA/jwpupil_segments.fits to the newer version
        expected_max = 40.7763  # nanometers, expected value for peak
        expected_rms = 11.24    # nm
    assert np.isclose(np.max(otelm.opd)/1e-9, expected_max, rtol=1e-3), "OPD max does not match expected value after 1 day slew."
    assert np.isclose(otelm.rms(), expected_rms, rtol=1e-3), "OPD rms does not match expected value after 1 day slew."


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

        
    # Test every DOF on A1-1 and SM and check the OTE state updated accordingly
    s = glob.glob(surdir+'/example_alldof_A1-SM_sur.xml')[0]
    print("Testing "+s)
    ote.reset()
    ote.move_sur(s)
    assert np.allclose(ote.segment_state[0],  [1, 2, 3, 4, 5, 6])
    assert np.allclose(ote.segment_state[-1], [1, 2, 3, 4, 5, 0])
    
    
    # Test moving one at a time. This test relies on specifics of what's in the image stacking SUR.
    s = glob.glob(surdir+'/example_image_stacking*sur.xml')[0]
    print("Testing moving one group at a time with "+s)
    ote.reset()
    sur = webbpsf.surs.SUR(s)

    ngroups = len(sur.groups)
    oldstate = ote.segment_state.copy()

    for igrp in range(1, ngroups+1):
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


def test_single_seg_psf(segmentid=1):
    """Test calculation of a single segment PSF, including options to remove piston/tip/tilt as used by MIRAGE

    """

    nrc = webbpsf.NIRCam()
    nrc.filter = 'F212N'
    nrc, ote = webbpsf.enable_adjustable_ote(nrc)
    ote.zero(zero_original=True)

    segname = webbpsf.constants.SEGNAMES_WSS_ORDER[segmentid-1][0:2]

    ote.move_seg_local(segname, xtilt=1, piston=-1)

    pupil = webbpsf.webbpsf_core.one_segment_pupil(segmentid)
    ote.amplitude = pupil[0].data


    psf = nrc.calc_psf(nlambda=1)

    ote.remove_piston = True
    ote.update_opd()
    psf_rm_piston = nrc.calc_psf(nlambda=1)
    assert np.allclose(psf[0].data, psf_rm_piston[0].data), "Piston removal should not affect the overall PSF"

    assert np.allclose( webbpsf.measure_centroid(psf), webbpsf.measure_centroid(psf_rm_piston)), "centroid should not shift"

    ote.remove_piston_tip_tilt = True
    ote.update_opd()
    psf_rm_ptt = nrc.calc_psf(nlambda=1)
    assert not np.allclose(psf[0].data, psf_rm_ptt[0].data), "Piston/Tip/Tip removal should shift the overall PSF"
    assert np.abs(webbpsf.measure_centroid(psf)[0] - webbpsf.measure_centroid(psf_rm_ptt)[0]) > 40, "centroid should shift susbtantially with/without tip/tilt removal"


def test_apply_field_dependence_model():
    ''' Test to make sure the field dependence model is giving sensible output'''

    # Get the OPD without any sort of field dependence
    ote = webbpsf.opds.OTE_Linear_Model_WSS(v2v3=None)
    opd_no_field_model = ote.opd.copy() * ote.get_transmission(0)

    # Get the OPD at the zero field point of v2 = 0, v3 = -468 arcsec
    # Center of NIRCAM fields, but not physically on a detector.
    ote.v2v3 = (0, -468) * u.arcsec
    ote._apply_field_dependence_model()
    opd_zero_field = ote.opd.copy() * ote.get_transmission(0)
    rms1 = np.sqrt(np.mean((opd_no_field_model - opd_zero_field) ** 2))

    # Get the OPD at some arbitrary nonzero field point
    ote.v2v3 = (1.8, -7) * u.arcmin
    ote._apply_field_dependence_model()
    opd_arb_field = ote.opd.copy() * ote.get_transmission(0)
    rms2 = np.sqrt(np.mean((opd_no_field_model - opd_arb_field) ** 2))

    assert(rms1 < 7e-9), "OPDs expected to match didn't, zero field"
    assert(rms2 > 7e-9), "OPDs expected to differ didn't"


def test_get_zernike_coeffs_from_smif():
    """ 
    Test that the OTE SM Influence function returns expected Hexike coefficients.
    """
    
    # Create an instance of the OTE linear model
    otelm = webbpsf.opds.OTE_Linear_Model_WSS()

    # Case 1: otelm.v2v3 is None, should return None
    otelm._apply_field_dependence_model()
    assert ( otelm._apply_field_dependence_model() is None)

    # Case 2: check coefficient at control point; should return zeros.
    assert( np.allclose(otelm._get_zernike_coeffs_from_smif(0., 0.), np.asarray([0.]*9) ))

    # Case 3: dx=1, dy=1, SM Poses all equal to 1 um
    telfer_zern = [-0.055279643, -0.037571947, -0.80840763, -0.035680581, -0.0036747300, 0.0033910640] # Taken from Telfer's tool
    # Convert Telfer's Zernikes to Hexikes:
    hexikes = [-telfer_zern[1], 
               2.*telfer_zern[0] - (60984./69531.)*telfer_zern[5], 
               telfer_zern[2], 
               (33./25)*telfer_zern[3], 
               (-33./25)*telfer_zern[4], 
               (1386./860.)*telfer_zern[5]]

    otelm.segment_state[-1, :] = 1.0
    
    assert (np.allclose(otelm._get_zernike_coeffs_from_smif(1.0, 1.0)[3:], hexikes, rtol=1e-3))

    # Case 4: test at MIRIM_FP1MIMF field point
    otelm.ote_ctrl_pt = pysiaf.Siaf('NIRCAM')['NRCA3_FP1'].reference_point('tel') *u.arcsec
    otelm.v2v3 = pysiaf.Siaf('MIRI')['MIRIM_FP1MIMF'].reference_point('tel') *u.arcsec
    telfer_zern_mirim_fp1mimf = np.asarray( [-0.25066019, 0.22840080, -0.53545999, -0.024227464, -0.0025191352, 0.00050082553]) # Taken from Telfer's tool
    # Convert Telfer's Zernikes to Hexikes:
    hexikes = hexikes = [-telfer_zern_mirim_fp1mimf[1], 
                         2.*telfer_zern_mirim_fp1mimf[0] - (60984./69531.)*telfer_zern_mirim_fp1mimf[5], 
                         telfer_zern_mirim_fp1mimf[2], 
                         (33./25)*telfer_zern_mirim_fp1mimf[3], 
                         (-33./25)*telfer_zern_mirim_fp1mimf[4], 
                         (1386./860.)*telfer_zern_mirim_fp1mimf[5]]
    
    otelm.segment_state[-1, :] = [300., 400., 100., 200., 5., 0.]
    dx =-(otelm.v2v3[0] - otelm.ote_ctrl_pt[0]).to(u.rad).value 
    dy = (otelm.v2v3[1] - otelm.ote_ctrl_pt[1]).to(u.rad).value

    assert (np.allclose(otelm._get_zernike_coeffs_from_smif(dx, dy)[3:], hexikes, rtol=1e-3))
    
