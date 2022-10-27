"""
Tests for opds.py
"""
import os

from astropy.io import fits
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pysiaf
import pytest
import webbpsf
import matplotlib.pyplot as plt

# Set up a pinned pysiaf version so as not to break tests with any pysiaf value updates
prd_data_dir = pysiaf.constants.JWST_PRD_DATA_ROOT.rsplit('PRD', 1)[0]
PRD34_NRC = os.path.join(prd_data_dir, 'PRDOPSSOC-034/SIAFXML/SIAFXML/NIRCam_SIAF.xml')
PRD34_MIRI = os.path.join(prd_data_dir, 'PRDOPSSOC-034/SIAFXML/SIAFXML/MIRI_SIAF.xml')

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
# Coefficients for A1 based on 1 day after maximum slew with no scaling predicted
#   using above file
COEFFS_A1 = np.array([-3.52633363e-09, -2.90050902e-09, 1.25432196e-09, -7.43319098e-12,
                      -5.82462948e-11, -1.27115922e-10, -1.91541104e-12, 3.64760396e-11,
                      4.97176630e-13])
# Coefficients for A4 based on 5 hours after maximum slew with no scaling,
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
    # Assert the coefficients
    for coeff, truth in zip(coeffs, coeff_truth):
        #assert np.round(coeff, decimals=4) == np.round(truth, decimals=4)
        coeff /= 1e-9 # Convert to nm so we are not dealing with such small numbers
        truth /= 1e-9 # Convert to nm so we are not dealing with such small numbers
        assert np.isclose(coeff, truth), "Coeffs do not match expected value after day slew."


def test_thermal_slew_partial_angle():
    """ total slew should give same total amplitude if broken into smaller jumps"""

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

def test_sur_basics():
    # test we can create a null SUR
    sur = webbpsf.surs.SUR()
    assert sur.ngroups==1
    assert sur.nmoves==0
    assert isinstance(sur.describe(), str)


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
    assert np.abs(webbpsf.measure_centroid(psf)[0] - webbpsf.measure_centroid(psf_rm_ptt)[0]) > 40, "centroid should shift substantially with/without tip/tilt removal"


def test_apply_field_dependence_model():
    ''' Test to make sure the field dependence model is giving sensible output

    Checks cases comparing master chief ray, center of NIRCam, and center of NIRISS.

    Note, the steps for performing this test are a little subtle. We want to disable the
    SI and OTE global field dependence terms, and enable ONLY the OTE nominal field
    dependence. Thus there are several calls to manually set only the nominal field dep to True

    '''
    rms = lambda array, mask: np.sqrt((array[mask]**2).mean())

    # Get the OPD without any sort of field dependence
    ote = webbpsf.opds.OTE_Linear_Model_WSS(v2v3=None)
    # By default, an OTE LOM with zero OPD will implicitly also disable the field dependence.
    # For this test we don't want that so we re-enable it here:
    ote._include_nominal_field_dep = True

    opd_no_field_model = ote.opd.copy()

    mask = ote.get_transmission(0) != 0

    # Get the OPD at the zero field point of v2 = 0, v3 = -468 arcsec
    # Center of NIRCAM fields, but not physically on a detector.
    ote.v2v3 = (0, -468) * u.arcsec
    ote._apply_field_dependence_model(assume_si_focus=False)  # Do this test directly on the OTE OPD, without any implicit focus adjustments
    opd_zero_field = ote.opd.copy() * ote.get_transmission(0)
    rms1 =  rms(opd_no_field_model - opd_zero_field, mask)

    assert(rms1 < 7e-9), "OPDs expected to match didn't, at center field (zero field dependence)"

    # Get the OPD at some arbitrary nonzero field point (Center of NRC A)
    ote.v2v3 = (1.4, -8.2) * u.arcmin
    ote._apply_field_dependence_model(assume_si_focus=False)
    opd_arb_field = ote.opd.copy() * ote.get_transmission(0)
    rms2 = rms(opd_no_field_model - opd_arb_field, mask)

    assert(rms2 > 7e-9), "OPDs expected to differ didn't"
    assert np.isclose(rms2, 26.1e-9, atol=1e-9), "field dep OPD at center of NIRCam A was not the expected value"

    # Now we invoke this via an SI class, to show that works too:
    # Get the OPD at the center of NIRISS
    nis = webbpsf.NIRISS()
    nis.pupilopd = None  # disable any global WFE, so we just look at the field dependent part
    nis.detector_position = (1024, 1024)
    nis, ote_nis = webbpsf.enable_adjustable_ote(nis)
    ote_nis._include_nominal_field_dep = True  # Same as above, need this for test with pupilopd=None

    # Test if we directly invoke the OTE model, in this case also disabling SI focus implicit optimization
    ote_nis._apply_field_dependence_model(assume_si_focus=False)
    opd_nis_cen = ote_nis.opd.copy()
    rms3 = rms(opd_nis_cen, mask)
    # The value used in the following test is derived from this model itself, so it's a bit circular;
    # but at least this test should suffice to detect any unintended significant change in the
    # outputs of this model
    assert np.isclose(rms3, 36.0e-9, atol=1e-9), "Field-dependent OTE WFE at selected field point (NIRISS center) didn't match expected value (test case: explicit call, assume_si_focus=False)"

    # Now test as usd in a webbpsf calculation, implicitly, and with the defocus backout ON
    # The WFE here is slightly less, due to the focus optimization
    nis = webbpsf.NIRISS()
    nis.pupilopd = None  # disable any global WFE, so we just look at the field dependent part
    nis.detector_position = (1024, 1024)
    osys = nis.get_optical_system()
    ote = osys.planes[0]
    ote._include_nominal_field_dep = True  # Same as above, need this for test with pupilopd=None
    ote.update_opd()
    opd_nis_cen_v2 = ote.opd
    rms4 = rms(opd_nis_cen_v2, mask)
    assert np.isclose(rms4, 28.0e-9, atol=1e-9), "Field-dependent OTE WFE at selected field point (NIRISS center) didn't match expected value(test case: implicit call, assume_si_focus=True."


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
    assert(np.allclose(otelm._get_hexike_coeffs_from_smif(0., 0.), np.asarray([0.] * 9)))

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
    
    assert (np.allclose(otelm._get_hexike_coeffs_from_smif(1.0, 1.0)[3:], hexikes, rtol=1e-3))

    # Case 4: test at MIRIM_FP1MIMF field point, with pinned pysiaf version
    otelm.ote_ctrl_pt = pysiaf.Siaf('NIRCAM', filename=PRD34_NRC)['NRCA3_FP1'].reference_point('tel') * u.arcsec
    otelm.v2v3 = pysiaf.Siaf('MIRI', filename=PRD34_MIRI)['MIRIM_FP1MIMF'].reference_point('tel') * u.arcsec
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

    assert (np.allclose(otelm._get_hexike_coeffs_from_smif(dx, dy)[3:], hexikes, rtol=1e-3))

def test_segment_tilt_signs(fov_pix = 50, plot=False, npix=1024):
    """Test that segments move in the direction expected when tilted.

    The local coordinate systems are non-obvious, to say the least. This verifies
    sign conventions and coordinates are consistent in the linear optical model and
    optical propagation code.

    """

    if plot:
        fig, axs = plt.subplots(3, 5, figsize=(14,9))#, sharex = True, sharey = True)

    nrc = webbpsf.NIRCam()

    ote = webbpsf.opds.OTE_Linear_Model_WSS(npix=npix)
    nrc.include_si_wfe = False # not relevant for this test

    tilt = 1.0

    # We aim for relatively minimalist PSF calcs, to reduce test runtime
    psf_kwargs = {'monochromatic': 2e-6,
                  'fov_pixels': fov_pix,
                  'oversample': 1,
                  'add_distortion': False}

    # Which way are things expected to move?
    #
    # A1:  +X rotation -> -Y pixels (DMS), +Y rotation -> -X pixels
    # B1: +X rotation -> +Y pixels, +Y rotation -> +X pixels
    # C1: +X rotation -> +X/+Y pixels, +Y rotation -> -Y/+X pixels
    # (for C1, A/B means A is the sqrt(3)/2 component, B is the 1/2 component)
    #
    # The above derived from Code V models by R. Telfer, subsequently cross checked by Perrin

    for i, iseg in enumerate(['A1', 'B1', 'C1']):
        ote.zero()

        pupil = webbpsf.webbpsf_core.one_segment_pupil(iseg, npix=npix)

        ote.amplitude = pupil[0].data
        nrc.pupil = ote

        # CENTERED PSF:
        psf = nrc.calc_psf(**psf_kwargs)
        cen_ref = webbpsf.measure_centroid(psf, boxsize=10, threshold=1)

        ote.move_seg_local(iseg, xtilt=tilt)
        # XTILT PSF:
        psfx = nrc.calc_psf(**psf_kwargs)
        cen_xtilt = webbpsf.measure_centroid(psfx, boxsize=10, threshold=1)

        if iseg.startswith("A"):
            assert cen_xtilt[0] < cen_ref[0], "Expected A1:  +X rotation -> -Y pixels (DMS coords)"
            assert np.isclose(cen_xtilt[1], cen_ref[1], atol=1), "Expected A1:  +X rotation -> no change in X"
        elif iseg.startswith("B"):
            assert cen_xtilt[0] > cen_ref[0], "Expected B1: +X rotation -> +Y pixels (DMS coords)"
            assert np.isclose(cen_xtilt[1], cen_ref[1], atol=1), "Expected B1:  +X rotation -> no change in Y"
        elif iseg.startswith("C"):
            assert cen_xtilt[0] > cen_ref[0], "Expected C1: +X rotation -> +X/+Y pixels"
            assert cen_xtilt[1] > cen_ref[1], "Expected C1: +X rotation -> +X/+Y pixels"

        if plot:
            axs[i, 0].imshow(psf[0].data, norm=matplotlib.colors.LogNorm(vmax=1e-2, vmin=1e-5), origin="lower")
            axs[i, 0].set_title(iseg+": centered")
            axs[i, 0].axhline(y=fov_pix/2)
            axs[i, 0].axvline(x=fov_pix/2)
            # PLOT RESULTING OPD:
            im = axs[i, 1].imshow(ote.opd, vmin=-4e-6, vmax=4e-6, origin="lower")
            axs[i, 1].set_title("OPD (yellow +)")
            axs[i, 2].imshow(psfx[0].data, norm=matplotlib.colors.LogNorm(vmax=1e-2, vmin=1e-5), origin="lower")
            axs[i, 2].set_title(iseg+": xtilt {} um".format(tilt))
            axs[i, 2].axhline(y=fov_pix/2)
            axs[i, 2].axvline(x=fov_pix/2)


        ote.zero()
        ote.move_seg_local(iseg, ytilt=tilt)
        # YTILT PSF:
        psfy = nrc.calc_psf(**psf_kwargs)
        cen_ytilt = webbpsf.measure_centroid(psfy, boxsize=10, threshold=1)

        if iseg.startswith("A"):
            assert cen_ytilt[1] < cen_ref[1], "Expected A1:  +Y rotation -> -X pixels (DMS coords)"
            assert np.isclose(cen_ytilt[0], cen_ref[0], atol=1), "Expected A1:  +Y rotation -> no change in Y"
        elif iseg.startswith("B"):
            assert cen_ytilt[0] > cen_ref[0], "Expected B1: +Y rotation -> +X pixels(DMS coords)"
            assert np.isclose(cen_ytilt[0], cen_ref[0], atol=1), "Expected B1:  +Y rotation -> no change in Y"
        elif iseg.startswith("C"):
            assert cen_ytilt[0] < cen_ref[0], "Expected C1: +Y rotation -> -Y/+X pixels"
            assert cen_ytilt[1] > cen_ref[1], "Expected C1: +Y rotation -> -Y/+X pixels"

        # PLOT RESULTING OPD:
        if plot:
            im = axs[i, 3].imshow(ote.opd, vmin=-4e-6, vmax=4e-6, origin="lower")
            axs[i, 3].set_title("OPD (yellow +)")
            axs[i, 4].imshow(psfy[0].data, norm=matplotlib.colors.LogNorm(vmax=1e-2, vmin=1e-5), origin="lower")
            axs[i, 4].set_title(iseg+": ytilt {} um".format(tilt))
            axs[i, 4].axhline(y=fov_pix/2)
            axs[i, 4].axvline(x=fov_pix/2)

def test_segment_tilt_signs_2048npix():
    """ Re-run same test as above, but with a different value for npix

    This verifies the LOM works as expected for a size other than 1024 pixels
    """
    test_segment_tilt_signs(npix=2048)

def test_changing_npix():
    '''
    Test that using different npix will result in same PSF
    '''
    # Create a NIRCam instance using the default npix=1024
    nircam_1024 = webbpsf.NIRCam()
    nircam_1024.pupilopd = None # Set to none so I don't have to worry about making new OPDs
    psf_1024 = nircam_1024.calc_psf(oversample=2, nlambda=1, add_distortion=False)

    # Create a NIRCam instance using npix=2048
    npix = 2048
    nircam_2048 = webbpsf.NIRCam()
    nircam_2048.pupil = os.path.join(webbpsf.utils.get_webbpsf_data_path(),
                                     f'jwst_pupil_RevW_npix{npix}.fits.gz')
    nircam_2048.pupilopd = None # Set to none so I don't have to worry about making new OPDs
    psf_2048 = nircam_2048.calc_psf(oversample=2, nlambda=1, add_distortion=False)

    # Let's check individual pixel values, at least where the PSF is not too dim.
    # Check all pixels which have > 1e-6 of the total flux (we can safely ignore pixels with very low intensity)
    mask = psf_1024[0].data>1e-6
    assert np.allclose(psf_1024[0].data[mask], psf_2048[0].data[mask], rtol=0.01), 'Pixel values differ by more than 1%'

    # Let's check that the total flux in the PSF does not change much.
    #  (A small amount is acceptable and not surprising, since higher resolution improves the fidelity at which
    #   we model light that is scattered by segment edges to very wide angles outside of the simulated PSF FOV)
    assert np.isclose(psf_1024[0].data.sum(), psf_2048[0].data.sum(), rtol=0.005), "PSF total flux should not change much"

    # Let's also check a derived property of the whole PSF: the FWHM.
    # The FWHM should be very close to identical for the two PSFs.
    assert np.isclose(webbpsf.measure_fwhm(psf_1024), webbpsf.measure_fwhm(psf_2048), rtol=0.0001), "PSF FWHM should not vary for different npix"

def test_pupilopd_none():
    """Test that setting pupilopd=None does in fact result in no WFE
    In particular, this tests that setting opd to None also results in
    disabling the field-dependent component of the OTE linear model,
    as well as setting the global WFE component to zero.
    """

    nrc = webbpsf.NIRCam()
    nrc.pupilopd = None
    nrc.include_si_wfe = False

    ote_lom = nrc.get_optical_system().planes[0]
    assert ote_lom.rms()==0, "RMS WFE should be strictly 0"

    psf_small = nrc.calc_psf(fov_pixels=50, monochromatic=2e-6, add_distortion=False)
    centroid = webbpsf.measure_centroid(psf_small, relativeto='center')
    assert np.abs(centroid[0]) < 1e-5, "Centroid should be (0,0)"
    assert np.abs(centroid[1]) < 1e-5, "Centroid should be (0,0)"

def test_get_rms_per_segment():
    nrc0 = webbpsf.NIRCam()
    nrc, ote = webbpsf.enable_adjustable_ote(nrc0)
    rms_per_seg = webbpsf.opds.get_rms_per_segment(ote.opd)

    assert len(rms_per_seg)==18, "Wrong number of elements in result. Must be 18!"
    for seg in rms_per_seg:
        assert isinstance(seg, str)
        assert len(seg)==2
        assert isinstance(rms_per_seg[seg], float)
        if seg != 'C3':
            assert 10 < rms_per_seg[seg] < 100


        assert np.isclose(rms_per_seg[seg], ote.rms(seg))
