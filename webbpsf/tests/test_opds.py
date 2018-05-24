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

