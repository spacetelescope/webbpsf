from webbpsf import wfirst


import pytest
# FIXME ensure this is no longer xfailed when https://github.com/mperrin/webbpsf/pull/75
# is merged
@pytest.mark.xfail
def test_wfirstimager_psf():
    """
    Just test that instantiating WFIRSTImager works and can compute a PSF without raising
    any exceptions
    """
    wi = wfirst.WFIRSTImager()
    wi.calcPSF()
