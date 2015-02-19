from webbpsf import wfirst

def test_wfirstimager_psf():
    """
    Just test that instantiating WFIRSTImager works and can compute a PSF without raising
    any exceptions
    """
    wi = wfirst.WFIRSTImager()
    wi.calcPSF()
