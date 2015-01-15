from webbpsf import wfirst

def test_wfirstimager_psf():
    """Just test that instantiating WFIRSTImager works and can compute a PSF"""
    wi = wfirst.WFIRSTImager()
    wi.calcPSF()
