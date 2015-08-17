from webbpsf import wfirst

def test_WFI_psf():
    """
    Just test that instantiating WFI works and can compute a PSF without raising
    any exceptions
    """
    wi = wfirst.WFI()
    wi.calcPSF()
