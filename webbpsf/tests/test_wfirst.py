from webbpsf import wfirst


try:
    import pytest
    _HAVE_PYTEST = True
except:
    _HAVE_PYTEST = False

if _HAVE_PYTEST: 

    @pytest.mark.xfail # Only actually expected to fail if missing pysynphot
    def test_wfirstimager_psf():
        """Just test that instantiating WFIRSTImager works and can compute a PSF"""
        wi = wfirst.WFIRSTImager()
        wi.calcPSF()
