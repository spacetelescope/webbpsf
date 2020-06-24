import webbpsf


def test_coarse_jitter():
    """ Test the jitter options for Coarse Point mode
    Simple test to verify functionality, and sanity check that more
    blur means a less sharp PSF.
    """
    nrc = webbpsf.NIRCam()
    nrc.include_si_wfe = False
    nrc.options['jitter']='PCS=Coarse'
    psf_coarse = nrc.calc_psf(nlambda=1, add_distortion=False, fov_pixels=30)

    nrc.options['jitter']='PCS=Coarse_Like_ITM'
    psf_coarse_like_itm = nrc.calc_psf(nlambda=1, add_distortion=False, fov_pixels=30)

    # These test values are handwaved without any particular rigor
    assert psf_coarse[0].header['JITRSTRL'] < 0.5, "Coarse point blurs the image less than expected"

    assert psf_coarse_like_itm[0].header['JITRSTRL'] < 0.05, "Coarse point (Like ITM) blurs the image less than expected"

    assert psf_coarse[0].header['JITRTYPE'].startswith('PCS Coarse'), "Header keyword not as expected"
    assert 'JITRCPV2' in psf_coarse[0].header, "Header doesn't contain the expected coarse point offset V2 keyword"
    assert 'JITRCPV3' in psf_coarse[0].header, "Header doesn't contain the expected coarse point offset V3 keyword"

