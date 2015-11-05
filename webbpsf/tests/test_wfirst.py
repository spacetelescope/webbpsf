from webbpsf import wfirst

def test_WFI_psf():
    """
    Just test that instantiating WFI works and can compute a PSF without raising
    any exceptions
    """
    wi = wfirst.WFI()
    wi.calcPSF()

def test_CGI_psf():
    """
    Just test that instantiating CGI works and can compute a PSF without raising
    any exceptions
    """
    charspc = wfirst.CGI()
    #wi.calcPSF()
    charspc.pupilopd = None
    charspc.filter = 'F770'
    charspc.apod_mask = 'CHARSPC'
    charspc.image_mask = 'CHARSPC_F770'
    charspc.pupil_mask = 'SPC26D88'

    print 'Reading instrument data from %s'%charspc._WebbPSF_basepath
    print 'Filter list:',charspc.filter_list

    monopsf = charspc.calcPSF(oversample=8, nlambda=1, display=True)
    wfirst.poppy.display_PSF(monopsf)

if __name__ == "__main__":
    """
    Just test that instantiating CGI works and can compute a PSF without raising
    any exceptions
    """

    charspc = wfirst.CGI()
    charspc.pupilopd = None
    charspc.filter = 'F770'
    charspc.apod_mask = 'CHARSPC'
    charspc.image_mask = 'CHARSPC_F770'
    charspc.pupil_mask = 'SPC26D88'

    print 'Reading instrument data from %s'%charspc._WebbPSF_basepath
    print 'Filter list:',charspc.filter_list

    charspc_monopsf = charspc.calcPSF(nlambda=1, fft_oversample=8, detector_oversample=2, fov_pixels=160, display=True)
    wfirst.poppy.display_PSF(charspc_monopsf)

    charspc_monopsf.writeto('charSPC_psf_test.fits', clobber=True)

    diskspc = wfirst.CGI()
    diskspc.pupilopd = None
    diskspc.filter = 'F465'
    diskspc.apod_mask = 'DISKSPC'
    diskspc.image_mask = 'DISKSPC_F465'
    diskspc.pupil_mask = 'SPC26D88'

    print 'Reading instrument data from %s'%diskspc._WebbPSF_basepath
    print 'Filter list:',diskspc.filter_list

    diskpsf = diskspc.calcPSF(nlambda=7, fft_oversample=8, detector_oversample=2, fov_pixels=250, display=True)
    wfirst.poppy.display_PSF(diskpsf)

    diskpsf.writeto('diskSPC_psf_test.fits', clobber=True)
