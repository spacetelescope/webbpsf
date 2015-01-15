import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits


import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())


from .. import webbpsf_core
import poppy

#poppy._log.setLevel(logging.INFO)


# The following functions are used in each of the test_<SI> files to
# test the individual SIs
def generic_output_test(iname):
    """ Basic test: Can we get PSFs of desired size and shape and sampling? 

    This is repeated for each SI (probably overkill but let's be thorough.)
    """

    _log.info("Testing image output sizes for %s " % iname)
    inst = webbpsf_core.Instrument(iname)
    pxscale = inst.pixelscale
    fov_arcsec = 5.0

    PSF = inst.calcPSF(nlambda=1, fov_pixels = 100, oversample=1)
    assert(PSF[0].data.shape[0] == 100)

    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    fov_pix = int(np.round(fov_arcsec / pxscale))
    assert(PSF[0].data.shape[0] == fov_pix)

    inst.options['parity'] = 'odd'
    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    assert( np.remainder(PSF[0].data.shape[0],2) == 1)

    inst.options['parity'] = 'even'
    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    assert( np.remainder(PSF[0].data.shape[0],2) == 0)

    # odd array, even oversampling = even
    inst.options['parity'] = 'odd'
    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=2)
    assert( np.remainder(PSF[0].data.shape[0],2) == 0)

    # odd array, odd oversampling = odd
    inst.options['parity'] = 'odd'
    PSF = inst.calcPSF(nlambda=1, fov_arcsec = fov_arcsec, oversample=3)
    assert( np.remainder(PSF[0].data.shape[0],2) == 1)



def do_test_source_offset(iname, distance=0.5,  nsteps=1, theta=0.0, tolerance=0.05, display=False):
    """ Test source offsets
    Does the star PSF center end up in the desired location?

    The tolerance threshold for success is by default 1/20th of a pixel 
    in the SI pixel units. But this can be adjusted by the calling function if needed.

    This is chosen somewhat arbitrarily as pretty good subpixel performance
    for most applications. Trying for greater accuracy would be limited by
    subpixel sampling in the simulations, as well as by the accuracy of the 
    centroid measuring function itself. 
    """
    _log.info("Calculating shifted image PSFs for "+iname)

    si = webbpsf_core.Instrument(iname)
    si.pupilopd=None

    oversample = 2

    # Calculations
    shift_req = []
    psfs = []

    # unshifted PSF
    #psfs.append(  nc.calcPSF(nlambda=1, oversample=oversample) )
    #shift_req.append(0)

    steps = np.linspace(0, distance, nsteps+1)
    for i, value in enumerate(steps):
        si.options['source_offset_r'] =  steps[i]
        si.options['source_offset_theta'] = theta
        #nc.options['source_offset_r'] = i*nc.pixelscale*5
        shift_req.append(si.options['source_offset_r'])
        psfs.append(  si.calcPSF(nlambda=1, oversample=oversample) )


    # Control case: an unshifted image
    cent0 = np.asarray(poppy.measure_centroid(psfs[0]))
    center_pix = (psfs[0][0].data.shape[0]-1)/2.0
    assert( abs(cent0[0] == center_pix) < 1e-3 )
    assert( abs(cent0[1] == center_pix) < 1e-3 )
    _log.info("Center of unshifted image: ({0:.3f}, {1:.3f}) pixels measured".format(*cent0))
    _log.info(" vs center of the array is ({0}, {0})".format(center_pix))

    if display:
        poppy.display_PSF(psfs[0])

    # Compare to shifted case(s)
    for i in range(1, nsteps+1):

        if display:
            poppy.display_PSF(psfs[i])

        cent = poppy.measure_centroid(psfs[i])
        rx = shift_req[i] * (-np.sin(theta*np.pi/180))
        ry = shift_req[i] * (np.cos(theta*np.pi/180))
        _log.info("   Shift_requested:\t(%10.3f, %10.3f) arcsec" % (rx, ry))
        shift = (cent-cent0) * (si.pixelscale/oversample)
        _log.info("   Shift_achieved: \t(%10.3f, %10.3f) arcsec" % (shift[1], shift[0]))
        assert( abs(rx -  shift[1]) <  (si.pixelscale*tolerance) )
        assert( abs(ry -  shift[0]) <  (si.pixelscale*tolerance) )



#------------------ generic infrastructure tests ----------------

def test_calcPSF_filter_arg():
    """ Tests the filter argument to the calcPSF function
    Can be used to set filter as same time as calculating a PSF
    (added for Pytest coverage completeness, even though is a minor bit of functionality)
    """
    nc = webbpsf_core.Instrument('NIRCam')
    nc.pupilopd=None

    nc.filter='F212N'
    psf1 = nc.calcPSF()

    nc.filter='F200W'
    psf2=nc.calcPSF(filter='F212N') # should override the filter setting just above

    assert(np.abs(psf1[0].data-psf2[0].data).max() < 1e-6)


def test_calcPSF_rectangular_FOV():
    """ Test that we can create rectangular FOVs """
    nc = webbpsf_core.Instrument('NIRCam')
    nc.pupilopd=None
    nc.filter='F212N'
 

    psf = nc.calcPSF(fov_arcsec=(2,4))
    assert(psf[0].data.shape[0]*2 == psf[0].data.shape[1])

    psf2 = nc.calcPSF(fov_pixels=(100,200), oversample=1)

    assert(psf2[0].data.shape==(100,200))


def test_cast_to_str():
    nc = webbpsf_core.NIRCam()

    assert str(nc)=='JWInstrument name=NIRCam'

def test_return_intermediates():
    import poppy
    import astropy.io.fits

    nc = webbpsf_core.NIRCam()
    nc.image_mask='maskswb'
    nc.pupil_mask='wedgelyot'

    psf, intermediates = nc.calcPSF(monochromatic=2e-6, return_intermediates=True)
    assert len(intermediates) == 4
    assert isinstance(intermediates[0], poppy.Wavefront)
    assert isinstance(psf, astropy.io.fits.HDUList)

def test_unicode_filter_names():
    """ See https://github.com/mperrin/webbpsf/issues/18
    Bug reported by Brian York in which unicode filternames made
    webbpsf 0.2.8 fail during atpy table lookup. Believed to actually be
    an atpy bug, now irrelevant since we're using astropy.table, but
    let's add an easy test case to be sure.
    """

    nc = webbpsf_core.NIRCam()
    nc.filter=unicode('f212n')
    psf_unicode = nc.calcPSF(nlambda=1)
    nc.filter='f212n'
    psf_str = nc.calcPSF(nlambda=1)

    assert np.array_equal(psf_unicode[0].data, psf_str[0].data)

#------------------    Utility Function Tests    ----------------------------


def test_instrument():
    nc = webbpsf_core.Instrument('NIRCam')


    try:
        import pytest
    except:
        _log.warning('Skipping last step in test_instrument because pytest is not installed.')
        return # We can't do this next test if we don't have the pytest.raises function.

    with pytest.raises(ValueError) as excinfo:
        tmp = webbpsf_core.Instrument('ACS')
    assert excinfo.value.message.startswith('Incorrect instrument name')


def test_calc_or_load_PSF(outputdir=None):
    if outputdir is None:
        import tempfile
        outputdir = tempfile.gettempdir()


    nc = webbpsf_core.NIRCam()


    filename =  os.path.join(outputdir, "test_calc_or_load_output.fits")
    if os.path.exists(filename): os.unlink(filename)

    webbpsf_core.calc_or_load_PSF(filename, nc, monochromatic=2e-6)

    assert os.path.exists(filename)

    #this one should not re-calc since the file already exists:
    # TODO - add some checking here of file modification date/times
    webbpsf_core.calc_or_load_PSF(filename, nc, monochromatic=2e-6)
    assert os.path.exists(filename)

    # this one should recalc since we explicitly ask it to
    webbpsf_core.calc_or_load_PSF(filename, nc, monochromatic=2e-6, clobber=True)
    assert os.path.exists(filename)

#--------------------------------------------------------------------------------

