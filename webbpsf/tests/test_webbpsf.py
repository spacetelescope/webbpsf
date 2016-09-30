from __future__ import division, print_function, absolute_import, unicode_literals
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits


import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())


from .. import webbpsf_core
import poppy
from .test_errorhandling import _exception_message_starts_with



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

    # fov in pixels
    PSF = inst.calc_psf(nlambda=1, fov_pixels = 100, oversample=1)
    assert(PSF[0].data.shape[0] == 100)
    # fov in arcsec
    PSF = inst.calc_psf(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    fov_pix = int(np.round(fov_arcsec / pxscale))
    assert(PSF[0].data.shape[0] == fov_pix)

    # even and odd array sizes, no oversampling
    inst.options['parity'] = 'odd'
    PSF = inst.calc_psf(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    assert( np.remainder(PSF[0].data.shape[0],2) == 1)

    inst.options['parity'] = 'even'
    PSF = inst.calc_psf(nlambda=1, fov_arcsec = fov_arcsec, oversample=1)
    assert( np.remainder(PSF[0].data.shape[0],2) == 0)

    # odd array, even oversampling = even
    inst.options['parity'] = 'odd'
    PSF = inst.calc_psf(nlambda=1, fov_arcsec = fov_arcsec, oversample=2)
    assert( np.remainder(PSF[0].data.shape[0],2) == 0)

    # odd array, odd oversampling = odd
    inst.options['parity'] = 'odd'
    PSF = inst.calc_psf(nlambda=1, fov_arcsec = fov_arcsec, oversample=3)
    assert( np.remainder(PSF[0].data.shape[0],2) == 1)

def do_test_source_offset(iname, distance=0.5,  nsteps=1, theta=0.0, tolerance=0.05, monochromatic=None, display=False):
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

    if iname=='NIRSpec':
        si.image_mask = None # remove default MSA since it overcomplicates this test.

    oversample = 2

    # Calculations
    shift_req = []
    psfs = []

    # unshifted PSF
    #psfs.append(  nc.calc_psf(nlambda=1, oversample=oversample) )
    #shift_req.append(0)

    steps = np.linspace(0, distance, nsteps+1)
    for i, value in enumerate(steps):
        si.options['source_offset_r'] =  steps[i]
        si.options['source_offset_theta'] = theta
        #nc.options['source_offset_r'] = i*nc.pixelscale*5
        shift_req.append(si.options['source_offset_r'])
        psfs.append(  si.calc_psf(nlambda=1, monochromatic=monochromatic, oversample=oversample) )


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

        deltax =  abs(rx -  shift[1])
        deltay =  abs(ry -  shift[0])
        _log.info("   X offset:\t{0:.3f}\t\tTolerance:\t{1:.3f}".format(deltax, (si.pixelscale*tolerance)))
        assert( deltax  <  (si.pixelscale*tolerance) )
        _log.info("   Y offset:\t{0:.3f}\t\tTolerance:\t{1:.3f}".format(deltay, (si.pixelscale*tolerance)))
        assert( deltay  <  (si.pixelscale*tolerance) )



#------------------ generic infrastructure tests ----------------

def test_opd_selected_by_default():
    """
    Regression test for https://github.com/mperrin/webbpsf/issues/73

    Ensure an OPD map is set by default when instantiating an instrument
    """
    instruments = [
        webbpsf_core.NIRCam,
        webbpsf_core.MIRI,
        webbpsf_core.NIRSpec,
        webbpsf_core.NIRISS,
        webbpsf_core.FGS
    ]
    for InstrumentClass in instruments:
        ins = InstrumentClass()
        assert ins.pupilopd is not None, "No pupilopd set for {}".format(InstrumentClass)


def test_calc_psf_rectangular_FOV():
    """ Test that we can create rectangular FOVs """
    nc = webbpsf_core.Instrument('NIRCam')
    nc.pupilopd=None
    nc.filter='F212N'

    side = round(2/nc.pixelscale) *nc.pixelscale
    # pick something that can be done in integer pixels given NIRCam's sampling


    psf = nc.calc_psf(fov_arcsec=(side, 2*side))
    assert(psf[0].data.shape[0]*2 == psf[0].data.shape[1])

    psf2 = nc.calc_psf(fov_pixels=(100,200), oversample=1)

    assert(psf2[0].data.shape==(100,200))


def test_cast_to_str():
    nc = webbpsf_core.NIRCam()

    assert str(nc)=='<JWST: NIRCam>'


def test_return_intermediates():
    import poppy
    import astropy.io.fits

    nc = webbpsf_core.NIRCam()
    nc.image_mask='maskswb'
    nc.pupil_mask='wedgelyot'

    osys = nc._getOpticalSystem()

    psf, intermediates = nc.calc_psf(monochromatic=2e-6, return_intermediates=True)
    assert len(intermediates) == len(osys.planes)
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
    psf_unicode = nc.calc_psf(nlambda=1)
    nc.filter='f212n'
    psf_str = nc.calc_psf(nlambda=1)

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
    assert _exception_message_starts_with(excinfo,'Incorrect instrument name')


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

