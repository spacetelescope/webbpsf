from __future__ import division, print_function, absolute_import, unicode_literals
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

try:
    import pytest
    _HAVE_PYTEST = True
except:
    _HAVE_PYTEST = False

import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())

from .test_errorhandling import _exception_message_starts_with

from .. import webbpsf_core
from .. import utils
from .. import conf


def test_logging_restart():
    """ Test turning off and on the logging, and then put it back the way it was."""
    level = conf.logging_level

    conf.logging_level = 'NONE'
    utils.restart_logging()

    conf.logging_level = 'INFO'
    utils.restart_logging()
 
    conf.logging_level = level
    utils.restart_logging()


def test_logging_setup():
    """ Test changing log config settings, and then put it back the way it was."""
    loglevel = conf.logging_level
    logfn = conf.logging_filename

    _log.debug("Setting logging to OFF")
    utils.setup_logging(level=None, filename=None)
    _log.debug("Setting logging to WARN, and writing to file")
    utils.setup_logging(level='WARN', filename='test_log_file.txt')
    _log.debug("Setting logging to previous settings: {0}, {1}".format(loglevel, logfn))
    utils.setup_logging(level=loglevel, filename=logfn)

    try:
        import pytest
    except:
        _log.warning('Skipping last step in test_logging_setup because pytest is not installed.')
        return # We can't do this next test if we don't have the pytest.raises function.

    with pytest.raises(TypeError) as excinfo:
        utils.setup_logging(level='some junk')
    assert _exception_message_starts_with(excinfo,'Provided value for configuration item logging_level not valid:')


def test_diagnostic():

    res = utils.system_diagnostic()
    assert("webbpsf version" in res)
    assert("poppy version" in res)


def test_measure_strehl():
    # default NIRCam 2 micron PSF
    # FIXME this test will need reworking with the move to separate OTE and SI OPDs
    # for now I just doubled the tolerance to 6% instead of 3.
    wave=2.12e-6

    nc = webbpsf_core.NIRCam()
    nc.filter='F212N'
    defpsf = nc.calc_psf(nlambda=1)
    meas_strehl = utils.measure_strehl(defpsf, display=False, verbose=False)
    assert meas_strehl <= 1.0, 'measured Strehl cannot be > 1'
    assert meas_strehl >  0.7, 'measured Strehl is implausibly low for NIRCam'


    #compare to answer from Marechal approx on OPD rms WFE
    opdfile = fits.open(os.path.join(utils.get_webbpsf_data_path(),'NIRCam','OPD', nc.pupilopd))
    wfe_rms = opdfile[0].header['WFE_RMS']  # nm

    marechal_strehl = np.exp( -((wfe_rms *1e-9)/wave*(2*np.pi))**2)
    assert np.abs(meas_strehl-marechal_strehl) < 0.06, 'measured Strehl for that OPD file is too discrepant from the expected value from Marechal appoximation.'


    # and test a perfect PSF too
    perfnc = webbpsf_core.NIRCam()
    perfnc.filter='F212N'
    perfnc.pupilopd = None
    perfnc.include_si_wfe = False
    perfpsf = perfnc.calc_psf(nlambda=1)
    meas_perf_strehl = utils.measure_strehl(perfpsf, display=False, verbose=False)
    assert np.abs(meas_perf_strehl-1.0) < 0.01, 'measured Strehl for perfect PSF is insufficiently close to 1.0: {}'.format(meas_perf_strehl)
    assert meas_perf_strehl <= 1.0, 'measured Strehl cannot be > 1, even for a perfect PSF'


