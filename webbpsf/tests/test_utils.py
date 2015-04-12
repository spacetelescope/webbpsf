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


# The following would need raw_input to be mocked in order
# to run successfully as a unit test. Not a priority for right now...
#def test_check_new_install():
#    utils.check_for_new_install(force=True)
