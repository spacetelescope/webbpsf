import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits


import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())


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
    """ Test turning off and on the logging, and then put it back the way it was."""
    loglevel = conf.logging_level
    logfn = conf.logging_filename

    utils.setup_logging(level=None, filename='')
    utils.setup_logging(level='WARN', filename='test_log_file.txt')
    utils.setup_logging(level=loglevel, filename=logfn)


def test_diagnostic():

    res = utils.system_diagnostic()
    assert("webbpsf version" in res)
    assert("poppy version" in res)


# The following would need raw_input to be mocked in order
# to run successfully as a unit test. Not a priority for right now...
#def test_check_new_install():
#    utils.check_for_new_install(force=True)
