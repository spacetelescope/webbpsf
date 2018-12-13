# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
WebbPSF: Simulated Point Spread Functions for the James Webb Space Telescope
----------------------------------------------------------------------------

WebbPSF produces simulated PSFs for the James Webb Space Telescope, NASA's next
flagship infrared space telescope. WebbPSF can simulate images for any of the
four science instruments plus the fine guidance sensor, including both direct
imaging and coronagraphic modes.

Developed by Marshall Perrin and collaborators at STScI, 2010-2018.

Documentation can be found online at https://webbpsf.readthedocs.io/
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
# make use of astropy affiliate framework to set __version__, __githash__, and
# add the test() helper function
from ._astropy_init import *
# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

__minimum_python_version__ = "3.5"

class UnsupportedPythonError(Exception):
    pass

if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):
    raise UnsupportedPythonError("webbpsf does not support Python < {}".format(__minimum_python_version__))


# This tuple gives the *minimum* version of the WebbPSF data package
# required. If changes to the code and data mean WebbPSF won't work
# properly with an old data package, increment this version number.
# (It's checked against $WEBBPSF_DATA/version.txt)
DATA_VERSION_MIN = (0, 7, 0)

import astropy
from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `webbpsf`.
    """

    default_oversampling = _config.ConfigItem(4, 'Default '+
            'oversampling factor: number of times more finely sampled than '+
            'an integer pixel for the grid spacing in the PSF calculation.')

    default_output_mode = _config.ConfigItem('both', "Should output include the oversampled PSF, a copy rebinned onto the integer detector spacing, or both? Options: 'oversampled','detector','both' ")
    default_fov_arcsec = _config.ConfigItem( 5.0, "Default field of view size, in arcseconds per side of the square ")

# Should be package settings:
    WEBBPSF_PATH = _config.ConfigItem('from_environment_variable','Directory path to data files required for WebbPSF calculations, such as OPDs and filter transmissions. This will be overridden by the environment variable $WEBBPSF_PATH, if present.')
    autoconfigure_logging = _config.ConfigItem(
        False,
        'Should WebbPSF configure logging for itself and POPPY? This adds handlers that report '
        'calculation progress and information'
    )
    logging_level = _config.ConfigItem(
        ['INFO', 'DEBUG', 'WARN', 'ERROR', 'CRITICAL', 'NONE'],
        # (the default value is the first item in the options list)
        'Desired logging level for WebbPSF optical calculations.'
    )
    logging_filename = _config.ConfigItem("none", "Desired filename to save log messages to.")
    logging_format_screen = _config.ConfigItem(
        '[%(name)7s] %(message)s',
        'Format for lines logged to the screen.'
    )
    logging_format_file = _config.ConfigItem(
        '%(asctime)s [%(name)s:%(levelname)s] %(filename)s:%(lineno)d: %(message)s',
        'Format for lines logged to a file.'
    )

conf = Conf()

def _save_config():
    """ Save package configuration variables using the Astropy.config system

    NOTE: The functionality for saving config was was deprecated as of astropy v0.4
    See http://astropy.readthedocs.org/en/latest/config/config_0_4_transition.html

    This code is an undocumented workaround as advised by mdboom for the specific
    purpose of saving webbpsf GUI state, logging state, and related.
    """

    from astropy.config import configuration
    configuration._save_config("webbpsf")


from . import utils
from .utils import setup_logging, restart_logging, system_diagnostic, measure_strehl

if not _ASTROPY_SETUP_:
    if conf.autoconfigure_logging:
        restart_logging(verbose=True)

from poppy import ( display_psf, display_psf_difference, display_ee, measure_ee, # current names
        display_profiles, radial_profile,
        measure_radial, measure_fwhm, measure_sharpness, measure_centroid,
        specFromSpectralType, fwcentroid)

from .webbpsf_core import (Instrument, JWInstrument, NIRCam, NIRISS, NIRSpec,
    MIRI, FGS)

from .opds import enable_adjustable_ote

from .wfirst import WFI, CGI

from .jupyter_gui import show_notebook_interface

try:
    from .wxgui import wxgui
    _HAVE_WX_GUI = True
except ImportError:
    _HAVE_WX_GUI = False

try:
    from .tkgui import tkgui
    _HAVE_TK_GUI = True
except ImportError:
    _HAVE_TK_GUI = False



#if (_HAVE_WX_GUI or _HAVE_TK_GUI):

    #import warnings
    #warnings.warn("Warning: Neither Tk nor wx GUIs could be imported. "
    #              "Graphical interface disabled")
#else:
def gui(preferred='wx'):
    """ Start the WebbPSF GUI with the selected interface

    Parameters
    -------------
    preferred : string
        either 'wx' or 'ttk' to indicate which GUI toolkit should be started.


    """
    if preferred == 'wx' and _HAVE_WX_GUI:
        wxgui()
        pass
    elif preferred=='ttk' or _HAVE_TK_GUI:
        tkgui()
    else:
        raise NotImplementedError("Neither TK nor WX GUI libraries are available. Cannot start GUI.")
