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

import sys
from astropy import config as _config

try:
    from .version import __version__
except ImportError:
    __version__ = ""

__minimum_python_version__ = "3.10"


class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple(
    (int(val) for val in __minimum_python_version__.split("."))
):
    raise UnsupportedPythonError(
        "webbpsf does not support Python < {}".format(__minimum_python_version__)
    )


# This tuple gives the *minimum* version of the WebbPSF data package
# required. If changes to the code and data mean WebbPSF won't work
# properly with an old data package, increment this version number.
# (It's checked against $WEBBPSF_DATA/version.txt)
DATA_VERSION_MIN = (1, 5, 0)


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `webbpsf`.
    """

    default_oversampling = _config.ConfigItem(
        4,
        "Default "
        + "oversampling factor: number of times more finely sampled than "
        + "an integer pixel for the grid spacing in the PSF calculation.",
    )

    default_output_mode = _config.ConfigItem(
        "both",
        "Should output include the oversampled PSF, a copy rebinned onto the integer detector spacing, or both?\
        Options: 'oversampled','detector','both' ",
    )
    default_fov_arcsec = _config.ConfigItem(
        5.0, "Default field of view size, in arcseconds per side of the square "
    )

    # Should be package settings:
    WEBBPSF_PATH = _config.ConfigItem(
        "from_environment_variable",
        "Directory path to data files required for WebbPSF calculations, such as OPDs and filter transmissions.\
        This will be overridden by the environment variable $WEBBPSF_PATH, if present.",
    )
    autoconfigure_logging = _config.ConfigItem(
        False,
        "Should WebbPSF configure logging for itself and POPPY? This adds handlers that report "
        "calculation progress and information",
    )
    logging_level = _config.ConfigItem(
        ["INFO", "DEBUG", "WARN", "ERROR", "CRITICAL", "NONE"],
        # (the default value is the first item in the options list)
        "Desired logging level for WebbPSF optical calculations.",
    )
    logging_filename = _config.ConfigItem(
        "none", "Desired filename to save log messages to."
    )
    logging_format_screen = _config.ConfigItem(
        "[%(name)7s] %(message)s", "Format for lines logged to the screen."
    )
    logging_format_file = _config.ConfigItem(
        "%(asctime)s [%(name)s:%(levelname)s] %(filename)s:%(lineno)d: %(message)s",
        "Format for lines logged to a file.",
    )


conf = Conf()

from . import utils  # noqa - must go after config
from . import trending  # noqa - must go after config
from .utils import setup_logging, restart_logging, system_diagnostic, measure_strehl  # noqa - must go after config

from poppy import (  # noqa
    display_psf,
    display_psf_difference,
    display_ee,
    measure_ee,
    display_profiles,
    radial_profile,
    measure_radial,
    measure_fwhm,
    measure_sharpness,
    measure_centroid,
    specFromSpectralType,
    fwcentroid,
)

from .webbpsf_core import (  # noqa
    instrument,
    SpaceTelescopeInstrument,
    JWInstrument,
    NIRCam,
    NIRISS,
    NIRSpec,
    MIRI,
    FGS,
)

from .opds import enable_adjustable_ote  # noqa

from .roman import WFI, RomanCoronagraph  # noqa

from .match_data import setup_sim_to_match_file  # noqa
