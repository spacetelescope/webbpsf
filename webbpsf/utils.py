from __future__ import print_function
import os, sys

import logging
_log = logging.getLogger('webbpsf')

from . import conf

_DISABLE_FILE_LOGGING_VALUE = 'none'


class FilterLevelRange(object):
    def __init__(self, min_level, max_level):
        self.min_level = min_level
        self.max_level = max_level
    def filter(self, record):
        if record.levelno >= self.min_level and record.levelno <= self.max_level:
            return 1
        else:
            return 0


def restart_logging(verbose=True):
    """Restart logging using the same settings as the last WebbPSF
    session, as stored in the configuration system.

    Parameters
    ----------
    verbose : boolean
        Should this function print the new logging targets to
        standard output?
    """

    level = str(conf.logging_level).upper()
    lognames = ['webbpsf', 'poppy']

    root_logger = logging.getLogger()
    root_logger.handlers = []

    if level in ['DEBUG', 'INFO', 'WARN', 'ERROR', 'CRITICAL']:
        level_id = getattr(logging, level)  # obtain one of the DEBUG, INFO, WARN,
                                            # or ERROR constants
        if verbose:
            print("WebbPSF log messages of level {0} and above will be shown.".format(level))
    elif level == 'NONE':
        root_logger.handlers = []  # n.b. this will clear any handlers other libs/users configured
        return
    else:
        raise ValueError("Invalid logging level: {}".format(level))

    for name in lognames:
        logger = logging.getLogger(name)
        logger.setLevel(level_id)


    # set up screen logging
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.addFilter(FilterLevelRange(
        min_level=logging.DEBUG,
        max_level=logging.INFO
    ))

    stderr_handler = logging.StreamHandler(stream=sys.stderr)
    stderr_handler.addFilter(FilterLevelRange(
        min_level=logging.WARNING,
        max_level=logging.CRITICAL
    ))
    formatter = logging.Formatter(conf.logging_format_screen)
    stderr_handler.setFormatter(formatter)
    stdout_handler.setFormatter(formatter)

    root_logger.addHandler(stdout_handler)
    root_logger.addHandler(stderr_handler)

    if verbose:
        print("WebbPSF log outputs will be directed to the screen.")

    # set up file logging
    filename = conf.logging_filename
    if filename is None or filename.strip().lower() != _DISABLE_FILE_LOGGING_VALUE:
        hdlr = logging.FileHandler(filename)

        formatter = logging.Formatter(conf.logging_format_file)
        hdlr.setFormatter(formatter)

        root_logger.addHandler(hdlr)

        if verbose:
            print("WebbPSF log outputs will also be saved to file {}".format(filename))


def setup_logging(level='INFO', filename=None):
    """Allows selection of logging detail and output locations
    (screen and/or file)

    This is a convenience wrapper to Python's built-in logging package,
    as used by `webbpsf` and `poppy`. By default, this sets up log
    messages to be written to the screen, but the user can also
    request logging to a file.

    Editing the WebbPSF config file to set `autoconfigure_logging = True`
    (and any of the logging settings you wish to persist) instructs
    WebbPSF to apply your settings on import. (This is not
    done by default in case you have configured `logging` yourself
    and don't wish to overwrite your configuration.)

    For more advanced log handling, see the Python logging module's
    own documentation.

    Parameters
    -------------
    level : str
        Name of log output to show. Defaults to 'INFO', set to 'DEBUG'
        for more extensive messages, or to 'WARN' or 'ERROR' for fewer.
    filename : str, optional
        Filename to write the log output to. If not set, output will
        just be displayed on screen. (Default: None)

    Examples
    -----------

    >>> webbpsf.setup_logging(filename='webbpsflog.txt')

    This will save all log messages to 'webbpsflog.txt' in the current
    directory. If you later start another copy of webbpsf in a
    different directory, that session will also write to
    'webbpsflog.txt' in *that* directory. Alternatively you can specify
    a fully qualified absolute path to save all your logs to one
    specific file.

    >>> webbpsf.setup_logging(level='WARN')

    This will show only WARNING or ERROR messages on screen, and not
    save any logs to files at all (since the filename argument is None)
    """

    # implementation note: All this function actually does is apply the
    # defaults into the configuration system, then calls restart_logging to
    # do the actual work.
    level = str(level).upper()

    # The astropy config system will enforce the limited set of values for the logging_level
    # parameter by raising a TypeError on this next line if we feed in an invalid string.
    conf.logging_level = level

    if filename is None:
        # Use the string 'none' as a sentinel value for astropy.config
        filename = _DISABLE_FILE_LOGGING_VALUE

    conf.logging_filename = filename
    restart_logging(verbose=True)


MISSING_WEBBPSF_DATA_MESSAGE = """
 ***********  ERROR  ******  ERROR  ******  ERROR  ******  ERROR  ***********
 *                                                                          *
 *  WebbPSF requires several data files to operate.                         *
 *  These files could not be located automatically at this time, or this    *
 *  version of the software requires a newer set of reference files than    *
 *  you have installed.  For more details see:                              *
 *                                                                          *
 *           http://pythonhosted.org/webbpsf/installation.html              *
 *                                                                          *
 *  under "Installing the Required Data Files".                             *
 *  WebbPSF will not be able to function properly until the appropriate     *
 *  reference files have been downloaded to your machine and installed.     *
 *                                                                          *
 ****************************************************************************
"""


def get_webbpsf_data_path(data_version_min=None):
    """Get the WebbPSF data path

    Simply checking an environment variable is not always enough, since
    for packaging this code as a Mac .app bundle, environment variables are
    not available since .apps run outside the Terminal or X11 environments.

    Therefore, check first the environment variable WEBBPSF_PATH, and secondly
    check the configuration file in the user's home directory.

    If data_version_min is supplied (as a 3-tuple of integers), it will be
    compared with the version number from version.txt in the WebbPSF data
    package.
    """
    import os
    path_from_config = conf.WEBBPSF_PATH  # read from astropy configuration
    if path_from_config == 'from_environment_variable':
        path = os.getenv('WEBBPSF_PATH')
        if path is None:
            raise EnvironmentError("Environment variable $WEBBPSF_PATH is not set!")
    else:
        path = path_from_config

    # at minimum, the path must be a valid directory
    if not os.path.isdir(path):
        raise IOError("WEBBPSF_PATH ({}) is not a valid directory path!".format(path))

    if data_version_min is not None:
        # Check if the data in WEBBPSF_PATH meet the minimum data version
        version_file_path = os.path.join(path, 'version.txt')
        try:
            with open(version_file_path) as f:
                contents = f.read().strip()
                # keep only first 3 elements for comparison (allows "0.3.4.dev" or similar)
                parts = contents.split('.')[:3]
            version_tuple = tuple(map(int, parts))
        except (IOError, ValueError):
            raise EnvironmentError(
                "Couldn't read the version number from {}. (Do you need to update the WebbPSF "
                "data? See http://pythonhosted.org/webbpsf/installation.html#data-install "
                "for a link to the latest version.)".format(version_file_path)
            )

        if not version_tuple >= data_version_min:
            raise EnvironmentError(
                "WebbPSF data package has version {cur}, but {min} is needed. "
                "See http://pythonhosted.org/webbpsf/installation.html#data-install "
                "for a link to the latest version.".format(
                    cur=contents,
                    min='{}.{}.{}'.format(*data_version_min)
                )
            )

    return path


DIAGNOSTIC_REPORT = """
OS: {os}
Python version: {python}
numpy version: {numpy}
poppy version: {poppy}
webbpsf version: {webbpsf}

tkinter version: {tkinter}
wxpython version: {wxpython}

astropy version: {astropy}
pysynphot version: {pysyn}
pyFFTW version: {pyfftw}

Floating point type information for numpy.float:
{finfo_float}
Floating point type information for numpy.complex:
{finfo_complex}
"""

def system_diagnostic():
    """ return various helpful/informative information about the
    current system. For instance versions of python & available packages.

    Mostly undocumented function...
    """

    # There is probably a more clever way to do the following via introspection?

    import platform
    import os
    import poppy
    import numpy
    from .version import version
    try:
        import ttk
        ttk_version = ttk.__version__
    except ImportError:
        ttk_version = 'not found'

    try:
        import wx
        wx_version = wx.__version__
    except ImportError:
        wx_version = 'not found'

    try:
        import pyfftw
        pyfftw_version = pyfftw.version
    except ImportError:
        pyfftw_version = 'not found'

    try:
        import pysynphot
        pysynphot_version = pysynphot.__version__
    except ImportError:
        pysynphot_version = 'not found'


    try:
        import astropy
        astropy_version = astropy.__version__
    except ImportError:
        astropy_version = 'not found'

    result = DIAGNOSTIC_REPORT.format(
        os=platform.platform(),
        numpy=numpy.__version__,
        python=sys.version.replace("\n", " "),
        poppy=poppy.__version__,
        webbpsf=version,
        tkinter=ttk_version,
        wxpython=wx_version,
        pyfftw=pyfftw_version,
        pysyn=pysynphot_version,
        astropy=astropy_version,
        finfo_float=numpy.finfo(numpy.float),
        finfo_complex=numpy.finfo(numpy.complex),
    )
    return result
