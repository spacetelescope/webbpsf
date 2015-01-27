import os, sys
import ConfigParser

import logging
_log = logging.getLogger('webbpsf')

from . import conf


def restart_logging(verbose=True):
    """ Restart logging using the same settings as the last WebbPSF session, as stored in the configuration system. """

    level = str(conf.logging_level).upper()
    lognames = ['webbpsf', 'poppy']


    if level =='NONE':
        # disable logging
        lev = logging.CRITICAL  # we don't generate any CRITICAL flagged log items, so
                                # setting the level to this is effectively the same as ignoring
                                # all log events. FIXME there's likely a cleaner way to do this.
        if verbose: print "No log messages will be shown from WebbPSF."
    elif level in ['DEBUG', 'INFO','WARN','ERROR']:
        lev = logging.__dict__[level] # obtain one of the DEBUG, INFO, WARN, or ERROR constants
        if verbose: print "WebbPSF log messages of level {0} and above will be shown.".format(level)
    else:
        raise ValueError("Invalid logging level: "+level)
        return

    for name in lognames:
        logging.getLogger(name).setLevel(lev)

    # set up screen logging
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
    if verbose: print("WebbPSF log outputs will be directed to the screen.")

    # set up file logging
    filename = conf.logging_filename
    if filename is None or filename.strip().lower() != 'none':
        hdlr = logging.FileHandler(filename)

        formatter = logging.Formatter('%(asctime)s %(name)-10s: %(levelname)-8s %(message)s')
        hdlr.setFormatter(formatter)

        for name in lognames:
            logging.getLogger(name).addHandler(hdlr)

        if verbose: print("WebbPSF log outputs will also be saved to file "+filename)


def setup_logging(level='INFO',  filename=None):
    """ Allows selection of logging detail and output locations (screen and/or file)

    This is a convenience wrapper to Python's built-in logging package, as used
    by webbpsf and poppy.  By default, this sets up log messages to be written
    to the screen, but the user can also request logging to a file. 

    The settings applied here are stored persistently between sessions using the
    astropy.config system. In many cases, you can just configure these once when you
    first install webbpsf and then logging will continue to operate as desired 
    without any additional intervention.


    For more advanced log handling, see the Python logging module's own documentation.
    
    Parameters
    -------------
    level : str
        name of log output to show. Defaults to 'INFO', set to 'DEBUG' for
        more extensive messages, or to "WARN" or "ERROR" for fewer. 
    filename : str, optional
        Filename to write the log output to. If not set, output will just 
        be displayed on screen. 


    Examples
    -----------

    >>> webbpsf.setup_logging(filename='webbpsflog.txt')

    This will save all log messages to 'webbpsflog.txt' in the current directory.
    If you later start another copy of webbpsf in a different directory, that session
    will also write to 'webbpsflog.txt' in *that* directory. Alternatively you can 
    specify a fully qualified absolute path to save all your logs to one specific file.


    >>> webbpsf.setup_logging(level='WARN')

    This will show only WARNING or ERROR messages on screen, and not save any logs to
    files at all (since the filename argument is None)

    
    """

    # implementation note: All this function actually does is apply the
    # defaults into the configuration system, then calls restart_logging to
    # do the actual work.

    level = str(level).upper()

    # The astropy config system will enforce the limited set of values for the logging_level parameter
    # by raising a TypeError on this next line if we feed in an invalid string.
    conf.logging_level = level


    if filename is None: filename='none' # must be a string to store into the config system
    conf.logging_filename = filename

    #conf.logging_level.save()
    #conf.logging_filename.save()

    restart_logging(verbose=True)


def get_webbpsf_data_path():
    """ Get webbpsf data path

    Simply checking an environment variable is not always enough, since 
    for packaging this code as a Mac .app bundle, environment variables are 
    not available since .apps run outside the Terminal or X11 environments.

    Therefore, check first the environment variable WEBBPSF_PATH, and secondly
    check the configuration file in the user's home directory.
    """
    import os
    path = os.getenv('WEBBPSF_PATH') #, default= os.path.dirname(os.path.dirname(os.path.abspath(poppy.__file__))) +os.sep+"data" )
    if path is None:
        path = conf.WEBBPSF_PATH # read from astropy configuration system

    if (path is None) or not os.path.isdir(path):
        #(path == 'unknown') or (path == 'from_environment_variable'):
        _log.critical("Fatal error: Unable to find required WebbPSF data files!")
        print """
 *********  ERROR  ******  ERROR  ******  ERROR  ******  ERROR  *************
 *                                                                          *
 *  WebbPSF requires several data files to operate.                         *
 *  These files could not be located automatically at this                  *
 *  time. Please download them to a location of your                        *
 *  choosing and either                                                     *
 *   - set the environment variable $WEBBPSF_PATH to point there, or        *
 *   - set the WEBBPSF_PATH variable in your webbpsf.cfg config file        *
 *     (probably in ~/.astropy/config/ or similar location)                 *
 *                                                                          *
 *  See the WebbPSF documentation for more details.                         *
 *  WebbPSF will not be able to function properly until this has been done. *
 *                                                                          *
 ****************************************************************************
    """
        raise IOError('Missing or invalid WEBBPSF_PATH to required data files')

    return path


def check_for_new_install(force=False):
    """ Check for a new installation, and if so
    print a hopefully helpful explanatory message.
    """

    from .version import version as __version__
    if conf.last_version_ran == '0.0' or force:

        from . import _save_config
        import astropy.config

        conf.last_version_ran = __version__
        _save_config()
        #conf.last_version_ran.save()
        #save_config('webbpsf') # save default values to text file

        print """
  ***************************************************
  *           WebbPSF Initialization & Setup         *
  ****************************************************

    This appears to be the first time you have used WebbPSF. 
    
    Just so you know, there is a mailing list for users of
    webbpsf to keep you informed about any announcements of
    updates or enhancements to this software. If you would like
    to subscribe, please email 
        majordomo@stsci.edu
    with the message 
        subscribe webbpsf-users

    
    WebbPSF has some options that can be set using a 
    configuration file. An example configuration file with
    default values has been created in 
            {0}/webbpsf.{1}.cfg
    (unless such a config file was already present there)
    You can examine that file and change settings if desired.
    See the WebbPSF documentation for more detail. """.format(astropy.config.get_config_dir(), __version__)

        # check for data dir?
        path_from_env_var = os.getenv('WEBBPSF_PATH') 

        if path_from_env_var is not None:
            print """

    WebbPSF's required data files appear to be 
    installed at a path given by $WEBBPSF_PATH :
    {0} """.format(path_from_env_var)
        else:
            # the following will automatically print an error message if
            # the path is unknown in the config file.
            path_from_config = conf.WEBBPSF_PATH

            if path_from_config != 'unknown':
                print """
    WebbPSF's required data files appear to be 
    installed at a path given in the config file:
    {0} """.format(path_from_config)

        print """

    This message will not be displayed again.
    Press [Enter] to continue
    """
        any_key = raw_input()

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
    except:
        ttk_version = 'not found'

    try:
        import wx
        wx_version = wx.__version__
    except:
        wx_version = 'not found'

    try:
        import pyfftw
        pyfftw_version = pyfftw.version
    except:
        pyfftw_version = 'not found'

    try:
        import pysynphot
        pysynphot_version = pysynphot.__version__
    except:
        pysynphot_version = 'not found'


    try:
        import astropy
        astropy_version = astropy.__version__
    except:
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
