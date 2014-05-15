#!/usr/bin/env python
import os, sys

import logging
_log = logging.getLogger('webbpsf')

import ConfigParser 


from astropy.config import ConfigurationItem, get_config_dir, save_config

from . import conf


def restart_logging(verbose=True):
    """ Restart logging using the same settings as the last WebbPSF session, as stored in the configuration system. """

    level = conf.logging_level()
    lognames = ['webbpsf', 'poppy']
    if level.upper() =='NONE':
        # disable logging
        lev = logging.CRITICAL  # we don't generate any CRITICAL flagged log items, so
                                # setting the level to this is effectively the same as ignoring
                                # all log events. FIXME there's likely a cleaner way to do this.
        if verbose: print "No log messages will be shown from WebbPSF."
    else:
        lev = logging.__dict__[level.upper()] # obtain one of the DEBUG, INFO, WARN, or ERROR constants
        if verbose: print "WebbPSF log messages of level {0} and above will be shown.".format(level)

    for name in lognames:
        logging.getLogger(name).setLevel(lev)

    # set up screen logging
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
    if verbose: print("WebbPSF log outputs will be directed to the screen.")

    # set up file logging
    filename = conf.logging_filename()
    if filename.strip().lower() != 'none':
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

    level = level.upper()


    conf.logging_level.set(level)


    if filename is None: filename='none' # must be a string to write into the config system
    conf.logging_filename.set(filename)

    conf.logging_level.save()
    conf.logging_filename.save()

    restart_logging(verbose=True)



def check_for_new_install(force=False):
    """ Check for a new installation, and if so
    print a hopefully helpful explanatory message.
    """

    from .version import version
    if conf.last_version_ran() == '0.0' or force:
        from astropy.config import save_config, get_config_dir

        conf.last_version_ran.set(version)
        conf.last_version_ran.save()
        save_config('webbpsf') # save default values to text file

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
            {0}/webbpsf.cfg
    (unless such a config file was already present there)
    You can examine that file and change settings if desired.
    See the WebbPSF documentation for more detail. """.format(get_config_dir())

        # check for data dir?
        path_from_env_var = os.getenv('WEBBPSF_PATH') 
        WEBBPSF_PATH = ConfigurationItem('webbpsf_data_path','unknown','Directory path to data files required for WebbPSF calculations, such as OPDs and filter transmissions.', module='webbpsf.webbpsf_core')

        if path_from_env_var is not None:
            print """
    WebbPSF's required data files appear to be 
    installed at a path given by $WEBBPSF_PATH :
    {0} """.format(path_from_env_var)
        else:
            # the following will automatically print an error message if
            # the path is unknown in the config file.
            path_from_config = WEBBPSF_PATH()

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
        fftw3_version = 'yes, present'  # does not appear to have a version # string? 
    except:
        fftw3_version = 'not found'

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

    try:
        import pyfits
        pyfits_version = pyfits.__version__
    except:
        pyfits_version = 'not found'


    result = """
    OS: {os}
    Python version: {python}
    numpy version: {numpy}
    poppy version: {poppy}
    webbpsf version: {webbpsf}

    tkinter version: {tkinter}
    wxpython version: {wxpython}

    astropy version: {astropy}
    pysynphot version: {pysyn}
    pyfits version: {pyfits}
    PyFFTW version: {fftw3} """.format( os=platform.platform(), 
            numpy = numpy.__version__,
            python=sys.version.replace("\n"," "), 
            poppy=poppy.__version__, 
            webbpsf=version,
            tkinter=ttk_version,
            wxpython=wx_version,
            fftw3=fftw3_version,
            pysyn=pysynphot_version,
            astropy=astropy_version, 
            pyfits=pyfits_version)
    return result

