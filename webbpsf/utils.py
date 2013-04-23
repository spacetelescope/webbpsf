#!/usr/bin/env python
import os, sys
from ._version import __version__

import logging
_log = logging.getLogger('webbpsf')

import ConfigParser 

#class WebbPSFConfig(object):
#    """ Wrapper interface to a configuration file stored using the standard INI
#    file format and accessed using ConfigParser
#
#
#    Such a configuration file could look like:
#
#    [webbpsf]
#    datapath :  '/some/path'  ; optional, only used if $WEBBPSF_DATA is not defined
#    registered : True         ; user has already registered, don't ask again
#
#    [poppy]
#    use_multiproc : boolean ; should we spawn jobs over multi-processors by default?
#    multiproc_nprocess : integer    ; how many processes max? 
#
#
#
#    """
#    def __init__(self):
#        self.configdir = _get_webbpsf_config_path()
#        self.configfile = os.path.join(self.configdir, "webbpsf_config.ini")
#        self.defaults = {'registered': 'False',
#                         'first_invocation': 'True',
#                            'use_multiproc': 'False',
#                            'multiproc_nprocess': '4'}
#        if os.path.exists(self.configfile):
#            self._load_config()
#
#        else:
#            print "No configuration file detectet for WebbPSF - creating a default config now."
#            self._initialize_config()
#
#        if self.getboolean('webbpsf','first_invocation'):
#            wants_to_register, registered_ok = _register()
#            self.set('webbpsf', 'first_invocation', 'False')
#            self.set('webbpsf', 'registration_desired', str(wants_to_register))
#            self.set('webbpsf', 'registered', str(registered_ok))
#            self.save()
#
#
#
#    def get(self, section, name):
#        return self.parser.get(section, name)
#    def getboolean(self, section, name):
#        return self.parser.getboolean(section, name)
#
# 
#    def set(self, section, name, value):
#        return self.parser.set(section, name, value)
#
#
#    def _load_config(self):
#        """ 
#        read config from a file 
#        """
#        self.parser = ConfigParser.SafeConfigParser(defaults=self.defaults)
#        self.parser.read(self.configfile)
#        _log.info("Configuration file read from "+self.configfile)
#
#    def save(self):
#        """ Save the configuration information to disk """
#
#        filehandle = open(self.configfile, "w")
#        self.parser.write(filehandle)
#        filehandle.close()
#        _log.info("Configuration information saved to "+self.configfile)
#
#
#    def _initialize_config(self):
#        """ Create a new blank webbpsf configuration file. 
#        Also, query the user for whether they want to register for updates
#        """
#        if not os.path.isdir(self.configdir):
#            try:
#                os.makedirs(self.configdir)
#            except:
#                raise IOError("Could not create path to store configuration information: "+self.configdir)
#
#        _log.info("No configuration file found; creating "+self.configfile)
#
#        self.parser = ConfigParser.SafeConfigParser(defaults=self.defaults)
#        self.parser.add_section('poppy')
#        self.parser.add_section('webbpsf')
#        self.parser.add_section('gui')
#
#        self.save()
#

#def _get_webbpsf_config_path():
#    """ Find a reasonable location to store path and configuration information on 
#    your current platform. 
#
#    Code taken from http://stackoverflow.com/questions/1084697/how-do-i-store-desktop-application-data-in-a-cross-platform-way-for-python
#    """
#    APPNAME = "WebbPSF"
#
#    import sys
#    if sys.platform == 'darwin':
#        try:
#            from AppKit import NSSearchPathForDirectoriesInDomains
#            # http://developer.apple.com/DOCUMENTATION/Cocoa/Reference/Foundation/Miscellaneous/Foundation_Functions/Reference/reference.html#//apple_ref/c/func/NSSearchPathForDirectoriesInDomains
#            # NSApplicationSupportDirectory = 14
#            # NSUserDomainMask = 1
#            # True for expanding the tilde into a fully qualified path
#            appdatapath = os.path.join(NSSearchPathForDirectoriesInDomains(14, 1, True)[0], APPNAME)
#        except:
#            appdatapath = os.path.expanduser(path.join("~", "." + APPNAME)) 
#    elif sys.platform == 'win32':
#        appdatapath = os.path.join(os.environ['APPDATA'], APPNAME)
#    else:
#        appdatapath = os.path.expanduser(path.join("~", "." + APPNAME)) 
#
#    return appdatapath
#

from astropy.config import ConfigurationItem, get_config_dir, save_config

LOGGING_LEVEL = ConfigurationItem('logging_level',['INFO','DEBUG','WARN','ERROR','NONE'],'Desired logging level for WebbPSF optical calculations.')
LOGGING_FILENAME = ConfigurationItem('logging_filename',"none", "Desired filename to save log messages to.")
LAST_VERSION_RAN = ConfigurationItem('last_version_ran','0.0', 'Most recently used version of WebbPSF on this computer. This is used for detecting new or upgraded installations and providing some additional information to users.')


def _restart_logging(verbose=True):
    """ Restart logging using the same settings as the last WebbPSF session, as stored in the configuration system. """

    level = LOGGING_LEVEL()
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
    if verbose: print(" WebbPSF log outputs will be directed to the screen.")

    # set up file logging
    filename = LOGGING_FILENAME()
    if filename.strip().lower() != 'none':
        hdlr = logging.FileHandler(filename)

        formatter = logging.Formatter('%(asctime)s %(name)-10s: %(levelname)-8s %(message)s')
        hdlr.setFormatter(formatter)

        for name in lognames:
            logging.getLogger(name).addHandler(hdlr)

        if verbose: print(" WebbPSF log outputs will also be saved to file "+filename)




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


    LOGGING_LEVEL.set(level)


    if filename is None: filename='none' # must be a string to write into the config system
    LOGGING_FILENAME.set(filename)

    LOGGING_LEVEL.save()
    LOGGING_FILENAME.save()

    _restart_logging(verbose=True)


#    lognames = ['webbpsf', 'poppy']
#    if level.upper() =='NONE':
#        # disable logging
#        lev = logging.CRITICAL  # we don't generate any CRITICAL flagged log items, so
#                                # setting the level to this is effectively the same as ignoring
#                                # all log events. FIXME there's likely a cleaner way to do this.
#        print "No log messages will be shown."
#    else:
#        lev = logging.__dict__[level.upper()] # obtain one of the DEBUG, INFO, WARN, or ERROR constants
#        print "Log messages of level {0} and above will be shown.".format(level)
#
#    for name in lognames:
#        logging.getLogger(name).setLevel(lev)
#
#    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
#    print("Log outputs will be directed to the screen.")
#
#    if filename is not None :
#        hdlr = logging.FileHandler(filename)
#
#        formatter = logging.Formatter('%(asctime)s %(name)-10s: %(levelname)-8s %(message)s')
#        hdlr.setFormatter(formatter)
#
#        for name in lognames:
#            logging.getLogger(name).addHandler(hdlr)
#
#        print("Log outputs will also be saved to file "+filename)
#
#


def _check_for_new_install():
    from ._version import __version__
    if LAST_VERSION_RAN() == '0.0':
        print """
    *********************************************
    *           WebbPSF Initialization          *
    *********************************************

    This appears to be the first time you have used WebbPSF. 
    
    Just so you know, there is a mailing list for users of
    webbpsf to keep you informed about any announcements of
    updates or enhancements to this software. If you would like
    to subscribe, please email 
        majordomo@stsci.edu
    with the message 
        subscribe webbpsf-users

    This message will not be displayed again.
    Press [Enter] to continue
    """

        any_key = raw_input()

        LAST_VERSION_RAN.set(__version__)
        LAST_VERSION_RAN.save()
    
    #result = query_yes_no("Register?", default='yes')


def _register():
    """ Submit registration email
    
    Returns 2 bools, for "registration desired?" "if desired, successful?"
    """


    print """
    *********************************************
    *           WebbPSF Initialization          *
    *********************************************

This appears to be the first time you have used WebbPSF. 

Would you like to register your email address to 
stay informed of future versions, updates, etc? 
This will also register some basic information about
your system (OS, Python version, WebbPSF version, etc.)
to help us better support this software. 
"""
    result = query_yes_no("Register?", default='yes')

    if result =='no':
        return (False, False)
    elif result =='yes':
        confirmed = 'no'
        while confirmed == 'no':
            email = raw_input("Please enter your email address: ")
            email = email.strip()
            confirmed = query_yes_no("You entered '%s'; is this correct? " % email, default='yes')

        else:
                        
            import smtplib
            import string
            import platform
            import poppy

            print "now sending registration email to mperrin@stsci.edu."
             
            SUBJECT = "WebbPSF user registration"
            TO = "mperrin@stsci.edu"
            FROM = "do-not-reply-automated-webbpsf@stsci.edu"
            HOST = 'smtp.stsci.edu'
            text = """The following user would like to be informed about
future versions of WebbPSF:
    email: {0}
""".format( email)
            text += _system_diagnostic()
    
            BODY = string.join((
                    "From: %s" % FROM,
                    "To: %s" % TO,
                    "Subject: %s" % SUBJECT ,
                    "",
                    text
                    ), "\r\n")
            server = smtplib.SMTP(HOST)
            try:
                server.sendmail(FROM, [TO], BODY)
                server.quit()
                print """
    Registration email sent.
    """
                return (True, True)
            except:
                print """ Error in sending registration email! Check your internet connection and try again later?"""
                return (True, False)

        
     


def _system_diagnostic():
    """ return various helpful/informative information about the
    current system. For instance versions of python & available packages.

    Mostly undocumented function...
    """

    # There is probably a more clever way to do the following via introspection?

    import platform
    import os 
    import poppy
    import numpy
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
        import fftw3
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
    FFTW3 version: {fftw3} """.format( os=platform.platform(), 
            numpy = numpy.__version__,
            python=sys.version.replace("\n"," "), 
            poppy=poppy.__version__, 
            webbpsf=__version__,
            tkinter=ttk_version,
            wxpython=wx_version,
            fftw3=fftw3_version,
            pysyn=pysynphot_version,
            astropy=astropy_version, 
            pyfits=pyfits_version)
    return result


## {{{ http://code.activestate.com/recipes/577058/ (r2)
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.
    
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":"yes",   "y":"yes",  "ye":"yes",
             "no":"no",     "n":"no"}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return default
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")
## end of http://code.activestate.com/recipes/577058/ }}}

