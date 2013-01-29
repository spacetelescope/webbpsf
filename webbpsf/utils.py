#!/usr/bin/env python
import os, sys
from ._version import __version__


def _get_webbpsf_config_path():
    """ Find a reasonable location to store path and configuration information on 
    your current platform. 

    Code taken from http://stackoverflow.com/questions/1084697/how-do-i-store-desktop-application-data-in-a-cross-platform-way-for-python
    """
    APPNAME = "WebbPSF"

    import sys
    if sys.platform == 'darwin':
        try:
            from AppKit import NSSearchPathForDirectoriesInDomains
            # http://developer.apple.com/DOCUMENTATION/Cocoa/Reference/Foundation/Miscellaneous/Foundation_Functions/Reference/reference.html#//apple_ref/c/func/NSSearchPathForDirectoriesInDomains
            # NSApplicationSupportDirectory = 14
            # NSUserDomainMask = 1
            # True for expanding the tilde into a fully qualified path
            appdatapath = os.path.join(NSSearchPathForDirectoriesInDomains(14, 1, True)[0], APPNAME)
        except:
            appdatapath = os.path.expanduser(path.join("~", "." + APPNAME)) 
    elif sys.platform == 'win32':
        appdatapath = os.path.join(os.environ['APPDATA'], APPNAME)
    else:
        appdatapath = os.path.expanduser(path.join("~", "." + APPNAME)) 

    return appdatapath


def _check_config():
    """ Check for the existence and values of a WebbPSF configuration file
    If none is found, call _initialize_config() to initialize it
    """


def _initialize_config():
    """ Create a new blank webbpsf configuration file. 
    Also, query the user for whether they want to register for updates
    """
    pass

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

             
            SUBJECT = "WebbPSF user registration"
            TO = "mperrin@stsci.edu"
            FROM = "do-not-reply-automated-webbpsf@stsci.edu"
            HOST = 'smtp.stsci.edu'
            text = """The following user would like to be informed about
future versions of WebbPSF:
    email: {0}
    OS: {1}
    Python version: {2}
    poppy version: {3}
    webbpsf version: {4}
""".format( email, platform.platform(), sys.version, poppy.__version__, __version__)
    
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

