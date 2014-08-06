#import astropy.config
from astropy import config as _config

# Package-global configuration items here.
# segregated into a file "conf" to ease migration to the revamped config system in astropy >= 0.4

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

    if (path is None) or (path == 'unknown'):
        _log.critical("Fatal error: Unable to find required WebbPSF data files!")
        print """
 ********* WARNING ****** WARNING ****** WARNING ****** WARNING *************
 *                                                                          *
 *  WebbPSF requires several data files to operate.                         *
 *  These files could not be located automatically at this                  *
 *  time. Please download them to a location of your                        *
 *  choosing and either                                                     *
 *   - set the environment variable $WEBBPSF_PATH to point there, or        *
 *   - set the webbpsf_data_path variable in the configuration file         *
 *     path given above.                                                    *
 *                                                                          *
 *  See the WebbPSF documentation for more details.                         *
 *  WebbPSF will not be able to function properly until this has been done. *
 *                                                                          *
 ****************************************************************************
    """

    return path



def save_config():
    """ Save package configuration variables using the Astropy.config system 
    
    NOTE: This functionality is available as of astropy v0.3, but is being
    deprecated or changes as of astropy v0.4 See 

    http://astropy.readthedocs.org/en/latest/config/config_0_4_transition.html
    """

    from astropy.config import configuration
    configuration._save_config("webbpsf")

