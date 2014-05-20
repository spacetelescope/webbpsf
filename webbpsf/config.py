#import astropy.config
from astropy import config as _config

# Package-global configuration items here.
# segregated into a file "conf" to ease migration to the revamped config system in astropy >= 0.4

class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `poppy`.
    """

    use_multiprocessing = _config.ConfigItem(False,
            'Should PSF calculations run in parallel using multiple processers'+
            'using the Python multiprocessing framework (if True; faster but '+
            'does not allow display of each wavelength) or run serially in a '+
            'single process(if False; slower but shows the calculation in '+
            'progress. Also a bit more robust.?)')


# Should probably be science state in astropy>=0.4 schema:

    default_oversampling = _config.ConfigItem(4, 'Default '+
            'oversampling factor: number of times more finely sampled than '+
            'an integer pixel for the grid spacing in the PSF calculation.')

    default_output_mode = _config.ConfigItem('both', "Should output include the oversampled PSF, a copy rebinned onto the integer detector spacing, or both? Options: 'oversampled','detector','both' ")
    default_fov_arcsec = _config.ConfigItem( 5.0, "Default field of view size, in arcseconds per side of the square ")

# Should be package settings:
    WEBBPSF_PATH = _config.ConfigItem('unknown','Directory path to data files required for WebbPSF calculations, such as OPDs and filter transmissions. This will be overridden by the environment variable $WEBBPSF_PATH, if present.')


# Settings cloned here from poppy
#   see _apply_settings_to_poppy below...
    use_multiprocessing= _config.ConfigItem( False, 'Should PSF calculations run in parallel using the Python multiprocessing framework (if True; faster but does not allow display of each wavelength) or run serially in a single process (if False; slower but shows the calculation in progress. Also a bit more robust.?)')
    n_processes= _config.ConfigItem(4, 'Maximum number of additional worker processes to spawn. PSF calculations are likely RAM limited more than CPU limited for higher N on modern machines, particularly for oversampling >=4. Set to 0 to have the computer attempt to choose an intelligent default based on available cores and RAM.')
    use_fftw = _config.ConfigItem(True, 'Use FFTW for FFTs (assuming it is available)?  Set to False to force numpy.fft always, True to try importing and using FFTW via PyFFTW.')


    # the default value is the first item in the options list:
    logging_level =  _config.ConfigItem(['INFO','DEBUG','WARN','ERROR','NONE'],'Desired logging level for WebbPSF optical calculations.')
    logging_filename =  _config.ConfigItem("none", "Desired filename to save log messages to.")
    last_version_ran =  _config.ConfigItem('0.0', 'Most recently used version of WebbPSF on this computer. This is used for detecting new or upgraded installations and providing some additional information to users.')



def _apply_settings_to_poppy():
    """Use webbpsf's settings to override any of the
    same settings in poppy. This is admittedly perhaps overbuilt to have identical
    settings in both packages, but the intent is to shield typical users of webbpsf
    from having to think about the existence of the underlying library. They can 
    just deal with one set of settings.
    """

    import poppy

    poppy.conf.use_multiprocessing =   conf.use_multiprocessing
    poppy.conf.n_processes =   conf.n_processes
    poppy.conf.use_fftw = conf.use_fftw
    poppy.conf.default_image_display_fov = conf.default_fov_arcsec
 
conf = Conf()

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

