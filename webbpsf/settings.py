import astropy.config

# Package-global configuration items here.
# segregated into this file because "namespaces are a honking great idea - let's do more of those!"



WEBBPSF_PATH = astropy.config.ConfigurationItem('webbpsf_data_path','unknown','Directory path to data files required for WebbPSF calculations, such as OPDs and filter transmissions. This will be overridden by the environment variable $WEBBPSF_PATH, if present.')

default_oversampling = astropy.config.ConfigurationItem('default_oversampling', 4, 'Default oversampling factor: number of times more finely sampled than an integer pixel for the grid spacing in the PSF calculation.')
default_output_mode = astropy.config.ConfigurationItem('default_output_mode', 'Both as FITS extensions', "Should output include the oversampled PSF, a copy rebinned onto the integer detector spacing, or both? Options: 'oversampled','detector','both' ")
default_fov_arcsec = astropy.config.ConfigurationItem('default_fov_arcsec', 5.0, "Default field of view size, in arcseconds per side of the square ")



# Settings cloned here from poppy
#   see _apply_settings_to_poppy below...
use_multiprocessing= astropy.config.ConfigurationItem('use_multiprocessing', False, 'Should PSF calculations run in parallel using the Python multiprocessing framework (if True; faster but does not allow display of each wavelength) or run serially in a single process (if False; slower but shows the calculation in progress. Also a bit more robust.?)')
n_processes= astropy.config.ConfigurationItem('n_processes', 4, 'Maximum number of additional worker processes to spawn. PSF calculations are likely RAM limited more than CPU limited for higher N on modern machines, particularly for oversampling >=4. Set to 0 to have the computer attempt to choose an intelligent default based on available cores and RAM.')
use_fftw = astropy.config.ConfigurationItem('use_fftw', True, 'Use FFTW for FFTs (assuming it is available)?  Set to False to force numpy.fft always, True to try importing and using FFTW via PyFFTW.')



def _apply_settings_to_poppy():
    """Use webbpsf's settings to override any of the
    same settings in poppy. This is admittedly perhaps overbuilt to have identical
    settings in both packages, but the intent is to shield typical users of webbpsf
    from having to think about the existence of the underlying library. They can 
    just deal with one set of settings.
    """

    import poppy

    poppy.settings.use_multiprocessing.set(  use_multiprocessing() )
    poppy.settings.n_processes.set(  n_processes() )
    poppy.settings.use_fftw.set(  use_fftw() )
    poppy.settings.default_image_display_fov.set (default_fov_arcsec() )
 

def get_webbpsf_data_path():
    """ Get webbpsf data path

    Simply checking an environment variable is not always enough, since 
    for packaging this code as a Mac .app bundle, environment variables are 
    not available since .apps run outside the Terminal or X11 environments.

    Therefore, check first the environment variable WEBBPSF_PATH, and secondly
    check a configuration file ~/.webbpsf in the user's home directory.
    """
    import os
    path = os.getenv('WEBBPSF_PATH') #, default= os.path.dirname(os.path.dirname(os.path.abspath(poppy.__file__))) +os.sep+"data" )
    if path is None:
        path = WEBBPSF_PATH() # read from astropy configuration system

        if path == 'unknown':
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
    """ Save package configuration variables using the Astropy.config system """
    astropy.config.save_config('webbpsf')

