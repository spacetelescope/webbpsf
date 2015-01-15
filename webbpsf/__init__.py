# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" 
WebbPSF: Simulated Point Spread Functions for the James Webb Space Telescope
-------------------------------------------------------------------------------

WebbPSf produces simulated PSFs for the James Webb Space Telescope, NASA's next flagship
infrared space telescope. WebbPSF can simulate images for any of the four science instruments plus the
fine guidance sensor, including both direct imaging and coronagraphic modes. 

Developed by Marshall Perrin at STScI, 2010-2012. 

Documentation can be found online at http://www.stsci.edu/jwst/software/webbpsf/

WebbPSF requires a large amount of input data for its simulations, including optical path difference (OPD) maps,
filter transmission curves, and coronagraph Lyot mask shapes. These data files are not included in this
source distribution available from PYPI. Please see the main WebbPSF web page, linked above, to download
the required data tar file.

This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
#if not _ASTROPY_SETUP_:
#    from example_mod import *
#import warnings

import astropy
from astropy import config as _config

class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `poppy`.
    """

    #use_multiprocessing = _config.ConfigItem(False,
    #        'Should PSF calculations run in parallel using multiple processers'+
    #        'using the Python multiprocessing framework (if True; faster but '+
    #        'does not allow display of each wavelength) or run serially in a '+
    #        'single process(if False; slower but shows the calculation in '+
    #        'progress. Also a bit more robust.?)')


# Should probably be science state in astropy>=0.4 schema:

    default_oversampling = _config.ConfigItem(4, 'Default '+
            'oversampling factor: number of times more finely sampled than '+
            'an integer pixel for the grid spacing in the PSF calculation.')

    default_output_mode = _config.ConfigItem('both', "Should output include the oversampled PSF, a copy rebinned onto the integer detector spacing, or both? Options: 'oversampled','detector','both' ")
    default_fov_arcsec = _config.ConfigItem( 5.0, "Default field of view size, in arcseconds per side of the square ")

# Should be package settings:
    WEBBPSF_PATH = _config.ConfigItem('from_environment_variable','Directory path to data files required for WebbPSF calculations, such as OPDs and filter transmissions. This will be overridden by the environment variable $WEBBPSF_PATH, if present.')


    # the default value is the first item in the options list:
    logging_level =  _config.ConfigItem(['INFO','DEBUG','WARN','ERROR','NONE'],'Desired logging level for WebbPSF optical calculations.')
    logging_filename =  _config.ConfigItem("none", "Desired filename to save log messages to.")
    last_version_ran =  _config.ConfigItem('0.0', 'Most recently used version of WebbPSF on this computer. This is used for detecting new or upgraded installations and providing some additional information to users.')

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


from poppy import (display_PSF, display_PSF_difference, display_EE, display_profiles, radial_profile,
        measure_EE, measure_radial, measure_fwhm, measure_sharpness, measure_centroid, measure_strehl,
        specFromSpectralType, fwcentroid)

from .webbpsf_core import (Instrument, JWInstrument, NIRCam, NIRISS, NIRSpec,
    MIRI, FGS)

from . import utils
from .utils import setup_logging #, _system_diagnostic, _check_for_new_install, _restart_logging

if not _ASTROPY_SETUP_:

    utils.check_for_new_install()    # display informative message if so.

    # this should display a warning to the user if they don't have WEBBPSF_PATH
    # defined in either the environment or in webbpsf.cfg
    _data_path = utils.get_webbpsf_data_path()


    utils.restart_logging()          # restart logging based on saved settings.



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



if not (_HAVE_WX_GUI or _HAVE_TK_GUI):
    warnings.warn("Warning: Neither Tk nor wx GUIs could be imported. "
                  "Graphical interface disabled")
else:
    def gui(preferred='wx'):
        """ Start the WebbPSF GUI with the selected interface

        Parameters
        -------------
        preferred : string
            either 'wx' or 'ttk' to indicate which GUI toolkit should be started.


        """
        if preferred == 'wx' and _HAVE_WX_GUI:
            wxgui()
            #try:
#            wxgui()
            #except:
                #raise ImportError("wxpython GUI for webbpsf not available ")
            pass
        elif preferred=='ttk' or _HAVE_TK_GUI:
            #try:
            tkgui()
            #except:
                #raise ImportError("ttk GUI for webbpsf not available")
        else:
            raise NotImplementedError("Neither TK nor WX GUI libraries are available. Cannot start GUI.")




