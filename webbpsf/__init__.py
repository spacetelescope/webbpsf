# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" WebbPSF
PSF simulations for the James Webb Space Telescope

"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

import warnings

import astropy


from .config import conf


from poppy import (display_PSF, display_PSF_difference, display_EE, display_profiles, radial_profile,
        measure_EE, measure_radial, measure_fwhm, measure_sharpness, measure_centroid, measure_strehl,
        specFromSpectralType, fwcentroid)


from .webbpsf_core import Instrument, JWInstrument, NIRCam, NIRISS, NIRSpec,MIRI,FGS

from . import utils
from .utils import setup_logging #, _system_diagnostic, _check_for_new_install, _restart_logging

utils.check_for_new_install()    # display informative message if so.

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

#def test( verbose=False ) :
#    import os, pytest
#
#    # find the directory where the test package lives
#    from . import tests
#    dir = os.path.dirname( tests.__file__ )
#
#    # assemble the py.test args
#    args = [ dir ]
#
#    # run py.test
#    try :
#        return pytest.main( args )
#    except SystemExit as e :
#        return e.code
#
