
from poppy import (display_PSF, display_PSF_difference, display_EE, display_profiles, radial_profile,
        measure_EE, measure_radial, measure_fwhm, measure_sharpness, measure_centroid, measure_strehl,
        measure_anisotropy, specFromSpectralType, fwcentroid)


from .webbpsf_core import Instrument, JWInstrument, NIRCam, NIRISS, NIRSpec,MIRI,FGS

from ._version import __version__


from .utils import _initialize_config, _register


def setup_logging(filename=None, level='info'):
    """ Provide a default handler for log messages, if the 
    user has not already specified one. 

    This is a convenience wrapper to Python's built-in
    logging package, as used by webbpsf and poppy. 
    By default, this sets up log messages to be written to the screen,
    but the user can also request logging to a file. 

    For more advanced log handling, see the Python logging module's own documentation.
    
    Parameters
    -------------
    filename : str, optional
        Filename to write the log output to. If not set, output will just 
        be displayed on screen. 
    level : str
        name of log output to show. Defaults to 'info', set to 'debug' for
        more extensive messages. 


    Examples
    -----------

    >>> webbpsf.setup_logging('my_log_output.txt')

    
    """
    import logging
    import os

    if level.lower() =='info':
        lev = logging.INFO
    else:
        lev = logging.DEBUG

    lognames = ['webbpsf', 'poppy']

    for name in lognames:
        logging.getLogger(name).setLevel(lev)

    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')

    if filename is not None :
        hdlr = logging.FileHandler(filename)

        formatter = logging.Formatter('%(asctime)s %(name)-10s: %(levelname)-8s %(message)s')
        hdlr.setFormatter(formatter)

        for name in lognames:
            logging.getLogger(name).addHandler(hdlr)

        print("Log outputs will be saved to file "+filename)
    else:
        print("Log outputs will only be shown on screen.")




try: 
    from .wxgui import wxgui  
except:
    pass

try: 
    from .tkgui import tkgui  
except:
    pass



if ('tkgui' not in dir()) and ('wxgui' not in dir()):
    print ("Warning: Neither Tk nor wx GUIs could be imported. Graphical interface disabled")
else:

    def gui(preferred='wx'):
        """ Start the WebbPSF GUI with the selected interface

        Parameters
        -------------
        preferred : string
            either 'wx' or 'ttk' to indicate which GUI toolkit should be started.


        """
        if preferred == 'wx':
            try:
                wxgui()
            except:
                raise ImportError("wxpython GUI for webbpsf not available")
        elif preferred=='ttk':
            try:
                tkgui()
            except:
                raise ImportError("ttk GUI for webbpsf not available")
