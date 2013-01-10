from poppy import (display_PSF, display_PSF_difference, display_EE, display_profiles, radial_profile,
        measure_EE, measure_radial, measure_fwhm, measure_sharpness, measure_centroid, measure_strehl,
        measure_anisotropy, specFromSpectralType, fwcentroid, __version__)


from .webbpsf_core import Instrument, JWInstrument, NIRCam, NIRISS, NIRSpec,MIRI,FGS




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
