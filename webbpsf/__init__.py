from .webbpsf import Instrument, NIRCam, NIRISS, NIRSpec,MIRI,FGS
from .webbpsfgui import gui  # quick function to start the GUI and run it instantly
from .webbpsfgui import WebbPSF_GUI  # actual GUI class 

from poppy import (display_PSF, display_PSF_difference, display_EE, display_profiles, radial_profile,
        measure_EE, measure_radial, measure_fwhm, measure_sharpness, measure_centroid, measure_strehl,
        measure_anisotropy, specFromSpectralType, fwcentroid)

from poppy._version import __version__
