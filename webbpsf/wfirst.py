import os.path
import poppy

import webbpsf_core

class WFIRSTInstrument(webbpsf_core.SpaceTelescopeInstrument):
    """
    WFIRSTInstrument contains data and functionality common to WFIRST
    instruments, such as setting the pupil shape
    """
    def __init__(self, *args, **kwargs):
        super(WFIRSTInstrument, self).__init__(*args, **kwargs)
        # the AFTA_symmetrical.fits pupil is 516x516, covering a 2.4m diameter
        pixelscale = 2.4 / 516.0
        self.pupil = poppy.FITSOpticalElement(
            transmission=os.path.join(self._WebbPSF_basepath, 'AFTA_symmetrical.fits'),
            pixelscale=pixelscale
        )
        self.pupilopd = None  # until we have some OPD maps and a FITS pupil of the right shape

class WFIRSTImager(WFIRSTInstrument):
    """
    WFIRSTImager represents to the to-be-named wide field imager
    for the WFIRST mission
    """
    def __init__(self):
        scale = 110e-3  # arcsec/px, WFIRST-AFTA SDT report v2 (p. 58)
        super(WFIRSTImager, self).__init__("WFIRSTImager", pixelscale=scale)
    def _validate_config(self):
        return True