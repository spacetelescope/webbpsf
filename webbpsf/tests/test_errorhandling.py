from __future__ import division, print_function, absolute_import, unicode_literals
# This file contains code for testing various error handlers and user interface edge cases,
# as opposed to testing the main body of functionality of the code.

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits


import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())



from .. import webbpsf_core
from .. import utils
import poppy

def _exception_message_starts_with(excinfo, message_body):
    if sys.version_info.major<3:
        return excinfo.value.message.startswith(message_body)
    else:
        return excinfo.value.args[0].startswith(message_body)


try:
    import pytest
    _HAVE_PYTEST = True
except:
    _HAVE_PYTEST = False

if _HAVE_PYTEST: 
    def test_calc_psf_catch_incompatible_oversampling():
        """ Test that we can create rectangular FOVs """
        nc = webbpsf_core.Instrument('NIRCam')
        nc.pupilopd=None
        nc.filter='F212N'

        with pytest.raises(ValueError) as excinfo:
            psf = nc.calc_psf(oversample=2, detector_oversample=10, fft_oversample=4)
        assert _exception_message_starts_with(excinfo,"You cannot specify simultaneously the oversample= option with the detector_oversample and fft_oversample options. Pick one or the other!")

# Not clear we need to test this; astropy logging framework will protect against invalid levels I think.
#    def test_invalid_logging_level():
#        with pytest.raises(ValueError) as excinfo:
#             utils.
#        assert excinfo.value.message == "You cannot specify simultaneously the oversample= option with the detector_oversample and fft_oversample options. Pick one or the other!"



    def test_invalid_masks():
        nc = webbpsf_core.NIRCam()

        # first, test case indepedencence. These are all converted to upper case internally & automatically
        nc.image_mask = 'maskswb'
        assert nc.image_mask =='MASKSWB'

        with pytest.raises(ValueError) as excinfo:
            nc.image_mask = 'junk'
        assert _exception_message_starts_with(excinfo,"Instrument NIRCam doesn't have an image mask called 'JUNK'.")

        # now repeat for pupil masks:
        nc.pupil_mask = 'circlyot'
        assert nc.pupil_mask =='CIRCLYOT'

        with pytest.raises(ValueError) as excinfo:
            nc.pupil_mask = 'junk'
        assert _exception_message_starts_with(excinfo, "Instrument NIRCam doesn't have a pupil mask called 'JUNK'.")




#--------------------------------------

# %% TODO %%  Why does including the following make the unit tests all fail?
if _HAVE_PYTEST and False: 
    # the following can only be run if we have pytest. This should more or less 
    # always be true since it's bundled with astropy
    def test_get_webbpsf_data_path_invalid():

        real_variable = os.getenv('WEBBPSF_PATH')
        real_setting =  conf.WEBBPSF_PATH

        # Test that we can overright the WEBBPSF_PATH setting here, by disabling/deleting the environment variable first if it is present
        conf.WEBBPSF_PATH = '/usr/'
        try:
            del os.environ['WEBBPSF_PATH']
        except:
            pass
        assert utils.get_webbpsf_data_path() == '/usr/'

        # Test that we get an error if we make the path invalid
        conf.WEBBPSF_PATH = 'some junk'
        with pytest.raises(IOError) as excinfo:
            tmp = utils.get_webbpsf_data_path()
        assert _exception_message_starts_with(excinfo,'Missing or invalid WEBBPSF_PATH to required data files')

        conf.WEBBPSF_PATH = real_setting
        os.putenv('WEBBPSF_PATH', real_variable)

        # and have we successfully put the environment back the way it was?
        assert (os.getenv('WEBBPSF_PATH') == real_variable)
        assert conf.WEBBPSF_PATH == real_setting


