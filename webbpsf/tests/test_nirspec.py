import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import astropy.units as u
import pysiaf

import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())

from .. import webbpsf_core
import poppy


#------------------    NIRSpec Tests    ----------------------------

from .test_webbpsf import generic_output_test, do_test_source_offset, do_test_set_position_from_siaf

test_nirspec= lambda : generic_output_test('NIRSpec')

# Use a larger than typical tolerance when testing NIRSpec offsets. The
# pixels are so undersampled (0.1 arcsec!) that it's unreasonable to try for
# better than 1/10th of a pixel precision using default settings.
test_nirspec_source_offset_00 = lambda : do_test_source_offset('NIRSpec', theta=0.0, tolerance=0.1, monochromatic=3e-6)
test_nirspec_source_offset_45 = lambda : do_test_source_offset('NIRSpec', theta=45.0, tolerance=0.1, monochromatic=3e-6)

test_nirspec_set_siaf = lambda : do_test_set_position_from_siaf('NIRSpec')

def test_nirspec_slit_apertures():
    """Test that we can use slit and aperture names that don't map to a specific detector
    Verify that the V2 and V3 coordinates are reported as expected.
    """
    nrs = webbpsf_core.NIRSpec()

    for apname in [ 'NRS_FULL_IFU', 'NRS_S200A1_SLIT']:
        nrs.set_position_from_aperture_name(apname)

        ap = pysiaf.Siaf('NIRSpec')[apname]

        assert np.isclose(nrs._tel_coords()[0].to_value(u.arcsec),ap.V2Ref)
        assert np.isclose(nrs._tel_coords()[1].to_value(u.arcsec),ap.V3Ref)

