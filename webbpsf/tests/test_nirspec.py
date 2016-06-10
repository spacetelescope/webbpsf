from __future__ import division, print_function, absolute_import, unicode_literals
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

import logging
_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())

from .. import webbpsf_core
import poppy


#------------------    NIRSpec Tests    ----------------------------

from .test_webbpsf import generic_output_test, do_test_source_offset

test_nirspec= lambda : generic_output_test('NIRSpec')

# Use a larger than typical tolerance when testing NIRSpec offsets. The
# pixels are so undersampled (0.1 arcsec!) that it's unreasonable to try for
# better than 1/10th of a pixel precision using default settings.
test_nirspec_source_offset_00 = lambda : do_test_source_offset('NIRSpec', theta=0.0, tolerance=0.1, monochromatic=3e-6)
test_nirspec_source_offset_45 = lambda : do_test_source_offset('NIRSpec', theta=45.0, tolerance=0.1, monochromatic=3e-6)
