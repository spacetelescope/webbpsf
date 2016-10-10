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


#------------------    MIRI Tests    ----------------------------

from .test_webbpsf import generic_output_test, do_test_source_offset

test_miri= lambda : generic_output_test('MIRI')
test_miri_source_offset_00 = lambda : do_test_source_offset('MIRI', theta=0.0, monochromatic=8e-6)
test_miri_source_offset_45 = lambda : do_test_source_offset('MIRI', theta=45.0, monochromatic=8e-6)


def do_test_miri_fqpm(nlambda=1, clobber=True, angle=0.0, offset=0.0, oversample=2, outputdir=None, display=False, save=False):
    miri = webbpsf_core.MIRI()
    miri.pupilopd = None
    miri.filter='F1065C'
    miri.image_mask = 'FQPM1065'
    miri.pupil_mask = 'MASKFQPM'
 
    #for offset in np.linspace(0.0, 1.0, nsteps):
    #miri.options['source_offset_theta'] = 0.0
    miri.options['source_offset_r'] = offset

    #for angle in [0,45]:
    miri.options['source_offset_theta'] = angle 
    psf = miri.calc_psf(oversample=oversample, nlambda=nlambda, save_intermediates=False, display=display)

    if save:
        if outputdir is None:
            import tempfile
            outputdir = tempfile.gettempdir()


        fn = os.path.join(outputdir, 'test_miri_fqpm_t{0}_r{1:.2f}.fits'.format(angle,offset))
        psf.writeto(fn, clobber=clobber)

    #FIXME - add some assertion tests here. 

def test_miri_fqpm_centered(*args, **kwargs):
    do_test_miri_fqpm(angle=0.0, offset=0.0)


def test_miri_fqpm_offset_00(*args, **kwargs):
    do_test_miri_fqpm(angle=0.0, offset=1.0)

def test_miri_fqpm_offset_45(*args, **kwargs):
    do_test_miri_fqpm(angle=45.0, offset=1.0)


