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

from test_webbpsf import generic_output_test, do_test_source_offset

test_miri= lambda : generic_output_test('MIRI')
test_miri_00 = lambda : do_test_source_offset('MIRI', theta=0.0)
test_miri_45 = lambda : do_test_source_offset('MIRI', theta=45.0)


def test_miri_fqpm(theta=0.0, nsteps=3, nlambda=1, clobber=True, outputdir=None, display=False):
    miri = webbpsf_core.MIRI()
    miri.pupilopd = None
    miri.filter='F1065C'
    miri.image_mask = 'FQPM1065'
    miri.pupil_mask = 'MASKFQPM'
    
    oversample=2

    if outputdir is None:
        import tempfile
        outputdir = tempfile.gettempdir()

    for offset in np.linspace(0.0, 1.0, nsteps):
        miri.options['source_offset_theta'] = 0.0
        miri.options['source_offset_r'] = offset

        for angle in [0,45]:
            fn = os.path.join(outputdir, 'test_miri_fqpm_t{0}_r{1:.2f}.fits'.format(angle,offset))
            if not os.path.exists(fn) or clobber:
                miri.options['source_offset_theta'] = angle 
                psf = miri.calcPSF(oversample=oversample, nlambda=nlambda, save_intermediates=False, display=display)
                psf.writeto(fn, clobber=clobber)

    #FIXME - add some assertion tests here. 


