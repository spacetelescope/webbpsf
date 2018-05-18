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

#------------------    NIRISS Tests    ----------------------------

from .test_webbpsf import generic_output_test, do_test_source_offset, do_test_set_position_from_siaf
test_niriss= lambda : generic_output_test('NIRISS')
test_niriss_source_offset_00 = lambda : do_test_source_offset('NIRISS', theta=0.0, monochromatic=3.0e-6)
test_niriss_source_offset_45 = lambda : do_test_source_offset('NIRISS', theta=45.0, monochromatic=3.0e-6)

test_niriss_set_siaf = lambda : do_test_set_position_from_siaf('NIRISS', 
        ['NIS_FP1MIMF', 'NIS_SUB64', 'NIS_SOSSFULL','NIS_SOSSTA','NIS_AMI1'])

def test_niriss_auto_pupil():
    """ Test switching between CLEAR and CLEARP
    depending on selected filter or wavelengths
    """

    niriss = webbpsf_core.NIRISS()
    assert niriss.pupil_mask is None

    niriss.filter='F277W'
    niriss.calc_psf(nlambda=1)
    assert niriss.pupil_mask == 'CLEARP'

    niriss.filter='F090W'
    niriss.calc_psf(nlambda=1)
    assert niriss.pupil_mask is None

    niriss.filter='F480M'
    niriss.calc_psf(nlambda=1)
    assert niriss.pupil_mask == 'CLEARP'

    niriss.filter='F200W'
    niriss.calc_psf(nlambda=1)
    assert niriss.pupil_mask is None

def test_niriss_gr700xd():
    '''
    Smoke-test calculations with the GR700XD custom optic
    present in the system. This is a regression test for
    https://github.com/mperrin/webbpsf/issues/148
    '''
    niriss = webbpsf_core.NIRISS()
    niriss.filter = 'CLEAR'
    niriss.pupil_mask = 'GR700XD'
    niriss.calc_psf(monochromatic=1e-6, fov_pixels=2)
