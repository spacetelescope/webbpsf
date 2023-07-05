from .. import webbpsf_core

# ------------------    MIRI/MRS Tests    ----------------------------

#from .test_webbpsf import generic_output_test, do_test_source_offset, do_test_set_position_from_siaf

def test_monochromatic_psf():
    mrs = webbpsf_core.MRS()
    out = mrs.calc_psf(monochromatic=5* 1e-6, display=True)
    assert out is not None