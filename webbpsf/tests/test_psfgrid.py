import astropy.convolution
from astropy.io import fits
import numpy as np
import pytest

from .. import gridded_library
from .. import webbpsf_core


# @pytest.mark.skip()
def test_compare_to_calc_psf_oversampled():
    """Check that the output PSF matches calc_psf and is saved in the correct location:
    for a distorted, oversampled case"""
    oversample = 2
    fov_pixels = 11

    fgs = webbpsf_core.FGS()
    fgs.detector = "FGS1"
    grid = fgs.psf_grid(all_detectors=False, num_psfs=4, oversample=oversample, fov_pixels=fov_pixels)

    psfnum = 1
    loc = grid[0].header["DET_YX{}".format(psfnum)]
    locy = int(loc.split()[0][1:-1])
    locx = int(loc.split()[1][:-1])
    gridpsf = grid[0].data[psfnum, :, :]

    fgs.detector_position = (locx, locy)
    calcpsf = fgs.calc_psf(oversample=oversample, fov_pixels=fov_pixels)["OVERDIST"].data
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    convpsf = astropy.convolution.convolve(calcpsf, kernel)

    assert gridpsf.shape == calcpsf.shape
    assert np.array_equal(gridpsf, convpsf)


# @pytest.mark.skip()
def test_comapre_to_calc_psf_detsampled():
    """Check that the output PSF matches calc_psf and is saved in the correct location:
    for an un-distorted, detector sampled case"""
    oversample = 2
    fov_pixels = 11

    nis = webbpsf_core.NIRISS()
    nis.filter = "F090W"
    nis.detector = "NIS"
    grid = nis.psf_grid(all_detectors=False, num_psfs=4, use_detsampled_psf=True, add_distortion=False,
                       oversample=oversample, fov_pixels=fov_pixels)

    psfnum = 1
    loc = grid[0].header["DET_YX{}".format(psfnum)]
    locy = int(loc.split()[0][1:-1])
    locx = int(loc.split()[1][:-1])
    gridpsf = grid[0].data[psfnum, :, :]

    nis.detector_position = (locx, locy)
    nis.options['output_mode'] = 'Detector Sampled Image'
    calcpsf = nis.calc_psf(oversample=oversample, fov_pixels=fov_pixels)["DET_SAMP"].data
    kernel = astropy.convolution.Box2DKernel(width=1)
    convpsf = astropy.convolution.convolve(calcpsf, kernel)

    assert gridpsf.shape == calcpsf.shape
    assert np.array_equal(gridpsf, convpsf)


# @pytest.mark.skip()
def test_setting_values():
    """Test the different ways to set filters and detectors"""
    oversample = 2
    fov_pixels = 1

    mir = webbpsf_core.MIRI()
    mir.filter = "F560W"

    # Method 1
    grid1 = mir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=fov_pixels, oversample=oversample)

    # Method 2
    mir.filter = "F560W"
    mir.detector = "MIRIM"
    grid2 = mir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=fov_pixels, oversample=oversample)

    # Check they are the same
    assert np.array_equal(grid1[0].data, grid2[0].data)


# @pytest.mark.skip()
def test_all():
    """Check that running all the detectors works"""
    fgs = webbpsf_core.FGS()
    grid = fgs.psf_grid(all_detectors=True, num_psfs=4, fov_pixels=10, oversample=2)

    # Check they shape
    assert len(grid) == 2
    assert grid[0][0].data.shape == (4, 20, 20)


# @pytest.mark.skip()
def test_one_psf():
    """Check that setting num_psfs = 1 produces the PSF in the right location"""
    oversample = 2
    fov_pixels = 11

    nis = webbpsf_core.NIRISS()
    nis.filter = "F140M"
    grid1 = nis.psf_grid(all_detectors=False, num_psfs=1, add_distortion=True, oversample=oversample,
                         fov_pixels=fov_pixels, single_psf_centered=True, use_detsampled_psf=False)

    nis.detector_position = (10, 0)  # it's set as (x,y)
    grid2 = nis.psf_grid(all_detectors=False, num_psfs=1, add_distortion=True, oversample=oversample,
                         fov_pixels=fov_pixels, single_psf_centered=False, use_detsampled_psf=False)

    calc = nis.calc_psf(add_distortion=True, oversample=2, fov_pixels=11)
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    conv = astropy.convolution.convolve(calc["OVERDIST"].data, kernel)

    assert grid1[0].header["DET_YX0"] == "(1024, 1024)"  # the default is the center of the NIS aperture
    assert grid2[0].header["DET_YX0"] == "(0, 10)"  # it's in (y,x)
    assert np.array_equal(conv, grid2[0].data[0, :, :])


# @pytest.mark.skip()
def test_nircam_errors():
    """Check that NIRCam has checks for incorrect value setting"""
    longfilt = "F250M"
    shortfilt = "F140M"
    longdet = "NRCB5"
    shortdet = "NRCA3"

    nir = webbpsf_core.NIRCam()

    # Shouldn't error
    nir.filter = longfilt
    nir.detector = longdet
    nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=1)  # no error

    nir.filter = shortfilt
    nir.detector = shortdet
    nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=1)  # no error

    # Should error
    with pytest.raises(ValueError) as excinfo:
        nir.filter =longfilt
        nir.detector =shortdet
        nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=1)  # error
    assert "ValueError" in str(excinfo)

    with pytest.raises(ValueError) as excinfo:
        nir.filter =shortfilt
        nir.detector =longdet
        nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=1)  # error
    assert "ValueError" in str(excinfo)

