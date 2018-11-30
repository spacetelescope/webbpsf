import os

import astropy.convolution
from astropy.io import fits
import numpy as np
import pytest

from .. import gridded_library
from .. import webbpsf_core
from .. import utils


def test_compare_to_calc_psf_oversampled():
    """
    Check that the output PSF matches calc_psf and is saved in the correct slice of the array:
    for a distorted, oversampled case

    This case also uses an even length array, so we'll need to subtract 0.5 from the detector
    position because psf_grid value in the meta data has been shifted during calc_psf to account
    for it being an even length array and this shift shouldn't happen 2x (ie again in calc_psf
    call below)
    """
    oversample = 2
    fov_pixels = 10

    # Create PSF grid
    fgs = webbpsf_core.FGS()
    fgs.detector = "FGS1"
    grid = fgs.psf_grid(all_detectors=False, num_psfs=4, oversample=oversample, fov_pixels=fov_pixels)

    # Pull one of the PSFs out of the grid
    psfnum = 1
    loc = grid.meta["grid_xypos"][psfnum]
    locy = int(float(loc[1]) - 0.5)
    locx = int(float(loc[0]) - 0.5)
    gridpsf = grid.data[psfnum, :, :]

    # Using meta data, create the expected same PSF via calc_psf
    fgs.detector_position = (locx, locy)
    calcpsf = fgs.calc_psf(oversample=oversample, fov_pixels=fov_pixels)["OVERDIST"].data
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    convpsf = astropy.convolution.convolve(calcpsf, kernel)

    # Compare to make sure they are in fact the same PSF
    assert gridpsf.shape == calcpsf.shape
    assert np.array_equal(gridpsf, convpsf)


def test_compare_to_calc_psf_detsampled():
    """
    Check that the output PSF matches calc_psf and is saved in the correct slice of the array:
    for an un-distorted, detector sampled case
    """
    oversample = 2
    fov_arcsec = 0.5

    # Create PSF grid
    mir = webbpsf_core.MIRI()
    mir.filter = "F560W"
    mir.detector = "MIRIM"
    grid = mir.psf_grid(all_detectors=False, num_psfs=4, use_detsampled_psf=True, add_distortion=False,
                        oversample=oversample, fov_arcsec=fov_arcsec)

    # Pull one of the PSFs out of the grid
    psfnum = 1
    loc = grid.meta["grid_xypos"][psfnum]
    locy = int(float(loc[1]))
    locx = int(float(loc[0]))
    gridpsf = grid.data[psfnum, :, :]

    # Using meta data, create the expected same PSF via calc_psf
    mir.detector_position = (locx, locy)
    mir.options['output_mode'] = 'Detector Sampled Image'
    calcpsf = mir.calc_psf(oversample=oversample, fov_arcsec=fov_arcsec)["DET_SAMP"].data
    kernel = astropy.convolution.Box2DKernel(width=1)
    convpsf = astropy.convolution.convolve(calcpsf, kernel)

    # Compare to make sure they are in fact the same PSF
    assert gridpsf.shape == calcpsf.shape
    assert np.array_equal(gridpsf, convpsf)


def test_all():
    """
    Check that running all the detectors works (ie setting all_detectors=True). In
    particular for NIRCam, test that the detectors pulled are correct
    (shortwave vs longwave) with respect to the filter
    """
    nir = webbpsf_core.NIRCam()
    longfilt = "F250M"
    shortfilt = "F140M"

    # Case 1: Shortwave -> check that only the SW detectors are applied for the SW filter
    nir.filter = shortfilt
    grid1 = nir.psf_grid(all_detectors=True, num_psfs=1, add_distortion=False, fov_pixels=1, oversample=2)
    det_list = []
    for hdu in grid1:
        det_list.append(hdu.meta["detector"][0])

    assert len(grid1) == len(gridded_library.CreatePSFLibrary.nrca_short_detectors)
    assert set(det_list) == set(gridded_library.CreatePSFLibrary.nrca_short_detectors)

    # Case 2: Longwave -> check that only the LW detectors are applied for the LW filter
    nir.filter = longfilt
    grid2 = nir.psf_grid(all_detectors=True, num_psfs=1, add_distortion=False, fov_pixels=1, oversample=2)
    det_list = []
    for hdu in grid2:
        det_list.append(hdu.meta["detector"][0])

    assert len(grid2) == len(gridded_library.CreatePSFLibrary.nrca_long_detectors)
    assert set(det_list) == set(gridded_library.CreatePSFLibrary.nrca_long_detectors)


def test_one_psf():
    """Check that setting num_psfs = 1 produces the PSF in the right location"""
    oversample = 2
    fov_pixels = 11

    nis = webbpsf_core.NIRISS()
    nis.filter = "F140M"

    # Case 1: The PSF is centered on the detector (with single_psf_centered=True)
    grid1 = nis.psf_grid(all_detectors=False, num_psfs=1, add_distortion=True, oversample=oversample,
                         fov_pixels=fov_pixels, single_psf_centered=True, use_detsampled_psf=False)

    # Case 2: The PSF is set to a specific position (with nis.detector_position = (10, 0))
    nis.detector_position = (10, 0)  # it's set as (x,y)
    grid2 = nis.psf_grid(all_detectors=False, num_psfs=1, add_distortion=True, oversample=oversample,
                         fov_pixels=fov_pixels, single_psf_centered=False, use_detsampled_psf=False)

    # Compare Case 2 to the calc_psf output to make sure it's placing the PSF in the right location
    calc = nis.calc_psf(add_distortion=True, oversample=2, fov_pixels=11)
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    convpsf = astropy.convolution.convolve(calc["OVERDIST"].data, kernel)

    assert grid1.meta["grid_xypos"] == [(1023, 1023)]  # the default is the center of the NIS aperture
    assert grid2.meta["grid_xypos"] == [(10, 0)]  # it's in (x,y)
    assert np.array_equal(convpsf, grid2.data[0, :, :])


def test_nircam_errors():
    """Check that there are checks for incorrect value setting - particularly with NIRCam"""
    longfilt = "F250M"
    shortfilt = "F140M"
    longdet = "NRCB5"
    shortdet = "NRCA3"

    nir = webbpsf_core.NIRCam()

    # Shouldn't error - applying SW to SW and LW to LW
    nir.filter = longfilt
    nir.detector = longdet
    nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=1, detector_oversample=2, fft_oversample=2)

    nir.filter = shortfilt
    nir.detector = shortdet
    nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=1, detector_oversample=2, fft_oversample=2)

    # Should error - Bad filter/detector combination (LW filt to SW det)
    with pytest.raises(RuntimeError) as excinfo:  # Errors inside calc_psf() call
        nir.filter = longfilt
        nir.detector = shortdet
        nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=1)  # error
    assert "RuntimeError" in str(excinfo)

    # Should error - Bad filter/detector combination (SW filt to LW det)
    with pytest.raises(RuntimeError) as excinfo:  # Errors inside calc_psf() call
        nir.filter = shortfilt
        nir.detector = longdet
        nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=1)  # error
    assert "RuntimeError" in str(excinfo)

    # Should error - Bad num_psfs entry (must be a square number)
    with pytest.raises(ValueError) as excinfo:
        nir.psf_grid(all_detectors=False, num_psfs=2, fov_pixels=1)  # error
    assert "ValueError" in str(excinfo)


def test_saving(tmpdir):
    """Test saving files works properly"""

    # Create a temp directory to place file in
    file = str(tmpdir.join("test1"))

    # Test using default calc_psf values
    fgs = webbpsf_core.FGS()
    fgs.filter = "FGS"
    fgs.detector = "FGS2"
    grid = fgs.psf_grid(all_detectors=False, num_psfs=4, save=True, outfile=file, overwrite=True)

    # Check that the saved file matches the returned file (and thus that the save worked through properly)
    with fits.open(os.path.join(file[:-5], "test1_fgs2_fgs.fits")) as infile:
        # Check data
        assert np.array_equal(infile[0].data, grid.data)

        # Check meta data
        model = utils.to_griddedpsfmodel(infile)
        assert model.meta.keys() == grid.meta.keys()
        assert model.meta["grid_xypos"] == grid.meta["grid_xypos"]
        assert model.meta["oversampling"] == grid.meta["oversampling"]

    # Remove temporary directory
    tmpdir.remove()
