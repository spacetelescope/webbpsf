import os

import astropy.convolution
from astropy.io import fits
import numpy as np
import pytest

from .. import gridded_library
from .. import webbpsf_core
from .. import roman
from .. import utils
from .. import detectors


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
    nlambda = 1

    # Create PSF grid
    fgs = webbpsf_core.FGS()
    fgs.detector = 'FGS1'
    grid = fgs.psf_grid(
        all_detectors=False, num_psfs=4, oversample=oversample, fov_pixels=fov_pixels, nlambda=nlambda, verbose=False
    )

    # Pull one of the PSFs out of the grid
    psfnum = 1
    loc = grid.grid_xypos[psfnum]
    locy = int(float(loc[1]) - 0.5)
    locx = int(float(loc[0]) - 0.5)
    gridpsf = grid.data[psfnum, :, :]

    # Using meta data, create the expected same PSF via calc_psf
    fgs.detector_position = (locx, locy)
    calcpsf = fgs.calc_psf(oversample=oversample, fov_pixels=fov_pixels, nlambda=nlambda)['OVERDIST'].data
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    convpsf = astropy.convolution.convolve(calcpsf, kernel)
    scalefactor = oversample**2  # normalization as used internally in GriddedPSFModel; see #302

    # Compare to make sure they are in fact the same PSF
    assert gridpsf.shape == calcpsf.shape, 'Shape mismatch'
    assert np.allclose(gridpsf, convpsf * scalefactor), 'Data values not as expected'


def test_compare_to_calc_psf_detsampled():
    """
    Check that the output PSF matches calc_psf and is saved in the correct slice of the array:
    for an un-distorted, detector sampled case
    """
    oversample = 2
    fov_arcsec = 0.5
    nlambda = 1

    # Create PSF grid
    mir = webbpsf_core.MIRI()
    mir.filter = 'F560W'
    mir.detector = 'MIRIM'
    grid = mir.psf_grid(
        all_detectors=False,
        num_psfs=4,
        use_detsampled_psf=True,
        add_distortion=False,
        oversample=oversample,
        fov_arcsec=fov_arcsec,
        nlambda=nlambda,
        verbose=False,
    )

    # Pull one of the PSFs out of the grid
    psfnum = 1
    loc = grid.grid_xypos[psfnum]
    locy = int(float(loc[1]))
    locx = int(float(loc[0]))
    gridpsf = grid.data[psfnum, :, :]

    # Using meta data, create the expected same PSF via calc_psf
    mir.detector_position = (locx, locy)
    mir.options['output_mode'] = 'Detector Sampled Image'
    calcpsf = mir.calc_psf(oversample=oversample, fov_arcsec=fov_arcsec, nlambda=nlambda)['DET_SAMP'].data
    kernel = astropy.convolution.Box2DKernel(width=1)
    convpsf = astropy.convolution.convolve(calcpsf, kernel)

    # Compare to make sure they are in fact the same PSF
    assert gridpsf.shape == calcpsf.shape
    assert np.array_equal(gridpsf, convpsf)


def test_all_detectors():
    """
    Check that running all the detectors works (ie setting all_detectors=True). In
    particular for NIRCam, test that the detectors pulled are correct
    (shortwave vs longwave) with respect to the filter
    """
    nir = webbpsf_core.NIRCam()
    longfilt = 'F250M'
    shortfilt = 'F140M'

    # Case 1: Shortwave -> check that only the SW detectors are applied for the SW filter
    nir.filter = shortfilt
    grid1 = nir.psf_grid(all_detectors=True, num_psfs=1, add_distortion=False, fov_pixels=4, oversample=2, verbose=False)
    det_list = []
    for hdu in grid1:
        det_list.append(hdu.meta['detector'][0])

    assert len(grid1) == len(gridded_library.CreatePSFLibrary.nrca_short_detectors)
    assert set(det_list) == set(gridded_library.CreatePSFLibrary.nrca_short_detectors)

    # Case 2: Longwave -> check that only the LW detectors are applied for the LW filter
    nir.filter = longfilt
    grid2 = nir.psf_grid(all_detectors=True, num_psfs=1, add_distortion=False, fov_pixels=4, oversample=2, verbose=False)
    det_list = []
    for hdu in grid2:
        det_list.append(hdu.meta['detector'][0])

    assert len(grid2) == len(gridded_library.CreatePSFLibrary.nrca_long_detectors)
    assert set(det_list) == set(gridded_library.CreatePSFLibrary.nrca_long_detectors)


def test_one_psf():
    """Check that setting num_psfs = 1 produces the PSF in the right location"""
    oversample = 2
    fov_pixels = 11

    nis = webbpsf_core.NIRISS()
    nis.filter = 'F140M'

    # Case 1: The PSF is centered on the detector (with single_psf_centered=True)
    grid1 = nis.psf_grid(
        all_detectors=False,
        num_psfs=1,
        add_distortion=True,
        oversample=oversample,
        fov_pixels=fov_pixels,
        single_psf_centered=True,
        use_detsampled_psf=False,
        verbose=False,
    )

    # Case 2: The PSF is set to a specific position (with nis.detector_position = (10, 0))
    nis.detector_position = (10, 0)  # it's set as (x,y)
    grid2 = nis.psf_grid(
        all_detectors=False,
        num_psfs=1,
        add_distortion=True,
        oversample=oversample,
        fov_pixels=fov_pixels,
        single_psf_centered=False,
        use_detsampled_psf=False,
        verbose=False,
    )

    # Compare Case 2 to the calc_psf output to make sure it's placing the PSF in the right location
    calc = nis.calc_psf(add_distortion=True, oversample=2, fov_pixels=11)
    # first, add IPC to the oversampled data from case 2, then convolve with detector binning
    detectors.apply_detector_ipc(calc, extname='OVERDIST')
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    convpsf = astropy.convolution.convolve(calc['OVERDIST'].data, kernel)
    scalefactor = oversample**2  # normalization as used internally in GriddedPSFModel; see #302

    assert np.all(
        grid1.grid_xypos == [[1023, 1023]]
    ), 'Center position not as expected'  # the default is the center of the NIS aperture
    assert np.all(grid2.grid_xypos == [(10, 0)]), 'Corner position not as expected'  # it's in (x,y)
    # check for near-equality of the PSFs computed both ways,
    #  but ignore the outer few pixels rows and columns, for which boundary wrapping leads to imperfect equality
    #  depending on the order of the convolutions. Ignore 2 rows/cols on either side based on oversample value
    assert np.allclose((convpsf * scalefactor)[2:-2, 2:-2], grid2.data[0, 2:-2, 2:-2]), 'PSF data values not as expected'


def test_nircam_errors():
    """Check that there are checks for incorrect value setting - particularly with NIRCam"""
    longfilt = 'F250M'
    shortfilt = 'F140M'
    longdet = 'NRCB5'
    shortdet = 'NRCA3'

    nir = webbpsf_core.NIRCam()

    # Shouldn't error - applying SW to SW and LW to LW
    nir.filter = longfilt
    nir.detector = longdet
    nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=4, detector_oversample=2, fft_oversample=2, verbose=False)

    nir.filter = shortfilt
    nir.detector = shortdet
    nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=4, detector_oversample=2, fft_oversample=2, verbose=False)

    # Should error - Bad filter/detector combination (LW filt to SW det)
    with pytest.raises(RuntimeError) as excinfo:  # Errors inside calc_psf() call
        nir.filter = longfilt
        nir.detector = shortdet
        nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=4, verbose=False)  # error
    assert 'RuntimeError' in str(excinfo)

    # Should error - Bad filter/detector combination (SW filt to LW det)
    with pytest.raises(RuntimeError) as excinfo:  # Errors inside calc_psf() call
        nir.filter = shortfilt
        nir.detector = longdet
        nir.psf_grid(all_detectors=False, num_psfs=1, fov_pixels=4, verbose=False)  # error
    assert 'RuntimeError' in str(excinfo)

    # Should error - Bad num_psfs entry (must be a square number)
    with pytest.raises(ValueError) as excinfo:
        nir.psf_grid(all_detectors=False, num_psfs=2, fov_pixels=4, verbose=False)  # error
    assert 'ValueError' in str(excinfo)


def test_saving(tmpdir):
    """Test saving files works properly"""

    # Create a temp directory to place file in
    directory = str(tmpdir)
    file = 'test1.fits'

    # Test using default calc_psf values
    fgs = webbpsf_core.FGS()
    fgs.filter = 'FGS'
    fgs.detector = 'FGS2'
    grid = fgs.psf_grid(
        all_detectors=False, num_psfs=4, nlambda=1, save=True, outdir=directory, outfile=file, overwrite=True
    )

    # Check that the saved file matches the returned file (and thus that the save worked through properly)
    with fits.open(os.path.join(directory, file[:-5] + '_fgs2.fits')) as infile:
        # Check data
        assert np.array_equal(infile[0].data, grid.data)

        # Check meta data
        model = utils.to_griddedpsfmodel(infile)
        assert model.meta.keys() == grid.meta.keys()
        assert np.all(model.grid_xypos == grid.grid_xypos)
        assert model.meta['oversampling'] == grid.meta['oversampling']

    # Remove temporary directory
    tmpdir.remove()


def test_2d_to_griddedpsfmodel():
    """Test that utils.to_griddedpsfmodel function works for a 2D HDUList input"""

    # Set up example 2D fits image
    data = np.ones((10, 10))
    primaryhdu = fits.PrimaryHDU(data)
    primaryhdu.header['DET_YX0'] = ('(1024, 1024)', 'The #0 PSFs (y,x) detector pixel position')
    primaryhdu.header['OVERSAMP'] = (5, 'oversampling value')
    hdu = fits.HDUList(primaryhdu)

    # Test that nothing errors when writing a GriddedPSFModel object
    model = utils.to_griddedpsfmodel(hdu)

    # Check the basic keywords are there
    assert 'det_yx0' in model.meta
    assert 'grid_xypos' in model.meta
    assert 'oversamp' in model.meta
    assert 'oversampling' in model.meta


def test_wfi():
    """Test that the psf_grid method works for the WFI class"""

    oversample = 2
    fov_pixels = 10
    nlambda = 1

    # Create PSF grid
    wfi = roman.WFI()
    grid = wfi.psf_grid(
        all_detectors=False,
        add_distortion=False,
        num_psfs=4,
        fov_pixels=fov_pixels,
        oversample=oversample,
        nlambda=nlambda,
        verbose=False,
    )

    # Pull one of the PSFs out of the grid
    psfnum = 1
    loc = grid.grid_xypos[psfnum]
    locy = int(float(loc[1]) - 0.5)
    locx = int(float(loc[0]) - 0.5)
    gridpsf = grid.data[psfnum, :, :]

    # Using meta data, create the expected same PSF via calc_psf
    wfi.detector_position = (locx, locy)
    calcpsf = wfi.calc_psf(oversample=oversample, fov_pixels=fov_pixels, nlambda=nlambda)['OVERSAMP'].data
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    convpsf = astropy.convolution.convolve(calcpsf, kernel)
    scalefactor = oversample**2  # normalization as used internally in GriddedPSFModel; see #302

    # Compare to make sure they are in fact the same PSF
    assert gridpsf.shape == calcpsf.shape, 'Shape mismatch'
    assert np.allclose(gridpsf, convpsf * scalefactor), 'Data values not as expected'
