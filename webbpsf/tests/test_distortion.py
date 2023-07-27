import pytest
import numpy as np
from astropy.io import fits

from .. import distortion
from .. import webbpsf_core


# @pytest.mark.skip()
def test_apply_distortion_skew():
    """
    Test the PSF axis is skewed appropriately by the apply_distortion function.

    Check that a rectangle of 1s run through the FGS1 distortion will be skewed (since FGS1 is heavily skewed in a
    known way). We'll check that the top left corner of the rectangle is higher up than the top right corner of the
    rectangle by checking the indices where the 1s begin/end.

    """

    # Create a baseline PSF to have shape/header keywords correct
    fgs = webbpsf_core.FGS()
    fgs.detector = "FGS1"
    fgs.options["output_mode"] = "Oversampled image"
    psf = fgs.calc_psf(add_distortion=False)

    # Set up new extensions (from webbpsf_core.JWInstrument._calc_psf_format_output)
    n_exts = len(psf)
    for ext in np.arange(n_exts):
        hdu_new = fits.ImageHDU(psf[ext].data, psf[ext].header)  # these will be the PSFs that are edited
        psf.append(hdu_new)
        ext_new = ext + n_exts
        psf[ext_new].header["EXTNAME"] = psf[ext].header["EXTNAME"][0:4] + "DIST"  # change extension name

    # Run data through the distortion function
    psf_siaf = distortion.apply_distortion(psf)

    # Rebin data to get 3rd extension
    fgs.options["output_mode"] = "Both extensions"
    fgs.options["detector_oversample"] = psf[0].header['DET_SAMP']
    webbpsf_core.SpaceTelescopeInstrument._calc_psf_format_output(fgs, result=psf_siaf, options=fgs.options)

    # Test the slope of the rectangle
    for ext in [2, 3]:
        left = psf_siaf[ext].data[:, 0]  # isolate the far left column
        right = psf_siaf[ext].data[:, -1]  # isolate the far right column

        indexes_left = [i for i, x in enumerate(left) if x != 0.]  # find the indices of the rectangle in the left col
        indexes_right = [i for i, x in enumerate(right) if x != 0.]  # find the indices of the rectangle in the right

        top_of_left = np.min(indexes_left)  # find the index of the top left corner of the rectangle
        top_of_right = np.min(indexes_right)  # find the index of the top right corner of the rectangle

        # Assert that the top of left > top of right due to FGS1 skew
        assert top_of_left > top_of_right, "FGS PSF does not have expected skew after distortion application"


# @pytest.mark.skip()
def test_apply_distortion_pixel_scale():
    """
    Test the pixel scale is changed by the apply_distortion function.

    Create a fake data set that has rows of constant value, so row 0 is all 0s, row 1 is all 1s, etc. Then distort it
    via apply_distortion(), which will change both the pixel scale and skew the data. If there was no skew in the
    data, there'd be a constant x and y pixel scale change, which wouldn't really affect the x direction since the
    function would be blending the same values together, and it would affect the y direction, but by the same amount
    (since again, the function is blending the same values at each index (ie blending 0 and 1 across the entire row).

    So subtract the linear shape caused by the skew out of the row and then check that across a specific row,
    the newly pixel-scale distorted ("blended") values are approximately equal.

    Use FGS1 for its large pixel scale change in this test

    """

    # Create a baseline PSF to have shape/header keywords correct
    fgs = webbpsf_core.FGS()
    fgs.detector = "FGS1"
    fgs.options["output_mode"] = "Oversampled image"
    psf = fgs.calc_psf(add_distortion=False)

    # Replace data with a fake image of row values equal to the row number
    data = np.zeros_like(psf[0].data)
    ny, nx = data.shape
    for i in np.arange(ny):
        data[i, :] = i

    # Replace the data in the PSF with the fake image
    psf[0].data = data

    # Set up new extensions (from webbpsf_core.JWInstrument._calc_psf_format_output)
    n_exts = len(psf)
    for ext in np.arange(n_exts):
        hdu_new = fits.ImageHDU(psf[ext].data, psf[ext].header)  # these will be the PSFs that are edited
        psf.append(hdu_new)
        ext_new = ext + n_exts
        psf[ext_new].header["EXTNAME"] = psf[ext].header["EXTNAME"][0:4] + "DIST"  # change extension name

    # Run data through the distortion function
    psf_siaf = distortion.apply_distortion(psf)

    # Rebin data to get 3rd extension (DET_DIST)
    fgs.options["output_mode"] = "Both extensions"
    fgs.options["detector_oversample"] = psf[0].header['DET_SAMP']
    webbpsf_core.SpaceTelescopeInstrument._calc_psf_format_output(fgs, result=psf_siaf, options=fgs.options)

    # Test that the change caused by the pixel distortion is approximately constant along the row
    # Choosing to check the 20th row.
    i = 20
    ext = 3

    # Crop off the edges due to skew / rotation that brings in 0s from beyond edge of detector
    psf_arr = psf_siaf[ext].data[5:-5, 5:-5]
    ncol = psf_arr.shape[1]
    inds = np.arange(ncol)

    # Model the skew with a basic linear function
    slope, intercept = np.polyfit(inds, psf_arr[i, :], 1)
    linear = (slope * inds) + intercept

    # Create a new 1D array that's your 20th row with the linear skew subtracted out
    final = psf_arr[i, :] - linear

    # Check the difference between adjacent values is the same to 1 decimal place
    diff = final[:-1] - final[1:]
    assert pytest.approx(diff, abs=0.1) == 0, \
        "FGS PSF does not have expected pixel scale distortion for adjacent pixels"

    # Check that the difference between the first and last value is also the same to 1 decimal
    assert pytest.approx(final[-1], abs=0.1) == final[0], "FGS PSF does not have expected pixel scale distortion in the " \
                                                      "entire row"


# @pytest.mark.skip()
def test_apply_rotation_error():
    """ Test that the apply_rotation function raises an error for NIRSpec and MIRI PSFs """

    # Create a PSF
    for inst in [webbpsf_core.NIRSpec(), webbpsf_core.MIRI()]:
        psf = inst.calc_psf(nlambda=1)  # done for speed

        # Test that running this function will raise a ValueError
        with pytest.raises(ValueError) as excinfo:
            distortion.apply_rotation(psf)
        assert "ValueError" in str(excinfo), "NIRSpec & MIRI PSFs should not be able to run through apply_rotation"


def test_distortion_with_custom_pixscale():
    """ Verifies the distortion model works properly even if the pixel scale is changed to
    a nonstandard value for the calculation. This tests/verifies the fix in PR 669:
        https://github.com/spacetelescope/webbpsf/pull/669
    """

    miri = webbpsf_core.MIRI()
    miri.pixelscale = 0.061
    psf = miri.calc_psf(fov_arcsec=2)

    # A symptom of the prior bug was the total sum of a distorted PSF would be very
    # discrepant from the sum of the undistorted PSF. So verif that symptom is not the case:

    assert np.isclose(psf[0].data.sum(), psf[3].data.sum(), rtol=0.001)
    assert np.isclose(psf[1].data.sum(), psf[3].data.sum(), rtol=0.001)
