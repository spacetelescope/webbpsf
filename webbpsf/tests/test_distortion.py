import pytest
import numpy as np

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
    psf = fgs.calc_psf()

    # Re-write PSF with fake data: An array of 0s with a horizontal rectangle of 1s spanning the entire x axis
    for ext in np.arange(4):
        y = psf[ext].data.shape[0]
        x = psf[ext].data.shape[1]
        psf[ext].data = np.zeros((y, x))
        psf[ext].data[10:30, :] = np.ones((20, x))

    # Run data through the distortion function
    psf_siaf = distortion.apply_distortion(psf)

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
    psf = fgs.calc_psf()

    # Re-write PSF with fake data: An array with rows of constant value (a row of 0s, row of 1s, row of 2s, etc)
    for ext in np.arange(4):
        size = psf[ext].data.shape[0]
        arr = np.zeros((size, size))
        for i in range(size):
            arr[i] = np.full((1, size), i)
        psf[ext].data = arr

    # Run data through the distortion function
    psf_siaf = distortion.apply_distortion(psf)

    # Test that the change caused by the pixel distortion is approximately constant along the row
    # Choosing to check the 20th row.
    i = 20
    ext = 3

    psf_arr = psf_siaf[ext].data
    size = psf_siaf[ext].data.shape[0]
    inds = np.arange(len(psf_arr))

    # Model the skew with a basic linear function
    yN = psf_arr[i, -1]
    y0 = psf_arr[i, 0]
    slope = (yN - y0) / len(psf_arr)
    linear = (slope * inds) + y0

    # Create a new 1D array that's your 20th row with the linear skew subtracted out
    # Add y0 because we want the values compared to 0th value, not subtracted down to 0 (just preference)
    final = psf_arr[i, :] - linear + y0

    # Check the difference between adjacent values is the same to 1 decimal place
    for i in range(size - 1):
        a = final[i]
        b = final[i + 1]

        # This is the same as assert round(a - b, 1) == 0
        assert pytest.approx(a, 0.1) == b, "FGS PSF does not have expected pixel scale distortion for adjacent pixels"

    # Check that the difference between the first and last value is also the same to 1 decimal
    assert pytest.approx(final[-1], 0.1) == final[0], "FGS PSF does not have expected pixel scale distortion in the " \
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


# @pytest.mark.skip()
def test_apply_miri_scattering_error():
    """ Test that the apply_miri_scattering function raises an error for non-MIRI PSFs """

    # Create a PSF
    nir = webbpsf_core.NIRCam()
    psf = nir.calc_psf()

    # Test that running this function will raise a ValueError
    with pytest.raises(ValueError) as excinfo:
        distortion.apply_miri_scattering(psf)
    assert "ValueError" in str(excinfo), "Non-MIRI PSFs should not be able to run through apply_miri_scattering"


# @pytest.mark.skip()
def test_apply_miri_scattering():
    """
    Test that a cross shape is added by the apply_miri_scattering function.

    Find the difference between the input and output PSF and check that the only non-zero values are
    along where the cross-shape would lie: i.e. lined up with the image's center.
    """

    # Create a PSF
    mir = webbpsf_core.MIRI()
    mir.filter = "F560W"  # this filter has a strong cross added
    psf = mir.calc_psf()

    # Because calc_psf automatically applies distortions to ext 2 and 3, we'll overwrite these with the undistorted PSFs
    psf[2].data = psf[0].data
    psf[3].data = psf[1].data

    # Run it through just the apply_miri_scattering function
    psf_cross = distortion.apply_miri_scattering(psf)

    for ext in [2, 3]:
        # Find the difference between the before and after PSF
        diff = psf_cross[ext].data - psf[ext].data

        # Test that the 4 corners of the box contain very small (close to 0) values
        ylen, xlen = diff.shape

        # Choose the start/stop points for these squares (each will take up 1/3 of the total array)
        first = 0
        second = int(0.33 * xlen)
        third = int(0.67 * xlen)
        fourth = xlen-1

        # Pull these squares out of the data
        square1 = diff[first:second, first:second]
        square2 = diff[first:second, third:fourth]
        square3 = diff[third:fourth, first:second]
        square4 = diff[third:fourth, third:fourth]

        # What value it is compared to depends on the sampling since the range varies by a factor of the oversampling
        if ext == 2:
            value = 5e-7
        else:
            value = 1e-6

        # Show that these corner squares contain very small values
        assert_statement = "should have lower values because the scattering shouldn't be adding much to this region." \
                           " It's too far away from where the cross is"
        assert np.all(square1 < value), "The LLCorner of the array {}".format(assert_statement)
        assert np.all(square2 < value), "The LRCorner of the array {}".format(assert_statement)
        assert np.all(square3 < value), "The ULCorner of the array {}".format(assert_statement)
        assert np.all(square4 < value), "The URCorner of the array {}".format(assert_statement)

        # Test that there is a cross in the box which has a higher value than the surrounding area
        xcen = int(xlen / 2)
        ycen = int(ylen / 2)

        # Pull 20 values along the cross in both the x and y direction to check
        # shift up 20 pixels to ignore 0s near center
        cross_values_list = []
        for i in np.arange(20) + 20:
            cross_values_list.append(diff[xcen + i, ycen])
            cross_values_list.append(diff[xcen, ycen + i])

        # Find the average value of the points on the cross and squares
        avg_cross = np.mean(cross_values_list)
        avg_edge = np.mean([square1, square2, square3, square4])

        # Show that the average cross value is greater than the average square value by a factor of >100
        assert avg_cross > avg_edge, "The avg value of the cross should be larger than the avg value of the surrounding"
        assert avg_cross / 100 > avg_edge, "The avg value of the cross should be larger than the avg value of the " \
                                           "surrounding by a factor of 100"


# @pytest.mark.skip()
def test_miri_conservation_energy():
    """
    Test that the miri scattering function follows conservation of energy, at least within a certain tolerance

    Compare the total sum of the pixels in a MIRI PSF before and after this scattering effect if applied. The PSFs
    should have almost the same intensity (almost, because there may be some light scattered off the detector at the
    edges in some cases, so we will add in a tolerance of 0.005).
    """

    # Create a PSF
    mir = webbpsf_core.MIRI()
    mir.filter = "F1000W"
    psf = mir.calc_psf()

    # Because calc_psf automatically applies distortions to ext 2 and 3, we'll overwrite these with the undistorted PSFs
    psf[2].data = psf[0].data
    psf[3].data = psf[1].data

    # Run it through just the apply_miri_scattering function
    psf_cross = distortion.apply_miri_scattering(psf)

    for ext in [2, 3]:
        psf_sum = np.sum(psf[ext].data.flatten())
        psf_cross_sum = np.sum(psf_cross[ext].data.flatten())

        assert pytest.approx(psf_sum, 0.005) == psf_cross_sum, "The energy conversation of the PSF before/after the " \
                                                               "scattering is added is greater than the tolerance of " \
                                                               "0.005"
