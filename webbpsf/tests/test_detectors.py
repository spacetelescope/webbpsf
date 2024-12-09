import pytest
import numpy as np
from astropy.io import fits

import poppy
import webbpsf.detectors as detectors
import webbpsf.webbpsf_core as webbpsf_core


# @pytest.mark.skip()
def test_apply_miri_scattering_error():
    """Test that the apply_miri_scattering function raises an error for non-MIRI PSFs"""

    # Create a PSF
    nir = webbpsf_core.NIRCam()
    psf = nir.calc_psf(nlambda=1, fov_pixels=5)

    # Test that running this function will raise a ValueError
    with pytest.raises(ValueError) as excinfo:
        detectors.apply_miri_scattering(psf)
    assert 'ValueError' in str(excinfo), 'Non-MIRI PSFs should not be able to run through apply_miri_scattering'


# @pytest.mark.skip()
def test_apply_miri_scattering():
    """
    Test that a cross shape is added by the apply_miri_scattering function.

    Find the difference between the input and output PSF and check that the only non-zero values are
    along where the cross-shape would lie: i.e. lined up with the image's center.
    """

    # Create a baseline PSF to have shape/header keywords correct
    mir = webbpsf_core.MIRI()
    mir.filter = 'F560W'  # this filter has a strong cross added
    mir.options['output_mode'] = 'Oversampled image'
    psf = mir.calc_psf(add_distortion=False, nlambda=1)

    # Set up new extensions (from webbpsf_core.JWInstrument._calc_psf_format_output)
    n_exts = len(psf)
    for ext in np.arange(n_exts):
        hdu_new = fits.ImageHDU(psf[ext].data, psf[ext].header)  # these will be the PSFs that are edited
        psf.append(hdu_new)
        ext_new = ext + n_exts
        psf[ext_new].header['EXTNAME'] = psf[ext].header['EXTNAME'][0:4] + 'DIST'  # change extension name

    # Run it through just the apply_miri_scattering function
    psf_cross = detectors.apply_miri_scattering(psf)

    # Rebin data to get 3rd extension
    mir.options['output_mode'] = 'Both extensions'
    mir.options['detector_oversample'] = 1
    webbpsf_core.SpaceTelescopeInstrument._calc_psf_format_output(mir, result=psf_cross, options=mir.options)

    # Test distortion function
    for ext in [2, 3]:
        # Find the difference between the before and after PSF
        diff = psf_cross[ext].data - psf_cross[ext - 2].data

        # Test that the 4 corners of the box contain very small (close to 0) values
        ylen, xlen = diff.shape

        # Choose the start/stop points for these squares (each will take up 1/10th of the total array per side)
        # Following the addition of the radial term to the cruciform model, we restrict to just the very corners
        # to avoid including much from that circular halo term.
        first = 0
        second = int(0.1 * xlen)
        third = int(0.9 * xlen)
        fourth = xlen - 1

        # Pull these squares out of the data
        square1 = diff[first:second, first:second]
        square2 = diff[first:second, third:fourth]
        square3 = diff[third:fourth, first:second]
        square4 = diff[third:fourth, third:fourth]

        # What value it is compared to depends on the sampling since the range varies by a factor of the oversampling
        if ext == 2:
            value = 5e-7
        else:
            value = 1.5e-6

        # Show that these corner squares contain very small values
        assert_statement = (
            "should have lower values because the scattering shouldn't be adding much to this region."
            " It's too far away from where the cross is"
        )
        assert np.all(square1 < value), 'The LLCorner of the array {}'.format(assert_statement)
        assert np.all(square2 < value), 'The LRCorner of the array {}'.format(assert_statement)
        assert np.all(square3 < value), 'The ULCorner of the array {}'.format(assert_statement)
        assert np.all(square4 < value), 'The URCorner of the array {}'.format(assert_statement)

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
        assert avg_cross > avg_edge, 'The avg value of the cross should be larger than the avg value of the surrounding'
        assert avg_cross / 100 > avg_edge, (
            'The avg value of the cross should be larger than the avg value of the ' 'surrounding by a factor of 100'
        )


# @pytest.mark.skip()
def test_miri_conservation_energy():
    """
    Test that the miri scattering function follows conservation of energy, at least within a certain tolerance

    Compare the total sum of the pixels in a MIRI PSF before and after this scattering effect if applied. The PSFs
    should have almost the same intensity (almost, because there may be some light scattered off the detector at the
    edges in some cases, so we will add in a tolerance of 0.005).
    """

    # Create a baseline PSF to have shape/header keywords correct
    mir = webbpsf_core.MIRI()
    mir.filter = 'F1000W'
    mir.options['output_mode'] = 'Oversampled image'
    psf = mir.calc_psf(add_distortion=False, nlambda=1)

    # Set up new extensions (from webbpsf_core.JWInstrument._calc_psf_format_output)
    n_exts = len(psf)
    for ext in np.arange(n_exts):
        hdu_new = fits.ImageHDU(psf[ext].data, psf[ext].header)  # these will be the PSFs that are edited
        psf.append(hdu_new)
        ext_new = ext + n_exts
        psf[ext_new].header['EXTNAME'] = psf[ext].header['EXTNAME'][0:4] + 'DIST'  # change extension name

    # Run it through just the apply_miri_scattering function
    psf_cross = detectors.apply_miri_scattering(psf)

    # Rebin data to get 3rd extension
    mir.options['output_mode'] = 'Both extensions'
    mir.options['detector_oversample'] = 1
    webbpsf_core.SpaceTelescopeInstrument._calc_psf_format_output(mir, result=psf_cross, options=mir.options)

    # Test distortion function
    for ext in [2, 3]:
        psf_sum = np.sum(psf_cross[ext - 2].data.flatten())
        psf_cross_sum = np.sum(psf_cross[ext].data.flatten())

        assert pytest.approx(psf_sum, 0.005) == psf_cross_sum, (
            'The energy conversation of the PSF before/after the '
            'scattering is added is greater than the tolerance of '
            '0.005'
        )


def test_ipc_oversampling_equivalence(oversamp=2):
    """Test that we can apply in either order the IPC model and binning to detector pixel scale,
    and get the same results independent of order of operations.

    This is necessary to verify the "intuitive" way of applying IPC to detector-sampled data,
    and the alternative way to apply it to higher resolution oversampled data, are equivalent.
    """
    nrc = webbpsf_core.NIRCam()

    testpsf = nrc.calc_psf(nlambda=1, oversample=oversamp, fov_pixels=5)

    # regular version, with IPC added after binning to detector sampling
    # this happens in normal calc_psf calls
    psf_detdist = testpsf['DET_DIST'].data.copy()  # Binned then has IPC added

    # apply IPC to oversampled extension, then bin
    # this happens in psf_grid calls
    detectors.apply_detector_ipc(testpsf, extname='OVERDIST')
    psf_detdist_v2 = poppy.utils.rebin_array(testpsf['OVERDIST'].data, (oversamp, oversamp))

    assert np.allclose(
        psf_detdist, psf_detdist_v2
    ), 'PSFs calculated should be equivalent for IPC convolution and binning in either order'


def test_ipc_basic_effect_on_psf_fwhm():
    """A basic test that the IPC model has the expected effect: making PSFs slightly broader.

    Tests that (a) there is no change to the first two extensions, which are 'pure optical PSF'
               (b) that the FWHM increases for the other two extensions, which are distortion+detector effects
    """
    nrc = webbpsf_core.NIRCam()
    psf_withipc = nrc.calc_psf(nlambda=1, fov_pixels=101)
    nrc.options['add_ipc'] = False
    psf_noipc = nrc.calc_psf(nlambda=1, fov_pixels=101)

    for extname in ['OVERSAMP', 'DET_SAMP']:
        fwhm_ipc = poppy.measure_fwhm(psf_withipc, ext=extname)
        fwhm_noipc = poppy.measure_fwhm(psf_noipc, ext=extname)
        assert np.isclose(fwhm_ipc, fwhm_noipc), f'Adding IPC should not have any effect on the {extname} data.'
        print(f'test ok for {extname}')

    for extname in ['OVERDIST', 'DET_DIST']:
        fwhm_ipc = poppy.measure_fwhm(psf_withipc, ext=extname)
        fwhm_noipc = poppy.measure_fwhm(psf_noipc, ext=extname)
        assert fwhm_ipc > fwhm_noipc, f'Adding IPC should not blur and increase the FWHM in the {extname} data.'
        print(f'test ok for {extname}')
