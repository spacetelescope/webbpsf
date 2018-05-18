from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
import six
import copy

import pysiaf
import poppy

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

from astropy.table import Table
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.filters import uniform_filter
from decimal import Decimal, ROUND_HALF_UP
from astropy.modeling.functional_models import Gaussian2D


def _get_default_SIAF(instrument, aper_name):
    """ Store the default SIAF values for distortion and rotation """

    # Create new naming because SIAF requires special capitalization
    if instrument == "NIRCAM":
        siaf_name = "NIRCam"
    elif instrument == "NIRSPEC":
        siaf_name = "NIRSpec"
    else:
        siaf_name = instrument

    # Select a single SIAF aperture
    siaf = pysiaf.Siaf(siaf_name)
    aper = siaf.apertures[aper_name]

    return aper


def apply_distortion(HDUlist_or_filename=None, fill_value=0):
    """
    Apply a distortion to the input PSF. The distortion comes from the SIAF 4-5 degree polynomial (depends on the
    instrument). This function pulls and applies the SIAF polynomial values using pysiaf, which ensures the most
    up-to-date values will be called.

    Parameters
    ----------

    HDUlist_or_filename :
        A PSF from WebbPSF, either as an HDUlist object or as a filename
    fill_value : float
        Value used to fill in any blank space by the skewed PSF. Default = 0

    """

    # Read in input PSF
    if isinstance(HDUlist_or_filename, six.string_types):
        hdu_list = fits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, fits.HDUList):
        hdu_list = HDUlist_or_filename
    else:
        raise ValueError("input must be a filename or HDUlist")

    # Create a copy of the PSF
    psf = copy.deepcopy(hdu_list)

    # Log instrument and detector names
    instrument = hdu_list[0].header["INSTRUME"].upper()
    aper_name = hdu_list[0].header["APERNAME"].upper()

    # Pull default values
    aper = _get_default_SIAF(instrument, aper_name)

    ext = 2  # edit the oversampled PSF, then bin it down to get the detector sampled PSF

    # Pull PSF header information
    pixelscale = psf[ext].header["PIXELSCL"]  # the pixel scale carries the over-sample value
    oversamp = psf[ext].header["OVERSAMP"]  # will be 1 for ext=1
    xpix_center = psf[ext].header["DET_X"]  # center x location in pixels
    ypix_center = psf[ext].header["DET_Y"]  # center y location in pixels
    len_y = psf[ext].shape[0]
    len_x = psf[ext].shape[1]

    # Convert the PSF center point from pixels to arcseconds using pysiaf
    xarc_center, yarc_center = aper.sci_to_idl(xpix_center, ypix_center)

    # ###############################################
    # Create an array of indices (in pixels) for where the PSF is located on the detector
    # 1) Set up blank indices (in pixels)
    ypix, xpix = np.indices((len_y, len_x), dtype=float)

    # 2) Shift indices to be centered on (0,0) (starting to transform into the Ideal frame)
    ypix -= (len_y - 1.) / 2.
    xpix -= (len_x - 1.) / 2.

    # 3) Convert these indices from pixels to arcseconds
    # Note: This also shifts the oversampled indices so they span the same region as the detector-sampled indices
    # but the oversampled array is still longer by a factor of the oversample
    yarc = ypix * pixelscale
    xarc = xpix * pixelscale

    # 4) Shift the indices so they match where on the detector the PSF is located
    yidl = yarc + yarc_center
    xidl = xarc + xarc_center

    # 5) Now that the indices are in the Ideal frame, convert them to the Science Frame using idl_to_sci
    # Going from Idl to Sci this way allows us to add in the distortion
    xsci, ysci = aper.idl_to_sci(xidl, yidl)

    # ###############################################
    # Create an array of indices (in pixels) that the final data will be interpolated on to
    # 1) Set up blank indices (in pixels)
    ynew, xnew = np.indices([len_y, len_x], dtype=float)

    # 2) Shift indices to be in the Ideal frame (centered on 0)
    xnew -= (len_x - 1.) / 2.
    ynew -= (len_y - 1.) / 2.

    # 3) Shift the oversampled indices so they span the same region as the detector-sampled indices
    # Note: the oversampled array is still longer by a factor of the oversample
    xnew /= oversamp
    ynew /= oversamp

    # 4) Shift the indices so they match where on the detector the PSF is located
    xnew += xpix_center
    ynew += ypix_center

    # ###############################################
    # Interpolate from the original indices (xsci, ysci) on to new indices (xnew, ynew)
    psf_new = griddata((xsci.flatten(), ysci.flatten()), psf[ext].data.flatten(), (xnew, ynew),
                       fill_value=fill_value)

    # Apply data to correct extensions
    psf[ext].data = psf_new

    # Now bin down over-sampled PSF to be detector-sampled and re-write ext=3
    detector_oversample = psf[ext].header["DET_SAMP"]
    psf[3].data = poppy.utils.rebin_array(psf_new, rc=(detector_oversample, detector_oversample))

    # Set new header keywords
    for ext in [2, 3]:
        psf[ext].header["DISTORT"] = ("True", "SIAF distortion coefficients applied")
        psf[ext].header["SIAF_VER"] = (pysiaf.JWST_PRD_VERSION, "SIAF PRD version used")

        degree = np.int(getattr(aper, 'Sci2IdlDeg'))
        number_of_coefficients = np.int((degree + 1) * (degree + 2) / 2)
        all_keys = aper.__dict__.keys()
        for axis in ['X', 'Y']:
            coeff_keys = np.sort(np.array([c for c in all_keys if 'Idl2Sci' + axis in c]))
            coeff = np.array([getattr(aper, c) for c in coeff_keys[0:number_of_coefficients]])
            for i in range(len(coeff)):
                key = "COEF_{}".format(coeff_keys[i][-3:])
                psf[ext].header[key] = (coeff[i], "SIAF distortion coefficient for {}".format(coeff_keys[i]))

    return psf

# #####################################################################################################################


def apply_rotation(HDUlist_or_filename=None, rotate_value=None, crop=True):
    """
    Apply the detector's rotation to the PSF. This is for NIRCam, NIRISS, and FGS. MIRI and NIRSpec's large rotation is
    already added inside WebbPSF's calculations.

    Parameters
    ----------

    HDUlist_or_filename :
        A PSF from WebbPSF, either as an HDUlist object or as a filename
    rotate_value : float
        Rotation in degrees that PSF needs to be. If set to None, function will pull the most up to date
        SIAF value. Default = None.
    crop : bool
        True or False to crop the PSF so it matches the size of the input PSF (e.g. so they could be more easily
        compared).

    """
    # Read in input PSF
    if isinstance(HDUlist_or_filename, six.string_types):
        hdu_list = fits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, fits.HDUList):
        hdu_list = HDUlist_or_filename
    else:
        raise ValueError("input must be a filename or HDUlist")

    # Create a copy of the PSF
    psf = copy.deepcopy(hdu_list)

    # Log instrument and detector names
    instrument = hdu_list[0].header["INSTRUME"].upper()
    aper_name = hdu_list[0].header["APERNAME"].upper()

    if instrument == "MIRI":
        raise ValueError("MIRI's rotation is already included in WebbPSF and shouldn't be added again.")

    if instrument == "NIRSPEC":
        raise ValueError("NIRSpec's rotation is already included in WebbPSF and shouldn't be added again.")

    # Set rotation value if not already set by a keyword argument
    if rotate_value is None:
        aper = _get_default_SIAF(instrument, aper_name)
        rotate_value = getattr(aper, "V3IdlYAngle")  # the angle to rotate the PSF in degrees

    # If crop = True, then reshape must be False - so invert this keyword
    reshape = np.invert(crop)

    ext = 2  # edit the oversampled PSF, then bin it down to get the detector sampled PSF

    psf_new = rotate(psf[ext].data, rotate_value, reshape=reshape)

    # Apply data to correct extensions
    psf[ext].data = psf_new

    # Now bin down over-sampled PSF to be detector-sampled and re-write ext=3
    detector_oversample = psf[ext].header["DET_SAMP"]
    psf[3].data = poppy.utils.rebin_array(psf_new, rc=(detector_oversample, detector_oversample))

    # Set new header keyword
    for ext in [2, 3]:
        psf[ext].header["ROTATION"] = (rotate_value, "PSF rotated to match detector rotation")

    return psf

# #####################################################################################################################


def _get_default_miri(filter):
    """
    Store the default values for the MIRI scattering cross artifact distortion transformation. Values come from
    MIRI-TN-00076-ATC_Imager_PSF_Issue_4.pdf
    """

    radius = 200  # radius of kernel profile (in MIRI detector pixels)

    aper = pysiaf.Siaf("MIRI").apertures["MIRIM_FULL"]
    rotate_value = getattr(aper, "V3IdlYAngle")  # = 4.4497  # rotation value pulled from most updated SIAF file

    filter_list = ['F560W', 'F770W', 'F1000W', 'F1130W', 'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W',
                   'FND', 'F1065C', 'F1140C', 'F1550C', 'F2300C']

    kernel_amp_list = [0.00220, 0.00139, 0.00034, 0.00007, 0.00011, 0.0, 0.0, 0.0, 0.0,
                       0.00087, 0.00010, 0.00007, 0.0, 0.0]  # detector scattering kernel amplitude

    # Set PSF values
    i_filter = filter_list.index(filter)
    kernel_amp = kernel_amp_list[i_filter]

    miri_scattering_default = {
        "radius": radius,
        "kernel_amp": kernel_amp,
        "rotate_value": rotate_value
    }

    return miri_scattering_default


def _make_kernel(amplitude, radius, nsamples):
    """
    Creates a detector scatter kernel function k(x-a) which defines the fraction of the signal in a pixel at pixel
    location 'a' which is scattered into a pixel at position 'x-a' in the along row or along column direction. For
    simplicity, we assume a simple exponential dependence.

    Code is from MIRI-TN-00076-ATC_Imager_PSF_Issue_4.pdf (originally in IDL).
    """

    # Update values based on oversampling of PSF
    fold = 25.0 * nsamples  # e-folding length is 25 MIRI detector pixels
    amplitude /= nsamples  # the signal is shared between n samples
    radius *= nsamples

    # Generate kernel functions for the halo v pixel distance
    distances_list = np.arange(2 * radius + 1, dtype=float)
    distances_list = np.abs(distances_list - distances_list[radius])
    kernel = amplitude * np.exp(-distances_list / fold)

    return kernel, amplitude, fold


def _apply_kernel(kernel, image, radius):
    """
    Applies the detector scattering kernel function created in _make_kernel function to an input image. Code is
    from MIRI-TN-00076-ATC_Imager_PSF_Issue_4.pdf (originally in IDL).

    While the current code form isn't as elegant as it could be, it will stay like this for the time being until
    we can find a convolution method that takes less time to execute than this for loop method.
    """

    shape = image.shape
    y_num = shape[0]
    x_num = shape[1]
    out_image = np.zeros([y_num, x_num])

    # Stepping through each index of the original image
    for y in np.arange(0, y_num):
        for x in np.arange(0, x_num):

            x_i1 = x - radius
            if x_i1 < 0:
                x_i1 = 0
            n_left = x - x_i1
            x_k1 = radius - n_left

            x_i2 = x + radius
            if x_i2 > (x_num - 1):
                x_i2 = x_num - 1
            n_right = x_i2 - x
            x_k2 = radius + n_right

            out_image[y, x_i1:x_i2] += image[y, x] * kernel[x_k1: x_k2]

    return out_image


def apply_miri_scattering(HDUlist_or_filename=None, radius=None, kernel_amp=None, rotate_value=None):
    """
    Apply a distortion caused by the MIRI scattering cross artifact effect. Description of distortion and code is
    adapted from MIRI-TN-00076-ATC_Imager_PSF_Issue_4.pdf (originally in IDL).

    Parameters
    ----------

    HDUlist_or_filename :
        A PSF from WebbPSF, either as an HDUlist object or as a filename
    radius: float
        Radius of kernel profile (MIRI pixels). If set to None, value will be set at 200 pixels. Default = None.
    kernel_amp: float
        Detector scattering kernel amplitude. If set to None, function will pull the value based on best fit analysis
        based on the input PSF's filter. Default = None.
    rotate_value: float
        The rotation of the MIRI detector in degrees

    """

    # Read in input PSF
    if isinstance(HDUlist_or_filename, six.string_types):
        hdu_list = fits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, fits.HDUList):
        hdu_list = HDUlist_or_filename
    else:
        raise ValueError("input must be a filename or HDUlist")

    # Create a copy of the PSF
    psf = copy.deepcopy(hdu_list)

    # Log instrument and detector names
    instrument = hdu_list[0].header["INSTRUME"].upper()
    filter = hdu_list[0].header["FILTER"].upper()

    if instrument != "MIRI":
        raise ValueError("MIRI's Scattering Effect should only be applied to MIRI PSFs")

    # Pull default values
    miri_scattering_default = _get_default_miri(filter)

    # Set values if not already set by a keyword argument
    if radius is None:
        radius = miri_scattering_default["radius"]
    if kernel_amp is None:
        kernel_amp = miri_scattering_default["kernel_amp"]
    if rotate_value is None:
        rotate_value = miri_scattering_default["rotate_value"]

    ext = 2

    # Set over-sample value
    cdp_samp = psf[ext].header["OVERSAMP"]  # the over-sample value for this ext. If det, it'll = 1 so no effect

    # Read in PSF
    in_psf = psf[ext].data

    # Make the kernel
    kernel, amplitude, fold = _make_kernel(kernel_amp, radius, cdp_samp)

    # Create scattering images by applying the kernel vertically/horizontally via transposing
    x_scattered_image = _apply_kernel(kernel, in_psf, radius)
    kernel[radius] = 0.0  # set this value to 0 for the 2nd application so you don't apply this value 2x to 1 point

    in_psf_tr = in_psf.T  # apply the kernel to the y-direction of the images
    y_scattered_image_tr = _apply_kernel(kernel, in_psf_tr, radius)  # but your output will still be 1D in x-dir
    y_scattered_image = y_scattered_image_tr.T  # so then make it vertical to be applied later

    # Rotate the scattering images (but keep same size) so they match the PSF
    x_scattered_image_rot = rotate(x_scattered_image, rotate_value, reshape=False)
    y_scattered_image_rot = rotate(y_scattered_image, rotate_value, reshape=False)

    # Add the vertical/horizontal scattering images to the PSF
    psf_new = in_psf + x_scattered_image_rot + y_scattered_image_rot

    # Apply data to correct extensions
    psf[ext].data = psf_new

    # Now bin down over-sampled PSF to be detector-sampled and re-write ext=3
    detector_oversample = psf[ext].header["DET_SAMP"]
    psf[3].data = poppy.utils.rebin_array(psf_new, rc=(detector_oversample, detector_oversample))

    for ext in [2, 3]:

        # Set new header keywords
        psf[ext].header["MIR_DIST"] = ("True", "MIRI detector scattering applied")
        psf[ext].header["KERN_AMP"] = (amplitude, "Kernel Amplitude used in kernel exponential")
        psf[ext].header["KERNFOLD"] = (fold, "e-folding length used in kernel exponential")
        psf[ext].header["KERN_RAD"] = (radius, "Radius of kernel profile (MIRI pixels)")

    return psf
