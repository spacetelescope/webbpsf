import copy

import astropy.convolution
import astropy.io.fits as fits
import numpy as np
import pysiaf
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import rotate


def _get_default_siaf(instrument, aper_name):
    """
    Create instance of pysiaf for the input instrument and aperture
    to be used later to pull SIAF values like distortion polynomial
    coefficients and rotation.

    Parameters
    ----------
    instrument : str
        The name of the instrument
    aper_name : str
        The name of the specific instrument aperture

    Returns
    -------
    aper : instance of pysiaf
    """

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


# Function for applying distortion from SIAF polynomials

def apply_distortion(hdulist_or_filename=None, fill_value=0):
    """
    Apply a distortion to the input PSF. The distortion comes from the SIAF 4-5 degree polynomial
    (depending on the instrument). This function pulls and applies the SIAF polynomial values
    using pysiaf package, which ensures the most up-to-date values will be called.

    Parameters
    ----------
    hdulist_or_filename :
        A PSF from WebbPSF, either as an HDUlist object or as a filename
    fill_value : float
        Value used to fill in any blank space by the skewed PSF. Default = 0

    Returns
    -------
    psf : HDUlist object
        PSF with distortion applied from SIAF polynomial
    """

    # Read in input PSF
    if isinstance(hdulist_or_filename, str):
        hdu_list = fits.open(hdulist_or_filename)
    elif isinstance(hdulist_or_filename, fits.HDUList):
        hdu_list = hdulist_or_filename
    else:
        raise ValueError("input must be a filename or HDUlist")

    # Create a copy of the PSF
    psf = copy.deepcopy(hdu_list)

    # Log instrument and detector names
    instrument = hdu_list[0].header["INSTRUME"].upper()
    aper_name = hdu_list[0].header["APERNAME"].upper()

    # Pull default values
    aper = _get_default_siaf(instrument, aper_name)

    ext = 1  # edit the oversampled PSF (OVERDIST extension)

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

    # 6) Shift the sci indices so they match the PSF's position again (moved slightly off from pysiaf calculation)
    xsci += xpix_center - np.median(xsci)
    ysci += ypix_center - np.median(ysci)

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

    # Set new header keywords
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


# Function for applying Rotation to NIRCam, NIRISS, and FGS

def apply_rotation(hdulist_or_filename=None, rotate_value=None, crop=True):
    """
    Apply the detector's rotation to the PSF. This is for NIRCam, NIRISS, and FGS.
    MIRI and NIRSpec's large rotation is already added inside WebbPSF's calculations.

    Parameters
    ----------
    hdulist_or_filename :
        A PSF from WebbPSF, either as an HDUlist object or as a filename
    rotate_value : float
        Rotation in degrees that PSF needs to be. If set to None, function
        will pull the most up to date SIAF value. Default = None.
    crop : bool
        True or False to crop the PSF so it matches the size of the input
        PSF (e.g. so they could be more easily compared).

    Returns
    -------
    psf : HDUlist object
        PSF with rotation applied from SIAF values
    """

    # Read in input PSF
    if isinstance(hdulist_or_filename, str):
        hdu_list = fits.open(hdulist_or_filename)
    elif isinstance(hdulist_or_filename, fits.HDUList):
        hdu_list = hdulist_or_filename
    else:
        raise ValueError("input must be a filename or HDUlist")

    # Create a copy of the PSF
    psf = copy.deepcopy(hdu_list)

    # Log instrument and detector names
    instrument = hdu_list[0].header["INSTRUME"].upper()
    aper_name = hdu_list[0].header["APERNAME"].upper()

    if instrument in ["MIRI", "NIRSPEC"]:
        raise ValueError("{}'s rotation is already included in WebbPSF and "
                         "shouldn't be added again.".format(instrument))

    # Set rotation value if not already set by a keyword argument
    if rotate_value is None:
        aper = _get_default_siaf(instrument, aper_name)
        rotate_value = getattr(aper, "V3IdlYAngle")  # the angle to rotate the PSF in degrees

    # If crop = True, then reshape must be False - so invert this keyword
    reshape = np.invert(crop)

    ext = 1  # edit the oversampled PSF (OVERDIST extension)

    psf_new = rotate(psf[ext].data, rotate_value, reshape=reshape)

    # Apply data to correct extensions
    psf[ext].data = psf_new

    # Set new header keyword
    psf[ext].header["ROTATION"] = (rotate_value, "PSF rotated to match detector rotation")

    return psf


# Functions for applying MIRI Detector Scattering Effect

def _make_miri_scattering_kernel(image, amplitude, nsamples):
    """
    Creates a detector scatter kernel function. For simplicity, we assume a
    simple exponential dependence. Code is adapted from
    MIRI-TN-00076-ATC_Imager_PSF_Issue_4.pdf (originally in IDL).

    Parameters
    ----------
    image : ndarray
        PSF array for which to make the kernel
    amplitude : float
        Amplitude of the kernel
    nsamples : int
        Amount by which the input PSF is oversampled

    Returns
    -------
    kernel_x : ndarray
        1D detector scattering kernel in the x direction
    """

    # Compute 1d indices
    x = np.arange(image.shape[1], dtype=float)
    x -= (image.shape[1]-1)/2
    x /= nsamples

    # Create 1d kernel
    kernel_x = amplitude * np.exp(-np.abs(x) / 25)

    # Reshape kernel to 2D image for use in convolution
    kernel_x.shape = (1, image.shape[1])

    return kernel_x


def _apply_miri_scattering_kernel(in_psf, kernel_x, oversample):
    """
    Applies the detector scattering kernel created in _make_miri_scattering_kernel
    function to an input image. Code is adapted from
    MIRI-TN-00076-ATC_Imager_PSF_Issue_4.pdf

    Parameters
    ----------
    in_psf : ndarray
        PSF array upon which to apply the kernel
    kernel_x : ndarray
        The 1D kernel in the x direction, output from _make_miri_scattering_kernel.
        This will be transposed to createt the kernel in the y direction.
    oversample : int
        Amount by which the input PSF is oversampled

    Returns
    -------
    im_conv_both : ndarray
        The input image convolved with the input kernel in both the x and
        y directions
    """

    # Apply the kernel via convolution in both the X and Y direction
    # Convolve the input PSF with the kernel for scattering in the X direction
    im_conv_x = astropy.convolution.convolve_fft(in_psf, kernel_x, boundary='fill', fill_value=0.0,
                                                 normalize_kernel=False, nan_treatment='fill')

    # Transpose to make a kernel for Y and convolve with that too
    im_conv_y = astropy.convolution.convolve_fft(in_psf, kernel_x.T, boundary='fill', fill_value=0.0,
                                                 normalize_kernel=False, nan_treatment='fill')

    # Sum together both the X and Y scattering.
    # Note, it appears we do need to correct the amplitude for the sampling factor. Might as well do that here.
    im_conv_both = (im_conv_x + im_conv_y)/(oversample**2)

    return im_conv_both


def apply_miri_scattering(hdulist_or_filename=None, kernel_amp=None):
    """
    Apply a distortion caused by the MIRI scattering cross artifact effect.
    In short we convolve a 2D exponentially decaying cross to the PSF where
    the amplitude of the exponential function is determined by the filter of
    the PSF. A full description of the distortion and the original code can
    be found in MIRI-TN-00076-ATC_Imager_PSF_Issue_4.pdf

    Note, this code **edits in place Extension 1 of the supplied HDUlist**. In the typical case where the
    input PSF is calculated as Extension 0, the calling function must put a copy of that into Extension 1
    which this will then modify. This happens in webbpsf_core.py/JWInstrument._calc_psf_format_output,
    which is where this is called from in the usual course of operation.

    Parameters
    ----------
    hdulist_or_filename :
        A PSF from WebbPSF, either as an HDUlist object or as a filename
    kernel_amp: float
        Detector scattering kernel amplitude. If set to None,
        function will pull the value based on best fit analysis
        using the input PSF's filter. Default = None.

    Returns
    -------
    psf : HDUlist object
        PSF with MIRI detector scattering effect applied
    """

    # Read in input PSF
    if isinstance(hdulist_or_filename, str):
        hdu_list = fits.open(hdulist_or_filename)
    elif isinstance(hdulist_or_filename, fits.HDUList):
        hdu_list = hdulist_or_filename
    else:
        raise ValueError("input must be a filename or HDUlist")

    # Create a copy of the PSF
    psf = copy.deepcopy(hdu_list)

    # Log instrument name and filter
    instrument = hdu_list[0].header["INSTRUME"].upper()
    filt = hdu_list[0].header["FILTER"].upper()

    if instrument != "MIRI":
        raise ValueError("MIRI's Scattering Effect should only be applied to MIRI PSFs")

    # Default kernel amplitude values from modeling in MIRI-TN-00076-ATC_Imager_PSF_Issue_4.pdf
    kernel_amp_dict = {'F560W': 0.00220, 'F770W': 0.00139, 'F1000W': 0.00034,
                       'F1130W': 0.00007, 'F1280W': 0.00011, 'F1500W': 0.0,
                       'F1800W': 0.0, 'F2100W': 0.0, 'F2550W': 0.0, 'FND': 0.00087,
                       'F1065C': 0.00010, 'F1140C': 0.00007, 'F1550C': 0.0,
                       'F2300C': 0.0}

    # Set values if not already set by a keyword argument
    if kernel_amp is None:
        kernel_amp = kernel_amp_dict[filt]

    ext = 1  # edit the oversampled PSF (OVERDIST extension)

    # Set over-sample value
    oversample = psf[ext].header["DET_SAMP"]

    # Read in PSF
    in_psf = psf[ext].data

    # Make the kernel
    kernel_x = _make_miri_scattering_kernel(in_psf, kernel_amp, oversample)

    # Apply the kernel via convolution in both the X and Y direction to produce a 2D output
    im_conv_both = _apply_miri_scattering_kernel(in_psf, kernel_x, oversample)

    # Add this 2D scattered light output to the PSF
    psf_new = in_psf + im_conv_both

    # To ensure conservation of intensity, normalize the psf
    psf_new *= in_psf.sum() / psf_new.sum()

    # Apply data to correct extensions
    psf[ext].data = psf_new

    # Set new header keywords
    psf[ext].header["MIR_DIST"] = ("True", "MIRI detector scattering applied")
    psf[ext].header["KERN_AMP"] = (kernel_amp, "Amplitude (A) in kernel function A*exp(-x/B)")
    psf[ext].header["KERNFOLD"] = (25, "e-folding length (B) in kernel func A*exp(-x/B)")

    return psf
