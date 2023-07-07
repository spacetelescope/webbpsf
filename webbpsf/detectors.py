import copy

import os
import numpy as np
import scipy
import webbpsf
from webbpsf import utils
from astropy.convolution.kernels import CustomKernel
from astropy.convolution import convolve
from astropy.io import fits
import astropy.convolution
import scipy.signal as signal




def apply_detector_ipc(psf_hdulist):
    """Apply a model for interpixel capacitance
    NIRCam: IPC and PPC values derived during ground I&T, primarily ISIM-CV3 from Jarron Leisenring
    these IPC/PPC kernels will be update after flight values are available.
    For NIRCam only PPC effects are also included, these are relatively small compared to the IPC contribution
    MIRI: Convolution kernels from JWST-STScI-002925 by Mike Engesser
    NIRISS: Convolution kernels and base code provided by Kevin Volk
    The IPC kernel files are derived from IPC measurements made from NIRISS commissioning dark ramps by Chris Willott.

    For NIRISS the user needs to have the right kernels under $WEBBPSF_PATH/NIRISS/IPC/
    These kernels should be available with webbpsf data > Version 1.1.1

    You can turn On/Off IPC effects as an option.
     For example: inst.option['add_ipc'] = False, where inst is the instrument class. Default is True.

       """

    inst = psf_hdulist[0].header['INSTRUME'].upper()
    if inst =='NIRCAM':



        det2sca = {
            'NRCA1': '481', 'NRCA2': '482', 'NRCA3': '483', 'NRCA4': '484', 'NRCA5': '485',
            'NRCB1': '486', 'NRCB2': '487', 'NRCB3': '488', 'NRCB4': '489', 'NRCB5': '490',
        }

        webbpsf.webbpsf_core._log.info("Detector IPC: NIRCam (added)")
        ext = 'DET_DIST'
        det = psf_hdulist[ext].header['DET_NAME'] #detector name
        # IPC effect
        # read the SCA extension for the detector
        sca_path = os.path.join(utils.get_webbpsf_data_path(), 'NIRCam', 'IPC', 'KERNEL_IPC_CUBE.fits')
        kernel = CustomKernel(fits.open(sca_path)[det2sca[det]].data[0]) # we read the first slide in the cube
        # apply to DET_DIST
        out_ipc  = convolve(psf_hdulist[ext].data, kernel)

        # PPC effect
        # read the SCA extension for the detector
        sca_path = os.path.join(utils.get_webbpsf_data_path(), 'NIRCam', 'IPC', 'KERNEL_PPC_CUBE.fits')
        kernel = CustomKernel(fits.open(sca_path)[det2sca[det]].data[0])  # we read the first slide in the cube
        # apply to DET_DIST after IPC effect
        out_ipc_ppc = convolve(out_ipc, kernel)

        psf_hdulist[ext].header['IPCINST'] = ('NIRCam', 'Interpixel capacitance (IPC)')
        psf_hdulist[ext].header['IPCTYPA'] = (det2sca[det], 'Post-Pixel Coupling (PPC)')
        psf_hdulist[ext].data = out_ipc_ppc

    if inst =='MIRI':
        a = webbpsf.constants.INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS[inst]
        webbpsf.webbpsf_core._log.info("Detector IPC: MIRI")
        alpha = webbpsf.constants.INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS[inst][0]
        beta = webbpsf.constants.INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS[inst][1]
        c = webbpsf.constants.INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS[inst][2]  # real observation noise adjustment
        miri_kernel = np.array([[c, beta, c],
                                [alpha, 1 - 2 * alpha - 2 * beta - 4 * c, alpha],
                                [c, beta, c]])
        kernel = CustomKernel(miri_kernel)
        # apply to DET_DIST
        ext = 'DET_DIST'
        out  = convolve(psf_hdulist[ext].data, kernel)
        psf_hdulist[ext].header['IPCINST'] = ('MIRI', 'Interpixel capacitance (IPC)')
        psf_hdulist[ext].header['IPCTYPA'] = (alpha, 'coupling coefficient alpha')
        psf_hdulist[ext].header['IPCTYPB'] = (beta, 'coupling coefficient beta')
        psf_hdulist[ext].data = out

    #
    if inst == 'NIRISS':
        ext = 'DET_DIST'
        # this set-up the input variables as required by Kevin Volk IPC code
        image = psf_hdulist[ext].data
        xposition = psf_hdulist[ext].header["DET_X"]
        yposition = psf_hdulist[ext].header["DET_Y"]
        # find the voidmask fits file
        voidmask10 = os.path.join(utils.get_webbpsf_data_path(), 'NIRISS', 'IPC', 'voidmask10.fits')

        if os.path.exists(voidmask10):
            maskimage = fits.getdata(voidmask10)
        else:
            maskimage = None
            webbpsf.webbpsf_core._log.info("Error reading the file voidmask10.fits.  Will assume a non-void position.")

        nchannel = int(yposition) // 512
        try:
            flag = maskimage[nchannel, int(xposition)]
        except:
            # This marks the pixel as non-void by default if the maskimage is not
            # read in properly
            flag = 0
        frag1 = ['A', 'B', 'C', 'D']
        frag2 = ['notvoid', 'void']

        ipcname = 'ipc5by5median_amp' + frag1[nchannel] + '_' + \
                      frag2[flag] + '.fits'
        ipc_file = os.path.join(utils.get_webbpsf_data_path(), 'NIRISS', 'IPC', ipcname)
        if os.path.exists(ipc_file):
            ipckernel = fits.getdata(ipc_file)
            newimage = signal.fftconvolve(image, ipckernel, mode='same')
            psf_hdulist[ext].header['IPCINST'] = ('NIRISS', 'Interpixel capacitance (IPC)')
            psf_hdulist[ext].header['IPCTYPA'] = (ipcname, 'kernel file used for IPC correction')
        else:
            newimage = numpy.copy(image)
            psf_hdulist[ext].header['IPCINST'] = ('NIRISS', 'Interpixel capacitance (IPC)')
            psf_hdulist[ext].header['IPCTYPA'] = ('NIRISS', 'No kernel file found')
            sf_hdulist[ext].header['IPCTYPB'] = ('NIRISS', 'No IPC correction applied')
            webbpsf.webbpsf_core._log.info("No IPC correction for NIRISS. Check kernel files.")

        psf_hdulist[ext].data = newimage

    if inst in ["FGS", "NIRSpec"]:
        ext = 'DET_DIST'
        psf_hdulist[ext].header['IPCINST'] = (inst, 'No IPC correction applied')
        webbpsf.webbpsf_core._log.info("IPC corrections are not implemented yet for {}".format(inst))


    return psf_hdulist


def apply_detector_charge_diffusion(psf_hdulist, options):
    """Apply a model for charge diffusion of photoelectrons within an H2RG
    This is a PLACEHOLDER, temporary heuristic

    """

    sigma = options.get('charge_diffusion_sigma')

    if sigma is None:
        # look up default from constants
        inst = psf_hdulist[0].header['INSTRUME'].upper()
        key = f"NIRCAM_{psf_hdulist[0].header['CHANNEL'][0]}W" if inst=='NIRCAM' else inst
        sigma = webbpsf.constants.INSTRUMENT_DETECTOR_CHARGE_DIFFUSION_DEFAULT_PARAMETERS[key]

    ext = 1  # Apply to the 'OVERDIST' extension

    webbpsf.webbpsf_core._log.info("Detector charge diffusion: Convolving with Gaussian with sigma={0:.3f} arcsec".format(sigma))
    out = scipy.ndimage.gaussian_filter(psf_hdulist[ext].data, sigma / psf_hdulist[0].header['PIXELSCL'])
    psf_hdulist[ext].header['CHDFTYPE'] = ('gaussian', 'Type of detector charge diffusion model')
    psf_hdulist[ext].header['CHDFSIGM'] = (sigma, '[arcsec] Gaussian sigma for charge diff model')
    psf_hdulist[ext].data = out

    return psf_hdulist

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
        This will be transposed to create the kernel in the y direction.
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
                                                 normalize_kernel=False, nan_treatment='fill', allow_huge = True)

    # Transpose to make a kernel for Y and convolve with that too
    im_conv_y = astropy.convolution.convolve_fft(in_psf, kernel_x.T, boundary='fill', fill_value=0.0,
                                                 normalize_kernel=False, nan_treatment='fill', allow_huge = True)

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

    # The above values are from that tech report, but empirically we need higher values to
    # better match the MIRI CDP PSFS. See e.g. MIRI_FM_MIRIMAGE_F560W_PSF_07.02.00.fits
    # and https://github.com/spacetelescope/webbpsf/issues/415
    kernel_amp_corrections = {'F560W': 4.05, 'F770W': 4.1, 'F1000W': 3.8,
                               'F1130W': 2.5, 'F1280W': 2.5, 'F1065C': 2.5, 'F1140C': 2.5}

    # Set values if not already set by a keyword argument
    if kernel_amp is None:
        kernel_amp = kernel_amp_dict[filt]

        if filt in kernel_amp_corrections:
            kernel_amp *= kernel_amp_corrections[filt]

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