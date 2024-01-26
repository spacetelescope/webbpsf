import copy

import os
import numpy as np
import scipy
import webbpsf
from webbpsf import utils, constants
from astropy.convolution.kernels import CustomKernel
from astropy.convolution import convolve
from astropy.io import fits
import astropy.convolution
import scipy.signal as signal


def get_detector_ipc_model(inst, header):
    """Retrieve detector interpixel capacitance model

    The details of the available calibration data vary per instrument.

    Parameters:
    -----------
    inst : string
        instrument name
    header : astropy.io.fits.Header
        FITS header

    Returns:
    --------
    kernel : numpy.ndarray
        Convolution kernel
    meta : dict
        Metadata about that kernel to be saved in FITS header

    """

    inst = inst.upper()
    det = header['DET_NAME'] #detector name

    meta = dict()

    if inst == 'NIRCAM':

        det2sca = {
            'NRCA1': '481', 'NRCA2': '482', 'NRCA3': '483', 'NRCA4': '484', 'NRCA5': '485',
            'NRCB1': '486', 'NRCB2': '487', 'NRCB3': '488', 'NRCB4': '489', 'NRCB5': '490',
        }

        webbpsf.webbpsf_core._log.info(f"Detector IPC: NIRCam {det} (added)")
        # IPC effect
        # read the SCA extension for the detector
        sca_path = os.path.join(utils.get_webbpsf_data_path(), 'NIRCam', 'IPC', 'KERNEL_IPC_CUBE.fits')
        kernel_ipc = CustomKernel(fits.open(sca_path)[det2sca[det]].data[0]) # we read the first slice in the cube

        # PPC effect
        # read the SCA extension for the detector
        ## TODO: This depends on detector coordinates, and which readout amplifier. if in subarray, then the PPC effect is always like in amplifier 1
        sca_path_ppc = os.path.join(utils.get_webbpsf_data_path(), 'NIRCam', 'IPC', 'KERNEL_PPC_CUBE.fits')
        kernel_ppc = CustomKernel(fits.open(sca_path_ppc)[det2sca[det]].data[0])  # we read the first slice in the cube

        kernel = (kernel_ipc, kernel_ppc)  # Return two distinct convolution kernels in this case

        meta['IPCINST'] = ('NIRCam', 'Interpixel capacitance (IPC)')
        meta['IPCTYPA'] = (det2sca[det], 'NRC SCA num used for IPC and PPC model')
        meta['IPCFILE'] = (os.path.basename(sca_path), 'IPC model source file')
        meta['PPCFILE'] = (os.path.basename(sca_path_ppc), 'PPC model source file')

    elif inst =='MIRI':
        webbpsf.webbpsf_core._log.info("Detector IPC: MIRI")

        a = webbpsf.constants.INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS[inst]
        alpha = webbpsf.constants.INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS[inst][0]
        beta = webbpsf.constants.INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS[inst][1]
        c = webbpsf.constants.INSTRUMENT_IPC_DEFAULT_KERNEL_PARAMETERS[inst][2]  # real observation noise adjustment
        miri_kernel = np.array([[c, beta, c],
                                [alpha, 1 - 2 * alpha - 2 * beta - 4 * c, alpha],
                                [c, beta, c]])
        kernel = CustomKernel(miri_kernel)

        meta['IPCINST'] = ('MIRI', 'Interpixel capacitance (IPC)')
        meta['IPCTYPA'] = (alpha, 'coupling coefficient alpha')
        meta['IPCTYPB'] = (beta, 'coupling coefficient beta')
        meta['IPCFILE'] = ('webbpsf.constants', 'IPC model source file')

    elif inst == 'NIRISS':
        # NIRISS IPC files distinguish between the 4 detector readout channels, and
        # whether or not the pixel is within the region of a large detector epoxy void
        # that is present in the NIRISS detector.

        # this set-up the input variables as required by Kevin Volk IPC code
        #image = psf_hdulist[ext].data
        xposition = header["DET_X"]
        yposition = header["DET_Y"]
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
            kernel = fits.getdata(ipc_file)
            #newimage = signal.fftconvolve(image, ipckernel, mode='same')
            meta['IPCINST'] = ('NIRISS', 'Interpixel capacitance (IPC)')
            meta['IPCTYPA'] = (ipcname, 'kernel file used for IPC correction')
            meta['IPCFILE'] = (os.path.basename(ipc_file), 'IPC model source file')
        else:
            kernel = None

            meta['IPCINST'] = ('NIRISS', 'Interpixel capacitance (IPC)')
            meta['IPCTYPA'] = ('NIRISS', 'No kernel file found')
            meta['IPCTYPB'] = ('NIRISS', 'No IPC correction applied')
            meta['IPCFILE'] = ('Not found', 'IPC model source file')
            webbpsf.webbpsf_core._log.info(f"NIRISS IPC kernel file {ipc_file} not found.")


    elif inst in ["FGS", "NIRSPEC", "WFI"]:
        kernel = None     # No IPC models yet implemented for these
        meta['IPCFILE'] = ('Not found', 'IPC model source file')

    return kernel, meta


def apply_detector_ipc(psf_hdulist, extname = 'DET_DIST'):
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


    Parameters
    ----------
    psf_hdulist : astropy.io.fits.HDUList
        A HDUList containing a webbpsf simulation result
    extname : string
        Which extension name to apply this to. This gets a bit tricky. In the normal calc_psf code path, this
        is applied to detector-sampled data, *after* binning the oversampled data to detector resolution. This
        is most intuitive, and in some sense better represents the actual physics of this effect. However in the
        psf_grid code path for making ePSFs, we need to be able to apply this model to oversampled PSFs.

    """

    # In cases for which the user has asked for the IPC to be applied to a not-present extension, we have nothing to add this to

    #if psf_hdulist is  None: return
    if extname not in psf_hdulist:
        webbpsf.webbpsf_core._log.debug(f"Skipping IPC simulation since ext {extname} is not found")
        return

    # This avoid applying IPC corrections twice
    keyword = 'IPCINST'
    if keyword in psf_hdulist[extname].header._keyword_indices: return

    inst = psf_hdulist[extname].header['INSTRUME'].upper()
    oversample = psf_hdulist[extname].header['OVERSAMP']

    kernel, meta = get_detector_ipc_model(inst, psf_hdulist[extname].header)
    if kernel is not None:

        if inst.upper()=='NIRCAM':
            # For NIRCam we have distinct models for IPC and PPC effects. Needs two convolutions.
            ipckernel, ppckernel = kernel

            if oversample !=1:
                ipckernel = oversample_ipc_model(ipckernel, oversample)
                ppckernel = oversample_ipc_model(ppckernel, oversample)

            out_ipc_0  = convolve(psf_hdulist[extname].data, ipckernel)
            out_ipc  = convolve(out_ipc_0, ppckernel)
        elif inst.upper()=='NIRISS':
            # the NIRISS code provided by Kevin Volk was developed for a different convolution function
            if oversample !=1:
                kernel = oversample_ipc_model(kernel, oversample)
            out_ipc = signal.fftconvolve(psf_hdulist[extname].data, kernel, mode='same')
        else:
            if oversample !=1:
                kernel = oversample_ipc_model(kernel, oversample)
            out_ipc  = convolve(psf_hdulist[extname].data, kernel)

        # apply kernel to DET_DIST
        psf_hdulist[extname].data = out_ipc

        # save metadata to header
        for key in meta:
            psf_hdulist[extname].header[key] = meta[key]
        psf_hdulist[extname].header.add_history("Applied detector interpixel capacitance (IPC) model")

    else:
        webbpsf.webbpsf_core._log.info("IPC corrections are not implemented yet for {}".format(inst))
        psf_hdulist[extname].header['IPCINST'] = (inst, 'No IPC correction applied')



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
    psf_hdulist[ext].header.add_history("Applied detector charge diffusion model.")
    psf_hdulist[ext].header['CHDFTYPE'] = ('gaussian', 'Type of detector charge diffusion model')
    psf_hdulist[ext].header['CHDFSIGM'] = (sigma, '[arcsec] Gaussian sigma for charge diff model')
    psf_hdulist[ext].data = out

    return psf_hdulist


def oversample_ipc_model(kernel, oversample):
    """ Transform an IPC model convolution kernel to be applied to oversampled data.

    The correct way to do this turns out to be to intersperse zeros into the array, turning it
    into a sparse comb function. This is because the IPC is a discrete effect that acts on pixels,
    rather than a continuous function.

    (This is non-intuitive but mathematically yields precisely consistent results for either order
    of binning then applying IPC, or applying IPC then binning).

    Parameters
    ----------
    kernel : numpy.ndarray
        Convolution kernel for IPC model
    oversample : int
        Oversampling factor

    Returns a version of the kernel resampled and padded for use on oversampled data, for instance an ePSF

    """



    oversampling_kernel = np.zeros((oversample, oversample))
    oversampling_kernel[(oversample-1)//2, (oversample-1)//2] = 1

    kernel_oversample = np.kron(kernel, oversampling_kernel)

    if oversample % 2 == 0:
        # pad with an extra row and column of zeros, to convert into a symmetrical and odd-sized kernel
        npix = kernel_oversample.shape[0]
        padded_kernel = np.zeros((npix+1, npix+1))
        padded_kernel[1:, 1:] = kernel_oversample
        kernel_oversample =  padded_kernel

    return kernel_oversample


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

    # reduce intensity in the inner part, since the cruciform is suppressed at small radii
    kernel_x[np.abs(x) < constants.MIRI_CRUCIFORM_INNER_RADIUS_PIX] *= 0.5

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
    # In-flight correction based on measured cycle 1 ePSFs, coarsely
    for k in kernel_amp_corrections:
        kernel_amp_corrections[k] *= 0.5

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
