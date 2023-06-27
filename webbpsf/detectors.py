import os
import numpy as np
import scipy
import webbpsf
from webbpsf import utils
from astropy.convolution.kernels import CustomKernel
from astropy.convolution import convolve
import scipy.signal as signal
from astropy.io import fits

def apply_detector_ipc(psf_hdulist):
    """Apply a model for interpixel capacitance
    NIRCam: IPC and PPC values derived during ground I&T, primarily ISIM-CV3 from Jarron Leisenring
    these IPC/PPC kernels will be update it after flight values are available.
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
            psf_hdulist[ext].header['IPCTYPA'] = ('NIRISS', 'no kernel file found')
            sf_hdulist[ext].header['IPCTYPB'] = ('NIRISS', 'no IPC correction applied')
            webbpsf.webbpsf_core._log.info("No IPC correction for NIRISS. Check kernel files.")

        psf_hdulist[ext].data = newimage


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
