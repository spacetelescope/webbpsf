import numpy as np
import scipy
import webbpsf



def apply_h2rg_charge_diffusion(psf_hdulist, options):
    """Apply a model for charge diffusion of photoelectrons within an H2RG
    This is a PLACEHOLDER, temporary heuristic

    """

    sigma = options.get('charge_diffusion_sigma', 0.)

    webbpsf.webbpsf_core._log.info("Calculating jitter using " + str(options['jitter']))

    webbpsf.webbpsf_core._log.info("Jitter: Convolving with Gaussian with sigma={0:.3f} arcsec".format(sigma))

    out = scipy.ndimage.gaussian_filter(psf_hdulist[0].data, sigma / psf_hdulist[0].header['PIXELSCL'])

    psf_hdulist[0].header['CHDFTYPE'] = ('gaussian', 'Type of model for detector charge diffusion applied')
    psf_hdulist[0].header['CHDFSIGM'] = ('gaussian', 'Gaussian sigma for detector `charge diffusion model [arcsec]')
    psf_hdulist[0].data = out

    return psf_hdulist