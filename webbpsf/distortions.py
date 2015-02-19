import os
from astropy.io import fits
import poppy
from poppy import zernike
import numpy as np

import logging
_log = logging.getLogger(__name__)
_log.setLevel(logging.WARNING)

from . import utils

def jwexike_basis(nterms=15, rho=None, theta=None, pupil_filename='pupil_revV.fits'):
    """Return a cube of distortion polynomials orthonormalized and evaluated
    over the JWST pupil.

    Parameters
    ----------
    nterms : int
        Number of JWexike terms to return
    rho, theta : array_like
        Image plane coordinates. `rho` should be 1 at the desired pixel radius,
        `theta` should be in radians
    """
    pupil_path = os.path.join(utils.get_webbpsf_data_path(), pupil_filename)
    pupil = fits.getdata(pupil_path)
    npix = pupil.shape[0]

    if rho is not None:
        if rho is None or theta is None:
            raise ValueError("You must supply both `theta` and `rho`, or neither.")
        npix = rho.shape[0]
        shape = rho.shape
        use_polar = True
    else:
        shape = (npix, npix)
        use_polar = False

def derive_jwexikes(nterms=11, pupil=None):
    """Derive a set of orthonormal polynomials over the JWST pupil, following the
    method of Mahajan and Dai 2006

    Parameters
    ----------
    nterms : int
        How many terms to compute, starting with the piston term for `nterms` = 1 and counting
        up analogous to the indexing convention for Zernikes in Noll et al. 1976
    pupil : array
        Square array with pupil transmittance values in the range [0, 1]. Output will be the same
        shape as `pupil`.

    Returns
    -------
    jwexikes : array of shape (nterms,) + pupil.shape
        The orthonormal polynomials evaluated over the pupil provided, in order according to
        the indexing convention of Noll et al. 1976
    """

    A = pupil.sum()

    # precompute zernikes
    shape = pupil.shape
    Z = [np.zeros(shape), ]
    for j in range(nterms+1):
        Z.append(zernike.zernike1(j + 1, npix=shape[0], outside=0., mask_outside=False))

    G = [np.zeros(shape), np.ones(shape)]  # array of G_i etc. intermediate fn
    H = [np.zeros(shape), np.ones(shape) * pupil]  # array of hexikes
    c = {}  # coefficients hash

    for j in np.arange(nterms - 1) + 1:  # can do one less since we already have the piston term
        _log.debug("  j = " + str(j))
        # Compute the j'th G, then H
        nextG = Z[j+1] * pupil
        for k in np.arange(j) + 1:
            #print "   k=%d\tbefore: %s" %(k, np.isfinite(nextG.sum()))
            c[(j+1, k)] = -1/A * (Z[j+1] * H[k] * pupil).sum()
            if c[(j+1, k)] != 0:
                nextG += c[(j+1, k)] * H[k]
            #print "    c[%s] = %f\tafter: %s" % (str((j+1, k)), c[(j+1, k)], np.isfinite(nextG.sum()) )
            _log.debug("    c[%s] = %f" % (str((j+1, k)), c[(j+1, k)] ))

        #print "   nextG integral: %f" % sqrt((nextG**2).sum() / A )
        nextH = nextG / np.sqrt((nextG ** 2).sum() / A)

        G.append(nextG)
        H.append(nextH)

        #TODO - contemplate whether the above algorithm is numerically stable
        # cf. modified gram-schmidt algorithm discussion on wikipedia.

    # drop the 0th null element, return the rest
    return H[1:], c, pupil.copy()


def plot_jwexikes(nterms=20, **kwargs):
    """ Test the jwexikes functions and display the results """
    plotny = int(np.floor(np.sqrt(nterms)))
    plotnx = int(nterms/plotny)

    from matplotlib import pyplot as plt
    fig = plt.gcf()
    fig.clf()

    H, c, ap = derive_jwexikes(nterms=nterms, **kwargs)

    wgood = np.where(ap != 0)
    ap[np.where(ap == 0)] = np.nan

    for j in np.arange(nterms):
        ax = fig.add_subplot(plotny,plotnx,j+1, frameon=False, xticks=[], yticks=[])

        n, m = zernike.noll_indices(j+1)

        ax.imshow(H[j] * ap , vmin=-3, vmax=3.0)
        #print "j = %d\tzmin = %f\tzmax=%f" % (j, np.nanmin(Z), np.nanmax(Z))
        npix = ap.shape[0]
        ax.text(npix*0.7, npix*0.1, "$JW_%d^{%d}$" % (n, m), fontsize=20)
        print "Term %d:   std dev is %f. (should be near 1)" % (j+1, H[j][wgood].std())

    plt.draw()