
import re
from copy import deepcopy

import numpy as np
from scipy.io import loadmat

from scipy import ndimage
import webbpsf



matfile = '/itar/jwst/tel/share/webbpsf/webbpsf-data-source/JWST_thermal_response_data05_31_2017.mat'
def load_matfile(matfile):
    """Load MATLAB file"""
    return loadmat(matfile, mat_dtype=True, struct_as_record=False)

def grab_opd_from_mat(mat, key):
    '''
    Given a MATLAB file with delta OPDs,
    return a properly-scaled OPD map
    in nd array (in microns).

    Follows directions provided by S. Knight.
    '''
    mat = load_matfile(matfile)

    # Grab the data from the requested key
    data = deepcopy(mat[key][0, 0])
    opd = data.MAP

    # Scale to microns
    opd *= data.WVL[0, 0] / data.SSZ[0, 0]

    return opd[::-1] #flip y-axis

def grab_dates_from_mat(matfile=matfile):
    """
    Get the list of dates of the thermal OPDs from the matfile
    Return in order of secs to days
    """
    mat = load_matfile(matfile)

    all_dates = [d for d in mat.keys() if 'day' in d]
    secs = sorted([s for s in all_dates if 'sec' in s], key=natural_key)
    days = sorted(set(all_dates) ^ set(secs), key=natural_key)
    return secs + days

def natural_key(string_):
    'sory list by natural numbers'
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_) if s]

def get_known_pmsa_centers():
    aperture = webbpsf.optics.WebbPrimaryAperture()
    return aperture.seg_centers

def vcoords_to_pixels(vx, vy, image_dim, pixelscale):
    pixx = vx / pixelscale
    pixy = vy / pixelscale
    center = ((image_dim[0] - 1) / 2., (image_dim[1] - 1) / 2.)
    return pixx + center[1], pixy + center[0]

def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)



def grab_opd_info(opd_data_str, npix=1024):
    """
    Get the information about the JWST OTE

    Parameters:
    ----------
    opd_data_str: str
        The name of the OPD being queried from JWST_thermal_response_data05_31_2017.mat
    npix: int
        Size of the image being created

    Returns
    -------
    opd: array
        Array of the OPD asked for
    labeldict: dictionary
        A nested dictionary containing masks (labels) for each JWST segment
    newlabels: array
        A combined label plot with segments properly combined
    newcentroids: dictionary
        The center of each segment in pixels(?)
    segradius: float
        The radius of each segment
    """
    # Label the segments
    opd = grab_opd_from_mat(matfile, opd_data_str)
    labels, nfeatures = ndimage.label(opd)
    yind, xind = np.indices(labels.shape)

    # Take mean position of each segment and then sort
    centroids = np.array([(np.mean(yind[labels == i]),
                           np.mean(xind[labels == i])) for i in np.arange(1, nfeatures+1)])

    # Get ideal centroids
    known_centers = get_known_pmsa_centers()

    # Handle obscured segments and construct segment label dictionary
    labeldict = {}
    for seg, (vx, vy) in known_centers.items():
        pixscale = 0.006526714427860697
        segradius = 1.5 / 2. / pixscale #pixels
        pixx, pixy = vcoords_to_pixels(vx, vy, (npix, npix), pixscale)
        alldist = np.array([distance(x, y, pixx, pixy) for (y, x) in centroids])
        matches = np.where(alldist < segradius)[0]
        labeldict[seg[:2]] = np.zeros((npix, npix))
        for m in matches:
            labeldict[seg[:2]][labels == (m+1)] = np.unique(labels[labels == (matches[0]+1)]) # assign val of first match

    newcentroids = {key:(np.mean(yind[im > 0]), np.mean(xind[im > 0 ])) for key, im in labeldict.items()}

    # Make new (combined) label plot with segments properly combined
    newlabels = np.zeros((npix, npix))
    for val in labeldict.values():
        newlabels += val

    return opd, labeldict, newlabels, newcentroids, segradius
