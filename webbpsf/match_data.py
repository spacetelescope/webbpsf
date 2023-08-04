## Functions to match or fit PSFs to observed JWST data
import astropy
import astropy.units as u
import astropy.io.fits as fits

import scipy.optimize

import webbpsf


def setup_sim_to_match_file(filename_or_HDUList, verbose=True, plot=False):
    """ Setup a webbpsf Instrument instance matched to a given 
    """
    if isinstance(filename_or_HDUList,str):
        if verbose:
            print(f"Setting up sim to match {filename_or_HDUList}")
        header = fits.getheader(filename_or_HDUList)
    else:
        header = filename_or_HDUList[0].header
        if verbose:
            print(f"Setting up sim to match provided FITS HDUList object")

    inst = webbpsf.instrument(header['INSTRUME'])

    if inst.name=='MIRI' and header['FILTER']=='P750L':
        # webbpsf doesn't model the MIRI LRS prism spectral response
        print("Please note, webbpsf does not currently model the LRS spectral response. Setting filter to F770W instead.")
        inst.filter='F770W'
    else:
        inst.filter=header['filter']
    inst.set_position_from_aperture_name(header['APERNAME'])

    dateobs = astropy.time.Time(header['DATE-OBS']+"T"+header['TIME-OBS'])
    inst.load_wss_opd_by_date(dateobs, verbose=verbose, plot=plot)


    # per-instrument specializations
    if inst.name == 'NIRCam':
        if header['PUPIL'].startswith('MASK'):
            inst.pupil_mask = header['PUPIL']
            inst.image_mask = header['CORONMSK'].replace('MASKA', 'MASK')  # note, have to modify the value slightly for
                                                                           # consistency with the labels used in webbpsf
    elif inst.name == 'MIRI':
        if inst.filter in ['F1065C', 'F1140C', 'F1550C']:
            inst.image_mask = 'FQPM'+inst.filter[1:5]
        elif inst.filter == 'F2300C':
            inst.image_mask = 'LYOT2300'
        elif header['FILTER'] == 'P750L':
            inst.pupil_mask = 'P750L'
            if header['APERNAME'] == 'MIRIM_SLIT':
                inst.image_mask = 'LRS slit'

    # TODO add other per-instrument keyword checks

    if verbose:
        print(f"""
Configured simulation instrument for:
    Instrument: {inst.name}
    Filter: {inst.filter}
    Detector: {inst.detector}
    Apername: {inst.aperturename}
    Det. Pos.: {inst.detector_position} {'in subarray' if "FULL" not in inst.aperturename else ""}
    Image plane mask: {inst.image_mask}
    Pupil plane mask: {inst.pupil_mask}
    """)

    return inst

