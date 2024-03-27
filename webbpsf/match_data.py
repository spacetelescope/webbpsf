## Functions to match or fit PSFs to observed JWST data
import astropy
import astropy.units as u
import astropy.io.fits as fits

import scipy.optimize

import webbpsf
import pysiaf


def setup_sim_to_match_file(filename_or_HDUList, verbose=True, plot=False, choice='closest'):
    """ Setup a webbpsf Instrument instance matched to a given dataset

    Parameters
    ----------
    filename_or_HDUlist : str or astropy.io.fits HDUList
        file to load
    verbose : bool
        be more verbose?
    plot : bool
        plot?
    choice : string
        Method to choose which OPD file to use, e.g. 'before', 'after', or 'closest'
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
    inst.load_wss_opd_by_date(dateobs, verbose=verbose, plot=plot, choice=choice)


    # per-instrument specializations
    if inst.name == 'NIRCam':
        if header['PUPIL'].startswith('MASK'):
            if header['PUPIL'] == 'MASKBAR':
                inst.pupil_mask = header['CORONMSK'].replace('MASKA', 'MASK')
            else:
                inst.pupil_mask = header['PUPIL']
            inst.image_mask = header['CORONMSK'].replace('MASKA', 'MASK')  # note, have to modify the value slightly for
                                                                           # consistency with the labels used in webbpsf

            # The apername keyword is not always correct for cases with dual-channel coronagraphy
            # in some such cases, APERNAME != PPS_APER. Let's ensure we have the proper apername for this channel:
            apername = get_coron_apname(header)
            inst.set_position_from_aperture_name(apername)


        elif header['PUPIL'].startswith('F'):
            inst.filter = header['PUPIL']
        else:
            inst.pupil_mask = header['PUPIL']


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



def get_coron_apname(input):
    """Get aperture name from header or data model

    Handles edge cases for dual-channel coronagraphy.

    By Jarron Leisenring originally in webbpsf_ext, copied here by permission

    Parameters
    ==========
    input : fits.header.Header or datamodels.DataModel
        Input header or data model
    """

    if isinstance(input, (fits.header.Header)):
        # Aperture names
        apname = input['APERNAME']
        apname_pps = input['PPS_APER']
        subarray = input['SUBARRAY']
    else:
        # Data model meta info
        meta = input.meta

        # Aperture names
        apname = meta.aperture.name
        apname_pps = meta.aperture.pps_name
        subarray = meta.subarray.name

    # print(apname, apname_pps, subarray)

    # No need to do anything if the aperture names are the same
    # Also skip if MASK not in apname_pps
    if ((apname==apname_pps) or ('MASK' not in apname_pps)) and ('400X256' not in subarray):
        apname_new = apname
    else:
        # Should only get here if coron mask and apname doesn't match PPS
        apname_str_split = apname.split('_')
        sca = apname_str_split[0]
        image_mask = get_mask_from_pps(apname_pps)

        # Get subarray info
        # Sometimes apname erroneously has 'FULL' in it
        # So, first for subarray info in apname_pps 
        if ('400X256' in apname_pps) or ('400X256' in subarray):
            apn0 = f'{sca}_400X256'
        elif ('FULL' in apname_pps):
            apn0 = f'{sca}_FULL'
        else:
            apn0 = sca

        apname_new = f'{apn0}_{image_mask}'

        # Append filter or NARROW if needed
        pps_str_arr = apname_pps.split('_')
        last_str = pps_str_arr[-1]
        # Look for filter specified in PPS aperture name
        if ('_F1' in apname_pps) or ('_F2' in apname_pps) or ('_F3' in apname_pps) or ('_F4' in apname_pps):
            # Find all instances of "_"
            inds = [pos for pos, char in enumerate(apname_pps) if char == '_']
            # Filter is always appended to end, but can have different string sizes (F322W2)
            filter = apname_pps[inds[-1]+1:]
            apname_new += f'_{filter}'
        elif last_str=='NARROW':
            apname_new += '_NARROW'
        elif ('TAMASK' in apname_pps) and ('WB' in apname_pps[-1]):
            apname_new += '_WEDGE_BAR'
        elif ('TAMASK' in apname_pps) and (apname_pps[-1]=='R'):
            apname_new += '_WEDGE_RND'

    # print(apname_new)

    # If apname_new doesn't exit, we need to fall back to apname
    # even if it may not fully make sense.
    if apname_new in pysiaf.Siaf('NIRCam').apernames:
        return apname_new
    else:
        return apname


def get_mask_from_pps(apname_pps):
    """Get mask name from PPS aperture name
    
    The PPS aperture name is of the form:
        NRC[A/B][1-5]_[FULL]_[MASK]
    where MASK is the name of the coronagraphic mask used.

    For target acquisition apertures the mask name can be
    prependend with "TA" (eg., TAMASK335R).

    Return '' if MASK not in input aperture name.
    """

    if 'MASK' not in apname_pps:
        return ''

    pps_str_arr = apname_pps.split('_')
    for s in pps_str_arr:
        if 'MASK' in s:
            image_mask = s
            break

    # Special case for TA apertures
    if 'TA' in image_mask:
        # Remove TA from mask name
        image_mask = image_mask.replace('TA', '')

        # Remove FS from mask name
        if 'FS' in image_mask:
            image_mask = image_mask.replace('FS', '')

        # Remove trailing S or L from LWB and SWB TA apertures
        if ('WB' in image_mask) and (image_mask[-1]=='S' or image_mask[-1]=='L'):
            image_mask = image_mask[:-1]

    return image_mask
