# Functions for interacting with the MAST archive and JWST measured data


import functools
import os

import astropy
import astropy.io.fits as fits
import astropy.time
import astropy.units as u
import numpy as np
from astropy.time import Time, TimeDelta
from astroquery.mast import Mast, Observations

import webbpsf.utils

# Login and authentication

service = 'Mast.Jwst.Filtered.Wss'


def mast_retrieve_opd(filename, output_path=None, verbose=False, redownload=False):
    """Download an OPD from MAST. Files are saved in the WebbPSF data folder.
    If file is already present locally, the download is skipped and the cached file is used.
    """

    if output_path is None:
        output_path = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'MAST_JWST_WSS_OPDs')
    else:
        output_path = output_path

    output_filename = os.path.join(output_path, filename)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    if os.path.exists(output_filename) and not redownload:
        if verbose:
            print(f'Found OPD file previously downloaded: {filename}')
        return output_filename

    data_uri = f'mast:JWST/product/{filename}'

    # Download the file
    url_path = Mast._portal_api_connection.MAST_DOWNLOAD_URL
    Mast._download_file(url_path + '?uri=' + data_uri, output_filename)
    return output_filename


def download_all_opds(opdtable, verbose=False):
    """Download all OPDs included in some table."""

    for row in opdtable:
        webbpsf.mast_wss.mast_retrieve_opd(row['fileName'], verbose=verbose)


# Functions for searching and retrieving OPDs based on time


def mast_wss_date_query(date, tdelta):
    """Search for OPDs within a specified range of a given date"""

    t_start, t_stop = date - tdelta, date + tdelta

    params = {
        'columns': '*',
        'filters': [
            {'paramName': 'ap_type', 'values': ['OPD']},
            {'paramName': 'date_obs_mjd', 'values': [{'min': t_start.mjd, 'max': t_stop.mjd}]},
        ],
    }

    return Mast.service_request(service, params)


def mast_wss_opds_around_date_query(date, verbose=True):
    """Retrieve OPDs preceding and following a given date

    returns: tuple of two URIs for the data before and after, followed by two fractional date offsets
        i.e. it's one tuple with two string filenames, follows by two floats for how many days before/after

    """

    if not isinstance(date, Time):
        date = Time(date)
        # raise ValueError("Please supply the date as an astropy.time.Time instance")

    # Set date range (units of days) for the query
    # Note: start with a small value (+-1.5 day) so the MAST query doesn't start off too large
    # This is consistent with expected WFS cadence
    tdelta = TimeDelta(1.5 * u.day, format='jd')

    # With a too-small date range, this initial query may return a "NoResultsWarning"
    obs_table = mast_wss_date_query(date, tdelta)

    nfound = 0
    # If the initial query:
    # - returns no results OR
    # - does not include an OPD that precedes the given date OR
    # - does not include an OPD that follows the given date
    # Run the query again with a larger date range
    while len(obs_table) < 1 or min(obs_table['date_obs_mjd']) > date.mjd or max(obs_table['date_obs_mjd']) < date.mjd:
        if tdelta >= 6 * u.day:
            if verbose:
                print(
                    'Could not find JWST OPDs both before and after the specified date.',
                    'Date outside of the available range of WFS data.'
                )

            if len(obs_table) == 0:
                error_message = (
                    'Cannot find ANY OPDs in MAST within a week before/after that date.'
                    'Date is likely outside the range of valid data.'
                )
                raise RuntimeError(error_message)
            elif max(obs_table['date_obs_mjd']) < date.mjd:
                # if len(obs_table) == 1 : #and min(obs_table['date_obs_mjd']) < date.mjd:
                if verbose:
                    print('Found at least one OPD before that date, but no OPD after that date.')
                closest = [np.argmin(np.abs(obs_table['date_obs_mjd'] - date.mjd))]
                obs_table = obs_table[closest]

                nfound = 1
                break

        tdelta *= 2
        if verbose:
            print(f'iterating query, tdelta={tdelta}')

        obs_table = mast_wss_date_query(date, tdelta)
        nfound = 2

    if verbose:
        print(f'\nMAST OPD query around UTC: {date}')
        print(f'                        MJD: {date.mjd}')

    # In case we only found one file within the searched date range. This will most often be the case if searching for
    # an OPD for the very most recent observations, for which there may not yet be any "after" measurement.
    if nfound == 1:
        current_opd = obs_table
        if verbose:
            print('\nOnly found one OPD file when searching  :')
            print(f'URI -- {current_opd[0]["dataURI"]}')
            print(f'Date (MJD) -- {current_opd[0]["date_obs_mjd"]}')

        return (current_opd[0]['fileName'], 'Not found', current_opd[0]['date_obs_mjd'] - date.mjd, np.nan)

    # In case you provide a datetime that exactly matches the datetime of an OPD file
    elif obs_table[date.mjd - obs_table['date_obs_mjd'] == 0]:
        current_opd = obs_table[
            obs_table['date_obs_mjd'] - date.mjd == 0
        ]  # Get files with date_obs_mjd == provided datetime
        if verbose:
            print('\nThe given datetime *exactly* matches the datetime of an OPD file:')
            print(f'URI -- {current_opd[0]["dataURI"]}')
            print(f'Date (MJD) -- {current_opd[0]["date_obs_mjd"]}')

        return (
            current_opd[0]['fileName'],
            current_opd[0]['fileName'],
            current_opd[0]['date_obs_mjd'] - date.mjd,
            current_opd[0]['date_obs_mjd'] - date.mjd,
        )

    # Otherwise, print some details about the immediately preceding and following OPD files
    else:
        temp_table = obs_table[obs_table['date_obs_mjd'] - date.mjd < 0]  # Get files with date_obs_mjd < provided datetime
        prev_opd = temp_table[[obs['date_obs_mjd'] == max(temp_table['date_obs_mjd']) for obs in temp_table]]
        temp_table = obs_table[obs_table['date_obs_mjd'] - date.mjd > 0]  # Get files with date_obs_mjd > provided datetime
        next_opd = temp_table[[obs['date_obs_mjd'] == min(temp_table['date_obs_mjd']) for obs in temp_table]]

        if verbose:
            print('\nOPD immediately preceding the given datetime:')
            print(f'\tURI:\t {prev_opd[0]["dataURI"]}')
            print(f'\tDate (MJD):\t {prev_opd[0]["date_obs_mjd"]:.4f}')
            print(f'\tDelta time:\t {prev_opd[0]["date_obs_mjd"] - date.mjd:.4f} days')

            print('\nOPD immediately following the given datetime:')
            print(f'\tURI:\t {next_opd[0]["dataURI"]}')
            print(f'\tDate (MJD):\t {next_opd[0]["date_obs_mjd"]:.4f}')
            print(f'\tDelta time:\t {next_opd[0]["date_obs_mjd"] - date.mjd:.4f} days')

    return (
        prev_opd[0]['fileName'],
        next_opd[0]['fileName'],
        prev_opd[0]['date_obs_mjd'] - date.mjd,
        next_opd[0]['date_obs_mjd'] - date.mjd,
    )


def get_opd_at_time(date, choice='closest', verbose=False, output_path=None):
    """Get an estimated OPD at a given time based on measured OPDs from JWST wavefront sensing monitoring

    Parameters
    ----------
    date : string or astropy.time.Time instance
        Datetime, either as a string giving ISO time in UTC, e.g. 2021-12-25T11:20:00, or an astropy Time instance.
    choice : string
        How to choose which OPD to use. Allowable values include the following:
        - 'before': use OPD measured shortly prior to the specified time
        - 'closest': use OPD measured closest to the specified time
        - 'after': use OPD measured shortly afterto the specified time
        - 'average': use a weighted average of the before and after OPDs, weighted based on proximity in time.

    Returns
    --------
    Filename for the retrieved OPD, or FITS HDUList instance if an average is to be used (TBC).
    """

    if isinstance(date, str):
        date = astropy.time.Time(date, format='isot')

    prev_opd_fn, post_opd_fn, prev_dtime, post_dtime = mast_wss_opds_around_date_query(date, verbose=verbose)

    if choice == 'before':
        if verbose:
            print(
                f'User requested choosing OPD before date {date}, which is {prev_opd_fn}, delta time {prev_dtime:.3f} days'
            )
        return mast_retrieve_opd(prev_opd_fn, output_path=output_path)
    elif choice == 'after':
        if verbose:
            print(f'User requested choosing OPD after date {date}, which is {post_opd_fn}, delta time {post_dtime:.3f} days')
        return mast_retrieve_opd(post_opd_fn, output_path=output_path)
    elif choice == 'average':
        if verbose:
            print(f'User requested calculating OPD time averaged around {date}')
        mast_retrieve_opd(prev_opd_fn, output_path=output_path)
        mast_retrieve_opd(post_opd_fn, output_path=output_path)
        raise NotImplementedError('Not yet implemented')
    elif choice == 'closest':
        closest_fn, closest_dt = (
            (post_opd_fn, post_dtime) if abs(post_dtime) < abs(prev_dtime) else (prev_opd_fn, prev_dtime)
        )
        if verbose:
            print(
                f'User requested choosing OPD time closest in time to {date},',
                f'which is {closest_fn}, delta time {closest_dt:.3f} days'
            )
        return mast_retrieve_opd(closest_fn, output_path=output_path)


# Functions for format conversion of OPDs


def import_wss_opd(filename, npix_out=1024, verbose=False):
    """Import an OPD produced by the JWST WSS, and convert to the right format for use with WebbPSF.

    This includes:
    - Rescale from the input size, probably 256x256, up to 1024x1024 (or any other requested size).
    - To reduce interpolation artifacts along edges, before interpolating, pad OPD values along segment edges
    - Update FITS header to include keywords as needed for WebbPSF
    - Copy OPD values from the 1st extension (where the WSS puts it) to the 0th extension (where webbpsf wants it)

    Note, this function does NOT subtract off the SI WFE portion; see "load_wss_opd" for that.

    Parameters
    ----------
    filename : string
        Filename for input OPD file, optionally including path. Should resolve to an accessible local path,
        i.e. this file should be already downloaded on disk.
    npix_out : int
        Number of pixels per side in the converted output OPD


    Returns
    -------
    astropy.fits.HDUList instance for the converted OPD.

    """

    wasopd = fits.open(filename)

    inputOPD = wasopd['RESULT_PHASE'].data
    npix_in = inputOPD.shape[0]

    wasopd[0].header.add_history('OPD file retrieved from MAST for use by WebbPSF.')

    wasopd[0].header.add_history('Converting input for use with WebbPSF:')
    wasopd[0].header.add_history(f'  Input file is: {filename:s} ')
    wasopd[0].header.add_history(f'  Need to rescale from {npix_in:d}x{npix_in:d} to {npix_out:d}x{npix_out:d} pixels.')

    if verbose:
        print(f'Converting {filename:s} from {npix_in:d}x{npix_in:d} to {npix_out:d}x{npix_out:d}')

    # First, pad/dilate the OPD to fill in invalid pixels (0 value) adjacent to valid pixels.
    #  We do this before interpolating to minimize edge effects from the
    #  initial coarse resolution on the segment gaps
    mask = inputOPD != 0
    paddedOPD = webbpsf.utils.border_extrapolate_pad(inputOPD, mask)
    wasopd[0].header.add_history('  Dilated OPD values to fill adjacent invalid pixels (i.e. fill in gaps)')

    # interpolate to larger size
    newopd = webbpsf.utils.rescale_interpolate_opd(paddedOPD, npix_out)
    wasopd[0].header.add_history(f'  Interpolated array to {npix_out:d}x{npix_out:d} pixels across')

    # Convert units from microns to meters (as expected by poppy)
    # WSS output files are in units of microns, though this is not given in the headers.
    newopd *= 1e-6
    wasopd[0].header.add_history('  Converted units from microns to meters')

    # Update FITS header
    wasopd[0].header['BUNIT'] = 'm'
    wasopd[0].header['PUPLDIAM'] = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER
    wasopd[0].header['PUPLSCAL'] = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER / npix_out

    # WebbPSF expects OPDs in the 0th extension, not 1st, so copy the values there too
    wasopd[0].data = newopd

    # The WSS puts many addtional phase retrieval products in later OPDs, which webbpsf doesn't need.
    while len(wasopd) > 2:
        del wasopd[2]

    return wasopd


# Functions for dealing with time series or entire set of OPDs


def infer_pre_or_post_correction(row):
    """Use the activity label to infer if a given WFS measurement OPD is from pre or post correction"""

    act = row['activity']

    # handle some special cases
    if row['visitId'] == 'V01163111001':
        # replacement for obs 10 which couldn't get into Track
        return 'post'  # WFE was too bad to guide, so we applied correction using WFSC Commissioning in Coarse,
        # Then took this weak lens data after.
    elif row['visitId'].startswith('V01445'):
        # thermal slew cold attitude WFS is all sensing-only
        return 'pre'

    # infer based on activity label
    lookup = {
        '02101': 'pre',  # wfsc only, any diversity, pre move (common, early)
        '02104': 'post',  # wfsc only, diversity PM8, post move
        '02106': 'F187N pre',  # wfsc only, diversity ALL+187N, pre move (rare, post MIMF correction)
        '02107': 'post',  # wfsc only, diversity ALL, post move
        '02109': 'post',  # wfsc only, diversity ALL+187N, post move
        '0210E': 'F187N post',  # wfsc only, diversity ALL+187N, post move, using F187N
        # visits with jitter sensing included:
        '03104': 'pre',  # visit with jitter sensing, any diversity, pre move
        '03107': 'post',  # visit with jitter sensing, diversity PM8, post move
        '03109': 'F187N pre',  # visit with jitter sensing, diversity ALL+187N, pre move, using F187N?
        '0310A': 'post',  # visit with jitter sensing, diversity ALL, post move
        '0310C': 'post',  # visit with jitter sensing, diversity ALL+187N, post move
        '0310H': 'F187N post',  # visit with jitter sensing, diversity ALL+187N, post move, using F187N
    }
    if act in lookup:
        return lookup[act]
    elif act.startswith('04'):
        # rare special case which has a genwaitmain, e.g. 1163:206
        # which increments the group by 1
        return lookup['03' + act[2:]]
    else:
        return 'UNKNOWN'


def retrieve_mast_opd_table(aperture_list=['NRCA3_FP1', 'NRCA1_FP6'], verbose=False):
    """Retrieve table of OPDs from MAST.

    Returns : Astropy table listing available OPDs and metadata such as dates and sensing type.

    """

    # Construct the query and execute the search to retrieve available OPDs in MAST
    params = {
        'columns': '*',
        'filters': [{'paramName': 'ap_type', 'values': ['OPD']}, {'paramName': 'apername', 'values': aperture_list}],
    }
    obs_table = Mast.service_request(service, params)

    if verbose:
        print('\n\nTotal products with apername = {}: {}'.format(aperture_list, len(obs_table)))

    # Now perform some manipulations on the result, to select and add additional columns
    colnames_we_want = get_colnames()
    columns_we_want = [obs_table[colname] for colname in colnames_we_want]
    # insert a column with times as ISO time strings, instead of MJD
    columns_we_want.insert(
        0, astropy.table.Column(astropy.time.Time(obs_table['date_obs_mjd'], format='mjd').isot, name='date')
    )
    # insert a column with just the OSS activity label, e.g. '02101'
    columns_we_want.insert(3, astropy.table.Column([a[-5:] for a in obs_table['obs_id']], name='activity'))

    opdtable = astropy.table.Table(columns_we_want)
    # Update the visit ID to all start with a prepended initial letter V
    opdtable['visitId'] = ['V' + vid for vid in opdtable['visitId']]
    opdtable.sort('date')

    opdtable = add_columns_to_track_corrections(opdtable)

    return opdtable


def get_colnames(updated=False):
    """Get list of column names from MAST that are desired to populate the OPD Astropy table
    Parameters
    ----------
    updated : boolean
        toggle whether or not to include the full colnames list with columns to track corrections
        these columns are added via add_columns_to_track_corrections
    Returns : list of colnames from MAST that will be used in the opdtable
    """
    colnames = ['date_obs_mjd', 'visitId', 'apername', 'corr_id', 'fileName', 'dataURI']
    if updated:
        colnames.append('wfs_measurement_type')
        colnames.append('is_post_correction')
        colnames.append('is_pre_correction')

    return colnames


def add_columns_to_track_corrections(opdtable):
    """Add useful columns which help track when there were corrections

    Parameters
    ----------
    opdtable : astropy.table.Table
        Table of available OPDs, Default None: as returned by retrieve_mast_opd_table()
    Returns : Updated Astropy table listing available OPDs and metadata such as dates and sensing type.

    Additional columns:
    wfs_measurement_type : pre | post  (mirror moves)
    is_post_correction : whether it is after a correction has been applied
    is_pre_correction : whether it is prior to a correction

    """
    #
    pre_or_post = []

    for row in opdtable:
        pre_or_post.append(infer_pre_or_post_correction(row))

    # Add column for is this WFS measurement made immediately after a correction
    opdtable['wfs_measurement_type'] = pre_or_post
    opdtable['is_post_correction'] = ['post' in a for a in pre_or_post]

    # add column for is this a WFS measurement made right before a correction,
    # which has a corresponding 2nd measurement right after
    has_correction = []
    for row in opdtable:
        if row['is_post_correction']:
            # if this is a measurement after a correction, it itself doesn't have another correction after it
            has_correction.append(False)
        else:
            # Find if there is a row with same visit ID, and which is post correction
            matching_post_correction = (opdtable['visitId'] == row['visitId']) & opdtable['is_post_correction']
            has_correction.append(matching_post_correction.sum() > 0)
    opdtable['is_pre_correction'] = has_correction

    return opdtable


def deduplicate_opd_table(opdtable, drop_f187n=True, verbose=False):
    """Filter out duplicates from the OPD table, to achieve a list of unique OPD measurements.

    This is needed because, in some cases, there are multiple analyses for the same input data, so
    the list of all OPDs has some redundancies. This function filters those out so there's at most
    one unique OPD per each input WFS dataset.

    For cases with multi-wavelength WFS, this also filters out the F187N results, similarly in order to
    result in one OPD per input sensing instance.

    Parameters
    -----------
    opdtable : astropy.table.Table instance
        OPD table returned from retrieve_mast_opd_table
    drop_f187n : bool
        Filter out the F187N data?
    verbose : bool
        Print more verbose output?


    """
    # Remove redundant APs
    # Also discard any which are F187N

    # Find where we should start; let's only look at APs after OTE alignment is complete
    # mimf2_corection_index = np.where([v.startswith('V01467') for v in opdtable['visitId']])[0].max()
    # if verbose:
    #    print("Retaining only OPDs after OTE alignment complete.")

    measurement_dates_encountered = set()
    redundant_aps = set()
    # Some APs are known to be not good, due to issues in the analyses
    invalid_aps = [
        'R2022083106',
    ]

    unique_indices = []

    # Some datetimes have multiple redundant APs
    # work through the list in forward order to figure out which such cases to ignore
    # Update: do it in backwards order, so that the last/latest AP for a particular observation
    # is considered the most authoritative.
    # all_row_indices = list(range(mimf2_corection_index, len(opdtable)))
    all_row_indices = list(range(len(opdtable)))
    all_row_indices.reverse()
    for row_index in all_row_indices:
        if drop_f187n and ('F187N' in opdtable[row_index]['wfs_measurement_type']):
            if verbose:
                print(f"{opdtable[row_index]['fileName']} is F187N. Ignoring it.")
            continue

        if opdtable[row_index]['corr_id'] in invalid_aps:
            if verbose:
                print(
                    f"{opdtable[row_index]['fileName']}",
                    "is flaggd as an invalid AP due to known analysis issue(s). Ignoring it."
                )
            continue
        elif opdtable[row_index]['date'] in measurement_dates_encountered:
            # sometimes we end up with multiple redundant APs per a given measurement
            # If so, only plot the first one, and ignore any redundant ones.
            if verbose:
                print(f"{opdtable[row_index]['fileName']} is redundant for a measurement in a prior AP. Ignoring it.")
            redundant_aps.add(opdtable[row_index]['fileName'])
            continue
        else:
            # We should use this AP
            unique_indices.append(row_index)
            measurement_dates_encountered.add(opdtable[row_index]['date'])

    # Because we went through the table in reverse order, un-reverse the list of indices
    unique_indices.reverse()

    return opdtable[unique_indices]


def filter_opd_table(opdtable, start_time=None, end_time=None, wfs_measurement_type=None):
    """Filter an OPD table, for instance to select a date range

    Parameters
    ----------
    opdtable : astropy.table.Table instance
        Table of OPDs, from retrieve_mast_opd_table
    start_time, end_time : strings or astropy.Time.Time instances
        Define start and/or end time for selecting a date range
    wfs_measurement_type : string
        WFS measurement type, e.g. "pre" or "post" for pre-correction or post correction, etc

    Returns a filtered copy of the input table.
    """
    valid_mask = np.array(len(opdtable), bool)
    if start_time:
        if isinstance(start_time, str):
            start_time = astropy.time.Time(start_time)
        valid_mask = valid_mask & (opdtable['date_obs_mjd'].value > start_time.mjd)
    if end_time:
        if isinstance(end_time, str):
            end_time = astropy.time.Time(end_time)
        valid_mask = valid_mask & (opdtable['date_obs_mjd'].value < end_time.mjd)
    if wfs_measurement_type:
        valid_mask = valid_mask & (opdtable['wfs_measurement_type'] == wfs_measurement_type)

    return opdtable[valid_mask]


def get_corrections(opdtable):
    """Identify all mirror corrections from a given set of OPDs

    Parameters
    ----------
    opdtable : astropy.table.Table instance
        Table of OPDs, as returned from

    Returns: Astropy table listing mirror corrections and the associated WFS measurements before/after each

    """

    # Iterate over the table to identify which WFS visits included corrections.

    # Match up pairs of WFS measurements on the same visit. These indicate visits with corrections.
    pre_correction_indices = []
    post_correction_indices = []

    corr_ids = []

    for i, row in enumerate(opdtable):
        if not row['is_pre_correction']:
            continue

        # Check if there is a matching post correction OPD from the same visit?
        w = (opdtable['visitId'] == row['visitId']) & (opdtable['is_post_correction'])
        if sum(w) > 0:
            # If we find something, yes this visit had a correction.
            pre_correction_indices.append(i)
            post_correction_indices.append(np.where(w)[0][0])

            # infer the correction ID
            # based on the assumption that the move probably came from the immediately prior WSS session
            prior_corr_id = opdtable[i - 1]['corr_id']
            corr_ids.append(prior_corr_id)

    correction_table = astropy.table.Table(
        [
            corr_ids,
            opdtable[pre_correction_indices]['visitId'],
            opdtable[pre_correction_indices]['date'],
            opdtable[pre_correction_indices]['fileName'],
            opdtable[post_correction_indices]['date'],
            opdtable[post_correction_indices]['fileName'],
        ],
        names=[
            'WFC ID',
            'Visit ID',
            'Pre Move Sensing Time',
            'Pre Move Sensing OPD Filename',
            'Post Move Sensing Time',
            'Post Move Sensing OPD Filename',
        ],
    )

    return correction_table


# Functions for retrieving images


def set_params(parameters):
    """Helper function for setting up the dicts needed to call Mast.service_request"""
    return [{'paramName': p, 'values': v} for p, v in parameters.items()]


def query_wfsc_images_by_program(prog, obs, detector='NRCA3', productlevel=3):
    """MAST query to get filenames of WFSC images from a given program and observation

    Returns table of MAST query results
    """
    # Query MAST for level-3 WFSC combined products
    keywords = {
        'program': [str(prog)],
        'observtn': [str(obs)],
        'detector': [detector],
        'productLevel': [str(productlevel)],
    }

    # Restructuring the keywords dictionary to the MAST syntax
    params = {'columns': '*', 'filters': set_params(keywords)}
    service = 'Mast.Jwst.Filtered.Nircam'
    tab = Mast.service_request(service, params)
    return tab


def query_wfsc_images_latest(detector='NRCA3', productlevel=3):
    """MAST query to get filenames of the most recent WFSC images

    Returns table of MAST query results
    """
    # Query MAST for the most recent level-3 WFSC combined products
    now = astropy.time.Time.now()
    keywords = {
        'template': ['WFSC NIRCam Fine Phasing'],
        'detector': [detector],
        'date_obs_mjd': [{'min': now.mjd - 7, 'max': now.mjd}],
        'productLevel': [str(productlevel)],
    }

    # Restructuring the keywords dictionary to the MAST syntax
    params = {
        'columns': 'filename, apername, program, observtn, visit_id, act_id,  date_obs_mjd, date_obs, productLevel',
        'filters': set_params(keywords),
    }
    service = 'Mast.Jwst.Filtered.Nircam'
    tab = Mast.service_request(service, params)
    tab.sort(keys=['date_obs_mjd'], reverse=True)

    # some sanity checking:
    if tab[0]['visit_id'] != tab[1]['visit_id']:
        raise RuntimeError('The latest two images in MAST are not from the same visit')
        # these are indeed two files from the same visi
    else:
        return tab[0:2]  # Just return that most recent pair


def download_wfsc_images(program=None, obs=None, verbose=False, **kwargs):
    """Download to current directory some WFSC weak lens images.
    Either specify a WFSC program and observation number, or else by default
    the most recent WFSC images available will be retrieved.

    Parameters
    ----------
    program, obs : int or str
        Program and observation numbers.
    verbose : bool
        Be more verbose?

    Other kwargs for detector and productlevel are passed through
    to the query functions.

    Returns list of image filenames
    """

    if program is None and obs is None:
        if verbose:
            print('Querying latest available WFSC images')
        filetable = query_wfsc_images_latest(**kwargs)
    else:
        if verbose:
            print(f'Querying WFSC images from program {program}, observation {obs}')
        filetable = query_wfsc_images_by_program(program, obs, **kwargs)

    # If we found > 0 available files for that visit, then we're done searching and can go on to download
    if verbose:
        date_obs = astropy.time.Time(filetable[0]['date_obs_mjd'], format='mjd')
        print(
            f"Found {len(filetable)} level 3 data products from",
            f"{filetable[0]['program']}:{filetable[0]['observtn']} around {date_obs.iso[0:16]}"
        )
    if len(filetable) > 0:
        file_list = []
        for row in filetable:
            file_list.append(row['filename'])
            data_uri = f"mast:JWST/product/{row['filename']}"
            Observations.download_file(data_uri)
        return file_list
    else:
        raise RuntimeError('Error, could not find any WFSC image data matching the specified search parameters.')


@functools.lru_cache
def get_visit_nrc_ta_image(visitid, verbose=True, kind='cal'):
    """Retrieve from MAST the NIRCam target acq image for a given visit.

    This retrieves an image from MAST and returns it as a HDUList variable
    without writing to disk.
    """

    keywords = {
        'visit_id': [visitid.lstrip('V')],  # note: drop the initial character 'V'
        'exp_type': ['NRC_TACQ'],
    }

    # Restructuring the keywords dictionary to the MAST syntax
    params = {'columns': '*', 'filters': set_params(keywords)}

    service = 'Mast.Jwst.Filtered.Nircam'
    t = Mast.service_request(service, params)
    filename = t[0]['filename']

    # If user manually specifies rate or uncal, retrieve that instead
    if kind == 'rate' or kind == 'uncal':
        filename = filename.replace('_cal.fits', f'_{kind}.fits')

    if verbose:
        print(f'TA filename: {filename}')
    mast_file_url = f'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/{filename}'
    import urllib

    try:
        ta_hdul = fits.open(mast_file_url)
    except urllib.error.HTTPError as err:
        if err.code == 401:  # Unauthorized
            # Use MAST API to allow retrieval of exclusive access data, if relevant
            import tempfile

            import astroquery

            mast_api_token = os.environ.get('MAST_API_TOKEN', None)
            mast_obs = astroquery.mast.ObservationsClass(mast_api_token)
            uri = f'mast:JWST/product/{filename}'
            with tempfile.NamedTemporaryFile() as temp:
                mast_obs.download_file(uri, local_path=temp.name, cache=False)
                ta_hdul = astropy.io.fits.open(temp.name)
        else:
            raise  # re-raise any errors other than 401 for permissions denied

    return ta_hdul


# Functions for retrieving metadata about science observations


def _query_program_visit_times_by_inst(program, instrument, verbose=False):
    """Get the start and end times of all completed visits in a program, per instrument.
    Not intended for general use; this is mostly a helper to query_program_visit_times.

    Getting the vststart_mjd and visitend_mjd fields requires using the instrument keywords
    interface, so one has to specify which instrument ahead of time.

    Parameters
    ----------
    program : int or str
        Program ID
    instrument : str
        instrument name
    verbose : bool
        be more verbose in output?

    returns list of (visitid, start, end) tuples.

    """

    from astroquery.mast import Mast

    svc_table = {
        'MIRI': 'Mast.Jwst.Filtered.Miri',
        'NIRCAM': 'Mast.Jwst.Filtered.NIRCam',
        'NIRSPEC': 'Mast.Jwst.Filtered.NIRSpec',
        'NIRISS': 'Mast.Jwst.Filtered.NIRISS',
    }

    service = svc_table[instrument.upper()]

    collist = 'filename, program, observtn, visit_id, vststart_mjd, visitend_mjd, bstrtime'
    all_columns = False

    keywords = {
        'program': [
            str(program),
        ]
    }
    parameters = {'columns': '*' if all_columns else collist, 'filters': set_params(keywords)}

    if verbose:
        print('MAST query parameters:')
        print(parameters)

    responsetable = Mast.service_request(service, parameters)
    responsetable.sort(keys='bstrtime')

    visit_times = []

    for row in responsetable:
        visit_times.append(('V' + row['visit_id'], row['vststart_mjd'], row['visitend_mjd']))

    visit_times = set(visit_times)
    return list(visit_times)


def query_program_visit_times(program, verbose=False):
    """Get the start and end times of all completed visits in a program.

    Parameters
    ----------
    program : int or str
        Program ID
    verbose : bool
        be more verbose in output?

    Returns astropy Table with columns for visit ID and start and end times.
    """

    from astroquery.mast import Observations

    obs = Observations.query_criteria(obs_collection=['JWST'], proposal_id=[program])
    # Annoyingly, that query interface doesn't return start/end times
    instruments = set([val.split('/')[0] for val in set(obs['instrument_name'])])

    visit_times = []
    for inst in instruments:
        if verbose:
            print(f'querying for visits using {inst}')
        visit_times += _query_program_visit_times_by_inst(program, inst)

    vids = [v[0] for v in visit_times]
    starts = astropy.time.Time([float(v[1]) for v in visit_times], format='mjd')
    ends = astropy.time.Time([float(v[2]) for v in visit_times], format='mjd')

    # visit_times = np.asarray(visit_times)
    return astropy.table.Table([vids, starts, ends], names=('visit_id', 'start_mjd', 'end_mjd'))
