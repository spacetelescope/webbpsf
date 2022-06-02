# Functions for interacting with the MAST archive and JWST measured data


import os
from astroquery.mast import Mast
from astropy import table
from astropy.time import Time, TimeDelta
import astropy.time
import astropy.io.fits as fits
import webbpsf.utils

### Login and authentication

service = 'Mast.Jwst.Filtered.Wss'
mast_login_ok = None

def mast_wss_login():
    """Login to MAST via API token, for file downloads

    """
    global mast_login_ok

    if mast_login_ok:
        return

    mast_api_token = os.environ.get('MAST_API_TOKEN')

    if mast_api_token is None:
        raise RuntimeError("You must define the MAST_API_TOKEN environment variable to use this function. See https://auth.mast.stsci.edu/info")
    Mast.login(token=mast_api_token)
    mast_login_ok = True


def mast_retrieve_opd(filename, verbose=False, redownload=False):
    """Download an OPD from MAST. Files are saved in the WebbPSF data folder.
    If file is already present locally, the download is skipped and the cached file is used.
    """

    output_path = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'MAST_JWST_WSS_OPDs')
    output_filename = os.path.join(output_path, filename)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    mast_wss_login()

    if os.path.exists(output_filename) and not redownload:
        if verbose:
            print(f"Found OPD file previously downloaded: {filename}")
        return output_filename

    data_uri = f'mast:JWST/product/{filename}'

    # Download the file
    url_path = Mast._portal_api_connection.MAST_DOWNLOAD_URL
    Mast._download_file(url_path + "?uri=" + data_uri, output_filename);

    return output_filename

### Functions for searching and retrieving OPDs based on time


def mast_wss_date_query(date, tdelta):
    """Search for OPDs within a specified range of a given date"""

    t_start, t_stop = date-tdelta,date+tdelta

    params = {"columns":"*",
              "filters":[{"paramName":"ap_type","values":["OPD"]},
                         {"paramName":"date_obs_mjd","values":[{"min":t_start.mjd,"max":t_stop.mjd}]}]}

    return Mast.service_request(service, params)

def mast_wss_opds_around_date_query(date, verbose=True):
    """Retrieve OPDs preceding and following a given date

    returns: tuple of two URIs for the data before and after.

    """


    if not isinstance(date, Time):
        raise ValueError("Please supply the date as an astropy.time.Time instance")


    # Set date range (units of days) for the query
    # Note: start with a small value (+-1.5 day) so the MAST query doesn't start off too large
    # This is consistent with expected WFS cadence
    tdelta = TimeDelta(1.5, format='jd')

    # With a too-small date range, this initial query may return a "NoResultsWarning"
    obs_table = mast_wss_date_query(date, tdelta)

    # If the initial query:
    # - returns no results OR
    # - does not include an OPD that precedes the given date OR
    # - does not include an OPD that follows the given date
    # Run the query again with a larger date range
    while len(obs_table) < 1 or min(obs_table['date_obs_mjd']) > date.mjd or max(obs_table['date_obs_mjd']) < date.mjd:
        tdelta *= 2
        if verbose:
            print(f"iterating query, tdelta={tdelta}")

        if tdelta > 14:
            print("Insufficient JWST OPD data can be found within 1 week of the specified date. This date is likely outside of the valid range")
            raise RuntimeError("Cannot find OPDs in MAST before/after that date. Date is likely outside the range of valid data.")
        obs_table = mast_wss_date_query(date, tdelta)

    if verbose:
        print(f'\nMAST OPD query around UTC: {date}')
        print(f'                        MJD: {date.mjd}')

    # In case you provide a datetime that exactly matches the datetime of an OPD file
    if obs_table[date.mjd-obs_table['date_obs_mjd'] == 0]:
        current_opd = obs_table[obs_table['date_obs_mjd']-date.mjd == 0] # Get files with date_obs_mjd == provided datetime
        if verbose:
            print('\nThe given datetime *exactly* matches the datetime of an OPD file:')
            print(f'URI -- {current_opd[0]["dataURI"]}')
            print(f'Date (MJD) -- {current_opd[0]["date_obs_mjd"]}')

        return (current_opd[0]["fileName"], current_opd[0]["fileName"])

    # Otherwise, print some details about the immediately preceding and following OPD files
    else:
        temp_table = obs_table[ obs_table['date_obs_mjd'] - date.mjd < 0]  # Get files with date_obs_mjd < provided datetime
        prev_opd = temp_table[[obs['date_obs_mjd'] == max(temp_table['date_obs_mjd']) for obs in temp_table]]
        temp_table = obs_table[ obs_table['date_obs_mjd'] - date.mjd > 0]  # Get files with date_obs_mjd > provided datetime
        next_opd = temp_table[[obs['date_obs_mjd'] == min(temp_table['date_obs_mjd']) for obs in temp_table]]

        if verbose:
            print('\nOPD immediately preceding the given datetime:')
            print(f'\tURI:\t {prev_opd[0]["dataURI"]}')
            print(f'\tDate (MJD):\t {prev_opd[0]["date_obs_mjd"]:.4f}')
            print(f'\tDelta time:\t {prev_opd[0]["date_obs_mjd"]-date.mjd:.4f} days')

            print('\nOPD immediately following the given datetime:')
            print(f'\tURI:\t {next_opd[0]["dataURI"]}')
            print(f'\tDate (MJD):\t {next_opd[0]["date_obs_mjd"]:.4f}')
            print(f'\tDelta time:\t {next_opd[0]["date_obs_mjd"]-date.mjd:.4f} days')

    return (prev_opd[0]["fileName"], next_opd[0]["fileName"])


def get_opd_at_time(date, choice='before', verbose=False):
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

    prev_opd_fn, post_opd_fn = mast_wss_opds_around_date_query(date, verbose=verbose)

    if choice== 'before':
        if verbose: print(f"User requested choosing OPD before date {date}, which is {prev_opd_fn}")
        fn = mast_retrieve_opd(prev_opd_fn)
        return fn
    elif choice== 'after':
        if verbose: print(f"User requested choosing OPD after date {date}, which is {post_opd_fn}")
        fn = mast_retrieve_opd(post_opd_fn)
        return fn
    elif choice== 'average':
        if verbose: print(f"User requested calculating OPD time averaged around {date}")
        fn_pre = mast_retrieve_opd(pre_opd_fn)
        fn_post = mast_retrieve_opd(post_opd_fn)
        raise NotImplementedError("Not yet implemented")
    elif choice== 'closest':
        if verbose: print(f"User requested choosing OPD time closest in time to {date}")
        raise NotImplementedError("Not yet implemented")



### Functions for format conversion of OPDs

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

    inputOPD = wasopd["RESULT_PHASE"].data
    npix_in = inputOPD.shape[0]

    wasopd[0].header.add_history("Converted for use with WebbPSF:")
    wasopd[0].header.add_history(f"  Converting input file {filename:s} ")
    wasopd[0].header.add_history(f"  from {npix_in:d}x{npix_in:d} to {npix_out:d}x{npix_out:d}.")

    if verbose:
        print(f"Converting {filename:s} from {npix_in:d}x{npix_in:d} to {npix_out:d}x{npix_out:d}")

    # First, pad/dilate the OPD to fill in invalid pixels (0 value) adjacent to valid pixels.
    #  We do this before interpolating to minimize edge effects from the
    #  initial coarse resolution on the segment gaps
    mask = inputOPD != 0
    paddedOPD = webbpsf.utils.border_extrapolate_pad(inputOPD, mask)
    wasopd[0].header.add_history(f"  Dilated OPD values to fill adjacent invalid pixels (i.e. fill in gaps)")

    # interpolate to larger size
    newopd = webbpsf.utils.rescale_interpolate_opd(paddedOPD, npix_out)
    wasopd[0].header.add_history(f"  Interpolated array to {npix_out:d}x{npix_out:d} pixels across")

    # Convert units from microns to meters (as expected by poppy)
    # WSS output files are in units of microns, though this is not given in the headers.
    newopd *= 1e-6
    wasopd[0].header.add_history("Converted units from microns to meters")

    # Update FITS header
    wasopd[0].header["BUNIT"] = 'm'
    wasopd[0].header["PUPLDIAM"] = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER
    wasopd[0].header["PUPLSCAL"] = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER / npix_out

    # WebbPSF expects OPDs in the 0th extension, not 1st, so copy the values there too
    wasopd[0].data = newopd

    # The WSS puts many addtional phase retrieval products in later OPDs, which webbpsf doesn't need.
    while (len(wasopd) > 2):
        del wasopd[2]

    return wasopd
