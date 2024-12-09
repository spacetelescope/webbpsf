import calendar
import functools
import os

import astropy
import astropy.io.fits as fits
import astropy.time
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import poppy
import scipy.interpolate
from astroquery.mast import Observations
from matplotlib.backends.backend_pdf import PdfPages

import webbpsf


def _read_opd(filename, auto_download=True):
    """Trivial utilty function to read OPD from a WSS-output FITS file
    If the file does not exist locally, try to retrieve it from MAST automatically."""
    full_file_path = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'MAST_JWST_WSS_OPDs', filename)
    if not os.path.exists(full_file_path) and auto_download:
        webbpsf.mast_wss.mast_retrieve_opd(filename)
    opdhdu = fits.open(full_file_path)
    opd = opdhdu[1].data.copy()
    return opd, opdhdu


def get_datetime_utc(opdhdul, return_as='string'):
    """Utilty function. Retrieve the date time either as a string with ISO date, or an astropy time object"""
    if return_as == 'string':
        return opdhdul[0].header['DATE-OBS'] + ' ' + opdhdul[0].header['TIME-OBS'][:-4]  # leave off milliseconds
    elif return_as == 'astropy':
        return astropy.time.Time(
            opdhdul[0].header['DATE-OBS'] + 'T' + opdhdul[0].header['TIME-OBS'], format='isot', scale='utc'
        )
    else:
        raise ValueError(f'Invalid value for return_as={return_as}')


def wavefront_time_series_plot(
    opdtable, start_date=None, end_date=None, ymin=0, ymax=250, threshold=80, label_visits=True, label_events=True
):
    """Make a time series plot of total WFS versus time

    Parameters
    ----------
    opdtable : astropy.table.Table
        OPD table, retrieved from MAST. See functons in mast_wss.py
    start_date, end_date : strings or astropy.time.Time objects
        Start and end dates for the plot time range, specified as strings giving the
        ISO format datetime, or astropy time.Time instances. Default is March 2022 to present.
    label_visits : bool
        Label program_id:visit_id for each WFS visit.
    label_events : bool
        Label events along the x-axis in time defined in the events {}

    Returns
    -------
    Nothing, but makes a plot

    """

    # Plot figure showing WFE vs time

    rmses = []

    if 'rms_wfe' in opdtable.colnames:
        rmses = opdtable['rms_wfe']

    dates = astropy.time.Time(opdtable['date'], format='isot')
    pre_or_post = []

    for row in opdtable:
        if os.path.isfile(row['fileName']) is False:
            full_file_path = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'MAST_JWST_WSS_OPDs', row['fileName'])
        else:
            full_file_path = row['fileName']
        if 'rms_wfe' not in opdtable.colnames:
            rmses.append(fits.getheader(full_file_path, ext=1)['RMS_WFE'])
        pre_or_post.append(webbpsf.mast_wss.infer_pre_or_post_correction(row))

    where_pre = ['pre' in a for a in pre_or_post]
    where_post = ['post' in a for a in pre_or_post]

    events = {
        '2022-03-22T22:55:00': ('Cooler State 4', 'red'),
        '2022-04-10T17:24:00': ('Cooler State 6', 'red'),
        '2022-04-07T20:13:00': ('Cooler State 5', 'red'),
        # '2022-04-15T06:33:00': ("NIRISS FPA Heater on (NIS-24)", 'orange'),
        # '2022-04-10T04:30':  ('NIRSpec FPA CCH to level 30','orange'),
        # '2022-04-10T01:00': ('NIRSpc bench & FPA heaters adjust', 'orange'),  # to level 10
        # '2022-04-11T16:00:00': ('FGS trim heater 0 to 4', 'orange'),
        '2022-04-23T06:30:00': ('OTE alignment complete', 'green'),
        '2022-05-12T18:30:00': ('Thermal slew cold attitude start', 'blue'),
        '2022-05-20T06:30:00': ('Thermal slew cold attitude end', 'blue'),
        '2022-05-23T00:00:00': ('Larger micrometeorite strike on C3', 'red'),
        '2022-06-19T18:00:00': ('Coarse move of C3 for astigmatism correction', 'red'),
        # '2022-06-27T00:00:00': ('NIRSpec safing, not in thermal control', 'orange'),
        '2022-07-12T00:00:00': ('Large outlier tilt event on B5+C5', 'orange'),
        '2024-02-25T20:00:00': ('Large outlier wing tilt event on -V2 wing', 'orange'),
    }

    plt.figure(figsize=(16, 12))

    rms_nm = np.asarray(rmses) * 1000

    routine_pids = [
        1163,
        2586,
        2724,
        2725,
        2726,
        4431,  # Commissioning OTE-26 and Cycle 1 maintenance
        4500,
        4501,
        4502,
        4503,
        4504,
        4505,
        4506,
        4507,
        4508,
        4509,
    ]  # Cycle 2
    routine_pids += list(range(6689, 6699))  # Cycle 3

    is_routine = np.asarray([int(v[1:6]) in routine_pids for v in opdtable[where_pre]['visitId']])

    # Plot all, with connecting line
    plt.plot_date(
        dates.plot_date,
        rms_nm,
        ls='-',
        color='gray',
    )
    # Plot maintenance visits
    plt.plot_date(
        dates[where_pre][is_routine].plot_date,
        rms_nm[where_pre][is_routine],
        'o',
        xdate=True,
        label='WF Maintenance Sensing',
        color='C0',
    )
    # Plot other visits (MIMF etc)
    plt.plot_date(
        dates[where_pre][~is_routine].plot_date,
        rms_nm[where_pre][~is_routine],
        'o',
        xdate=True,
        label='Other Sensing (MIMF etc)',
        color='purple',
    )
    # Plot corrections.
    plt.plot_date(
        dates[where_post].plot_date,
        rms_nm[where_post],
        'v',
        xdate=True,
        label='Wavefront Corrections',
        color='C2',
    )

    ax = plt.gca()
    ax.set_ylabel('Observatory WFE (OTE+NRC)\n[nm rms]', fontweight='bold', fontsize=15)
    ax.set_xlabel('Date, UTC', fontweight='bold', fontsize=15)

    # ymin, ymax = 40, 200
    # ymin, ymax = 0, 250
    ax.set_ylim(ymin, ymax)
    ax.axhline(threshold - 10, ls=':', color='gray')
    ax.axhline(threshold, ls=':', color='orange')
    ax.axhline(threshold + 20, ls=':', color='gray')

    if start_date is None:
        start_date = astropy.time.Time('2022-03-20 00:00:00')
    if end_date is None:
        end_date = astropy.time.Time.now() + 7 * u.day
    start_date = astropy.time.Time(start_date)  # cast any input date string to Time
    end_date = astropy.time.Time(end_date)

    ax.set_xlim(start_date.plot_date, end_date.plot_date)

    if end_date - start_date < 365:
        ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(interval=1))
        ax.xaxis.set_minor_locator(matplotlib.dates.DayLocator())
    else:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(interval=1))
        ax.xaxis.set_minor_locator(matplotlib.dates.WeekdayLocator(interval=1))

    ax.tick_params('x', length=10)
    for tick in ax.get_xticklabels():
        tick.set_rotation(75)

    # label events
    if label_events:
        for timestamp, (event, color) in events.items():
            d = astropy.time.Time(timestamp, format='isot')
            # Limit events that happened after start_date:
            if d >= start_date and d <= end_date:
                plt.axvline(d.plot_date, color=color, ls=':', alpha=0.5)
                ax.text(d.plot_date + 0.25, ymax * 0.95, event, color=color, rotation=90, verticalalignment='top',
                        alpha=0.7, clip_on=True)

    # Connect measurements on the same visit
    for row, rms in zip(opdtable[where_post], rms_nm[where_post]):
        w = opdtable[where_pre]['visitId'] == row['visitId']
        if sum(w) > 0:
            d_post = astropy.time.Time(row['date'], format='isot')
            d_pre = astropy.time.Time(opdtable[where_pre][w]['date'][0], format='isot')
            r_post = rms
            r_pre = rms_nm[where_pre][w][0]
            plt.plot([d_pre.plot_date, d_post.plot_date], [r_pre, r_post], color='green')

    # Label visit IDs for routine WFS
    if label_visits:
        for row, rms in zip(opdtable[where_pre][is_routine], rms_nm[where_pre][is_routine]):
            d_pre = astropy.time.Time(row['date'], format='isot')
            visit = row['visitId']
            short_visit = f'  {int(visit[1:6])}:{int(visit[6:9])}'
            if d_pre.datetime >= start_date and d_pre.datetime <= end_date and rms < ymax:
                plt.text(d_pre.plot_date, rms, short_visit, rotation=65, color='C0', fontsize=10, clip_on=True)

    plt.plot(dates.plot_date[where_pre][is_routine], rms_nm[where_pre][is_routine], ls='none', alpha=0.5)
    plt.plot(dates.plot_date[where_post], rms_nm[where_post], ls='none', color='C2', alpha=0.5)

    plt.title('Observatory WFE over time at NIRCam Field Point 1', fontweight='bold', fontsize=15)

    plt.legend(loc='upper right')


def wfe_histogram_plot(
    opdtable,
    start_date=None,
    end_date=None,
    thresh=None,
    pid=None,
    download_opds=True,
    mark_corrections='lines',
    ote_only=False,
    min_wfe=60,
    max_wfe=None
):
    """Plot histogram and cumulative histogram of WFE over some time range.

    Parameters
    ----------
    opdtable : astropy.table.Table
        OPD table, retrieved from MAST. See functons in mast_wss.py
    start_date, end_date : astropy.time.Time objects, or strings
        Start and end dates for the plot time range. Default is March 2022 to present.
    thresh : int
        threshold to filter the RMS WFE
    download_opds : bool
        toggle downloading of OPDs from MAST
    pid : int, optional
        Program ID for which dates of observations are to be overplotted
    min_wfe, max_wfe : float or None
        if provided, minimum and maximum WFE values to limit plot display at, in nanometers RMS. These will set the
        Y axis range for the top time series plot, and the X axis range for the histogram plot. Note, any histogram
        values outside of this range are clipped such that they will appear as an 'extra pile up' of data at the
        extreme values. I.e. any fraction of time above max_wfe will show up in the histogram as if extra time right
        at the max_wfe value. This is intentional such that the fraction of time shown in the histogram always sums
        to 100%.

    Returns
    -------
    Nothing, but makes a plot

    """

    if start_date is None:
        start_date = astropy.time.Time('2022-07-16')
    if end_date is None:
        end_date = astropy.time.Time.now()
    start_date = astropy.time.Time(start_date)  # cast any input string to Date object
    end_date = astropy.time.Time(end_date)  # cast any input string to Date object

    # Look up wavefront sensing and mirror move corrections for that month
    if pid:
        pid_dates = get_dates_for_pid(pid)
        if pid_dates is not None:
            pid_dates = pid_dates[(pid_dates >= start_date) & (pid_dates <= end_date)]
        else:
            pid = None

    opdtable0 = webbpsf.mast_wss.deduplicate_opd_table(opdtable)
    opdtable1 = webbpsf.mast_wss.filter_opd_table(opdtable0, start_time=start_date, end_time=end_date)

    if download_opds:
        webbpsf.mast_wss.download_all_opds(opdtable1)

    # Retrieve all RMSes, from the FITS headers.
    # These are observatory WFE (OTE + NIRCam), at the WFS sensing field point
    rmses = []

    if 'rms_wfe' in opdtable1.colnames:
        rmses = opdtable1['rms_wfe']

    mjds = []
    pre_or_post = []
    for row in opdtable1:
        if download_opds:
            full_file_path = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'MAST_JWST_WSS_OPDs', row['fileName'])
        else:
            full_file_path = row['fileName']
        if 'rms_wfe' not in opdtable1.colnames:
            if ote_only is False:
                rmses.append(fits.getheader(full_file_path, ext=1)['RMS_WFE'])
            elif ote_only is True:
                opd_data = fits.getdata(full_file_path, ext=1)
                mask = opd_data != 0
                sensing_apername = fits.getheader(full_file_path, ext=0)['APERNAME']

                # Get WSS Target Phase Map for the sensing aperture
                # Note that the sensing maintenance program changed field point from NRC A3 to A1 around Dec 2024.
                was_targ_file = webbpsf.utils.get_target_phase_map_filename(sensing_apername)


                target_1024 = astropy.io.fits.getdata(was_targ_file)
                target_256 = poppy.utils.krebin(target_1024, (256, 256)) / 16
                wf_si = target_256 * mask  # Nircam target phase map at FP1

                rmses.append(webbpsf.utils.rms(opd_data - wf_si, mask=mask))

        mjds = opdtable1['date_obs_mjd']
        pre_or_post.append(webbpsf.mast_wss.infer_pre_or_post_correction(row))

    where_post = ['post' in a for a in pre_or_post]

    dates = astropy.time.Time(opdtable1['date'], format='isot')

    # Interpolate those RMSes into an even grid over time
    interp_fn = scipy.interpolate.interp1d(mjds, rmses, kind='nearest')

    mjdrange = np.linspace(np.min(mjds), np.max(mjds), 2048)
    interp_rmses = interp_fn(mjdrange)

    if max_wfe is not None:
        interp_rmses = np.clip(interp_rmses, min_wfe / 1e3, max_wfe / 1e3)
        # note the interp_wfe array is in microns, so have to rescale min and max_wfe values before clipping

    # Plot
    hspace = 0.3
    nrows = 2
    fig, axes = plt.subplots(figsize=(16, 10), nrows=nrows, gridspec_kw={'hspace': hspace})

    ms = 14  # markersize

    sensing_markers, = axes[0].plot_date(dates.plot_date, np.asarray(rmses) * 1e3, '.', ms=ms, ls='-',
                                         label='Sensing visit')
    axes[0].xaxis.set_major_locator(matplotlib.dates.DayLocator(bymonthday=[1]))
    axes[0].xaxis.set_minor_locator(matplotlib.dates.DayLocator(interval=1))
    axes[0].tick_params('x', length=10, rotation=30)

    if mark_corrections == 'lines':
        # Add vertical lines for corrections
        icorr = 0
        for i, idate in enumerate(where_post):
            if idate is True:
                plot = axes[0].axvline(dates[i].plot_date, ymin=0, ymax=500, linestyle='dashed', color='red')
                if icorr == 0:
                    plot.set_label('Corrections')
                    icorr += 1
        axes[0].legend()
    elif mark_corrections == 'triangles':
        yval = (np.asarray(rmses) * 1e3).max() * 1.01
        if max_wfe:
            yval = min(yval, max_wfe * 0.98)
        axes[0].scatter(
            dates[where_post].plot_date,
            np.ones(np.sum(where_post)) * yval,
            marker='v',
            s=100,
            color='limegreen',
            label='Corrections',
        )
        axes[0].legend()
    elif mark_corrections == 'arrows':
        rms_nm = np.asarray(rmses) * 1e3

        for i, idate in enumerate(where_post):
            if idate:
                axes[0].annotate('', xy=(dates[i].plot_date, rms_nm[i]),
                                 xytext=(dates[i - 1].plot_date, rms_nm[i - 1] if
                                         (max_wfe is None or rms_nm[i-1] < max_wfe) else max_wfe),
                                 color='limegreen', arrowprops=dict(arrowstyle='-|>', color='limegreen', linewidth=2,
                                                                    mutation_scale=20),
                                 xycoords='data'
                )
        arrow_placeholder = matplotlib.lines.Line2D([], [], color='limegreen', marker='v', ls='none',
                                                    markersize=10, label='Corrections')
        axes[0].legend(handles=[sensing_markers, arrow_placeholder])

    if pid:
        axes[0].set_ylim(0.975 * axes[0].get_ylim()[0], 1.025 * axes[0].get_ylim()[1])
    if max_wfe:
        axes[0].set_ylim(min_wfe, max_wfe)

    fig_title = 'OTE' if ote_only else 'Observatory'
    ylabel = 'OTE-only' if ote_only else 'OTE+NIRCam'

    axes[0].set_xlabel('Date')
    axes[0].set_ylabel(f'RMS WFE\n({ylabel})', fontweight='bold')

    axes[0].set_title(
        f'{fig_title} WFE from {start_date.isot[0:10]} to {end_date.isot[0:10]}', fontsize=14, fontweight='bold'
    )

    if thresh:
        axes[0].axhline(thresh, color='C2', label='OTE Correction threshold', linestyle='dashed')

    axes[0].tick_params(right=True, which='both', direction='in')

    nbins = 100
    binwidth = 1
    minbin = np.round(np.min(interp_rmses * 1e3) - binwidth)
    maxbin = np.round(np.max(interp_rmses * 1e3) + binwidth)

    axes[1].set_title(
        f'{fig_title} WFE Histogram from {start_date.isot[0:10]} to {end_date.isot[0:10]}', fontsize=14, fontweight='bold'
    )

    hist_values = axes[1].hist(
        interp_rmses * 1e3, density=True, bins=np.arange(minbin, maxbin, binwidth), color='#1f77b4', rwidth=0.95
    )
    axes[1].set_ylabel('Fraction of time with this WFE', fontweight='bold', color='#1f77b4')
    axes[1].set_xlabel('RMS Wavefront Error [nm]')
    axes[1].minorticks_on()

    ax_right = axes[1].twinx()
    ax_right.hist(interp_rmses * 1e3, density=True, bins=nbins, cumulative=1, histtype='step', lw=3, color='C1')
    ax_right.set_ylabel('Cumulative fraction of time\nwith this WFE or better', color='C1', fontweight='bold')
    ax_right.minorticks_on()

    xmin = min_wfe
    xmax = interp_rmses.max() * 1e3
    ymax = hist_values[0].max()
    axes[1].set_xticks(np.arange(xmin, xmax, 1), minor=True)

    ax_right.set_xlim(xmin, xmax)
    axes[1].set_xlim(xmin, xmax)
    ax_right.set_ylim(0, 1)

    if thresh:
        for i in [1]:
            if thresh <= xmax:
                axes[i].axvline(thresh, color='C2', linestyle='dashed')
            fractime = (interp_rmses * 1e3 < thresh).sum() / len(interp_rmses)
            axes[i].text(
                xmin + 0.68 * (xmax - xmin),
                0.65 * ymax,
                f'{fractime*100:.1f}% of the time has \nmeasured {ylabel} WFE < {thresh}',
                color='C2',
                fontweight='bold',
                fontsize=14,
            )

    # Add vertical lines for dates of PID observations
    if pid:
        for i, obs_date in enumerate(pid_dates):
            axes[0].axvline(obs_date.plot_date, color='darkgrey', ls='--', alpha=0.5, zorder=1)
            y_star = axes[0].get_ylim()[0] + 0.10 * np.diff(axes[0].get_ylim())
            if i == 0:
                label = 'PID {:d} obs. date(s)'.format(pid)
            else:
                label = None
            axes[0].scatter(obs_date.plot_date, y_star, marker='*', s=200, color='darkgrey', label=label)



# Wavefront Drifts Plot #####


def show_opd_image(
    array,
    vmax=0.300,
    ax=None,
    title=None,
    colorbar=False,
    labelrms=True,
    mask=None,
    deltatime_hrs=None,
    fontsize=12,
    maskc3=None,
):
    """Show and annotate an OPD

    Note, assumes the input OPD array is in units of microns, consistent with WSS outputs
    """
    from webbpsf.optical_budget import rms

    if ax is None:
        plt.figure(figsize=(10, 8))
        im = plt.imshow(array, cmap=matplotlib.cm.RdBu_r, vmin=-vmax, vmax=vmax, origin='lower')
        ax = im.axes
    else:
        im = ax.imshow(array, cmap=matplotlib.cm.RdBu_r, vmin=-vmax, vmax=vmax, origin='lower')

    ax.set_facecolor('gray')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    if labelrms:
        ax.text(
            0.02,
            0.02,
            f'{rms(array, mask) * 1000:.1f} nm rms',
            color='yellow',
            fontweight='bold',
            fontsize=fontsize,
            verticalalignment='bottom',
            transform=ax.transAxes,
        )
        if maskc3 is not None:
            ax.text(
                array.shape[0] - 5,
                5,
                f'ex. C3:\n{rms(array, maskc3)*1000:.1f} nm rms',
                color='lightgray',
                fontweight='bold',
                fontsize=fontsize,
                horizontalalignment='right',
            )

    if deltatime_hrs is not None:
        # ax.text(array.shape[0], 5,
        ax.text(
            0.995,
            0.02,
            f'{rms(array, mask) * 1000 / deltatime_hrs:.2f} nm/hr ',
            color='cyan',
            fontweight='bold',
            fontsize=fontsize,
            horizontalalignment='right',
            verticalalignment='bottom',
            transform=ax.transAxes,
        )
        if maskc3 is not None:
            ax.text(
                array.shape[0] - 5,
                -5,
                f'{rms(array, maskc3)/deltatime_hrs*1000:.2f} nm/hr ',
                color='lightgray',
                fontweight='bold',
                fontsize=fontsize,
                horizontalalignment='right',
                verticalalignment='top',
            )

    if colorbar:
        plt.colorbar(mappable=im, label='nanometers rms')

    if title is not None:
        ax.set_title(
            title,
            fontsize=16,
            fontweight='bold',
        )

    return ax


# Single measurement trending plot


@functools.lru_cache
def _get_mimf2_focus_offset_model(npix=256):
    # Precompute a model for the MIMF2 focus adjustment
    ote = webbpsf.opds.OTE_Linear_Model_WSS(npix=npix)
    ote.move_sm_local(piston=-1.789125, trans_unit='micron')  # Piston from M2022042201_sur.xml
    mimf2_focus_offset = ote.opd * 1e6  # convert meters to microns
    return mimf2_focus_offset


def single_measurement_trending_plot(
    opdtable, row_index=-1, reference=None, verbose=True, vmax=0.15, ignore_missing=True, subtract_target=True
):
    """Wavefront trending plot for a single measurement

    Parameters
    ----------

    opdtable : astropy.table.Table
        Table of available OPDs, as returned by retrieve_mast_opd_table()
    row_index : int
        Index into that table. Which row to make a plot for? Default: latest OPD.
    verbose: bool
        be more verbose in output?
    vmax : float
        Image display scale max for OPD, in microns. Defaults to 0.15 microns = 150 nanometers
    reference: str, or None
        Reference OPD to use for comparison, where the string is the date format (e.g. 2022-10-31).
        Will select closest match from the opdtable. Default: MIMF2 meas. + focus offset.


    """

    mimf2_focus_offset = _get_mimf2_focus_offset_model()

    filename = opdtable[row_index]['fileName']
    if verbose:
        print(f'Generating trending plot for {filename}')

    # Read current measurement
    opd, opdhdu = _read_opd(filename)
    mask = opd != 0
    opd[~mask] = np.nan
    sensing_apername = opdhdu[0].header['APERNAME']

    # hdr_rmswfe = opdhdu[1].header['RMS_WFE']
    visit = opdhdu[0].header['OBS_ID'][0:12]

    nanmask = np.ones_like(opd)
    nanmask[~mask] = np.nan

    mask_segs = webbpsf.utils.get_pupil_mask(opd.shape[0], label_segments=True)
    mask_without_C3 = (mask_segs != 0) & (mask_segs != 12)

    # Read prior WFS measurement
    prior_opd_index = row_index - 1
    if 'F187N' in opdtable[prior_opd_index]['wfs_measurement_type']:
        if verbose:
            print('     Immediately prior OPD file appears to be in F187N; skipping that and looking at next earlier OPD.')
            prior_opd_index = prior_opd_index - 1
    prev_filename = opdtable[prior_opd_index]['fileName']
    if verbose:
        print('     Previous OPD is:', prev_filename)

    prev_opd, prev_opd_hdu = _read_opd(prev_filename)
    prev_sensing_apername = prev_opd_hdu[0].header['APERNAME']

    if subtract_target:
        if verbose:
            print('     Subtracting NIRCam SI WFE target phase map')

        # Get WSS Target Phase Map for the sensing aperture
        # Note that the sensing maintenance program changed field point from NRC A3 to A1 around Dec 2024.
        was_targ_file = webbpsf.utils.get_target_phase_map_filename(sensing_apername)
        prev_was_targ_file = webbpsf.utils.get_target_phase_map_filename(prev_sensing_apername)


        target_1024 = astropy.io.fits.getdata(was_targ_file)
        target_256 = poppy.utils.krebin(target_1024, (256, 256)) / 16  # scale factor for rebinning w/out increasing values

        if prev_was_targ_file != was_targ_file:
            prev_target_1024 = astropy.io.fits.getdata(prev_was_targ_file)
            prev_target_256 = poppy.utils.krebin(prev_target_1024,
                                                 (256, 256)) / 16  # scale factor for rebinning w/out increasing values
        else:
            prev_target_256 = target_256


        opd -= target_256
        prev_opd -= prev_target_256

    # Compute deltas and decompose
    deltatime = get_datetime_utc(opdhdu, return_as='astropy') - get_datetime_utc(prev_opd_hdu, return_as='astropy')
    deltatime_hrs = deltatime.to_value(u.hour)
    delta_opd = opd - prev_opd
    fit, coeffs = webbpsf.opds.decompose_opd_segment_PTT(delta_opd)

    # Compare to reference nominal OPD
    # For this we use the ~ best correction, immediately prior to MIMF-2 measurement
    # And we add in a model for the MIMF2 focus offset, applied just after that.

    if reference is None:
        ref_row = opdtable[opdtable['visitId'] == 'V01163030001']
        _, ref_opd_hdu = _read_opd(ref_row['fileName'][0])
        ref_opd = ref_opd_hdu[1].data + mimf2_focus_offset
        ref_label = '(from MIMF2)'
    else:
        actual_ref = min(astropy.time.Time(opdtable['date']), key=lambda x: abs(x - astropy.time.Time(reference))).value
        ref_row = opdtable[opdtable['date'] == actual_ref]
        _, ref_opd_hdu = _read_opd(ref_row['fileName'][0])
        ref_opd = ref_opd_hdu[1].data
        ref_label = get_datetime_utc(ref_opd_hdu)

    if subtract_target:
        ref_opd -= target_256

    # Read associated post-correction measurement, if present
    if opdtable[row_index]['is_pre_correction']:
        show_post_move = True
        meas_title = 'Measurement (Pre Move)'
        match = opdtable[(opdtable['visitId'] == opdtable[row_index]['visitId']) & opdtable['is_post_correction']][0]
        post_opd, post_opd_hdu = _read_opd(match['fileName'])

        if subtract_target:
            post_opd -= target_256

        # Compare the post-move sensing to the reference OPD
        diff2_title = 'post-move - reference'
        delta_opd2 = post_opd - ref_opd
        delta_opd2[~mask] = np.nan

        # infer the correction ID
        # based on the assumption that the move probably came from the immediately prior WSS session
        prior_corr_id = opdtable[row_index - 1]['corr_id']
        if verbose:
            print('   That measurement has a correction, therefore reading in post-move sensing data and SUR.')
            print(f'   Inferred prior correction ID is {prior_corr_id}')
        sur_fn = f'{prior_corr_id}_sur.xml'
    else:
        show_post_move = False
        meas_title = 'Measurement'
        # Compare this latest sensing measurement to the reference OPD
        diff2_title = 'latest - reference'
        delta_opd2 = opd - ref_opd

    fit2, coeffs = webbpsf.opds.decompose_opd_segment_PTT(delta_opd2)

    # If we have info on the WAS' suggested corrections, we can plot that too.

    if 'EXPECTED' in opdhdu:
        show_correction = True
        correction = opdhdu['EXPECTED'].data - opdhdu['RESULT_PHASE'].data
        correction_mask = np.ones_like(correction, float)
        correction[correction == 0] = np.nan
    else:
        show_correction = False

    # Plotting ###########
    # Plot setup
    # fig, axes = plt.subplots(figsize=(8.5,11), nrows=3, ncols=4)
    fig, axes = plt.subplots(
        figsize=(16, 12),
        nrows=3,
        ncols=4,
        gridspec_kw={'top': 0.9, 'bottom': 0.01, 'hspace': 0.2, 'left': 0.05, 'right': 0.87},
    )

    title = f'{filename}    {visit}     {get_datetime_utc(opdhdu)}'
    if subtract_target:
        title += '\nNIRCam FP1 Target Phase Map Subtracted'
    plt.suptitle(title, fontweight='bold', fontsize=18)

    # Row 1: Latest measurement, and correction if present #####
    fontsize = 11
    # Panel 1: latest OPD
    iax = axes[0, 0]
    show_opd_image(opd, ax=iax, title=None, vmax=vmax, mask=mask, maskc3=mask_without_C3, fontsize=fontsize)
    iax.set_title(f'{meas_title}\n{get_datetime_utc(opdhdu)}', fontsize=fontsize * 1.2, fontweight='bold')

    if show_post_move:
        iax = axes[0, 1]
        show_opd_image(
            post_opd * nanmask, ax=iax, title=None, vmax=vmax, mask=mask, maskc3=mask_without_C3, fontsize=fontsize
        )
        iax.set_title(
            f'Measurement (Post Move)\n{get_datetime_utc(post_opd_hdu)}', fontsize=fontsize * 1.2, fontweight='bold'
        )

        sur_opd = webbpsf.opds.sur_to_opd(sur_fn, ignore_missing=ignore_missing)
        # cosmetic: handle masking slightly differently here to accomodate slightly different edge pixels
        surnanmask = nanmask.copy()
        surnanmask[sur_opd == 0] = np.nan

        iax = axes[0, 2]
        show_opd_image(sur_opd * surnanmask, ax=iax, title=None, vmax=vmax, mask=mask)
        iax.set_title(f'Mirror Move Commanded\nSUR {sur_fn}', fontsize=fontsize)

        if np.all(np.isnan(sur_opd)):
            iax.text(
                0.5,
                0.5,
                f'SUR file not found:\n{sur_fn}\n(Download of these is \nnot automated yet)',
                color='white',
                transform=iax.transAxes,
            )
        # iax.text(128, 128, f"SUR retrieval\nNot yet implemented", color='black', fontweight='bold',
        #        horizontalalignment='center')

        iax = axes[0, 3]
        show_opd_image((post_opd - opd) * nanmask, ax=iax, title=None, vmax=vmax, mask=mask, maskc3=mask_without_C3)
        iax.set_title('Mirror Move Measured\nDelta WFE', fontsize=fontsize)

    else:
        for ax in axes[0, 1:4]:
            ax.set_visible(False)
        fig.text(
            0.55, 0.77, 'Sensing-only visit. No mirror moves.', alpha=0.3, horizontalalignment='center', fontsize=fontsize
        )
    # Row 2 #####
    # Compare to immediate prior OPD

    # Panel 2-1: prior OPD
    iax = axes[1, 0]
    show_opd_image(prev_opd * nanmask, ax=iax, title='Prior WFS', vmax=vmax, maskc3=mask_without_C3, fontsize=fontsize)
    iax.set_title(f'Prior Measurement\n{get_datetime_utc(prev_opd_hdu)}', fontsize=fontsize * 1.2, fontweight='bold')

    # Panel 2-2: difference
    iax = axes[1, 1]
    show_opd_image(delta_opd, ax=iax, vmax=vmax, deltatime_hrs=deltatime_hrs, fontsize=fontsize)
    iax.set_title(f'Difference\n{deltatime.to(astropy.units.hour):.2f}', fontsize=fontsize * 1.1)

    # Panel 2-3: proposed correction
    iax = axes[1, 2]
    show_opd_image(fit, ax=iax, vmax=vmax, deltatime_hrs=deltatime_hrs, fontsize=fontsize)
    iax.set_title('Controllable Modes\nin difference', fontsize=fontsize * 1.1)

    # Panel 2-4:
    iax = axes[1, 3]
    show_opd_image(delta_opd - fit, ax=iax, vmax=vmax, fontsize=fontsize)
    iax.set_title('High order WFE\nin difference', fontsize=fontsize * 1.1)

    # Row 3 #####

    # Panel 3-1: ref OPD
    iax = axes[2, 0]
    show_opd_image(ref_opd * nanmask, ax=iax, title='Prior WFS', vmax=vmax, mask=mask, fontsize=fontsize)
    iax.set_title(f'Reference Measurement\n{ref_label}', fontsize=fontsize * 1.2, fontweight='bold')

    # Panel 3-2: difference
    iax = axes[2, 1]
    show_opd_image(delta_opd2 * nanmask, ax=iax, vmax=vmax, mask=mask, fontsize=fontsize)
    iax.set_title(f'Difference\n{diff2_title}', fontsize=fontsize * 1.1)

    # Panel 3-3: proposed correction
    iax = axes[2, 2]

    if show_correction:
        # Note, the WSS does not remove piston from the correction so it's not zero mean.
        # Fix that here before display. This also affects the derived rms value displayed.
        correction -= np.nanmean(correction[correction_mask == 1])

        show_opd_image(-correction, ax=iax, vmax=vmax, mask=correction_mask, fontsize=fontsize)
        iax.set_title('Controllable modes\nfrom WSS proposed correction', fontsize=fontsize * 1.1)

    #        show_opd_image(fit2, ax=axes[2,3], vmax=vmax, mask=mask, fontsize=fontsize)

    else:
        show_opd_image(fit2, ax=iax, vmax=vmax, mask=mask, fontsize=fontsize)
        iax.set_title('Controllable Modes\nin difference', fontsize=fontsize * 1.1)

    # Panel 3-4:
    iax = axes[2, 3]
    show_opd_image(delta_opd2 - fit2, ax=iax, vmax=vmax, maskc3=mask_without_C3, fontsize=fontsize)
    iax.set_title('High order WFE\nin difference', fontsize=fontsize * 1.1)

    cax = fig.add_axes([0.91, 0.41, 0.01, 0.4])

    cmap = matplotlib.cm.RdBu_r
    mappable_nm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(-vmax * 1000, vmax * 1000))
    fig.colorbar(mappable_nm, cax=cax, label='WFE [nm]')

    # Draw a bar plot showing potential correction and green/yellow/red thresholds
    cax = fig.add_axes([0.90, 0.03, 0.05, 0.2])
    norm = matplotlib.colors.Normalize(0, 0.15)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('thresholds', ['LightGreen', 'Yellow', 'Salmon'], N=3)
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax)
    cb.set_ticks([0.05, 0.10])
    cax.axvline(0.33, color='black', alpha=0.5)
    cax.axvline(0.66, color='black', alpha=0.5)
    cax.set_title(' Piston\nXtilt  \nYtilt  ', rotation=90, fontsize=14)
    ptts = np.abs(coeffs.reshape(18, 3))
    cax.plot(np.arange(18) / (64) + 0.02, ptts[:, 0], marker='+', ls='none', color='black')
    cax.plot(np.arange(18) / (64) + 0.34, ptts[:, 1], marker='+', ls='none', color='navy')
    cax.plot(np.arange(18) / (64) + 0.67, ptts[:, 2], marker='+', ls='none', color='black')
    cax.set_ylim(0, 0.15)
    fig.text(0.89, 0.30, 'Potential\nCorrections:', fontsize=13, fontweight='bold', horizontalalignment='left')


def series_of_measurement_trending_plots(opdtable, ignore_missing=False, start_date=None):
    """Generate the wavefront trending plot for all NRCA3 FP1 OPDs since the completion of OTE alignment

    Calls wavefront_trending_plot for all rows in the table, except any duplicates are ignored.

    Note, for this to work correctly, opdtable must be a ~complete table including back to the MIMF measurement.

    """

    opdtable = webbpsf.mast_wss.deduplicate_opd_table(opdtable)

    if start_date:
        start_time = astropy.time.Time(start_date)

    # Now that we have a clean list, iterate over all:
    year_yday = astropy.time.Time.now().yday[0:8]
    filename = f'wavefront_trending_{year_yday}.pdf'

    with PdfPages(filename) as pdf:
        for row_index in reversed(range(len(opdtable))):
            if start_date:
                if astropy.time.Time(opdtable[row_index]['date']) < start_time:
                    continue

            if opdtable[row_index]['is_post_correction']:
                # If this row is a post correction measurment, it doesn't get its own plot.
                # it gets folded in with the pre correction measurement plot
                continue

            single_measurement_trending_plot(opdtable, row_index=row_index, ignore_missing=ignore_missing)

            pdf.savefig()  # saves the current figure into a pdf page

        wavefront_time_series_plot(opdtable)
        pdf.savefig()
    return filename


# Some reference data : label which corrections in early cycle 1 had sigificant wing tilts
cids_with_left_wing_tilts = ['R20220523', 'R20220713', 'R20220715', 'O20220802', 'R20220802', 'R20221019']
cids_with_right_wing_tilts = ['R20220606', 'R20220627']


def cleanup_global_tt_for_wing_tilt(delta_opd, wing='left'):
    """Special case: cleanup global TT for the 16 segments not affected by a particularly large wing tilt event
    This helps remove systematics from the WFS results, which aren't actual motions of the mirrors.

    """
    npix = delta_opd.shape[0]

    segmask = webbpsf.utils.get_pupil_mask(npix=npix, label_segments=True)
    # Fit global PTT to everything that isn't on the affected wing
    if wing == 'left':  # -V2 wing
        bigtiltmask = (segmask != 0) & (segmask != 17) & (segmask != 16) & (segmask != 15)
    elif wing == 'right':  # +V2 wing
        bigtiltmask = (segmask != 0) & (segmask != 9) & (segmask != 10) & (segmask != 11)
    else:
        raise ValueError("Invalid value for wing parameter. Must be 'left' or 'right'. ")
    ptt_coeffs = poppy.zernike.decompose_opd_nonorthonormal_basis(
        delta_opd, aperture=(delta_opd != 0) & bigtiltmask, nterms=3
    )

    ptt_fit = poppy.zernike.compose_opd_from_basis(
        ptt_coeffs, poppy.zernike.zernike_basis, npix=npix, aperture=webbpsf.utils.get_pupil_mask(npix=npix)
    )
    # Subtract that from the delta OPD
    return delta_opd - ptt_fit


def wavefront_drift_plots(
    opdtable, start_time, end_time, verbose=False, vmax=100, n_per_row=9, label_cid=False, label_visit=False
):
    """Generate plots of wavefront drifts over time.

    This plots the natural drift in wavefront of the observatory
    (i.e. what it did on its own, not showing corrections).

    Parameters
    -----------
    start_time, end_time : Astropy Time instances
        Start and end times to define the range of OPDs to plot for
    verbose : bool
        Print more verbose output
    vmax : int
         Maximum value in nanometers for the color stretch of OPDs.
    n_per_row : int
        Number of plots per row
    label_cid : bool
        Add label for WSS correction ID to each plot
    label_visit : bool
        Add label for APT/OSS program and observation number to each plot
    """

    def vprint(*text):
        if verbose:
            print(*text)

    date_range_mask = (opdtable['date_obs_mjd'].value > start_time.mjd) & (opdtable['date_obs_mjd'].value < end_time.mjd)

    which_opds_mask = date_range_mask  # & (~ redundant_aps_mask) & ap_type_mask

    sum(which_opds_mask)

    # Iterate over all relevant OPDs to retrieve the OPD, and calc the WFE RMS
    # Setup arrays to store results from the iteration:
    n = np.sum(which_opds_mask)
    rms_obs = np.zeros(n)
    rms_ote = np.zeros(n)
    wf_obs = np.zeros((n, 256, 256))
    wf_ote = np.zeros((n, 256, 256))

    opd, opdhdu = _read_opd(opdtable[which_opds_mask][0]['fileName'])
    mask = opd != 0
    sensing_apername = opdhdu[0].header['APERNAME']

    # Get WSS Target Phase Map for the sensing aperture
    # Note that the sensing maintenance program changed field point from NRC A3 to A1 around Dec 2024.

    was_targ_file = webbpsf.utils.get_target_phase_map_filename(sensing_apername)

    target_1024 = astropy.io.fits.getdata(was_targ_file)
    target_256 = poppy.utils.krebin(target_1024, (256, 256))

    wf_si = target_256 * mask

    dates = []

    # Iterate over all selected OPDs
    for i, row in enumerate(opdtable[which_opds_mask]):
        vprint(f"Loading OPD from {row['fileName']}")

        opd, opdhdu = _read_opd(row['fileName'])

        wf_obs[i] = opd * 1e-6  # convert from microns to meters
        wf_ote[i] = (opd - wf_si) * 1e-6

        rms_obs[i] = webbpsf.utils.rms(wf_obs[i], mask)
        rms_ote[i] = webbpsf.utils.rms(wf_ote[i], mask)

        dates.append(row['date'])

    # Compute all the delta WFEs, which can come either from actual mirror drifts, or corrections
    wf_deltas = wf_ote[1:] - wf_ote[0:-1]

    labelfontsize = 7

    # Set up plot axes etc
    n_to_plot = sum(opdtable[which_opds_mask]['wfs_measurement_type'] == 'pre') - 1
    nrows = int(np.ceil(n_to_plot / n_per_row))

    fig, axes = plt.subplots(
        figsize=(16, nrows * 2 + 1),
        nrows=nrows,
        ncols=n_per_row,
        gridspec_kw={'hspace': 0.3, 'wspace': 0.01, 'left': 0.01, 'right': 0.93, 'bottom': 0.01, 'top': 0.92},
    )
    axes_f = axes.flat

    nanmask = np.zeros_like(wf_ote[0]) + np.nan
    nanmask[mask] = 1
    last_date_obs = np.nan

    deltas_shown = []

    # Iterate to make plots
    i = 0  # axes counter
    for ir, row in enumerate(opdtable[which_opds_mask]):
        if ir == 0:
            # Skip first row, for which we can't compute a delta to the prior
            last_date_obs = row['date_obs_mjd']
            continue
        if row['wfs_measurement_type'] == 'post':
            # don't plot corrections, but do record the time for use in the subsequent delta
            last_date_obs = row['date_obs_mjd']
            continue

        vprint(row['fileName'], row['date'], row['wfs_measurement_type'])
        date = row['date']

        delta_opd = wf_deltas[ir - 1]
        cid = row['fileName'][0:11]
        deltat = row['date_obs_mjd'] - last_date_obs

        if cid[0:9] in cids_with_left_wing_tilts:
            vprint('  Adjusting delta OPD to exclude left wing (-V2) from global PTT')
            delta_opd = cleanup_global_tt_for_wing_tilt(delta_opd)
        elif cid[0:9] in cids_with_right_wing_tilts:
            vprint('  Adjusting delta OPD to exclude right wing (+V2) from global PTT')
            delta_opd = cleanup_global_tt_for_wing_tilt(delta_opd, wing='right')

        deltas_shown.append(delta_opd)

        title = f'{date[0:10]}  {date[11:16]}\n$\\Delta T =${deltat * 24:.1f} hr'
        if label_cid:
            title = f'{cid}\n' + title
        if label_visit:
            visit = row['visitId']
            title = f'{int(visit[1:6])}:{int(visit[6:9])}\n' + title

        show_opd_image(
            delta_opd * 1e6 * nanmask,
            mask=mask,
            ax=axes_f[i],
            # deltatime_hrs=deltat*24,
            vmax=vmax / 1e3,
            fontsize=labelfontsize,
        )
        axes_f[i].set_title(title, fontsize=9)

        last_date_obs = row['date_obs_mjd']
        i += 1  # increment plot counter for next plot

    # Finishing details: hide unused axes, add colorbar and supertitle
    for j in range(i, len(axes_f)):
        axes_f[j].set_visible(False)
        pass

    cax = fig.add_axes([0.94, 0.26, 0.01, 0.4])
    mappable_nm = matplotlib.cm.ScalarMappable(cmap=matplotlib.cm.RdBu_r, norm=matplotlib.colors.Normalize(-vmax, vmax))
    fig.colorbar(mappable_nm, cax=cax, label='WFE [nm]')

    plt.suptitle(
        f'Telescope Wavefront Drifts: {start_time.isot[0:10]} to {end_time.isot[0:10]}\n', fontsize=20, fontweight='bold'
    )

    plt.savefig('wf_drifts.pdf')


# Monthly Trending Plots, including OPDs, RMS WFE and PSF EE


def get_month_start_end(year, month):
    _, ndays = calendar.monthrange(year, month)
    start_date_str = f'{year:04d}-{month:02d}-01'
    end_date_str = f'{year:04d}-{month:02d}-{ndays:02d}'

    start_date = astropy.time.Time(start_date_str)
    end_date = astropy.time.Time(end_date_str)

    return start_date, end_date


def filter_opdtable_for_daterange(start_date, end_date, opdtable):
    """Filter existing opdtable for a given time range
    This includes the last measurement in the prior time range too (if applicable), so we can compute a delta
    to the first one
    """
    # Start a little early, such that we are going to have at least 1 WFS before the start date
    pre_start_date = astropy.time.Time(start_date) - astropy.time.TimeDelta(4 * u.day)
    opdtable = webbpsf.mast_wss.filter_opd_table(opdtable, start_time=pre_start_date, end_time=end_date)
    if len(opdtable) == 0:
        raise ValueError('The opdtable is empty for this date range.')

    # Trim the table to have 1 and only 1 precursor measurement -
    # we'll use this to compute the drift for the first WFS in the time period
    is_pre = [astropy.time.Time(row['date']) < start_date for row in opdtable]
    opdtable['is_pre'] = is_pre
    opdtable = opdtable[np.sum(is_pre) - 1:]

    return opdtable


def filter_opdtable_for_month(year, mon, opdtable, shift_dates_offset=None):
    """Filter existing opdtable for a given month
    This includes the last measurement in the prior month too (if applicable), so we can compute a delta
    to the first one
    """
    start_date, end_date = get_month_start_end(year, mon)
    if shift_dates_offset:
        start_date += shift_dates_offset * u.day
        end_date += shift_dates_offset * u.day

    return filter_opdtable_for_daterange(start_date, end_date, opdtable)


def get_opdtable_for_daterange(start_date, end_date):
    """Return table of OPD measurements for date range.

    This includes the last measurement preceding this date range, too, so we
    can compute the first delta at the start of this range.
    """
    # Retrieve full OPD table, then trim to the selected time period
    opdtable0 = webbpsf.mast_wss.retrieve_mast_opd_table()
    opdtable0 = webbpsf.mast_wss.deduplicate_opd_table(opdtable0)

    opdtable = filter_opdtable_for_daterange(start_date, end_date, opdtable0)
    return opdtable


def get_opdtable_for_month(year, mon, **kwargs):
    """Return table of OPD measurements for a given month.
    retrieve the opd table and filter according to month/year designation
    """
    # Retrieve full OPD table, then trim to the selected time period
    opdtable0 = webbpsf.mast_wss.retrieve_mast_opd_table()
    opdtable0 = webbpsf.mast_wss.deduplicate_opd_table(opdtable0)

    opdtable = filter_opdtable_for_month(year, mon, opdtable0, **kwargs)
    return opdtable


def check_colnames(opdtable):
    """Check for expected column names in the OPD Table
    Parameters
    -----------
    opdtable: astropy.table.Table
        Table of available OPDs
    """
    colnames = webbpsf.mast_wss.get_colnames(True)
    if all(colname in colnames for colname in opdtable.colnames) is False:
        colnames = webbpsf.mast_wss.get_colnames()
        if all(colname in opdtable.colnames for colname in colnames) is False:
            raise KeyError('Expected the following column names in the opdtable: ' + str(colnames))
        else:
            opdtable = webbpsf.mast_wss.add_columns_to_track_corrections(opdtable)
    return opdtable


def get_dates_for_pid(pid, project='jwst'):
    """Check the archive for the start dates of each observation of the specified PID

    Parameters
    -----------
    pid: int
         Program ID for which the observation start dates are to be overplotted

    Returns:
    start_times: list of astropy.time.Time dates

    """

    try:
        obs_table = Observations.query_criteria(project='jwst', proposal_id=pid)

        cond_calib = obs_table['calib_level'] > 1
        obs_table = obs_table[cond_calib]

        obs_table['obs_num'] = [x['obs_id'][7:10] if x['calib_level'] == 2 else x['obs_id'][8:13] for x in obs_table]
        obs_table['obs_num'] = [x.strip('_') for x in obs_table['obs_num']]

        obs_by_num = obs_table.group_by('obs_num')

        start_times = []
        for key, group in zip(obs_by_num.groups.keys, obs_by_num.groups):
            start_times.append(group['t_min'].min())

        start_times = astropy.time.Time(start_times, format='mjd')
        start_times = start_times.to_value('fits')  # , subfmt='date')
        start_times = astropy.time.Time(np.unique(start_times))

        return start_times
    except:  # noqa - adequate errors can be raised all needing this error
        print('No access to data for PID {:d}'.format(pid))
        return


def monthly_trending_plot(year, month, verbose=True, instrument='NIRCam', filter='F200W', vmax=200, pid=None, opdtable=None,
                          shift_dates_offset=None):
    """Make monthly trending plot showing OPDs, mirror moves, RMS WFE, and the resulting PSF EEs

    year, month : integers
        Self explanatory
    pid : int, optional
        Program ID for which dates of observations are to be overplotted
    verbose : bool
        Print more verbose text output
    vmax : float
        Image display vmax for OPDs, given here in units of nanometers.
    opdtable : astropy.table.Table
        Table of available OPDs, Default None: as returned by retrieve_mast_opd_table()
    shift_dates_offset: int
        number of dates to shift later the time period, for showing trends around the month divisions
    """

    def vprint(*text):
        if verbose:
            print(*text)

    start_date, end_date = get_month_start_end(year, month)

    if shift_dates_offset:
        start_date += shift_dates_offset * u.day
        end_date += shift_dates_offset * u.day

    # Look up wavefront sensing and mirror move corrections for that month
    if pid:
        pid_dates = get_dates_for_pid(pid)
        if pid_dates is not None:
            pid_dates = pid_dates[(pid_dates >= start_date) & (pid_dates <= end_date)]
        else:
            pid = None

    try:
        if opdtable is None:
            vprint('obtaining opdtable for month')
            opdtable = get_opdtable_for_month(year, month, shift_dates_offset=shift_dates_offset)
        else:
            vprint('filtering opdtable for month')
            opdtable = check_colnames(opdtable)
            opdtable = filter_opdtable_for_month(year, month, opdtable, shift_dates_offset=shift_dates_offset)
    except ValueError as e:
        print(e)
        raise

    corrections_table = webbpsf.mast_wss.get_corrections(opdtable)

    inst = webbpsf.instrument(instrument)
    inst.filter = filter

    apmask = webbpsf.utils.get_pupil_mask()

    vprint(f'Computing PSFs for {len(opdtable)} OPD measurements in {year}-{month:02d}')
    vprint(f'Calculating for {instrument}, {filter}')

    psfs = []
    dates = []
    mjds = []
    rms_obs = []
    rms_ote = []
    wfes_ote = []

    # for row_index in reversed(range(mimf2_corection_index, len(opdtable))):
    for row_index in range(len(opdtable)):
        opd_fn = opdtable[row_index]['fileName']
        vprint(opdtable[row_index]['date'], opdtable[row_index]['fileName'])

        inst.load_wss_opd(opd_fn, verbose=verbose)

        vprint(f'Calculating PSF for {opd_fn}')
        psf = inst.calc_psf(fov_pixels=101, nlambda=1, add_distortion=False)

        psfs.append(psf)
        dates.append(opdtable[row_index]['date'])
        mjds.append(opdtable[row_index]['date_obs_mjd'])

        # Reconstruct some particular wavefronts so we can compute their RMSes
        # Here note that load_wss_opd stores some data in FITS extensions, precisely to enable this sort of thing
        # The primary WFE output from that function is the OTE WFE at the master chief ray (between NIRCam modules).
        wfe_ote = inst.pupilopd['PRIMARY'].data + inst.pupilopd['SENSING_OTE_FD_WFE'].data

        cid = opdtable[row_index]['fileName'][0:9]
        if cid in cids_with_left_wing_tilts:
            vprint('  Adjusting delta OPD to exclude left wing (-V2) from global PTT, for wing tilt')
            wfe_ote = cleanup_global_tt_for_wing_tilt(wfe_ote)
        elif cid in cids_with_right_wing_tilts:
            vprint('  Adjusting delta OPD to exclude right wing (+V2) from global PTT, for wing tilt')
            wfe_ote = cleanup_global_tt_for_wing_tilt(wfe_ote, wing='right')

        wfe_obs = wfe_ote + inst.pupilopd['SENSING_SI_WFE'].data

        wfes_ote.append(wfe_ote)

        rms_ote.append(webbpsf.utils.rms(wfe_ote, mask=apmask))
        rms_obs.append(webbpsf.utils.rms(wfe_obs, mask=apmask))

    # Compute some derived quantities based on the PSFs
    dates_array = astropy.time.Time(dates, format='isot')

    ees = []
    for psf in psfs:
        ees.append(webbpsf.measure_ee(psf))

    rms_obs = np.asarray(rms_obs)
    rms_ote = np.asarray(rms_ote)

    sensing_tds = astropy.time.Time(corrections_table['Post Move Sensing Time'], format='isot') - astropy.time.Time(
        corrections_table['Pre Move Sensing Time'], format='isot'
    )
    correction_times = astropy.time.Time(corrections_table['Pre Move Sensing Time'], format='isot') + sensing_tds / 2

    npoints = sum(opdtable[1:]['wfs_measurement_type'] == 'pre')

    # Prepare for plotting

    cmap = matplotlib.cm.RdBu_r.copy()
    cmap.set_bad('0.4')

    def basic_show_image(image, ax, vmax=0.3, nanmask=1):
        ax.imshow(image * nanmask, cmap=cmap, vmin=-vmax, vmax=vmax, origin='lower')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_yticks([])

    # Make the plots
    fig = plt.figure(
        constrained_layout=False,
        figsize=(16, 10),
    )

    subfigs = fig.subfigures(
        2,
        1,
        hspace=0.02,
        height_ratios=[2, 2],
    )
    for sf in subfigs:
        sf.set_facecolor('none')  # work around display issue that was hiding the connection patches

    min_n_im_axes = 12
    im_axes = subfigs[0].subplots(
        3,
        max(npoints, min_n_im_axes),
        gridspec_kw={'wspace': 0.001, 'left': 0.06, 'right': 0.98, 'bottom': 0.02, 'top': 0.8},
    )  # , gridspec_kw={'wspace':1.000, 'left': 0})

    axes = subfigs[1].subplots(2, 1, gridspec_kw={'left': 0.07, 'right': 0.98, 'top': 0.97})

    fs = 14  # Font size for axes labels

    title_extra = f", plus {shift_dates_offset} days" if shift_dates_offset else ""
    fig.suptitle(f'WF Trending for {year}-{month:02d}{title_extra}', fontsize=fs * 1.5, fontweight='bold')

    # Plot 1: Wavefront Error

    axes[0].plot_date(dates_array.plot_date, rms_obs * 1e9, color='C1', ls='-', label='Observatory WFE at NIRCam NRCA3')
    axes[0].plot_date(dates_array.plot_date, rms_ote * 1e9, color='C0', ls='-', label='Telescope WFE')

    for ax in axes:
        for corr_date in correction_times:
            ax.axvline(corr_date.plot_date, color='darkgreen', ls='--', alpha=0.5)

    # axes[0].axhline(59, ls=":", color='gray')
    # axes[0].axhline(80, ls=":", color='gray')
    axes[0].fill_between(
        [start_date.plot_date - 0.5, end_date.plot_date + 0.5],
        [59, 59],
        [80, 80],
        color='blue',
        alpha=0.08,
        label='Wavefront control target range',
    )

    axes[0].set_ylim(0.8 * rms_ote.min() * 1e9, 1.2 * rms_obs.max() * 1e9)
    axes[0].set_ylabel('Wavefront Error\n[nm rms]', fontsize=fs, fontweight='bold')
    axes[0].set_xticklabels([])

    # Plot 2: Encircled Energies

    ee_ax_ylim = 0.04
    ee_measurements = {}
    for i, ee_npix in enumerate([10, 2.5]):
        color = f'C{i*2+2}'
        ee_rad = inst.pixelscale * ee_npix
        ees_at_rad = np.asarray([e(ee_rad) for e in ees])
        ee_measurements[ee_npix] = ees_at_rad  # save for later

        median_ee = np.median(ees_at_rad)
        ee_ax_ylim = np.max(
            [ee_ax_ylim, np.abs((ees_at_rad - median_ee) / median_ee).max() * 1.1]
        )  # display tweak: adjust the plot Y scale sensibly to its contents

        axes[1].plot_date(
            dates_array.plot_date,
            (ees_at_rad - median_ee) / median_ee,
            ls='-',
            color=color,
            label=f'$\\Delta$EE within {ee_rad:.2f} arcsec ({ee_npix} pix)',
        )

        axes[1].text(
            0.01,
            0.75 - i * 0.12,
            f'Median EE within {ee_rad:.2f} arcsec = {median_ee:.3f}',
            color=color,
            fontweight='bold',
            transform=axes[1].transAxes,
        )

    axes[1].fill_between(
        [start_date.plot_date - 0.5, end_date.plot_date + 0.5],
        -0.03,
        0.03,
        color='gray',
        alpha=0.1,
        label='±3% change (stability requirement)',
    )
    axes[1].set_xlabel('Date', fontsize=fs, fontweight='bold')
    axes[1].set_ylabel(f'% Change in \nEncircled Energy\n{instrument} {filter}', fontsize=fs, fontweight='bold')
    axes[1].set_ylim(0.5, 1.0)
    axes[1].axhline(0, ls=':', color='gray')
    axes[1].set_ylim(-ee_ax_ylim, ee_ax_ylim)

    # Configure Axes for the time series plots
    for ax in axes[0:2]:
        ax.set_xlim(
            start_date.plot_date - 0.5,
            end_date.plot_date + 0.5,
        )
        ax.xaxis.set_minor_locator(matplotlib.dates.DayLocator())
        ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(byweekday=matplotlib.dates.MONDAY, interval=1))

    # Add vertical lines for dates of PID observations
    if pid:
        for iax, ax in enumerate(axes):
            for i, obs_date in enumerate(pid_dates):
                ax.axvline(obs_date.plot_date, color='darkgrey', ls='--', alpha=0.5, zorder=1)
                y_star = ax.get_ylim()[0] + 0.20 * np.diff(ax.get_ylim())
                if iax == 0 and i == 0:
                    label = 'PID {:d} obs. date(s)'.format(pid)
                else:
                    label = None
                ax.scatter(obs_date.plot_date, y_star, marker='*', s=200, color='darkgrey', label=label)

    axes[0].legend(loc='best', fontsize=9, ncol=2)
    axes[1].legend(loc='upper right', fontsize=9)

    nanmask = np.zeros_like(apmask) + np.nan
    nanmask[apmask == 1] = 1
    plot_index = -1
    from matplotlib.patches import ConnectionPatch

    ax0ymax = axes[0].get_ylim()[1]  # for use in connection patches below
    for i, row in enumerate(opdtable):
        if verbose:
            print(
                row['fileName'],
                row['date'],
                'Sensing' if row['wfs_measurement_type'] == 'pre' else 'Post Mirror Move',
                f'Obs WFE: {rms_obs[i] * 1e9:.1f} nm\tOTE WFE: {rms_ote[i]*1e9:.1f} nm',
            )
        date = row['date']
        if astropy.time.Time(date) < start_date:
            continue

        delta_opd = wfes_ote[i] - wfes_ote[i - 1] if i > 0 else np.zeros((1024, 1024))
        vmax_micron = vmax / 1000
        rms_label = None
        if row['wfs_measurement_type'] == 'pre':
            plot_index += 1
            # Plot WFS in row 1
            basic_show_image(
                wfes_ote[i] * 1e6, ax=im_axes[0, plot_index], nanmask=nanmask, vmax=vmax_micron
            )  # , title=None)
            rms_label = im_axes[0, plot_index].text(
                20, 20, f'{webbpsf.utils.rms(wfes_ote[i], mask=apmask)*1e9:.1f}', color='yellow', fontsize=fs * 0.6
            )
            # Plot drift since last measurement in row 2
            basic_show_image(delta_opd * 1e6, ax=im_axes[1, plot_index], nanmask=nanmask, vmax=vmax_micron)
            im_axes[1, plot_index].text(
                20, 20, f'{webbpsf.utils.rms(delta_opd, mask=apmask)*1e9:.1f}', color='yellow', fontsize=fs * 0.6
            )
            # mark row 3 as blank (default)
            basic_show_image(delta_opd + np.nan, ax=im_axes[2, plot_index], nanmask=nanmask, vmax=vmax_micron)
            im_axes[0, plot_index].set_title(row['date'][0:10])

            # Make a fancy connector line from the OPD plots to the correct time in the time series plots
            # Draw this for all sensing instances pre-move, but no need to repeat for any post-correction sensing
            cp = ConnectionPatch(
                [0.5, 0],
                (dates_array[i].plot_date, ax0ymax),
                coordsA='axes fraction',
                coordsB='data',
                axesA=im_axes[2, plot_index],
                axesB=axes[0],
                color='darkgreen',
                ls='--',
                alpha=0.5,
            )
            fig.add_artist(cp)

        else:
            # Update row 1 to show post-mirror-move WFS
            basic_show_image(
                wfes_ote[i] * 1e6, ax=im_axes[0, plot_index], nanmask=nanmask, vmax=vmax_micron
            )  # , title=None)
            if rms_label is not None:
                rms_label.set_visible(False)
                del rms_label  # delete the previously-written one for the pre-move sensing

            im_axes[0, plot_index].text(
                20, 20, f'{webbpsf.utils.rms(wfes_ote[i], mask=apmask)*1e9:.1f}', color='yellow', fontsize=fs * 0.6
            )
            # Plot correction in row 3
            basic_show_image(delta_opd * 1e6, ax=im_axes[2, plot_index], nanmask=nanmask, vmax=vmax_micron)
            im_axes[2, plot_index].text(
                20, 20, f'{webbpsf.utils.rms(delta_opd, mask=apmask)*1e9:.1f}', color='yellow', fontsize=fs * 0.6
            )

    for i, label in enumerate(['Measured\nWFE', 'Drifts', 'Mirror\nCorrections']):
        im_axes[i, 0].yaxis.set_visible(True)
        im_axes[i, 0].set_ylabel(label + '\n\n', fontsize=fs, fontweight='bold')

    for j in range(npoints, min_n_im_axes):
        for i in range(3):
            im_axes[i, j].set_visible(False)

    if shift_dates_offset:
        outname = f'wf_trending_{year}-{month:02d}+{shift_dates_offset}days.pdf'
    else:
        outname = f'wf_trending_{year}-{month:02d}.pdf'
    plt.savefig(outname, dpi=200, bbox_inches='tight')
    vprint(f'Saved to {outname}')

    wfs_type = [('Sensing' if row['wfs_measurement_type'] == 'pre' else 'Post Mirror Move') for row in opdtable]
    result_table = astropy.table.QTable(
        [
            opdtable['date'],
            opdtable['fileName'],
            wfs_type,
            rms_obs * 1e9 * u.nm,
            rms_ote * 1e9 * u.nm,
            ee_measurements[2.5],
            ee_measurements[10],
        ],
        names=['Date', 'Filename', 'WFS Type', 'RMS WFE (OTE+SI)', 'RMS WFE (OTE only)', 'EE(2.5 pix)', 'EE(10pix)'],
    )
    return result_table


def plot_phase_retrieval_crosscheck(fn, vmax_fraction=1.0):
    """Make a plot which displays the quality of the phase retrievals,
    comparing the WSS-generated retrieved PSF with the actual measured PSF,
    for both positive and negative defocus.

    This can be used to investigate particular phase retrieval results,
    among other things offering a quick check against the occasional
    previously-undetected binary star which may crop up in the WFS target pool.

    Parameters
    ----------
    fn : string
        Filename of an OPD file retrieved from MAST.
    vmax_fraction : float
        Scale factor for setting the plot log scale vmax relative to the
        image peak pixel. vmax_fraction = 1.0 sets the scale vmax to the
        image max value, vmax_fraction = 0.1 sets it to 1/10th that, etc.
    """
    from skimage.registration import phase_cross_correlation

    _ = webbpsf.mast_wss.mast_retrieve_opd(fn)

    opd, hdul = webbpsf.trending._read_opd(fn)

    # measured
    wlm8_m = hdul[5].data
    wlp8_m = hdul[10].data
    # computed
    wlm8_c = poppy.utils.pad_or_crop_to_shape(hdul[6].data, wlm8_m.shape)
    wlp8_c = poppy.utils.pad_or_crop_to_shape(hdul[11].data, wlm8_m.shape)

    cm = matplotlib.cm.inferno.copy()
    cm.set_bad(cm(0))

    fscale = wlm8_m.sum() / wlm8_c.sum()

    norm = matplotlib.colors.LogNorm(wlm8_m.max() / 1e4 * vmax_fraction, wlm8_m.max() * vmax_fraction)

    # Shift to register
    shift, error, diffphase = phase_cross_correlation(wlm8_m, wlm8_c * fscale)
    wlm8_cs = np.roll(wlm8_c, tuple(np.array(shift, int)), axis=(0, 1))

    shift, error, diffphase = phase_cross_correlation(wlp8_m, wlp8_c * fscale)
    wlp8_cs = np.roll(wlp8_c, tuple(np.array(shift, int)), axis=(0, 1))

    plt.subplots_adjust(wspace=0.01, hspace=0.01, left=0.2, right=0.8)
    fig, axes = plt.subplots(
        figsize=(16, 9),
        ncols=3,
        nrows=2,
        gridspec_kw={'left': 0.10, 'right': 0.90, 'wspace': 0.01, 'bottom': 0.02, 'top': 0.85},
    )

    axes[0, 0].imshow(wlm8_m, norm=norm, cmap=cm)
    axes[1, 0].imshow(wlp8_m, norm=norm, cmap=cm)

    axes[0, 1].imshow(wlm8_cs * fscale, norm=norm, cmap=cm)
    axes[1, 1].imshow(wlp8_cs * fscale, norm=norm, cmap=cm)

    axes[0, 2].imshow(wlm8_m - wlm8_cs * fscale, norm=norm, cmap=cm)
    axes[1, 2].imshow(wlp8_m - wlp8_cs * fscale, norm=norm, cmap=cm)

    nomask = np.ones_like(wlp8_m)
    rmserr_m1 = webbpsf.utils.rms(wlm8_m - wlm8_c * fscale, nomask)
    rmserr_m2 = webbpsf.utils.rms(wlp8_m - wlp8_c * fscale, nomask)
    axes[0, 2].text(10, 5, f'RMS counts: {rmserr_m1:.4}', color='white', fontsize=10)
    axes[1, 2].text(10, 5, f'RMS counts: {rmserr_m2:.4}', color='white', fontsize=10)

    axes[0, 0].set_ylabel('WLM8', fontsize=18)
    axes[1, 0].set_ylabel('WLP8', fontsize=18)
    axes[0, 0].set_title('Measured', fontsize=18)
    axes[0, 1].set_title('WAS PR Calc', fontsize=18)
    axes[0, 2].set_title('Difference\n(Measured - WAS PR)', fontsize=18)

    fig.suptitle(f"{hdul[0].header['CORR_ID']}     {hdul[0].header['TSTAMP']}     {fn}", fontsize=18, fontweight='bold')

    cax = plt.axes([0.92, 0.02, 0.02, 0.83])
    plt.colorbar(mappable=axes[0, 0].images[0], cax=cax, label='Surface Brightness [MJy/sr]')

    for ax in axes.flat:
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

    return fig


def plot_wfs_obs_delta(fn1, fn2, vmax_fraction=1.0, download_opds=True):
    """Display comparison of two weak lens observations

    This compares the actual measured WL data, not the derived wavefronts.
    Useful for quick checks of e.g. anything odd in the data, presence of a binary, etc.

    See also plot_phase_retrieval_crosscheck

    Parameters
    ----------
    fn1, fn2 : string
        Filenames of two OPD files retrieved from MAST.
        (Does not need to include full path, if files were downloaded by webbpsf)
    vmax_fraction : float
        Scale factor for setting the plot log scale vmax relative to the
        image peak pixel. vmax_fraction = 1.0 sets the scale vmax to the
        image max value, vmax_fraction = 0.1 sets it to 1/10th that, etc.

    """
    from skimage.registration import phase_cross_correlation

    if download_opds:
        _ = webbpsf.mast_wss.mast_retrieve_opd(fn1)
        _ = webbpsf.mast_wss.mast_retrieve_opd(fn2)

    opd, hdul1 = webbpsf.trending._read_opd(fn1)

    wlm8_m1 = hdul1[5].data
    wlp8_m1 = hdul1[10].data

    opd, hdul2 = webbpsf.trending._read_opd(fn2)

    wlm8_m2 = hdul2[5].data
    wlp8_m2 = hdul2[10].data

    cm = matplotlib.cm.inferno
    cm.set_bad(cm(0))

    fscale = wlm8_m1.sum() / wlm8_m2.sum()

    norm = matplotlib.colors.LogNorm(wlm8_m1.max() / 1e4 * vmax_fraction, wlm8_m1.max() * vmax_fraction)

    # Shift to register
    shift, error, diffphase = phase_cross_correlation(wlm8_m1, wlm8_m2 * fscale)
    wlm8_m2s = np.roll(wlm8_m2, tuple(np.array(shift, int)), axis=(0, 1))

    shift, error, diffphase = phase_cross_correlation(wlp8_m1, wlp8_m2 * fscale)
    wlp8_m2s = np.roll(wlp8_m2, tuple(np.array(shift, int)), axis=(0, 1))

    plt.subplots_adjust(wspace=0.01, hspace=0.01, left=0.2, right=0.8)
    fig, axes = plt.subplots(
        figsize=(16, 9),
        ncols=3,
        nrows=2,
        gridspec_kw={'left': 0.1, 'right': 0.9, 'wspace': 0.01, 'bottom': 0.05, 'top': 0.85},
    )

    axes[0, 0].imshow(wlm8_m1, norm=norm, cmap=cm)
    axes[1, 0].imshow(wlp8_m1, norm=norm, cmap=cm)

    axes[0, 1].imshow(wlm8_m2 * fscale, norm=norm, cmap=cm)
    axes[1, 1].imshow(wlp8_m2 * fscale, norm=norm, cmap=cm)

    axes[0, 2].imshow(wlm8_m1 - wlm8_m2s * fscale, norm=norm, cmap=cm)
    axes[1, 2].imshow(wlp8_m1 - wlp8_m2s * fscale, norm=norm, cmap=cm)

    nomask = np.ones_like(wlp8_m1)
    rmserr_m1 = webbpsf.utils.rms(wlm8_m1 - wlm8_m2s * fscale, nomask)
    rmserr_m2 = webbpsf.utils.rms(wlp8_m1 - wlp8_m2s * fscale, nomask)
    axes[0, 2].text(0, 0, f'RMS counts: {rmserr_m1:.4}', color='white', fontsize=12)
    axes[1, 2].text(0, 0, f'RMS counts: {rmserr_m2:.4}', color='white', fontsize=12)

    axes[0, 0].set_ylabel('WLM8', fontsize=18)
    axes[1, 0].set_ylabel('WLP8', fontsize=18)
    axes[0, 0].set_title(f'Measured: \n{os.path.basename(fn1)}', fontsize=18)
    axes[0, 1].set_title(f'Measured: \n{os.path.basename(fn2)}', fontsize=18)
    axes[0, 2].set_title('Difference\n ', fontsize=18)

    figure_title = (
        f"{hdul1[0].header['CORR_ID']}, {hdul1[0].header['TSTAMP'][:-3]}   vs.   "
        f"{hdul2[0].header['CORR_ID']}, {hdul2[0].header['TSTAMP'][:-3]}"
    )
    fig.suptitle(
        figure_title,
        fontsize=20,
        fontweight='bold',
    )

    cax = plt.axes([0.92, 0.02, 0.02, 0.83])
    plt.colorbar(mappable=axes[0, 0].images[0], cax=cax, label='Surface Brightness [MJy/sr]')

    for ax in axes.flat:
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

    return fig


def show_wfs_around_obs(filename, verbose='True'):
    """Make a helpful plot showing available WFS before and after some given science
    observation. This can be used to help inform how much WFE variability there was around that time.

    See also show_wfs_during_program.

    Parameters
    ----------
    filename : str
        A filename of some JWST data

    """

    header = fits.getheader(filename)

    def get_datetime_from_header(header):
        return astropy.time.Time(header['DATE-OBS'] + 'T' + header['TIME-OBS'])

    def vprint(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    # Retrieve header info, and WFS data before and after
    inst = webbpsf.instrument(header['INSTRUME'])
    inst.filter = header['filter']
    inst.set_position_from_aperture_name(header['APERNAME'])

    dateobs = get_datetime_from_header(header)
    vprint(f'File {filename} observed at {dateobs}')

    vprint('Retrieving WFS before that obs...', end='')
    inst.load_wss_opd_by_date(dateobs, choice='before', verbose=False)
    wfe_before = inst.get_wfe('total')
    wfe_before_dateobs = get_datetime_from_header(inst.pupilopd[0].header)
    vprint(f' WFS at {wfe_before_dateobs}')

    vprint('Retrieving WFS after that obs...', end='')
    inst.load_wss_opd_by_date(dateobs, choice='after', verbose=False)
    wfe_after = inst.get_wfe('total')
    wfe_after_dateobs = get_datetime_from_header(inst.pupilopd[0].header)

    vprint(f' WFS at {wfe_after_dateobs}')

    fnbase = os.path.basename(filename)

    # Setup axes
    fig = plt.figure(figsize=(16, 9), constrained_layout=False)

    gs = matplotlib.gridspec.GridSpec(2, 4, figure=fig, hspace=0.3, height_ratios=[1, 2])
    ax_t = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])
    ax3 = fig.add_subplot(gs[1, 2])
    ax4 = fig.add_subplot(gs[1, 3])

    # Plot and annotate timeline at top
    ax_t.plot_date([wfe_before_dateobs.plot_date, dateobs.plot_date, wfe_after_dateobs.plot_date], [0, 0, 0])
    ax_t.axhline(0, ls=':')

    ax_t.xaxis.set_major_locator(matplotlib.dates.DayLocator(interval=1))
    ax_t.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))
    ax_t.xaxis.set_minor_locator(matplotlib.dates.HourLocator(interval=1))
    ax_t.tick_params(which='major', width=2, length=5)

    ax_t.yaxis.set_visible(False)
    ax_t.set_ylim(-0.5, 1)
    ax_t.margins(0.1)
    ax_t.text(wfe_before_dateobs.plot_date, 0.2, 'WFS before', rotation=70, color='C0')
    ax_t.text(wfe_after_dateobs.plot_date, 0.2, 'WFS after', rotation=70, color='C0')
    ax_t.scatter(dateobs.plot_date, 0, marker='s', color='C2', zorder=10)
    ax_t.set_title(f'Wavefront Sensing around Observation\n {fnbase}', fontweight='bold', fontsize=18)

    ax_t.text(
        (wfe_before_dateobs + (dateobs - wfe_before_dateobs) / 2).plot_date,
        -0.2,
        f'{(dateobs-wfe_before_dateobs).to_value(u.day):.3f} days',
        horizontalalignment='center',
        color='black',
    )
    ax_t.text(
        (dateobs + (wfe_after_dateobs - dateobs) / 2).plot_date,
        -0.2,
        f'{(wfe_after_dateobs-dateobs).to_value(u.day):.3f} days',
        horizontalalignment='center',
        color='black',
    )

    # Not really an axis label for the top plot, but this is a convenient/easy way to
    # get this text into the middle of the figure:
    ax_t.set_xlabel(
        f" \nShowing WFE for {header['APERNAME']}  (inferred from WFS at NRCA3 FP1)", fontweight='bold', fontsize=18
    )

    # Compute linear weighted interpolated estimate WFS at time obs
    #  This is sort of a dirty trick...
    wfs_deltat = wfe_after_dateobs - wfe_before_dateobs
    weight_before = (wfe_after_dateobs - dateobs) / wfs_deltat
    weight_after = (dateobs - wfe_before_dateobs) / wfs_deltat
    wfe_weighted = wfe_before * weight_before + wfe_after * weight_after

    # one more annotation for above plot, using weights to be clever about spacing
    ax_t.text(
        dateobs.plot_date,
        0.2,
        f'Science obs:\n{fnbase}',
        rotation=0,
        color='C2',
        horizontalalignment='right' if weight_after > weight_before else 'left',
    )

    # Retrieve ap mask for overplotting OPDS with the borders nicely grayed out
    apmask = webbpsf.utils.get_pupil_mask()
    nanmask = np.ones_like(apmask)
    nanmask[apmask == 0] = np.nan

    # Plot the OPDs
    vmax = 0.15
    webbpsf.trending.show_opd_image(wfe_before * nanmask * 1e6, ax=ax1, vmax=vmax, fontsize=10)
    ax1.set_title('WFS Before\n ', color='C0', fontweight='bold')
    webbpsf.trending.show_opd_image(wfe_weighted * nanmask * 1e6, ax=ax2, vmax=vmax, fontsize=10)
    ax2.set_title('Time-weighted Linear \nEstimate at Obs Time', color='C2', fontweight='bold')
    webbpsf.trending.show_opd_image(wfe_after * nanmask * 1e6, ax=ax3, vmax=vmax, fontsize=10)
    ax3.set_title('WFS After\n ', color='C0', fontweight='bold')

    webbpsf.trending.show_opd_image((wfe_after - wfe_before) * nanmask * 1e6, ax=ax4, vmax=vmax, fontsize=10)
    ax4.set_title('Delta WFE\nAfter-Before', color='C1', fontweight='bold')


def show_wfs_during_program(
    program, verbose=False, ax=None, ref_wavefront_date=None, ref_wavefront_visit=None, start_date=None, end_date=None
):
    """Show WFS data for the entire time interval in which a given program was observed.

    Plots time series of the WFS measured RMS WFE as seen in NIRCam, the start and end times
    of all the observations in that program, and the delta wavefront RMS relative to either
    the median wavefront during that time period, or to a specifed date.

    See also show_wfs_around_obs.

    Parameters
    ----------
    program : str or int
        Program ID number
    verbose : str
        be more verbose in output?
    ax : matplotlib.Axes instance, or None
        an existing Axes to plot into, or else a new one will be created.
    ref_wavefront_date : date-like str or None
        Optional, to specify which date's wavefront sensing should be used as
        the reference wavefront for computation of the plotted delta WFE RMS.
        The closest wavefront in time to the specified date will be used.
        If not set, the median wavefront over the entire time period will be used.
    ref_wavefront_visit : str
        Another way to specify which date's wavefront sensing should be used as the
        reference. Provide a science visit ID from the program (e.g. "V01234002001")
        and the closest wavefront to the observation start time of that visit will be used.
    start_date, end_date : strings or astropy.time.Time
        Start and/or end dates for the time period to display. If not set,
        reasonable default values will be computed that include all observations
        for that science program plus some padding time on either side.
    """

    # Query mast for when the observations took place
    science_visit_table = webbpsf.mast_wss.query_program_visit_times(program, verbose=verbose)
    science_visit_table.sort(keys=['start_mjd'])
    if verbose:
        for row in science_visit_table:
            print(f" Found {row['visit_id']} starting at {row['start_mjd'].iso}")

    # Figure out reasonable start and end dates for the time interval to display,
    # or use values provided by the user
    sci_duration = science_visit_table['start_mjd'].max() - science_visit_table['start_mjd'].min()
    plot_padding_time_range = max(sci_duration.to(u.day) * 0.2, 4 * u.day)
    if start_date is None:
        start_date = science_visit_table['start_mjd'].min() - plot_padding_time_range
    else:
        start_date = astropy.time.Time(start_date)
    if end_date is None:
        end_date = science_visit_table['end_mjd'].max() + plot_padding_time_range
    else:
        end_date = astropy.time.Time(end_date)

    # Look up wavefront sensing and mirror move corrections for that range
    opdtable = get_opdtable_for_daterange(start_date, end_date)

    # Iterate over the WFS measurements to retrieve the OPDs and RMS WFE
    wfs_dates = []
    rms_obs = []
    opds = []
    for row_index in range(len(opdtable)):
        opd_fn = opdtable[row_index]['fileName']
        try:
            opd, opd_hdul = webbpsf.trending._read_opd(opd_fn)
        except FileNotFoundError:
            webbpsf.mast_wss.mast_retrieve_opd(opd_fn, verbose=verbose)
            opd, opd_hdul = webbpsf.trending._read_opd(opd_fn)

        opds.append(opd)
        wfs_dates.append(opdtable[row_index]['date'])
        rms_obs.append(opd_hdul[1].header['RMS_WFE'] * 1000)
    wfs_dates_array = astropy.time.Time(wfs_dates, format='isot')
    opd_array = np.asarray(opds)

    # Compute delta WFE, either relative to median or a specifed date
    if ref_wavefront_visit is not None:
        for row in science_visit_table:
            if row['visit_id'] == ref_wavefront_visit:
                ref_wavefront_date = row['start_mjd'].iso
                if verbose:
                    print(f" Will compare wavefronts relative to {row['visit_id']} starting at {row['start_mjd'].iso}")

    if ref_wavefront_date:
        refdate = astropy.time.Time(ref_wavefront_date)
        delta_times = wfs_dates_array - refdate
        closest = np.argmin(np.abs(delta_times))
        reference_opd = opds[closest]
        ref_label = 'WFS near ' + ref_wavefront_date
        if verbose:
            print(f"Computing delta OPDs relative to date {ref_wavefront_date}, using opd {opdtable[closest]['fileName']}")
    else:
        median_opd = np.median(opd_array, axis=0)
        reference_opd = median_opd
        ref_label = 'median wavefront'
        if verbose:
            print('Computing delta OPDs relative to median OPD over that time period.')

    delta_opds = opd_array - reference_opd
    mask = opd_array.sum(axis=0) != 0
    delta_rmses = [webbpsf.utils.rms(d, mask=mask) * 1000 for d in delta_opds]

    # Plot!
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), ncols=1, nrows=1)

    ax.plot_date(
        wfs_dates_array.plot_date, rms_obs, '+', color='C1', ls='-', label='Measured RMS Wavefront Error at NIRCam NRCA3'
    )
    ax.plot_date(
        wfs_dates_array.plot_date,
        delta_rmses,
        'none',
        color='C0',
        ls='--',
        label=f'RMS Delta Wavefront relative to {ref_label}',
    )

    plot_sci_y = 25
    ax.scatter(
        science_visit_table['start_mjd'].plot_date,
        [plot_sci_y] * len(science_visit_table),
        s=150,
        marker='*',
        color='black',
        label=f'Program {program} observations',
    )
    for row in science_visit_table:
        ax.fill_betweenx([0, 120], row['start_mjd'].plot_date, (row['end_mjd']).plot_date, color='gray', alpha=0.2)
        ax.text(row['start_mjd'].plot_date, plot_sci_y + 4, row['visit_id'], rotation=45, fontsize='small', clip_on=True)

    ax.set_ylim(0, 120)
    ax.legend(framealpha=0.99)
    ax.xaxis.set_major_formatter(matplotlib.dates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
    ax.set_ylabel('Measured Wavefront Error RMS [nm]')
    ax.set_xlim(start_date.plot_date, end_date.plot_date)


def delta_wfe_around_time(datetime, plot=True, ax=None, vmax=0.05, return_filenames=False):
    """Return the delta OPD in the period around some specified time

    Inferred from the available WFS measurements before and after that time.

    Parameters
    ----------
    datetime : string
        date time in ISO8601 format, or similar
    plot : bool
        make plot?
    vmax : float
        vmax for plot, in microns.

    Returns the delta OPD, in units of microns.

    """

    prev_opd_fn, post_opd_fn, prev_delta_t, post_delta_t = webbpsf.mast_wss.mast_wss_opds_around_date_query(
        datetime, verbose=False
    )

    prev_opd, prev_hdul = _read_opd(prev_opd_fn)
    post_opd, post_hdul = _read_opd(post_opd_fn)

    delta_opd = post_opd - prev_opd
    mask = prev_opd != 0

    nanmask = np.ones_like(prev_opd)
    nanmask[~mask] = np.nan

    if plot:
        if ax is None:
            ax = plt.gca()
        show_opd_image(delta_opd * nanmask, ax=ax, vmax=vmax)
        plt.colorbar(mappable=ax.images[0], label='WFE [microns]')

        ax.set_title(f'$\\Delta$WFE in the {post_delta_t - prev_delta_t:.2f} d around {datetime}')
        ax.set_xlabel(f'{post_opd_fn} - \n{prev_opd_fn}')
        ax.set_xticks([])
        ax.xaxis.set_visible(True)

    return delta_opd


# Functions for image comparisons
def show_nrc_ta_img(visitid, ax=None, return_handles=False):
    """Retrieve and display a NIRCam target acq image"""

    hdul = webbpsf.mast_wss.get_visit_nrc_ta_image(visitid)

    ta_img = hdul['SCI'].data
    mask = np.isfinite(ta_img)
    rmean, rmedian, rsig = astropy.stats.sigma_clipped_stats(ta_img[mask])
    bglevel = rmedian

    vmax = np.nanmax(ta_img) - bglevel
    cmap = matplotlib.cm.viridis.copy()
    cmap.set_bad('orange')

    norm = matplotlib.colors.AsinhNorm(
        linear_width=vmax * 0.003,
        vmax=vmax,  # vmin=0)
        vmin=-1 * rsig,
    )

    if ax is None:
        ax = plt.gca()
    ax.imshow(ta_img - bglevel, norm=norm, cmap=cmap, origin='lower')
    ax.set_title(f"NIRCam TA on {visitid}\n{hdul[0].header['DATE-OBS']}")
    ax.set_ylabel('[Pixels]')
    ax.text(0.05, 0.9, hdul[0].header['TARGPROP'], color='white', transform=ax.transAxes)

    if return_handles:
        return hdul, ax, norm, cmap, bglevel


def nrc_ta_image_comparison(visitid, verbose=False, show_centroids=False):
    """Retrieve a NIRCam target acq image and compare to a simulation

    Parameters:
    -----------
    visitid : string
        Visit ID. Should be one of the WFSC visits, starting with a NIRCam target acq, or at least
        some other sort of visit that begins with a NIRCam target acquisition.
    """
    from skimage.registration import phase_cross_correlation

    fig, axes = plt.subplots(figsize=(10, 5), ncols=3)

    # Get and plot the observed TA image
    hdul, ax, norm, cmap, bglevel = show_nrc_ta_img(visitid, ax=axes[0], return_handles=True)
    im_obs = hdul['SCI'].data.copy()
    im_obs_err = hdul['ERR'].data
    im_obs_dq = hdul['DQ'].data

    im_obs_clean = im_obs.copy()
    im_obs_clean[im_obs_dq & 1] = np.nan  # Mask out any DO_NOT_USE pixels.
    im_obs_clean = astropy.convolution.interpolate_replace_nans(im_obs_clean, kernel=np.ones((5, 5)))

    # Make a matching sim
    nrc = webbpsf.setup_sim_to_match_file(hdul, verbose=False)
    opdname = nrc.pupilopd[0].header['CORR_ID'] + '-NRCA3_FP1-1.fits'
    if verbose:
        print('Calculating PSF to match that TA image...')
    psf = nrc.calc_psf(fov_pixels=im_obs.shape[0])

    # Align and Shift:
    im_sim = psf['DET_DIST'].data  # Use the extension including distortion and IPC

    shift, _, _ = phase_cross_correlation(im_obs_clean, im_sim, upsample_factor=32)
    if verbose:
        print(f'Shift to register sim to data: {shift} pix')
    im_sim_shifted = scipy.ndimage.shift(im_sim, shift, order=5)

    # figure out the background level and scale factor
    scalefactor = np.nanmax(im_obs_clean) / im_sim.max()
    if verbose:
        print(f'Scale factor to match sim to data: {scalefactor}')

    im_sim_scaled_aligned = im_sim_shifted * scalefactor

    # Plot
    if show_centroids:
        # OSS CENTROIDS ###
        # First, see if we can retrieve the on-board TA centroid measurment from the OSS engineering DB in MAST
        try:
            import misc_jwst  # Optional dependency, including engineering database access tools

            # If we have that, retrieve the log for this visit, extract the OSS centroids, and convert to same
            # coordinate frame as used here:
            osslog = misc_jwst.engdb.get_ictm_event_log(hdul[0].header['VSTSTART'], hdul[0].header['VISITEND'])
            try:
                oss_cen = misc_jwst.engdb.extract_oss_TA_centroids(osslog, 'V' + hdul[0].header['VISIT_ID'])
                # Convert from full-frame (as used by OSS) to detector subarray coords:
                oss_cen_sci = nrc._detector_geom_info.aperture.det_to_sci(*oss_cen)
                oss_cen_sci_pythonic = np.asarray(oss_cen_sci) - 1  # convert from 1-based pixel indexing to 0-based
                oss_centroid_text = f'\n   OSS centroid: {oss_cen_sci_pythonic[0]:.2f}, {oss_cen_sci_pythonic[1]:.2f}'
                axes[0].scatter(oss_cen_sci_pythonic[0], oss_cen_sci_pythonic[1], color='0.5', marker='+', s=50)
                axes[0].text(
                    oss_cen_sci_pythonic[0],
                    oss_cen_sci_pythonic[1],
                    'OSS  ',
                    color='0.9',
                    verticalalignment='center',
                    horizontalalignment='right',
                )
                if verbose:
                    print(f'OSS centroid on board:  {oss_cen}  (full det coord frame, 1-based)')
                    print(
                        f'OSS centroid converted: {oss_cen_sci_pythonic}',
                        f'(sci frame in {nrc._detector_geom_info.aperture.AperName}, 0-based)'
                    )
                    full_ap = nrc.siaf[nrc._detector_geom_info.aperture.AperName[0:5] + '_FULL']
                    oss_cen_full_sci = np.asarray(full_ap.det_to_sci(*oss_cen)) - 1
                    print(f'OSS centroid converted: {oss_cen_full_sci}  (sci frame in {full_ap.AperName}, 0-based)')

            except RuntimeError:
                if verbose:
                    print('Could not parse TA coordinates from log. TA may have failed?')
                oss_cen_sci_pythonic = None

            # WCS COORDINATES ###
            import jwst.datamodels

            model = jwst.datamodels.open(hdul)
            targ_coords = astropy.coordinates.SkyCoord(model.meta.target.ra, model.meta.target.dec, frame='icrs', unit=u.deg)
            targ_coords_pix = model.meta.wcs.world_to_pixel(targ_coords)  # returns x, y
            axes[0].scatter(targ_coords_pix[0], targ_coords_pix[1], color='magenta', marker='+', s=50)
            axes[0].text(
                targ_coords_pix[0],
                targ_coords_pix[1] + 2,
                'WCS',
                color='magenta',
                verticalalignment='bottom',
                horizontalalignment='center',
            )
            axes[0].text(
                0.95,
                0.18,
                f'From WCS: {targ_coords_pix[0]:.2f}, {targ_coords_pix[1]:.2f}',
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=axes[0].transAxes,
                color='white',
            )

            if verbose:
                print(f'Star coords from WCS: {targ_coords_pix}')
                if oss_cen_sci_pythonic is not None:
                    print(f'WCS offset =  {np.asarray(targ_coords_pix) - oss_cen_sci_pythonic} pix')

        except ImportError:
            oss_centroid_text = ''

        # WEBBPSF CENTROIDS ###
        cen = webbpsf.fwcentroid.fwcentroid(im_obs_clean)
        axes[0].scatter(cen[1], cen[0], color='red', marker='+', s=50)
        axes[0].text(cen[1], cen[0], '  webbpsf', color='red', verticalalignment='center')
        axes[0].text(
            0.95,
            0.05,
            f' webbpsf Centroid: {cen[1]:.2f}, {cen[0]:.2f}' + oss_centroid_text,
            horizontalalignment='right',
            verticalalignment='bottom',
            transform=axes[0].transAxes,
            color='white',
        )

    axes[1].imshow(im_sim_scaled_aligned, norm=norm, cmap=cmap, origin='lower')
    axes[1].set_title(f'Simulated PSF in F212N\nusing {opdname}')

    diffim = im_obs - bglevel - im_sim_scaled_aligned

    dofs = np.isfinite(diffim).sum() - 4  # 4 estimated parameters: X and Y offsets, flux scaling, background level
    reduced_chisq = np.nansum(((diffim / im_obs_err) ** 2)) / dofs

    axes[2].imshow(diffim, cmap=cmap, norm=norm, origin='lower')
    axes[2].set_title('Difference image\nafter alignment and scaling')
    axes[2].text(
        0.05,
        0.9,
        f'$\\chi^2_r$ = {reduced_chisq:.3g}' + ('  Alert, not a good fit!' if (reduced_chisq > 1.5) else ''),
        transform=axes[2].transAxes,
        color='white' if (reduced_chisq < 1.5) else 'yellow',
    )

    for ax in axes:
        fig.colorbar(ax.images[0], ax=ax, orientation='horizontal', label=hdul['SCI'].header['BUNIT'])

    plt.tight_layout()

    outname = f'nrc_ta_comparison_{visitid}.pdf'
    plt.savefig(outname)
    print(f' => {outname}')
