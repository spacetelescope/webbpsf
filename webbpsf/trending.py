import calendar
import datetime
import functools
import os

import astropy
import astropy.io.fits as fits
import astropy.time
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scipy.interpolate

import poppy
import webbpsf

def _read_opd(filename):
    """Trivial utilty function to read OPD from a WSS-output FITS file"""
    full_file_path = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'MAST_JWST_WSS_OPDs', filename)
    opdhdu = fits.open(full_file_path)
    opd = opdhdu[1].data
    return opd, opdhdu

def get_datetime_utc(opdhdul, return_as='string'):
    """Utilty function. Retrieve the date time either as a string with ISO date, or an astropy time object"""
    if return_as == 'string':
        return opdhdul[0].header['DATE-OBS'] + " " + opdhdul[0].header['TIME-OBS'][:-4] # leave off milliseconds
    elif return_as == 'astropy':
        return astropy.time.Time(opdhdul[0].header['DATE-OBS'] + "T" + opdhdul[0].header['TIME-OBS'],
                                format='isot', scale='utc')
    else:
        raise ValueError(f"Invalid value for return_as={return_as}")


def wavefront_time_series_plot(opdtable, start_date=None, end_date=None, label_visits=True, label_events=True):
    """ Make a time series plot of total WFS versus time

    Parameters
    ----------
    opdtable : astropy.table.Table
        OPD table, retrieved from MAST. See functons in mast_wss.py
    start_date, end_date : datetime.datetime objects
        Start and end dates for the plot time range. Default is March 2022 to present.
    label_visits : bool
        Label program_id:visit_id for each WFS visit.

    Returns
    -------
    Nothing, but makes a plot

    """

    # Plot figure showing WFE vs time

    rmses = []

    dates = astropy.time.Time(opdtable['date'], format='isot')
    pre_or_post = []

    for row in opdtable:
        full_file_path = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'MAST_JWST_WSS_OPDs', row['fileName'])
        rmses.append(fits.getheader(full_file_path, ext=1)['RMS_WFE'])
        pre_or_post.append(webbpsf.mast_wss.infer_pre_or_post_correction(row))

    where_pre = ['pre' in a for a in pre_or_post]
    where_post = ['post' in a for a in pre_or_post]

    events = {'2022-03-22T22:55:00': ('Cooler State 4', 'red'),
              '2022-04-10T17:24:00': ('Cooler State 6', 'red'),
              '2022-04-07T20:13:00': ('Cooler State 5', 'red'),
              #'2022-04-15T06:33:00': ("NIRISS FPA Heater on (NIS-24)", 'orange'),
              # '2022-04-10T04:30':  ('NIRSpec FPA CCH to level 30','orange'),
              #'2022-04-10T01:00': ('NIRSpc bench & FPA heaters adjust', 'orange'),  # to level 10
              #'2022-04-11T16:00:00': ('FGS trim heater 0 to 4', 'orange'),
              '2022-04-23T06:30:00': ('OTE alignment complete', 'green'),
              '2022-05-12T18:30:00': ('Thermal slew cold attitude start', 'blue'),
              '2022-05-20T06:30:00': ('Thermal slew cold attitude end', 'blue'),
              '2022-05-23T00:00:00': ('Larger micrometeorite strike on C3', 'red'),
              '2022-06-19T18:00:00': ('Coarse move of C3 for astigmatism correction', 'red'),
              #'2022-06-27T00:00:00': ('NIRSpec safing, not in thermal control', 'orange'),
              '2022-07-12T00:00:00': ('Large outlier tilt event on B5+C5', 'orange')
              }

    plt.figure(figsize=(16, 12))

    rms_nm = np.asarray(rmses) * 1000

    routine_pids = [1163, 2586, 2724, 2725, 2726]  # commissioning OTE-26 and cycle 1 maintenance

    is_routine = np.asarray([int(v[1:6]) in routine_pids for v in opdtable[where_pre]['visitId']])

    # Plot all, with connecting line
    plt.plot_date(dates.plot_date, rms_nm, ls='-', color='gray', )
    # Plot maintenance visits
    plt.plot_date(dates[where_pre][is_routine].plot_date, rms_nm[where_pre][is_routine], 'o', xdate=True,
                  label="WF Maintenance Sensing",
                  color='C0')
    # Plot other visits (MIMF etc)
    plt.plot_date(dates[where_pre][~is_routine].plot_date, rms_nm[where_pre][~is_routine], 'o', xdate=True,
                  label="Other Sensing (MIMF etc)",
                  color='purple')
    # Plot corrections.
    plt.plot_date(dates[where_post].plot_date, rms_nm[where_post], 'v', xdate=True, label="Wavefront Corrections",
                  color='C2', )

    ax = plt.gca()
    ax.set_ylabel('Observatory WFE (OTE+NRC)\n[nm rms]', fontweight='bold', fontsize=15)
    ax.set_xlabel("Date, UTC", fontweight='bold', fontsize=15)

    # ymin, ymax = 40, 200
    ymin, ymax = 0, 250
    ax.set_ylim(ymin, ymax)
    ax.axhline(70, ls=":", color='gray')
    ax.axhline(80, ls=":", color='orange')
    ax.axhline(100, ls=":", color='gray')

    if start_date is None:
        start_date = datetime.datetime(2022, 3, 20, 0)
    if end_date is None:
        end_date = datetime.datetime.now() + datetime.timedelta(days=7)
    # start_date = datetime.datetime(2022, 2, 20, 0)
    ax.set_xlim(start_date, end_date)

    ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(interval=1))
    ax.xaxis.set_minor_locator(matplotlib.dates.DayLocator())
    ax.tick_params('x', length=10)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    # label events
    if label_events:
        for timestamp, (event, color) in events.items():
            d = astropy.time.Time(timestamp, format='isot')
            plt.axvline(d.plot_date, color=color, ls=':', alpha=0.5)
            ax.text(d.plot_date + 0.25, ymax * 0.95, event, color=color, rotation=90, verticalalignment='top', alpha=0.7)

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
            short_visit = f"  {int(visit[1:6])}:{int(visit[6:9])}"
            if d_pre.datetime > start_date and rms < ymax:
                plt.text(d_pre.plot_date, rms, short_visit, rotation=65, color='C0', fontsize=10)

    plt.plot(dates.plot_date[where_pre][is_routine], rms_nm[where_pre][is_routine], ls='none', alpha=0.5)
    plt.plot(dates.plot_date[where_post], rms_nm[where_post], ls='none', color='C2', alpha=0.5)

    plt.title("Observatory WFE over time at NIRCam Field Point 1", fontweight='bold', fontsize=15)

    plt.legend(loc='upper right')


def wfe_histogram_plot(opdtable, start_date=None, end_date=None, thresh=None):
    """ Plot histogram and cumulative histogram of WFE over some time range.
    """

    if start_date is None:
        start_date = astropy.time.Time('2022-07-16')
    if end_date is None:
        end_date = astropy.time.Time.now()



    opdtable0 = webbpsf.mast_wss.deduplicate_opd_table(opdtable)
    opdtable1 = webbpsf.mast_wss.filter_opd_table(opdtable0, start_time=start_date, end_time=end_date)

    webbpsf.mast_wss.download_all_opds(opdtable1)

    # Retrieve all RMSes, from the FITS headers. 
    # These are observatory WFE (OTE + NIRCam), at the WFS sensing field point
    rmses=[]
    mjds = []
    for row in opdtable1:
        full_file_path = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'MAST_JWST_WSS_OPDs', row['fileName'])
        rmses.append(fits.getheader(full_file_path, ext=1)['RMS_WFE'])
        mjds = opdtable1['date_obs_mjd']

    dates = astropy.time.Time(opdtable1['date'], format='isot')

    # Interpolate those RMSes into an even grid over time
    interp_fn = scipy.interpolate.interp1d(mjds, rmses, kind='linear')

    mjdrange = np.linspace(np.min(mjds), np.max(mjds), 2048)
    interp_rmses = interp_fn(mjdrange)

    # Plot
    fig, axes = plt.subplots(figsize=(16,12), nrows=3, gridspec_kw = {'hspace':0.3})

    axes[0].plot_date(dates.plot_date, np.asarray(rmses)*1e3, 'o', ls='-')


    #plt.plot(mjdrange, interp_rmses, marker='+')
    axes[0].set_xlabel("Date")
    axes[0].set_ylabel("RMS WFE\n(OTE+NIRCam)", fontweight='bold')
    axes[0].set_title(f"Observatory WFE from {start_date.isot[0:10]} to {end_date.isot[0:10]}",
                     fontsize=14, fontweight='bold')

    if thresh:
        axes[0].axhline(thresh, color='C2')

    nbins=100

    axes[1].set_title(f"Observatory WFE Histogram from {start_date.isot[0:10]} to {end_date.isot[0:10]}",
                     fontsize=14, fontweight='bold')

    axes[1].hist(interp_rmses*1e3, density=True, bins=nbins)
    axes[1].set_ylabel("Fraction of time with this WFE", fontweight='bold')

    axes[2].hist(interp_rmses*1e3, density=True, bins=nbins, cumulative=1, histtype='step', lw=3, color='C1');
    axes[2].set_ylabel("Cumulative fraction of time\nwith this WFE or better", fontweight='bold')
    axes[2].set_title(f"Observatory WFE Cumulative Histogram from {start_date.isot[0:10]} to {end_date.isot[0:10]}",
                     fontsize=14, fontweight='bold')

    axes[2].set_xlabel("RMS Wavefront Error [nm]")
    axes[2].set_ylim(0,1)
    axes[2].set_xlim(60, interp_rmses.max()*1e3)

    if thresh: 
        for i in [1,2]:
            axes[i].axvline(thresh, color='C2')    
        fractime = (interp_rmses*1e3 < thresh).sum()/len(interp_rmses)
        axes[2].text(thresh+0.5, 0.1, 
                     f"{fractime*100:.1f}% of the time has WFE < {thresh}", color='C2',
                    fontweight='bold', fontsize=14)    



##### Wavefront Drifts Plot #####

def show_opd_image(array, vmax=0.300, ax=None, title=None, colorbar=False, labelrms=True, mask=None,
                   deltatime_hrs=None, fontsize=12, maskc3=None):
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
        ax.text(0.02, 0.02, f"{rms(array, mask) * 1000:.1f} nm rms", color='yellow', fontweight='bold',
                fontsize=fontsize, verticalalignment='bottom',
                transform=ax.transAxes)
        if maskc3 is not None:
            ax.text(array.shape[0]-5, 5,
                    f"ex. C3:\n{rms(array, maskc3)*1000:.1f} nm rms", color='lightgray', fontweight='bold',
                    fontsize=fontsize,
                    horizontalalignment='right')

    if deltatime_hrs is not None:
        # ax.text(array.shape[0], 5,
        ax.text(0.995, 0.02,
                f"{rms(array, mask) * 1000 / deltatime_hrs:.2f} nm/hr ",
                color='cyan', fontweight='bold', fontsize=fontsize,
                horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
        if maskc3 is not None:
            ax.text(array.shape[0]-5, -5,
                   f"{rms(array, maskc3)/deltatime_hrs*1000:.2f} nm/hr ",
                   color='lightgray', fontweight='bold', fontsize=fontsize,
                   horizontalalignment='right', verticalalignment='top')


    if colorbar:
        plt.colorbar(mappable=im, label='nanometers rms')

    if title is not None:
        ax.set_title(title, fontsize=16, fontweight='bold', )

    return ax


# Single measurement trending plot

@functools.lru_cache
def _get_mimf2_focus_offset_model(npix=256):
    # Precompute a model for the MIMF2 focus adjustment
    ote = webbpsf.opds.OTE_Linear_Model_WSS(npix=npix)
    ote.move_sm_local(piston=-1.789125, trans_unit='micron')  # Piston from M2022042201_sur.xml
    mimf2_focus_offset = ote.opd * 1e6  # convert meters to microns
    return mimf2_focus_offset


def single_measurement_trending_plot(opdtable, row_index=-1, verbose=True, vmax=0.15, ignore_missing=True, subtract_target=True):
    """Wavefront trending plot for a single measurement

    Parameters
    ----------

    opdtable : astropy.table.Table
        Table of available OPDs, as returned by retrieve_mast_opd_table()
    row_index : int
        Index into that table. Which row to make a plot for?
    verbose: bool
        be more verbose in output?
    vmax : float
        Image display scale max for OPD, in microns. Defaults to 0.15 microns = 150 nanometers


    """

    mimf2_focus_offset = _get_mimf2_focus_offset_model()

    filename = opdtable[row_index]['fileName']
    if verbose:
        print(f"Generating trending plot for {filename}")

    # Read current measurement
    opd, opdhdu = _read_opd(filename)
    mask = opd != 0
    opd[~ mask] = np.nan

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
            print(
                "     Immediately prior OPD file appears to be in F187N; skipping that and looking at next earlier OPD.")
            prior_opd_index = prior_opd_index - 1
    prev_filename = opdtable[prior_opd_index]['fileName']
    if verbose:
        print("     Previous OPD is:", prev_filename)

    prev_opd, prev_opd_hdu = _read_opd(prev_filename)

    if subtract_target:
        if verbose:
            print("     Subtracting NIRCam SI WFE target phase map")

        # Get WSS Target Phase Map
        was_targ_file = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'NIRCam', 'OPD', 'wss_target_phase.fits')
        target_1024 = astropy.io.fits.getdata(was_targ_file)
        target_256 = poppy.utils.krebin(target_1024, (256, 256)) /16   # scale factor for rebinning w/out increasing values

        opd -= target_256
        prev_opd -= target_256


    # Compute deltas and decompose
    deltatime = get_datetime_utc(opdhdu, return_as='astropy') - get_datetime_utc(prev_opd_hdu, return_as='astropy')
    deltatime_hrs = deltatime.to_value(u.hour)
    delta_opd = (opd - prev_opd)
    fit, coeffs = webbpsf.opds.decompose_opd_segment_PTT(delta_opd)

    # Compare to reference nominal OPD
    # For this we use the ~ best correction, immediately prior to MIMF-2 measurement
    # And we add in a model for the MIMF2 focus offset, applied just after that.

    ref_row = opdtable[opdtable['visitId'] == 'V01163030001']
    _, ref_opd_hdu = _read_opd(ref_row['fileName'][0])
    ref_opd = ref_opd_hdu[1].data + mimf2_focus_offset

    if subtract_target:
        ref_opd -= target_256

    # Read associated post-correction measurment, if present
    if opdtable[row_index]['is_pre_correction']:
        show_post_move = True
        meas_title = 'Measurement (Pre Move)'
        match = opdtable[(opdtable['visitId'] == opdtable[row_index]['visitId']) & opdtable['is_post_correction']][0]
        post_opd, post_opd_hdu = _read_opd(match['fileName'])

        if subtract_target:
            post_opd -= target_256

        # Compare the post-move sensing to the reference OPD
        diff2_title = "post-move - reference"
        delta_opd2 = post_opd - ref_opd
        delta_opd2[~ mask] = np.nan

        # infer the correction ID
        # based on the assumption that the move probably came from the immediately prior WSS session
        prior_corr_id = opdtable[row_index - 1]['corr_id']
        if verbose:
            print("   That measurement has a correction, therefore reading in post-move sensing data and SUR.")
            print(f"   Inferred prior correction ID is {prior_corr_id}")
        sur_fn = f'{prior_corr_id}_sur.xml'
    else:
        show_post_move = False
        meas_title = 'Measurement'
        # Compare this latest sensing measurement to the reference OPD
        diff2_title = "latest - reference"
        delta_opd2 = opd - ref_opd

    fit2, coeffs = webbpsf.opds.decompose_opd_segment_PTT(delta_opd2)
    ############## Plotting
    # Plot setup
    # fig, axes = plt.subplots(figsize=(8.5,11), nrows=3, ncols=4)
    fig, axes = plt.subplots(figsize=(16, 12), nrows=3, ncols=4,
                             gridspec_kw={'top': 0.9, 'bottom': 0.01, 'hspace': 0.2,
                                          'left': 0.05, 'right': 0.87},
                             )

    title = f"{filename}    {visit}     {get_datetime_utc(opdhdu)}"
    if subtract_target:
        title += "\nNIRCam FP1 Target Phase Map Subtracted"
    plt.suptitle(title, fontweight='bold', fontsize=18)


    ####### Row 1: Latest measurement, and correction if present
    fontsize=11
    # Panel 1: latest OPD
    iax = axes[0, 0]
    show_opd_image(opd, ax=iax, title=None, vmax=vmax, mask=mask, maskc3=mask_without_C3, fontsize=fontsize)
    iax.set_title(f"{meas_title}\n{get_datetime_utc(opdhdu)}", fontsize=fontsize*1.2, fontweight='bold')

    if show_post_move:
        iax = axes[0, 1]
        show_opd_image(post_opd * nanmask, ax=iax, title=None, vmax=vmax, mask=mask, maskc3=mask_without_C3, fontsize=fontsize)
        iax.set_title(f"Measurement (Post Move)\n{get_datetime_utc(post_opd_hdu)}", fontsize=fontsize*1.2, fontweight='bold')

        sur_opd = webbpsf.opds.sur_to_opd(sur_fn, ignore_missing=ignore_missing)
        # cosmetic: handle masking slightly differently here to accomodate slightly different edge pixels
        sur_mask = sur_opd != 0
        surnanmask = nanmask.copy()
        surnanmask[sur_opd == 0] = np.nan

        iax = axes[0, 2]
        show_opd_image(sur_opd * surnanmask, ax=iax, title=None, vmax=vmax, mask=mask)
        iax.set_title(f"Mirror Move Commanded\nSUR {sur_fn}", fontsize=fontsize)

        if np.all(np.isnan(sur_opd)):
            iax.text(0.5, 0.5, f"SUR file not found:\n{sur_fn}\n(Download of these is \nnot automated yet)", color='white', transform=iax.transAxes)
        # iax.text(128, 128, f"SUR retrieval\nNot yet implemented", color='black', fontweight='bold',
        #        horizontalalignment='center')

        iax = axes[0, 3]
        show_opd_image((post_opd - opd) * nanmask, ax=iax, title=None, vmax=vmax, mask=mask, maskc3=mask_without_C3)
        iax.set_title(f"Mirror Move Measured\nDelta WFE", fontsize=fontsize)

    else:
        for ax in axes[0, 1:4]:
            ax.set_visible(False)
        fig.text(0.55, 0.77, "Sensing-only visit. No mirror moves.", alpha=0.3,
                 horizontalalignment='center', fontsize=fontsize)
    ####### Row 2
    # Compare to immedaite prior OPD

    # Panel 2-1: prior OPD
    iax = axes[1, 0]
    show_opd_image(prev_opd * nanmask, ax=iax, title='Prior WFS', vmax=vmax, maskc3=mask_without_C3, fontsize=fontsize)
    iax.set_title(f"Prior Measurement\n{get_datetime_utc(prev_opd_hdu)}", fontsize=fontsize*1.2, fontweight='bold')

    # Panel 2-2: difference
    iax = axes[1, 1]
    show_opd_image(delta_opd, ax=iax, vmax=vmax, deltatime_hrs=deltatime_hrs, fontsize=fontsize)
    iax.set_title(f"Difference\n{deltatime.to(astropy.units.hour):.2f}", fontsize=fontsize*1.1)

    # Panel 2-3: proposed correction
    iax = axes[1, 2]
    show_opd_image(fit, ax=iax, vmax=vmax, deltatime_hrs=deltatime_hrs, fontsize=fontsize)
    iax.set_title(f"Controllable Modes\nin difference", fontsize=fontsize*1.1)

    # Panel 2-4:
    iax = axes[1, 3]
    show_opd_image(delta_opd - fit, ax=iax, vmax=vmax, fontsize=fontsize)
    iax.set_title(f"High order WFE\nin difference", fontsize=fontsize*1.1)

    ####### Row 3

    # Panel 3-1: ref OPD
    iax = axes[2, 0]
    show_opd_image(ref_opd * nanmask, ax=iax, title='Prior WFS', vmax=vmax, mask=mask, fontsize=fontsize)
    iax.set_title(f"Reference Measurement\n(from MIMF2)", fontsize=fontsize*1.1)

    # Panel 3-2: difference
    iax = axes[2, 1]
    show_opd_image(delta_opd2 * nanmask, ax=iax, vmax=vmax, mask=mask, fontsize=fontsize)
    iax.set_title(f"Difference\n{diff2_title}", fontsize=fontsize*1.1)

    # Panel 3-3: proposed correction
    iax = axes[2, 2]
    show_opd_image(fit2, ax=iax, vmax=vmax, mask=mask, fontsize=fontsize)
    iax.set_title(f"Controllable Modes\nin difference", fontsize=fontsize*1.1)

    # Panel 3-4:
    iax = axes[2, 3]
    show_opd_image(delta_opd2 - fit2, ax=iax, vmax=vmax, maskc3=mask_without_C3, fontsize=fontsize)
    iax.set_title(f"High order WFE\nin difference", fontsize=fontsize*1.1)

    cax = fig.add_axes([0.91, 0.41, 0.01, 0.4])

    cmap = matplotlib.cm.RdBu_r
    mappable_nm = matplotlib.cm.ScalarMappable(cmap=cmap,
                                               norm=matplotlib.colors.Normalize(-vmax * 1000, vmax * 1000))
    fig.colorbar(mappable_nm, cax=cax, label='WFE [nm]')

    # Draw a bar plot showing potential correction and green/yellow/red thresholds
    cax = fig.add_axes([0.90, 0.03, 0.05, 0.2])
    norm = matplotlib.colors.Normalize(0, 0.15)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('thresholds', ['LightGreen', 'Yellow', 'Salmon'], N=3)
    cb = fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=cax)
    cb.set_ticks([0.05, 0.10])
    cax.axvline(0.33, color='black', alpha=0.5)
    cax.axvline(0.66, color='black', alpha=0.5)
    cax.set_title(" Piston\nXtilt  \nYtilt  ", rotation=90, fontsize=14)
    ptts = np.abs(coeffs.reshape(18, 3))
    cax.plot(np.arange(18) / (64) + 0.02, ptts[:, 0], marker='+', ls='none', color='black')
    cax.plot(np.arange(18) / (64) + 0.34, ptts[:, 1], marker='+', ls='none', color='navy')
    cax.plot(np.arange(18) / (64) + 0.67, ptts[:, 2], marker='+', ls='none', color='black')
    cax.set_ylim(0, 0.15)
    fig.text(0.89, 0.30, "Potential\nCorrections:", fontsize=13, fontweight='bold',
             horizontalalignment='left')


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
    ptt_coeffs = poppy.zernike.decompose_opd_nonorthonormal_basis(delta_opd,
                                                                  aperture=(delta_opd != 0) & bigtiltmask,
                                                                  nterms=3)

    ptt_fit = poppy.zernike.compose_opd_from_basis(ptt_coeffs, poppy.zernike.zernike_basis,
                                                   npix=npix,
                                                   aperture=webbpsf.utils.get_pupil_mask(npix=npix))
    # Subtract that from the delta OPD
    return delta_opd - ptt_fit


def wavefront_drift_plots(opdtable, start_time, end_time, verbose=False,
                          vmax=100, n_per_row=9,
                          label_cid=False, label_visit=False):
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
        if verbose: print(*text)


    date_range_mask = ((opdtable['date_obs_mjd'].value > start_time.mjd) &
                       (opdtable['date_obs_mjd'].value < end_time.mjd))

    which_opds_mask = date_range_mask  # & (~ redundant_aps_mask) & ap_type_mask

    sum(which_opds_mask)

    ### Iterate over all relevant OPDs to retrieve the OPD, and calc the WFE RMS
    # Setup arrays to store results from the iteration:
    n = np.sum(which_opds_mask)
    rms_obs = np.zeros(n)
    rms_ote = np.zeros(n)
    wf_obs = np.zeros((n, 256, 256))
    wf_ote = np.zeros((n, 256, 256))

    opd, opdhdu = _read_opd(opdtable[which_opds_mask][0]['fileName'])
    mask = opd != 0

    # Get WSS Target Phase Map
    was_targ_file = os.path.join(webbpsf.utils.get_webbpsf_data_path(), 'NIRCam', 'OPD', 'wss_target_phase.fits')
    target_1024 = astropy.io.fits.getdata(was_targ_file)
    target_256 = poppy.utils.krebin(target_1024, (256, 256))

    wf_si = target_256 * mask

    dates = []
    last_visit = ''

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

    fig, axes = plt.subplots(figsize=(16, 10), nrows=nrows, ncols=n_per_row,
                             gridspec_kw={'hspace': 0.3, 'wspace': 0.01,
                                          'left': 0.01, 'right': 0.93, 'bottom': 0.01, 'top': 0.92})
    axes_f = axes.flat

    is_correction = np.zeros(n, bool)
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

        title = f"{date[0:10]}  {date[11:16]}\n$\Delta T =${deltat * 24:.1f} hr"
        if label_cid:
            title = f"{cid}\n" + title
        if label_visit:
            visit = row['visitId']
            title = f"{int(visit[1:6])}:{int(visit[6:9])}\n" + title

        show_opd_image(delta_opd * 1e6 * nanmask, mask=mask, ax=axes_f[i],
                       # deltatime_hrs=deltat*24,
                       vmax=vmax / 1e3, fontsize=labelfontsize,
                       )
        axes_f[i].set_title(title, fontsize=9)

        last_date_obs = row['date_obs_mjd']
        i += 1  # increment plot counter for next plot

    # Finishing details: hide unused axes, add colorbar and supertitle
    for j in range(i, len(axes_f)):
        axes_f[j].set_visible(False)
        pass

    cax = fig.add_axes([0.94, 0.26, 0.01, 0.4])
    mappable_nm = matplotlib.cm.ScalarMappable(cmap=matplotlib.cm.RdBu_r,
                                               norm=matplotlib.colors.Normalize(-vmax, vmax))
    fig.colorbar(mappable_nm, cax=cax, label='WFE [nm]')

    plt.suptitle(f"Telescope Wavefront Drifts: {start_time.isot[0:10]} to {end_time.isot[0:10]}\n",
                 fontsize=20, fontweight='bold')

    plt.savefig('wf_drifts.pdf')


##### Monthly Trending Plots, including OPDs, RMS WFE and PSF EE

def get_month_start_end(year, month):
    _, ndays = calendar.monthrange(year, month)
    start_date_str = f"{year:04d}-{month:02d}-01"
    end_date_str = f"{year:04d}-{month:02d}-{ndays:02d}"

    start_date = astropy.time.Time(start_date_str)
    end_date = astropy.time.Time(end_date_str)

    return start_date, end_date


def get_opdtable_for_month(year, mon):
    """Return table of OPD measurements for a given month.

    This includes the last measurement in the prior month too, so we can compute a delta
    to the first one
    """
    start_date, end_date = get_month_start_end(year, mon)

    # Start a little early, such that we are going to have at least 1 WFS before the start date
    pre_start_date = astropy.time.Time(start_date) - astropy.time.TimeDelta(4 * u.day)

    # Retrieve full OPD table, then trim to the selected time period
    opdtable0 = webbpsf.mast_wss.retrieve_mast_opd_table()
    opdtable0 = webbpsf.mast_wss.deduplicate_opd_table(opdtable0)
    opdtable = webbpsf.mast_wss.filter_opd_table(opdtable0, start_time=pre_start_date, end_time=end_date)

    # Trim the table to have 1 and only 1 precursor measurement -
    # we'll use this to compute the drift for the first WFS in the time period
    is_pre = [astropy.time.Time(row['date']) < start_date for row in opdtable]
    opdtable['is_pre'] = is_pre
    opdtable = opdtable[np.sum(is_pre) - 1:]

    return opdtable


def monthly_trending_plot(year, month, verbose=True, instrument='NIRCam', filter='F200W', vmax=200):
    """Make monthly trending plot showing OPDs, mirror moves, RMS WFE, and the resulting PSF EEs

    year, month : integers
        Self explanatory
    verbose : bool
        Print more verbose text output
    vmax : float
        Image display vmax for OPDs, given here in units of nanometers.
    """

    def vprint(*text):
        if verbose: print(*text)

    # Look up wavefront sensing and mirror move corrections for that month
    start_date, end_date = get_month_start_end(year, month)
    opdtable = get_opdtable_for_month(year, month)
    corrections_table = webbpsf.mast_wss.get_corrections(opdtable)

    inst = webbpsf.instrument(instrument)
    inst.filter = filter

    apmask = webbpsf.utils.get_pupil_mask()

    vprint(f"Computing PSFs for {len(opdtable)} OPD measurements in {year}-{month:02d}")
    vprint(f"Calculating for {instrument}, {filter}")

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

        vprint(f"Calculating PSF for {opd_fn}")
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

    sensing_tds = (astropy.time.Time(corrections_table['Post Move Sensing Time']) -
                   astropy.time.Time(corrections_table['Pre Move Sensing Time']))
    correction_times = astropy.time.Time(corrections_table['Pre Move Sensing Time']) + sensing_tds / 2

    npoints = sum(opdtable[1:]['wfs_measurement_type'] == 'pre')

    # Prepare for plotting

    cmap = matplotlib.cm.RdBu_r.copy()
    cmap.set_bad('0.4')

    def basic_show_image(image, ax, vmax=.3, nanmask=1):
        ax.imshow(image * nanmask, cmap=cmap, vmin=-vmax, vmax=vmax, origin='lower')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_yticks([])

    # Make the plots
    fig = plt.figure(constrained_layout=False, figsize=(16, 9), )

    subfigs = fig.subfigures(2, 1, hspace=0.02, height_ratios=[2, 2], )

    min_n_im_axes = 12
    im_axes = subfigs[0].subplots(3, max(npoints,min_n_im_axes),
                                  gridspec_kw={'wspace': 0.001, 'left': 0.06, 'right': 0.98,
                                               'bottom': 0.02,
                                               'top': 0.8})  # , gridspec_kw={'wspace':1.000, 'left': 0})

    axes = subfigs[1].subplots(2, 1, gridspec_kw={'left': 0.07, 'right': 0.98, 'top': 0.97})

    fs = 14  # Font size for axes labels

    fig.suptitle(f'WF Trending for {year}-{month:02d}', fontsize=fs * 1.5, fontweight='bold')

    #### Plot 1: Wavefront Error

    axes[0].plot_date(dates_array.plot_date, rms_obs * 1e9,
                      color='C1', ls='-', label='Observatory WFE at NIRCam NRCA3')
    axes[0].plot_date(dates_array.plot_date, rms_ote * 1e9,
                      color='C0', ls='-', label='Telescope WFE')
    for ax in axes:
        for corr_date in correction_times:
            ax.axvline(corr_date.plot_date, color='darkgreen', ls='--', alpha=0.5)

    #axes[0].axhline(59, ls=":", color='gray')
    #axes[0].axhline(80, ls=":", color='gray')
    axes[0].fill_between( [start_date.plot_date - 0.5, end_date.plot_date + 0.5],
                          [59,59], [80, 80], color='blue', alpha=0.08, label='Wavefront control target range')
    axes[0].legend(loc='lower right', fontsize=9)

    axes[0].set_ylim(0, 150)
    axes[0].set_ylabel("Wavefront Error\n[nm rms]", fontsize=fs, fontweight='bold')
    axes[0].set_xticklabels([])

    #### Plot 2: Encircled Energies

    ee_ax_ylim = 0.04
    ee_measurements = {}
    for i, ee_npix in enumerate([10, 2.5]):
        color = f'C{i*2+2}'
        ee_rad = inst.pixelscale * ee_npix
        ees_at_rad = np.asarray([e(ee_rad) for e in ees])
        ee_measurements[ee_npix] = ees_at_rad  # save for later

        median_ee = np.median(ees_at_rad)
        ee_ax_ylim = np.max([ee_ax_ylim, np.abs(ees_at_rad-median_ee).max()*1.1]) # display tweak: adjust the plot Y scale sensibly to its contents

        axes[1].plot_date(dates_array.plot_date, ees_at_rad-median_ee, ls='-', color=color,
                          label=f"$\Delta$EE within {ee_rad:.2f} arcsec ({ee_npix} pix)")

        axes[1].text(0.01, 0.75-i*0.12, f'Median EE within {ee_rad:.2f} arcsec = {median_ee:.3f}', color=color,
                     fontweight='bold',
                     transform=axes[1].transAxes)

    axes[1].fill_between( [start_date.plot_date - 0.5, end_date.plot_date + 0.5], -0.03, 0.03, color='gray', alpha=0.1, label="±3% change (stability requirement)")
    axes[1].set_xlabel("Date", fontsize=fs, fontweight='bold')
    axes[1].set_ylabel(f"Change in \nEncircled Energy\n{instrument} {filter}", fontsize=fs, fontweight='bold')
    axes[1].set_ylim(0.5, 1.0)
    axes[1].axhline(0, ls=":", color='gray')
    axes[1].set_ylim(-ee_ax_ylim, ee_ax_ylim)
    axes[1].legend(loc='upper right', fontsize=9)

    # Configure Axes for the time series plots
    for ax in axes[0:2]:
        ax.set_xlim(start_date.plot_date - 0.5,
                    end_date.plot_date + 0.5, )
        ax.xaxis.set_minor_locator(matplotlib.dates.DayLocator())
        ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(byweekday=matplotlib.dates.MONDAY,
                                                                   interval=1))

    nanmask = np.zeros_like(apmask) + np.nan
    nanmask[apmask == 1] = 1
    plot_index = -1
    from matplotlib.patches import ConnectionPatch

    for i, row in enumerate(opdtable):
        if verbose:
            print(row['fileName'], row['date'],
                  "Sensing" if row['wfs_measurement_type'] == 'pre' else "Post Mirror Move",
                  f"Obs WFE: {rms_obs[i] * 1e9:.1f} nm\tOTE WFE: {rms_ote[i]*1e9:.1f} nm")
        date = row['date']
        if astropy.time.Time(date) < start_date:
            continue

        delta_opd = wfes_ote[i] - wfes_ote[i - 1] if i > 0 else np.zeros((1024, 1024))

        vmax_micron = vmax / 1000
        if row['wfs_measurement_type'] == 'pre':
            plot_index += 1
            # Plot WFS in row 1
            basic_show_image(wfes_ote[i] * 1e6, ax=im_axes[0, plot_index], nanmask=nanmask,
                             vmax=vmax_micron)  # , title=None)
            rms_label = im_axes[0, plot_index].text(20, 20, f"{webbpsf.utils.rms(wfes_ote[i], mask=apmask)*1e9:.1f}", color='yellow', fontsize=fs*0.6)
            # Plot drift since last measurement in row 2
            basic_show_image(delta_opd * 1e6, ax=im_axes[1, plot_index], nanmask=nanmask, vmax=vmax_micron)
            im_axes[1, plot_index].text(20, 20, f"{webbpsf.utils.rms(delta_opd, mask=apmask)*1e9:.1f}", color='yellow', fontsize=fs*0.6)
            # mark row 3 as blank (default)
            basic_show_image(delta_opd + np.nan, ax=im_axes[2, plot_index], nanmask=nanmask, vmax=vmax_micron)
            im_axes[0, plot_index].set_title(row['date'][0:10])

            # Make a fancy connector line from the OPD plots to the correct time in the time series plots
            # Draw this for all sensing instances pre-move, but no need to repeat for any post-correction sensing
            cp = ConnectionPatch([0.5, 0], (dates_array[i].plot_date, 150),
                                 coordsA='axes fraction', coordsB='data',
                                 axesA=im_axes[2, plot_index], axesB=axes[0],
                                 color='darkgreen', ls='--', alpha=0.5
                                 )
            fig.add_artist(cp)

        else:
            # Update row 1 to show post-mirror-move WFS
            basic_show_image(wfes_ote[i] * 1e6, ax=im_axes[0, plot_index], nanmask=nanmask,
                             vmax=vmax_micron)  # , title=None)
            rms_label.set_visible(False)
            del rms_label # delete the previously-written one for the pre-move sensing
            im_axes[0, plot_index].text(20, 20, f"{webbpsf.utils.rms(wfes_ote[i], mask=apmask)*1e9:.1f}", color='yellow', fontsize=fs*0.6)
            # Plot correction in row 3
            basic_show_image(delta_opd * 1e6, ax=im_axes[2, plot_index], nanmask=nanmask, vmax=vmax_micron)
            im_axes[2, plot_index].text(20, 20, f"{webbpsf.utils.rms(delta_opd, mask=apmask)*1e9:.1f}", color='yellow', fontsize=fs*0.6)

    for i, l in enumerate(['WF Sensing', "Drifts", 'Corrections']):
        im_axes[i, 0].yaxis.set_visible(True)
        im_axes[i, 0].set_ylabel(l + "\n\n", fontsize=fs, fontweight='bold')

    for j in range(npoints, min_n_im_axes):
        for i in range(3):
            im_axes[i,j].set_visible(False)

    outname = f'wf_trending_{year}-{month:02d}.pdf'
    plt.savefig(outname, dpi=200)
    vprint(f"Saved to {outname}")

    wfs_type = [("Sensing" if row['wfs_measurement_type'] == 'pre' else "Post Mirror Move") for row in opdtable]
    result_table = astropy.table.QTable([opdtable['date'], opdtable['fileName'], wfs_type,
                                        rms_obs*1e9*u.nm, rms_ote*1e9*u.nm, ee_measurements[2.5], ee_measurements[10]],
                                       names=['Date', 'Filename', 'WFS Type', 'RMS WFE (OTE+SI)', 'RMS WFE (OTE only)',
                                              'EE(2.5 pix)', 'EE(10pix)'])
    return result_table
