import os
import numpy as np
import matplotlib.pyplot as plt
import webbpsf
import poppy
import astropy.table as table
import astropy.units as u
from webbpsf.utils import rms


### JWST Optical Budgets Information
# This module makes extensive use of information from the JWST Optical Budget
# by Paul Lightsey et al.
# See Lightsey's 'Guide to the Optical Budget'
#
# Note, WebbPSF does *not* parameterize the optical parameters in strictly the same way
# as the budget, intentionally. That's a statistical document, whereas this software attempts to
# model specific instances. We also combine some of the terms differently here when presenting,
# however the overall sums are precisely consistent.
#
# See Section 7 of the Guide to the Optical Budget which describes the breakdown into
#  wfe residuals after WFSC + image motion equivalent + stability

wfe_budget_info = None


def load_wfe_budget_info():
    """ Load optical budget info, if that hasn't been done already.
    (THis caches to avoid multiple relaods)
    """
    global wfe_budget_info
    if wfe_budget_info is None:
        # Load information from the formal optical budget
        wfe_budget_filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                           'otelm',
                                           'jwst_wfe_summary_from_optical_budget.csv')
        wfe_budget_info = table.Table.read(wfe_budget_filename, header_start=1)
    return wfe_budget_info


def get_budget_info(instname, param_name):
    """ Return required and predicted from optical budget for a given quantity

    The optical budget is read from disk in summary form, and cached for reuse on subsequent calls.

    """
    global wfe_budget_info
    load_wfe_budget_info()
    row = wfe_budget_info[wfe_budget_info['Value']==param_name]
    return row[f"{instname}_req"].data[0], row[f"{instname}_pred"].data[0]


def get_optical_budget_ref_wavelength(inst):
    """Reference wavelengths for Strehl calculations used in the JWST optical budget

    Parameters
    ----------
    inst : webbpsf.Instrument instance

    """
    # see optical budget page 'Syst' cell V5 formula
    ref = {'NIRCam': 2.0,
           'NIRSpec': 3.17,
           'MIRI': 5.6,
           'NIRISS': 2.4,
           'FGS': 2.47}
    return ref[inst.name]

@poppy.utils.quantity_input(wavelength=u.micron, rms_wfe=u.nanometer)
def rms2strehl(rms_wfe, wavelength):
    """ Compute Strehl ratio given RMS wavefront error and wavelength
    via the Marechal approximation

    Parameters
    ----------
    rms_wfe : float
        RMS WFE, in nanometers by default unless other astropy.Unit is specified
    wavelength
        Wavelength in micron by default unless other astropy.Unit is specified

    TODO units

    """
    # Lightsey guide to the optical budget eq 1
    return np.exp(- (2*np.pi*rms_wfe / wavelength)**2 )


@poppy.utils.quantity_input(wavelength=u.micron)
def imagemotion2strehl(rms_jitter_per_axis, wavelength):
    """Compute Strehl Ratio given RMS image motion per axis

    """
    # Lightsey guide to the optical budget eq 2

    Dtel = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER* u.m

    return 1./(1+ 0.25*(np.pi*Dtel/wavelength)**2 * ( 2* rms_jitter_per_axis**2))

@poppy.utils.quantity_input(rms_jitter_per_axis=u.milliarcsecond, wavelength=u.micron)
def imagemotion2equiv_wfe(rms_jitter_per_axis, wavelength):
    """ Compute equivalent WFE for a given RMS image motion per axis

    This is a rough approximation, of course: equivalent in terms of impact on Strehl
    but not at all an accurate model of the actual optical physics

    """
    # Lightsey guide to the optical budget eq 3
    # However note the as-written equation in the PDF is incorrect, and inconsistent with the formula
    # used in the actual spreadsheet. The version used here was back-engineered from the Excel file
    # and is precisely consistent with the value used there

    Dtel = 6.59 * u.m   # match precisely the value used in Lightsey's optical budget
    rms_jitter_radians = rms_jitter_per_axis.to_value(u.radian)

    result = (wavelength/(2*np.pi)) * np.sqrt(np.log(1.0 + 0.5 * (np.pi * rms_jitter_radians * Dtel/wavelength)**2 ))
    return result.to(u.nanometer)




def show_opd(opd, mask = None, ax=None, vmax=200, title=None, annotate_budget=None, instname=None,
            titlebold=False):
    """ Display an OPD; helper function used in visual error budget plots

    Shows the OPD, and annotates the RMS WFE within that aperture

    """
    if ax is None:
        plt.figure()
        ax = plt.gca()

    opd_nm = opd * 1e9
    if mask is not None:
        opd_nm[mask==0] = np.nan

    cm = plt.get_cmap(poppy.conf.cmap_diverging)
    cm.set_bad('0.75', alpha=0)

    ax.patch.set_facecolor('0.8')
    ax.imshow(opd_nm, vmin=-vmax, vmax=vmax, cmap = cm, origin='lower')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    if title is not None:
        ax.set_title(title,
            fontweight='bold' if titlebold else None)


    rmswfe = rms(opd_nm, mask)
    ax.text(0.04, 0.04, f"{rmswfe:.1f} nm rms ", transform=ax.transAxes)

    if annotate_budget:
        req, pred = get_budget_info(instname, annotate_budget)
        ax.text(0.96, 0.11, f"Req:  {req:3d}", transform=ax.transAxes,
           horizontalalignment='right', color='maroon')
        ax.text(0.96, 0.04, f"Pred:  {pred:3d}", transform=ax.transAxes,
           horizontalalignment='right', color='navy')

    return ax

def extract_si_wfe(inst, wave):
    """ Retrieve the SI WFE portion of an instrument's total WFE model
    """
    aberration = inst._get_aberrations()
    return aberration.get_opd(wave)


def get_dynamic_vibe(rms_nm=4):
    """ WAG to make an OPD for the vibe term"""
    ote = webbpsf.opds.OTE_Linear_Model_WSS()

    amp=rms_nm/1000
    for seg in ote.segnames[0:18]:
        ote.move_seg_local(seg, xtilt=amp, ytilt=amp, delay_update=True)

    ote.update_opd()
    return ote.opd.copy()

def visualize_wfe_budget(inst, slew_delta_time=14 * u.day, slew_case='EOL', ptt_only=False, verbose=True):
    """Display a visual WFE budget showing the various terms that sum into the overall WFE for a given instrument

    Compares a WebbPSF instrument instance with the JWST optical budget for that instrument

    Parameters
    ----------
    inst : webbpsf.JWInstrument
        A JWST instrument instance
    slew_delta_time : astropy.Quantity time
        Time duration for thermal slew model
    slew_case : basestring
        'BOL' or 'EOL' for beginning of life or end of life thermal slew model. EOL is about 3x higher amplitude
    ptt_only : bool
        When decomposing wavefront into controllable modes, use a PTT-only basis? The default is to use all
        controllable pose modes. (This is mostly a leftover debug option at this point, not likely useful in general)
    verbose : bool
        Be more verbose
    """

    vprint = lambda x: print(x) if verbose else None

    loc = f'{inst.detector} at {inst.detector_position}'

    vprint("generating optical models")
    osys = inst.get_optical_system()

    wavelen = get_optical_budget_ref_wavelength(inst) *u.micron
    wave = osys.input_wavefront(wavelen)

    #-------------------------------------------
    # Gather and compute the wavefront map terms
    ote = osys.planes[0]

    # Get entrance pupil aperture
    aperture = ote.get_transmission(wave)

    #--- OTE static

    vprint("inferring OTE static WFE terms")
    # Get OTE total WFE
    wfe_ote = ote.get_opd(wave).copy()

    # Figure out the field dependent part and factor that out
    if ote._include_nominal_field_dep:
        wfe_ote_field_dep_nominal = ote._get_field_dependence_nominal_ote(ote.v2v3)
    else:  # pragma: no cover
        wfe_ote_field_dep_nominal = np.zeros_like(wfe_ote)
    wfe_ote_field_dep_mimf = ote._get_field_dependence_secondary_mirror(ote.v2v3)
    wfe_ote_field_dep = wfe_ote_field_dep_nominal + wfe_ote_field_dep_mimf

    wfe_ote -= wfe_ote_field_dep  # pull the field dep part out of the total OTE WFE to get the control-point WFE

    # Decompose OTE control-point total WFE into different spatial frequency bins
    vprint(" decomposing WFE into controllable and uncontrollable spatial frequencies")
    if ptt_only:
        basis = webbpsf.opds.JWST_WAS_PTT_Basis()
        ndof = 3
    else:
        basis = webbpsf.opds.JWST_WAS_Full_Basis()
        ndof = 6
    seg_coeffs = poppy.zernike.opd_expand_segments(wfe_ote, aperture=aperture,  nterms=18*ndof,
                                  basis=basis )
    vprint(" modeling controllable and uncontrollable spatial frequencies")
    wfe_ote_controllable = poppy.zernike.opd_from_zernikes(seg_coeffs, basis=basis )
    wfe_ote_uncontrollable = wfe_ote - wfe_ote_controllable

    # A note on reference frames: We display the OTE OPD in entrance pupil orientation,
    # to match the standard convention adopted by the JWST WSS and other tools.
    # SI WFE and other terms are therefore flipped below into this orientation for display.
    wfe_ote_static = wfe_ote

    #--- OTE dynamic
    vprint("inferring OTE dynamic WFE terms")
    # get LOS jitter, and infer a defocus-like WF map representation based on that.
    # This isn't real WFE, of course, but is how the optical budgets represent this.
    # Note, Lightsey's optical budget rev Y "predicted" for this is 3.2 mas, which is optimistic...
    jitter = inst.options.get('jitter_sigma', 0.005) * u.arcsecond
    los_as_wfe = imagemotion2equiv_wfe(jitter, wavelen)
    vprint(f"los jitter {jitter}, as wfe {los_as_wfe}")
    wfe_for_imagemotion = poppy.zernike.opd_from_zernikes([0,0,0,los_as_wfe.to_value(u.meter)], aperture=aperture, npix=1024, outside=0)

    # estimate wavefront thermal drifts + frill + IEC
    wfe_ote_dynamic_thermal = ote._get_dynamic_opd(case=slew_case, delta_time=slew_delta_time)
    wfe_ote_dynamic_vibe = get_dynamic_vibe()


    wfe_ote_dynamic = wfe_ote_dynamic_thermal + wfe_ote_dynamic_vibe + wfe_for_imagemotion

    #--- SI and ISIM wfe
    vprint("inferring ISIM + SI WFE terms")
    # Get SI WFE
    # (this is already in SI exit pupil orientation --- therefore flip to match OTE pupil for display)
    wfe_si = extract_si_wfe(inst, wave)
    wfe_si = np.flipud(wfe_si)

    # very rough WAG model - treat all the budgeted SI dynamic WFE as defocus
    # (here we use 5.8 nm rms in a circular aperture to get the budgeted 5 nm in the tricontagon)
    wfe_si_dynamic = poppy.zernike.opd_from_zernikes([0,0,0,5.8e-9], npix=wfe_si.shape[0])

    # Sum OTE and SI WFE to total WFE in pupil

    wfe_total = wfe_ote_static + wfe_ote_dynamic + wfe_si


    #-----------------
    # Plot wavefronts
    #
    vprint("displaying plots")
    fig, axes = plt.subplots(figsize=(16,16), nrows = 4, ncols = 4,
                             gridspec_kw={'hspace': 0.25})

    show_opd(wfe_total, aperture, title=f'Total WFE for OTE+{loc}', ax=axes[0,0],
             annotate_budget='System Performance', instname=inst.name,
             titlebold=True)
             # See also the System Total including uncertainty reserve

    # OTE static
    show_opd(wfe_ote_static, aperture, title = f'OTE total static wavefront', ax=axes[1,0],
             annotate_budget='OTE total static', instname=inst.name,
             titlebold=True)

    show_opd(wfe_ote_controllable, aperture, title = f'OTE controllable mode residuals\n(low+mid s.f.)', ax=axes[1,1],
                          annotate_budget='OTE residual controllable modes (mid freq)', instname=inst.name)
    show_opd(wfe_ote_uncontrollable, aperture, title = f'OTE uncontrollable WFE\n(high s.f.)', ax=axes[1,2],
                          annotate_budget='OTE uncontrollable high freq', instname=inst.name)
    show_opd(wfe_ote_field_dep, aperture, title = f'OTE field-dependent WFE\n(low s.f.)', ax=axes[1,3],
                          annotate_budget='OTE residual low freq (field dep)', instname=inst.name)
    # OTE dynamic
    show_opd(wfe_ote_dynamic, aperture, title = f'OTE total dynamic wavefront', ax=axes[2,0],
             annotate_budget='OTE total dynamic', instname=inst.name,
             titlebold=True)
    ax = show_opd(wfe_ote_dynamic_thermal, aperture, title = f'OTE thermal drifts', ax=axes[2,1],
            annotate_budget='OTE stability', instname=inst.name)
    ax.text(0.04, 0.92, slew_case, transform=ax.transAxes)
    ax.text(0.96, 0.92, f"45$^\circ$, {slew_delta_time}", transform=ax.transAxes, horizontalalignment='right',)
    show_opd(wfe_ote_dynamic_vibe, aperture, title = f'OTE vibe', ax=axes[2,2],
            annotate_budget='OTE vibe', instname=inst.name)
    ax = show_opd(wfe_for_imagemotion, aperture, title = f'image motion$^*$ (as equiv. WFE)', ax=axes[2,3],
            annotate_budget='Image motion (as equiv. WFE)', instname=inst.name)
    ax.text(0.04, 0.92, f'LOS jitter: {jitter}, 1$\sigma$/axis', transform=ax.transAxes)


    # ISIM and SI
    show_opd(wfe_si, aperture, title = f'ISIM+SI total ', ax=axes[3,0],
             annotate_budget='ISIM+SI total', instname=inst.name,
             titlebold=True)
    show_opd(wfe_si, aperture, title = f'{inst.name} internal WFE at {inst.detector}, {inst.detector_position}',
             ax=axes[3,1],
             annotate_budget='SI internal WFE', instname=inst.name)
    show_opd(wfe_si*0, aperture, title = f'ISIM struct. align.', ax=axes[3,2],
             annotate_budget='ISIM structural', instname=inst.name)
    show_opd(wfe_si_dynamic, aperture, title = f'ISIM+SI instability', ax=axes[3,3],
             annotate_budget='ISIM+SI instability', instname=inst.name)



    for ax in [axes[0, 1], axes[0, 2], axes[0, 3], ]:
        ax.set_visible(False)
    tc = inst._tel_coords()


    # Annotate rows to show WFE summation
    def bounds(axis):
        return axis.get_position().get_points()

    for irow in [1,2,3]:
        between_col0_col1 = (bounds(axes[irow,1])[0,0] + bounds(axes[irow,0])[1,0])/2
        between_col1_col2 = (bounds(axes[irow,2])[0,0] + bounds(axes[irow,1])[1,0])/2
        between_col2_col3 = (bounds(axes[irow,3])[0,0] + bounds(axes[irow,2])[1,0])/2


        meany_row = bounds(axes[irow,0])[:,1].mean()
        fig.text(between_col0_col1, meany_row, "=", color='purple', horizontalalignment='center',
                 fontweight='bold', fontsize=18)
        fig.text(between_col1_col2, meany_row, "+", color='purple', horizontalalignment='center',
                 fontweight='bold', fontsize=18)
        fig.text(between_col2_col3, meany_row, "+", color='purple', horizontalalignment='center',
             fontweight='bold', fontsize=18)

    # Annotate first column to show WFE summation
    between_row1_row2 = (bounds(axes[1,0])[0,1] + bounds(axes[2,0])[1,1])/2
    between_row2_row3 = (bounds(axes[2,0])[0,1] + bounds(axes[3,0])[1,1])/2

    meanx_col0 = bounds(axes[0,0])[:,0].mean()
    for irow, symbol in enumerate(['=', '+', '+']):
        between_rows = (bounds(axes[irow, 0])[0, 1] + bounds(axes[irow+1, 0])[1, 1]) / 2
        fig.text(meanx_col0, between_rows, symbol, color='royalblue', horizontalalignment='center',
                 fontweight='bold', fontsize=18)


    meany_row0 = bounds(axes[0,0])[:,1].mean()
    fig.text(bounds(axes[irow,1])[0,0], meany_row0+0.05, f"Optical model for {inst.name} wavefront error",
        fontweight='bold', fontsize=15)

    fig.text(bounds(axes[irow,1])[0,0], meany_row0, f"OTE OPD model: {inst.pupilopd}\n\n"
             f"{inst.name}, detector {inst.detector} at {inst.detector_position}, aperture = {inst.aperturename}\n"
             f"(V2, V3):  ({tc[0].value:.3f}, {tc[1].value:.3f}) arcmin\n\n"
            "All WFE shown as projected to OTE entrance pupil orientation.",
            verticalalignment='center')




