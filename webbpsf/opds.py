###############################################################################
#
#           JWST Linear Optical Models for Mirror Motions
#
# This file implements two independent linear optical models for
# how the JWST OTE wavefront will change in response to commanded mirror
# motions. There are two models for historical reasons; both give fairly
# consistent answers within the uncertainties set by the inherent limits of
# linearizing a nonlinear optical problem.  One optical model is derived from
# the OTE influence functions in Elliott 2011, JWST-STScI-02356, "Sensitivity
# of wavefront error to segment motions in the JWST primary mirror."
# The second is derived from OTE influence functions consistent with those used
# in the JWST Wavefront Analysis Software by Ball Aerospace.
# All code and algorithms in this file are by Marshall Perrin
#
# This information is *not* considered sensitive under ITAR or EAR:
#   See GSFC Export Control Checklist "JWST Performance and Calibration
#   Technical Data Assessment (Jan. 2010)" which describes technical data
#   on the JWST program that is not considered ITAR controlled, in particular
#   items 1.1.2: "models and measures of wavefront error on both the telescope
#   and in the instrument" and item 4.7: "the effects of the OTE's optical
#   elements on a PSF."
#   See also the STScI Mission Office June 2015 determination that the
#   Elliott paper JWST-STScI-02356 is not ITAR sensitive, in part because
#   its contents are derivable from the publicly available OTE optical
#   prescription, as published for instance in Lightsey et al. 2012 Opt. Eng.
#
###############################################################################
import functools
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.special as sp
import scipy
import astropy.table
import astropy.io.fits as fits
import astropy.units as u
import logging
from collections import OrderedDict
import warnings
from packaging.version import Version
import copy

import poppy
import poppy.zernike as zernike

import pysiaf

from . import constants
from . import surs
from . import utils
import webbpsf

_log = logging.getLogger('webbpsf')

__doc__ = """

Linear Models for JWST segment motions

JWST_OTE_LOM_Elliot : OTE Linear optical model based on
    the influence functions in Elliott 2011, JWST-STScI-02356.

JWST_OTE_LOM_WSS: OTE linear optical model based on the
    influence functions used in the JWST WFSC Software System

These give very similar results, though not precisely consistent
because of small differences in the parameter ranges used when
deriving the linear approximations.

"""
__location__ = os.path.dirname(os.path.abspath(__file__)) + os.sep


################################################################################

class OPD(poppy.FITSOpticalElement):
    """ Base class for JWST Optical Path Difference (OPD) files

    Provides methods for analyzing and displaying OPDS.

    If you want to manipulate one using the linear optical model, see the OPDbender class.


    This class is implemented as a child of poppy.FITSOpticalElement; as such it
    follows all the same behavior as that, including having all OPDs internally
    in units of meters.

    Key class attributes include
    .OPD            the OPD data (as a fits.PrimaryHDU)
    .pupilmask      a mask for the pupil (as a fits.ImageHDU)
    .pixelscale     Scale in meters/pixel
    .header         reference to the OPD data's header


    """

    def __init__(self, name='unnamed OPD', opd=None, opd_index=0, transmission=None,
                 segment_mask_file=None, npix=1024, **kwargs):
        """
        Parameters
        ----------
        opd : string, path to FITS file.
            FITS file to load an OPD from. The OPD must be specified in microns. This FITS file must be
            compatible with the format expected by poppy.FITSOpticalElement.
        transmission: str
            FITS file for pupil mask, with throughput from 0-1. If not explicitly provided, will be inferred from
            wherever is nonzero in the OPD file
        npix: int
            Number of pixels per side for the OPD. Can be 1024, 2048, 4096 with current data files.
        segment_mask_file : string, or None
            if None, will infer based on npix.
        ext : int, optional
            FITS extension to load OPD from
        opd_index : int, optional
            slice of a datacube to load OPD from, if the selected extension contains a datacube.

        """

        self.npix = npix

        if opd is None and transmission is None:
            _log.debug('Neither a pupil mask nor OPD were specified. Using the default JWST pupil.')
            transmission = os.path.join(utils.get_webbpsf_data_path(), f"jwst_pupil_RevW_npix{self.npix}.fits.gz")

        super(OPD, self).__init__(name=name,
                                  opd=opd, transmission=transmission,
                                  opd_index=opd_index, transmission_index=0,
                                  planetype=poppy.poppy_core.PlaneType.pupil, **kwargs)

        # Check that the shape of the OPD that has been passed, matches the npix parameters
        if self.opd is not None:
            if not self.opd.shape == (self.npix, self.npix):
                raise ValueError(f"npix value of {self.npix} does not match shape of OPD provided: {self.opd.shape}")

        if self.opd_header is None:
            self.opd_header = self.amplitude_header.copy()

        self.segnames = np.asarray([a[0:2] for a in constants.SEGNAMES_WSS_ORDER])

        if segment_mask_file is None:
            segment_mask_file = f'JWpupil_segments_RevW_npix{self.npix}.fits.gz'

        full_seg_mask_file = os.path.join(utils.get_webbpsf_data_path(), segment_mask_file)
        if not os.path.exists(full_seg_mask_file):
            # try without .gz
            full_seg_mask_file = os.path.join(utils.get_webbpsf_data_path(), f"JWpupil_segments_RevW_npix{npix}.fits")

        self._segment_masks = fits.getdata(full_seg_mask_file)
        if not self._segment_masks.shape[0] == self.npix:
            raise ValueError(f"The shape of the segment mask file {self._segment_masks.shape} does not match the shape expect: ({self.npix}, {self.npix})")

        self._segment_masks_version = fits.getheader(full_seg_mask_file)['VERSION']

        # Where are the centers of each segment?  From OTE design geometry
        self._seg_centers_m = {seg[0:2]: np.asarray(cen)
                               for seg, cen in constants.JWST_PRIMARY_SEGMENT_CENTERS}
        # convert the center of each segment to pixels for the current array sampling:
        self._seg_centers_pixels = {seg[0:2]: self.shape[0] / 2 + np.asarray(cen) / self.pixelscale.value
                                    for seg, cen in constants.JWST_PRIMARY_SEGMENT_CENTERS}

        # And what are the angles of the local control coordinate systems?
        self._rotations = {}
        self._control_xaxis_rotations = {}
        self._control_xaxis_rot_base = {'A': 180,  # Rotations of local control coord
                                        'B': 0,  # X axes, CCW relative to V2 axis
                                        'C': 60}  # for A,B,C1

        for i in range(18):
            seg = self.segnames[i]
            self._rotations[seg] = (int(seg[1]) - 1) * 60  # seg_rotangles[i]
            self._control_xaxis_rotations[seg] = (self._control_xaxis_rot_base[seg[0]] +
                                                  -1 * self._rotations[seg])

        self._seg_tilt_angles = {'A': -4.7644,  # Tilts of local normal vector
                                 'B': 9.4210,  # relative to the V1 axis, around
                                 'C': -8.1919}  # the local X axis, or Y axis for Cs

        self.name = name
        self.header = self.opd_header  # convenience reference

    def copy(self):
        """ Make a copy of a wavefront object """
        from copy import deepcopy
        return deepcopy(self)

    # We can add and subtract OPD objects using the following.
    def __add__(self, x):
        """ Operator add """

        output = self.copy()
        if isinstance(x, OPD):
            output.opd += x.data
            output.name += " + " + x.name
        else:
            output.opd += x
            output.name += " + " + str(x)
        return output

    def __radd__(self, x):
        return self.__add__(x)

    def __sub__(self, x):
        """ Operator subtract """

        output = self.copy()
        if isinstance(x, OPD):
            output.opd -= x.opd
            output.name += " - " + x.name
        else:
            output.opd -= x
            output.name += " - " + str(x)
        return output

    def as_fits(self, include_pupil=True):
        """ Return the OPD as a fits.HDUList object

        Parameters
        -----------
        include_pupil : bool
            Include the pupil mask as a FITS extension?
        """

        output = fits.HDUList([fits.PrimaryHDU(self.opd, self.opd_header)])
        output[0].header['EXTNAME'] = 'OPD'
        output[0].header['BUNIT'] = 'meter'  # Rescaled to meters in poppy_core

        if include_pupil:
            self.amplitude_header['EXTNAME'] = 'PUPIL'
            output.append(fits.ImageHDU(self.amplitude, self.amplitude_header))

        return output

    def writeto(self, outname, overwrite=True, **kwargs):
        """ Write OPD to a FITS file on disk """
        self.as_fits(**kwargs).writeto(outname, overwrite=overwrite)

    # ---- display and analysis
    def powerspectrum(self, max_cycles=50, sampling=5, vmax=100, iterate=False):
        """
        Compute the spatial power spectrum of an aberrated wavefront.

        Produces nice plots on screen.

        Returns an array [low, mid, high] giving the RMS spatial frequencies in the
        different JWST-defined spatial frequency bins:
            low:   <=5 cycles/aperture
            mid:   5 < cycles <= 30
            high:  30 < cycles

        """

        # import SFT

        def tvcircle(radius=1, xcen=0, ycen=0, center=None, **kwargs):
            """
                draw a circle on an image.

                    radius
                    xcen
                    ycen
                    center=     tuple in (Y,X) order.
            """
            if center is not None:
                xcen = center[1]
                ycen = center[0]
            t = np.arange(0, np.pi * 2.0, 0.01)
            t = t.reshape((len(t), 1))
            x = radius * np.cos(t) + xcen
            y = radius * np.sin(t) + ycen
            plt.plot(x, y, **kwargs)

        cmap = copy.copy(matplotlib.cm.get_cmap(poppy.conf.cmap_diverging))
        cmap.set_bad('0.3')

        plt.clf()
        plt.subplot(231)
        self.display(title='full wavefront', clear=False, colorbar=False, vmax=vmax)

        ps_pixel_size = 1. / sampling  # how many cycles per pixel
        trans = SFT.SFT3(self.data, max_cycles * 2, max_cycles * 2 * sampling)

        abstrans = np.abs(trans)

        extent = [-max_cycles, max_cycles, -max_cycles, max_cycles]

        plt.subplot(233)
        plt.imshow(abstrans, extent=extent, origin='lower')
        plt.title("Power Spectrum of the phase")
        plt.ylabel("cycles/aperture")
        tvcircle(radius=5, color='k', linewidth=1)  # , ls='--')
        tvcircle(radius=30, color='k', linewidth=1)  # 2, ls='--')
        plt.gca().set_xbound(-max_cycles, max_cycles)
        plt.gca().set_ybound(-max_cycles, max_cycles)

        y, x = np.indices(abstrans.shape)
        y -= abstrans.shape[0] / 2.
        x -= abstrans.shape[1] / 2.
        r = np.sqrt(x ** 2 + y ** 2) * ps_pixel_size

        mask = np.ones_like(self.data)
        mask[np.where(self.amplitude == 0)] = np.nan
        wgood = np.where(self.amplitude != 0)

        components = []
        for i, label in enumerate(['low', 'mid', 'high']):
            plt.subplot(2, 3, i + 4)
            if label == 'low':
                condition = r <= 5
            elif label == 'mid':
                condition = (r > 5) & (r <= 30)
            else:
                condition = (r > 30)
            filtered = trans * condition

            inverse = SFT.SFT3(filtered, max_cycles * 2, self.opd.shape[0], inverse=True)
            inverse = inverse[::-1, ::-1]  # I thought SFT did this but apparently this is necessary to get the high freqs right...

            plt.imshow(inverse.real * mask, vmin=(-vmax) / 1000., vmax=vmax / 1000,
                       cmap=cmap, origin='lower')  # vmax is in nm, but WFE is in microns, so convert
            plt.title(label + " spatial frequencies")
            rms = (np.sqrt((inverse.real[wgood] ** 2).mean()) * 1000)

            components.append(rms)
            plt.xlabel("%.3f nm RMS WFE" % rms)

        return np.asarray(components)

    def display_opd(self, ax=None, labelsegs=True, vmax=150., colorbar=True, clear=False, title=None, unit='nm',
                    pupil_orientation='entrance_pupil',
                    cbpad=None, colorbar_orientation='vertical',
                    show_axes=False, show_rms=True, show_v2v3=False,
                    cmap=None):
        """ Draw on screen the perturbed OPD

        Parameters
        -----------
        ax : matplotlib.Axes
            axes instance to display into.
        labelsegs : bool
            draw segment name labels on each segment? default True.
        show_axes: bool
            Draw local control axes per each segment
        show_rms : bool
            Annotate the RMS wavefront value
        show_v2v3:
            Draw the observatory V2V3 coordinate axes
        pupil_orientation : string
            either 'entrance_pupil' or 'exit_pupil', for which orientation we should display the OPD in.
        clear : bool
            Clear plot window before display? default true
        unit : str
            Unit for WFE. default is 'nm'

        """

        if unit == 'nm':
            scalefact = 1e9
        elif unit == 'micron':
            scalefact = 1e6
        elif unit == 'm' or unit == 'meter':
            scalefact = 1
        else:
            raise ValueError("unknown unit keyword")

        # The actual OPD data is stored in "entrance pupil" orientation, looking in to JWST,
        # by convention and for consistency with the JWST WSS.
        # In the case of JWST it turns out the "exit pupil" orientation is in effect a vertical
        # flip. This is because of (1) inversion in both X and Y axes from the physical entrance pupil
        # at the primary to the real image of that pupil on the FSM which is the OTE exit pupil, plus
        # (2) the change in viewing convention from "in front of OTE looking in" to "at the instruments
        # looking outward". That results in a flip in the X axis. Thus overall we can just flip the
        # Y axis to get to the exit pupil orientation. In this display function we can do that using the
        # origin parameter to matplotlib.imshow. In the actual optical propagation we do that with a
        # poppy coordinate transform instance.
        if pupil_orientation=='entrance_pupil':
            origin='lower'
        elif pupil_orientation == 'exit_pupil':
            origin = 'upper'
        else:
            raise ValueError("pupil_orientation must be 'entrance_pupil' or 'exit_pupil'.")

        if clear:
            if ax is not None:
                raise RuntimeError("'clear=True' is incompatible with passing in an Axes instance.")
            plt.clf()
        if cmap is None:
            cmap = copy.copy(matplotlib.cm.get_cmap(poppy.conf.cmap_diverging))
        cmap.set_bad('0.3')

        mask = np.ones_like(self.opd)
        mask[np.where(self.amplitude == 0)] = np.nan

        try:
            pupilscale = self.opd_header['PUPLSCAL']
            s = self.opd.shape
            extent = [a * pupilscale for a in [-s[0] / 2, s[0] / 2, -s[1] / 2, s[1] / 2]]
        except KeyError:
            extent = None

        if ax is None:
            ax = plt.gca()

        plot = ax.imshow(self.opd * mask * scalefact, vmin=-vmax, vmax=vmax, cmap=cmap, extent=extent, origin=origin)

        _log.debug("Displaying OPD. Vmax is %f, data max is %f " % (vmax, self.opd.max()))

        if title is None:
            title = self.name
        ax.set_title(title)
        if show_rms:
            ax.set_xlabel("RMS WFE = %.1f nm" % self.rms())
        if show_v2v3:
            utils.annotate_ote_pupil_coords(None, ax, orientation=pupil_orientation)

        if labelsegs:
            for seg in self.segnames[0:18]:
                self.label_seg(seg, ax=ax, show_axes=show_axes, pupil_orientation=pupil_orientation)
        if colorbar:
            if cbpad is None:
                cbpad = 0.05 if colorbar_orientation == 'vertical' else 0.15
            cb = plt.colorbar(plot, ax=ax, pad=cbpad, orientation=colorbar_orientation)
            cb.set_label("WFE [%s]" % unit)
        else:
            cb = None
        plt.draw()

        return ax, cb

    def label_seg(self, segment, ax=None, show_axes=False, color='black', pupil_orientation='entrance_pupil'):
        """ Annotate a plot with a text label for a particular segment """
        cx, cy = self._seg_centers_m[segment]

        # See note about pupil orientations in OPD.display_opd() for explanation
        # of the Y axis signs here
        if pupil_orientation == 'entrance_pupil':
            ysign = 1
        elif pupil_orientation == 'exit_pupil':
            ysign = -1
        else:
            raise ValueError("pupil_orientation must be 'entrance_pupil' or 'exit_pupil'.")

        offset = 0.2 if show_axes else 0

        if ax is None:
            ax = plt.gca()
        label = ax.text(cx + offset, (cy + offset)*ysign, segment, color=color, horizontalalignment='center', verticalalignment='center')

        if show_axes:
            ax_arrow_len = .3
            # if not ('C' in segment ):
            if True:

                for i, color, label in zip([0, 1, 2], ['green', 'blue', 'red'], ['x', 'y', 'z']):
                    vec = np.matrix([0, 0, 0]).transpose()  # xyz order
                    vec[i] = 1
                    b = self._rot_matrix_local_to_global(segment) * vec
                    b = np.asarray(b).flatten()  # Inelegant but it works

                    ax.arrow(cx, cy*ysign, ax_arrow_len * b[0], (ax_arrow_len * b[1])*ysign, color=color,
                             # width=ax,
                             head_width=.050, head_length=.080)  # in units of mm

                    xoffset = 0.1 if i == 2 else 0
                    ax.text(cx + ax_arrow_len * b[0] * 1.5 + xoffset, (cy + ax_arrow_len * b[1] * 1.5)*ysign, label,
                            color=color, fontsize=8,
                            horizontalalignment='center', verticalalignment='center'
                            )

        ax.get_figure().canvas.draw()
        return label

    def _rot_matrix_local_to_global(self, segname):
        """ Rotation matrix from Local to Global coordinates

        Inverse of _rot_matrix_global_to_local
        """
        from . import geometry
        tilt = self._seg_tilt_angles[segname[0]]
        xaxis_rot = self._control_xaxis_rotations[segname]
        if 'C' in segname:  # Cs are tilted about Y, ABs tilted about X
            return geometry.rot_matrix_z(xaxis_rot) * geometry.rot_matrix_y(tilt)
        else:
            return geometry.rot_matrix_z(xaxis_rot) * geometry.rot_matrix_x(tilt)

    def _rot_matrix_global_to_local(self, segname):
        """ Rotation matrix from Global to Local coordinates
        """
        from . import geometry
        tilt = self._seg_tilt_angles[segname[0]]
        xaxis_rot = self._control_xaxis_rotations[segname]
        if 'C' in segname:  # Cs are tilted about Y, ABs tilted about X
            return geometry.rot_matrix_y(-tilt) * geometry.rot_matrix_z(-xaxis_rot)
        else:
            return geometry.rot_matrix_x(-tilt) * geometry.rot_matrix_z(-xaxis_rot)

    def zern_seg(self, segment, vmax=150, unit='nm'):
        """ Show the Zernike terms applied to a given segment"""

        # the actual wavefront is always in units of microns.
        if unit == 'nm':
            scalefact = 1000.
        elif unit == 'micron':
            scalefact = 1.0
        else:
            raise ValueError("unknown unit keyword")

        nzerns = 11
        title = "Zernike components for " + segment
        zerns = np.zeros(nzerns)
        if segment + "-decenter" in self.state.keys():
            zerns += self._move(segment, type='decenter', vector=self.state[segment + "-decenter"], return_zernikes=True)
            title += ", displacement=[%.2e, %.2e, %.2e] um" % tuple(self.state[segment + "-decenter"])
        if segment + "-tilt" in self.state.keys():
            zerns += self._move(segment, type='tilt', vector=self.state[segment + "-tilt"], return_zernikes=True)
            title += ", tilts =[%.5f, %.5f, %.5f] urad" % tuple(self.state[segment + "-tilt"])

        _log.info("Zerns: " + str(zerns))
        fig = plt.gcf()
        fig.clf()

        npix = 200
        hexap = zernike.hex_aperture(npix)
        hexap[np.where(hexap == 0)] = np.nan
        cmap = copy.copy(matplotlib.cm.get_cmap('jet'))
        cmap.set_bad('0.5', alpha=0.0)

        for j in np.arange(nzerns) + 1:
            ax = fig.add_subplot(3, 4, j, frameon=False, xticks=[], yticks=[])

            # n, m = zernike.noll_indices(j)
            Z = zernike.zernike1(j, npix=npix)
            ax.imshow(Z * zerns[j - 1] * hexap * scalefact, vmin=-1 * vmax, vmax=vmax, cmap=cmap, origin='lower')
            ax.text(npix * 0.95, npix * 0.8, "$Z{:d}$".format(j), fontsize=20, horizontalalignment='right')
            ax.text(npix * 0.95, npix * 0.1, "{:.2e}".format(zerns[j - 1]), fontsize=15, horizontalalignment='right')

        fig.text(0.5, 0.95, title, horizontalalignment='center', fontsize=15)
        fig.text(0.95, 0.15, 'Segment RMS WFE = {:.2f} {} '.format(self.rms(segment) * (scalefact / 1000), unit), fontsize=10,
                 horizontalalignment='right')
        # fig.text(0.95, 0.05, 'OPDs scaled from %.2e - %.2e um' % (-vmax, vmax), horizontalalignment='right', fontsize=10)
        fig.text(0.95, 0.05, 'OPDs scaled from {:.2f} - {:.2f} {}'.format(-vmax, vmax, unit), horizontalalignment='right', fontsize=10)

        plt.draw()

    def rms(self, segment=None):
        """ Return RMS WFE in nanometers

        Parameters
        ----------
        segment : string
            Segment name, to compute RMS for a single segment. Leave unspecified (None) to compute
            RMS WFE for the entire aperture. Segments are identified by name: A1-A6, B1-B6, C1-C6
        """
        if segment is None:
            # RMS for whole aperture
            wgood = np.where((self.amplitude != 0) & (np.isfinite(self.opd)))
        else:
            assert (segment in self.segnames)
            iseg = np.where(self.segnames == segment)[0][0] + 1  # segment index from 1 - 18
            wseg = np.where((self._segment_masks == iseg) & (np.isfinite(self.opd)))
            wgood = wseg

        return np.sqrt((self.opd[wgood] ** 2).mean()) * 1e9

    def ptv(self, segment=None):
        """return peak to valley WFE in nanometers

        Parameters
        ----------
        segment : string
            Segment name, to compute RMS for a single segment. Leave unspecified (None) to compute
            for the entire aperture.
        """
        if segment is None:
            # RMS for whole aperture
            wgood = np.where(self.amplitude != 0)
        else:
            assert (segment in self.segnames)
            iseg = np.where(self.segnames == segment)[0][0] + 1  # segment index from 1 - 18
            wseg = np.where(self._segment_masks == iseg)
            wgood = wseg

        peak = self.opd[wgood].max()
        valley = self.opd[wgood].min()
        return (peak - valley) * 1000

    def estimated_Strehl(self, wavelength, verbose=True):
        """ Compute an estimated Strehl given a wavelength in meters

        Parameters
        -----------
        wavelength : float
            in meters
        verbose : bool
            should I print out an informative message?
            """
        # rms = sqrt((-1.0)*alog(strehl))*wavelength/(2*!pi)

        rms = self.rms() * 1e-9
        strehl = np.exp(-(rms * 2 * np.pi / wavelength) ** 2)
        if verbose:
            print("For wavelength = {0} meters, estimated Strehl = {1} for {2} nm rms WFE".format(wavelength, strehl, rms * 1e9))
        return strehl


################################################################################

class OTE_Linear_Model_Elliott(OPD):
    """ Perturb an existing JWST OTE wavefront OPD file, by applying changes in WFE
    derived from Erin Elliott's linear optical model.  See Elliott 20122, JWST-STScI-02356

    Note that this model is strictly accurate only for a single NIRCam field
    point, but we can apply it elsewhere and it should not be too far off for small motions.

    """

    def __init__(self, opd=None, opd_index=0, transmission=None, rm_ptt=False, rm_piston=False, zero=False):
        """
        Parameters
        ----------
        opdfile : str or fits.HDUList
            FITS file to load an OPD from. The OPD must be specified in microns.
        ext : int, optional
            FITS extension to load OPD from
        opd_index : int, optional
            slice of a datacube to load OPD from, if the selected extension contains a datacube.
        pupilfile : str
            FITS file for pupil mask, with throughput from 0-1. If not explicitly provided, will be inferred from
            wherever is nonzero in the OPD file.
        rm_ptt : bool
            Remove piston, tip and tilt from all segments if set. Default is False.
        rm_piston : bool
            Remove piston only from all segments if set. Default is False
        zero : bool
            If set, reate an OPD which is precisely zero in all locations. Default is False.


        """

        OPD.__init__(self, opd=opd, opd_index=opd_index, transmission=transmission)

        self._sensitivities = astropy.table.Table.read(os.path.join(__location__, 'otelm', 'seg_sens.txt'), format='ascii', delimiter='\t')
        self.state = {}
        self.remove_piston_tip_tilt = rm_ptt
        self.remove_piston_only = rm_piston

        self._opd_original = self.opd.copy()
        if zero:
            self.zero()

    # ---- overall state manipulation

    def reset(self):
        """ Reset an OPD to the state it was loaded from disk.

        i.e. undo all segment moves.
        """
        self.opd = self._opd_original.copy()
        _log.info("Reset to unperturbed OPD")

    def zero(self):
        self.opd *= 0
        self.state = {}
        _log.info("Set OPD to zero WFE!")

    def load_state(self, newstate):
        _log.info("Loading new state.")
        self.zero()
        for k in newstate.keys():
            segname, motion = k.split('-')
            self._move(segname, motion, newstate[k])

    def print_state_v1(self):
        keys = self.state.keys()
        keys.sort()
        if keys is None:
            print("No perturbations")
        else:
            print("state = {")
            for k in keys:
                print("    '%s' : %s," % (k, repr(self.state[k])))
            print("    }")

    def print_state(self, type='local'):
        keys = self.state.keys()

        if type == 'local':  # control coords
            print(
                "Segment poses in local Control coordinates: "
                "(microns for decenter & piston, microradians for tilts, milliradians for clocking):")
            print("  \t %10s %10s %10s %10s %10s %10s" % ("X dec", "Y dec", "Z piston", "X tilt", "Y tilt", "Clocking"))
            for segment in self.segnames[0:18]:
                if segment + "-tilt" in keys:
                    tilts = self.state[segment + "-tilt"].tolist()
                else:
                    tilts = [0, 0, 0]
                if segment + "-decenter" in keys:
                    decenters = self.state[segment + "-decenter"].tolist()
                else:
                    decenters = [0, 0, 0]

                print("%2s\t %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" % tuple([segment] + decenters + tilts))

        elif type == 'Report':
            raise NotImplementedError("Coord conversions need work")
            print("Segment positions in Report coordinates: "
                  "(microns for decenter, microradians for alpha & beta, milliradians for gamma):")
            print("  \t %10s %10s %10s %10s %10s %10s" % ("X dec", "Y dec", "Z dec", "alpha", "beta", "gamma"))
            for segment in self.segnames:
                if segment + "-tilt" in keys:
                    tilts = self.state[segment + "-tilt"].tolist()
                else:
                    tilts = [0, 0, 0]
                if segment + "-decenter" in keys:
                    decenters = self.state[segment + "-decenter"].tolist()
                else:
                    decenters = [0, 0, 0]

                print("%2s\t %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" % tuple([segment] + decenters + tilts))
        elif type == 'Matlab':
            print("Segment positions in Matlab coordinates: (millimeters and degrees)")
            print("  \t %12s %12s %12s %12s %12s %12s" % ("X dec", "Y dec", "Z dec", "alpha", "beta", "gamma"))
            for segment in self.segnames:
                if segment + "-tilt" in keys:
                    tilts = self.state[segment + "-tilt"].tolist()
                else:
                    tilts = [0, 0, 0]
                if segment + "-decenter" in keys:
                    decenters = self.state[segment + "-decenter"].tolist()
                else:
                    decenters = [0, 0, 0]

                decenters = (np.array([decenters[1], decenters[0], -1 * decenters[2]]) / 1000).tolist()
                tilts[2] *= 1000  # convert millirad to microrad for consistency
                tilts = (np.array([tilts[1], tilts[0], -1 * tilts[2]]) * 1e-6 * 180 / np.pi).tolist()

                print("%2s\t %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e" % tuple([segment] + decenters + tilts))

    # ---- segment manipulation via linear model

    def _get_seg_sensitivities(self, segment='A1', type='decenter'):
        assert (segment in self.segnames)
        refseg = 2 if segment[0] == 'C' else 1
        fields = ['%s_%s_%s%d' % (axis, type, segment[0], refseg) for axis in ['X', 'Y', 'Z']]
        # divide by 1e6 since sensitivities in units of microns and we want meters

        return np.asarray([self._sensitivities[f] for f in fields]) / 1e6

    def _record(self, segment='A1', type='decenter', values=np.zeros(3)):
        """ update the state structure with the current segment positions """
        key = "%s-%s" % (segment, type)
        if key in self.state.keys():
            self.state[key] += values
        else:
            self.state[key] = values

    def print_(self):
        keys = self.state.keys()
        keys.sort()
        if keys is None:
            print("No perturbations")
        else:
            for k in keys:
                print("%s\t%s" % (k, str(self.state[k])))

    def _apply_zernikes_to_seg(self, segment, zernike_coeffs, coordsys='local', debug=False):
        """ Apply Zernike perturbations to a given segment

        Parameters
        ----------
        segment : string
            name of segment, A1 through C6
        zernike_coeffs : iterable of floats
            Zernike coefficients, in units of microns
        coordsys : string
            Coordinate system to apply the Zernikes in, either "global" or "local".
            Local is the BATC-defined "Control" coordinates for each segment.
        """
        assert (segment in self.segnames)

        iseg = np.where(self.segnames == segment)[0][0] + 1  # segment index from 1 - 18
        wseg = np.where(self._segment_masks == iseg)

        # determine the X and Y zernike coefficients for each segment
        # determine the center of each segment, as the mean of the min and max X and Y values
        # FIXME this should just use the BATC-provided center coordinates; see jwst_ote3d.py

        Y, X = np.indices(self.opd.shape, dtype=float)
        cx = np.mean([X[wseg].min(), X[wseg].max()])
        cy = np.mean([Y[wseg].min(), Y[wseg].max()])

        Y -= cy
        X -= cx

        # We have to compute the segment radius in pixels exactly here, to make sure we don't
        # have any pixels that are > that radius, based on the discrete pixelization and the
        # varying illuminated areas of some of the segments. FIXME this could be done better once
        # the segment centers are update as noted above. I think.
        # Hmm, and needs slight rounding up so add a ceil() call?
        R = np.sqrt(X[wseg] ** 2 + Y[wseg] ** 2)
        seg_radius = np.ceil(R.max())

        _log.debug("Segment {} is centered at pixel loc ({:.1f}, {:.1f}) with radius {:.1f} pix".format(segment, cx, cy, seg_radius))

        # Ah hell, this already does the rotations here??
        ang = -self._rotations[segment]  # This definitely has to be a positive sign,
        # to get the right X, Y for locally-defined zernikes
        # 2016-06: No that seems incorrect. Add negative sign;
        # aha it depends on which kind of segment! This gets things
        # working correctly for the A segments.
        if 'B' in segment:
            ang += 180  # As are opposite Bs
        elif 'C' in segment:
            ang += 60  # Cs are 60 deg back

        ang = np.deg2rad(ang)

        Xr = X * np.cos(ang) + Y * np.sin(ang)
        Yr = X * -np.sin(ang) + Y * np.cos(ang)

        # These are the BATC "Control" coordinates for each segment
        Xrc = Xr[wseg]
        Yrc = Yr[wseg]

        # But Erin Elliott's linear optical model was computed on a slightly
        # different coordinate system, which has a different sign convention,
        # and which depends on segment type

        if 'A' in segment:
            Xrc_elliott = -Xrc
            Yrc_elliott = Yrc
        else:
            Xrc_elliott = Xrc
            Yrc_elliott = -Yrc

        theta = np.arctan2(Yrc_elliott, Xrc_elliott)
        Rw = np.sqrt(
            Xrc_elliott ** 2 + Yrc_elliott ** 2) / seg_radius
        # This normalizes the Zernike input radial coord to maximum 1 on the segment aperture.

        if debug:
            plt.subplot(121)
            plt.imshow(Xr * (self._segment_masks == iseg), origin='lower')
            plt.title("Local X_Control for " + segment)
            plt.subplot(122)
            plt.imshow(Yr * (self._segment_masks == iseg), origin='lower')
            plt.title("Local Y_Control for" + segment)
            plt.draw()

        if self.remove_piston_tip_tilt:
            zernike_coeffs[0:3] = 0
        elif self.remove_piston_only:
            zernike_coeffs[0] = 0

        for i in range(len(zernike_coeffs)):
            zern = zernike.zernike1(i + 1, rho=Rw, theta=theta) * zernike_coeffs[i]
            self.opd[wseg] += zern

        outtxt = "Zs=[" + ", ".join(['%.1e' % z for z in zernike_coeffs]) + "]"
        _log.debug("     " + outtxt)

    def _move(self, segment, type='tilt', vector=None, display=False, return_zernikes=False):
        """
        Internal function that actually modifies the OPD by applying the linear optical model.
        Don't use this directly, use tilt() or displace() instead.

        The provided vector must be in **control** coordinates.

        """

        # print "Segment rotation for %s is %f" % (segment, self._rotations[segment])
        # ang = self._rotations[segment] * np.pi/180
        # local_coordX = vector[0] * np.cos(ang) + vector[1] * np.sin(ang)
        # local_coordY = vector[0] *-np.sin(ang) + vector[1] * np.cos(ang)
        # local_vector = np.array([local_coordX, local_coordY, vector[2]])
        local_vector = vector
        if type == 'tilt':
            local_vector[
                2] /= 1000  # convert Z tilt to milliradians instead of microradians because that is what the sensitivity tables use
            units = 'microradians for tip/tilt, milliradians for clocking'
        else:
            units = 'microns'

        sensitivities = self._get_seg_sensitivities(segment, type=type)

        zernike_coeffs_2d = sensitivities * local_vector[:, np.newaxis]
        zernike_coeffs = zernike_coeffs_2d.sum(axis=0)

        if return_zernikes:
            return zernike_coeffs
        else:
            _log.info("Segment %s requested %s: %s in %s" % (segment, type, str(vector), units))
            _log.debug("    local %s: %s" % (type, str(local_vector)))
            self._record(segment, type, vector)
            self._apply_zernikes_to_seg(segment, zernike_coeffs)

            if display:
                self.display()

    def tilt(self, segment, tiltX=0.0, tiltY=0.0, tiltZ=0.0, unit='urad', display=False, coordsys='elliott'):
        """ Tilt/rotate a segment some angle around X, Y, or Z.

        Parameters
        -----------
        segment : str
            Segment name, e.g. 'A1'
        tiltX, tiltY, tiltZ : floats
            Tilt angle, in microradians by default.
        unit : str
            Unit for displacements. Can be 'urad', 'radian', 'milliradian', 'arcsec', 'arcmin', 'milliarcsec'
        display : bool
            Display after moving?
        coordsys : string
            Name of coordinate system the input tilts are with respect to. Can be 'local', 'global', 'control', or 'elliott'

        """

        if coordsys == 'control':
            # flip signs to match Ball's control coordinates convention relative to Erin's LOM
            # don't use *= to avoid editing the value in the calling function
            if 'A' in segment:
                tiltX = -1 * tiltX
            elif 'B' in segment:
                tiltY = -1 * tiltY
            elif 'C' in segment:
                tiltX = -1 * tiltX
            tiltZ = -1 * tiltZ

        tilts = np.array([tiltX, tiltY, tiltZ]).astype(float)

        self.opd_header.add_history('Rotation: %s %s' % (str(tuple(tilts)), unit))

        # new sensitivity matrices are in urad for  alpha and beta, mrad for gamma.
        # first just convert all to urad.
        if unit.endswith('s'):
            unit = unit[:-1]
        unit = unit.lower()
        if unit == 'urad':
            pass
        elif unit == 'milliarcsec':
            tilts *= (1e6 * np.pi / (180. * 60 * 60 * 1000))
        elif unit == 'arcsec':
            tilts *= (1e6 * np.pi / (180. * 60 * 60))
        elif unit == 'arcmin':
            tilts *= (1e6 * np.pi / (180. * 60))
        elif unit == 'radian' or unit == 'rad':
            tilts *= 1e6
        elif unit == 'milliradian' or unit == 'mrad':
            tilts *= 1e3
        else:
            raise NotImplementedError('unknown unit')

        self._move(segment, 'tilt', tilts, display=display)

    def displace(self, segment, distX=0.0, distY=0.0, distZ=0.0, unit='micron', display=False, coordsys='elliott'):
        """ Move a segment some distance in X, Y, and Z.

        Parameters
        -----------
        segment : str
            Segment name, e.g. 'A1'
        distX, distY, distZ : floats
            Displacement distance, in microns by default.
        unit : str
            Unit for displacements. Can be 'micron', 'millimeter','nanometer', 'mm', 'nm', 'um'
        display : bool
            Display after moving?
        coordsys : string
            Name of coordinate system the input tilts are with respect to. Can be 'control', 'global' or 'elliott'
        """

        # Convert the provided coordinates into CONTROL coords

        if coordsys == 'elliott':
            pass  # no conversion needed
        elif coordsys == 'control':
            # flip signs of displacements as needed, due to differing conventions and
            # handedness. Depends on segment type.

            if 'A' in segment:
                # distY = -distY
                distZ = -distZ
            else:
                # distX = -distX
                distZ = -distZ
        elif coordsys == 'global':
            raise NotImplementedError('global not yet implemented')
        else:
            raise ValueError("invalid coordsys param")

        vector = np.array([distX, distY, distZ])

        self.opd_header.add_history('Displacement: %s %s' % (str(tuple(vector)), unit))
        if unit.endswith('s'):
            unit = unit[:-1]
        unit = unit.lower()
        if unit == 'micron' or unit == 'um':
            pass
        elif unit == 'millimeters' or unit == 'millimeter' or unit == 'mm':
            vector *= 1000
        elif unit == 'nm' or unit == 'nanometer' or unit == 'nanometers':
            vector /= 1000
        else:
            raise ValueError("Unknown unit for length: %s" % unit)

        self._move(segment, 'decenter', vector, display=display)

    def sinewave(self, cyclesX=0.0, cyclesY=0.0, amplitude=100., display=True):
        """ Add a (nonphysical) sine-wave phase of a given spatial frequency in X and Y, specified
            in terms of cycles per aperture.


            Those cycles for spatial frequency are defined in accordance with the JWST OTE spec,
            i.e. cycles per 6.605 m circumscribing circle, rather than the actual aperture shape.

            cyclesX, Y : float
                cycles / aperture
            amplitude : float
                amplitude in nm

        """

        if cyclesX == 0 and cyclesY == 0:
            _log.info("Must specify either X or Y cycles - doing nothing since neither was given")
            return

        ref_ap_diam = 6.605
        Y, X = np.indices(self.opd.shape, dtype=float)

        Y -= Y.mean()  # center
        X -= X.mean()  # center
        # Y *= 2*np.pi* cyclesY/(ref_ap_diam/self.pixelscale)
        # X *= 2*np.pi* cyclesX/(ref_ap_diam/self.pixelscale)

        wsegs = np.where(self._segment_masks != 0)

        # self.opd[wsegs] += (np.cos(X)[wsegs] + np.cos(Y)[wsegs]) * (amplitude/1000.)  # convert amplitude to microns

        self.opd[wsegs] += np.cos(2 * np.pi * (X * cyclesX + Y * cyclesY) / (ref_ap_diam / self.pixelscale))[wsegs] * (
                amplitude / 1000.)  # convert amplitude to microns

        # self.opd[wsegs] = 0.2 # debug

        _log.info("added sine wave: (%.2f, %.2f, %.2f)" % (cyclesX, cyclesY, amplitude))

        if display:
            self.display()

    def perturb_all(self, display=True, verbose=True, **kwargs):
        """ Randomly perturb all segments


        There is no good easy way to dial in a desired level of WFE right now.
        Best/easiest is to try iteratively adjusting the multiplier keyword and
        repeat until you get the desired level of RMS WFE.


        """
        for seg in self.segnames:
            self.perturb(seg, **kwargs)

        if display:
            plt.clf()
            self.display(colorbar=True, vmax=1)
        if verbose:
            print("")
            self.print_state(type='Report')
            print("")
            self.print_state(type='Matlab')

    def perturb(self, segment, max=False, multiplier=1.0):
        """ Randomly perturb a segment


        There is no good easy way to dial in a desired level of WFE right now.
        Best/easiest is to try iteratively adjusting the multiplier keyword and
        repeat until you get the desired level of RMS WFE.

        Parameters
        ----------
        segment : str
            segment name
        multiplier : float
            Scale factor for making larger or smaller motions.
        max : bool
            force the motion to the largest possible motion consistent with the linear optical
            model bounds.

        """
        bound = [13, 13, 0.2, 0.9, 0.9, 200]  # bounds on linear regime, the region over which the
        # linear model can be considered trustworthy.
        stepsize = [1, 1, 0.1, 0.1, 0.1, 1]  # what relative size motion should we put in the various terms?

        if max:
            steps = bound

        else:
            # stick in a random step between -1 and 1 times the desired step size
            steps = (np.random.random_sample(6) - 0.5) * 2 * np.array(stepsize) * multiplier

        self.displace(segment, steps[0], steps[1], steps[2])
        self.tilt(segment, steps[3], steps[4], steps[5])

        rms = self.rms(segment)
        print("After perturbation, segment %s has RMS WFE = %.1f nm" % (segment, rms))


class OTE_Linear_Model_WSS(OPD):
    """ Perturb an existing wavefront OPD file, by applying changes in WFE
    based on a linear optical model that is algorithmically consistent with the
    JWST wavefront analysis software.

    Note, the coordinate system conventions for this model are somewhat non-obvious,
    because they follow the local control coordinates for each hexapod mechanism.

    """

    def __init__(self, name='Unnamed OPD', opd=None, opd_index=0, transmission=None,
                 segment_mask_file=None, zero=False, rm_ptt=False,
                 rm_piston=False, v2v3=None, control_point_fieldpoint='nrca3_full',
                 npix=1024, include_nominal_field_dependence=True):
        """
        Parameters
        ----------
        opdfile : str or fits.HDUList
            FITS file to load an OPD from. The OPD must be specified in microns.
        opd_index : int, optional
            FITS extension to load OPD from
        transmission: str or fits.HDUList
            FITS file to load aperture transmission from.
        opd_index : int, optional
            slice of a datacube to load OPD from, if the selected extension contains a datacube.
        segment_mask_file : str
            FITS file for pupil mask, with throughput from 0-1. If not explicitly provided, will
            use JWpupil_segments_RevW_npix1024.fits, or equivalent for other value of npix
        zero: bool
            Load a perfectly zero OPD, overriding anything present in the opdfile parameter.
        rm_ptt : bool
            Remove piston, tip, and tilt terms, per segment? This is mostly for visualizing the higher order parts of
            the LOM.
        v2v3 : tuple of 2 astropy.Quantities
            Tuple giving V2,v3 coordinates as quantities, typically in arcminutes, or None to default to
            the master chief ray location between the two NIRCam modules.
        include_nominal_field_dependence : bool
            Include the Zernike polynomial model for OTE field dependence for the nominal OTE.
            Note, if OPD is None, then this will be ignored and the nominal field dependence will be disabled.
        control_point_fieldpoint: str
            A parameter used in the field dependence model for a misaligned secondary mirror.
            Name of the field point where the OTE MIMF control point is located, on instrument defined by "control_point_instr".
            Default: 'nrca3_full'.
            The OTE control point is the field point to which the OTE has been aligned and defines the field angles
            for the field-dependent SM pose aberrations.
        npix : int
            Size of OPD: npix x npix

        """

        OPD.__init__(self, name=name, opd=opd, opd_index=opd_index, transmission=transmission,
                     segment_mask_file=segment_mask_file, npix=npix)
        self.v2v3 = v2v3

        # load influence function table:
        self._influence_fns = astropy.table.Table.read(os.path.join(__location__, 'otelm', 'JWST_influence_functions_control_with_sm.fits'))

        # With updated sign convention in poppy 1.0.0, the WSS influence function values can be used in WebbPSF directly,
        # with no change in sign.
        if Version(poppy.__version__) < Version('1.0'):
            # For earlier poppy versions, fix IFM sign convention for consistency to WSS
            cnames = self._influence_fns.colnames
            for icol in cnames[3:]:
                self._influence_fns[icol] *= -1

        # WFTP10 hotfix for RoC sign inconsistency relative to everything else, due to outdated version of WAS IFM used in table construction.
        # FIXME update the IFM file on disk and then delete the next three lines
        roc_rows = self._influence_fns['control_mode']=='ROC'
        for icol in self._influence_fns.colnames[3:]:
            self._influence_fns[icol][roc_rows] *= -1

        self._control_modes = ['Xtilt', 'Ytilt', 'Piston', 'Clocking', 'Radial', 'ROC']
        self._sm_control_modes = ['Xtilt', 'Ytilt', 'Xtrans', 'Ytrans', 'Piston']
        # controllable modes in WAS order; yes it's not an obvious ordering but that's the order of the
        # WAS influence function matrix for historical reasons.
        self.state = {}
        self.segment_state = np.zeros((19, 6), dtype=float)  # 18 segs, 6 controllable DOF each, plus SM
        self.segnames = np.asarray(list(self.segnames) + ['SM'])  # this model, unlike the above, knows about the SM.

        self._opd_original = self.opd.copy()  # make a separate copy
        self.remove_piston_tip_tilt = rm_ptt
        self.remove_piston_only = rm_piston
        # Arbitrary additional perturbations can be added, ad hoc, as either Zernike or Hexike coefficients over the
        # whole primary. Use move_global_zernikes for a convenient interface to this.
        self._number_global_zernikes = 9
        self._global_zernike_coeffs = np.zeros(self._number_global_zernikes)
        self._global_hexike_coeffs = np.zeros(self._number_global_zernikes)

        # Field dependence model data
        # Note, if the OTE is set to None, we disable this automatically. This is to enable modeling the ideal case with
        # truly NO WFE for opd=None.
        self._include_nominal_field_dep = include_nominal_field_dependence if opd else False
        self._field_dep_file = None
        self._field_dep_hdr = None
        self._field_dep_data = None

        # Thermal OPD parameters
        self.delta_time = 0.0
        self.start_angle = 0.0
        self.end_angle = 0.0
        self.scaling = None
        self._thermal_model = OteThermalModel() # Initialize thermal model object

        self._thermal_wfe_amplitude = 0.0
        self._frill_wfe_amplitude = 0.0
        self._iec_wfe_amplitude = 0.0

        self.meta = OrderedDict()  # container for arbitrary extra metadata

        # DETERMINE INSTRUMENT BASED ON CONTROL FIELD POINT NAME:
        self.control_point_fieldpoint = control_point_fieldpoint
        control_point_detector = self.control_point_fieldpoint.split("_")[0]
        if 'nrc' in control_point_detector:
            control_point_instr = 'nircam'
        elif 'miri' in control_point_detector:
                control_point_instr = 'miri'
        elif 'nis' in control_point_detector:
                control_point_instr = 'niriss'
        elif 'nrs' in control_point_detector:
                control_point_instr = 'nirspec'

        self.ote_control_point = pysiaf.Siaf(control_point_instr)[self.control_point_fieldpoint.upper()].reference_point('tel')*u.arcsec
        
        if zero:
            self.zero()
        else:
            if self.v2v3 is not None:
                # Run an update immediately after initialization, to apply the field dependence model automatically
                self.update_opd()

    # ---- overall state manipulation

    def reset(self):
        """ Reset an OPD to the state it was loaded from disk.

        i.e. undo all segment moves.
        """
        self.opd = self._opd_original.copy()
        self.segment_state *= 0
        _log.info("Reset to unperturbed OPD")

    def zero(self, zero_original=False):
        """ Reset an OPD to precisely zero everywhere.

        Parameters
        ----------
        zero_original : bool
            should we zero out the stored copy of the original OPD?
            If so, then even using the reset() function won't set this
            back to the original value.

        """
        self.opd *= 0
        self.segment_state *= 0
        self._thermal_wfe_amplitude = 0
        self._frill_wfe_amplitude = 0
        self._iec_wfe_amplitude = 0
        if zero_original:
            self._opd_original *= 0
        self.name = "Null OPD"
        _log.info("Set OPD to zero WFE!")

    def print_state(self):
        keys = self.state.keys()

        print("Segment poses in Control coordinates: (microns for decenter & piston, microradians for tilts and clocking):")
        print("  \t %10s %10s %10s %10s %10s %10s" % tuple(self._control_modes))
        for i, segment in enumerate(self.segnames[0:18]):
            thatsegment = self.segment_state[i]

            print("%2s\t %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" % tuple([segment] + thatsegment.tolist()))
        if len(self.segnames) == 19:  # SM is present
            print("Secondary Mirror Pose in Control coordinates: ")
            print("  \t %10s %10s %10s %10s %10s     n/a" % tuple(self._sm_control_modes))
            segment = 18
            thatsegment = self.segment_state[18]
            print("%2s\t %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" % tuple([segment] + thatsegment.tolist()))

    # ---- segment manipulation via linear model

    def _get_seg_sensitivities(self, segment='A1'):
        assert (segment in self.segnames)

        # Get rows matching PMSA moves
        wseg = np.where((np.asarray(self._influence_fns['segment_moved'], dtype=str) == segment) &
                        (np.asarray(self._influence_fns['segment_affected'], dtype=str) == segment))
        table = self._influence_fns[wseg]

        coeffs = np.zeros((6, 9))
        # this code requires ordering is the same as expected (but verifies that)
        nhexike = 9

        assert len(table) == len(self._control_modes), 'Got wrong number of expected records from the table'

        for i, label in enumerate(self._control_modes):
            if table[i]['control_mode'] != self._control_modes[i]:
                raise RuntimeError("Influence function table has unexpected ordering")
            for h in range(nhexike):
                coeffs[i, h] = table[i]['Hexike_{}'.format(h)]

        # the coefficients are in the table natively in units of microns,
        # as preferred by most Ball code. WebbPSF works natively in meters
        # for wavefront, so we have to convert from microns to meters here:
        return coeffs * 1e-6

    def _get_seg_sensitivities_from_sm(self, segment='A1'):
        assert (segment in self.segnames)

        # Get rows matching PMSA moves
        wseg = np.where((np.asarray(self._influence_fns['segment_moved'], dtype=str) == 'SM') &
                        (np.asarray(self._influence_fns['segment_affected'], dtype=str) == segment))
        table = self._influence_fns[wseg]

        coeffs = np.zeros((5, 9))
        # this code requires ordering is the same as expected (but verifies that)
        nhexike = 9

        assert len(table) == len(self._sm_control_modes), 'Got wrong number of expected records from the table'

        for i, label in enumerate(self._sm_control_modes):
            if table[i]['control_mode'] != self._sm_control_modes[i]:
                raise RuntimeError("Influence function table has unexpected ordering")
            for h in range(nhexike):
                coeffs[i, h] = table[i]['Hexike_{}'.format(h)]

        # the coefficients are in the table natively in units of microns,
        # as preferred by most Ball code. WebbPSF works natively in meters
        # for wavefront, so we have to convert from microns to meters here:
        return coeffs * 1e-6

    def _apply_hexikes_to_seg(self, segment, hexike_coeffs, debug=False):
        """ Apply Hexike perturbations to a given segment, using the
        same hexike ordering as Ball does.

        Parameters
        ----------
        segment : string
            name of segment, A1 through C6
        hexike_coeffs : iterable of floats
            Zernike coefficients, in units of meters
        """
        assert (segment in self.segnames)

        iseg = np.where(self.segnames == segment)[0][0] + 1  # segment index from 1 - 18
        wseg = np.where(self._segment_masks == iseg)

        # Note that the influence function matrix values already take into
        # account the rotations between segment coordinate systems.
        # so here we can just work in unrotated coordinates.

        # determine the X and Y hexike coordinates for each segment
        # determine the center of each segment using the BATC-provided center coordinates; see constants.py

        Y, X = np.indices(self.opd.shape, dtype=float)
        cx, cy = self._seg_centers_pixels[segment]

        seg_radius = (X[wseg].max() - X[wseg].min()) / 2.0

        _log.debug("Segment %s is centered at pixel loc (%.1f, %.1f) with radius %.1f pix" % (segment, cx, cy, seg_radius))

        # These are the BATC "Control" coordinates for each segment
        Y = (Y - cy) / seg_radius
        X = (X - cx) / seg_radius
        # cut out just the values of interest inside this segment.
        Xc = X[wseg]
        Yc = Y[wseg]

        apmask = np.ones_like(Xc)  # by construction, we're only evaluating this for the good pixels

        hexikes = zernike.hexike_basis_wss(x=Xc, y=Yc, nterms=len(hexike_coeffs),
                                           aperture=apmask)
        # returns a list of hexike array values each with the same shape as Xc

        if self.remove_piston_tip_tilt:
            # Save the values of the PTT we are removing, for optional reference elsewhere
            self.meta[f'S{iseg:02d}PISTN'] = (hexike_coeffs[0], f"[m] Hexike piston coeff for segment {segment}")
            self.meta[f'S{iseg:02d}XTILT'] = (hexike_coeffs[1], f"[m] Hexike X tilt coeff for segment {segment}")
            self.meta[f'S{iseg:02d}YTILT'] = (hexike_coeffs[2], f"[m] Hexike Y tilt coeff for segment {segment}")

            hexike_coeffs[0:3] = 0
        elif self.remove_piston_only:
            self.meta[f'S{iseg:02d}PISTN'] = (hexike_coeffs[0], f"[m] Hexike piston coeff for segment {segment}")
            hexike_coeffs[0] = 0
        else:
            # If remove_piston_tip_tilt is off, we shouldn't output those extra keywords, so
            # delete any prior saved values for those.
            try:
                del self.meta[f'S{iseg:02d}PISTN']
                del self.meta[f'S{iseg:02d}XTILT']
                del self.meta[f'S{iseg:02d}YTILT']
            except:
                pass

        for i in range(len(hexike_coeffs)):
            if hexike_coeffs[i] == 0:
                continue
            self.opd[wseg] += hexikes[i] * hexike_coeffs[i]


    def _apply_global_zernikes(self):
        """ Apply Zernike perturbations to the whole primary

        """

        perturbation = poppy.zernike.opd_from_zernikes(self._global_zernike_coeffs,
                                                       npix=self.npix,
                                                       outside=0,
                                                       basis=poppy.zernike.zernike_basis_faster)
        # Add perturbation to the opd
        self.opd += perturbation

    def _apply_global_hexikes(self, coefficients=None):
        """ Apply Hexike perturbations to the whole primary

        Parameters
        ----------
        coefficients : ndarray. By default this applies the self._global_hexike_coeffs values, but
            you can override that by optionally providing different coefficients to this function.
            In particular this is used behind the scenes for adding the thermal slew global terms 
            onto any other global hexikes already present. See update_opd()
        """

        if coefficients is None:
            coefficients = self._global_hexike_coeffs

        def _get_basis(*args, **kwargs):
            """ Convenience function to make basis callable """
            return basis
        # Define aperture as the full OTE
        aperture = self._segment_masks != 0
        basis = poppy.zernike.hexike_basis_wss(nterms=self._number_global_zernikes, npix=self.npix, aperture=aperture > 0.)
        # Use the Hexike basis to reconstruct the global terms
        perturbation = poppy.zernike.opd_from_zernikes(coefficients,
                                                       basis=_get_basis,
                                                       aperture=aperture > 0.,
                                                       outside=0)
        perturbation[~np.isfinite(perturbation)] = 0.0
        # Add perturbation to the opd
        self.opd += perturbation

    #---- OTE field dependence is implemented across the next several functions ----
    def _apply_field_dependence_model(self, reference='global', **kwargs):
        """Apply field dependence model for OTE wavefront error spatial variation.

        Includes SM field-dependent model for OTE wavefront error as a function of field angle and SM pose amplitude,
        and model for field-dependence of a nominal (perfectly aligned) OTE.

        Updates self.opd based on V2V3 coordinates and amplitude of SM pose.

        """
        if self.v2v3 is None:
            # if this OTE linear model instance does not have field coordinates set,
            # we cannot apply field dependence.
            return

        # Model field dependence for the nominal OTE, based on Code V model coefficients
        if self._include_nominal_field_dep:
            field_dep_nominal = self._get_field_dependence_nominal_ote(self.v2v3, reference=reference, **kwargs)
            self.opd += field_dep_nominal
            self.opd_header['HISTORY'] = 'Applied OTE nominal field-dependent aberrations'

        # Model field dependence from a misaligned secondary mirror
        # (only do so if the SM is in fact misaligned)
        if not np.allclose(self.segment_state[-1], 0):
            field_dep_sm_perturbation = self._get_field_dependence_secondary_mirror(self.v2v3)
            self.opd += field_dep_sm_perturbation
            self.opd_header['HISTORY'] = 'Applied OTE SM alignment field-dependent aberrations (MIMF)'


    def _get_hexike_coeffs_from_smif(self, dx=0., dy=0.):
        '''     
        Apply Secondary Mirror field dependence based on SM pose and field angle,
        using coefficients from the Secondary Mirror Influence Functions

        NOTES:
        Wavefront Phi(x,y) = Sum[ m_j * Sum[ (alpha_jk*dx + beta_jk*dy)*Z_k ] ]
        where m_j is the SM mode error (0: X-trans, 1: Y-trans, 2: X-tilt, 3: Y-tilt);
        dx/dy are the field separation from the control point, in radians;
        Z_k are the Zernike (3: defocus, 4: 0-Degree Astig, 5: 45-Degree Astig);
        alpha and beta are taken from the SMIF_Matrix/hexike.csv calculated 
        by Randal Telfer.
        
        See Also Ball SER 2288152 "JWST WFSC MIMF Control Algorithm Design".

        Parameters
        ----------
        dx, dy = field angles, in radians, from the OTE control point (nominally, NRCA3_FULL).

        Returns:
        ----------
        Hexike coefficients, in microns, for the first nine (9) coefficients, with PTT explicitly set to zero.

        '''
        
        # GET SM POSE:
        # RE-ORDER SUCH THAT (1) X-TRANS, (2) Y-TRANS, (3) X-TILT, (4) Y-TILT.
        sm_errors = [self.segment_state[-1, 2],    # X-TRANS
                     self.segment_state[-1, 3],    # Y-TRANS
                     self.segment_state[-1, 0],    # X-TILT
                     self.segment_state[-1, 1],    # Y-TILT
                     self.segment_state[-1, 4]]    # PISTON

        # GET SM INFLUENCE MATRIX:
        # HEXIKE PROJECTED ONTO ENTRANCE PUPIL (WHAT WEBBPSF NEEDS):
        smif = astropy.table.Table.read(os.path.join(__location__, 'otelm', "SMIF_hexike.csv"), header_start=5)

        alphas = smif[smif['Type'] == 'alpha']
        betas  = smif[smif['Type'] == 'beta']

        # NOW, WE SUM UP THE ALPHAS AND BETAS FOR EACH HEXIKE,
        # AS DESCRIBED IN THE DOC STRING ABOVE.
        coeffs = np.zeros((6))
        for i, icol in enumerate(smif.colnames[2:]):
            coeffs[i] = np.sum( (alphas[icol]*dx + betas[icol]*dy )*sm_errors)

        # PAD IN ZEROS FOR PISTON, TIP, TILT:
        coeffs = np.insert(coeffs, 0, [0., 0., 0.])
    
        return coeffs

    def _get_field_dependence_secondary_mirror(self, v2v3):
        """Calculate field dependence model for OTE with misaligned secondary mirror,
        as is to be corrected by MIMF (Multi Instrument Multi Field) sensing.

        Returns OPD based on V2V3 coordinates using model for secondary mirror influence functions.

        Parameters
        ----------
        v2v3 : tuple
            JWST focal plane coordinates (V2,V3) as a tuple of astropy Quantities with angular dimension
        """
        # Model field dependence from any misalignment of the secondary mirror
        dx = -(v2v3[0] - self.ote_control_point[0]).to( u.rad).value
        # NEGATIVE SIGN IN THE ABOVE B/C TELFER'S FIELD ANGLE COORD. SYSTEM IS (X,Y) = (-V2,V3)
        dy = (v2v3[1] - self.ote_control_point[1]).to(u.rad).value
        z_coeffs = self._get_hexike_coeffs_from_smif(dx, dy)

        if np.any(z_coeffs != 0):
            perturbation = poppy.zernike.opd_from_zernikes(z_coeffs, npix=self.npix,
                                                           basis=poppy.zernike.hexike_basis_wss, aperture=self.amplitude,
                                                           outside=0)
        else:
            # shortcut for the typical case where the SM is well aligned
            perturbation = np.zeros((self.npix, self.npix), float)

        if Version(poppy.__version__) < Version('1.0'):
            wfe_sign = -1  # In earlier poppy versions, fix sign convention for consistency with WSS
        else:
            wfe_sign = 1

        for i in range(3,9):
            self.opd_header[f'SMIF_H{i}'] = (z_coeffs[i], f"Hexike coeff from S.M. influence fn model")
        self.opd_header['HISTORY'] = ('Field point (x,y): ({})'.format(v2v3))
        self.opd_header['HISTORY'] = (
            'Control point: {} ({})'.format(self.control_point_fieldpoint.upper(), self.ote_control_point))
        self.opd_header['HISTORY'] = ('Delta x/y: {} {} (radians))'.format(dx, dy))

        return perturbation * (1e-6 * wfe_sign)


    def _get_field_dependence_nominal_ote(self, v2v3, reference='global',
                                          zern_num=78, legendre_num=None, assume_si_focus=True):
        """Calculate field dependence model for OTE nominal wavefront error spatial variation,

        Returns OPD based on V2V3 coordinates using model(s) for spatial variations.

        zern_num and legendre_num allow for fewer terms to be used in calculating the wavefront at a field angle
        than are specified in the calibration files.  This is to reduce required computation if the increase in
        accuracy isn't needed.

        The precomputed data files support up to 136 Zernikes, but using higher numbers provides diminishing
        returns for the increased computation time. Tests have shown that using 78 Zernikes yields within
        7-10 nm of the result of the full model, and 36 Zernikes yields within 15 nm.

        Parameters
        ----------
        v2v3 : tuple
           (v2,v3) coordinates
        reference : str
           should always be "global" in the current implementation.
        zern_num : int
            Number of Zernikes to use to represent the WFE at each field point
        legendre_num : int
        assume_si_focus : bool
            Assume that the SIs have been independently refocused, i.e. take out the net average
            defocus per each SI from the overall focal plane curvature
        """
        # The code for this is now split up into several sub-functions for improved clarity.

        if v2v3 is None:  # pragma: no cover
            return 0

        instrument = utils.determine_inst_name_from_v2v3(v2v3)
        if not self._load_ote_field_dep_data(instrument, reference):  # pragma: no cover
            return 0

        field_coeff_order, f_ang_unit, opd_to_meters, zern_num, legendre_num = self._validate_ote_field_dep_data(instrument,
                zern_num=zern_num, legendre_num=legendre_num)

        zernike_coeffs = self._get_zernikes_for_ote_field_dep(v2v3, instrument, field_coeff_order, f_ang_unit,
                zern_num, legendre_num, assume_si_focus=assume_si_focus)

        # Apply perturbation to OPD according to Zernike coefficients calculated above.
        perturbation = poppy.zernike.opd_from_zernikes(zernike_coeffs * opd_to_meters,
                                                       npix=self.npix,
                                                       basis=poppy.zernike.zernike_basis_faster,
                                                       outside=0)

        return perturbation

    def _load_ote_field_dep_data(self, instrument, reference):
        """Load tables of Zernikes vs field position

        Returns True if successful (file loaded OK, or was already loaded), False for failure.
        """

        base_path = utils.get_webbpsf_data_path()
        field_dep_file = os.path.join(base_path, f'{instrument}/OPD/field_dep_table_{instrument.lower()}.fits')

        # For efficiency, load from disk only if needed. And, for back-compatibility, fail gracefully if file not found
        if self._field_dep_file != field_dep_file:
            self._field_dep_file = field_dep_file
            _log.info(f'Loading field dependent model parameters from {self._field_dep_file}')

            try:
                self._field_dep_hdr = fits.getheader(field_dep_file)
                # Read in the data table with the coefficients for our model
                #   hdu[1] ==> local reference point
                #   hdu[2] ==> # global reference point
                if reference == 'global':
                    ext = 2
                elif reference == 'local':  # pragma: no cover
                    ext = 1
                    # we don't need both options for this. Simpler to only support one (even though the data files have two)
                    raise ValueError("Local field dependent OTE coordinates discouraged; use global")
                else:  # pragma: no cover
                    raise ValueError('Invalid wavefront reference')
                self._field_dep_data = fits.getdata(field_dep_file, ext=ext)

            except FileNotFoundError:  # pragma: no cover
                warnings.warn(f"Could not load {self._field_dep_file}; OTE field dependence model disabled")
                return False
        return True


    def _validate_ote_field_dep_data(self, instrument, zern_num=None, legendre_num=None):
        """Sanity check inputs after loading OTE field dep data. 

        Returns several parameter values used subsequently in the OTE OPD calculation
        """
        # Pull useful parameters from header
        hdr = self._field_dep_hdr

        # Check to make sure that we've got the right instrument
        if hdr['instr'].lower() != instrument.lower():  # pragma: no cover
            ValueError('Instrument inconsistent with field dependence file')

        # Check to make sure that the file has the right type of data in it and throw exception if not
        if hdr['wfbasis'] != 'Noll Zernikes':  # pragma: no cover
            raise ValueError('Data file contains data with unsupported wavefront polynomial expansion')

        if hdr['fiebasis'] != 'Legendre Polynomials':  # pragma: no cover
            raise ValueError('Data file contains data with unsupported field polynomial expansion')

        num_wavefront_coeffs = hdr['ncoefwf']
        if zern_num is None:
            zern_num = num_wavefront_coeffs
        if (zern_num > num_wavefront_coeffs):  # pragma: no cover
            raise ValueError('Data file contains fewer wavefront coefficients than specified')

        num_field_coeffs = hdr['ncoeffie']
        if legendre_num is None:
            legendre_num = num_field_coeffs
        if (legendre_num > num_field_coeffs):  # pragma: no cover
            raise ValueError('Data file contains fewer field coefficients than specified')

        field_coeff_order = hdr['fieorder']
        # If we aren't using all the Legendre terms in the table we can update the order of
        # Legendre polynomial required so we don't have to calculate all of them.
        if legendre_num is not None:
            field_coeff_order = int(np.ceil((-3 + np.sqrt(1 + 8 * legendre_num)) / 2))

        # Check the field angle units in the input file
        if hdr['fangunit'] == 'degrees':
            f_ang_unit = u.deg
        elif hdr['fangunit'] == 'arcmin':
            f_ang_unit = u.arcmin
        elif hdr['fangunit'] == 'arscec':
            f_ang_unit = u.arcsec
        else:  # pragma: no cover
            raise ValueError('Field angle unit specified in file is not supported')

        # Check the OPD units in the input file
        if hdr['opdunit'] == 'nm':
            opd_to_meters = 1e-9
        elif hdr['opdunit'] == 'um':
            opd_to_meters = 1e-6
        elif hdr['opdunit'] == 'm':
            opd_to_meters = 1
        elif hdr['opdunit'] == 'pm':
            opd_to_meters = 1e-12
        else:  # pragma: no cover
            ValueError('OPD unit specified in file is not supported')

        return field_coeff_order, f_ang_unit, opd_to_meters, zern_num, legendre_num

    def _get_zernikes_for_ote_field_dep(self, v2v3, instrument, field_coeff_order, f_ang_unit, 
            zern_num, legendre_num , assume_si_focus=True):
        """ Calculate the Zernike coeffs from the lookup table data

        This is the main numerical piece of the algorithm.
        """
        hdr = self._field_dep_hdr
        # Extent of box in field over which the data is defined
        min_x_field = hdr['MINXFIE'] * u.arcmin
        max_x_field = hdr['MAXXFIE'] * u.arcmin
        min_y_field = hdr['MINYFIE'] * u.arcmin
        max_y_field = hdr['MAXYFIE'] * u.arcmin

        # Calculate field angle for our model from the V2/V3 coordinates
        x_field_pt = hdr['v2sign'] * (v2v3[0] - hdr['v2origin'] * f_ang_unit)
        y_field_pt = hdr['v3sign'] * (v2v3[1] - hdr['v3origin'] * f_ang_unit)

        _log.info(f'Calculating field-dependent OTE OPD at v2 = {v2v3[0]:.3f}, v3 = {v2v3[1]:.3f}')
        _log.debug(f'Calculating field-dependent OTE OPD at CodeV X field = {x_field_pt:.3f}, Y field= {y_field_pt:.3f}')

        # Confirm that the calculated field point is within our model's range
        if ((x_field_pt < min_x_field) or (x_field_pt > max_x_field) or
                (y_field_pt < min_y_field) or (y_field_pt > max_y_field)):
            # If not within the valid region, find closest point that is
            x_field_pt0 = x_field_pt*1
            y_field_pt0 = y_field_pt*1
            x_field_pt = np.clip(x_field_pt0, min_x_field, max_x_field)
            y_field_pt = np.clip(y_field_pt0, min_y_field, max_y_field)

            clip_dist = np.sqrt((x_field_pt-x_field_pt0)**2 + (y_field_pt-y_field_pt0)**2)
            if clip_dist > 0.1*u.arcsec:
                # warn the user we're making an adjustment here (but no need to do so if the distance is trivially small)
                warnings.warn(f'For (V2,V3) = {v2v3}, Field point {x_field_pt}, {y_field_pt} not within valid region for field dependence model of OTE WFE for {instrument}: {min_x_field}-{max_x_field}, {min_y_field}-{max_y_field}. Clipping to closest available valid location, {clip_dist} away from the requested coordinates.')

        # Get value of Legendre Polynomials at desired field point.  Need to implement model in G. Brady's prototype
        # polynomial basis code, independent of that code for now.  Perhaps at some point in the future this model
        # can become more tightly coupled with WebbPSF/Poppy and we just call it here instead.
        # Calculate value of Legendre at all orders at our field point of interest.

        field_center_x = (max_x_field + min_x_field) / 2
        field_center_y = (max_y_field + min_y_field) / 2

        x_field_pt_norm = float((x_field_pt - field_center_x) / ((max_x_field - min_x_field) / 2))
        y_field_pt_norm = float((y_field_pt - field_center_y) / ((max_y_field - min_y_field) / 2))

        _log.debug(f'{instrument} max_x={max_x_field} min_x={min_x_field} max_y={max_y_field} min_y={min_y_field}')
        _log.debug(f'Normalized field point {x_field_pt_norm}, {y_field_pt_norm}')
        poly_x1d = np.zeros(field_coeff_order + 1)
        poly_y1d = np.zeros(field_coeff_order + 1)
        for index in range(0, field_coeff_order + 1):
            leg_poly1d = sp.legendre(index)
            poly_y1d[index] = leg_poly1d(y_field_pt_norm)
            poly_x1d[index] = leg_poly1d(x_field_pt_norm)

        #Calculate product of x and y Legendre value for all combinations of orders
        poly_val_2d = np.einsum('i,j', poly_y1d, poly_x1d)

        #Reorder and rearrange values to correspond to the single index ordering the input coefficients assume
        map1 = []
        map2 = []
        for index_i in range(0, field_coeff_order + 1):
            map1 += list(range(index_i, -1, -1))
            map2 += list(range(0, index_i + 1))
        poly_vals = poly_val_2d[map1, map2]
        poly_vals = poly_vals[0:legendre_num]

        # poly_vals now has the value of all of the Legendre polynomials at our field point of interest.  So now we
        # need to multiply each value there with the value of the Legendre coefficients in each column and sum.  That
        # sum will be the value of a Zernike, so we loop to repeat that for each Zernike coefficient

        zernike_coeffs = np.zeros(zern_num)
        for z_index in range(0, zern_num):
            cur_legendre = self._field_dep_data.field(z_index)
            zernike_coeffs[z_index] = np.einsum('i, i->', cur_legendre[0:legendre_num], poly_vals)

        zernike_coeffs[0:3] = 0  # ignore piston/tip/tilt

        if assume_si_focus:
            # Assume SIs have been refocused, so there is no net defocus from the OTE field dep.
            # Take out the average focus per each SI, using precomputed values
            # These values computed using 'OTE As-Built Field Dependence Test.ipynb'
            avg_si_defocus_values = {'NIRCam':  -11.399957391866572,
                'NIRISS':  -31.412664805859624,
                'MIRI':  -32.074808824153386,
                'FGS':  -7.271566403275971,
                'NIRSpec':  -31.469227175636576,}  # these are in nm
            si_defocus_zern = avg_si_defocus_values[instrument]

            zernike_coeffs[3] -= si_defocus_zern

        return zernike_coeffs


    def move_seg_local(self, segment, xtilt=0.0, ytilt=0.0, clocking=0.0, rot_unit='urad',
                       radial=None, xtrans=None, ytrans=None, piston=0.0, roc=0.0, trans_unit='micron', display=False,
                       delay_update=False, absolute=False):
        """ Move a segment in pose and/or ROC, using segment-local control coordinates.

        These motions are always commanded in the segment-local "Control" coordinate systems,
        which are distinct for each segment.

        Parameters
        -----------
        segment : str
            Segment name, e.g. 'A1'.
        xtilt, ytilt, clocking : floats
            Tilt angle, in microradians by default. 'xtilt' means tilt *around the X axis*, and similarly for ytilt.
        radial, ytrans,xtrans : floats
            Displacement distance, in microns by default. Note the 'radial' and 'xtrans', 'ytrans' are redundant and
            included for convenience; the Ball WAS linear optical model uses radial translation as the control DoF, but
            physically that maps to either x or y translation depending on whether A, B, or C segment, and the Ball MCS
            algorithms expect X and Y translations.  We support both ways of describing this here.
        piston : float
            Displacement distance for piston.
        roc : float
            radius of curvature mechanism adjustment, in microns.
        trans_unit : str
            Unit for translations. Can be 'meter', 'micron', 'millimeter', 'nanometer', 'm', 'mm', 'nm', 'um'
        rot_unit : str
            Unit for rotations. Can be 'urad', 'radian', 'milliradian', 'arcsec', 'arcmin', 'milliarcsec'
        absolute : bool
            Same meaning as for JWST SURs: if true, move the segment to exactly this position. Otherwise moves
            are treated as incremental relative moves from the current position.
        display : bool
            Display after moving?
        delay_update : bool
            hold off on computing the WFE change? This is useful for computational efficiency if you're
            moving a whole bunch of segments at once. Incompatible with display=True.

        """

        if segment == 'SM':
            raise ValueError("SM not supported by move_seg_local. Use move_sm_local instead.")

        # Handle tilts and clocking
        tilts = np.array([xtilt, ytilt, clocking], dtype=float)
        if np.abs(tilts).sum() > 0:
            self.opd_header.add_history('Rotation: %s %s' % (str(tuple(tilts)), rot_unit))

            # sensitivity matrices are in microns per microradian
            # so convert all to urad.
            if rot_unit.endswith('s'):
                rot_unit = rot_unit[:-1]
            rot_unit = rot_unit.lower()
            if rot_unit == 'urad' or rot_unit == 'microrad' or rot_unit == 'microradian':
                pass
            elif rot_unit == 'milliarcsec':
                tilts *= (1e6 * np.pi / (180. * 60 * 60 * 1000))
            elif rot_unit == 'arcsec':
                tilts *= (1e6 * np.pi / (180. * 60 * 60))
            elif rot_unit == 'arcmin':
                tilts *= (1e6 * np.pi / (180. * 60))
            elif rot_unit == 'radian' or rot_unit == 'rad':
                tilts *= 1e6
            elif rot_unit == 'milliradian' or rot_unit == 'mrad':
                tilts *= 1e3
            else:
                raise NotImplementedError('unknown rot_unit')

        # Handle displacements. First special handling for radial vs x/y trans to allow both ways.
        if radial is None:
            if xtrans is None and ytrans is None:
                radial = 0
            else:
                if 'A' in segment:
                    radial = -1 * ytrans
                elif 'B' in segment:
                    radial = 1 * ytrans
                else:
                    radial = 1 * xtrans
        else:
            if xtrans is not None or ytrans is not None:
                raise RuntimeError("Cannot specify x/ytrans and radial at the same time.")

        vector = np.asarray([piston, radial], dtype=float)
        if np.abs(vector).sum() > 0:
            # influence functions are in microns WFE per micron, so convert all to microns
            if trans_unit.endswith('s'):
                trans_unit = trans_unit[:-1]
            trans_unit = trans_unit.lower()
            if trans_unit == 'micron' or trans_unit == 'um':
                pass
            elif trans_unit == 'mm' or trans_unit == 'millimeter':
                vector *= 1000
            elif trans_unit == 'nm' or trans_unit == 'nanometer' or trans_unit == 'nanometers':
                vector /= 1000
            elif trans_unit == 'm' or trans_unit == 'meter':
                vector *= 1e6
            else:
                raise ValueError("Unknown trans_unit for length: %s" % trans_unit)

            self.opd_header.add_history('Displacement: %s %s' % (str(tuple(vector)), trans_unit))

        # Handle ROC. ROC doesn't support any units conversions
        if np.abs(roc) != 0:
            if trans_unit == 'micron' or trans_unit == 'um':
                pass
            elif trans_unit == 'sag':
                roc *= 1e6
            self.opd_header.add_history('ROC: %s %s' % (roc, 'micron'))

        iseg = np.where(self.segnames == segment)[0][0]

        # Ordering = xtilt, ytilt, piston, clocking, rad trans, roc
        update_vector = [tilts[0], tilts[1], vector[0], tilts[2], vector[1], roc]
        if absolute:
            self.segment_state[iseg][:] = update_vector
        else:
            self.segment_state[iseg][:] += update_vector

        if not delay_update:
            self.update_opd(display=display)

    def move_seg_global(self, segment, xtilt=0.0, ytilt=0.0, clocking=0.0, rot_unit='urad',
                        radial=None, xtrans=0.0, ytrans=0.0, piston=0.0, roc=0.0, trans_unit='micron', display=False,
                        delay_update=False, absolute=False):
        """ Move a segment in pose and/or ROC, using PM global V coordinates..

        These motions are converted into the segment-local "Control" coordinate systems,
        which are distinct for each segment.

        Parameters
        -----------
        segment : str
            Segment name, e.g. 'A1'. Use 'SM' for the secondary mirror.
        xtilt, ytilt, clocking : floats
            Tilt angle, in microradians by default. 'xtilt' means tilt *around the X axis*, and similarly for ytilt.
        radial, ytrans,xtrans : floats
            Displacement distance, in microns by default. Note the 'radial' and 'xtrans', 'ytrans' are redundant and
            included for convenience; the Ball WAS linear optical model uses radial translation as the control DoF, but
            physically that maps to either x or y translation depending on whether A, B, or C segment, and the Ball MCS
            algorithms expect X and Y translations.  We support both ways of describing this here.
        piston : float
            Displacement distance for piston.
        roc : float
            radius of curvature mechanism adjustment, in microns.
        trans_unit : str
            Unit for translations. Can be 'micron', 'millimeter','nanometer', 'mm', 'nm', 'um'
        rot_unit : str
            Unit for rotations. Can be 'urad', 'radian', 'milliradian', 'arcsec', 'arcmin', 'milliarcsec'
        absolute : bool
            Same meaning as for JWST SURs: if true, move the segment to exactly this position. Otherwise moves
            are treated as incremental relative moves from the current position.
        display : bool
            Display after moving?
        delay_update : bool
            hold off on computing the WFE change? This is useful for computational efficiency if you're
            moving a whole bunch of segments at once. Incompatible with display=True.

        """

        # Handle tilts and clocking
        tilts = np.array([xtilt, ytilt, clocking], dtype=float)

        # Handle displacements. First special handling for radial vs x/y trans to allow both ways.
        # FIXME if radial provided, assign that value to either ytrans or xtrans
        if radial is not None:
            if 'A' in segment:
                ytrans = -1 * radial
            elif 'B' in segment:
                ytrans = 1 * radial
            else:
                xtrans = 1 * radial

        vector = np.asarray([xtrans, ytrans, piston], dtype=float)

        # convert to local coords

        local_tilt = tilts * 1.0

        self.move_seg_global(segment,
                             xtilt=local_tilt[0],
                             ytilt=local_tilt[1],
                             clocking=local_tilt[2], )
        # FIXME
        # FIXME - need to complete this fn!

        if np.abs(tilts).sum() > 0:
            self.opd_header.add_history('Rotation: %s %s' % (str(tuple(tilts)), rot_unit))

            # sensitivity matrices are in microns per microradian
            # so convert all to urad.
            if rot_unit.endswith('s'):
                rot_unit = rot_unit[:-1]
            rot_unit = rot_unit.lower()
            if rot_unit == 'urad':
                pass
            elif rot_unit == 'milliarcsec':
                tilts *= (1e6 * np.pi / (180. * 60 * 60 * 1000))
            elif rot_unit == 'arcsec':
                tilts *= (1e6 * np.pi / (180. * 60 * 60))
            elif rot_unit == 'arcmin':
                tilts *= (1e6 * np.pi / (180. * 60))
            elif rot_unit == 'radian' or rot_unit == 'rad':
                tilts *= 1e6
            elif rot_unit == 'milliradian' or rot_unit == 'mrad':
                tilts *= 1e3
            else:
                raise NotImplementedError('unknown rot_unit')

        # Handle displacements. First special handling for radial vs x/y trans to allow both ways.
        if radial is None:
            if xtrans is None and ytrans is None:
                radial = 0
            else:
                if 'A' in segment:
                    radial = -1 * ytrans
                elif 'B' in segment:
                    radial = 1 * ytrans
                else:
                    radial = 1 * xtrans
        else:
            if xtrans is not None or ytrans is not None:
                raise RuntimeError("Cannot specify x/ytrans and radial at the same time.")

        vector = np.asarray([piston, radial], dtype=float)
        if np.abs(vector).sum() > 0:
            # influence functions are in microns WFE per micron, so convert all to microns
            if trans_unit.endswith('s'):
                trans_unit = trans_unit[:-1]
            trans_unit = trans_unit.lower()
            if trans_unit == 'micron' or trans_unit == 'um':
                pass
            elif trans_unit == 'mm' or trans_unit == 'millimeter':
                vector *= 1000
            elif trans_unit == 'nm' or trans_unit == 'nanometer' or trans_unit == 'nanometers':
                vector /= 1000
            elif trans_unit == 'm' or trans_unit == 'meter':
                vector *= 1e6
            else:
                raise ValueError("Unknown trans_unit for length: %s" % trans_unit)

            self.opd_header.add_history('Displacement: %s %s' % (str(tuple(vector)), trans_unit))

        # Handle ROC. ROC doesn't support any units conversions
        if np.abs(roc) != 0:
            self.opd_header.add_history('ROC: %s %s' % (roc, 'micron'))

        iseg = np.where(self.segnames == segment)[0][0]

        # Ordering = xtilt, ytilt, piston, clocking, rad trans, roc
        update_vector = [tilts[0], tilts[1], vector[0], tilts[2], vector[1], roc]
        if absolute:
            self.segment_state[iseg][:] = update_vector
        else:
            self.segment_state[iseg][:] += update_vector

        if not delay_update:
            self.update_opd(display=display)

    def move_sm_local(self, xtilt=0.0, ytilt=0.0, rot_unit='urad',
                      xtrans=0.0, ytrans=0.0, piston=0.0, trans_unit='micron', display=False,
                      delay_update=False):
        """ Move the secondary mirror in pose, using segment-local control coordinates.

        These motions are always commanded in the segment-local "Control" coordinate systems,
        which are distinct for each segment. The SM is also handled a bit differently than all the
        PMSAs in terms of its local coordinate system, which this function handles behind the scenes.

        Parameters
        -----------
        segment : str
            Segment name, e.g. 'A1'. Use 'SM' for the secondary mirror.
        xtilt, ytilt, clocking : floats
            Tilt angle, in microradians by default. 'xtilt' means tilt *around the X axis*, and similarly for ytilt.
        xtrans, ytrans, piston : floats
            Displacement distance, in microns by default.
        trans_unit : str
            Unit for translations. Can be 'micron', 'millimeter','nanometer', 'mm', 'nm', 'um'
        rot_unit : str
            Unit for rotations. Can be 'urad', 'radian', 'milliradian', 'arcsec', 'arcmin', 'milliarcsec'
        display : bool
            Display after moving?
        delay_update : bool
            hold off on computing the WFE change? This is useful for computational efficiency if you're
            moving a whole bunch of segments at once. Incompatible with display=True.
        """

        # Handle tilts and clocking (which is always 0)
        tilts = np.array([xtilt, ytilt, 0]).astype(float)
        if np.abs(tilts).sum() > 0:
            self.opd_header.add_history('Rotation: %s %s' % (str(tuple(tilts)), rot_unit))

            # sensitivity matrices are in microns per microradian
            # so convert all to urad.
            if rot_unit.endswith('s'):
                rot_unit = rot_unit[:-1]
            rot_unit = rot_unit.lower()
            if rot_unit == 'urad':
                pass
            elif rot_unit == 'milliarcsec':
                tilts *= (1e6 * np.pi / (180. * 60 * 60 * 1000))
            elif rot_unit == 'arcsec':
                tilts *= (1e6 * np.pi / (180. * 60 * 60))
            elif rot_unit == 'arcmin':
                tilts *= (1e6 * np.pi / (180. * 60))
            elif rot_unit == 'radian' or rot_unit == 'rad':
                tilts *= 1e6
            elif rot_unit == 'milliradian' or rot_unit == 'mrad':
                tilts *= 1e3
            else:
                raise NotImplementedError('unknown rot_unit')

        # Handle displacements
        vector = np.asarray([xtrans, ytrans, piston])
        if np.abs(vector).sum() > 0:

            if trans_unit.endswith('s'):
                trans_unit = trans_unit[:-1]
            trans_unit = trans_unit.lower()
            if trans_unit == 'micron' or trans_unit == 'um':
                pass
            elif trans_unit == 'millimeters' or trans_unit == 'millimeter' or trans_unit == 'mm':
                vector *= 1000
            elif trans_unit == 'nm' or trans_unit == 'nanometer' or trans_unit == 'nanometers':
                vector /= 1000
            elif trans_unit == 'meter':
                vector *= 1000000
            else:
                raise ValueError("Unknown trans_unit for length: %s" % trans_unit)

            self.opd_header.add_history('Displacement: %s %s' % (str(tuple(vector)), trans_unit))

        iseg = 18

        self.segment_state[iseg][0] += tilts[0]  # xtilt
        self.segment_state[iseg][1] += tilts[1]  # ytilt
        self.segment_state[iseg][2] += vector[0]  # xtrans
        self.segment_state[iseg][3] += vector[1]  # ytrans
        self.segment_state[iseg][4] += vector[2]  # piston

        if not delay_update:
            self.update_opd(display=display)

    def move_global_zernikes(self, zvector, unit='micron',
                             absolute=False, delay_update=False, display=False):
        """ Add one or more aberrations specified arbitrarily as Zernike polynomials.
        This assumes no particular physics for the mirror motions, and allows adding
        any arbitrary WFE.

        Parameters
        ----------
        zvector : list or ndarray
            Zernike coefficients


        Note that the Zernikes are interpreted as being with respect to the
        *CIRCUMSCRIBING* circle.
        """

        if len(zvector) > len(self._global_zernike_coeffs):
            raise RuntimeError("Too many Zernike coefficients supplied. " +
                               "Need to increase length of global zernike coeffs vector in OTE_Linear_model_WSS.__init__")

        vector = np.asarray(zvector)
        # Convert to meters, since that's what the OPDs are in
        if unit.endswith('s'):
            unit = unit[:-1]
        unit = unit.lower()
        if unit == 'micron' or unit == 'um':
            vector *= 1e-6
        elif unit == 'mm' or unit == 'millimeter':
            vector *= 1e-3
        elif unit == 'nm' or unit == 'nanometer' or unit == 'nanometers':
            vector *= 1e-9
        elif unit == 'm' or unit == 'meter':
            pass
        else:
            raise ValueError("Unknown unit for Zernike wavefront RMS: %s" % unit)

        if absolute:
            self._global_zernike_coeffs *= 0
        for i in range(len(zvector)):  # for loop is inelegant but easily handles len(zvector) < len(coeffs).
            self._global_zernike_coeffs[i] += zvector[i]

        if not delay_update:
            self.update_opd(display=display)

    def move_sur(self, sur_file, group=None, verbose=False, reverse=False):
        """
        Move using a JWST Segment Update Request file

        Parameters
        ----------
        sur_file : file name, or SUR object instance
            Path to SUR XML file, or a webbpsf.surs.SUR object
        group : one-based int index
            Index to a single group to run. Default is to run all groups. Note,
            this index counts up from 1 (not 0) for consistency with group indexing
            in the SUR files themselves.

        verbose : bool
            Flag controlling whether moves are printed.
        reverse : bool
            Run this SUR "backwards", i.e. in opposite order of all groups and
            flipping the sign of all moves. (This can be useful for certain
            testing and mock data generation scenarios.)

        Returns
        -------

        """

        if isinstance(sur_file, surs.SUR):
            sur = sur_file
        else:
            sur = surs.SUR(sur_file)

        if group is not None:
            if group == 0:
                raise ValueError("Group indices start at 1, not 0.")
            groups = [sur.groups[group-1]]
        else:
            groups = sur.groups
        groupnum = list(np.arange(len(groups))+1)

        sign = -1 if reverse else 1
        if reverse:
            if verbose:
                print("Applying SUR in reverse: flipping group order and sign of all moves.")
            groups = reversed(groups)
            groupnum = reversed(groupnum)

        for igrp, grp in zip(groupnum, groups):
            if verbose:
                print("Moving segments for group {}".format(igrp))
            for update in grp:
                if verbose:
                    print("Move seg {} by {}".format(update.segment, str(update)))
                if update.type == 'pose':
                    if update.coord != 'local':
                        raise NotImplementedError("Only local moves supported!")

                    # FIXME - consider whether we should check for
                    # heterogeneous sets of units here...
                    rot_unit = update.units['X_TILT']
                    trans_unit = update.units['X_TRANS']

                    if update.segment == 'SM':
                        self.move_sm_local(xtilt=update.moves['X_TILT']*sign,
                                           ytilt=update.moves['Y_TILT']*sign,
                                           xtrans=update.moves['X_TRANS']*sign,
                                           ytrans=update.moves['Y_TRANS']*sign,
                                           piston=update.moves['PISTON']*sign,
                                           rot_unit=rot_unit,
                                           trans_unit=trans_unit,
                                           delay_update=True)
                    else:
                        self.move_seg_local(update.segment[0:2],
                                            xtilt=update.moves['X_TILT']*sign,
                                            ytilt=update.moves['Y_TILT']*sign,
                                            xtrans=update.moves['X_TRANS']*sign,
                                            ytrans=update.moves['Y_TRANS']*sign,
                                            piston=update.moves['PISTON']*sign,
                                            clocking=update.moves['CLOCK']*sign,
                                            absolute=update.absolute,
                                            rot_unit=rot_unit,
                                            trans_unit=trans_unit,
                                            delay_update=True)

                elif update.type == 'roc':
                    self.move_seg_local(update.segment[0:2],
                                        roc=update.moves['ROC']*sign,
                                        absolute=update.absolute,
                                        trans_unit=update.units['ROC'],
                                        delay_update=True)

                else:
                    raise NotImplementedError("Only moves of type='pose' or 'roc' are supported.")
        self.update_opd()

    #---- OPD perturbations from drifts are implemented in the next several functions ----
    def thermal_slew(self, delta_time, start_angle=-5,end_angle=45,
                     scaling=None, display=False, case='EOL', delay_update=False):
        """ Update the OPD based on presence of a pitch angle change between
        observations.

        Use a delta slew time along with the beginning and ending angles of the
        observatory relative to the sun (or the user can define a scaling factor)
        to determine the expected WFE caused by thermal variations.
        Note: The start_angle and end_angle are used together, but will be ignored
        if the scaling variable is set to something other than "None".

        The maximum HOT to COLD pitch angles are -5 to 45 degrees. With regards
        to this, we make some assumptions:
        1. A COLD to HOT slew is just the negative of the HOT to COLD slew
        2. The scaling factor can be simplified to a simple ratio of angles (this is
           a gross over-simplification due to lack of a better model)

        The HOT to COLD vs COLD to HOT nature of the slew is determined by the start
        and end angles

        Note, multiple calls in a row to this function are NOT cumulative; rather, the
        model internally resets to the initial starting OPD each time, and calculates
        a single slew. This is intentional to be more reproducible and well defined,
        with less hidden history state. If you need a more complex time evolution,
        build that yourself by summing individual delta OPDs.

        Parameters
        ----------
        delta_time: astropy.units quantity object
            The time between observations. Default units: "hour"
        start_angle: float
            The starting sun pitch angle, in degrees between -5 and +45
        end_angle: float
            The ending sun pitch angle, in degrees between -5 and +45
        scaling: float between 0 and 1
            Scaling factor that can be used instead of the start_angle
            and end_angle parameters. This directly sets the amplitude
            of the drift and overrides the angles and case settings.
        display: bool
            Display the updated OPD
        case : string
            either "BOL" for current best estimate at beginning of life, or
            "EOL" for more conservative prediction at end of life. The amplitude
            of the thermal drift is roughly 3x lower for BOL (13 nm after 14 days)
            versus EOL (43 nm after 14 days).
        delay_update: bool
            Users typically only need to call this directly if they have set the
            "delay_update" parameter to True in some function call to move mirrors.

        """
        # Check values
        if (start_angle < -5) or (end_angle > 45):
            raise ValueError("Start or end angle is outside of acceptable range of -5 to 45 degrees.")

        # Convert Delta time to units of days
        delta_time = convert_quantity(delta_time, to_units=u.day) #this returns astropy units quantity
        self.delta_time = delta_time.value
        self.start_angle = start_angle
        self.end_angle = end_angle
        self.scaling = scaling
        self._thermal_model.case = case

        # Update the header info
        self.opd_header['BUNIT'] = 'meter'
        self.opd_header['DELTA_T'] = (self.delta_time, "Delta time after slew [d]")
        self.opd_header['STARTANG'] = (self.start_angle, "Starting sun pitch angle [deg]")
        self.opd_header['ENDANG'] = (self.end_angle, "Ending sun pitch angle [deg]")
        self.opd_header['THRMCASE'] = (self._thermal_model.case, "Thermal model case, beginning or end of life")
        if scaling:
            self.opd_header['SCALING'] = (self.scaling, 'Scaling factor for delta slew')

        if not delay_update:
            self.update_opd(display=display)


    def _get_thermal_slew_coeffs(self, segid):
        """
        Get the WSS Hexike coefficients for the OPD describing the changes that have been
        caused by a change in pitch angle between two observations.

        Parameters:
        -----------
        segid: str
            Segment to be fit. 'SM' will fit the global focus term. Any other
            segment name will fits 9 Hexikes to that segment
        """
        if not self.scaling:
            num = np.sin(np.radians(self.end_angle)) - np.sin(np.radians(self.start_angle))
            den = np.sin(np.radians(45.)) - np.sin(np.radians(-5.))
            scaling = num / den

        else:
            scaling = self.scaling

        coeffs = self._thermal_model.get_coeffs(segid, self.delta_time)
        return scaling*coeffs


    def update_opd(self, display=False, verbose=False):
        """ Update the OPD based on the current linear model values.

        Users typically only need to call this directly if they have set the
        "delay_update" parameter to True in some function call to move mirrors.

        """

        # always start from the input OPD, then apply perturbations
        self.opd = self._opd_original.copy()
        self.opd_header.add_history('')
        self.opd_header.add_history('OTE linear model: updated OPD based on input segment poses')

        sm = 18
        sm_pose_coeffs = self.segment_state[sm].copy()[0:5]  # 6th row is n/a for SM
        sm_pose_coeffs.shape = (5, 1)  # to allow broadcasting below


        total_segment_state = self.segment_state + self._get_frill_drift_poses() + self._get_iec_drift_poses()

        for iseg, segname in enumerate(self.segnames[0:18]):
            pose_coeffs = total_segment_state[iseg].copy()
            if np.all(pose_coeffs == 0) and np.all(sm_pose_coeffs == 0) and self.delta_time==0:
                continue
            else:

                sensitivities = self._get_seg_sensitivities(segname)
                pose_coeffs.shape = (6, 1)  # to allow broadcasting in next line

                hexike_coeffs = sensitivities * pose_coeffs  # this will be a 6,9 array
                hexike_coeffs = hexike_coeffs.sum(axis=0)  # sum across all controlled modes
                # to yield a set of 9 overall coeffs

                sm_sensitivities = self._get_seg_sensitivities_from_sm(segname)
                hexike_coeffs_from_sm = sm_sensitivities * sm_pose_coeffs  # will be 5,9 array
                hexike_coeffs_from_sm = hexike_coeffs_from_sm.sum(axis=0)  # sum to get 9
                hexike_coeffs_from_thermal = self._get_thermal_slew_coeffs(segname)
                hexike_coeffs_combined = hexike_coeffs + hexike_coeffs_from_sm + hexike_coeffs_from_thermal

                if verbose:
                    print("Need to move segment {} by {} ".format(segname, pose_coeffs.flatten()))
                    print("plus SM moved by {} ".format(sm_pose_coeffs.flatten()))
                    print("plus segment moved by {} due to thermal contribution".format(hexike_coeffs_from_thermal))
                    print("   Hexike coeffs for {}: {}".format(segname, hexike_coeffs))

                self._apply_hexikes_to_seg(segname, hexike_coeffs_combined)

        # The thermal slew model for the SM global defocus is implemented as a global hexike.
        # So we have to combine that with the _global_hexikes array here
        global_hexike_coeffs_combined = self._global_hexike_coeffs.copy()
        if self.delta_time != 0.0:
            global_hexike_coeffs_combined[4] += self._get_thermal_slew_coeffs('SM')

        # Apply Global Zernikes, and/or hexikes
        if not np.all(global_hexike_coeffs_combined == 0):
            self._apply_global_hexikes(global_hexike_coeffs_combined)
        if not np.all(self._global_zernike_coeffs == 0):
            self._apply_global_zernikes()
            
        self._apply_field_dependence_model()

        if display:
            self.display()

    def _get_dynamic_opd(self, delta_time=14*u.day, case='EOL'):
        """Return ONLY the dynamic portion of the OPD, including thermal slew + frill + IEC
        Normally these are all folded up as part of the update_opd method, but it can be
        useful in certain circumstances to extract just this part, e.g. for plotting or analysis

        This is therefore a modified version of the update_opd function.
        This is a bit of a hack, as internal OPD state gets modified and then reset as a
        side effect of this.

        """
        # always start from the input OPD, then apply perturbations
        # self.opd = self._opd_original.copy()
        # self.opd_header.add_history('')
        # self.opd_header.add_history('OTE linear model: updated OPD based on input segment poses')
        self.opd *= 0  # start from a zero OPD and just add perturbations; we will undo this later!

        # Change a bunch of parameters (we will also undo this later!)
        thermal_slew_params = (self.delta_time, self.start_angle, self.end_angle, self._thermal_model.case)
        # default is max slew and 2 weeks
        self.start_angle= -5
        self.start_angle=  45
        self._thermal_model.case = case
        self.delta_time = delta_time.to_value(u.day)

        # sm = 18
        # sm_pose_coeffs = self.segment_state[sm].copy()[0:5]  # 6th row is n/a for SM
        # sm_pose_coeffs.shape = (5, 1)  # to allow broadcasting below
        sm_pose_coeffs = np.zeros( (5,1) )

        drift_segment_state = self._get_frill_drift_poses() + self._get_iec_drift_poses()

        for iseg, segname in enumerate(self.segnames[0:18]):
            pose_coeffs = drift_segment_state[iseg].copy()
            if np.all(pose_coeffs == 0) and np.all(sm_pose_coeffs == 0) and delta_time == 0:
                continue
            else:

                sensitivities = self._get_seg_sensitivities(segname)
                pose_coeffs.shape = (6, 1)  # to allow broadcasting in next line

                hexike_coeffs = sensitivities * pose_coeffs  # this will be a 6,9 array
                hexike_coeffs = hexike_coeffs.sum(axis=0)  # sum across all controlled modes
                # to yield a set of 9 overall coeffs

                sm_sensitivities = self._get_seg_sensitivities_from_sm(segname)
                hexike_coeffs_from_sm = sm_sensitivities * sm_pose_coeffs  # will be 5,9 array
                hexike_coeffs_from_sm = hexike_coeffs_from_sm.sum(axis=0)  # sum to get 9
                hexike_coeffs_from_thermal = self._get_thermal_slew_coeffs(segname)
                hexike_coeffs_combined = hexike_coeffs + hexike_coeffs_from_sm + hexike_coeffs_from_thermal

                self._apply_hexikes_to_seg(segname, hexike_coeffs_combined)

        # The thermal slew model for the SM global defocus is implemented as a global hexike.
        # So we have to combine that with the _global_hexikes array here
        global_hexike_coeffs_combined = self._global_hexike_coeffs.copy()
        if delta_time != 0.0:
            global_hexike_coeffs_combined[4] += self._get_thermal_slew_coeffs('SM')

        # Apply Global Zernikes, and/or hexikes
        if not np.all(global_hexike_coeffs_combined == 0):
            self._apply_global_hexikes(global_hexike_coeffs_combined)
        # if not np.all(self._global_zernike_coeffs == 0):
        #    self._apply_global_zernikes()

        #self._apply_field_dependence_model()

        drift_opd = self.opd.copy()  # save this, to return below

        # Important: Undo the overwriting of self.opd we just did
        self.delta_time, self.start_angle, self.end_angle, self._thermal_model.case = thermal_slew_params
        self.update_opd()

        return drift_opd

    def apply_frill_drift(self, amplitude=None, random=False, case='BOL', delay_update=False):
        """ Apply model of segment PTT motions for the frill-induced drift.

        This is additive with other WFE terms.

        Parameters
        ----------
        amplitude : float
            Amplitude of drift in nm rms to apply
        random : bool
            if True, choose a random amplitude from within the expected range
            for either the BOL or EOL cases. The assumed model is a uniform
            distribution between 0 and a maximum amplitude of 8.6 or 18.4 nm rms
            respectively.
        case : string
            either "BOL" for current best estimate at beginning of life, or
            "EOL" for more conservative prediction at end of life. Only relevant
            if random=True.
        delay_update : bool
            hold off on computing the WFE change? This is useful for computational efficiency if you're
            making many changes at once.
        """

        if random:
            if case.upper() == 'BOL':
                max_amp = 8.6
            elif case.upper() == 'EOL':
                max_amp = 18.4
            else:
                raise ValueError(f"Unknown value for parameter case: {case}.")

            amplitude = np.random.uniform(0, max_amp)
            _log.info(f"Applying random frill drift with amplitude {amplitude} nm rms (out of max {max_amp} nm rms).")
        elif amplitude is None:
            raise ValueError("if random=False, you must provide a value for the amplitude.")
        else:
            _log.info(f"Applying frill drift with amplitude {amplitude} nm rms.")
        self._frill_wfe_amplitude = amplitude

        if not delay_update:
            self.update_opd()

    def _get_frill_drift_poses(self):
        """ Return segment poses for current frill drift state
        """
        # These segment piston/tip/tilt misalignments are normalized to give 1 nm rms.
        # This segment state approximates the OTE in-flight prediction from John Johnston / Joe Howard.
        # Developed in "Generate Mock Flight Predicts.ipynb" by Perrin.
        ote_seg_motions_frill = np.array(
             [ [-0.00728 ,  0.      ,  0.00273 ,  0.      ,  0.      ,  0.      ],
               [ 0.00364 ,  0.00364 ,  0.00091 ,  0.      ,  0.      ,  0.      ],
               [-0.00455 ,  0.00182 , -0.00455 ,  0.      ,  0.      ,  0.      ],
               [ 0.      ,  0.      , -0.00273 ,  0.      ,  0.      ,  0.      ],
               [ 0.001092, -0.00364 , -0.0091  ,  0.      ,  0.      ,  0.      ],
               [ 0.00455 , -0.00455 ,  0.00455 ,  0.      ,  0.      ,  0.      ],
               [-0.00728 ,  0.      ,  0.00273 ,  0.      ,  0.      ,  0.      ],
               [ 0.00273 ,  0.      ,  0.00273 ,  0.      ,  0.      ,  0.      ],
               [-0.002275, -0.00455 ,  0.00728 ,  0.      ,  0.      ,  0.      ],
               [-0.00273 ,  0.      ,  0.009555,  0.      ,  0.      ,  0.      ],
               [-0.00273 , -0.00091 ,  0.01183 ,  0.      ,  0.      ,  0.      ],
               [ 0.      , -0.0091  , -0.00091 ,  0.      ,  0.      ,  0.      ],
               [ 0.0091  , -0.00182 ,  0.0182  ,  0.      ,  0.      ,  0.      ],
               [ 0.      ,  0.      , -0.00455 ,  0.      ,  0.      ,  0.      ],
               [ 0.002275,  0.00455 ,  0.01183 ,  0.      ,  0.      ,  0.      ],
               [ 0.00273 ,  0.      ,  0.009555,  0.      ,  0.      ,  0.      ],
               [-0.00182 ,  0.00364 ,  0.00728 ,  0.      ,  0.      ,  0.      ],
               [-0.00273 ,  0.      ,  0.00273 ,  0.      ,  0.      ,  0.      ],
               [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])/16

        return ote_seg_motions_frill * self._frill_wfe_amplitude


    def apply_iec_drift(self, amplitude=None, random=False, case='BOL', delay_update=False):
        """ Apply model of segment PTT motions for the drift seen at OTIS induced by the IEC
        (Instrument Electronics Compartment) heater resistors. This effect was in part due to
        non-flight-like ground support equipment mountings, and is not expected in flight at
        the same levels it was seen at JSC. We model it anyway, at an amplitude consistent with
        upper limits for flight.

        This is additive with other WFE terms.

        Parameters
        ----------
        amplitude : float
            Amplitude of drift in nm rms to apply
        random : bool
            if True, choose a random amplitude from within the expected range
            for either the BOL or EOL cases. The assumed model is a sinusoidal drift
            between 0 and 3.5 nm, i.e. a random variate from the arcsine distribution
            times 3.5.
        case : string
            either "BOL" for current best estimate at beginning of life, or
            "EOL" for more conservative prediction at end of life. Only relevant
            if random=True. (Note, for IEC drift the amplitude is the same regardless.)
        delay_update : bool
            hold off on computing the WFE change? This is useful for computational efficiency if you're
            making many changes at once.
        """

        # These segment piston/tip/tilt misalignments are normalized to give 1 nm rms.
        if random:
            max_amp = 3.5  # regardless of EOL or BOL

            amplitude = scipy.stats.arcsine().rvs() * max_amp
            _log.info(f"Applying random IEC-induced drift with amplitude {amplitude} nm rms (out of max {max_amp} nm rms).")
        elif amplitude is None:
            raise ValueError("if random=False, you must provide a value for the amplitude.")
        else:
            _log.info(f"Applying IEC-induced drift with amplitude {amplitude} nm rms.")
        self._iec_wfe_amplitude = amplitude

        if not delay_update:
            self.update_opd()

    def _get_iec_drift_poses(self):

        # These segment piston/tip/tilt motions approximate the OTE observed
        # oscillation seen at JSC OTIS cryo vac test, due to IEC heater thermal loading
        # through GSE paths.

        # Developed in the PFR190 study by Perrin & Telfer as
        # oscillation_opd_case1_var2_ptt.fits.gz; decomposed into segment PTT by Perrin
        # in "Make PSF Sim Library for JWST cycle 1 props - Development.ipynb"

        ote_seg_motions_iec = np.array(
              [[ 0.5211,  0.2925, -0.3567,  0.    ,  0.    ,  0.    ],
               [-0.0652,  0.2788,  0.1896,  0.    ,  0.    ,  0.    ],
               [ 0.2491, -0.2203,  0.1522,  0.    ,  0.    ,  0.    ],
               [ 0.4078,  0.1386,  0.0301,  0.    ,  0.    ,  0.    ],
               [-0.0663, -0.0702,  0.2119,  0.    ,  0.    ,  0.    ],
               [ 0.3488, -0.3986, -0.2216,  0.    ,  0.    ,  0.    ],
               [ 0.4363, -0.4034, -0.5762,  0.    ,  0.    ,  0.    ],
               [-0.584 , -0.4198, -0.0271,  0.    ,  0.    ,  0.    ],
               [ 1.1319, -0.195 ,  0.8383,  0.    ,  0.    ,  0.    ],
               [ 0.3116, -0.4376,  0.4631,  0.    ,  0.    ,  0.    ],
               [-0.0052,  0.6276, -0.1083,  0.    ,  0.    ,  0.    ],
               [ 0.1727,  0.6529, -0.5457,  0.    ,  0.    ,  0.    ],
               [-0.3807, -0.4189, -0.602 ,  0.    ,  0.    ,  0.    ],
               [-0.4924, -0.1094,  0.0996,  0.    ,  0.    ,  0.    ],
               [ 1.0861, -0.1953,  0.8371,  0.    ,  0.    ,  0.    ],
               [ 0.4202, -0.6961,  0.3543,  0.    ,  0.    ,  0.    ],
               [ 0.7644,  0.6466, -0.0861,  0.    ,  0.    ,  0.    ],
               [ 0.1887,  0.0825, -0.7851,  0.    ,  0.    ,  0.    ],
               [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ]])/1000
        return ote_seg_motions_iec * self._iec_wfe_amplitude

    def header_keywords(self):
        """ Return info we would like to save in FITS header of output PSFs
        """
        keywords = OrderedDict()
        keywords['OTETHMDL'] = (self._thermal_model.case, "OTE Thermal slew model case")
        keywords['OTETHSTA'] = (self.start_angle, "OTE Starting pitch angle for thermal slew model")
        keywords['OTETHEND'] = (self.start_angle, "OTE Ending pitch angle for thermal slew model")
        keywords['OTETHRDT'] = (self.delta_time, "OTE Thermal slew model delta time after slew")
        keywords['OTETHRWF'] = (self._thermal_wfe_amplitude, "OTE WFE amplitude from 'thermal slew' term")
        keywords['OTEFRLWF'] = (self._frill_wfe_amplitude, "OTE WFE amplitude from 'frill tension' term")
        keywords['OTEIECWF'] = (self._iec_wfe_amplitude, "OTE WFE amplitude from 'IEC thermal cycling' term")
        keywords.update(self.meta)
        return keywords


################################################################################


def enable_adjustable_ote(instr):
    """
    Set up a WebbPSF instrument instance to have a modifiable OTE
    wavefront error OPD via an OTE linear optical model (LOM).

    Parameters
    ----------
    inst : WebbPSF Instrument instance
        an instance of one of the WebbPSF instrument classes.

    Returns
    --------
    a modified copy of that instrument set up to use the LOM, and
    the associated instance of the LOM.

    """

    # if the instrument already is set up for an adjustable OTE model, then no op and return that
    if isinstance(instr.pupilopd, OTE_Linear_Model_WSS):
        return instr, instr.pupilopd

    import copy
    instcopy = copy.copy(instr)
    if instr.pupilopd is None:
        opdpath = None
    elif isinstance(instr.pupilopd, fits.HDUList):
        opdpath = instr.pupilopd
    else:
        opdpath = instr.get_opd_file_full_path(instr.pupilopd)

    pupilpath = instr.pupil

    name = "Modified OPD from " + str(instr.pupilopd)
    opd = OTE_Linear_Model_WSS(name=name,
                               opd=opdpath, transmission=pupilpath)

    opd.v2v3 = instr._tel_coords()  # copy field location, for use in field dependence models.
    instcopy.pupilopd = opd
    instcopy.pupil = opd

    return instcopy, opd


def setup_image_array(ote, radius=1, size=None, inverted=False, reset=False, verbose=False,
                      acfs_only=False,
                      guide_seg=None, guide_radius=10.0):
    """
    Apply tilts to put the segments in an image array configuration.

    Parameters
    -----------
    ote : OTE_Linear_Model_WSS instance
        The telescope model to be adjusted. This will be modified by this function.
    radius : float
        Desired radius to B segments, in arcseconds as seen on the focal plane.
    size : string, optional
        Another way of specifying the image array size. Use one of 'small', 'medium', or
        'large' for the standard sizes used in OTE commissioning. Or 'cmimf' for the size
        used in Coarse MIMF, which is in between small and medium.
    guide_seg : string
        relevant mostly for coarse MIMF and image stacking. Kick out a segment to guide?
        The segment will be offset during coarse MIMF to its position desired during GA2, i.e.
        to its position in the large array.
    acfs_only : bool
        Only tilt the ACFs (applicable to JSC OTIS only)
    inverted : bool
        Invert the array
    reset : bool
        Discard any previously applied mirror moves if true
    verbose : bool
        Print more text about moves.
    """

    assert isinstance(ote, OTE_Linear_Model_WSS), "First argument has to be a linear optical model instance."

    nircam_pixelscale = 0.0311
    standard_sizes = {'small': 80 * nircam_pixelscale,
                      'cmimf': 120 * nircam_pixelscale,
                      'medium': 300 * nircam_pixelscale,
                      'large': 812 * nircam_pixelscale
                      }

    if size is not None:
        radius = standard_sizes[size]

    # how many microradians of segment tilt per arcsecond of PSF motion?
    # note factor of 2 since reflection
    arcsec_urad = (1 * u.arcsec).to(u.urad).value / 2

    if reset:
        ote.reset()

    # Image Arrays used in flight
    size = radius * -1 if inverted else radius
    for i in range(1, 7):
        ote.move_seg_local('A' + str(i), xtilt=size * arcsec_urad / 2, delay_update=True)
        ote.move_seg_local('B' + str(i), xtilt=-size * arcsec_urad, delay_update=True)
        ote.move_seg_local('C' + str(i), ytilt=size * np.sqrt(3) / 2 * arcsec_urad, delay_update=True)

    if guide_seg is not None:
        # Undo the regular tilt for this segment, and then move it to
        # the side.
        if 'A' in guide_seg:
            xtilt = (-size + standard_sizes['large']) / 2 * arcsec_urad
            ytilt = 0
        elif 'B' in guide_seg:
            xtilt = -size * arcsec_urad
            ytilt = guide_radius * arcsec_urad
        else:
            ytilt = (size * np.sqrt(3) / 2 + guide_radius) * arcsec_urad
            xtilt = 0
        ote.move_seg_local(guide_seg, xtilt=xtilt, ytilt=ytilt, delay_update=True)

    ote.update_opd(verbose=verbose)


def random_unstack(ote, radius=1, verbose=False):
    """ Unstack the segments by randomly perturbing them in tip and tilt.


    Parameters
    ----------
    radius : float
        Scale the perturbations to have this value for 1-sigma per axis
    """
    assert isinstance(ote, OTE_Linear_Model_WSS), "First argument has to be a linear optical model instance."

    # how many microradians of segment tilt per arcsecond of PSF motion?
    # note factor of 2 since reflection
    arcsec_urad = (1 * u.arcsec).to(u.urad).value / 2

    xoffsets = np.random.randn(18) * radius
    yoffsets = np.random.randn(18) * radius

    for i, seg in enumerate(constants.SEGNAMES):
        ote.move_seg_local(seg, xtilt=xoffsets[i] * arcsec_urad,
                           ytilt=yoffsets[i] * arcsec_urad, delay_update=True)

    ote.update_opd(verbose=verbose)

# --------------------------------------------------------------------------------

def test_OPDbender():
    plt.figure(1)
    tel = OPDbender()
    tel.displace('A1', 1, 0, 0, display=False)
    tel.displace('A2', 0, 1, 0, display=False)
    tel.displace('A3', 0, 0, .03, display=False)
    tel.displace('A4', 0, -1, 0, display=False)
    tel.displace('A5', 1, -1, 0, display=False)

    tel.tilt('B1', .1, 0, 0)
    tel.tilt('B2', 0, .1, 0)
    tel.tilt('B3', 0, 0, 100)

    tel.display()

    plt.figure(2)
    tel.zern_seg('B3')

    print("")
    print("")
    print("RMS WFE is ", tel.rms())

    tel.print_state()


def test2_OPDbender(filename='OPD_RevV_nircam_132.fits'):
    orig = OPDbender(filename)

    plot_kwargs = {'colorbar_orientation': 'horizontal', 'clear': False}

    plt.clf()
    plt.subplot(131)
    orig.draw(title='Input OPD from \n' + filename, **plot_kwargs)

    perturbed = orig.copy()
    perturbed.perturb_all(multiplier=0.2, draw=False)

    plt.subplot(132)
    perturbed.draw(title='OPD after small random perturbation', **plot_kwargs)

    plt.subplot(133)
    diff = (perturbed - orig)
    diff.draw(title='Difference ({0:.1f} nm rms)'.format(diff.rms()), **plot_kwargs)


#-------------------------------------------------------------------------------
# Thermal


class OteThermalModel(object):
    """
    Create an object for a delta_time that predicts the WSS Hexike coefficients
    for an OPD that represents the impact of thermal variation caused by a change
    in pitch angle relative to the sun.

    Given a time in units of seconds, minutes, hours, or days.

    Parameters:
    -----------
    delta_time: tuple, (number, astropy.units quantity object)
        Include the number and units for the delta time between observations.
    case : string, 'BOL' or 'EOL'
        Model case for drift amplitudes. As of the 2020 PSR modeling cycle, the beginning of life
        has 13 nm rms after about a week, EOL has about 43 nm rms after a week.

    Returns:
    --------
    coeffs: array-like
        WSS Hexike Coefficients for JWST OPD based on thermal variations over
        delta_time

    """
    def __init__(self, case='BOL'):
        """
        Set up the object such that it can be used for any time, delta_time
        """
        self.nterms = 9
        # Load fitting values table:
        mypath = os.path.dirname(os.path.abspath( __file__ ))+os.sep
        # This table is in units of microns
        self._fit_file = os.path.join(mypath, 'otelm', 'thermal_OPD_fitting_parameters_9H_um.fits')
        self._fit_data = fits.getdata(self._fit_file)
        self.case = case


    @staticmethod
    def second_order_thermal_response_function(x, tau_1, gn_1, tau_2, gn_2):
        """ Second order thermal response function """
        return gn_1 * (1 - np.exp(-1 * (x / tau_1))) + gn_2 * (1 - np.exp(-1 * (x / tau_2)))


    def check_units(self, coeffs):
        """
        Make sure that the coefficients are in the correct units - meters.
        (Adapted from poppy.poppy_core.FITSOpticalElement)
        """
        header = fits.getheader(self._fit_file, ext=1)
        opdunits = header['BUNIT']
        # normalize and drop any trailing 's'
        opdunits = opdunits.lower()
        if opdunits.endswith('s'):
            opdunits = opdunits[:-1]
        if opdunits in ('meter', 'm'):
            pass
        elif opdunits in ('micron', 'um', 'micrometer'):
            coeffs *= 1e-6
        elif opdunits in ('nanometer', 'nm'):
            coeffs *= 1e-9
        return coeffs


    def get_coeffs(self, segid, delta_time):
        """ Given the segid name (either 'SM' or any of the segment names under
        constants.SEGNAMES), return the global or local (to each segment) Hexike
        coefficients

        Assume that delta_time is a float in units of days.
        """
        if delta_time == 0.0:
            if segid == 'SM':
                return 0.0
            else:
                return np.zeros(self.nterms)
        else:
            coeffs = OteThermalModel.second_order_thermal_response_function(delta_time,
                                         self._fit_data[self._fit_data['segs'] == segid]['tau1'],
                                         self._fit_data[self._fit_data['segs'] == segid]['Gn1'],
                                         self._fit_data[self._fit_data['segs'] == segid]['tau2'],
                                         self._fit_data[self._fit_data['segs'] == segid]['Gn2'])
            if len(coeffs) == 0:
                _log.warning("Invalid segment ID. No coefficients returned")
                coeffs = 0.0
            elif segid == 'SM':
                coeffs = self.check_units(coeffs[0])
            else:
                coeffs = self.check_units(coeffs)

            if self.case.upper()=='BOL':
                # Beginning of life predictions as of 2020 have much lower amplitude WFE drift
                # than the EOL model that was fit to produce the coefficients here.
                coeffs *= 0.35

            return coeffs


def convert_quantity(input_quantity, from_units=None, to_units=u.day):
    """
    Convert an input quantity (expecting an astropy units quantity), to a
    specified output quantity.

    (This defaults to units of time but can be used for any quantity as long
    as both from_units and to_units are set. If the from_units is not set, it
    will assume units of hours.)

    Parameters:
    -----------
    input_quantity: astropy units quantity, int/float
        Give an input quantity as an astropy units quantity or int/float.
        If the user passes in an int or float, units of *HOURS* will be assumed.
    from_units: astropy unit
        If input_quantity is not an astropy units quantity, this needs to be
        set, otherwise a unit of hours is assumed, regardless of to_units parameter
    to_units: astropy unit
        Default: u.day
        Set the astropy unit the convert to.

    Returns:
    --------
    output_quantity: astropy units quantity
        Return a an astropy units quantity in units set by to_units
    """
    try:
        output_quantity = input_quantity.to(to_units)
    except AttributeError:
        if from_units:
            input_quantity *= from_units
        else:
            input_quantity *= u.hour
        output_quantity = input_quantity.to(to_units)

    return output_quantity


#--------------------------------------------------------------------------------
# WFE decomposition


class JWST_WAS_PTT_Basis(object):
    def __init__(self):
        """ Segment piston/tip/tilt basis using the same conventions as JWST WAS
        i.e. local mechanical control coordinates per each segment and its local
        orientation.

        Similar to poppy.zernike.Segment_PTT_Basis, but:
            (a) specifically matches the JWST aperture geometry exactly, and
            (b) matches the local control coordinates for JWST segment controls.

        Useful for decomposing WFE maps into segment piston, tip, tilts.
        See poppy.zernike.opd_expand_segments()
        and coeffs_to_seg_state() in this file.

        See also JWST_WAS_Full_Basis, which includes the other three degrees of freedom

        """

        # Internally this is implemented as a wrapper on an OTE Linear WFE model;
        # we use the degrees of freedom of that model directly as the basis functions here

        self.ote = OTE_Linear_Model_WSS()
        self.nsegments=18

    def aperture(self):
        """ Return the overall aperture across all segments """
        return self.ote.amplitude

    def __call__(self, nterms=None, npix=1024, outside=np.nan):
        """ Generate PTT basis ndarray for the specified aperture

        Parameters
        ----------
        nterms : int
            Number of terms. Set to 3x the number of segments.
        npix : int
            Size, in pixels, of the aperture array.
        outside : float
            Value for pixels outside the specified aperture.
            Default is `np.nan`, but you may also find it useful for this to
            be 0.0 sometimes.

        """
        if nterms is None:
            nterms = 3*self.nsegments
        elif nterms > 3*self.nsegments:
            raise ValueError("nterms must be <= {} for the specified segment aperture.".format(3*self.nsegments))

        # Re-use the machinery inside the OTE Linear model class class to set up the
        # arrays defining the segment and zernike geometry.

        # If multiple calls to this function use different values for npix, we may have to 
        # update/reset the OTE LOM instance here.
        if self.ote.opd.shape[0] != npix:
            self.ote = OTE_Linear_Model_WSS(npix=npix)

        # For simplicity we always generate the basis for all the segments
        # even if for some reason the user has set a smaller nterms.
        basis = np.zeros((self.nsegments*3, npix, npix))
        basis[:] = outside
        for i, segname in enumerate(self.ote.segnames[0:18]):
            # We do these intentionally with the base units, though those result in unphysically large moves

            iseg = i+1
            wseg = np.where(self.ote._segment_masks==iseg)

            # Piston
            self.ote.zero()
            self.ote.move_seg_local(segname, piston=1, trans_unit='meter')
            basis[i*3][wseg] = self.ote.opd[wseg]

            # Tip
            self.ote.zero()
            self.ote.move_seg_local(segname, xtilt=1, rot_unit='radian')
            basis[i*3+1][wseg] = self.ote.opd[wseg]

            #Tilt
            self.ote.zero()
            self.ote.move_seg_local(segname, ytilt=1, rot_unit='radian')
            basis[i*3+2][wseg] = self.ote.opd[wseg]

        return basis[0:nterms]


class JWST_WAS_Full_Basis(object):
    def __init__(self):
        """ Segment pose full basis using the same conventions as JWST WAS
        i.e. local mechanical control coordinates per each segment and its local
        orientation.

        Similar to JWST_WAS_PTT_Basis, but:
            - includes clocking, radial translation, and radius of curvature degrees of freedom too

            (Note, azimuthal translation is intentionally not controlled in the JWST
            alignment schema, but is left as a redundant degree of freedom given the azimuthal
            rotational symmetry of the observatory.)

        Useful for decomposing WFE maps into segment piston, tip, tilts, translations, RoC.
        See poppy.zernike.opd_expand_segments()
        and coeffs_to_seg_state() in this file.

        See also JWST_WAS_PTT_Basis, which includes just the piston, tip, tilt degrees of freedom
        """

        # Internally this is implemented as a wrapper on an OTE Linear WFE model;
        # we use the degrees of freedom of that model directly as the basis functions here

        self.ote = OTE_Linear_Model_WSS()
        self.nsegments=18

    def aperture(self):
        """ Return the overall aperture across all segments """
        return self.ote.amplitude

    def __call__(self, nterms=None, npix=1024, outside=np.nan):
        """ Generate basis ndarray for the specified aperture

        Parameters
        ----------
        nterms : int
            Number of terms. Set to 6x the number of segments.
        npix : int
            Size, in pixels, of the aperture array.
        outside : float
            Value for pixels outside the specified aperture.
            Default is `np.nan`, but you may also find it useful for this to
            be 0.0 sometimes.

        """
        ndof = 6

        if nterms is None:
            nterms = ndof*self.nsegments
        elif nterms > ndof*self.nsegments:
            raise ValueError("nterms must be <= {} for the specified segment aperture.".format(3*self.nsegments))

        # Re-use the machinery inside the OTE Linear model class class to set up the
        # arrays defining the segment and zernike geometry.

        # For simplicity we always generate the basis for all the segments
        # even if for some reason the user has set a smaller nterms.
        basis = np.zeros((self.nsegments*6, npix, npix))
        basis[:] = outside
        for i, segname in enumerate(self.ote.segnames[0:18]):
            # We do these intentionally with the base units, though those result in unphysically large moves

            iseg = i+1
            wseg = np.where(self.ote._segment_masks==iseg)

            # Piston
            self.ote.zero()
            self.ote.move_seg_local(segname, piston=1, trans_unit='meter')
            basis[i*ndof][wseg] = self.ote.opd[wseg]

            # Tip
            self.ote.zero()
            self.ote.move_seg_local(segname, xtilt=1, rot_unit='radian')
            basis[i*ndof+1][wseg] = self.ote.opd[wseg]

            # Tilt
            self.ote.zero()
            self.ote.move_seg_local(segname, ytilt=1, rot_unit='radian')
            basis[i*ndof+2][wseg] = self.ote.opd[wseg]

            # Clocking
            self.ote.zero()
            self.ote.move_seg_local(segname, clocking=1, rot_unit='radian')
            basis[i * ndof + 3][wseg] = self.ote.opd[wseg]

            # Radial Translation
            self.ote.zero()
            self.ote.move_seg_local(segname, radial=1, trans_unit='meter')
            basis[i * ndof + 4][wseg] = self.ote.opd[wseg]

            # Radius of Curvature
            self.ote.zero()
            self.ote.move_seg_local(segname, roc=1, trans_unit='meter')
            basis[i * ndof + 5][wseg] = self.ote.opd[wseg]
        return basis[0:nterms]


def coeffs_to_seg_state(coeffs):
    """ Convert coefficients from Zernike fit to OTE linear model segment state

    Unit conversion and axis index reordering.

    Example usage:

    coeffs = poppy.zernike.opd_expand_segments(some_opd, aperture=ote.amplitude, basis=jw_ptt_basis, nterms=54)
    ote.segment_state = coeffs_to_seg_state(coeffs)


    """
    seg_state = np.zeros((18,6))
    coeffs_tab = coeffs.reshape(18,3)
    seg_state[:,2] = coeffs_tab[:,0]  # piston is 3rd column
    seg_state[:,0] = coeffs_tab[:,1]  # tip in 1st column
    seg_state[:,1] = coeffs_tab[:,2]  # tilt in 2nd column
    return seg_state*1e6   # convert from meters & radians to micro units


@functools.lru_cache
def _get_lom(npix):
    """Initialize a few things that are computationally expensive so we cache for multiple reuses"""

    ote_lom = OTE_Linear_Model_WSS(npix=npix)
    jw_ptt_basis = JWST_WAS_PTT_Basis()
    return ote_lom, jw_ptt_basis


def decompose_opd_segment_PTT(opd, plot=False, plot_vmax=None):
    """Phase decomposition of an OPD into PMSA piston, tip, tilt modes

    Parameters
    ----------
    opd : 2d ndarray
        OPD array
    plot : bool
        Display diagnostic plots in addition to doing the fit
    plot_vmax : float
        If plot is used, this allows you to adjust the vmin/vmax in the plot.


    Returns
    -------

    fit_opd : 2d ndarray
        Projection of the input OPD into JWST PTT modes
    coeffs : float array
        Coefficients per mode, in order corresponding to the webbpsf.opds.JWST_WAS_PTT_Basis class,
        which is {piston, xtilt, ytilt} repeated per segments in order

    """
    npix = opd.shape[0]

    ote_lom, jw_ptt_basis = _get_lom(npix)

    if ote_lom.opd.shape[0] != npix:
        ote_lom = OTE_Linear_Model_WSS(npix=npix)

    combined_mask = ((ote_lom.amplitude != 0) & np.isfinite(opd))


    coeffs = poppy.zernike.opd_expand_segments(opd, aperture=combined_mask,
                                               basis=jw_ptt_basis, nterms=54,iterations=4)
    fit = poppy.zernike.compose_opd_from_basis(coeffs, basis=jw_ptt_basis, npix=npix)

    fit[~combined_mask]=np.nan
    if plot:
        fig, ax = plt.subplots(figsize=(16,4), nrows=1, ncols=1)
        if not plot_vmax:
            plot_vmax = np.abs(opd).max()
        masked_opd = opd.copy()
        masked_opd[~combined_mask]=np.nan
        webbpsf.trending.show_opd_image(np.hstack((masked_opd, fit, opd-fit)), ax=ax, vmax=plot_vmax, labelrms=False )
        plt.colorbar(mappable=ax.images[0])
        plt.title('OPD, fit to segment PTT terms, and residuals')

    return fit, coeffs


def sur_to_opd(sur_filename, ignore_missing=False, npix=256):
    """Utilty function to load a SUR and compute delta OPD"""
    ote = OTE_Linear_Model_WSS(npix=npix)

    if not os.path.exists(sur_filename):
        if not ignore_missing:
            raise FileNotFoundError(f"Missing SUR: {sur_filename}. Download of these should eventually be automated; for now, manually retrieve from WSSTAS at https://wsstas.stsci.edu/wsstas/staticPage/showContent/RecentSURs?primary=master.png")
        else:
            return np.zeros((npix,npix), float)
    ote.move_sur(sur_filename)
    return ote.opd * 1e6 # convert from meters to microns



#--------------------------------------------------------------------------------
# Coarse track pointing (for early commissioning simulations)

def get_coarse_blur_parameters(t0, duration, pixelscale, plot=False, case=1,):
    """ Extract coarse blur center offset and convolution kernel from the Coarse Point sim time series

    Parameters
    -----------

    pcsmodel : astropy.Table
        High resolution time series data from JWST ACS sims
    t0 : float
        Start time, in seconds, for the time period of interest (typically the exposure start time.)
    duration : float
        Exposure duration, in seconds
    pixelscale : float
        Pixelscale in arcsec/pix for the output convolution kernel. Typically the NIRCam
        detector pixel scale, or an oversampled version thereof.
    case : int
        Which model output case to use.
        Model 1 has lower drift rate and yields approximately Gaussian jitter with 1 sigma = 0.15 per axis
        Model 2 has a larger linear drift/trend of about 2 arcsec over the 2 hours.

    Returns
    -------
    cen : 2-tuple of floats
        Mean offset in V2,V3 during the exposure
    kernel : 2D ndarray
        Convolution kernel to pass to WebbPSF, generated from the LOS model during the observation
        sampled/rasterized into the specified pixel scale.
    """

    pcsmodel = astropy.table.Table.read(os.path.join(__location__, 'otelm', f'coarse_track{case}_sim_pointing.fits'))

    wt = (t0 < pcsmodel['time']) & (pcsmodel['time'] < t0+duration)
    ns = wt.sum()

    # Extract coordinates for the requested time period
    coords = np.zeros((2,wt.sum()), float )
    coords[0] = pcsmodel['deltaV2'][wt]
    coords[1] = pcsmodel['deltaV3'][wt]

    cen = coords.mean(axis=1)       # Center

    dc = (coords-cen.reshape(2,1) )   # differential coords, in arcsec

    # Set up box to raster the curve into
    halfbox = np.ceil(np.abs(dc).max()/pixelscale)
    boxsize = int(2*halfbox+1) # must be an odd number for astropy convolution
    kernel = np.zeros((boxsize, boxsize))

    # Compute coords relative to lower right corner of raster box
    lowerright = -halfbox*pixelscale
    dc_pixels = np.array(np.round((dc - lowerright)/pixelscale), int)

    # Raster the curve into the array
    # have to do this via for loop rather than array indexing, to handle repeated indices
    for x, y in dc_pixels.transpose():
        kernel[y, x] += 1

    # optional display
    if plot:
        plt.figure()
        plt.plot(pcsmodel['deltaV2'], pcsmodel['deltaV3'], label='every 0.25 s')
        plt.plot(coords[0], coords[1], label='during exposure')

        plt.plot(pcsmodel['deltaV2'][0], pcsmodel['deltaV3'][0], marker='*', color='cyan', label='start of time series')
        plt.plot(cen[0], cen[1], marker='*', color='black', label='mean in exposure')

        plt.gca().set_aspect('equal')
        plt.legend(fontsize=7)
        plt.xlabel("Delta V2 [mas]")
        plt.ylabel("Delta V3 [mas]")


        plt.figure()
        plt.plot(dc[0], dc[1], color='C1', label='during exposure' )
        plt.plot((dc_pixels[0]-halfbox)*pixelscale, (dc_pixels[1]-halfbox)*pixelscale, color='C2', label='rounded to pixels')
        plt.gca().set_aspect('equal')
        plt.legend(fontsize=7)
        plt.xlabel("Delta V2 [mas]")
        plt.ylabel("Delta V3 [mas]")


        plt.figure()
        plt.imshow(kernel, cmap = matplotlib.cm.gray, origin='lower')
        plt.title(f"Convolution kernel at t={t0}, d={duration} s\nOffset={cen} arcsec", fontsize=10)
        plt.ylabel('Delta V3 [pixels]')

    return cen, kernel

#--------------------------------------------------------------------------------
# Wavefront decomposition and related


def get_rms_per_segment(opd, plot=False):
    """ Calculate RMS WFE per segment

    Parameters
    ----------
    opd : 2d float ndarray
        OPD. Assumed to be in units of meters
    plot : bool
        Plot images for diagnostics?

    Returns
    -------
    rms_per_seg : dict
        Dict of the form {"A1": 75.2, [etc]} mapping string segment names to
        RMS in nanometers per each segment.

    """

    npix=opd.shape[0]

    segmap_fn = os.path.join(utils.get_webbpsf_data_path(), f"JWpupil_segments_RevW_npix{npix}.fits.gz")
    segmap = fits.getdata(segmap_fn)

    rms_per_seg = dict()

    for longsegname in constants.SEGNAMES_WSS:
        segid, segnum = longsegname.split('-')

        # Calculate RMS per segment. Convert to nanometers
        segmask = segmap==int(segnum)
        rms_per_seg[segid] = utils.rms(opd, mask=segmask)*1e9

        if plot:
            plt.figure()
            plt.imshow(opd * segmask)
            plt.title(f"{segid}: {rms_per_seg[segid]:.2f} nm rms")

    return rms_per_seg
