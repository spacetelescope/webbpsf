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

from __future__ import division

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy
import astropy.table
import astropy.io.fits as fits
import astropy.units as u
import logging
import six

import poppy
import poppy.zernike as zernike
from . import constants
from . import utils

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
                 segment_mask_file='JWpupil_segments.fits',
                 **kwargs):
        """
        Parameters
        ----------
        opd : string, path to FITS file.
            FITS file to load an OPD from. The OPD must be specified in microns. This FITS file must be
            compatible with the format expected by poppy.FITSOpticalElement.
        transmission: str
            FITS file for pupil mask, with throughput from 0-1. If not explicitly provided, will be inferred from
            wherever is nonzero in the OPD file.


        ext : int, optional
            FITS extension to load OPD from
        slice : int, optional
            slice of a datacube to load OPD from, if the selected extension contains a datacube.

        """
        mypath = os.path.dirname(os.path.abspath(__file__)) + os.sep
        if opd is None and transmission is None:
            _log.debug('Neither a pupil mask nor OPD were specified. Using the default JWST pupil.')
            transmission = os.path.join(utils.get_webbpsf_data_path(), "jwst_pupil_revW_npix1024.fits.gz")

        super(OPD, self).__init__(name='Modified OPD',
                                  opd=opd, transmission=transmission,
                                  opd_index=opd_index, transmission_index=0,
                                  planetype=poppy.poppy_core.PlaneType.pupil, **kwargs)

        if self.opd_header is None:
            self.opd_header = self.amplitude_header.copy()

        self.segnames = np.asarray([a[0:2] for a in constants.SEGNAMES_WSS_ORDER])

        full_seg_mask_file = os.path.join(utils.get_webbpsf_data_path(), segment_mask_file)
        self._segment_masks = fits.getdata(full_seg_mask_file)

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
        """ Return an OPD as a fits.HDUList object

        Parameters
        -----------
        include_pupil : bool
            Include the pupil mask as a FITS extension?
        """

        output = fits.HDUList([self._opdHDU])
        output[0].header.update('EXTNAME', 'OPD')
        if include_pupil:
            # puphdu= fits.ImageHDU(self._pupilHDU, name='PUPIL')
            self.amplitude_header.update('EXTNAME', 'PUPIL')
            output.append(fits.ImageHDU(self.amplitude, self.amplitude_header))
            # hdus.append(puphdu)
            # hdus[0].header.update("EXTEND",'True', 'File may contain extensions')
        return output

    def writeto(self, outname, clobber=True, **kwargs):
        """ Write OPD to a FITS file on disk """
        self.as_fits(**kwargs).writeto(outname, clobber=clobber)

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

        cmap = matplotlib.cm.get_cmap(poppy.conf.cmap_diverging)
        cmap.set_bad('0.3')

        plt.clf()
        plt.subplot(231)
        self.display(title='full wavefront', clear=False, colorbar=False, vmax=vmax)

        ps_pixel_size = 1. / sampling  # how many cycles per pixel
        trans = SFT.SFT3(self.data, max_cycles * 2, max_cycles * 2 * sampling)

        abstrans = np.abs(trans)

        extent = [-max_cycles, max_cycles, -max_cycles, max_cycles]

        plt.subplot(233)
        plt.imshow(abstrans, extent=extent)
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
                       cmap=cmap)  # vmax is in nm, but WFE is in microns, so convert
            plt.title(label + " spatial frequencies")
            rms = (np.sqrt((inverse.real[wgood] ** 2).mean()) * 1000)

            components.append(rms)
            plt.xlabel("%.3f nm RMS WFE" % rms)

        return np.asarray(components)

    def display_opd(self, ax=None, labelsegs=True, vmax=150., colorbar=True, clear=False, title=None, unit='nm',
                    cbpad=None, colorbar_orientation='vertical',
                    show_axes=False, show_rms=True,
                    cmap=None):
        """ Draw on screen the perturbed OPD

        Parameters
        -----------
        ax : matplotlib.Axes
            axes instance to display into.
        labelsegs : bool
            draw segment name labels on each segment? default True.

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

        if clear:
            if ax is not None:
                raise RuntimeError("'clear=True' is incompatible with passing in an Axes instance.")
            plt.clf()
        if cmap is None:
            cmap = matplotlib.cm.get_cmap(poppy.conf.cmap_diverging)
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

        plot = ax.imshow(self.opd * mask * scalefact, vmin=-vmax, vmax=vmax, cmap=cmap, extent=extent)

        _log.debug("Displaying OPD. Vmax is %f, data max is %f " % (vmax, self.opd.max()))

        if title is None:
            title = self.name
        ax.set_title(title)
        if show_rms:
            ax.set_xlabel("RMS WFE = %.1f nm" % self.rms())

        if labelsegs:
            for seg in self.segnames[0:18]:
                self.label_seg(seg, ax=ax, show_axes=show_axes)
        if colorbar:
            if cbpad is None:
                cbpad = 0.05 if colorbar_orientation == 'vertical' else 0.15
            cb = plt.colorbar(plot, ax=ax, pad=cbpad, orientation=colorbar_orientation)
            cb.set_label("WFE [%s]" % unit)
        else:
            cb = None
        plt.draw()

        return ax, cb

    def label_seg(self, segment, ax=None, show_axes=False, color='black'):
        # Y, X = np.indices(self.opd.shape)

        # base = {'A':0, 'B':6,'C':12}
        # iseg = np.where(self.segnames == segment)[0][0] + 1  # segment index from 1 - 18
        # iseg = base[segment.upper()[0]]+int(segment[1])
        # wseg = np.where(self._segment_masks == iseg)
        # pupilscale = self.opd_header['PUPLSCAL']
        cx, cy = self._seg_centers_m[segment]
        # cx = (np.mean([X[wseg].min(), X[wseg].max()])  -512) * pupilscale
        # cy = (np.mean([Y[wseg].min(), Y[wseg].max()])  -512) * pupilscale

        offset = 0.2 if show_axes else 0

        if ax is None: ax = plt.gca()
        label = ax.text(cx + offset, cy + offset, segment, color=color, horizontalalignment='center', verticalalignment='center')

        if show_axes:
            ax_arrow_len = .3
            # if not ('C' in segment ):
            if True:

                for i, color, label in zip([0, 1, 2], ['green', 'blue', 'red'], ['x', 'y', 'z']):
                    vec = np.matrix([0, 0, 0]).transpose()  # xyz order
                    vec[i] = 1
                    b = self._rot_matrix_local_to_global(segment) * vec
                    b = np.asarray(b).flatten()  # Inelegant but it works

                    ax.arrow(cx, cy, ax_arrow_len * b[0], ax_arrow_len * b[1], color=color,
                             # width=ax,
                             head_width=.050, head_length=.080)  # in units of mm

                    xoffset = 0.1 if i == 2 else 0
                    ax.text(cx + ax_arrow_len * b[0] * 1.5 + xoffset, cy + ax_arrow_len * b[1] * 1.5, label,
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
        cmap = matplotlib.cm.jet
        cmap.set_bad('0.5', alpha=0.0)

        for j in np.arange(nzerns) + 1:
            ax = fig.add_subplot(3, 4, j, frameon=False, xticks=[], yticks=[])

            # n, m = zernike.noll_indices(j)
            Z = zernike.zernike1(j, npix=npix)
            ax.imshow(Z * zerns[j - 1] * hexap * scalefact, vmin=-1 * vmax, vmax=vmax, cmap=cmap)
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

    def __init__(self, opd=None, opd_index=0, transmission=None, rm_ptt=False, zero=False):
        """
        Parameters
        ----------
        opdfile : str or fits.HDUList
            FITS file to load an OPD from. The OPD must be specified in microns.
        ext : int, optional
            FITS extension to load OPD from
        slice : int, optional
            slice of a datacube to load OPD from, if the selected extension contains a datacube.
        pupilfile : str
            FITS file for pupil mask, with throughput from 0-1. If not explicitly provided, will be inferred from
            wherever is nonzero in the OPD file.
        rm_ptt : bool
            Remove piston, tip and tilt from all segments if set. Default is False.
        zero : bool
            If set, reate an OPD which is precisely zero in all locations. Default is False.


        """

        OPD.__init__(self, opd=opd, opd_index=opd_index, transmission=transmission)

        mypath = os.path.dirname(os.path.abspath(__file__)) + os.sep
        self._sensitivities = astropy.table.Table.read(os.path.join(mypath, 'otelm', 'seg_sens.txt'), format='ascii', delimiter='\t')
        self.state = {}
        self.remove_piston_tip_tilt = rm_ptt

        self._opd_original = self.opd.copy()
        if zero: self.zero()

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
            plt.imshow(Xr * (self._segment_masks == iseg))
            plt.title("Local X_Control for " + segment)
            plt.subplot(122)
            plt.imshow(Yr * (self._segment_masks == iseg))
            plt.title("Local Y_Control for" + segment)
            plt.draw()

        if self.remove_piston_tip_tilt:
            zernike_coeffs[0:3] = 0
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

            if display: self.display()

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
        if unit.endswith('s'): unit = unit[:-1]
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
        if unit.endswith('s'): unit = unit[:-1]
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

        if display: self.display()

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

    def __init__(self, name='Unnamed OPD', opd=None, opd_index=0, transmission=None, segment_mask_file='JWpupil_segments.fits',
                 zero=False, jsc=False, rm_ptt=False):
        """
        Parameters
        ----------
        opdfile : str or fits.HDUList
            FITS file to load an OPD from. The OPD must be specified in microns.
        ext : int, optional
            FITS extension to load OPD from
        slice : int, optional
            slice of a datacube to load OPD from, if the selected extension contains a datacube.
        pupilfile : str
            FITS file for pupil mask, with throughput from 0-1. If not explicitly provided, will be inferred from
            wherever is nonzero in the OPD file.
        jsc : bool
            Enable JSC OTIS test specific options?
        rm_ptt : bool
            Remove piston, tip, and tilt? This is mostly for visualizing the higher order parts of
            the LOM.

        """

        OPD.__init__(self, name=name, opd=opd, opd_index=opd_index, transmission=transmission, segment_mask_file=segment_mask_file)

        mypath = os.path.dirname(os.path.abspath(__file__)) + os.sep
        # load influence function table:
        self._influence_fns = astropy.table.Table.read(os.path.join(mypath, 'otelm', 'JWST_influence_functions_control_with_sm.fits'))
        self._control_modes = ['Xtilt', 'Ytilt', 'Piston', 'Clocking', 'Radial', 'ROC']
        self._sm_control_modes = ['Xtilt', 'Ytilt', 'Xtrans', 'Ytrans', 'Piston']
        # controllable modes in WAS order; yes it's not an obvious ordering but that's the order of the
        # WAS influence function matrix for historical reasons.
        self.state = {}
        self.segment_state = np.zeros((19, 6), dtype=float)  # 18 segs, 6 controllable DOF each, plus SM
        self.segnames = np.asarray(list(self.segnames) + ['SM'])  # this model, unlike the above, knows about the SM.

        self._opd_original = self.opd.copy()  # make a separate copy
        self._jsc = jsc
        self.remove_piston_tip_tilt = rm_ptt
        self._global_zernike_coeffs = np.zeros(15)
        if self._jsc:
            self._jsc_acf_tilts = np.zeros((3, 2))  # only for JSC sims. Tilts in microradians.

            # helper diagram taped to WSS monitor says:
            # 'ACF1' = ABC4, 'ACF2' = ABC6, 'ACF3' = ABC2.

            self._jsc_acf_cens = np.array([[1.90478, 0.6781127],  # Segs ABC2  V2, V3
                                           [-0.382703, -1.96286],  # Secs ABC4
                                           [-1.52316, 1.342146]])  # Segs ABC6
            self._jsc_acf_centers_pixels = self.shape[0] / 2 + self._jsc_acf_cens / self.pixelscale.value
        if zero: self.zero()

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
        if self._jsc:
            print("JSC Autocollimating flat tilts: ")
            print("  \t %10s %10s " % ("Xtilt", "Ytilt"))
            for i in range(3):
                print("  \t %10s %10s " % tuple(self._jsc_acf_tilts[i]))

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
            if table[i]['control_mode'] != self._control_modes[i]: raise RuntimeError("Influence function table has unexpected ordering")
            for h in range(nhexike):
                coeffs[i, h] = table[i]['Hexike_{}'.format(h)]

        if self._jsc:
            coeffs *= 2  # double pass on the SM for PAAH config

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
            if table[i]['control_mode'] != self._sm_control_modes[i]: raise RuntimeError("Influence function table has unexpected ordering")
            for h in range(nhexike):
                coeffs[i, h] = table[i]['Hexike_{}'.format(h)]

        if self._jsc:
            coeffs *= 2  # double pass on the SM for PAAH config

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
            Zernike coefficients, in units of microns
        """
        assert (segment in self.segnames)

        iseg = np.where(self.segnames == segment)[0][0] + 1  # segment index from 1 - 18
        wseg = np.where(self._segment_masks == iseg)

        # Note that the influence function matrix values already take into
        # account the rotations between segment coordinate systems.
        # so here we can just work in unrotated coordinates.

        # determine the X and Y hexike coordinates for each segment
        # determine the center of each segment, as the mean of the min and max X and Y values
        # FIXME this should just use the BATC-provided center coordinates; see jwst_ote3d.py

        Y, X = np.indices(self.opd.shape, dtype=float)
        cx, cy = self._seg_centers_pixels[segment]
        # cx = np.mean([X[wseg].min(), X[wseg].max()])
        # cy = np.mean([Y[wseg].min(), Y[wseg].max()])

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

        # returns a list of hexikes each with the same shape as Xc
        if self.remove_piston_tip_tilt:
            hexike_coeffs[0:3] = 0

        for i in range(len(hexike_coeffs)):
            if hexike_coeffs[i] == 0:
                continue
            self.opd[wseg] += hexikes[i] * hexike_coeffs[i]

        # outtxt="Hs=["+", ".join(['%.1e'%z for z in hexike_coeffs])+"]"
        # _log.debug("     "+outtxt)

    def _apply_acf_tilt(self, acfnum):
        """ Apply tilt moves for one of the JSC autocollimating flats

        """

        coeffs = self._jsc_acf_tilts[acfnum]  # Xtilt, Ytilt

        Y, X = np.indices(self.opd.shape, dtype=float)
        cx, cy = self._jsc_acf_centers_pixels[acfnum]
        acf_radius = 1.52 / 2 / self.pixelscale.value
        Y = (Y - cy) / acf_radius
        X = (X - cx) / acf_radius
        R = np.sqrt(X ** 2 + Y ** 2)
        wacf = np.where(R <= 1)
        Xc = X[wacf]
        Yc = Y[wacf]

        zern_xtilt = Yc * 2e-6  # remember, "Xtilt" means tilt around the X axis
        zern_ytilt = Xc * 2e-6  # Times 1e-6 to convert from microradians of tilt to meters of WFE
        # Times 2 since double pass, negative since facing the other way

        self.opd[wacf] += coeffs[0] * zern_xtilt + coeffs[1] * zern_ytilt

    def _apply_global_zernikes(self):
        """ Apply Zernike perturbations to the whole primary

        """

        if not self.opd.shape == (1024, 1024):
            raise NotImplementedError("Code need to be generalized for OPD sizes other than 1024**2")

        perturbation = poppy.opd_from_zernikes(self._global_zernike_coeffs,
                                               npix=1024, basis=poppy.zernike.zernike_basis_faster)

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

        if segment == 'SM':
            raise ValueError("SM not supported by move_seg_local. Use move_sm_local instead.")

        # Handle tilts and clocking
        tilts = np.array([xtilt, ytilt, clocking], dtype=float)
        if np.abs(tilts).sum() > 0:
            self.opd_header.add_history('Rotation: %s %s' % (str(tuple(tilts)), rot_unit))

            # sensitivity matrices are in microns per microradian
            # so convert all to urad.
            if rot_unit.endswith('s'): rot_unit = rot_unit[:-1]
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
            if trans_unit.endswith('s'): trans_unit = trans_unit[:-1]
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
            if rot_unit.endswith('s'): rot_unit = rot_unit[:-1]
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
            if trans_unit.endswith('s'): trans_unit = trans_unit[:-1]
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
            if rot_unit.endswith('s'): rot_unit = rot_unit[:-1]
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

            if trans_unit.endswith('s'): trans_unit = trans_unit[:-1]
            trans_unit = trans_unit.lower()
            if trans_unit == 'micron' or trans_unit == 'um':
                pass
            elif trans_unit == 'millimeters' or trans_unit == 'millimeter' or trans_unit == 'mm':
                vector *= 1000
            elif trans_unit == 'nm' or trans_unit == 'nanometer' or trans_unit == 'nanometers':
                vector /= 1000
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

    def move_jsc_acf(self, acfnum, xtilt=0.0, ytilt=0.0, unit='urad', absolute=False,
                     delay_update=False, display=False):
        """ Move autocollimating flats at JSC.
            NOTE - THIS IS ONLY APPLICABLE TO JSC OTIS CRYO - NOT FLIGHT!

        """
        if not self._jsc:
            raise RuntimeError("This instance of the linear model is not configured for JSC.")
        if absolute:
            self._jsc_acf_tilts[acfnum] = [xtilt, ytilt]
        else:
            self._jsc_acf_tilts[acfnum] += [xtilt, ytilt]

        if not delay_update:
            self.update_opd(display=display)

    def move_global_zernikes(self, zvector, unit='micron',
                             absolute=False):
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
        if unit.endswith('s'): unit = unit[:-1]
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

    def move_sur(self, sur_file, group=None, verbose=False):
        """ Move using a JWST Segment Update Request file

        """
        import jwxml

        sur = jwxml.SUR(sur_file)
        for grp in sur.groups:
            for update in grp:
                if verbose:
                    print("Move seg {} by {}".format(update.segment, str(update)))
                if update.type == 'pose':
                    if update.coord != 'local':
                        raise NotImplementedError("Only local moves supported!")

                    # FIXME - consider whether we should check for
                    # heterogenous sets of units here...
                    rot_unit = update.units['X_TILT']
                    trans_unit = update.units['X_TRANS']

                    self.move_seg_local(update.segment[0:2],
                                        xtilt=update.moves['X_TILT'],
                                        ytilt=update.moves['Y_TILT'],
                                        xtrans=update.moves['X_TRANS'],
                                        ytrans=update.moves['Y_TRANS'],
                                        piston=update.moves['PISTON'],
                                        absolute=update.absolute,
                                        rot_unit=rot_unit,
                                        trans_unit=trans_unit,
                                        delay_update=True)

                elif update.type == 'roc':
                    self.move_seg_local(update.segment[0:2],
                                        roc=update.moves['ROC'],
                                        absolute=update.absolute,
                                        delay_update=True)

                else:
                    raise NotImplementedError("Only moves of type='pose' or 'roc' are supported.")
        self.update_opd()

    def update_opd(self, display=False, verbose=False):
        """ Update the OPD based on the current linear model values.

        Users typically only need to call this directly if they have set the
        "delay_update" parameter to True in some function call to move mirrors.

        """

        # always start from the input OPD, then apply perturbations
        self.opd = self._opd_original.copy()

        sm = 18
        sm_pose_coeffs = self.segment_state[sm].copy()[0:5]  # 6th row is n/a for SM
        sm_pose_coeffs.shape = (5, 1)  # to allow broadcasting below

        for iseg, segname in enumerate(self.segnames[0:18]):
            pose_coeffs = self.segment_state[iseg].copy()
            if np.all(pose_coeffs == 0) and np.all(sm_pose_coeffs == 0):
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
                hexike_coeffs_combined = hexike_coeffs + hexike_coeffs_from_sm

                if verbose:
                    print("Need to move segment {} by {} ".format(segname, pose_coeffs.flatten()))
                    print("plus SM moved by {} ".format(sm_pose_coeffs.flatten()))
                    print("   Hexike coeffs: {}".format(hexike_coeffs))

                self._apply_hexikes_to_seg(segname, hexike_coeffs_combined)

        # Apply Global Zernikes
        if not np.all(self._global_zernike_coeffs == 0):
            self._apply_global_zernikes()

        # Apply NASA JSC OTIS test ACF tilts (not relevant in flight)
        if self._jsc and np.any(self._jsc_acf_tilts != 0):
            for iacf in range(3):
                self._apply_acf_tilt(iacf)
                if verbose:
                    print("Tilted JSC ACF {} by {}".format(iacf, self._jsc_acf_tilts[iacf]))

        if display:
            self.display()


################################################################################


def enable_adjustable_ote(instr, jsc=False):
    """
    Set up a WebbPSF instrument instance to have a modifiable OTE
    wavefront error OPD via an OTE linear optical model (LOM).

    Paramters
    ---------
    inst : WebbPSF Instrument instance
        an instance of one of the WebbPSF instrument classes.
    jsc : bool
        Use ACF pupil for JSC pass and a half test configuration

    Returns
    --------
    a modified copy of that instrumet set up to use the LOM, and
    the associated instance of the LOM.

    """
    import copy
    instcopy = copy.copy(instr)
    if instr.pupilopd is None:
        opdpath = None
    elif isinstance(instr.pupilopd, fits.HDUList):
        opdpath = instr.pupilopd
    else:
        # assume it is a string and try to use as filename
        # either an absolute or relative directory path if that works,
        # or else infer that it's a filename in the WebbPSF data directory.
        if not os.path.exists(instr.pupilopd):
            opdpath = os.path.join(instr._datapath, 'OPD', instr.pupilopd)
        else:
            opdpath = instr.pupilopd

    if jsc:
        pupilpath = os.path.join(utils.get_webbpsf_data_path(), "jwst_pupil_JSC_OTIS_Cryo.fits")
    else:
        pupilpath = instr.pupil

    name = "Modified OPD from " + str(instr.pupilopd)
    opd = OTE_Linear_Model_WSS(name=name,
                               opd=opdpath, transmission=pupilpath, jsc=jsc)

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
        'large' for the standard sizes used in OTE commissioning.
    guide_seg : string
        relevant mostly for coarse MIMF and image stacking. Kick out a segment to guide?
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

    jsc = ((size == 'jsc') or (size == 'jsc_compact') or (size == 'jsc_inverted'))
    if size is not None:
        if not jsc:
            nircam_pixelscale = 0.0311
            standard_sizes = {'small': 80 * nircam_pixelscale,
                              'medium': 300 * nircam_pixelscale,
                              'large': 812 * nircam_pixelscale
                              }
            radius = standard_sizes[size]

    # how many microradians of segment tilt per arcsecond of PSF motion?
    # note factor of 2 since reflection
    arcsec_urad = (1 * u.arcsec).to(u.urad).value / 2

    if reset:
        ote.reset()

    if not jsc:
        # Image Arrays used in flight
        size = radius * -1 if inverted else radius
        for i in range(1, 7):
            ote.move_seg_local('A' + str(i), xtilt=-size * arcsec_urad / 2, delay_update=True)
            ote.move_seg_local('B' + str(i), xtilt=size * arcsec_urad, delay_update=True)
            ote.move_seg_local('C' + str(i), ytilt=-size * np.sqrt(3) / 2 * arcsec_urad, delay_update=True)
    else:
        if not acfs_only:
            # Image Arrays used for JSC OTIS Cryo
            # 6 umicradian tilt of each of ABC 2,4,6
            # plus 50 microradian tilt of the ACFs
            # Standard tilts used at JSC are as follows. See BATC SER 2508696 by K. Smith and L. Coyle.
            xt = -5.1961524228
            yt = 3

            if size == 'jsc_inverted':
                xt *= -1
                yt *= -1
            for i in [2, 4, 6]:
                ote.move_seg_local('A' + str(i), xtilt=xt, ytilt=yt, delay_update=True)
                ote.move_seg_local('B' + str(i), xtilt=xt, ytilt=-yt, delay_update=True)
                ote.move_seg_local('C' + str(i), xtilt=-xt, ytilt=yt, delay_update=True)

        # ACF tilts. Also see BATC SER 2508696

        if size == 'jsc':  # regular "radial array"
            # acftilts = [[  6.868, -20.901],
            # [-21.647,   3.926],
            # [ 14.228,  16.780] ]
            acftilts = [  # rV2       rV3
                [6.868, -20.901],  # ABC2
                [-21.647, 3.926],  # ABC4
                [14.228, 16.780]]  # ABC6
        elif size == 'jsc_inverted':  # inverted "radial array"
            acftilts = [  # rV2       rV3
                [-14.228, -16.780],  # ABC2
                [21.647, -3.926],  # ABC4  CORRECT
                [-6.868, 20.901]]  # ABC6
        elif size == 'jsc_compact':
            acftilts = [[20.901, 6.868],  # ACF1 = ABC2
                        [-3.926, -21.647],  # ACF2 = ABC4
                        [-16.780, 14.228]]  # ACF3 = ABC6
            # acftilts = [[-16.780,  14.228],      # ACF1 = ABC4
            # [ -3.926, -21.647],      # ACF2 = ABC6
            # [ 20.901,   6.868] ]     # ACF3 = ABC2

        else:
            raise ValueError("Unknown array configuration.")

        for i in range(3):
            ote.move_jsc_acf(i, xtilt=acftilts[i][0], ytilt=acftilts[i][1], absolute=True, delay_update=True)

    if guide_seg is not None:
        # Undo the regular tilt for this segment, and then move it to
        # the side.
        if 'A' in guide_seg:
            xtilt = size * arcsec_urad / 2
            ytilt = guide_radius * arcsec_urad
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

def segment_primary(infile='JWpupil.fits'):
    """ Given a FITS file containing the JWST aperture, create an array
    with the PMSAs numbered from 1 - 18.

    This is used to create the 'JWpupil_segments.fits' file used by the OPDbender code.
    You should not need to run this routine unless you want to change the model of the
    JWST pupil.

    The basic algorithm is general, but right now this routine contains
    a bunch of hard-coded values tuned for the 1024x1024 Rev V pupil,
    so probably you will need to edit it for any other pupil.

    M. P. 2011-02-15


    """

    def renumber_array(array):
        uniq, inds = np.unique1d(array, return_inverse=True)
        # res = np.zeros_like(res2.shape)
        newuniq = range(len(uniq))
        res = np.array([newuniq[i] for i in inds]).astype(np.int16)
        res.shape = array.shape
        return res

    im = fits.getdata(infile)
    markers = np.zeros_like(im).astype(np.int16)
    xm, ym = np.ogrid[0:1024:102, 0:1024:100]
    markers[xm, ym] = np.arange(xm.size * ym.size).reshape((xm.size, ym.size))
    res2 = scipy.ndimage.watershed_ift(im.astype(np.uint8), markers)
    res2[xm, ym] = res2[xm - 1, ym - 1]  # remove the isolate seeds

    res = renumber_array(res2)

    # split and merge segments as needed - hard coded for Rev V pupil

    res3 = res.copy()  # np.zeros_like(res2.shape)
    merge_inds = [[0, 3], [2, 4], [11, 15], [18, 19]]
    for pair in merge_inds:
        res3[np.where(res == pair[1])] = pair[0]

    Y, X = np.indices(res.shape)
    w10h = np.where((res == 10) & (Y >= 512))
    res3[w10h] = 31
    w10h = np.where((res == 12) & (Y >= 512))
    res3[w10h] = 32

    res3 = renumber_array(res3)

    plt.clf()
    # plt.subplot(121)
    # plt.imshow(res)
    plt.subplot(121)
    plt.imshow(res3)

    for i in range(19):
        w = np.where(res3 == i)
        mx = X[w].mean()
        my = Y[w].mean()
        plt.text(mx, my, str(i), color='k')

    segs = ['A' + str(i + 1) for i in range(6)]
    for i in range(6): segs.append('B' + str(i + 1))
    for i in range(6): segs.append('C' + str(i + 1))
    seg_inds = [9, 18, 10, 4, 8, 17, 15, 14, 5, 1, 3, 11, 13, 7, 2, 0, 6, 12]

    result = np.zeros_like(res3)
    mxs = []
    mys = []
    for i in range(18):
        w = np.where(res3 == seg_inds[i])
        result[w] = i + 1
        mxs.append(X[w].mean())
        mys.append(Y[w].mean())

    plt.subplot(122)
    plt.imshow(result)
    for i in range(18): plt.text(mxs[i], mys[i], segs[i], color='k', horizontalalignment='center', verticalalignment='center')

    hdu = fits.PrimaryHDU((result * im).astype(np.uint8))
    for i in range(18):
        hdu.header.update('PMSA_' + str(i + 1), segs[i])

    # TODO copy relevant keywords and history from input FITS file header
    hdu.writeto("JWpupil_segments.fits", clobber=True)


def create_jsc_pupil(plot=False):
    """Create a pupil mask for the JSC Pass-and-a-half (PAAH)
    configuration with the 3 ACFs"""
    import webbpsf
    # Infer properties of the default pupil used with WebbPSF
    nc = webbpsf.NIRCam()
    defaultpupil = fits.open(nc.pupil)
    defpupilsize = defaultpupil[0].header['PUPLDIAM']

    # Create masks for the ACFs
    rad = 0.76  # Diam 1.52 per Randal Telfer

    acfs = []

    xcs = [-.382703, -1.52316, 1.90478]
    ycs = [-1.96286, 1.342146, 0.6781127]
    for x1, y1 in zip(xcs, ycs):
        acf = poppy.CircularAperture(radius=rad, shift_x=x1, shift_y=y1)
        acfs.append(acf)

    acfs = poppy.CompoundAnalyticOptic(acfs, mergemode='or')
    acfs.pupil_diam = defpupilsize

    # Create a FITS Optical Element instance with this mask applied:
    acfpupil = poppy.FITSOpticalElement(transmission=nc.pupil)
    acfmask = acfs.sample(npix=1024)
    acfpupil.amplitude *= acfmask
    acfpupil.name = "JWST Pupil for JSC PAAH"

    acfpupil.amplitude_header.add_history(" ")
    acfpupil.amplitude_header.add_history("**Modified for JSC Cryo Pass and a Half**")
    for x1, y1 in zip(xcs, ycs):
        acfpupil.amplitude_header.add_history("  Added JSC ACF: r={:.2f}, center=({:.4f},{:.4f}) m V2/V3".format(
            rad, x1, y1))
    acfpupil.amplitude_header.add_history("   Coords from Code V models via R. Telfer")
    acfpupil.amplitude_header['CONTENTS'] = "JWST Pupil for JSC OTIS Cryo"
    if plot:
        plt.figure(figsize=(10, 10))
        plt.imshow(acfs.sample(npix=1024) + defaultpupil[0].data)

    return acfpupil


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
