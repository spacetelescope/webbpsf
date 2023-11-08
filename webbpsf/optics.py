import os
import poppy
import poppy.utils
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from astropy.table import Table
import astropy.io.fits as fits
import astropy.units as units

from scipy.interpolate import griddata, RegularGridInterpolator
from scipy.ndimage import rotate

from . import utils
from . import constants

import logging

_log = logging.getLogger('webbpsf')

import pysiaf


#######  Classes for modeling aspects of JWST's segmented active primary #####

def segment_zernike_basis(segnum=1, nterms=15, npix=512, outside=np.nan):
    """ Basis set in the style of poppy.zernike.zernike_basis for segment-level
    Zernike polynomials for one segment at a time in JWST's aperture.

    Parameters
    ------------
    segnum : integer
        1 to 18, number of JWST segment. Uses same numbering convention as the WSS.
    nterms : integer
        Number of Zernike polynomial terms to return
    npix : integer
        Number of pixels per side of the array
    outside : float
        Value to fill the array with outside of the valid segment.

    """
    from .webbpsf_core import segname

    aper = WebbPrimaryAperture(label_segments=True)
    w = poppy.Wavefront(
        npix=npix,
        diam=constants.JWST_CIRCUMSCRIBED_DIAMETER
    )
    segmask = aper.get_transmission(w)

    segname = segname(segnum)
    cenx, ceny = aper.seg_centers[segname]

    # nominal point to point diam for A and B segments;
    # ignoring slight departures from ideal hexes for now.
    seg_radius = constants.JWST_SEGMENT_RADIUS

    y, x = w.coordinates()

    r = np.sqrt((y - ceny) ** 2 + (x - cenx) ** 2) / seg_radius
    theta = np.arctan2((y - ceny) / seg_radius, (x - cenx) / seg_radius)
    r[segmask != segnum] = np.nan
    theta[segmask != segnum] = np.nan

    wg = np.where(segmask == segnum)
    outzerns = np.full((nterms, npix, npix), outside, dtype=float)
    outzerns_tmp = poppy.zernike.zernike_basis(
        nterms=nterms,
        rho=r[wg],
        theta=theta[wg],
        outside=outside
    )
    for iz in range(nterms):
        outzerns[iz][wg] = outzerns_tmp[iz]

    return outzerns


class WebbPrimaryAperture(poppy.AnalyticOpticalElement):
    """ The JWST telescope primary mirror geometry, in all its
    hexagonal obscured complexity. Note this has **just the aperture shape**
    and not any wavefront error terms.

    JWST design pupil geometry and segment coordinates taken
    from Paul Lightsey's spreadsheet: "2010.03.16 Transmission X Area Budget.xls".
    That document was in turn based on Ball Aerospace drawing 2220169 Rev B,
    and the OTE Cryogenic Optics ICD, BATC doc # C327693.


    This class has no wavefront errors, it's just the pupil geometry including
    the segments (which are not quite perfect hexagons for manufacturing reasons
    related to trying to tile a curved surface with hexagons while maintaining
    uniform clearance between adjacent segments) and the secondary mirror support
    struts, including the bumps on the +V3 strut for the mid boom hinge and mag dampers.


    .. warning:: At high sampling factors, PSF calculations become a LOT slower.

    By default, this produces an aperture with values 0 and 1 for the transmission.
    By setting the parameter label_segments=True, you can instead have it generate a map of
    which segment number is in which location.
    """

    def __init__(self, name="WebbPrimaryAperture", label_segments=False, **kwargs):
        super(WebbPrimaryAperture, self).__init__(name=name, **kwargs)
        self.label_segments = label_segments
        self.segdata = constants.JWST_PRIMARY_SEGMENTS
        self.strutdata = constants.JWST_PRIMARY_STRUTS
        self.seg_centers = dict(constants.JWST_PRIMARY_SEGMENT_CENTERS)

    def get_transmission(self, wave):

        segpaths = {}
        strutpaths = []
        for segname, vertices in self.segdata:
            segpaths[segname] = matplotlib.path.Path(vertices)
        for strutname, vertices in self.strutdata:
            strutpaths.append(matplotlib.path.Path(vertices))

        y, x = wave.coordinates()
        pts = np.asarray([a for a in zip(x.flat, y.flat)])
        npix = wave.shape[0]
        out = np.zeros((npix, npix))

        # paint the segments 1 but leave out the SMSS struts
        for segname, p in segpaths.items():
            res = p.contains_points(pts)
            res.shape = (npix, npix)
            out[res] = 1 if not self.label_segments else int(segname.split('-')[1])
        for p in strutpaths:
            res = p.contains_points(pts)
            res.shape = (npix, npix)
            out[res] = 0
        return out


# Note - the following is **NOT USED YET **
# This will be finished up and used in a subsequent release to
# apply the OTE field dependence. For now just the fixed per SI stuff
# is there.
class WebbOTEPupil(poppy.FITSOpticalElement):
    """The complex OTE pupil, including:
        1) the aperture geometry, based on the cryo ICD detailed coordinates
        2) high spatial frequency WFE from the as-built mirrors in Rev G optical model
        3) mid frequencies from Rev W optical budget
        4) low frequency field-dependent WFE from the Rev G optical model.

        Parameters
        -----------
        level : '
    """

    def __init__(self, instrument=None, level='requirements', opd_index=0, **kwargs):

        if instrument is not None:
            self.instrument = instrument
            self.instr_name = instrument.name
            self.tel_coords = instrument._tel_coords()
        else:
            self.instrument = None
            self.instr_name = "NIRCam"
            # TODO figure out default V2V3 coords here
            self.tel_coords = (0, 0)  # ? TODO

        # determine filename for pupil amplitude array
        aperture_file = 'jwst_pupil_revW_npix1024.fits.gz'
        aperture_file = os.path.abspath(os.path.join(
            utils.get_webbpsf_data_path(), aperture_file
        ))

        # determine filename for the OPD array
        #   This should contain a precomputed combination of
        #   Rev G high spatial frequencies and
        #   Rev W mid spatial frequencies
        # Depends on what the 'level' parameter is.

        if level == 'perfect':
            opd_file = os.path.join(
                utils.get_webbpsf_data_path(),
                'OPD_jwst_ote_perfectly_aligned.fits'
            )
        elif level in ('predicted', 'requirements'):
            opd_file = os.path.join(
                utils.get_webbpsf_data_path(),
                self.instr_name,
                'OPD',
                'OPD_RevW_ote_for_{}_{}.fits'.format(self.instr_name, level)
            )
        else:
            raise ValueError("Invalid/unknown wavefront error level")

        super(WebbOTEPupil, self).__init__(name='JWST Primary',
                                           transmission=aperture_file,
                                           opd=opd_file,
                                           **kwargs)

        if self.instrument is not None:
            # we need a field point to be able to use this so
            # just skip it if we don't have one.

            # determine Zernike coeffs for field dependent error
            # based on Rev G field dependence model.

            coeffs = np.zeros(22)
            self.zernike_coeffs = coeffs

            # TODO apply that to as a modification to the OPD array.


#######  Custom Optics used in JWInstrument classes  #####


class NIRSpec_three_MSA_shutters(poppy.AnalyticOpticalElement):
    """ Three NIRSpec MSA shutters, adjacent vertically."""

    def get_transmission(self, wave):
        """ Compute the transmission inside/outside of the field stop.

        The area of an open shutter is 0.2 x 0.45, while the shutter pitch is 0.26x0.51
        The walls separating adjacent shutters are 0.06 arcsec wide.
        """

        msa_width = 0.2
        msa_height = 0.45
        msa_wall = 0.06

        if not isinstance(wave, poppy.Wavefront):
            raise ValueError("get_transmission must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._IMAGE)

        y, x = wave.coordinates()

        self.transmission = np.zeros(wave.shape)
        # get the innermost shutter than spans the Y axis
        w_inside_1 = np.where((abs(y) < (msa_height / 2)) & (abs(x) < (msa_width / 2)))
        self.transmission[w_inside_1] = 1
        # get the adjacent shutters one above and one below.
        w_inside_2 = np.where((abs(y) > (msa_height / 2) + msa_wall) & (abs(y) < msa_height * 1.5 + msa_wall) & (abs(x) < (msa_width / 2)))
        self.transmission[w_inside_2] = 1

        return self.transmission


class NIRSpec_MSA_open_grid(poppy.AnalyticOpticalElement):
    """ An infinite repeating region of the NIRSpec MSA grid"""

    def get_transmission(self, wave):
        """ Compute the transmission inside/outside of the field stop.

        The area of an open shutter is 0.2 x 0.45, while the shutter pitch is 0.26x0.51
        The walls separating adjacent shutters are 0.06 arcsec wide.
        """

        msa_width = 0.2
        msa_height = 0.45
        msa_wall = 0.06
        msa_x_pitch = 0.26
        msa_y_pitch = 0.51

        if not isinstance(wave, poppy.Wavefront):
            raise ValueError("get_transmission must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._IMAGE)

        y, x = wave.coordinates()
        # xnew =  x*np.cos(np.deg2rad(self.angle)) + y*np.sin(np.deg2rad(self.angle))
        # ynew = -x*np.sin(np.deg2rad(self.angle)) + y*np.cos(np.deg2rad(self.angle))
        # x,y = xnew, ynew

        mask_vert_walls = np.abs(np.mod(np.abs(x), msa_x_pitch) - (msa_x_pitch / 2)) < msa_wall / 2
        mask_horz_walls = np.abs(np.mod(np.abs(y), msa_y_pitch) - (msa_y_pitch / 2)) < msa_wall / 2

        self.transmission = np.ones(wave.shape)
        self.transmission[mask_vert_walls] = 0
        self.transmission[mask_horz_walls] = 0

        return self.transmission


class NIRISS_GR700XD_Grism(poppy.AnalyticOpticalElement):
    """ Custom optic class to model the NIRISS SOSS grim GR700XD

    This includes both the pupil mask file and the cylindrical lens

    Based on inputs from Loic Albert, Anand Sivaramakrishnan, and Andre Martel
    In particular see FGS_TFI_UdM_035_RevD for details of the NIRISS GR700XD
    measurement, and JWST-STScI-003338 for detector orientation and layout.

    GRISM DESIGN:

    The grism (and cylinder) are per design rotated by 2 degrees so as to be able
    to sample an emission line across different pixel position along the spatial
    direction (kind of resampling the line and not be limited by intra pixel
    response).

    From Loic Albert's NIRISS technical report:

        * surface sag for the cylinder: 3.994 micron peak
        * limited to 3.968 microns for the 26 mm FOV mask

    From Loic Albert's email to Marshall 2013-07-18:

            I do have an update concerning the geometry of the GR700XD pupil
            mask. It turns out that they clocked the grism by about 2.25 degrees wrt the
            OTE system of coordinates. However, the field mask did not follow and is still
            aligned along the OTE s.o.c. That was a mistake that fortunately does have much
            impact.

            Comdev is in the process of modelling a new mask for the
            Spare grism. Remember that we will swap the current FLight GR700XD for
            its Spare which offers much improved performances. The mask size will
            be a little different (rectangular) and this time will be clocked 2.25
            degrees along with the grism.

            The sign of the rotation of the grism will have to be
            devised by trying the 2 possibilities an looking at the resulting tilt
            of the monochromatic PSF and the position of that PSF on the detector.
            Attached is a simulation of what we expect based on my own PSF
            generator.

            The cylinder lens has a well characterized power (actually radius of curvature). The values are:
                current Flight: 22.85 meters
                Spare: 22.39 meters

            Prism physical size: pupil is 26 mm on a side for the current prism, will be 28 mm for the spare

    From Loic Albert's email to Marshall 2013-09-19:

        The latest news on this front are:

        1 - The current Flight mask is attached. It is 26x26 mm. The mask and grism are
            *not* aligned along the same coordinates. That was a mistake. I'll forward you
            a message from Michael M., our optics expert at CSA.

        2 - The Spare mask (likely the one which will fly) is not built yet. The mask
            will be aligned along the grism coordinate and both will be clocked 2.2 deg wrt
            the OTE.

        3 - A ghost analysis showed that the current grism clocking will suffer from
            large ghosts. So we are studying how to clock the Spare grism in its cell to
            minimize ghosts. Likely a 90 degrees rotation will be applied to baseline of
            point 2.

        From Michael.Maszkiewicz@asc-csa.gc.ca:

            As far as I understand now, we have two rotations in the as-built
            GR700. One rotation is for the prism-grism combo by 2 deg CCW, looking along
            the local +z axis, and the second rotation is for  the mask by 3.05 deg  but
            CW. As a result there is total 5.05 deg rotation between grism and its mask.
            See my annotations to your drawing attached.

    From Loic Albert's email to Marshall 2014-05-20:

        I should have pointed that the power assumed in my simulations for the cylindrical lens
        was off. It was one of the conclusions of CV1RR. The actual radius of curvature of the
        cylinder is 25.3 meters (rather than the smaller figure I used before).

     ORIENTATION:

        See Figure 2 of JWST-STScI-003338
        In "DMS" coordinates, as projected looking outwards onto the sky,
        The GR700XD grating trace is near the extreme right edge of the detector
        with long wavelengths closest to (2048,2048) and short wavelengths nearest (2048,0)
        (The raw detector coordinates are very different from this due to a 180 degree rotation)

        **PLEASE NOTE** that the DMS when processing spectral data performs an additional transformation:
            For spectral data, the science X-axis is aligned with the detector
            dispersion direction and the science frame Y-axis is at a right angle
            to the X-axis in a right-handed coordinate system (Swade 2003)

        We choose here to ignore that complication; WebbPSF simulates the 2D sky projected
        image in "Sci" coordinates in the terminology for SIAF from Lallo et al.
        In this coordinate system, the dispersion from the cylinder lens is aligned
        almost along V2 and the longer wavelengths are oriented toward +V3.




    Parameters
    ----------
    which : string
        'initial' or 'spare'. Properties are hard coded.
    """

    #
    #    transmission : string filename
    #        file for the pupil transmission function
    #    cylinder_sag_mm : float
    #        physical thickness of the cylindrical lens, in millimeters
    #    rotation_angle : float
    #        degrees clockwise for the orientation of the cylinder's dispersing axis. Default
    #        of 92.25 should be consistent with initial NIRISS girsm and spare, except for
    #        sign ambiguity.
    #    rotate_mask : bool
    #        should the field mask be rotated along with the cylinder? False for first gen initial
    #        prism, true for expected spare replacement.

    def __init__(self, name='GR700XD', which='Bach',
                 # cylinder_radius=22.85,  cylinder_sag_mm=4.0, rotation_angle=92.25, rotate_mask=False, transmission=None,
                 **kwargs):
        # Initialize the base optical element with the pupil transmission and zero OPD

        if which == 'LLNL':
            raise NotImplementedError("Rotated field mask for LLNL grism not yet implemented!")
        elif which == 'Bach':
            transmission = os.path.join(utils.get_webbpsf_data_path(), "NIRISS/optics/MASKGR700XD.fits.gz")
        else:
            raise NotImplementedError("Unknown grating name:" + which)

        poppy.AnalyticOpticalElement.__init__(self, name=name, planetype=poppy.poppy_core._PUPIL, **kwargs)

        # UPDATED NUMBERS 2013-07:
        # See Document FGS_TFI_UdM_035_RevD

        _log.debug("Computing properties for {0} grism".format(which))
        if which == 'Bach':
            # ---- Phase properties ---------------
            # 3.994 microns P-V over 27.02 mm measured (Loic's email)
            # This is **surface sag**, corresponding to P-V of 6.311 waves at lambda=632.8 nm.
            # should correspond to 3.698 microns over 26 mm clear aperture.
            self.prism_size = 0.02702  # 27.02 millimeters for the physical prism
            self.prism_clear_aperture = 0.0260  # 26 mm clear aperture for the prism + mount
            self.cylinder_rotation_angle = 2  # was 2.25

            # self.cylinder_radius = 22.85 # radius of curvature  ; Nominal
            # but they discarded that and used 25.3 instead
            # From Lafreniere's wfe_cylindricallens.pro:
            #  "OVERRIDE PREVIOUS CASES AFTER CV1RR RESULTS:"
            self.cylinder_radius = 25.3  # radius of curvature

            # ---- Amplitude Transmission / Pupil shape ---------------
            self.pupil_size_mm = 26.0
            # Note that the IDL code says 26 mm is 683.75 pixels using the assumed demagnification
            self.pupil_rotation_angle = 2.0

        else:
            # 5.8 microns P-V over 32.15 mm (Loic's email)
            # should correspond to 4.38 microns over 28 mm clear aperture
            self.cylinder_radius = 22.39  # radius of curvature
            self.prism_size = 0.03215  # millimeters for the physical prism
            self.prism_clear_aperture = 0.0280  # clear aperture for the prism + mount
            self.cylinder_rotation_angle = 2.25

        # We need to know the magnification scale of the NIRISS reimaged pupil
        # in order to compute the curvature in the full-pupil units that POPPY uses
        # internally

        # pupil magnification computed from 22 mm clear aperture reported =
        # 857-169 pixels = 699 pixels in the 2D array which has scale =.00645604
        # = 4.44175 meters projected on the primary

        # 2014-05-21 but wait, that's actually 26 mm!
        # so the 699 pixels at 0.00645604 m/pixel = 4.512 meters implies the magnificationa 173 not 170
        # but, double wait, it's actually more like 687 pixels across rather than 699 so that makes it 170 again.

        # therefore the magnification is 0.1708 meters projected on the primary / mm in the NIRISS pupil
        # self.pupil_demagnification =  170.8367 # meters on the primary / meters in the NIRISS pupil
        # self.pupil_demagnification =  173.56 # meters on the primary / meters in the NIRISS pupil

        # Anand says:
        #  nominally the circumscribing circle at the PW of NIRISS is ~40mm.  I use 39mm for the nrm, but it's slightly field-dependent.  Compare that to the 6.6... PM circle?
        self.pupil_demagnification = 6.6 / 0.040  # about 165

        # perform an initial population of the OPD array for display etc.
        tmp = self.get_phasor(poppy.Wavefront(2e-6))

    def get_opd(self, wave):
        """ Make an OPD array corresponding to the cylindrical weak lens
        used for defocusing the spectrum in the perpendicular-to-dispersion direction.

        Parameters
        ----------
        wave : float or obj
            either a scalar wavelength (meters) or a Wavefront object
        """

        # wave should be a Wavefront object
        # wavelength is an astropy.units type
        if isinstance(wave, poppy.Wavefront):
            wavelength = wave.wavelength
        else:
            wave = poppy.Wavefront(wavelength=float(wave))
            wavelength = wave.wavelength

        # compute indices in pixels, relative to center of plane, with rotation
        # units of these are meters
        y, x = wave.coordinates()

        ang = np.deg2rad(self.cylinder_rotation_angle)
        x = np.cos(ang) * x - np.sin(ang) * y
        y = np.sin(ang) * x + np.cos(ang) * y

        _log.debug(" Rotating local grism axes by {0} degrees".format(self.cylinder_rotation_angle))

        # From IDL code by David Lafreniere:
        #  ;the cylindrical defocus
        # x=(dindgen(pupdim)-pupdim/2)#replicate(1,pupdim)
        # y0=(rpuppix^2+sag[s]^2)/(2*sag[s])
        # wfe1=y0-sqrt(y0^2-x^2)
        # if sag[s] lt 1.e-5 then wfe1=0.d0

        # Here I will just translate that to Python exactly, making use of the
        # variables here:

        # rpuppix = radius of pupil in pixels
        # rpuppix = self.amplitude_header['DIAM'] / self.amplitude_header['PUPLSCAL'] / 2
        # Calculate the radius of curvature of the cylinder, based on
        # the chord length and height

        # In this case we're assuming the cylinder is precisely as wide as the projected
        # telescope pupil. This doesn't seem guaranteed:
        #  * actual chord length across cylinder: 27.02 mm.
        #  * projected primary scale at NIRISS = ?

        _log.debug(" Computing GR700XD cylinder based on RoC: {0:.3g} meters".format(self.cylinder_radius))
        _log.debug(
            " Computing GR700XD cylinder based on pupil demagnification: {0:.3g} primary to grism".format(self.pupil_demagnification))

        # Compute the overall sag of the cylinder lens at its outer edge. This is not actually used, it's
        # just for cross-check of the values
        # the sag will depend on half the pupil size since that's the offset from center to edge
        sag0 = np.sqrt(self.cylinder_radius ** 2 - (self.prism_size / 2) ** 2) - self.cylinder_radius
        _log.debug(" Computed GR700XD cylinder sag at lens outer edge (for cross check only): {0:.3g} meters".format(sag0))

        # now compute the spatially dependent sag of the cylinder, as projected onto the primary

        # what is the pupil scale at the *reimaged pupil* of the grism?
        pupil_scale_m_per_pix = 38.0255e-6  # Based on UdeM info in wfe_cylindricallens.pro
        # sag = np.sqrt(self.cylinder_radius**2 - (x*self.amplitude_header['PUPLSCAL']/self.pupil_demagnification)**2) - self.cylinder_radius
        sag = np.sqrt(self.cylinder_radius ** 2 - (x / self.pupil_demagnification) ** 2) - self.cylinder_radius
        # sag = self.cylinder_radius -  np.sqrt(self.cylinder_radius**2 - (x * pupil_scale_m_per_pix )**2 )

        # what we really want to do is take the physical properties of the as-built optic, and interpolate into that
        # to compute the OPD after remapping based on the pupil scale (and distortion?)
        # y0=(rpuppix**2+self.cylinder_sag**2)/(2*self.cylinder_sag)
        # wfe1=y0-np.sqrt(y0**2-x**2)

        _log.debug(" Cylinder P-V: {0:.4g} meters physical sag across full array".format(sag.max() - sag.min()))

        # no OPD in opaque regions (makes no difference in propagation but improves display)
        if self._transmission.shape != sag.shape:
            tmp = self.get_transmission()  # Update the ._transmission attribute
        sag[self._transmission == 0] = 0
        wnz = np.where(self._transmission != 0)  # use this just for display of the log messages:
        _log.debug(" Cylinder P-V: {0:.4g} meters physical sag across clear aperture".format(sag[wnz].max() - sag[wnz].min()))

        # scale for index of refraction
        index = self.ZnS_index(wavelength)
        opd = sag * (index - 1)
        lambda_micron = wavelength.to(units.micron).value
        _log.debug(" Scaling for ZnS index of refraction {0} at {1:.3g} microns".format(index, lambda_micron))
        _log.debug(
            " Cylinder P-V: {0:.4g} meters optical sag at {1:.3g} microns across clear aperture".format(opd[wnz].max() - opd[wnz].min(),
                                                                                                        lambda_micron))
        return opd

    def get_transmission(self, wave):
        """ Make array for the pupil obscuration appropriate to the grism
        """

        if isinstance(wave, poppy.Wavefront):
            wavelength = wave.wavelength
        else:
            wave = poppy.Wavefront(wavelength=float(wave))
            wavelength = wave.wavelength

        y, x = wave.coordinates()
        ang = np.deg2rad(self.pupil_rotation_angle)
        x = np.cos(ang) * x - np.sin(ang) * y
        y = np.sin(ang) * x + np.cos(ang) * y

        _log.debug("Rotating local pupil mask axes by {0} degrees".format(self.cylinder_rotation_angle))

        pupil_halfsize_m = self.pupil_size_mm / 2 / 1000 * self.pupil_demagnification
        pupilmask = np.ones_like(x)
        pupilmask[np.abs(x) > pupil_halfsize_m] = 0
        pupilmask[np.abs(y) > pupil_halfsize_m] = 0

        self._transmission = pupilmask

        return pupilmask

    @poppy.utils.quantity_input(wavelength=units.meter)
    def ZnS_index(self, wavelength, temperature=40):
        """ Return cryogenic index of refraction of ZnS (Cleartran)

        Based on IDL function index_cleartran provided by Loic Albert at U de Montreal
        Which was in turn based on Leviton and Fray 2013
        http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1744938
        doi:10.1117/12.2024817

        """
        lambda_micron = wavelength.to(units.micron).value

        # Sellmeier dispersion model
        # From Leviton & Frey measurements (SPIE preprint) (assumes lambda in microns)
        S_1 = np.asarray([[3.35933, -5.12262e-4, 1.01086e-5, -4.14798e-8, 6.91051e-11]])
        S_2 = np.asarray([[0.706131, 4.89603e-4, -8.91159e-6, 3.81621e-8, -6.54805e-11]])
        S_3 = np.asarray([[4.02154, -2.93193e-2, 2.31080e-4, -7.57289e-07, 8.31188e-10]])
        S_ij = np.concatenate((S_1, S_2, S_3), axis=0)
        lambda_1 = np.array([[0.161151, -8.93057E-06, 2.73286E-07, -1.23408E-09, 2.29917E-12]])
        lambda_2 = np.array([[0.282427, -4.66636E-05, 7.55906E-07, -2.77513E-09, 4.35237E-12]])
        lambda_3 = np.array([[41.1590, -0.161010, 1.23906E-03, -3.95895E-06, 4.16370E-09]])
        lambda_ij = np.concatenate((lambda_1, lambda_2, lambda_3))

        n2minus1 = 0.0
        T = temperature
        for i in range(3):
            S_i = S_ij[i, 0] + S_ij[i, 1] * T + S_ij[i, 2] * T ** 2.0 + S_ij[i, 3] * T ** 3.0 + S_ij[i, 4] * T ** 4.0
            lambda_i = lambda_ij[i, 0] + lambda_ij[i, 1] * T + lambda_ij[i, 2] * T ** 2.0 + lambda_ij[i, 3] * T ** 3.0 + lambda_ij[
                i, 4] * T ** 4.0
            n2minus1 += S_i * lambda_micron ** 2.0 / (lambda_micron ** 2.0 - lambda_i ** 2.0)

        cleartran_index = np.sqrt(1.0 + n2minus1)
        return cleartran_index

    def display(self, opd_vmax=6e-6, *args, **kwargs):
        "Same as regular display for any other optical element, except opd_vmax default changed"
        poppy.AnalyticOpticalElement.display(self, *args, opd_vmax=opd_vmax, **kwargs)


class NIRISS_CLEARP(poppy.CompoundAnalyticOptic):
    """NIRISS 'CLEARP' pupil, including PAR obscuration

        **CAUTIONARY NOTE** TODO: This class represents this
        optic as having a circular outer edge; in reality the
        hardware has a 4% oversized tricontagon mask around the
        JWST pupil image. However as the primary mirror should
        serve as the pupil stop, in practice this model
        simplification should not affect output PSFs in imaging
        modes. This simplification may be removed in a future
        version of WebbPSF.
        See https://github.com/mperrin/webbpsf/issues/71


        CLEARP pupil info from:
           MODIFIED CALIBRATION OPTIC HOLDER - NIRISS
           DRAWING NO 196847  REV 0  COMDEV
           Design file name 196847Rev0.pdf sent by Loic Albert

        Properties:
          39 mm outer diam, corresponds to the circumscribing pupil of JWST
          2.0 mm vane width
          6.0 mm radius for central obstruction

        Note the circumscribing pupil of JWST is 6603.464 mm in diameter
        (Ball SER on geometric optics model: BALL-JWST-SYST-05-003)
        and therefore the NIRISS pupil magnification is 6.603464/39.0
        = 0.1693 meters (JWST primary) per mm (NIRISS internal pupil)

        Pupil distortions are not included in this model.

    """

    def __init__(self, *args, **kwargs):
        # CLEARP pupil info from:
        #   MODIFIED CALIBRATION OPTIC HOLDER - NIRISS
        #   DRAWING NO 196847  REV 0  COMDEV
        #   Design file name 196847Rev0.pdf sent by Loic Albert
        # Properties:
        #  39 mm outer diam, corresponds to the circumscribing pupil of JWST
        #  2.0 mm vane width
        #  6.0 mm radius for central obstruction
        # Note the circumscribing pupil of JWST is 6603.464 mm in diameter
        #  (Ball SER on geometric optics model: BALL-JWST-SYST-05-003)

        pupil_mag = 6.603464 / 39.0
        poppy.CompoundAnalyticOptic.__init__(self, (
            poppy.SecondaryObscuration(secondary_radius=6.0 * pupil_mag,
                                       support_width=2.0 * pupil_mag,
                                       n_supports=3,
                                       support_angle_offset=90 + 180,  # align first support with +V2 axis
                                       # but invert to match OTE exit pupil
                                       *args, **kwargs),
            poppy.CircularAperture(radius=39 * pupil_mag / 2,
                                   *args, **kwargs)),
                                             name='CLEARP')


class NIRCam_BandLimitedCoron(poppy.BandLimitedCoron):
    """ Band Limited Coronagraph

    Parameters
    ----------
    name : string
        Descriptive name. Must be one of the defined NIRCam coronagraphic mask names.
    module : string
        A or B
    nd_squares : bool
        Include the ND squares in the mask simulation? (Not an option in the real instrument;
        solely for certain simulation checks.)
    bar_offset : float
        Offset along coronagraphic bar (wedge) occulter, in arcseconds.
        Used for computing a PSF at a different position along the wedge, while
        keeping the convention that the target star has zero tip/tilt.
        This option is used to MANUALLY specify a specific position along the bar;
        see also the following option auto_offset.
    auto_offset : string or None
        Set to a NIRCam filter name to automatically offset to the nominal
        position along the bar for that filter. See bar_offset if you want to set
        to some arbitrary position.
    shift_x, shift_y : floats or None
        X and Y offset shifts applied to the occulter, via the standard mechanism for
        poppy.AnalyticOpticalElements. Like bar_offset but allows for 2D offsets, and
        applies to both bar and wedge coronagraphs.  This is IN ADDITION TO any offset
        from bar_offset.
    """
    allowable_kinds = ['nircamcircular', 'nircamwedge']
    """ Allowable types of BLC supported by this class"""

    def __init__(self, name="unnamed BLC", kind='nircamcircular', module='A', nd_squares=True,
                 bar_offset=None, auto_offset=None, **kwargs):
        super(NIRCam_BandLimitedCoron, self).__init__(name=name, kind=kind, **kwargs)
        if module not in ['A', 'B']:
            raise ValueError("module parameter must be 'A' or 'B'.")
        self.module = module
        self.nd_squares = nd_squares

        if self.name == 'MASK210R':
            self.sigma = 5.253
            self.kind = 'nircamcircular'
        elif self.name == 'MASK335R':
            self.sigma = 3.2927866
            self.kind = 'nircamcircular'
        elif self.name == 'MASK430R':
            self.sigma = 2.58832
            self.kind = 'nircamcircular'
        elif self.name == 'MASKSWB':
            self.kind = 'nircamwedge'
            # coeffs set in lookup table inside getPhasor
        elif self.name == 'MASKLWB':
            self.kind = 'nircamwedge'
            # coeffs set in lookup table inside getPhasor
        else:
            raise NotImplementedError("invalid name for NIRCam occulter: " + self.name)

        # EDIT: updated on 8 Dec 2021 to grab offsets directly from pySIAF
        self.siaf = pysiaf.Siaf('NIRCAM')
        self.offset_swb = {filt: self.get_bar_offset_from_siaf(filt, channel='SW')
                           for filt in ["F182M", "F187N", "F210M", "F212N", "F200W", 'narrow']}
        self.offset_lwb = {filt: self.get_bar_offset_from_siaf(filt, channel='LW')
                           for filt in ["F250M", "F300M", "F277W", "F335M", "F360M", "F356W", "F410M", "F430M", "F460M", "F480M", "F444W", 'narrow']}

        if bar_offset is None and auto_offset is not None:
            offsets = self.offset_swb if self.name.lower() == 'maskswb' else self.offset_lwb
            try:
                bar_offset = offsets[auto_offset]
                _log.debug("Set bar offset to {} based on requested filter {} on {}.".format(bar_offset, auto_offset, self.name))
            except:
                raise ValueError("Filter {} does not have a defined nominal offset position along {}".format(auto_offset, self.name))

        if bar_offset is not None:
            if self.kind == 'nircamcircular':
                raise ValueError("bar_offset option only makes sense with the bar occulters.")
            self.bar_offset = float(bar_offset)
            _log.debug("Set offset along {} to {} arcsec.".format(self.name, self.bar_offset))
        else:
            self.bar_offset = None

    def get_bar_offset_from_siaf(self, filt, channel='LW'):
        """ Get bar offset directly from SIAF.
        
        """
        
        if channel == 'SW':
            refapername = 'NRCA4_MASKSWB'
            apername = 'NRCA4_MASKSWB_' + filt.upper()
        else: # otherwise default to LW
            refapername = 'NRCA5_MASKLWB'
            apername = 'NRCA5_MASKLWB_' + filt.upper()
        offset_arcsec = np.sqrt((self.siaf.apertures[refapername].V2Ref - self.siaf.apertures[apername].V2Ref)**2 + (self.siaf.apertures[refapername].V3Ref - self.siaf.apertures[apername].V3Ref)**2)
        sign = np.sign(self.siaf.apertures[refapername].V2Ref - self.siaf.apertures[apername].V2Ref)
        
        return sign * offset_arcsec

    def get_transmission(self, wave):
        """ Compute the amplitude transmission appropriate for a BLC for some given pixel spacing
        corresponding to the supplied Wavefront.

        Based on the Krist et al. SPIE paper on NIRCam coronagraph design

        Note that the equations in Krist et al specify the intensity transmission of the occulter,
        but what we want to return here is the amplitude transmittance. That is the square root
        of the intensity, of course, so the equations as implemented here all differ from those
        written in Krist's SPIE paper by lacking an exponential factor of 2. Thanks to John Krist
        for pointing this out.

        """
        import scipy.special
        if not isinstance(wave, poppy.Wavefront):  # pragma: no cover
            raise ValueError("BLC getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == poppy.poppy_core._IMAGE)

        y, x = self.get_coordinates(wave)

        if self.bar_offset is not None:
            x += float(self.bar_offset)

        if self.kind == 'nircamcircular':
            r = poppy.accel_math._r(x, y)
            sigmar = self.sigma * r

            # clip sigma: The minimum is to avoid divide by zero
            #             the maximum truncates after the first sidelobe to match the hardware
            bessel_j1_zero2 = scipy.special.jn_zeros(1, 2)[1]
            sigmar.clip(np.finfo(sigmar.dtype).tiny, bessel_j1_zero2, out=sigmar)  # avoid divide by zero -> NaNs
            if poppy.accel_math._USE_NUMEXPR:
                import numexpr as ne
                jn1 = scipy.special.j1(sigmar)
                self.transmission = ne.evaluate("(1 - (2 * jn1 / sigmar) ** 2)")
            else:
                self.transmission = (1 - (2 * scipy.special.j1(sigmar) / sigmar) ** 2)
            self.transmission[r == 0] = 0  # special case center point (value based on L'Hopital's rule)

        elif self.kind == 'nircamwedge':
            # This is hard-coded to the wedge-plus-flat-regions shape for NIRCAM

            # the scale fact should depend on X coord in arcsec, scaling across a 20 arcsec FOV.
            # map flat regions to 2.5 arcsec each
            # map -7.5 to 2, +7.5 to 6. slope is 4/15, offset is +9.5
            wedgesign = 1 if self.name == 'MASKSWB' else -1  # wide ends opposite for SW and LW

            scalefact = (2 + (x * wedgesign + 7.5) * 4 / 15).clip(2, 6)

            # Working out the sigma parameter vs. wavelength to get that wedge pattern is non trivial
            # This is NOT a linear relationship. See calc_blc_wedge helper fn below.

            if self.name == 'MASKSWB':  # np.abs(self.wavelength - 2.1e-6) < 0.1e-6:
                polyfitcoeffs = np.array([2.01210737e-04, -7.18758337e-03, 1.12381516e-01,
                                          -1.00877701e+00, 5.72538509e+00, -2.12943497e+01,
                                          5.18745152e+01, -7.97815606e+01, 7.02728734e+01])
            elif self.name == 'MASKLWB':  # elif np.abs(self.wavelength - 4.6e-6) < 0.1e-6:
                polyfitcoeffs = np.array([9.16195583e-05, -3.27354831e-03, 5.11960734e-02,
                                          -4.59674047e-01, 2.60963397e+00, -9.70881273e+00,
                                          2.36585911e+01, -3.63978587e+01, 3.20703511e+01])
            else:
                raise NotImplementedError("invalid name for NIRCam wedge occulter")

            sigmas = np.poly1d(polyfitcoeffs)(scalefact)

            sigmar = sigmas * np.abs(y)
            # clip sigma: The minimum is to avoid divide by zero
            #             the maximum truncates after the first sidelobe to match the hardware
            sigmar.clip(min=np.finfo(sigmar.dtype).tiny, max=2 * np.pi, out=sigmar)
            self.transmission = (1 - (np.sin(sigmar) / sigmar) ** 2)
            self.transmission[y == 0] = 0  # special case center point (value based on L'Hopital's rule)
            # the bar should truncate at +- 10 arcsec:

            woutside = np.where(np.abs(x) > 10)
            self.transmission[woutside] = 1.0

        if self.nd_squares:
            # add in the ND squares. Note the positions are not exactly the same in the two wedges.
            # See the figures  in Krist et al. of how the 6 ND squares are spaced among the 5
            # coronagraph regions
            # Note: 180 deg rotation needed relative to Krist's figures for the flight SCI orientation:

            if ((self.module == 'A' and self.name == 'MASKLWB') or
                (self.module == 'B' and self.name == 'MASK210R')):
                # left edge:
                # has one fully in the corner and one half in the other corner, half outside the 10x10 box
                wnd_5 = np.where(
                    ((y < -5) & (y > -10)) &
                    (
                            ((x > 5) & (x < 10)) |
                            ((x < -7.5) & (x > -12.5))
                    )
                )
                wnd_2 = np.where(
                    ((y < 10) & (y > 8)) &
                    (
                            ((x > 8) & (x < 10)) |
                            ((x < -9) & (x > -11))
                    )
                )
            elif ((self.module == 'A' and self.name == 'MASK210R') or
                  (self.module == 'B' and self.name == 'MASKSWB')):
                # right edge
                wnd_5 = np.where(
                    ((y < -5) & (y > -10)) &
                    (
                            ((x < 12.5) & (x > 7.5)) |
                            ((x < -5) & (x > -10))
                    )
                )
                wnd_2 = np.where(
                    ((y < 10) & (y > 8)) &
                    (
                            ((x < 11) & (x > 9)) |
                            ((x < -8) & (x > -10))
                    )
                )
            else:
                # the others have two, one in each corner, both halfway out of the 10x10 box.
                wnd_5 = np.where(
                    ((y < -5) & (y > -10)) &
                    (np.abs(x) > 7.5) &
                    (np.abs(x) < 12.5)
                )
                wnd_2 = np.where(
                    ((y < 10) & (y > 8)) &
                    (np.abs(x) > 9) &
                    (np.abs(x) < 11)
                )

            self.transmission[wnd_5] = np.sqrt(1e-3)
            self.transmission[wnd_2] = np.sqrt(1e-3)

            # Add in the opaque border of the coronagraph mask holder.
            if ((self.module == 'A' and self.name == 'MASKLWB') or
                (self.module == 'B' and self.name == 'MASK210R')):
                # left edge
                woutside = np.where((x > 10) & (y > -11.5))
                self.transmission[woutside] = 0.0
            elif ((self.module == 'A' and self.name == 'MASK210R') or
                  (self.module == 'B' and self.name == 'MASKSWB')):
                # right edge
                woutside = np.where((x < -10) & (y > -11.5))
                self.transmission[woutside] = 0.0
            # mask holder edge
            woutside = np.where(y > 10)
            self.transmission[woutside] = 0.0

            # edge of mask itself
            # TODO the mask edge is complex and partially opaque based on CV3 images?
            # edge of glass plate rather than opaque mask I believe. To do later.
            # The following is just a temporary placeholder with no quantitative accuracy.
            # but this is outside the coronagraph FOV so that's fine - this only would matter in
            # modeling atypical/nonstandard calibration exposures.
            wedge = np.where((y < -11.5) & (y > -13))
            self.transmission[wedge] = 0.7

        if not np.isfinite(self.transmission.sum()):
            # stop()
            _log.warn("There are NaNs in the BLC mask - correcting to zero. (DEBUG LATER?)")
            self.transmission[np.where(np.isfinite(self.transmission) == False)] = 0
        return self.transmission

    def display(self, annotate=False, annotate_color='cyan', annotate_text_color=None, grid_size=20, *args, **kwargs):
        """Same as regular display for any other optical element, except adds annotate option
        for the LWB offsets """
        poppy.AnalyticOpticalElement.display(self, grid_size=grid_size, *args, **kwargs)
        if annotate:

            shift_dx = getattr(self, 'shift_x', 0) - getattr(self, 'bar_offset', 0)
            shift_dy = getattr(self, 'shift_y', 0)

            if annotate_text_color is None:
                annotate_text_color = annotate_color
            if self.name.lower() == 'maskswb' or self.name.lower() == 'masklwb':
                offset = self.offset_swb if self.name.lower() == 'maskswb' else self.offset_lwb
                for filt, offset in offset.items():
                    if 'W' in filt:
                        horiz, vert, voffset = 'right', 'top', -0.5
                    else:
                        horiz, vert, voffset = 'left', 'bottom', +0.5
                    matplotlib.pyplot.plot(offset + shift_dx, shift_dy, marker='+', color=annotate_color, clip_on=True)
                    matplotlib.pyplot.text(offset + shift_dx, voffset + shift_dy, filt, color=annotate_text_color, rotation=75,
                                           horizontalalignment=horiz, verticalalignment=vert, clip_on=True)
            ax = matplotlib.pyplot.gca()
            # Fix the axis scaling if any of the overplots exceeded it
            ax.set_xlim(-grid_size / 2, grid_size / 2)
            ax.set_ylim(-grid_size / 2, grid_size / 2)


# Helper functions for NIRcam occulters.
# The following are no longer used in practice, but were used to derive the
# table of polynomial coefficients that is now hard-coded inside
# the NIRCam_BandLimitedCoron case for the nircam wedge occulters.


def _width_blc(desired_width, approx=None, plot=False):
    """ The calculation of sigma parameters for the wedge BLC function is not straightforward.

    This function numerically solves the relevant equation to determine the sigma required to
    achieve a given HWHM.

    It uses recursion to iterate to a higher precision level.
    """

    loc = desired_width

    if approx is None:
        sigma = np.linspace(0, 20, 5000)
    else:
        sigma = np.linspace(approx * 0.9, approx * 1.1, 100000.)
    lhs = loc * np.sqrt(1 - np.sqrt(0.5))
    rhs = np.sin(sigma * loc) / sigma
    diff = np.abs(lhs - rhs)
    wmin = np.where(diff == np.nanmin(diff))
    sig_ans = sigma[wmin][0]

    if approx:
        return sig_ans
    else:
        # use recursion
        sig_ans = _width_blc(loc, sig_ans)

    if plot:
        check = (1 - (np.sin(sig_ans * loc) / sig_ans / loc) ** 2) ** 2
        # plt.plot(sigma, lhs)
        plt.clf()
        plt.plot(sigma, rhs)
        plt.axhline(lhs)

        print("sigma = %f implies HWHM = %f" % (sig_ans, loc))
        print(" check: 0.5 == %f" % (check))
    return sig_ans


def _calc_blc_wedge(deg=4, wavelength=2.1e-6):
    """ This function determines the desired sigma coefficients required to
    achieve a wedge from 2 to 6 lam/D.

    It returns the coefficients of a polynomial fit that maps from
    nlambda/D to sigma.

    """
    import scipy
    r = np.linspace(2, 6, 161)
    difflim = wavelength / 6.5 * 180. * 60 * 60 / np.pi
    sigs = [_width_blc(difflim * ri) for ri in r]

    pcs = scipy.polyfit(r, sigs, deg)
    p = np.poly1d(pcs)
    plt.plot(r, sigs, 'b')
    plt.plot(r, p(r), "r--")
    diffs = (sigs - p(r))
    print("Poly fit:" + repr(pcs))
    print("  fit rms: " + str(diffs.std()))


def _trim_nan_image(xgrid, ygrid, zgrid):
    """NaN Trimming of Image
    
    Remove rows/cols with NaN's while trying to preserve
    the maximum footprint of real data.
    """
    
    xgrid2, ygrid2, zgrid2 = xgrid, ygrid, zgrid
    
    # Create a mask of NaN'ed values
    nan_mask = np.isnan(zgrid2)
    nrows, ncols = nan_mask.shape
    # Determine number of NaN's along each row and col
    num_nans_cols = nan_mask.sum(axis=0)
    num_nans_rows = nan_mask.sum(axis=1)
    
    # First, crop all rows/cols that are only NaN's
    xind_good = np.where(num_nans_cols < nrows)[0]
    yind_good = np.where(num_nans_rows < ncols)[0]
    # get border limits
    x1, x2 = (xind_good.min(), xind_good.max()+1)
    y1, y2 = (yind_good.min(), yind_good.max()+1)
    # Trim of NaN borders
    xgrid2 = xgrid2[x1:x2]
    ygrid2 = ygrid2[y1:y2]
    zgrid2 = zgrid2[y1:y2,x1:x2]
    
    # Find a optimal rectangle subsection free of NaN's
    # Iterative cropping
    ndiff = 5
    while np.isnan(zgrid2.sum()):
        # Make sure ndiff is not negative
        if ndiff<0:
            break

        npix = zgrid2.size

        # Create a mask of NaN'ed values
        nan_mask = np.isnan(zgrid2)
        nrows, ncols = nan_mask.shape
        # Determine number of NaN's along each row and col
        num_nans_cols = nan_mask.sum(axis=0)
        num_nans_rows = nan_mask.sum(axis=1)

        # Look for any appreciable diff row-to-row/col-to-col
        col_diff = num_nans_cols - np.roll(num_nans_cols,-1) 
        row_diff = num_nans_rows - np.roll(num_nans_rows,-1)
        # For edge wrapping, just use last minus previous
        col_diff[-1] = col_diff[-2]
        row_diff[-1] = row_diff[-2]
        
        # Keep rows/cols composed mostly of real data 
        # and where number of NaN's don't change dramatically
        xind_good = np.where( ( np.abs(col_diff) <= ndiff  ) & 
                              ( num_nans_cols < 0.5*nrows ) )[0]
        yind_good = np.where( ( np.abs(row_diff) <= ndiff  ) & 
                              ( num_nans_rows < 0.5*ncols ) )[0]
        # get border limits
        x1, x2 = (xind_good.min(), xind_good.max()+1)
        y1, y2 = (yind_good.min(), yind_good.max()+1)
    
        # Trim of NaN borders
        xgrid2 = xgrid2[x1:x2]
        ygrid2 = ygrid2[y1:y2]
        zgrid2 = zgrid2[y1:y2,x1:x2]
        
        # Check for convergence
        # If we've converged, reduce 
        if npix==zgrid2.size:
            ndiff -= 1
                
    # Last ditch effort in case there are still NaNs
    # If so, remove rows/cols 1 by 1 until no NaNs
    while np.isnan(zgrid2.sum()):
        xgrid2 = xgrid2[1:-1]
        ygrid2 = ygrid2[1:-1]
        zgrid2 = zgrid2[1:-1,1:-1]
            
    return xgrid2, ygrid2, zgrid2


def _fix_zgrid_NaNs(xgrid, ygrid, zgrid, rot_ang=0):
    """Fix NaN's in Zernike Grid
    
    We trim NaN's within `zgrid`, then generate an extrapolation function
    using `RegularGridInterpolator`. A rotation angle can also be specified
    to maximize the number of remaining data points due to irregular
    polygons of the real `zgrid` data.
    
    Returns `zgrid` with the NaN's fixed using the extrapolation function.
    
    Parameter
    =========
    xgrid : ndarray
        1D V2 regular grid information
    ygrid : ndarray
        1D V3 regular grid information
    zgrid : ndarray
        2D Zernike grid
    rot_ang : float
        Option to rotate grid data for more optimal
        trimming of NaN's.
    """
            
    # Rotate zgrid
    if rot_ang != 0:
        zgrid = rotate(zgrid, rot_ang, reshape=False, order=1, cval=np.nan)
        
    # There will be some NaN's along the border that need to be replaced
    ind_nan = np.isnan(zgrid)
    # Remove rows/cols with NaN's
    xgrid2, ygrid2, zgrid2 = _trim_nan_image(xgrid, ygrid, zgrid)
    
    # Create regular grid interpolator function for extrapolation of NaN's
    func = RegularGridInterpolator((ygrid2,xgrid2), zgrid2, method='linear',
                                   bounds_error=False, fill_value=None)

    # Replace NaNs
    X, Y = np.meshgrid(xgrid,ygrid)
    pts = np.array([Y[ind_nan], X[ind_nan]]).transpose()
    zgrid[ind_nan] = func(pts)

    # De-rotate clipped zgrid image and redo RegularGridInterpolator
    if rot_ang != 0:
        # De-rotate
        zgrid = rotate(zgrid, -rot_ang, reshape=False, order=1, cval=np.nan)
        # There will be some NaNs along the border that need to be replaced
        ind_nan = np.isnan(zgrid)
        # Remove rows/cols 1 by 1 until no NaNs
        xgrid2, ygrid2, zgrid2 = _trim_nan_image(xgrid, ygrid, zgrid)

        # Create regular grid interpolator function for extrapolation of NaN's
        func = RegularGridInterpolator((ygrid2,xgrid2), zgrid2, method='linear',
                                       bounds_error=False, fill_value=None)

        # Replace NaNs
        pts = np.array([Y[ind_nan], X[ind_nan]]).transpose()
        zgrid[ind_nan] = func(pts)
        
    return zgrid


def _get_initial_pupil_sampling(instrument):
    """Utility function to retrieve the sampling of the first plane in some optical system.

    Returns: npix, pixelscale
    """
    # Determine the pupil sampling of the first aperture in the
    # instrument's optical system
    if isinstance(instrument.pupil, poppy.OpticalElement):
        # This branch needed to handle the OTE Linear Model case
        npix = instrument.pupil.shape[0]
        pixelscale = instrument.pupil.pixelscale
    else:
        # these branches to handle FITS files, by name or as an object
        if isinstance(instrument.pupil, fits.HDUList):
            pupilheader = instrument.pupil[0].header
        else:
            pupilfile = os.path.join(instrument._datapath, "OPD", instrument.pupil)
            pupilheader = fits.getheader(pupilfile)

        npix = pupilheader['NAXIS1']
        pixelscale = pupilheader['PUPLSCAL'] * units.meter / units.pixel
    return npix, pixelscale


# Field dependent aberration class for JWST instruments
class WebbFieldDependentAberration(poppy.OpticalElement):
    """ Field dependent aberration generated from Zernikes measured in ISIM CV testing

    Parameters
    -----------
    include_oversize : bool
        Explicitly model the 4% oversize for pupil tolerance

    """

    def __init__(self, instrument, include_oversize=False, **kwargs):
        super(WebbFieldDependentAberration, self).__init__(
            name="Aberrations",
            **kwargs
        )

        self.instrument = instrument
        self.instr_name = instrument.name

        # work out which name to index into the CV results with, if for NIRCam
        is_nrc_coron = False  # Define NRC coronagraph variable for conciseness 
        if instrument.name == 'NIRCam':
            channel = instrument.channel[0].upper()
            lookup_name = "NIRCam{channel}W{module}".format(
                channel=channel,
                module=instrument.module
            )
            # Check for coronagraphy; Set is_ncr_coron to True for Lyot pupil mask
            pupil_mask = self.instrument._pupil_mask
            is_nrc_coron = (pupil_mask is not None) and ( ('LYOT' in pupil_mask.upper()) or ('MASK' in pupil_mask.upper()) )
        elif instrument.name == 'FGS':
            # 'GUIDER1' or 'GUIDER2'
            assert instrument.detector in ('FGS1', 'FGS2')
            lookup_name = 'Guider' + instrument.detector[3]
        else:
            lookup_name = instrument.name
        _log.debug("Retrieving Zernike coefficients for " + lookup_name)

        self.tel_coords = instrument._tel_coords()

        # load the Zernikes table here         
        zfile = "si_zernikes_isim_cv3.fits"
        # Check special case NIRCam coronagraphy
        if is_nrc_coron:
            zfile = "si_zernikes_coron_wfe.fits"
        zernike_file = os.path.join(utils.get_webbpsf_data_path(), zfile)

        if not os.path.exists(zernike_file):
            raise RuntimeError("Could not find Zernike coefficients file {} \
                               in WebbPSF data directory".format(zfile))
        else:
            self.ztable_full = Table.read(zernike_file)

        npix, self.pixelscale = _get_initial_pupil_sampling(self.instrument)

        self.ztable = self.ztable_full[self.ztable_full['instrument'] == lookup_name]

        # Figure out the closest field point

        telcoords_am = self.tel_coords.to(units.arcmin).value
        v2 = self.ztable['V2']
        v3 = self.ztable['V3']
        r = np.sqrt((telcoords_am[0] - v2) ** 2 + (telcoords_am[1] - v3) ** 2)
        closest = np.argmin(r)

        # Save closest ISIM CV3 WFE measured field point for reference
        self.row = self.ztable[closest]

        self.name = "{instrument} internal WFE at V2V3=({v2:.2f},{v3:.2f})', near {field_point}".format(
            instrument=lookup_name,
            field_point=self.row['field_point_name'],
            v2=telcoords_am[0], v3=telcoords_am[1]
        )
        self.si_wfe_type = ("Interpolated", 
                "SI WFE was interpolated between available meas.")

        # Retrieve those Zernike coeffs
        # Field point interpolation
        v2_tel, v3_tel = telcoords_am
        coeffs = []
        for i in range(1, 37):
            zkey = 'Zernike_{}'.format(i)
            zvals = self.ztable[zkey]

            # Cubic interpolation of of non-uniform 2D grid
            cf = griddata((v2, v3), zvals, (v2_tel, v3_tel), method='cubic').tolist()
            
            # Want to perform extrapolation if field point outside of bounds
            if np.isnan(cf):
                if i==1:
                    self.si_wfe_type = ("Extrapolated", 
                            "SI WFE was extrapolated outside available meas.")

                # To extrapolate outside the measured field points, we proceed 
                # in two steps.  This first creates a fine-meshed cubic fit 
                # over the known field points, fixes any NaN's using 
                # RegularGridInterpolator, then again uses RegularGridInterpolator 
                # on the fixed data to extrapolate the requested field point.

                # In principle, the first call of RegularGridInterpolator can be 
                # used to extrapolate the requested field point to eliminate 
                # the intermediate step, but this method enables use of all the 
                # real data rather than the trimmed data set. RGI is a rather
                # quick process, so added overheads should be negligible.

                # Full field V2/V3 limits for each instrument.
                # Produces better initial extrapolation with fewer 
                # interpolation artifacts in RGI.
                if lookup_name == 'Guider1':
                    v2_min, v2_max, v3_min, v3_max = (2.2, 4.7, -12.9, -10.4)
                elif lookup_name == 'Guider2':
                    v2_min, v2_max, v3_min, v3_max = (-0.8, 1.6, -12.9, -10.4)
                elif lookup_name == 'NIRISS': 
                    v2_min, v2_max, v3_min, v3_max = (-6.0, -3.6, -12.9, -10.4)
                elif lookup_name == 'MIRI': 
                    v2_min, v2_max, v3_min, v3_max = (-8.3, -6.1, -7.3, -5.2)
                elif lookup_name == 'NIRSpec': 
                    v2_min, v2_max, v3_min, v3_max = (3.7, 9.0, -9.8, -4.5)
                elif (lookup_name == 'NIRCamLWA') or (lookup_name == 'NIRCamSWA'): 
                    v2_min, v2_max, v3_min, v3_max = (0.2, 2.7, -9.5, -7.0)
                elif (lookup_name == 'NIRCamLWB') or (lookup_name == 'NIRCamSWB'): 
                    v2_min, v2_max, v3_min, v3_max = (-2.7, -0.2, -9.5, -7.0)
                else:
                    v2_min, v2_max, v3_min, v3_max = (v2.min(), v2.max(), v3.min(), v3.max())

                # For NIRCam coronagraphy, add 50" to V3 limits
                if is_nrc_coron:
                    v3_min += 50. / 60.
                    v3_max += 50. / 60.

                # Create fine mesh grid
                dstep = 1. / 60. # 1" steps
                xgrid = np.arange(v2_min, v2_max+dstep, dstep)
                ygrid = np.arange(v3_min, v3_max+dstep, dstep)
                X, Y = np.meshgrid(xgrid,ygrid)

                # Cubic interpolation of all points
                # Will produce a number of NaN's that need to be extrapolated over
                zgrid = griddata((v2, v3), zvals, (X, Y), method='cubic')

                # Want to rotate zgrid image of some SIs to minimize NaN clipping
                if 'NIRSpec' in lookup_name:
                    rot_ang = 43
                elif 'MIRI' in lookup_name:
                    rot_ang = -5
                elif 'NIRISS' in lookup_name:
                    rot_ang = 2
                else:
                    rot_ang = 0

                # Fix the NaN's within zgrid array
                # Perform specified rotation for certain SIs
                # Trim rows/cols
                zgrid = _fix_zgrid_NaNs(xgrid, ygrid, zgrid, rot_ang=rot_ang)

                # Create final function for extrapolation
                func = RegularGridInterpolator((ygrid,xgrid), zgrid, method='linear',
                                               bounds_error=False, fill_value=None)
                # Extrapolate at requested (V2,V3) coordinates
                cf = func( (v3_tel, v2_tel) ).tolist()

            coeffs.append(cf)

        self.zernike_coeffs = coeffs

        # Generate an OPD on the same sampling as the input wavefront -
        # but implicitly inverted in coordinate system
        # to match the OTE exit pupil orientation

        if include_oversize:
            # Try to model the oversized gaps around the internal pupils.
            # This is only relevant if you are trying to model pupil shear or rotations,
            # and in general we don't have good WFE data outside the nominal pupil anyway
            # so let's leave this detail off by default.

            # internal pupils for NIRISS and MIRI instruments are 4 percent
            # oversized tricontagons
            if self.instrument.name == "NIRISS":
                self.amplitude = fits.getdata(os.path.join(
                    utils.get_webbpsf_data_path(),
                    'tricontagon_oversized_4pct.fits.gz')
                )
                # cut out central region to match the OPD, which is hard coded
                # to 1024
                self.amplitude = self.amplitude[256:256 + 1024, 256:256 + 1024]
            elif self.instrument.name == "MIRI":
                self.amplitude = fits.getdata(os.path.join(
                    utils.get_webbpsf_data_path(),
                    'MIRI',
                    'optics',
                    'MIRI_tricontagon_oversized_rotated.fits.gz')
                )

            else:
                # internal pupil is a 4 percent oversized circumscribing circle?
                # For NIRCam:
                # John stansberry 2016-09-07 reports "It is definitely oversized, but isn't really
                # circular... Kinda vaguely 6-sided I guess. [...] I can dig up
                # a drawing and/or some images that show the pupil stop."
                y, x = np.indices((npix, npix), dtype=float)
                y -= (npix - 1) / 2.0
                x -= (npix - 1) / 2.0
                r = np.sqrt(y ** 2 + x ** 2)
                self.amplitude = (r < (npix - 1) / 2.0 * 1.04).astype(int)

            self.opd = poppy.zernike.opd_from_zernikes(
                coeffs,
                npix=npix,
                aperture=self.amplitude,
                outside=0
            )
        else:
            self.opd = poppy.zernike.opd_from_zernikes(
                coeffs,
                npix=npix,
                outside=0
            )
            self.amplitude = (self.opd != 0).astype(int)
    
    def header_keywords(self):
        """ Return info we would like to save in FITS header of output PSFs
        """
        from collections import OrderedDict
        keywords = OrderedDict()
        keywords['SIWFETYP'] = self.si_wfe_type
        keywords['SIWFEFPT'] = (self.row['field_point_name'], "Closest ISIM CV3 WFE meas. field point")
        for i in range(1, 36):
            keywords['SIZERN{}'.format(i)] = (self.zernike_coeffs[i - 1], "[m] SI WFE coeff for Zernike term {}".format(i))
        return keywords

    # wrapper just to change default vmax
    def display(self, *args, **kwargs):
        if 'opd_vmax' not in kwargs:
            kwargs.update({'opd_vmax': 2.5e-7})

        return super(WebbFieldDependentAberration, self).display(*args, **kwargs)


class NIRSpecFieldDependentAberration(WebbFieldDependentAberration):
    """ Subclass that adds to the above the division into fore-optics
    and spectrograph optics for NIRSpec.

    The available end-to-end optical test data for NIRSpec from ISIM CV3
    do not allow distinguishing which optical planes have which amounts of
    aberration. However, the NIRSpec team performed extensive metrology
    during the assembly of NIRSpec FM2, both of individual components and
    of the assembled system using a shack hartmann WFS temporarily placed within
    the optical system.   [should add document number here to the report with those data!]

    Based on those data, Maurice Te Plate recommended to Marshall Perrin that the
    CV3 WFE should be apportioned 1/3 to the fore-optics and 2/3 to the
    spectrograph optics (collimator & camera). Given the uncertainties and
    available data that seems sufficiently precise for current purposes.

    """

    def __init__(self, instrument, where='fore', **kwargs):
        super(NIRSpecFieldDependentAberration, self).__init__(instrument, **kwargs)

        if where == 'fore':
            self.name = 'NIRSpec fore-optics WFE, near {}'.format(self.row['field_point_name'])
            self.scalefactor = 1. / 3
        else:
            self.name = 'NIRSpec spectrograph WFE, near {}'.format(self.row['field_point_name'])
            self.scalefactor = 2. / 3

        # apply scale factor to split up the OPD, and that's all we need to do.
        self.opd *= self.scalefactor


class NIRCamFieldAndWavelengthDependentAberration(WebbFieldDependentAberration):
    """ Subclass that adds to the above the wavelength dependent variation in defocus for
    NIRCam.

    The model for this is based on NIRCam models and ISIM CV2 test data, as
    provided by Randal Telfer to Marshall Perrin. It uses a combination of
    model design predictions continuously at all wavelengths based on the
    properties of the glasses in the refractive optical design, plus some small
    tweaks to achieve better agreement with the CV test measurements of defocus
    at a small subset of wavelengths.

    """

    def __init__(self, instrument, **kwargs):
        super(
            NIRCamFieldAndWavelengthDependentAberration,
            self).__init__(
            instrument,
            **kwargs)

        # Polynomial equations fit to defocus model. Wavelength-dependent focus 
        # results should correspond to Zernike coefficients in meters.
        # Fits were performed to the SW and LW optical design focus model 
        # as provided by Randal Telfer. 
        # See plot at https://github.com/spacetelescope/webbpsf/issues/179
        # The relative wavelength dependence of these focus models are very
        # similar for coronagraphic mode in the Zemax optical prescription,
        # so we opt to use the same focus model in both imaging and coronagraphy.
        defocus_to_rmswfe = -1.09746e7 # convert from mm defocus to meters (WFE)
        sw_focus_cf = np.array([-5.169185169, 50.62919436, -201.5444129, 415.9031962,  
                                -465.9818413, 265.843112, -59.64330811]) / defocus_to_rmswfe
        lw_focus_cf = np.array([0.175718713, -1.100964635, 0.986462016, 1.641692934]) / defocus_to_rmswfe
        self.fm_short = np.poly1d(sw_focus_cf)
        self.fm_long  = np.poly1d(lw_focus_cf)

        # Coronagraphic tilt (`ctilt`) offset model
        # Primarily effects the LW channel (approximately a 0.031mm diff from 3.5um to 5.0um).
        # SW module is small compared to LW, but we include it for completeness.
        # Values have been determined using the Zernike offsets as reported in the 
        # NIRCam Zemax models. The center reference positions will correspond to the 
        # NIRCam target acquisition filters (3.35um for LW and 2.1um for SW)
        sw_ctilt_cf = np.array([125.849834, -289.018704]) / 1e9
        lw_ctilt_cf = np.array([146.827501, -2000.965222, 8385.546158, -11101.658322]) / 1e9
        self.ctilt_short = np.poly1d(sw_ctilt_cf)
        self.ctilt_long  = np.poly1d(lw_ctilt_cf)

        # Get the representation of focus in the same Zernike basis as used for
        # making the OPD. While it looks like this does more work here than needed
        # by making a whole basis set, in fact because of caching behind the scenes
        # this is actually quick
        basis = poppy.zernike.zernike_basis_faster(
            nterms=len(self.zernike_coeffs),
            npix=self.opd.shape[0],
            outside=0
        )
        self.defocus_zern = basis[3]
        self.tilt_zern = basis[2]

        
    def get_opd(self, wave):
        """
        Parameters
        ----------
        wave : float or obj
            either a scalar wavelength (meters) or a Wavefront object
        """

        if isinstance(wave, poppy.Wavefront):
            wavelength = wave.wavelength
        else:
            wave = poppy.Wavefront(wavelength=float(wave))
            wavelength = wave.wavelength

        # Check for coronagraphy
        pupil_mask = self.instrument._pupil_mask
        is_nrc_coron = (pupil_mask is not None) and ( ('LYOT' in pupil_mask.upper()) or ('MASK' in pupil_mask.upper()) )

        # Which wavelength was used to generate the OPD map we have already
        # created from zernikes?
        if self.instrument.channel.upper() == 'SHORT':
            focusmodel = self.fm_short
            opd_ref_wave = 2.12
            opd_ref_focus = focusmodel(opd_ref_wave)
        else:
            focusmodel = self.fm_long
            opd_ref_wave = 3.23
            # All LW WFE measurements were made using F323N,
            # which has it's own focus that deviates from focusmodel().
            # But only do this for direct imaging SI WFE values,
            # because coronagraph WFE was measured in Zemax (no additional focus power).
            if is_nrc_coron:
                opd_ref_focus = focusmodel(opd_ref_wave)
            else:
                opd_ref_focus = 1.206e-7 # Not coronagraphy (e.g., imaging)

        # If F323N or F212N, then no focus offset necessary
        wave_um = wavelength.to(units.micron).value
        if ('F323N' in self.instrument.filter) or ('F212N' in self.instrument.filter):
            deltafocus = 0
        else:
            deltafocus = focusmodel(wave_um) - opd_ref_focus

        _log.info("  Applying OPD focus adjustment based on NIRCam focus vs wavelength model")
        _log.info("  Modified focus from {} to {} um: {:.3f} nm wfe".format(
            opd_ref_wave, wave_um, -deltafocus * 1e9)
        )

        # Apply defocus
        mod_opd = self.opd - deltafocus * self.defocus_zern

        # Apply wavelength-dependent tilt offset for coronagraphy
        # We want the reference wavelength to be that of the target acq filter
        # Final offset will position TA ref wave at the OPD ref wave location
        #   (wave_um - opd_ref_wave) - (ta_ref_wave - opd_ref_wave) = wave_um - ta_ref_wave
        if is_nrc_coron:
            if self.instrument.channel.upper() == 'SHORT':
                ctilt_model = self.ctilt_short
                ta_ref_wave = 2.10
            else: 
                ctilt_model = self.ctilt_long
                ta_ref_wave = 3.35

            tilt_offset = ctilt_model(wave_um) - ctilt_model(ta_ref_wave)
            _log.info("  Applying OPD tilt adjustment based on NIRCam tilt vs wavelength model")
            _log.info("  Modified tilt from {} to {} um: {:.3f} nm wfe".format(
                ta_ref_wave, wave_um, tilt_offset * 1e9)
            )

            # Apply tilt offset
            mod_opd = mod_opd + tilt_offset * self.tilt_zern

        rms = np.sqrt((mod_opd[mod_opd != 0] ** 2).mean())
        _log.info("  Resulting OPD has {:.3f} nm rms".format(rms * 1e9))

        return mod_opd


class MIRIFieldDependentAberrationAndObscuration(WebbFieldDependentAberration):
    """ Subclass that adds to the above the field dependent obscuration
    from the MIRI internal calibration source pickoff mirror.

    The model for this was derived by Randal Telfer based on the optical model file
    OTE_MIRI_20150223.seq, provided by Scott Rohrbach.

    In this case we do turn on by default the tricontagon outline since we have to worry about
    pupil shape anyway.

    """

    def __init__(self, instrument, include_oversize=True, **kwargs):
        super(MIRIFieldDependentAberrationAndObscuration, self).__init__(
            instrument,
            include_oversize=include_oversize,
            **kwargs
        )

        # figure out the XAN, YAN coordinates in degrees,
        # since that is what Randal's linear model expects

        xanyan = instrument._xan_yan_coords().to(units.degree)

        xan = xanyan[0].value
        yan = xanyan[1].value

        #    Telfer:
        #    Here is the matrix that reproduces the projection of the
        #    obscuration on the primary mirror, V2-V3 coordinates and
        #    radius in mm, as a function of XAN,YAN in degrees.
        #    So, V2 is:
        #
        #           V2 = -20882.636 * XAN -680.661 * YAN  - 1451.682.
        #
        #                     XAN            YAN         Const
        #           V2     -20882.636     -680.661    -1451.682
        #           V3        815.955    26395.552    -2414.406
        #           Rad       176.864     -392.545      626.920

        # we implement the above here, and convert the outputs to meters:
        self.obsc_v2 = (-20882.636 * xan - 680.661 * yan - 1451.682) * 0.001
        self.obsc_v3 = (815.955 * xan + 26395.552 * yan - 2414.406) * 0.001
        self.obsc_r = (176.864 * xan - 392.545 * yan + 626.920) * 0.001

        # generate coordinates. N.B. this assumed hard-coded pixel scale and
        # array size.
        pixel_scale = constants.JWST_CIRCUMSCRIBED_DIAMETER / 1024
        y, x = poppy.Wavefront.pupil_coordinates((1024, 1024), pixel_scale)

        # Now, the v2 and v3 coordinates calculated above are as projected back to
        # the OTE entrance pupil
        # But the OTE exit pupil as seen at the MIRI internal pupil is rotated by
        # 5 degrees with respect to that, and flipped in handedness as well
        # (but only in V3, given webbpsf axes conventions relative to the definition of the V frame)
        # Therefore we must transform the v2 and v3 to match the wavefront coords at the
        # intermediate plane.

        angle = np.deg2rad(instrument._rotation)
        proj_v2 = np.cos(angle) * self.obsc_v2 - np.sin(angle) * self.obsc_v3
        proj_v3 = -np.sin(angle) * self.obsc_v2 + np.cos(angle) * self.obsc_v3

        # handle V3 flip from OTE entrance to exit pupils
        # no flip needed for V2 since that's already implicitly done between
        # the V frame looking "in" to the OTE vs WebbPSF simulations looking
        # "out" from the detector toward the sky.
        proj_v3 *= -1

        mask = np.sqrt((y - proj_v3) ** 2 + (x - proj_v2) ** 2) < self.obsc_r
        self.amplitude[mask] = 0

        # No need to subclass any of the methods; it's sufficient to set the custom
        # amplitude mask attribute value.


# Alternative implementation that just reads OPDs from some file
class LookupTableFieldDependentAberration(poppy.OpticalElement):
    """ Retrieve OPDs from a lookup table over many field points.
    This is pretty much a hack, hard-coded for a specific data delivery from Ball!
    Intended for OTE team WFR4 and MIMF KDP Practice data prep, not generalized beyond that.

    Parameters
    -----------
    add_niriss_defocus: bool
        add 0.8 microns PTV defocus to NIRISS only (for WFR4 test)
    rm_ptt: bool
        Remove piston, tip, tilt
    rm_center_ptt : bool
        If rm_ptt, use the center value for each detector rather than per field point
    nwaves: float
        Number of waves to defocus SM, if add_sm_defocus_pos or add_sm_defocus_neg is True.
    add_sm_defocus: bool
        If True, add "nwaves" of SM defocus, measured at a reference wavelength of 2.0 microns.


    Usage:
    ------
        inst = webbpsf.NIRCam()  # or any other SI
        inst._si_wfe_class = LookupTableFieldDependentAberration()

    """

    def __init__(self, instrument, field_points_file=None, phasemap_file=None,
                 which_exercise='MIMF_KDP_2',
                 add_niriss_defocus=None, rm_ptt=None, rm_center_ptt=None,
                 add_mimf_defocus=False, add_sm_defocus=False, nwaves=None, **kwargs):
        super().__init__(
            name="Aberrations",
            **kwargs
        )
        import warnings

        self.instrument = instrument
        self.instr_name = instrument.name

        self.instrument.pupilopd=None  # Do not add in the usual telescope WFE on top of this;
                                       # This model provided by Ball includes both the telescope and the SI WFE combined.

        self.rm_ptt = rm_ptt

        self.which_exercise = which_exercise

        if self.which_exercise == 'WFR4':
            add_niriss_defocus=True
            rm_ptt = True
            rm_center_ptt = True
        elif self.which_exercise == 'MIMF_KDP':
            add_niriss_defocus=False
            rm_ptt = False
            rm_center_ptt = False
        elif self.which_exercise == 'LRE4' or self.which_exercise == 'LRE4-OTE26':
            add_niriss_defocus=False
            rm_ptt = False
            rm_center_ptt = False
        elif self.which_exercise == 'MIMF_KDP_2':
            add_niriss_defocus=False
            rm_ptt = False
            rm_center_ptt = False
 

        if self.instr_name =='NIRCam':
            self.instr_name += " "+self.instrument.module
        elif self.instr_name == 'FGS':
            self.instr_name = self.instrument.detector

        self.tel_coords = instrument._tel_coords()

        # load the OPD lookup map table (datacube) here

        import webbpsf.constants
        if self.which_exercise == 'WFR4':
            fp_path = '/ifs/jwst/tel/wfr4_mirage_sims/phase_maps_from_ball/'
            if field_points_file is None:
                field_points_file = fp_path + 'The_Field_Coordinates.txt'
            if phasemap_file is None:
                phasemap_file = fp_path + 'phase_maps.fits'

            self.phasemap_file = phasemap_file

            self.table = Table.read(field_points_file, format='ascii', names=('V2', 'V3'))

            self.yoffset = -7.8

            self.table['V3'] += self.yoffset # Correct from -YAN to actual V3

            self.phasemaps = fits.getdata(phasemap_file)

            self.phasemap_pixelscale = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER/256 * units.meter / units.pixel
            resample = True

        elif self.which_exercise == 'MIMF_KDP':
            fp_path = '/ifs/jwst/tel/MIMF_KDP_Practice/Ball_Phase_Maps/'
            field_points_file = fp_path + 'coordinates.txt'
            self.table = Table.read(field_points_file, format='ascii.basic', names=('XWAS', 'YWAS'))

            # Convert coordinate table to V2V3 in arcminutes
            self.table['V2'] = -self.table['XWAS']
            self.table['V3'] = self.table['YWAS'] - 468/60

            phasemap_file = fp_path + 'all_26Feb2021.fits'
            self.phasemaps = fits.getdata(phasemap_file)
            self.phasemaps = self.phasemaps.reshape(7*11*11, 256, 256)
            self.phasemap_pixelscale = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER / 256 * units.meter / units.pixel
            resample = True
        elif self.which_exercise == 'LRE4' or self.which_exercise == 'LRE4-OTE26':
            fp_path = '/ifs/jwst/tel/LRE4/from_ball/'
            if self.which_exercise == 'LRE4':
                field_points_file = fp_path + 'coordinates.ecsv'
                phasemap_file = fp_path + 'rescaled_opds_for_OTE-25.2.fits'
            elif self.which_exercise == 'LRE4-OTE26':
                field_points_file = fp_path + 'coordinates-ote26.ecsv'
                phasemap_file = fp_path + 'rescaled_opds_for_OTE-26.fits'

            self.table = Table.read(field_points_file)
            self.phasemaps = fits.getdata(phasemap_file)
            # Phase maps have been pre-zoomed in this case by the import notebook
            resample = False
            self.phasemap_pixelscale = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER / 1024 * units.meter / units.pixel
        elif self.which_exercise == 'MIMF_KDP_2':
            fp_path = '/ifs/jwst/tel/MIMF_KDP_Practice_Sept2021/Ball_Phase_Maps/'

            # Convert coordinate table to V2V3 in arcminutes
            xcoords = fits.getdata(fp_path+"xcor.fits")
            ycoords = fits.getdata(fp_path+"ycor.fits")
            V2 = -xcoords.flatten()
            V3 =  ycoords.flatten() - 468/60
            self.table = Table([V2,V3], names=['V2','V3'])

            phasemap_file = fp_path + 'complete_wf.fits'
            self.phasemaps = fits.getdata(phasemap_file)
            self.phasemaps = self.phasemaps.reshape(7*11*11, 256, 256)
            self.phasemap_pixelscale = webbpsf.constants.JWST_CIRCUMSCRIBED_DIAMETER / 256 * units.meter / units.pixel
            resample = True
 
        self.phasemap_file = phasemap_file


        # Determine the pupil sampling of the first aperture in the
        # instrument's optical system
        if isinstance(instrument.pupil, poppy.OpticalElement):
            # This branch needed to handle the OTE Linear Model case
            self.pixelscale = instrument.pupil.pixelscale
        else:
            # these branches to handle FITS files, by name or as an object
            if isinstance(instrument.pupil, fits.HDUList):
                pupilheader = instrument.pupil[0].header
            else:
                pupilfile = os.path.join(instrument._datapath, "OPD", instrument.pupil)
                pupilheader = fits.getheader(pupilfile)

            self.pixelscale = pupilheader['PUPLSCAL'] * units.meter / units.pixel

        # Figure out the closest field point

        telcoords_am = self.tel_coords.to(units.arcmin).value
        print(f"Requested field point has coord {telcoords_am}")

        v2 = self.table['V2']
        v3 = self.table['V3']
        r = np.sqrt((telcoords_am[0] - v2) ** 2 + (telcoords_am[1] - v3) ** 2)
        closest = np.argmin(r)   # if there are two field points with identical coords or equal distance just one is returned

        print(f"Closest field point is row {closest}:\n{self.table[closest]}")

        # Save closest ISIM CV3 WFE measured field point for reference
        self.row = self.table[closest]


        self.name = "{instrument} at V2V3=({v2:.2f},{v3:.2f}) Lookup table WFE from ({v2t:.2f},{v3t:.2f})".format(
            instrument=self.instr_name,
            v2=telcoords_am[0], v3=telcoords_am[1],
            v2t=self.row['V2'], v3t=self.row['V3']

        )
        self.si_wfe_type = ("Lookup Table",
                "SI + OTE WFE from supplied lookup table of phase maps.")

        # Retrieve the phase map

        phasemap = self.phasemaps[closest]

        # The phase maps are provided in OTE entrance pupil orientation, however we need it to be
        # in exit pupil orientation, so flip it vertically here.
        phasemap = phasemap[::-1]
        print("Flipped input phase map vertically into exit pupil orientation.")

        if resample:
            if phasemap.shape[0] != 256:
                raise NotImplementedError("Hard coded for Ball delivery of 256 pixel phase maps")
    
            # Resample to 1024 across, by replicating each pixel into a 4x4 block
            resample_factor = 4
            phasemap_big = np.kron(phasemap, np.ones((resample_factor,resample_factor)))
        else:
            # no resampling / zooming needed
            phasemap_big = phasemap

        self.opd = phasemap_big * 1e-6   # Convert from microns to meters
        self.amplitude = np.ones_like(self.opd)

        if rm_ptt:
            apmask = self.opd != 0
            if rm_center_ptt:
                # Remove the PTT values at the center of each instrument, rather than per field point. This
                # leaves in the field dependence but takes out the bulk offset
                # These values are just a precomputed lookup table of the coefficients returned by the
                # opd_expand_nonorthonormal call just below, for the center field point on each.

               coeffs_per_si = {"NIRCam A": [-3.50046880e-10, -7.29120639e-08, -1.39751567e-08],
                                "NIRCam B": [-2.45093780e-09, -2.51804001e-07, -2.64821753e-07],
                                "NIRISS": [-1.49297771e-09, -2.11111038e-06, -3.99881993e-07],
                                "FGS1": [ 9.86180620e-09, -5.94041500e-07,  1.18953161e-06],
                                "FGS2": [ 4.84327424e-09, -8.24285481e-07,  5.09791593e-07],
                                "MIRI": [-8.75766849e-09, -1.27850277e-06, -1.03467567e-06],}
               coeffs = coeffs_per_si[self.instr_name]
            else:
                coeffs = poppy.zernike.opd_expand_nonorthonormal(self.opd, aperture=apmask, nterms=3)
            ptt_only = poppy.zernike.opd_from_zernikes(coeffs, aperture=apmask, npix=self.opd.shape[0], outside=0)
            self.opd -= ptt_only
            print(f"Removing piston, tip, tilt from the input wavefront. Coeffs for  {self.instr_name}: {coeffs},")

        if add_mimf_defocus:
            self.instrument.options['defocus_waves'] = 0.8
            self.instrument.options['defocus_wavelength'] = 1e-6  # Add 0.8 microns PTV defocus

        if add_niriss_defocus and self.instr_name=='NIRISS':
            # The Ball delivery was supposed to have defocused NIRISS for rehearsal purposes, but didn't.
            # So fix that here.
            self.instrument.options['defocus_waves'] = 0.8
            self.instrument.options['defocus_wavelength'] = 1e-6  # Add 0.8 microns PTV defocus
            warnings.warn("Adding defocus=0.8 waves for NIRISS!")

        if add_sm_defocus:
            if nwaves:
                print("ADDING DEFOCUS {:4.1f} WAVES at 2.0 microns".format(nwaves))
                self.instrument.options['defocus_waves'] = nwaves
                self.instrument.options['defocus_wavelength'] = 2.0e-6
            else:
                print("Not adding any defocus; set nwaves")


    def header_keywords(self):
        """ Return info we would like to save in FITS header of output PSFs
        """
        from collections import OrderedDict
        keywords = OrderedDict()
        keywords['SIWFETYP'] = self.si_wfe_type
        keywords['SIWFEFPT'] = ( f"{self.row['V2']:.3f}, {self.row['V3']:.3f}", "Closest lookup table meas. field point")
        keywords['SIWFEFIL'] = self.phasemap_file
        return keywords

    # wrapper just to change default vmax
    def display(self, *args, **kwargs):
        if 'opd_vmax' not in kwargs:
            kwargs.update({'opd_vmax': 2.5e-7})

        return super().display(*args, **kwargs)


class NIRCamFieldDependentWeakLens(poppy.OpticalElement):
    """Higher-fidelity model of NIRCam weak lens(es), based on calibrated as-built performance
    and field dependence.

    Includes field-dependent variations in defocus power, and in astigmatism. Includes variation of the
    +4 lens' effective OPD when used in a pair with either the +8 or -8 lens.

    These are modeled as the specific values from the nearest neighbor ISIM CV calibration point,
    with no interpolation between them included at this time.

    See R. Telfer, 'NIRCam Weak Lens Characterization and Performance', JWST-REF-046515

    Parameters
    -----------
    name : str
        WLP8, WLM8, WLP4, WLM4, WLP12.

    center_fp_only : bool
        For debugging; override to set no field dependence and just use the average center field point power

    include_power, include_astigmatism : bool
        Can be used to selectively enable/disable parts of the optical model. Intended for debugging; should no
        need to be set by users in general.

    """

    def __init__(self, name='WLP8', instrument=None, center_fp_only=False, verbose=False, include_power=True,
                 include_astigmatism=True, **kwargs):
        super().__init__(name=name)

        self.ref_wavelength = 2.12e-6  # reference wavelength for defocus

        self.verbose = verbose
        if instrument is None:
            self.module = 'A'
            self.v2v3_coords = (0, -468 / 60)
            npix = 1024
        else:
            self.module = instrument.module
            self.v2v3_coords = instrument._tel_coords()
            npix, pixelscale = _get_initial_pupil_sampling(instrument)

        self.ztable_full = None

        ## REFERENCE:
        # NIRCam weak lenses, values from WSS config file, PRDOPSFLT-027
        #                  A         B
        # WLP4_diversity =   8.27309     8.3443         diversity in microns
        # WLP8_diversity =  16.4554     16.5932
        # WLM8_diversity = -16.4143    -16.5593
        # WL_wavelength =    2.12                       Wavelength, in microns

        if center_fp_only or instrument is None:
            # use the center field point power only. No field dependence

            # Power in P-V waves at center field point in optical model
            # JWST-REF-046515, table 2      Mod A:    Mod B:
            power_at_center_fp = {'WLM8': (-8.0188, -7.9521),
                                  'WLM4': (-4.0285, -3.9766),
                                  'WLP4': (3.9797, 3.9665),
                                  'WLP8': (8.0292, 7.9675),
                                  'WLP12': (12.0010, 11.9275)}

            power_pv = power_at_center_fp[self.name][0 if self.module == 'A' else 1]
            astig0 = 0
            astig45 = 0

        else:
            closest_fp = self.find_closest_isim_fp_name(instrument)
            if verbose: print(closest_fp)
            power_pv, astig0, astig45 = self.lookup_empirical_lens_power(name, closest_fp)

        self.power_pv_waves = power_pv
        pv2rms_norm = self.ref_wavelength / (2 * np.sqrt(3))  # convert desired PV waves to RMS microns for power
        # since the below function wants inputs in RMS

        self.power_rms_microns = power_pv * pv2rms_norm

        zernike_coefficients = np.zeros(6)
        if include_power:
            zernike_coefficients[3] = self.power_rms_microns
        if include_astigmatism:
            zernike_coefficients[4] = astig0
            zernike_coefficients[5] = astig45
        self.zernike_coefficients = zernike_coefficients

        self.opd = poppy.zernike.opd_from_zernikes(
            zernike_coefficients,
            npix=npix,
            outside=0
        )
        self.amplitude = np.ones_like(self.opd)

    def find_closest_isim_fp_name(self, instr):
        """Find the closest ISIM CV field point to a given instrument object,
        i.e. the field point closest to the configured detector and coordinates
        """

        if self.ztable_full is None:
            zernike_file = os.path.join(utils.get_webbpsf_data_path(), "si_zernikes_isim_cv3.fits")
            self.ztable_full = Table.read(zernike_file)

        lookup_name = f"NIRCam{instr.channel.upper()[0]}W{instr.module}"
        ztable = self.ztable_full[self.ztable_full['instrument'] == lookup_name]

        self._ztable = ztable
        self._instr = instr
        telcoords_am = instr._tel_coords().to(units.arcmin).value
        if self.verbose: print(telcoords_am)
        r = np.sqrt((telcoords_am[0] - ztable['V2']) ** 2 + (telcoords_am[1] - ztable['V3']) ** 2)
        # Save closest ISIM CV3 WFE measured field point for reference
        row = ztable[r == r.min()]
        return row['field_point_name']

    def lookup_empirical_lens_power(self, lens_name, field_point_name):
        """ Lookup lens power and astigmatism versus field position, from empirical calibrations from ISIM CV testing

        """
        mypath = os.path.dirname(os.path.abspath(__file__)) + os.sep
        wl_data_file = os.path.join(mypath, 'otelm', 'NIRCam_WL_Empirical_Power.csv')
        wl_data = Table.read(wl_data_file, comment='#', header_start=1)

        field_point_row = wl_data[wl_data['Field'] == field_point_name]
        if self.verbose: print(field_point_row)

        defocus_name = lens_name[2:]

        power = field_point_row[defocus_name].data[0]
        # Fringe zernike coefficients, from Telfer's table
        z5 = field_point_row[defocus_name+"_Z5"].data[0]
        z6 = field_point_row[defocus_name + "_Z6"].data[0]

       # Have to convert Zernike normalization and order from fringe to noll, and nanometers to meters
        astig0 = z6 / np.sqrt(6)*1e-9
        astig45 = z5 / np.sqrt(6)*1e-9

        if self.verbose: print(power)
        return power, astig0, astig45
