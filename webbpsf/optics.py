from __future__ import division, print_function, absolute_import, unicode_literals
import os
import poppy
import numpy as np

import astropy.table
import astropy.io.fits as fits

from . import utils
from .webbpsf_core import _log


#######  Custom Optics used in JWInstrument classes  #####


class NIRSpec_three_MSA_shutters(poppy.AnalyticOpticalElement):
    """ Three NIRSpec MSA shutters, adjacent vertically."""

    def get_transmission(self,wave):
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

        y, x= wave.coordinates()

        self.transmission = np.zeros(wave.shape)
        # get the innermost shutter than spans the Y axis
        w_inside_1 = np.where( (abs(y) < (msa_height/2))  & (abs(x) < (msa_width/2)))
        self.transmission[w_inside_1] = 1
        # get the adjacent shutters one above and one below.
        w_inside_2 = np.where( (abs(y) > (msa_height/2)+msa_wall) & (abs(y) < msa_height*1.5+msa_wall)  & (abs(x) < (msa_width/2)))
        self.transmission[w_inside_2] = 1

        return self.transmission


class NIRSpec_MSA_open_grid(poppy.AnalyticOpticalElement):
    """ An infinite repeating region of the NIRSpec MSA grid"""

    def get_transmission(self,wave):
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

        y, x= wave.coordinates()
        #xnew =  x*np.cos(np.deg2rad(self.angle)) + y*np.sin(np.deg2rad(self.angle))
        #ynew = -x*np.sin(np.deg2rad(self.angle)) + y*np.cos(np.deg2rad(self.angle))
        #x,y = xnew, ynew

        mask_vert_walls  = np.abs(np.mod(np.abs(x), msa_x_pitch) - (msa_x_pitch/2)) < msa_wall/2
        mask_horz_walls  = np.abs(np.mod(np.abs(y), msa_y_pitch) - (msa_y_pitch/2)) < msa_wall/2



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

        I should have pointed that the power assumed in my simulations for the cylindrical lens was off. It was one of the conclusions of CV1RR. The actual radius of curvature of the cylinder is 25.3 meters (rather than the smaller figure I used before).

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
            #cylinder_radius=22.85,  cylinder_sag_mm=4.0, rotation_angle=92.25, rotate_mask=False, transmission=None,
            shift=None):
        # Initialize the base optical element with the pupil transmission and zero OPD


        if which=='LLNL':
            raise NotImplementedError("Rotated field mask for LLNL grism not yet implemented!")
        elif which=='Bach':
            transmission=os.path.join( utils.get_webbpsf_data_path(), "NIRISS/optics/MASKGR700XD.fits.gz")
        else:
            raise NotImplementedError("Unknown grating name:"+which)

        self.shift=shift
        poppy.AnalyticOpticalElement.__init__(self, name=name, planetype=poppy.poppy_core._PUPIL)

        # UPDATED NUMBERS 2013-07:
        # See Document FGS_TFI_UdM_035_RevD

        _log.debug("Computing properties for {0} grism".format(which))
        if which =='Bach':
            #---- Phase properties ---------------
            # 3.994 microns P-V over 27.02 mm measured (Loic's email)
            # This is **surface sag**, corresponding to P-V of 6.311 waves at lambda=632.8 nm.
            # should correspond to 3.698 microns over 26 mm clear aperture.
            self.prism_size = 0.02702 # 27.02 millimeters for the physical prism
            self.prism_clear_aperture = 0.0260 # 26 mm clear aperture for the prism + mount
            self.cylinder_rotation_angle = 2 # was 2.25

            #self.cylinder_radius = 22.85 # radius of curvature  ; Nominal
            # but they discarded that and used 25.3 instead
            # From Lafreniere's wfe_cylindricallens.pro:
            #  "OVERRIDE PREVIOUS CASES AFTER CV1RR RESULTS:"
            self.cylinder_radius = 25.3 # radius of curvature

            #---- Amplitude Transmission / Pupil shape ---------------
            self.pupil_size_mm = 26.0
            # Note that the IDL code says 26 mm is 683.75 pixels using the assumed demagnification
            self.pupil_rotation_angle = 2.0

        else:
            # 5.8 microns P-V over 32.15 mm (Loic's email)
            # should correspond to 4.38 microns over 28 mm clear aperture
            self.cylinder_radius = 22.39 # radius of curvature
            self.prism_size = 0.03215 # millimeters for the physical prism
            self.prism_clear_aperture = 0.0280 # clear aperture for the prism + mount
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
        #self.pupil_demagnification =  170.8367 # meters on the primary / meters in the NIRISS pupil
        #self.pupil_demagnification =  173.56 # meters on the primary / meters in the NIRISS pupil

        # Anand says:
        #  nominally the circumscribing circle at the PW of NIRISS is ~40mm.  I use 39mm for the nrm, but it's slightly field-dependent.  Compare that to the 6.6... PM circle?
        self.pupil_demagnification = 6.6 / 0.040 # about 165

        # perform an initial population of the OPD array for display etc.
        tmp = self.get_phasor(poppy.Wavefront(2e-6))

    def get_opd(self, wave):
        """ Make an OPD array corresponding to the cylindrical weak lens
        used for defocusing the spectrum in the perpendicular-to-dispersion direction.
        """

        if isinstance(wave, poppy.Wavefront):
            wavelength=wave.wavelength
        else:
            wavelength=float(wave)
            wave =  poppy.Wavefront(wavelength=wave)

        # compute indices in pixels, relative to center of plane, with rotation
        # units of these are meters
        y, x = wave.coordinates()

        ang = np.deg2rad(self.cylinder_rotation_angle )
        x = np.cos(ang)*x - np.sin(ang)*y
        y = np.sin(ang)*x + np.cos(ang)*y

        _log.debug(" Rotating local grism axes by {0} degrees".format(self.cylinder_rotation_angle))

        # From IDL code by David Lafreniere:
        #  ;the cylindrical defocus
        #x=(dindgen(pupdim)-pupdim/2)#replicate(1,pupdim)
        #y0=(rpuppix^2+sag[s]^2)/(2*sag[s])
        #wfe1=y0-sqrt(y0^2-x^2)
        #if sag[s] lt 1.e-5 then wfe1=0.d0

        # Here I will just translate that to Python exactly, making use of the
        # variables here:

        # rpuppix = radius of pupil in pixels
        #rpuppix = self.amplitude_header['DIAM'] / self.amplitude_header['PUPLSCAL'] / 2
        # Calculate the radius of curvature of the cylinder, bsaed on
        # the chord length and height

        # In this case we're assuming the cylinder is precisely as wide as the projected
        # telescope pupil. This doesn't seem guaranteed:
        #  * actual chord length across cylinder: 27.02 mm.
        #  * projected primary scale at NIRISS = ?

        _log.debug(" Computing GR700XD cylinder based on RoC: {0:.3g} meters".format(self.cylinder_radius))
        _log.debug(" Computing GR700XD cylinder based on pupil demagnification: {0:.3g} primary to grism".format(self.pupil_demagnification))


        # Compute the overall sag of the cylinder lens at its outer edge. This is not actually used, it's
        # just for cross-check of the values
        # the sag will depend on half the pupil size since that's the offset from center to edge
        sag0 = np.sqrt(self.cylinder_radius**2 - (self.prism_size/2)**2) - self.cylinder_radius
        _log.debug(" Computed GR700XD cylinder sag at lens outer edge (for cross check only): {0:.3g} meters".format(sag0))

        # now compute the spatially dependent sag of the cylinder, as projected onto the primary

        # what is the pupil scale at the *reimaged pupil* of the grism?
        pupil_scale_m_per_pix = 38.0255e-6 # Based on UdeM info in wfe_cylindricallens.pro
        #sag = np.sqrt(self.cylinder_radius**2 - (x*self.amplitude_header['PUPLSCAL']/self.pupil_demagnification)**2) - self.cylinder_radius
        sag = np.sqrt(self.cylinder_radius**2 - (x/self.pupil_demagnification)**2) - self.cylinder_radius
        #sag = self.cylinder_radius -  np.sqrt(self.cylinder_radius**2 - (x * pupil_scale_m_per_pix )**2 )


        # what we really want to do is take the physical properties of the as-built optic, and interpolate into that
        # to compute the OPD after remapping based on the pupil scale (and distortion?)
        #y0=(rpuppix**2+self.cylinder_sag**2)/(2*self.cylinder_sag)
        #wfe1=y0-np.sqrt(y0**2-x**2)

        _log.debug(" Cylinder P-V: {0:.4g} meters physical sag across full array".format(sag.max()-sag.min()) )

        #fits.writeto('py_opd.fits', sag*(self.ZnSe_index(wavelength) -1), clobber=True)
        # remove piston offset
        #sag -= sag[wnz].min()   # normalize to 0 at the minimum
        #sag -= sag[wnz].mean()    # normalize around the mean
        sag[self.amplitude == 0] = 0 # no OPD in opaque regions (makes no difference in propagation but improves display)
        wnz = np.where(self.amplitude != 0) # this is just for display of the log messages:
        _log.debug(" Cylinder P-V: {0:.4g} meters physical sag across clear aperture".format(sag[wnz].max()-sag[wnz].min()) )


        # scale for index of refraction
        index = self.ZnS_index(wavelength)
        opd = sag *  ( index-1)
        _log.debug(" Scaling for ZnS index of refraction {0} at {1:.3g} microns".format(index, wavelength*1e6))
        _log.debug(" Cylinder P-V: {0:.4g} meters optical sag at {1:.3g} microns across clear aperture".format(opd[wnz].max()-opd[wnz].min(), wavelength*1e6) )
        return opd

    def get_transmission(self, wave):
        """ Make array for the pupil obscuration appropriate to the grism
        """

        if isinstance(wave, poppy.Wavefront):
            wavelength=wave.wavelength
        else:
            wavelength=float(wave)
            wave =  poppy.Wavefront(wavelength=wave)
        y, x = wave.coordinates()
        ang = np.deg2rad(self.pupil_rotation_angle )
        x = np.cos(ang)*x - np.sin(ang)*y
        y = np.sin(ang)*x + np.cos(ang)*y

        _log.debug("Rotating local pupil mask axes by {0} degrees".format(self.cylinder_rotation_angle))

        pupil_halfsize_m = self.pupil_size_mm / 2  / 1000 * self.pupil_demagnification
        pupilmask = np.ones_like(x)
        pupilmask[np.abs(x) > pupil_halfsize_m] = 0
        pupilmask[np.abs(y) > pupil_halfsize_m] = 0

        return pupilmask


    def ZnS_index(self, wavelength, temperature=40):
        """ Return cryogenic index of refraction of ZnS (Cleartran)

        Based on IDL function index_cleartran provided by Loic Albert at U de Montreal
        Which was in turn based on Leviton and Fray 2013
        http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1744938
        doi:10.1117/12.2024817

        """
        lambda_micron = wavelength*1e6

        # Sellmeier dispersion model
        #From Leviton & Frey measurements (SPIE preprint) (assumes lambda in microns)
        S_1 = np.asarray([[3.35933,-5.12262e-4,1.01086e-5,-4.14798e-8,6.91051e-11]])
        S_2 = np.asarray([[0.706131,4.89603e-4,-8.91159e-6,3.81621e-8,-6.54805e-11]])
        S_3 = np.asarray([[4.02154,-2.93193e-2,2.31080e-4,-7.57289e-07,8.31188e-10]])
        S_ij = np.concatenate( (S_1,S_2,S_3), axis=0)
        lambda_1 = np.array([[0.161151, -8.93057E-06, 2.73286E-07, -1.23408E-09, 2.29917E-12]])
        lambda_2 = np.array([[0.282427, -4.66636E-05, 7.55906E-07, -2.77513E-09, 4.35237E-12]])
        lambda_3 = np.array([[41.1590, -0.161010, 1.23906E-03, -3.95895E-06, 4.16370E-09]])
        lambda_ij = np.concatenate((lambda_1,lambda_2,lambda_3))


        n2minus1 = 0.0
        T=temperature
        for i in range(3):
            S_i =       S_ij[i,0]      + S_ij[i,1]     *T + S_ij[i,2]     *T**2.0 + S_ij[i,3]     *T**3.0 + S_ij[i,4]     *T**4.0
            lambda_i =  lambda_ij[i,0] + lambda_ij[i,1]*T + lambda_ij[i,2]*T**2.0 + lambda_ij[i,3]*T**3.0 + lambda_ij[i,4]*T**4.0
            n2minus1 += S_i*lambda_micron**2.0/(lambda_micron**2.0 - lambda_i**2.0)

        cleartran_index = np.sqrt(1.0 + n2minus1)
        return cleartran_index


    def display(self, opd_vmax=6e-6, *args, **kwargs):
        "Same as regular display for any other optical element, except opd_vmax default changed"
        poppy.AnalyticOpticalElement.display(self,*args, opd_vmax=opd_vmax, **kwargs)


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
    def __init__(self):
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

        pupil_mag = 6.603464/39.0
        poppy.CompoundAnalyticOptic.__init__( self, (
                poppy.SecondaryObscuration( secondary_radius = 6.0*pupil_mag,
                                                      support_width = 2.0*pupil_mag,
                                                      n_supports = 3,
                                                      support_angle_offset=90+180), # align first support with +V2 axis
                                                                                # but invert to match OTE exit pupil
                poppy.CircularAperture( radius = 39 * pupil_mag /2) ), name = 'CLEARP')



class NIRCam_BandLimitedCoron(poppy.BandLimitedCoron):
    """ Band Limited Coronagraph """
    allowable_kinds = ['nircamcircular', 'nircamwedge']
    """ Allowable types of BLC supported by this class"""

    def __init__(self, name="unnamed BLC", kind='nircamcircular',  module='A',
            **kwargs):
        super(NIRCam_BandLimitedCoron, self).__init__(name=name, kind=kind, **kwargs)
        if module not in ['A','B']:
            raise ValueError("module parameter must be 'A' or 'B'.")
        self.module=module

        if self.name=='MASK210R':
            self.sigma = 5.253
            self.kind = 'nircamcircular'
        elif self.name=='MASK335R':
            self.sigma=3.2927866
            self.kind = 'nircamcircular'
        elif self.name=='MASK430R':
            self.sigma=2.58832
            self.kind = 'nircamcircular'
        elif self.name == 'MASKSWB':
            self.kind = 'nircamwedge'
            # coeffs set in lookup table inside getPhasor
        elif self.name == 'MASKLWB':
            self.kind = 'nircamwedge'
            # coeffs set in lookup table inside getPhasor
        else:
            raise NotImplementedError("invalid name for NIRCam occulter: "+self.name)



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
        if self.kind == 'nircamcircular':
            r = np.sqrt(x ** 2 + y ** 2)
            sigmar = self.sigma * r
            # clip sigma: The minimum is to avoid divide by zero
            #             the maximum truncates after the first sidelobe to match the hardware
            sigmar.clip(np.finfo(sigmar.dtype).tiny, 2*np.pi, out=sigmar)  # avoid divide by zero -> NaNs
            self.transmission = (1 - (2 * scipy.special.jn(1, sigmar) / sigmar) ** 2)
            self.transmission[r==0] = 0   # special case center point (value based on L'Hopital's rule)

        elif self.kind == 'nircamwedge':
            # This is hard-coded to the wedge-plus-flat-regions shape for NIRCAM

            # the scale fact should depend on X coord in arcsec, scaling across a 20 arcsec FOV.
            # map flat regions to 2.5 arcsec each
            # map -7.5 to 2, +7.5 to 6. slope is 4/15, offset is +9.5
            scalefact = (2 + (-x + 7.5) * 4 / 15).clip(2, 6)

            # Working out the sigma parameter vs. wavelength to get that wedge pattern is non trivial
            # This is NOT a linear relationship. See calc_blc_wedge helper fn below.

            if self.name == 'MASKSWB': #np.abs(self.wavelength - 2.1e-6) < 0.1e-6:
                polyfitcoeffs = np.array([2.01210737e-04, -7.18758337e-03, 1.12381516e-01,
                                          -1.00877701e+00, 5.72538509e+00, -2.12943497e+01,
                                          5.18745152e+01, -7.97815606e+01, 7.02728734e+01])
                scalefact = scalefact[:, ::-1] # flip orientation left/right for SWB mask
            elif self.name == 'MASKLWB': #elif np.abs(self.wavelength - 4.6e-6) < 0.1e-6:
                polyfitcoeffs = np.array([9.16195583e-05, -3.27354831e-03, 5.11960734e-02,
                                          -4.59674047e-01, 2.60963397e+00, -9.70881273e+00,
                                          2.36585911e+01, -3.63978587e+01, 3.20703511e+01])
            else:
                raise NotImplementedError("invalid name for NIRCam wedge occulter")

            sigmas = scipy.poly1d(polyfitcoeffs)(scalefact)

            sigmar = sigmas * np.abs(y)
            # clip sigma: The minimum is to avoid divide by zero
            #             the maximum truncates after the first sidelobe to match the hardware
            sigmar.clip(min=np.finfo(sigmar.dtype).tiny, max=2*np.pi, out=sigmar)
            self.transmission = (1 - (np.sin(sigmar) / sigmar) ** 2)
            # TODO pattern should be truncated past first sidelobe
            self.transmission[x==0] = 0   # special case center point (value based on L'Hopital's rule)
            # the bar should truncate at +- 10 arcsec:
            woutside = np.where(np.abs(x) > 10)
            self.transmission[woutside] = 1.0


        # add in the ND squares. Note the positions are not exactly the same in the two wedges.
        # See the figures  in Krist et al. of how the 6 ND squares are spaced among the 5
        # corongraph regions
        # Note: 180 deg rotation needed relative to Krist's figures for the flight SCI orientation:
        x = x[::-1, ::-1]
        y = y[::-1, ::-1]
        if ((self.module=='A' and self.name=='MASKLWB') or
            (self.module=='B' and self.name=='MASK210R')):
            # left edge:
            # has one fully in the corner and one half in the other corner, half outside the 10x10 box
            wnd_5 = np.where(
                ((y > 5)&(y<10)) &
                (
                    ((x < -5) & (x > -10)) |
                    ((x > 7.5) & (x < 12.5))
                )
            )
            wnd_2 = np.where(
                ((y > -10)&(y<-8)) &
                (
                    ((x < -8) & (x > -10)) |
                    ((x > 9) & (x < 11))
                )
            )
        elif ((self.module=='A' and self.name=='MASK210R') or
              (self.module=='B' and self.name=='MASKSWB')):
            # right edge
            wnd_5 = np.where(
                ((y > 5)&(y<10)) &
                (
                    ((x > -12.5) & (x < -7.5)) |
                    ((x > 5) & (x <10))
                )
            )
            wnd_2 = np.where(
                ((y > -10)&(y<-8)) &
                (
                    ((x > -11) & (x < -9)) |
                    ((x > 8) & (x<10))
                )
            )
        else:
            # the others have two, one in each corner, both halfway out of the 10x10 box.
            wnd_5 = np.where(
                ((y > 5)&(y<10)) &
                (np.abs(x) > 7.5) &
                (np.abs(x) < 12.5)
            )
            wnd_2 = np.where(
                ((y > -10)&(y<-8)) &
                (np.abs(x) > 9) &
                (np.abs(x) < 11)
            )

        self.transmission[wnd_5] = np.sqrt(1e-3)
        self.transmission[wnd_2] = np.sqrt(1e-3)



        # Add in the opaque border of the coronagraph mask holder.
        if ((self.module=='A' and self.name=='MASKLWB') or
            (self.module=='B' and self.name=='MASK210R')):
            # left edge
            woutside = np.where((x < -10) & (y < 11.5 ))
            self.transmission[woutside] = 0.0
        elif ((self.module=='A' and self.name=='MASK210R') or
              (self.module=='B' and self.name=='MASKSWB')):
            # right edge
            woutside = np.where((x > 10) & (y < 11.5))
            self.transmission[woutside] = 0.0
        # mask holder edge
        woutside = np.where(y < -10)
        self.transmission[woutside] = 0.0

        # edge of mask itself
        # TODO the mask edge is complex and partially opaque based on CV3 images?
        # edge of glass plate rather than opaque mask I believe. To do later.
        # The following is just a temporary placeholder with no quantitative accuracy.
        # but this is outside the coronagraph FOV so that's fine - this only would matter in
        # modeling atypical/nonstandard calibration exposures.

        wedge = np.where(( y > 11.5) & (y < 13))
        self.transmission[wedge] = 0.7







        if not np.isfinite(self.transmission.sum()):
            #stop()
            _log.warn("There are NaNs in the BLC mask - correcting to zero. (DEBUG LATER?)")
            self.transmission[np.where(np.isfinite(self.transmission) == False)] = 0
        return self.transmission



# Helper functions for NIRcam occulters.
# The following are no longer used in practice, but were used to derive the
# table of polynomial coefficients that is now hard-coded inside
# the NIRCam_BandLimitedCoron case for the nircam wedge occulters.


def _width_blc(desired_width, approx=None, plot=False):
    """ The calculation of sigma parameters for the wedge BLC function is not straightforward.

    This function numerically solves the relevant equation to determine the sigma required to
    acheive a given HWHM.

    It uses recursion to iterate to a higher precision level.
    """

    loc = desired_width

    if approx is None:
        sigma = np.linspace(0, 20, 5000)
    else:
        sigma = np.linspace(approx*0.9, approx*1.1, 100000.)
    lhs = loc* np.sqrt(1 - np.sqrt(0.5))
    rhs = np.sin(sigma * loc) / sigma
    diff = np.abs(lhs - rhs)
    wmin = np.where(diff == np.nanmin(diff))
    sig_ans = sigma[wmin][0]

    if approx:
        return sig_ans
    else:
        # use recursion
        sig_ans = width_blc(loc, sig_ans)

    if plot:
        check =  (1-  (np.sin(sig_ans * loc)/sig_ans/loc)**2)**2
        #plt.plot(sigma, lhs)
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
    difflim = wavelen / 6.5 * 180.*60*60/np.pi
    sigs = [_width_blc(difflim * ri) for ri in r]

    pcs = scipy.polyfit(r, sigs, deg)
    p = scipy.poly1d(pcs)
    plt.plot(r, sigs, 'b')
    plt.plot(r, p(r), "r--")
    diffs = (sigs - p(r))
    print("Poly fit:" +repr(pcs))
    print("  fit rms: "+str(diffs.std()))



# Field dependent aberration class for JWST instruments

class JWST_Field_Dependent_Aberration(poppy.OpticalElement):
    """ Field dependent aberration generated from Zernikes measured in ISIM CV testing


    """
    def __init__(self, instrument, **kwargs):
        ""
        super(JWST_Field_Dependent_Aberration, self).__init__(name="Aberrations", **kwargs)

        self.instrument = instrument
        self.instr_name = instrument.name

        #work out which name to index into the CV results with, if for NIRCam
        if instrument.name == 'NIRCam':
            lookup_name = ("NIRCAM"+instrument.channel[0]+"W"+instrument.module).upper()
        elif instrument.name == 'FGS':
            lookup_name = 'GUIDER'+instrument.detector[3] # 'GUIDER1' or 'GUIDER2'
        else:
            lookup_name = instrument.name.upper()
        _log.debug("Retrieving zernike coeffs for "+lookup_name)

        self.tel_coords = instrument._tel_coords()

        # load the Zernikes table here

        self.ztable_full = astropy.table.Table.read(
            os.path.join(utils.get_webbpsf_data_path(),'zernikes_isim_cv2.fits'))
        # Determine the pupil sampling of the first aperture in the instrument's optical system
        npix=1024 # hard code for now


        self.ztable=self.ztable_full[self.ztable_full['instrument']==lookup_name]

        #Figure out the closest field point
        v2 = self.ztable['V2']
        v3 = self.ztable['V3']
        r = np.sqrt((self.tel_coords[0]-v2)**2+(self.tel_coords[1]-v3)**2)
        closest = np.argmin(r)

        self.row = self.ztable[closest]

        #self.name = '{} near {} ({:.3f}, {:.3f})'.format(lookup_name,self.row['field_point_name'],*self.tel_coords)
        self.name = '{} near {}'.format(lookup_name,self.row['field_point_name'])
        # Retrieve those Zernike coeffs (no interpolation for now)
        coeffs = [self.row['Zernike_{}'.format(i)] for i in range(1,36)]

        # Generate an OPD on the same sampling as the input wavefront -
        # but implicitly inverted in coordinate system
        # to match the OTE exit pupil orientation

        self.opd = poppy.zernike.opd_from_zernikes(coeffs, npix=npix,
                outside=0)  #*1e6 # convert to microns

        self.amplitude=fits.getdata(
            os.path.join(utils.get_webbpsf_data_path(),'tricontagon.fits'))
        #?   No the SI internal clear aperture is larger in general...
        # oversized tricontagon could be better
        #self.amplitude = np.ones_like(self.opd, dtype=float)

    def display(self, opd_vmax=2.5e-7, *args, **kwargs):
        # wrapper just to change default vmax
        super(JWST_Field_Dependent_Aberration, self).display(opd_vmax=opd_vmax, *args, **kwargs)



