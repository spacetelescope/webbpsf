"""
==============================
WFIRST Instruments (pre-alpha)
==============================

WARNING: This model has not yet been validated against other PSF
         simulations, and uses several approximations (e.g. for
         mirror polishing errors, which are taken from HST).
"""

import os.path
import poppy
import numpy as np
from . import webbpsf_core
from scipy.interpolate import griddata
from astropy.io import fits
import logging
_log = logging.getLogger('webbpsf').setLevel('DEBUG')

class WavelengthDependenceInterpolator(object):
    """WavelengthDependenceInterpolator can be configured with
    `n_zernikes` worth of Zernike coefficients at up to `n_wavelengths`
    wavelengths, and will let you `get_aberration_terms` for any
    wavelength in range interpolated linearly between measured/known
    points
    """
    def __init__(self, n_wavelengths=16, n_zernikes=22):
        self._n_wavelengths = n_wavelengths
        self._n_zernikes = n_zernikes
        self._aberration_terms = np.zeros((n_wavelengths, n_zernikes), dtype=np.float64)
        self._wavelengths = []

    def set_aberration_terms(self, wavelength, zernike_array):
        """Supply a reference `wavelength` and a `zernike_array`
        (of length `n_zernikes`) where the aberration is known
        """
        n_wavelengths_set = len(self._wavelengths)
        if wavelength not in self._wavelengths and n_wavelengths_set < self._n_wavelengths:
            self._wavelengths.append(wavelength)
            aberration_row_idx = n_wavelengths_set  # which is now index of last row
        elif wavelength in self._wavelengths:
            aberration_row_idx = self._wavelengths.index(wavelength)
        else:
            # can't add more wavelengths without allocating new _aberration_terms array
            raise ValueError("Already have information at {} wavelengths "
                             "(pass larger n_wavelengths to __init__?)".format(self._n_wavelengths))
        if len(zernike_array) != self._n_zernikes:
            raise ValueError("Expected {} aberration terms (pass different "
                             "n_zernikes to __init__?)".format(self._n_zernikes))
        self._aberration_terms[aberration_row_idx] = zernike_array

    def get_aberration_terms(self, wavelength):
        """Return the Zernike coefficients as interpolated for this
        `wavelength`"""
        # return array of length n_zernikes interpolated for this wavelength
        if wavelength in self._wavelengths:
            # aberration known exactly for this wavelength
            aberration_row_idx = self._wavelengths.index(wavelength)
            return self._aberration_terms[aberration_row_idx]
        else:
            # we have to interpolate @ this wavelength
            return griddata(self._wavelengths, self._aberration_terms, wavelength, method='linear')

    def get_zernike_coefficient(self, wavelength, zernike_subscript):
        """Return a single Zernike coefficient based on `wavelength`
        and `zernike_subscript` (counting from 1 in the Noll indexing
        convention)"""
        zernike_subscript -= 1  # 1-indexed to 0-indexed
        if zernike_subscript >= self._n_wavelengths:
            raise ValueError("No information at Zernike term Z_{}".format(zernike_subscript + 1))
        if wavelength in self._wavelengths:
            # aberration known exactly for this wavelength
            aberration_row_idx = self._wavelengths.index(wavelength)
            return self._aberration_terms[aberration_row_idx][zernike_subscript]
        else:
            return griddata(
                self._wavelengths,
                self._aberration_terms[:,zernike_subscript],
                wavelength,
                method='linear'
            )

class FieldDependentAberration(poppy.ZernikeWFE):
    """FieldDependentAberration incorporates aberrations that
    are interpolated in wavelength, x, and y pixel positions by
    computing the Zernike coefficients for a particular wavelength
    and position.
    """

    """By default, `get_aberration_terms` will zero out Z1, Z2, and Z3
    (piston, tip, and tilt) as they are not meaningful for telescope
    PSF calculations (the former is irrelevant, the latter two would
    be handled by a distortion solution). Change
    `_omit_piston_tip_tilt` to False to include the Z1-3 terms."""
    _omit_piston_tip_tilt = True

    def __init__(self, pixel_width, pixel_height,
                 name="Field-dependent Aberration", radius=1.0, oversample=1, interp_order=3):
        self.pixel_width, self.pixel_height = pixel_width, pixel_height
        self.x_pixel, self.y_pixel = pixel_width // 2, pixel_height // 2
        self._wavelength_interpolators = {}
        super(FieldDependentAberration, self).__init__(
            name=name,
            verbose=True,
            radius=radius,
            oversample=oversample,
            interp_order=interp_order
        )

    def getPhasor(self, wave):
        """Set the Zernike coefficients (for ZernikeWFE.getPhasor) based
        on the wavelength of the incoming wavefront and the pixel
        position
        """
        if not isinstance(wave, poppy.Wavefront):
            wavelength = wave
        else:
            wavelength = wave.wavelength
        self.coefficients = wavelength * self.get_aberration_terms(wavelength)
        return super(FieldDependentAberration, self).getPhasor(wave)

    def set_field_position(self, x_pixel, y_pixel):
        """Set the x and y pixel position on the detector for which to
        interpolate aberrations"""
        if x_pixel > self.pixel_width or x_pixel < 0:
            raise ValueError("Requested pixel_x position lies outside "
                             "the detector width ({})".format(x_pixel))
        if y_pixel > self.pixel_height or y_pixel < 0:
            raise ValueError("Requested pixel_y position lies outside "
                             "the detector height ({})".format(y_pixel))

        self.x_pixel, self.y_pixel = x_pixel, y_pixel

    def add_field_point(self, x_pixel, y_pixel, interpolator):
        """Supply a wavelength-space interpolator for a pixel position
        on the detector"""
        self._wavelength_interpolators[(x_pixel, y_pixel)] = interpolator

    def get_aberration_terms(self, wavelength):
        """Supply the Zernike coefficients for the aberration based on
        the wavelength and pixel position on the detector"""
        if (self.x_pixel, self.y_pixel) in self._wavelength_interpolators:
            # short path: this is a known point
            interpolator = self._wavelength_interpolators[(self.x_pixel, self.y_pixel)]
            coefficients = interpolator.get_aberration_terms(wavelength)
        else:
            # get aberrations at all field points
            field_points, aberration_terms = [], []
            for field_point_coords, point_interpolator in self._wavelength_interpolators.iteritems():
                field_points.append(field_point_coords)
                aberration_terms.append(point_interpolator.get_aberration_terms(wavelength))
            aberration_array = np.asarray(aberration_terms)
            assert len(aberration_array.shape) == 2, "computed aberration array is not 2D " \
                                                     "(inconsistent number of Zernike terms " \
                                                     "at each point?)"
            coefficients = griddata(
                np.asarray(field_points),
                np.asarray(aberration_terms),
                (self.x_pixel, self.y_pixel),
                method='linear'
            )
        if self._omit_piston_tip_tilt:
            _log.info("Omitting piston/tip/tilt")
            coefficients[:3] = 0.0  # omit piston, tip, and tilt Zernikes
        return coefficients

def _load_wfi_detector_aberrations(filename):
    from astropy.io import ascii
    zernike_table = ascii.read(filename)
    detectors = {}

    def build_detector_from_table(number, zernike_table):
        """Build a FieldDependentAberration optic for a detector using
        Zernikes Z1-Z22 at various wavelengths and field points"""
        single_detector_info = zernike_table[zernike_table['Cnfg#'] == number]
        field_points = set(single_detector_info['Field'])
        interpolators = {}
        detector = FieldDependentAberration(
            4096,
            4096,
            radius=WFIRSTInstrument.PUPIL_RADIUS,
            name="Field Dependent Aberration (SCA{:02})".format(number)
        )
        for field_id in field_points:
            field_point_rows = single_detector_info[single_detector_info['Field'] == field_id]
            local_x, local_y = field_point_rows[0]['Local_x'], field_point_rows[0]['Local_y']
            interpolator = build_wavelength_dependence(field_point_rows)

            midpoint_pixel = 4096 / 2
            # (local_x in mm / 10 um pixel size) -> * 1e2
            # local_x and _y range from -20.44 to +20.44, so adding to the midpoint pixel
            # makes sense to place (-20.44, -20.44) at (4, 4)
            pixx, pixy = int(midpoint_pixel + local_x * 1e2), int(midpoint_pixel + local_y * 1e2)

            detector.add_field_point(pixx, pixy, interpolator)
        return detector

    def build_wavelength_dependence(rows):
        """Build an interpolator object that interpolates Z1-Z22 in
        wavelength space"""
        wavelengths = set(rows['Wave(um)'])
        interpolator = WavelengthDependenceInterpolator(n_wavelengths=len(wavelengths),
                                                        n_zernikes=22)
        for row in rows:
            z = np.zeros(22)
            for idx in xrange(22):
                z[idx] = row['Z{}'.format(idx + 1)]
            interpolator.set_aberration_terms(row['Wave(um)'] * 1e-6, z)

        return interpolator

    detector_ids = set(zernike_table['Cnfg#'])
    for detid in detector_ids:
        detectors["SCA{:02}".format(detid)] = build_detector_from_table(detid, zernike_table)

    return detectors

class WFIRSTInstrument(webbpsf_core.SpaceTelescopeInstrument):
    PUPIL_RADIUS = 2.4 / 2.0
    """
    WFIRSTInstrument contains data and functionality common to WFIRST
    instruments, such as setting the pupil shape
    """
    telescope = "WFIRST"
    _detectors = {}
    """
    Dictionary mapping detector names to FieldDependentAberration optics.

    (Subclasses must populate this at `__init__`.)
    """
    _selected_detector = None
    """
    The name of the currently selected detector. Must be a key in _detectors, as validated by the
    `detector` property setter.

    (Subclasses must populate this at `__init__`.)
    """
    _detector_position = (2048.0, 2048.0)
    """
    The pixel position, accessed through the `detector_position` property
    and validated to be in-range on the selected detector.
    """

    def __init__(self, *args, **kwargs):
        super(WFIRSTInstrument, self).__init__(*args, **kwargs)
        self.pupil = os.path.join(self._WebbPSF_basepath, 'AFTA_WFC_C5_Pupil_Shortwave_Norm_2048px.fits')
        self.pupilopd = os.path.join(self._WebbPSF_basepath, 'upscaled_HST_OPD.fits')

        # n.b. WFIRSTInstrument subclasses must set these
        self._detectors = {}
        self._selected_detector = None

    @property
    def detector(self):
        """Detector selected for simulated PSF

        Used in calculation of field-dependent aberrations. Must be
        selected from detectors in the `detector_list` attribute.
        """
        return self._selected_detector

    @detector.setter
    def detector(self, value):
        if value.upper() not in self.detector_list:
            raise ValueError("Invalid detector. Valid detector names are: {}".format(', '.join(self.detector_list)))
        self._selected_detector = value.upper()

    @property
    def detector_list(self):
        """Detectors on which the simulated PSF could lie"""
        return sorted(self._detectors.keys())

    @property
    def detector_position(self):
        """The pixel position in (X, Y) on the detector"""
        return self._detector_position

    @detector_position.setter
    def detector_position(self, position):
        x_pixel, y_pixel = position
        detector = self._detectors[self._selected_detector]
        detector.set_field_position(x_pixel, y_pixel)

    def _get_aberrations(self):
        """Get the OpticalElement that applies the field-dependent
        optical aberrations. (Called in _getOpticalSystem.)"""
        return self._detectors[self._selected_detector]

    def _getFITSHeader(self, result, options):
        """Populate FITS Header keywords"""
        super(WFIRSTInstrument, self)._getFITSHeader(result, options)
        result[0].header['DETXPIXL'] = (self.detector_position[0], 'X pixel position (for field dependent aberrations)')
        result[0].header['DETYPIXL'] = (self.detector_position[1], 'Y pixel position (for field dependent aberrations)')
        result[0].header['DETECTOR'] = (self.detector, 'Detector selected')

class WFI(WFIRSTInstrument):
    """
    WFI represents to the to-be-named wide field imager
    for the WFIRST mission

    WARNING: This model has not yet been validated against other PSF
             simulations, and uses several approximations (e.g. for
             mirror polishing errors, which are taken from HST).
    """
    # "The H158, F184 and W149 filters and the grism are mounted with proximate cold pupil masks"
    # from the final draft of the SDT report, page 92, table 3-2
    UNMASKED_PUPIL_WAVELENGTH_MIN, UNMASKED_PUPIL_WAVELENGTH_MAX = 0.760e-6, 1.454e-6
    MASKED_PUPIL_WAVELENGTH_MIN, MASKED_PUPIL_WAVELENGTH_MAX = 1.380e-6, 2.000e-6

    def __init__(self):
        pixelscale = 110e-3  # arcsec/px, WFIRST-AFTA SDT report final version (p. 91)
        super(WFI, self).__init__("WFI", pixelscale=pixelscale)

        self._detectors = _load_wfi_detector_aberrations(os.path.join(self._datapath, 'zernikes.csv'))
        assert len(self._detectors.keys()) > 0
        self.detector = 'SCA01'
        self.detector_position = (2048.0, 2048.0)

        # Paths to the two possible pupils. The correct one is selected based on requested
        # wavelengths in _validate_config()
        self._unmasked_pupil_path = os.path.join(self._WebbPSF_basepath, 'AFTA_WFC_C5_Pupil_Shortwave_Norm_2048px.fits')
        self._masked_pupil_path = os.path.join(self._WebbPSF_basepath, 'AFTA_WFC_C5_Pupil_Mask_Norm_2048px.fits')

    def _validateConfig(self, **kwargs):
        """Validates that the WFI is configured sensibly

        This mainly consists of selecting the masked or unmasked pupil
        appropriately based on the wavelengths requested.
        """
        if self.pupil in (self._unmasked_pupil_path, self._masked_pupil_path):
            # Does the wavelength range fit entirely in an unmasked filter?
            wavelengths = np.array(kwargs['wavelengths'])
            wl_min, wl_max = np.min(wavelengths), np.max(wavelengths)
            # test shorter filters first; if wl range fits entirely in one of them, it's not going
            # to be the (masked) wideband filter
            if wl_min >= self.UNMASKED_PUPIL_WAVELENGTH_MIN and wl_max <= self.UNMASKED_PUPIL_WAVELENGTH_MAX:
                # use unmasked pupil optic
                self.pupil = self._unmasked_pupil_path
                _log.info("Using the unmasked WFI pupil shape based on wavelengths requested")
            elif wl_min >= self.MASKED_PUPIL_WAVELENGTH_MIN and wl_max <= self.MASKED_PUPIL_WAVELENGTH_MAX:
                # use masked pupil optic
                self.pupil = self._masked_pupil_path
                _log.info("Using the masked WFI pupil shape based on wavelengths requested")
            else:
                raise RuntimeError("Cannot figure out which WFI pupil shape to use from the "
                                   "wavelengths requested! (If you know which one you want, "
                                   "instantiate a poppy.FITSOpticalElement and assign it to the "
                                   "'pupil' attribute.)")
        else:
            # If the user has set the pupil to a custom value, let them worry about the
            # correct shape it should have
            pass
        super(WFI, self)._validateConfig(**kwargs)

class CGI(WFIRSTInstrument):
    """
    WFIRST Coronagraph Instrument

    """
    def __init__(self):
        pixelscale_BB = 0.01
        pixelscale_IFS = 0.01 
        pixelscale = 0.01 # arcsec/pix for broadband imager
        super(CGI, self).__init__("CGI", pixelscale=pixelscale)

        self.pupil = os.path.join(self._WebbPSF_basepath, 'AFTA_CGI_C5_Pupil_onax_1000px.fits')
        self.mode_list = ['CHARSPC', 'DISKSPC']
        self.image_mask_list = ['CHARSPC_F660', 'CHARSPC_F770', 'CHARSPC_F890', 'DISKSPC_F465', 'DISKSPC_F565', 'DISKSPC_F835', 'DISKSPC_F885']
        self.pupil_mask_list = ['SPC26D88']
        self.pupilopd = None
        self.aberration_optic = None
        self.options = {'force_coron':True}

    def _validate_config(self):
        if self.mode in self.mode_list:
            return True

    def _addAdditionalOptics(self, optsys, oversample=4):
        """Add coronagraphic or spectrographic optics for WFIRST CGI."""
        
        trySAM = False
        char_fpmres = 8
        disk_fpmres = 4

        if ('pupil_shift_x' in self.options and self.options['pupil_shift_x'] != 0) or \
           ('pupil_shift_y' in self.options and self.options['pupil_shift_y'] != 0):
            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        else: shift = None

        if self.mode == 'CHARSPC':
            optsys.addPupil(transmission=self._datapath+"optics/CHARSPC_SP_1000pix.fits.gz", name=self.mode, shift=None)
        elif self.mode == 'DISKSPC':
            optsys.addPupil(transmission=self._datapath+"optics/DISKSPC_SP_1000pix.fits.gz", name=self.mode, shift=None)

        if self.filter == 'F660':
            self.image_mask = 'CHARSPC_F660'
            optsys.addImage(transmission=self._datapath+"optics/CHARSPC_FPM_25WA90_2x65deg_-_FP1res%d_evensamp_D%d_F660.fits.gz"%(char_fpmres, 2*9*char_fpmres))
        elif self.filter == 'F770':
            self.image_mask = 'CHARSPC_F770'
            optsys.addImage(transmission=self._datapath+"optics/CHARSPC_FPM_25WA90_2x65deg_-_FP1res%d_evensamp_D%d_F770.fits.gz"%(char_fpmres, 2*9*char_fpmres))
        elif self.filter == 'F890':
            self.image_mask = 'CHARSPC_F890'
            optsys.addImage(transmission=self._datapath+"optics/CHARSPC_FPM_25WA90_2x65deg_-_FP1res%d_evensamp_D%d_F890.fits.gz"%(char_fpmres, 2*9*char_fpmres))

        elif self.filter == 'F465':
            self.image_mask = 'DISKSPC_F465'
            optsys.addImage(transmission=self._datapath+"optics/DISKSPC_FPM_65WA200_360deg_-_FP1res%d_evensamp_D%d_F465.fits.gz"%(disk_fpmres, 2*20*disk_fpmres))
        elif self.filter == 'F565':
            self.image_mask = 'DISKSPC_F565'
            optsys.addImage(transmission=self._datapath+"optics/DISKSPC_FPM_65WA200_360deg_-_FP1res%d_evensamp_D%d_F565.fits.gz"%(disk_fpmres, 2*20*disk_fpmres))
        elif self.filter == 'F835':
            self.image_mask = 'DISKSPC_F835'
            optsys.addImage(transmission=self._datapath+"optics/DISKSPC_FPM_65WA200_360deg_-_FP1res%d_evensamp_D%d_F835.fits.gz"%(disk_fpmres, 2*20*disk_fpmres))
        elif self.filter == 'F885':
            self.image_mask = 'DISKSPC_F885'
            optsys.addImage(transmission=self._datapath+"optics/DISKSPC_FPM_65WA200_360deg_-_FP1res%d_evensamp_D%d_F885.fits.gz"%(disk_fpmres, 2*20*disk_fpmres))
            
        if self.mode == 'CHARSPC' or self.mode == 'DISKSPC':
            self.pupil_mask = 'SPC26D88'
            optsys.addPupil(transmission=self._datapath+"optics/SPC_LS_26DS88_1000pix.fits.gz", name=self.pupil_mask, shift=shift)

        #if self.image_mask == 'CHARSPC_F660':
        #    optsys.addImage(transmission=self._datapath+"optics/CHARSPC_FPM_25WA90_2x65deg_-_FP1res%d_evensamp_D%d_F660.fits.gz"%(fpmres, 2*9*fpmres))
        #elif self.image_mask == 'CHARSPC_F770':
        #    optsys.addImage(transmission=self._datapath+"optics/CHARSPC_FPM_25WA90_2x65deg_-_FP1res%d_evensamp_D%d_F770.fits.gz"%(fpmres, 2*9*fpmres))
        #elif self.image_mask == 'CHARSPC_F890':
        #    optsys.addImage(transmission=self._datapath+"optics/CHARSPC_FPM_25WA90_2x65deg_-_FP1res%d_evensamp_D%d_F890.fits.gz"%(fpmres, 2*9*fpmres))

        #elif self.image_mask == 'DISKSPC_F465':
        #    optsys.addImage(transmission=self._datapath+"optics/DISKSPC_FPM_65WA200_360deg_-_FP1res%d_evensamp_D%d_F465.fits.gz"%(fpmres, 2*20*fpmres))
        #elif self.image_mask == 'DISKSPC_F565':
        #    optsys.addImage(transmission=self._datapath+"optics/DISKSPC_FPM_65WA200_360deg_-_FP1res%d_evensamp_D%d_F565.fits.gz"%(fpmres, 2*20*fpmres))
        #elif self.image_mask == 'DISKSPC_F835':
        #    optsys.addImage(transmission=self._datapath+"optics/DISKSPC_FPM_65WA200_360deg_-_FP1res%d_evensamp_D%d_F835.fits.gz"%(fpmres, 2*20*fpmres))
        #elif self.image_mask == 'DISKSPC_F885':
        #    optsys.addImage(transmission=self._datapath+"optics/DISKSPC_FPM_65WA200_360deg_-_FP1res%d_evensamp_D%d_F885.fits.gz"%(fpmres, 2*20*fpmres))

        #optsys.addPupil(name='RELAY')

        occ_box_size = 1.
        mft_optsys = poppy.MatrixFTCoronagraph(optsys, oversample=oversample, occulter_box=occ_box_size)

        return (mft_optsys, trySAM, occ_box_size)

    def _get_aberrations(self):
        """Get the OpticalElement that applies the field-dependent
        optical aberrations. (Called in _getOpticalSystem.)"""
        return None

    def _getFITSHeader(self, result, options):
        """Populate FITS Header keywords"""
        super(WFIRSTInstrument, self)._getFITSHeader(result, options)
