from __future__ import division, print_function, absolute_import, unicode_literals
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
_log = logging.getLogger('webbpsf')

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
            aberration_terms = griddata(self._wavelengths, self._aberration_terms, wavelength, method='linear')
            if np.any(np.isnan(aberration_terms)):
                raise RuntimeError("Attempted to get aberrations at wavelength "
                                   "outside the range of the reference data")
            return aberration_terms

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
    _field_position = None

    def __init__(self, pixel_width, pixel_height,
                 name="Field-dependent Aberration", radius=1.0, oversample=1, interp_order=3):
        self.pixel_width, self.pixel_height = pixel_width, pixel_height
        self.field_position = pixel_width // 2, pixel_height // 2
        self._wavelength_interpolators = {}
        self.pupil_diam = radius * 2.0
        super(FieldDependentAberration, self).__init__(
            name=name,
            verbose=True,
            radius=radius,
            oversample=oversample,
            interp_order=interp_order
        )

    def get_opd(self, wave, units='meters'):
        """Set the Zernike coefficients (for ZernikeWFE.getOPD) based
        on the wavelength of the incoming wavefront and the pixel
        position
        """
        if not isinstance(wave, poppy.Wavefront):
            wavelength = wave
        else:
            wavelength = wave.wavelength
        self.coefficients = wavelength * self.get_aberration_terms(wavelength)
        return super(FieldDependentAberration, self).get_opd(wave, units=units)

    @property
    def field_position(self):
        return self._field_position

    @field_position.setter
    def field_position(self, position):
        """Set the x and y pixel position on the detector for which to
        interpolate aberrations"""
        x_pixel, y_pixel = position
        if x_pixel > self.pixel_width or x_pixel < 0:
            raise ValueError("Requested pixel_x position lies outside "
                             "the detector width ({})".format(x_pixel))
        if y_pixel > self.pixel_height or y_pixel < 0:
            raise ValueError("Requested pixel_y position lies outside "
                             "the detector height ({})".format(y_pixel))

        self._field_position = x_pixel, y_pixel

    def add_field_point(self, x_pixel, y_pixel, interpolator):
        """Supply a wavelength-space interpolator for a pixel position
        on the detector"""
        self._wavelength_interpolators[(x_pixel, y_pixel)] = interpolator

    def get_aberration_terms(self, wavelength):
        """Supply the Zernike coefficients for the aberration based on
        the wavelength and pixel position on the detector"""
        if self.field_position in self._wavelength_interpolators:
            # short path: this is a known point
            interpolator = self._wavelength_interpolators[self.field_position]
            coefficients = interpolator.get_aberration_terms(wavelength)
        else:
            # get aberrations at all field points
            field_points, aberration_terms = [], []
            for field_point_coords, point_interpolator in self._wavelength_interpolators.items():
                field_points.append(field_point_coords)
                aberration_terms.append(point_interpolator.get_aberration_terms(wavelength))
            aberration_array = np.asarray(aberration_terms)
            assert len(aberration_array.shape) == 2, "computed aberration array is not 2D " \
                                                     "(inconsistent number of Zernike terms " \
                                                     "at each point?)"
            coefficients = griddata(
                np.asarray(field_points),
                np.asarray(aberration_terms),
                self.field_position,
                method='linear'
            )
            if np.any(np.isnan(coefficients)):
                raise RuntimeError("Attempted to get aberrations for an out-of-bounds field point")
        if self._omit_piston_tip_tilt:
            _log.debug("Omitting piston/tip/tilt")
            coefficients[:3] = 0.0  # omit piston, tip, and tilt Zernikes
        return coefficients

def _load_wfi_detector_aberrations(filename):
    from astropy.io import ascii
    zernike_table = ascii.read(filename)
    detectors = {}

    def build_detector_from_table(number, zernike_table):
        """Build a FieldDependentAberration optic for a detector using
        Zernikes Z1-Z22 at various wavelengths and field points"""
        single_detector_info = zernike_table[zernike_table['sca'] == number]
        field_points = set(single_detector_info['field_point'])
        interpolators = {}
        detector = FieldDependentAberration(
            4096,
            4096,
            radius=WFIRSTInstrument.PUPIL_RADIUS,
            name="Field Dependent Aberration (SCA{:02})".format(number)
        )
        for field_id in field_points:
            field_point_rows = single_detector_info[single_detector_info['field_point'] == field_id]
            local_x, local_y = field_point_rows[0]['local_x'], field_point_rows[0]['local_y']
            interpolator = build_wavelength_dependence(field_point_rows)

            midpoint_pixel = 4096 / 2
            # (local_x in mm / 10 um pixel size) -> * 1e2
            # local_x and _y range from -20.44 to +20.44, so adding to the midpoint pixel
            # makes sense to place (-20.44, -20.44) at (4, 4)
            pixx, pixy = (round(midpoint_pixel + local_x * 1e2),
                          round(midpoint_pixel + local_y * 1e2))

            detector.add_field_point(pixx, pixy, interpolator)
        return detector

    def build_wavelength_dependence(rows):
        """Build an interpolator object that interpolates Z1-Z22 in
        wavelength space"""
        wavelengths = set(rows['wavelength'])
        interpolator = WavelengthDependenceInterpolator(n_wavelengths=len(wavelengths),
                                                        n_zernikes=22)
        for row in rows:
            z = np.zeros(22)
            for idx in range(22):
                z[idx] = row['Z{}'.format(idx + 1)]
            interpolator.set_aberration_terms(row['wavelength'] * 1e-6, z)

        return interpolator

    detector_ids = set(zernike_table['sca'])
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
    def __init__(self, *args, **kwargs):
        super(WFIRSTInstrument, self).__init__(*args, **kwargs)


    # slightly different versions of the following two functions
    # from the parent superclass
    # in order to interface with the FieldDependentAberration class
    @property
    def detector_position(self):
        """The pixel position in (X, Y) on the detector"""
        return self._detectors[self._detector].field_position

    @detector_position.setter
    def detector_position(self, position):
        # exact copy of superclass function except we save the
        # into a different location.
        try:
            x, y = map(int, position)
        except ValueError:
            raise ValueError("Detector pixel coordinates must be pairs of nonnegative numbers, not {}".format(position))
        if x < 0 or y < 0:
            raise ValueError("Detector pixel coordinates must be nonnegative integers")
        if x > self._detector_npixels - 1 or y > self._detector_npixels - 1:
            raise ValueError("The maximum allowed detector pixel coordinate value is {}".format(self._detector_npixels - 1))

        self._detectors[self._detector].field_position = (int(position[0]),int(position[1]))

    def _get_aberrations(self):
        """Get the OpticalElement that applies the field-dependent
        optical aberrations. (Called in _getOpticalSystem.)"""
        return self._detectors[self._detector]

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

        self._detector_npixels=4096
        self._detectors = _load_wfi_detector_aberrations(os.path.join(self._datapath, 'wim_zernikes_cycle6.csv'))
        assert len(self._detectors.keys()) > 0
        self.detector = 'SCA01'

        # Paths to the two possible pupils. The correct one is selected based on requested
        # wavelengths in _validateConfig()
        self._unmasked_pupil_path = os.path.join(self._WebbPSF_basepath, 'wfc_pupil_rev_mcr.fits')
        self._masked_pupil_path = os.path.join(self._WebbPSF_basepath, 'wfc_pupil_masked_rev_mcr.fits')

        # Flag to en-/disable automatic selection of the appropriate pupil_mask
        self.auto_pupil = True
        self.pupil = self._unmasked_pupil_path
        self.pupilopd = os.path.join(self._WebbPSF_basepath, 'upscaled_HST_OPD.fits')

    def _validateConfig(self, **kwargs):
        """Validates that the WFI is configured sensibly

        This mainly consists of selecting the masked or unmasked pupil
        appropriately based on the wavelengths requested.
        """
        if self.auto_pupil and self.pupil in (self._unmasked_pupil_path, self._masked_pupil_path):
            # Does the wavelength range fit entirely in an unmasked filter?
            wavelengths = np.array(kwargs['wavelengths'])
            wl_min, wl_max = np.min(wavelengths), np.max(wavelengths)
            # test shorter filters first; if wl range fits entirely in one of them, it's not going
            # to be the (masked) wideband filter
            if wl_max <= self.UNMASKED_PUPIL_WAVELENGTH_MAX:
                # use unmasked pupil optic
                self.pupil = self._unmasked_pupil_path
                _log.info("Using the unmasked WFI pupil shape based on wavelengths requested")
            else:
                if wl_max > self.MASKED_PUPIL_WAVELENGTH_MAX:
                    _log.warn("Requested wavelength is > 2e-6 m, defaulting to masked pupil shape")
                # use masked pupil optic
                self.pupil = self._masked_pupil_path
                _log.info("Using the masked WFI pupil shape based on wavelengths requested")
        else:
            # If the user has set the pupil to a custom value, let them worry about the
            # correct shape it should have
            pass
        super(WFI, self)._validateConfig(**kwargs)



