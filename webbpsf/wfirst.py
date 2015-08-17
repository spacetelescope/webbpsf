"""
==============================
WFIRST Instruments (pre-alpha)
==============================

WARNING: No realistic wavefront error map was available for WFIRST at release time.
         This assumes a perfect telescope!
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
    def __init__(self, n_wavelengths=16, n_zernikes=22):
        self._n_wavelengths = n_wavelengths
        self._n_zernikes = n_zernikes
        self._aberration_terms = np.zeros((n_wavelengths, n_zernikes), dtype=np.float64)
        self._wavelengths = []

    def set_aberration_terms(self, wavelength, zernike_array):
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
        # return array of length n_zernikes interpolated for this wavelength
        if wavelength in self._wavelengths:
            # aberration known exactly for this wavelength
            aberration_row_idx = self._wavelengths.index(wavelength)
            return self._aberration_terms[aberration_row_idx]
        else:
            # we have to interpolate @ this wavelength
            return griddata(self._wavelengths, self._aberration_terms, wavelength, method='linear')

    def get_zernike_coefficient(self, wavelength, zernike_subscript):
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
    def __init__(self, pixel_width, pixel_height,
                 name="Field-dependent Aberration", radius=1.0, oversample=1, interp_order=3,
                 omit_piston_tip_tilt=True):
        self.pixel_width, self.pixel_height = pixel_width, pixel_height
        self.omit_piston_tip_tilt = omit_piston_tip_tilt
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
        if not isinstance(wave, poppy.Wavefront):
            wavelength = wave
        else:
            wavelength = wave.wavelength
        self.coefficients = wavelength * self.get_aberration_terms(wavelength)
        if self.omit_piston_tip_tilt:
            self.coefficients[:3] = 0  # zero out Z1-3 since they don't add anything useful
        return super(FieldDependentAberration, self).getPhasor(wave)

    def set_field_position(self, x_pixel, y_pixel):
        if x_pixel > self.pixel_width or x_pixel < 0:
            raise ValueError("Requested pixel_x position lies outside "
                             "the aperture width ({})".format(x_pixel))
        if y_pixel > self.pixel_height or y_pixel < 0:
            raise ValueError("Requested pixel_y position lies outside "
                             "the aperture height ({})".format(y_pixel))

        self.x_pixel, self.y_pixel = x_pixel, y_pixel

    def add_field_point(self, x_pixel, y_pixel, interpolator):
        self._wavelength_interpolators[(x_pixel, y_pixel)] = interpolator

    def get_aberration_terms(self, wavelength):
        # short path: this is a known point
        if (self.x_pixel, self.y_pixel) in self._wavelength_interpolators:
            interpolator = self._wavelength_interpolators[(self.x_pixel, self.y_pixel)]
            return interpolator.get_aberration_terms(wavelength)
        # get aberrations at all field points
        field_points, aberration_terms = [], []
        for field_point_coords, point_interpolator in self._wavelength_interpolators.iteritems():
            field_points.append(field_point_coords)
            aberration_terms.append(point_interpolator.get_aberration_terms(wavelength))
        aberration_array = np.asarray(aberration_terms)
        assert len(aberration_array.shape) == 2, "computed aberration array is not 2D " \
                                                 "(inconsistent number of Zernike terms " \
                                                 "at each point?)"
        computed_aberration = griddata(
            np.asarray(field_points),
            np.asarray(aberration_terms),
            (self.x_pixel, self.y_pixel),
            method='linear'
        )
        return computed_aberration

def _load_wfi_aberration_apertures(filename):
    from astropy.io import ascii
    zernike_table = ascii.read(filename)
    apertures = {}

    def build_detector_from_table(number, zernike_table):
        single_detector_info = zernike_table[zernike_table['Cnfg#'] == number]
        field_points = set(single_detector_info['Field'])
        interpolators = {}
        aperture = FieldDependentAberration(4096, 4096, radius=2.36)
        for field_id in field_points:
            field_point_rows = single_detector_info[single_detector_info['Field'] == field_id]
            local_x, local_y = field_point_rows[0]['Local_x'], field_point_rows[0]['Local_y']
            interpolator = build_wavelength_dependence(field_point_rows)

            midpoint_pixel = 4096 / 2
            # (local_x in mm / 10 um pixel size) -> * 1e2
            pixx, pixy = midpoint_pixel + local_x * 1e2, midpoint_pixel + local_y * 1e2

            aperture.add_field_point(pixx, pixy, interpolator)
        return aperture

    def build_wavelength_dependence(rows):
        wavelengths = set(rows['Wave(um)'])
        interpolator = WavelengthDependenceInterpolator(n_wavelengths=len(wavelengths),
                                                        n_zernikes=22)
        for row in rows:
            z = np.zeros(22)
            for idx in xrange(3, 22):  # omit piston, tip, tilt for WFI
                z[idx] = row['Z{}'.format(idx + 1)]
            interpolator.set_aberration_terms(row['Wave(um)'] * 1e-6, z)

        return interpolator

    aperture_ids = set(zernike_table['Cnfg#'])
    for apid in aperture_ids:
        apertures["SCA{:02}".format(apid)] = build_detector_from_table(apid, zernike_table)

    return apertures

class WFIRSTInstrument(webbpsf_core.SpaceTelescopeInstrument):
    """
    WFIRSTInstrument contains data and functionality common to WFIRST
    instruments, such as setting the pupil shape

    WARNING: No realistic wavefront error map was available for WFIRST at release time.
             This assumes a perfect telescope!
    """
    telescope = "WFIRST"
    _apertures = {}
    """
    Dictionary mapping aperture names to FieldDependentAberration optics.

    (Subclasses must populate this at `__init__`.)
    """
    _aperture_name = None
    """
    The name of the currently selected aperture. Must be a key in _apertures, as validated by the
    `aperture` property setter.

    (Subclasses must populate this at `__init__`.)
    """

    def __init__(self, *args, **kwargs):
        super(WFIRSTInstrument, self).__init__(*args, **kwargs)
        self.pupil = os.path.join(self._WebbPSF_basepath, 'AFTA_WFC_C5_Pupil_Shortwave_Norm_2048px.fits')
        self.pupilopd = os.path.join(self._WebbPSF_basepath, 'upscaled_HST_OPD.fits')

        # n.b. WFIRSTInstrument subclasses must set these
        self._apertures = {}
        self._aperture_name = None

    @property
    def aperture(self):
        return self._aperture_name

    @aperture.setter
    def aperture(self, value):
        if value not in self.aperture_list:
            raise ValueError("Invalid aperture. Valid aperture names are: {}".format(', '.join(self.aperture_list)))
        self._aperture_name = value

    @property
    def aperture_list(self):
        return sorted(self._apertures.keys())

    def get_aberrations(self):
        if 'pixel_position' in self.options:
            x_pixel, y_pixel = self.options['pixel_position']
            aperture_interpolator = self._apertures[self._aperture_name]
            aperture_interpolator.set_field_position(x_pixel, y_pixel)
            _log.info("Setting field position to ({}, {}) on {}".format(
                x_pixel, y_pixel, self._aperture_name
            ))
            return aperture_interpolator
        else:
            return None

class WFIRSTImager(WFIRSTInstrument):
    """
    WFIRSTImager represents to the to-be-named wide field imager
    for the WFIRST mission

    WARNING: No realistic wavefront error map was available for WFIRST at release time.
             This assumes a perfect telescope!
    """
    UNMASKED_PUPIL_WAVELENGTH_MIN, UNMASKED_PUPIL_WAVELENGTH_MAX = 0.760e-6, 1.454e-6
    MASKED_PUPIL_WAVELENGTH_MIN, MASKED_PUPIL_WAVELENGTH_MAX = 1.380e-6, 2.000e-6
    def __init__(self):
        scale = 110e-3  # arcsec/px, WFIRST-AFTA SDT report v2 (p. 58)
        super(WFIRSTImager, self).__init__("WFIRSTImager", pixelscale=scale)
        self._apertures = _load_wfi_aberration_apertures(os.path.join(self._datapath, 'zernikes.csv'))
        self._aperture_name = self.aperture_list[0]
        assert len(self._apertures.keys()) > 0

        # Paths to the two possible pupils. The correct one is selected based on requested
        # wavelengths in _validate_config()
        self._unmasked_pupil_path = os.path.join(self._WebbPSF_basepath, 'AFTA_WFC_C5_Pupil_Shortwave_Norm_2048px.fits')
        self._masked_pupil_path = os.path.join(self._WebbPSF_basepath, 'AFTA_WFC_C5_Pupil_Mask_Norm_2048px.fits')

    def _validate_config(self, **kwargs):
        if self.pupil in (self._unmasked_pupil_path, self._masked_pupil_path):
            # Does the wavelength range fit entirely in an unmasked filter?
            wavelengths = np.array(kwargs['wavelengths'])
            wl_min, wl_max = np.min(wavelengths), np.max(wavelengths)
            # test shorter filters first; if wl range fits entirely in one of them, it's not going
            # to be the (masked) wideband filter
            if wl_min >= UNMASKED_PUPIL_WAVELENGTH_MIN and wl_max <= UNMASKED_PUPIL_WAVELENGTH_MAX:
                # use unmasked pupil optic
                self.pupil = self._unmasked_pupil_path
                _log.info("Using the unmasked WFI pupil shape based on wavelengths requested")
            elif wl_min >= MASKED_PUPIL_WAVELENGTH_MIN and wl_max <= MASKED_PUPIL_WAVELENGTH_MAX:
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
        return super(WFIRSTImager, self)._validateConfig(**kwargs)
