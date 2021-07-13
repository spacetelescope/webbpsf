"""
=================
Roman Instruments
=================

WARNING: This model has not yet been validated against other PSF
         simulations, and uses several approximations (e.g. for
         mirror polishing errors, which are taken from HST).
"""

import os.path
import poppy
import numpy as np

from scipy.interpolate import griddata, RegularGridInterpolator
from astropy.io import fits
import astropy.units as u
import logging

from . import webbpsf_core
from .optics import _fix_zgrid_NaNs


_log = logging.getLogger('webbpsf')
import pprint

GRISM_FILTER = 'G150'
PRISM_FILTER = 'P120'

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
                if isinstance(wavelength, u.Quantity):
                    wavelength = wavelength.to(u.m).value
                wavelength_closest = np.clip(wavelength, np.min(self._wavelengths), np.max(self._wavelengths))
                _log.warn("Attempted to get aberrations at wavelength {:.2g} "
                          "outside the range of the reference data; clipping to closest wavelength {:.2g}".format(
                    wavelength, wavelength_closest))

                aberration_terms = griddata(self._wavelengths, self._aberration_terms, wavelength_closest,
                                            method='linear')
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

    def get_opd(self, wave):
        """Set the Zernike coefficients (for ZernikeWFE.getOPD) based
        on the wavelength of the incoming wavefront and the pixel
        position
        """
        if not isinstance(wave, poppy.Wavefront):
            wavelength = wave
        else:
            wavelength = wave.wavelength
        self.coefficients = wavelength * self.get_aberration_terms(wavelength)
        return super(FieldDependentAberration, self).get_opd(wave)

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
            field_position = tuple(self.field_position)
            coefficients = griddata(
                np.asarray(field_points),
                np.asarray(aberration_terms),
                field_position,
                method='linear'
            )
            if np.any(np.isnan(coefficients)):
                # FIND TWO CLOSEST INPUT GRID POINTS:
                dist = []
                corners = field_points[1:]  # use only the corner points
                for i, ip in enumerate(corners):
                    dist.append(np.sqrt(((ip[0] - field_position[0]) ** 2) + ((ip[1] - field_position[1]) ** 2)))
                min_dist_indx = np.argsort(dist)[:2]  # keep two closest points
                # DEFINE LINE B/W TWO POINTS, FIND ORTHOGONAL LINE AT POINT OF INTEREST,
                # AND FIND INTERSECTION OF THESE TWO LINES.
                x1, y1 = corners[min_dist_indx[0]]
                x2, y2 = corners[min_dist_indx[1]]
                dx = x2 - x1
                dy = y2 - y1
                a = (dy * (field_position[1] - y1) + dx * (field_position[0] - x1)) / (dx * dx + dy * dy)
                closest_interp_point = (x1 + a * dx, y1 + a * dy)
                # INTERPOLATE ABERRATIONS TO CLOSEST INTERPOLATED POINT:
                coefficients = griddata(
                    np.asarray(field_points),
                    np.asarray(aberration_terms),
                    closest_interp_point,
                    method='linear')
                # IF CLOSEST INTERPOLATED POINT IS STILL OUTSIDE THE INPUT GRID,
                # THEN USE NEAREST GRID POINT INSTEAD:
                if np.any(np.isnan(coefficients)):
                    coefficients = aberration_terms[min_dist_indx[0] + 1]
                    _log.warn("Attempted to get aberrations at field point {} which is outside the range "
                              "of the reference data; approximating to nearest input grid point".format(field_position))
                else:
                    _log.warn("Attempted to get aberrations at field point {} which is outside the range "
                              "of the reference data; approximating to nearest interpolated point {}".format(
                        field_position, closest_interp_point))
                assert not np.any(np.isnan(coefficients)), "Could not compute aberration " \
                                                           "at field point {}".format(field_position)
        if self._omit_piston_tip_tilt:
            _log.debug("Omitting piston/tip/tilt")
            coefficients[:3] = 0.0  # omit piston, tip, and tilt Zernikes
        return coefficients


def _load_wfi_detector_aberrations(filename):
    from astropy.io import ascii
    zernike_table = ascii.read(filename, encoding='utf-8-sig')
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
            radius=RomanInstrument.PUPIL_RADIUS,
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


class RomanInstrument(webbpsf_core.SpaceTelescopeInstrument):
    PUPIL_RADIUS = 2.4 / 2.0
    """
    RomanInstrument contains data and functionality common to Roman
    instruments, such as setting the pupil shape
    """
    telescope = "Roman"

    def __init__(self, *args, **kwargs):
        super(RomanInstrument, self).__init__(*args, **kwargs)
        self.options['jitter'] = 'gaussian'
        self.options['jitter_sigma'] = 0.014   # See https://roman.ipac.caltech.edu/sims/Param_db.html#telescope

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
            raise ValueError("Detector pixel coordinates must be pairs of nonnegative numbers, "
                             "not {}".format(position))
        if x < 0 or y < 0:
            raise ValueError("Detector pixel coordinates must be nonnegative integers")
        if x > self._detector_npixels - 1 or y > self._detector_npixels - 1:
            raise ValueError("The maximum allowed detector pixel "
                             "coordinate value is {}".format(self._detector_npixels - 1))

        self._detectors[self._detector].field_position = (int(position[0]), int(position[1]))

    def _get_aberrations(self):
        """Get the OpticalElement that applies the field-dependent
        optical aberrations. (Called in get_optical_system.)"""
        return self._detectors[self._detector]

    def _get_fits_header(self, result, options):
        """Populate FITS Header keywords"""
        super(RomanInstrument, self)._get_fits_header(result, options)
        result[0].header['DETXPIXL'] = (self.detector_position[0],
                                        'X pixel position (for field dependent aberrations)')
        result[0].header['DETYPIXL'] = (self.detector_position[1],
                                        'Y pixel position (for field dependent aberrations)')
        result[0].header['DETECTOR'] = (self.detector, 'Detector selected')


class WFIPupilController:
    """
    This is a helper class for the WFI and is used to swap in
    the correct pupil each time the detector is changed.
    The pupil depends on the selected detector, filter and the
    pupil_mask flag provided by the user.
    The user should not interact with this class directly, only
    through the API provided through the WFI class.
    """

    def __init__(self):
        self._datapath = None
        self._pupil_basepath = None

        self._pupil = None
        self._pupil_mask = None # new for webbpsf 1.0

        # Flag to en-/disable automatic selection of the appropriate pupil_mask
        self._auto_pupil = True

        # Flag to en-/disable automatic selection of the appropriate pupil file
        self._auto_pupil_mask = True

        self.pupil_path_formatters = {
            'skinny': 'RST_WIM_Filter_skinny_{0}.fits.gz',
            'wide': 'RST_WIM_Filter_F184_{0}.fits.gz',
            'grism': 'RST_WSM_Grism_grism_{0}.fits.gz',
            'prism': 'RST_WSM_Prism_prism_{0}.fits.gz'}

    @property
    def pupil(self):
        return self._pupil

    @pupil.setter
    def pupil(self, value):
        raise AttributeError('Pupil cannot be directly specified. '
                             'Use lock_pupil() instead.')

    @property
    def pupil_mask(self):
        """
        pupil_mask types:
        - "skinny"
        - "wide"
        - "prism"
        - "grism"
        """
        return self._pupil_mask

    @pupil_mask.setter
    def pupil_mask(self, name):
        """
        """

        raise AttributeError('Pupil mask cannot be directly specified. '
                             'Use lock_pupil_mask() instead.')

    def _get_filter_mask(self, wfi_filter):
        # FOR NEW PUPIL/MASK EDITS
        wfi_filter = wfi_filter.upper()
        if wfi_filter == GRISM_FILTER:
            return 'grism'
        elif wfi_filter == PRISM_FILTER:
            return 'prism'
        elif wfi_filter in ['F184', 'F213']:
            return 'wide'
        else:
            # this method should only be called after WFI.filter was validated,
            # so we assume all inputs are valid and direct those that don't pass
            # preceding cases to skinny
            return 'skinny'

    def set_base_path(self, datapath):
        """
        Sets the path to the WebbPSF data files.
        This should be set before this class is used.
        Parameters
        ----------
        datapath : string
            Path to WebbPSF-WFI data files
        """
        self._datapath = datapath
        self._pupil_basepath = os.path.join(self._datapath, "pupils")

    def update_pupil(self, filter, detector):
        """
        """
        if not self._auto_pupil:
            _log.info('Automatic pupil selection was locked; '
                      'using user-provided pupil.')
            return

        #if set_path and self._pupil_basepath is None:
        if self._pupil_basepath is None:
           raise Exception('update_pupil called before setting pupil file path')

        pupil_mask = (self._get_filter_mask(filter) if self._auto_pupil_mask
                      else self.pupil_mask)
        # log the else case?
        path_formatter = self.pupil_path_formatters[pupil_mask]

        pupil = os.path.join(self._pupil_basepath,
                             path_formatter.format(detector))

        self._pupil = pupil
        self._pupil_mask = pupil_mask

        _log.info(f"Using {'' if self._auto_pupil_mask else 'locked '}"
                  f"pupil mask '{pupil_mask}' and detector '{detector}'.")

    def lock_pupil(self, pupil_path):
        self._pupil_mask = None
        self._pupil = pupil_path
        self._auto_pupil = False

    def unlock_pupil(self):
        self._auto_pupil = True

    def lock_pupil_mask(self, pupil_mask):
        '''
        only skinny, wide, prism, grism
        '''
        if pupil_mask not in self.pupil_path_formatters.keys():
            raise Exception('invalid pupil mask')
        else:
            self._pupil_mask = pupil_mask
            self._auto_pupil_mask = False

    def unlock_pupil_mask(self):
        # how do we know which mask to reset to?
        self._auto_pupil_mask = True


class WFI(RomanInstrument):
    """
    WFI represents the Roman wide field imager
    for the Roman mission

    WARNING: This model has not yet been validated against other PSF
             simulations, and uses several approximations (e.g. for
             mirror polishing errors, which are taken from HST).
    """

    def __init__(self):
        """
        Initiate WFI
        """
        # pixel scale is from Roman-AFTA SDT report final version (p. 91)
        # https://roman.ipac.caltech.edu/sims/Param_db.html
        pixelscale = 110e-3 # arcsec/px # IS THIS STILL CORRECT?

        # Initialize the pupil controller
        self._pupil_controller = WFIPupilController()

        # Initialize the aberrations for super().__init__
        self._aberration_files = {}
        self._is_custom_aberration = False
        self._current_aberration_file = ""

        super(WFI, self).__init__("WFI", pixelscale=pixelscale)

        self._pupil_controller.set_base_path(self._datapath)

        self.pupil_mask_list = self._pupil_controller.pupil_path_formatters.keys()

        # Define default aberration files for WFI modes
        self._aberration_files = {
            'imaging': os.path.join(self._datapath, 'wim_zernikes_cycle9.csv'),
            'prism': os.path.join(self._datapath,
                                  'wsm_prism_zernikes_cycle9.csv'),
            'grism': os.path.join(self._datapath,
                                  'wsm_grism_zernikes_cycle9.csv'),
            'custom': None}

        # Load default detector from aberration file
        self._detector_npixels = 4096
        self._load_detector_aberrations(self._aberration_files[self.mode])
        self.detector = 'SCA01'

        self.opd_list = [os.path.join(self._WebbPSF_basepath,
                                      'upscaled_HST_OPD.fits')]
        self.pupilopd = self.opd_list[-1]

    def _addAdditionalOptics(self, optsys, **kwargs):
        return optsys, False, None

    def _load_detector_aberrations(self, path):
        """
        Helper function that, given a path to a file containing detector aberrations, loads the Zernike values and
        populates the class' dictator list with `FieldDependentAberration` detectors. This function achieves this by
        calling the `webbpsf.roman._load_wfi_detector_aberrations` function.

        Users should use the `override_aberrations` function to override current aberrations.

        Parameters
        ----------
        path : string
            Path to file containing detector aberrations
        """
        detectors = _load_wfi_detector_aberrations(path)
        assert len(detectors.keys()) > 0

        self._detectors = detectors
        self._current_aberration_file = path

    def _validate_config(self, **kwargs):
        """Validates that the WFI is configured sensibly

        This mainly consists of selecting the masked or unmasked pupil
        appropriately based on the wavelengths requested.
        """
        assert self.filter is not None, 'filter is None'
        assert self.detector is not None, 'detector is None'
        self._update_pupil()

        assert self.pupil is not None, 'pupil is None'
        super(WFI, self)._validate_config(**kwargs)

    def _get_filter_mode(self, wfi_filter):
        """
        Given a filter name, return the WFI mode

        Parameters
        ----------
        wfi_filter : string
            Name of WFI filter

        Returns
        -------
        mode : string
            Returns 'imaging', 'grism' or 'prism' depending on filter.

        Raises
        ------
        ValueError
            If the input filter is not found in the WFI filter list
        """

        wfi_filter = wfi_filter.upper()
        if wfi_filter == GRISM_FILTER:
            return 'grism'
        elif wfi_filter == PRISM_FILTER:
            return 'prism'
        elif wfi_filter in self.filter_list:
            return 'imaging'
        else:
            raise ValueError(f"Instrument {self.name} doesn't have a filter "
                             f"called {wfi_filter}.")

    def _update_pupil(self, filter=None, detector=None):
        if detector is None:
            detector = self.detector
        if filter is None:
            filter = self.filter

        if detector is not None and filter is not None:
            self._pupil_controller.update_pupil(filter=filter,detector=detector)

    @RomanInstrument.detector.setter
    def detector(self, value):
        if value.upper() not in self.detector_list:
            raise ValueError("Invalid detector. Valid detector names are: {}".format(', '.join(self.detector_list)))

        self._detector = value.upper()
        if self._detector is not None:
            self._update_pupil(detector=self._detector)

    @RomanInstrument.filter.setter
    def filter(self, value):

        # Update Filter
        # -------------
        value = value.upper()

        if value not in self.filter_list:
            raise ValueError("Instrument %s doesn't have a filter called %s." % (self.name, value))

        self._filter = value

        # Update Aberrations
        # ------------------
        # Check if _aberration_files has been initiated (not empty) and if aberrations are locked by user
        if self._aberration_files and not self._is_custom_aberration:

            # Identify aberration file for new mode
            mode = self._get_filter_mode(self._filter)
            aberration_file = self._aberration_files[mode]

            # If aberrations are not already loaded for the new mode,
            # load and replace detectors using the new mode's aberration file.
            if not os.path.samefile(self._current_aberration_file,
                                    aberration_file):
                self._load_detector_aberrations(aberration_file)

        # Update Pupil
        # ------------
        if self.detector is not None:
            self._update_pupil(filter=self._filter)

    @property
    def pupil(self):
        return self._pupil_controller.pupil

    @pupil.setter
    def pupil(self, value):
        if self._aberration_files: # check if we've been through super() in init (could be more sophisticated...)
            raise AttributeError('Pupil cannot be directly specified. '
                                 'Use lock_pupil() instead.')
        # ASK ROBEL:
        # else: # SHOULD WE JUST IGNORE super's self.pupil ASSIGNMENT??
        #     self._pupil_controller.pupil = value

    @property
    def pupil_mask(self):
        return self._pupil_controller.pupil_mask

    @pupil_mask.setter
    def pupil_mask(self, name):
        """
        Set the pupil mask

        Parameters
        ----------
        name : string
            Name of setting.
            Settings:
                - "skinny"
                - "wide"
                - "prism"
                - "grism"
        """
        raise AttributeError('Pupil mask cannot be directly specified. '
                             'Use lock_pupil_mask() instead.')

    @property
    def mode(self):
        """Current WFI mode"""
        return self._get_filter_mode(self.filter)

    @mode.setter
    def mode(self, value):
        """Mode is set by changing filters"""
        raise AttributeError("WFI mode cannot be directly specified; "
                             "WFI mode is set by changing filters.")

    def lock_aberrations(self, aberration_path):
        """
        This function loads user provided aberrations from a file and locks this instrument
        to only use the provided aberrations (even if the filter or mode change).
        To release the lock and load the default aberrations, use the `reset_override_aberrations` function.
        To load new user provided aberrations, simply call this function with the new path.

        To load custom aberrations, please provide a csv file containing the detector names,
        field point positions and Zernike values. The file should contain the following column names/values
        (comments in parentheses should not be included):
            - sca (Detector number)
            - wavelength (µm)
            - field_point (filed point number/id for SCA and wavelength, starts with 1)
            - local_x (mm, local detector coords)
            - local_y (mm, local detector coords)
            - global_x (mm, global instrument coords)
            - global_y (mm, global instrument coords)
            - axis_local_angle_x (XAN)
            - axis_local_angle_y (YAN)
            - wfe_rms_waves (nm)
            - wfe_pv_waves (waves)
            - Z1 (Zernike phase NOLL coefficients)
            - Z2 (Zernike phase NOLL coefficients)
            - Z3 (Zernike phase NOLL coefficients)
            - Z4 (Zernike phase NOLL coefficients)
              .
              .
              .

        Please refer to the default aberrations files for examples. If you have the WebbPSF data installed and defined,
        you can get the path to that file by running the following:
            >>> from webbpsf import roman
            >>> wfi = roman.WFI()
            >>> print(wfi._aberration_files["imaging"])

        Warning: You should not edit the default files!
        """
        self._load_detector_aberrations(aberration_path)
        self._aberration_files['custom'] = aberration_path
        self._is_custom_aberration = True

    def unlock_aberrations(self):
         """Release detector aberrations override and load defaults"""
         aberration_path = self._aberration_files[self.mode]
         self._load_detector_aberrations(aberration_path)
         self._aberration_files['custom'] = None
         self._is_custom_aberration = False

    def lock_pupil(self, pupil_path):
        if os.path.isfile(pupil_path):
            self._pupil_controller.lock_pupil(pupil_path)
        else:
            raise FileNotFoundError(f"{pupil_path} not found.")

        _log.warning("you're on your own...")

    def unlock_pupil(self):
        # QUESTION: how do we know which pupil to reset to? make it clear in the docstring that you're unlocking *and* resetting to default.
        # QUESTION: if you unlock_pupil_mask, lock_pupil, unlock_pupil, does the pupil mask also reset?
        self._pupil_controller.unlock_pupil()

        self._update_pupil() # reset pupil

    def lock_pupil_mask(self, pupil_mask):
        '''
        only skinny, wide, prism, grism
        '''
        self._pupil_controller.lock_pupil_mask(pupil_mask)

    def unlock_pupil_mask(self):
        # QUESTION: how do we know which mask to reset to? make it clear in the docstring that you're unlocking *and* resetting to default.
        self._pupil_controller.unlock_pupil_mask()

        self._update_pupil() # reset pupil mask


class CGI(RomanInstrument):
    """
    Roman Coronagraph Instrument

    Simulates the PSF of the Roman coronagraph.

    Current functionality is limited to the Shaped Pupil Coronagraph (SPC)
    observing modes, and these modes are only simulated with static, unaberrated
    wavefronts, without relay optics and without DM control. The design
    respresented here is an approximation to a baseline concept, and will be
    subject to change based on trades studies and technology development.

    Parameters
    ----------
    mode : str
        CGI observing mode. If not specified, the __init__ function
        will set this to a default mode 'CHARSPC_F660'
    pixelscale : float
        Detector pixelscale. If not specified, the pixelscale will default to
        0.02 arcsec for configurations usint the IMAGER camera and 0.025 arcsec
        for the IFS.
    fov_arcsec : float
        Field of view in arcseconds. If not specified, the field of view will
        default to 3.20 arcsec for the IMAGER camera and 1.76 arcsec for the IFS.

    """

    camera_list = ['IMAGER', 'IFS']
    filter_list = ['F660', 'F721', 'F770', 'F890']
    apodizer_list = ['CHARSPC', 'DISKSPC']
    fpm_list = ['CHARSPC_F660_BOWTIE', 'CHARSPC_F770_BOWTIE', 'CHARSPC_F890_BOWTIE', 'DISKSPC_F721_ANNULUS']
    lyotstop_list = ['LS30D88']

    _mode_table = {  # MODE             CAMERA    FILTER  APODIZER   FPM             LYOT STOP
        'CHARSPC_F660': ('IFS', 'F660', 'CHARSPC', 'CHARSPC_F660_BOWTIE', 'LS30D88'),
        'CHARSPC_F770': ('IFS', 'F770', 'CHARSPC', 'CHARSPC_F770_BOWTIE', 'LS30D88'),
        'CHARSPC_F890': ('IFS', 'F890', 'CHARSPC', 'CHARSPC_F890_BOWTIE', 'LS30D88'),
        'DISKSPC_F721': ('IMAGER', 'F721', 'DISKSPC', 'DISKSPC_F721_ANNULUS', 'LS30D88')}

    def __init__(self, mode=None, pixelscale=None, fov_arcsec=None, apply_static_opd=False):
        super(CGI, self).__init__("CGI", pixelscale=pixelscale)

        self._detector_npixels = 1024
        self._detectors = {camera: 'placeholder' for camera in self.camera_list}

        self.pupil_mask_list = self.lyotstop_list  # alias for use in webbpsf_core
        self.image_mask_list = self.fpm_list  # alias for use in webbpsf_core
        self.pupil = os.path.join(self._WebbPSF_basepath, 'AFTA_CGI_C5_Pupil_onax_256px_flip.fits')
        if apply_static_opd:
            self.pupilopd = os.path.join(self._WebbPSF_basepath, 'CGI', 'OPD', 'CGI_static_OPD.fits')
        else:
            self.pupilopd = None
        self.aberration_optic = None
        self.options = {'force_coron': True}
        # Allow the user to pre-emptively override the default instrument FoV and pixel scale
        if fov_arcsec is not None:
            self.fov_arcsec = fov_arcsec
            self._override_fov = True
        else:
            self._override_fov = False
        if pixelscale is not None:
            self._pixelscale = pixelscale
            self._override_pixelscale = True
        else:
            self._override_pixelscale = False

        if mode is None:
            self.print_mode_table()
            _log.info("Since the mode was not specified at instantiation, defaulting to CHARSPC_F660")
            self.mode = 'CHARSPC_F660'
        else:
            self.mode = mode

    @property
    def camera(self):
        """Currently selected camera name"""
        return self._camera

    @camera.setter
    def camera(self, value):
        value = value.upper()  # force to uppercase
        if value not in self.camera_list:
            raise ValueError("Instrument {0} doesn't have a camera called {1}.".format(self.name, value))
        self._camera = value
        if value == 'IMAGER':
            if not hasattr(self, 'fov_arcsec') or not self._override_fov:
                self.fov_arcsec = 3.2
            if not hasattr(self, 'pixelscale') or not self._override_pixelscale:
                self.pixelscale = 0.020  # Nyquist at 465 nm
        else:  # default to 'IFS'
            if not hasattr(self, 'fov_arcsec') or not self._override_fov:
                self.fov_arcsec = 2 * 0.82  # 2015 SDT report, Section 3.4.1.1.1:
                                            # IFS has 76 lenslets across the (2 x 0.82) arcsec FoV.
            if not hasattr(self, 'pixelscale') or not self._override_pixelscale:
                self.pixelscale = 0.025  # Nyquist at 600 nm

    # for CGI, there is one detector per camera and it should be set automatically.
    @property
    def detector(self):
        return self.camera

    @detector.setter
    def detector(self, value):
        raise RuntimeError("Can't set detector directly for CGI; set camera instead.")

    @property
    def filter(self):
        """Currently selected filter name"""
        return self._filter

    @filter.setter
    def filter(self, value):
        value = value.upper()  # force to uppercase
        if value not in self.filter_list:
            raise ValueError("Instrument {0} doesn't have a filter called {1}.".format(self.name, value))
        self._filter = value

    @property
    def apodizer(self):
        """Currently selected apodizer name"""
        return self._apodizer

    @apodizer.setter
    def apodizer(self, value):
        value = value.upper()  # force to uppercase
        if value not in self.apodizer_list:
            raise ValueError("Instrument {0} doesn't have a apodizer called {1}.".format(self.name, value))
        self._apodizer = value
        if value == 'DISKSPC':
            self._apodizer_fname = \
                os.path.join(self._datapath, "optics/DISKSPC_SP_256pix.fits.gz")
        else:  # for now, default to CHARSPC
            self._apodizer_fname = \
                os.path.join(self._datapath, "optics/CHARSPC_SP_256pix.fits.gz")

    @property
    def fpm(self):
        """Currently selected FPM name"""
        return self._fpm

    @fpm.setter
    def fpm(self, value):
        value = value.upper()  # force to uppercase
        if value not in self.fpm_list:
            raise ValueError("Instrument {0} doesn't have a FPM called {1}.".format(self.name, value))
        self._fpm = value
        if value.startswith('DISKSPC'):
            self._fpmres = 3
            self._owa = 20.
            self._Mfpm = int(np.ceil(self._fpmres * self._owa))
            self._fpm_fname = \
                os.path.join(self._datapath,
                             "optics/DISKSPC_FPM_65WA200_360deg_-_FP1res{0:d}_evensamp_D{1:03d}_{2:s}.fits.gz".format(
                                 self._fpmres, 2 * self._Mfpm, self.filter))
        else:
            self._fpmres = 4
            self._owa = 9.
            self._Mfpm = int(np.ceil(self._fpmres * self._owa))
            self._fpm_fname = \
                os.path.join(self._datapath,
                             "optics/CHARSPC_FPM_25WA90_2x65deg_-_FP1res{0:d}_evensamp_D{1:03d}_{2:s}.fits.gz".format(
                                 self._fpmres, 2 * self._Mfpm, self.filter))

    @property
    def lyotstop(self):
        """Currently selected Lyot stop name"""
        return self._lyotstop

    @lyotstop.setter
    def lyotstop(self, value):
        # preserve case for this one since we're used to that with the lyot mask names
        if value not in self.lyotstop_list:
            raise ValueError("Instrument {0} doesn't have a Lyot mask called {1}.".format(self.name, value))
        self._lyotstop = value
        self._lyotstop_fname = os.path.join(self._datapath, "optics/SPC_LS_30D88_256pix.fits.gz")

    @property
    def mode_list(self):
        """Available Observation Modes"""
        keys = self._mode_table.keys()
        keys = sorted(keys)
        return keys

    # mode works differently since it's a meta-property that affects the other ones:
    @property
    def mode(self):
        """Currently selected mode name"""
        for modename, settings in self._mode_table.items():
            if (self.camera == settings[0].upper() and self.filter == settings[1].upper() and
                    self.apodizer == settings[2].upper() and self.fpm == settings[3].upper() and
                    self.lyotstop == settings[4]):
                return modename
        return 'Custom'

    @mode.setter
    def mode(self, value):
        if value not in self.mode_list:
            raise ValueError("Instrument {0} doesn't have a mode called {1}.".format(self.name, value))
        settings = self._mode_table[value]
        self.camera = settings[0]
        self.filter = settings[1]
        self.apodizer = settings[2]
        self.fpm = settings[3]
        self.lyotstop = settings[4]
        _log.info('Set the following optical configuration:')
        _log.info('camera = {0}, filter = {1}, apodizer = {2}, fpm = {3}, lyotstop = {4}'.format(\
                  self.camera, self.filter, self.apodizer, self.fpm, self.lyotstop))

    def print_mode_table(self):
        """Print the table of observing mode options and their associated optical configuration"""
        _log.info("Printing the table of Roman CGI observing modes supported by WebbPSF.")
        _log.info("Each is defined by a combo of camera, filter, apodizer, "
                  "focal plane mask (FPM), and Lyot stop settings:")
        _log.info(pprint.pformat(self._mode_table))

    @property
    def detector_position(self):
        """The pixel position in (X, Y) on the detector"""
        return 512, 512

    @detector_position.setter
    def detector_position(self, position):
        raise RuntimeError("Detector position not adjustable for CGI")

    def _validate_config(self, **kwargs):
        super(CGI, self)._validate_config(**kwargs)

    def _addAdditionalOptics(self, optsys, oversample=4):
        """Add coronagraphic or spectrographic optics for Roman CGI."""

        trySAM = False

        if ('pupil_shift_x' in self.options and self.options['pupil_shift_x'] != 0) or \
                ('pupil_shift_y' in self.options and self.options['pupil_shift_y'] != 0):
            shift = (self.options['pupil_shift_x'], self.options['pupil_shift_y'])
        else:
            shift = None

        # Add the shaped pupil apodizer
        optsys.add_pupil(transmission=self._apodizer_fname, name=self.apodizer, shift=None)

        # Add the FPM
        optsys.add_image(transmission=self._fpm_fname, name=self.fpm)

        # Add Lyot stop
        self.pupil_mask = self.lyotstop
        optsys.add_pupil(transmission=self._lyotstop_fname, name=self.lyotstop, shift=shift)

        # Cast as MatrixFTCoronagraph; this configures the detector
        occ_box_size = 1.
        mft_optsys = poppy.MatrixFTCoronagraph(optsys, oversample=oversample, occulter_box=occ_box_size)

        return mft_optsys, trySAM, occ_box_size

    def _get_aberrations(self):
        """Get the OpticalElement that applies the field-dependent
        optical aberrations. (Called in get_optical_system.)"""
        return None

    def _get_fits_header(self, result, options):
        """Populate FITS Header keywords"""
        super(RomanInstrument, self)._get_fits_header(result, options)
        pupil_hdr = fits.getheader(self.pupil)
        apodizer_hdr = fits.getheader(self._apodizer_fname)
        fpm_hdr = fits.getheader(self._fpm_fname)
        lyotstop_hdr = fits.getheader(self._lyotstop_fname)

        result[0].header.set('MODE', self.mode, comment='Observing mode')
        result[0].header.set('CAMERA', self.camera, comment='Imager or IFS')
        result[0].header.set('APODIZER', self.apodizer, comment='Apodizer')
        result[0].header.set('APODTRAN', os.path.basename(self._apodizer_fname),
                             comment='Apodizer transmission')
        result[0].header.set('PUPLSCAL', apodizer_hdr['PUPLSCAL'],
                             comment='Apodizer pixel scale in m/pixel')
        result[0].header.set('PUPLDIAM', apodizer_hdr['PUPLDIAM'],
                             comment='Full apodizer array size, incl padding.')
        result[0].header.set('FPM', self.fpm, comment='Focal plane mask')
        result[0].header.set('FPMTRAN', os.path.basename(self._fpm_fname),
                             comment='FPM transmission')
        result[0].header.set('FPMSCAL', fpm_hdr['PIXSCALE'], comment='FPM spatial sampling, arcsec/pix')
        result[0].header.set('LYOTSTOP', self.lyotstop, comment='Lyot stop')
        result[0].header.set('LSTRAN', os.path.basename(self._lyotstop_fname),
                             comment='Lyot stop transmission')
        result[0].header.set('PUPLSCAL', lyotstop_hdr['PUPLSCAL'],
                             comment='Lyot stop pixel scale in m/pixel')
        result[0].header.set('PUPLDIAM', lyotstop_hdr['PUPLDIAM'],
                             comment='Lyot stop array size, incl padding.')
