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
            pixx, pixy = (round(midpoint_pixel + local_x * 1e2),
                          round(midpoint_pixel + local_y * 1e2))

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
            for idx in range(22):
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

    def __init__(self, *args, **kwargs):
        super(WFIRSTInstrument, self).__init__(*args, **kwargs)

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
        return self._detectors[self._selected_detector].field_position

    @detector_position.setter
    def detector_position(self, position):
        detector = self._detectors[self._selected_detector]
        detector.field_position = position

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

        # Paths to the two possible pupils. The correct one is selected based on requested
        # wavelengths in _validateConfig()
        self._unmasked_pupil_path = os.path.join(self._WebbPSF_basepath, 'AFTA_WFC_C5_Pupil_Shortwave_Norm_2048px.fits')
        self._masked_pupil_path = os.path.join(self._WebbPSF_basepath, 'AFTA_WFC_C5_Pupil_Mask_Norm_2048px.fits')

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


def show_notebook_interface(instrument):
    # Widget related imports.
    # (Currently not a hard dependency for the full webbpsf package, so we import
    # within the function.)
    import ipywidgets as widgets
    from IPython.display import display, clear_output
    from matplotlib import pyplot as plt

    try:
        import pysynphot
    except ImportError:
        raise ImportError("For now, PySynphot must be installed to use the notebook interface")

    # Clean up some warnings we know about so as not to scare the users
    import warnings
    from matplotlib.cbook import MatplotlibDeprecationWarning
    warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
    warnings.simplefilter('ignore', fits.verify.VerifyWarning)

    def make_binding_for_attribute(attribute):
        def callback(trait_name, new_value):
            setattr(instrument, attribute, new_value)
        return callback

    filter_selection = widgets.ToggleButtons(
        options=instrument.filter_list,
        value=instrument.filter,
        description='Filter:'
    )
    filter_selection.on_trait_change(
        make_binding_for_attribute('filter'),
        name='selected_label'
    )
    display(filter_selection)

    monochromatic_wavelength = widgets.BoundedFloatText(
        value=0.76,
        min=0.6,
        max=2.0,
    )
    monochromatic_wavelength.disabled = True
    monochromatic_toggle = widgets.Checkbox(description='Monochromatic calculation?')

    def update_monochromatic(trait_name, new_value):
        filter_selection.disabled = new_value
        monochromatic_wavelength.disabled = not new_value

    monochromatic_toggle.on_trait_change(update_monochromatic, name='value')

    display(widgets.HTML(value='''<p style="padding: 1em 0;">
    <span style="font-style:italic; font-size:1.0em">
    Monochromatic calculations can be performed for any wavelength in the 0.6 to 2.0 &micro;m range.
    </span></p>'''))  # kludge
    monochromatic_controls = widgets.HBox(children=(
            monochromatic_toggle,
            widgets.HTML(value='<span style="display: inline-block; width: 0.6em;"></span>'),
            monochromatic_wavelength,
            widgets.HTML(value='<span style="display: inline-block; width: 0.25em;"></span> &micro;m '),
    ))
    display(monochromatic_controls)

    display(widgets.HTML(value="<hr>"))

    source_selection = widgets.Select(
        options=poppy.specFromSpectralType('', return_list=True),
        value='G0V',
        description="Source spectrum"
    )
    display(source_selection)
    display(widgets.HTML(value="<hr>"))

    sca_selection = widgets.Dropdown(
        options=instrument.detector_list,
        value=instrument.detector,
        description='Detector:'
    )
    sca_selection.on_trait_change(
        make_binding_for_attribute('detector'),
        name='selected_label'
    )
    display(sca_selection)

    detector_field_points = [
        ('Top left', (4.0, 4092.0)),
        ('Bottom left', (4.0, 4.0)),
        ('Center', (2048.0, 2048.0)),
        ('Top right', (4092.0, 4092.0)),
        ('Bottom right', (4092.0, 4.0)),
    ]
    # enforce ordering of buttons
    detector_field_point_labels = [a[0] for a in detector_field_points]
    detector_field_points = dict(detector_field_points)

    def set_field_position(trait_name, new_value):
        instrument.detector_position = detector_field_points[new_value]

    field_position = widgets.ToggleButtons(options=detector_field_point_labels, value='Center', description='Detector field point:')
    field_position.on_trait_change(set_field_position, name='selected_label')
    display(field_position)

    calculate_button = widgets.Button(
        description="Calculate PSF",
        width='10em',
        color='white',
        background_color='#00c403',
        border_color='#318732'
    )
    display_osys_button = widgets.Button(
        description="Display Optical System",
        width='13em',
        color='white',
        background_color='#005fc4',
        border_color='#224A75'
    )
    clear_button = widgets.Button(
        description="Clear Output",
        width='10em',
        color='white',
        background_color='#ed4747',
        border_color='#911C1C'
    )
    progress = widgets.HTML(value='<progress>')

    OUTPUT_FILENAME = 'psf.fits'
    DOWNLOAD_BUTTON_HTML = """
    <a class="btn btn-info" href="files/{}" target="_blank">
        Download FITS image from last calculation
    </a>
    """
    download_link = widgets.HTML(value=DOWNLOAD_BUTTON_HTML.format(OUTPUT_FILENAME))

    def disp(*args):
        progress.visible = True
        plt.figure(figsize=(12, 8))
        instrument.display()
        progress.visible = None

    def calc(*args):
        progress.visible = True
        if monochromatic_toggle.value is True:
            psf = instrument.calcPSF(
                monochromatic=monochromatic_wavelength.value * 1e-6,
                display=True,
                outfile=OUTPUT_FILENAME,
                clobber=True
            )
        else:
            source = poppy.specFromSpectralType(source_selection.value)
            _log.debug("Got source type {}: {}".format(source_selection.value, source))
            psf = instrument.calcPSF(
                source=source,
                display=True,
                outfile=OUTPUT_FILENAME,
                clobber=True
            )
        fig, (ax_oversamp, ax_detsamp) = plt.subplots(1, 2)
        poppy.display_PSF(psf, ax=ax_oversamp)
        poppy.display_PSF(psf, ax=ax_detsamp, ext='DET_SAMP')
        progress.visible = None
        download_link.visible = True

    def clear(*args):
        clear_output()
        progress.visible = None
        download_link.visible = None

    calculate_button.on_click(calc)
    display_osys_button.on_click(disp)
    clear_button.on_click(clear)
    display(widgets.HTML(value="<br/>"))  # kludge
    buttons = widgets.HBox(children=[calculate_button, display_osys_button, clear_button])
    display(buttons)

    # Insert the progress bar, hidden by default
    display(progress)
    progress.visible = None
    # and the download link
    display(download_link)
    download_link.visible = None
