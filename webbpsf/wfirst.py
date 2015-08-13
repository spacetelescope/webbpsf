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
        if not isinstance(wave, poppy.Wavefront):
            wavelength = wave
        else:
            wavelength = wave.wavelength
        self.coefficients = wavelength * self.get_aberration_terms(wavelength)
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
            for idx in xrange(22):
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
    def __init__(self, *args, **kwargs):
        super(WFIRSTInstrument, self).__init__(*args, **kwargs)
        self.pupil = poppy.FITSOpticalElement(
            transmission=os.path.join(self._WebbPSF_basepath, 'AFTA-WFS_C5_Pupil_Mask.fits')
        )
        self.pupilopd = None  # until we have some OPD maps and a FITS pupil of the right shape


class WFIRSTImager(WFIRSTInstrument):
    """
    WFIRSTImager represents to the to-be-named wide field imager
    for the WFIRST mission

    WARNING: No realistic wavefront error map was available for WFIRST at release time.
             This assumes a perfect telescope!
    """
    def __init__(self):
        scale = 110e-3  # arcsec/px, WFIRST-AFTA SDT report v2 (p. 58)
        super(WFIRSTImager, self).__init__("WFIRSTImager", pixelscale=scale)
        self._apertures = _load_wfi_aberration_apertures(os.path.join(self._datapath, 'zernikes.csv'))
        self._aperture_name = self.aperture_list[0]
        assert len(self._apertures.keys()) > 0

    def _validate_config(self):
        return True

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

    # TODO jlong de-duplicate this functionality by providing a better hook in SpaceTelescopeInstrument to add the field-dependent aberrations
    def _getOpticalSystem(self, fft_oversample=2, detector_oversample=None, fov_arcsec=2, fov_pixels=None, options=None):
        """ Return an OpticalSystem instance corresponding to the instrument as currently configured.

        When creating such an OpticalSystem, you must specify the parameters needed to define the
        desired sampling, specifically the oversampling and field of view.

        Parameters
        ----------

        fft_oversample : int
            Oversampling factor for intermediate plane calculations. Default is 2
        detector_oversample: int, optional
            By default the detector oversampling is equal to the intermediate calculation oversampling.
            If you wish to use a different value for the detector, set this parameter.
            Note that if you just want images at detector pixel resolution you will achieve higher fidelity
            by still using some oversampling (i.e. *not* setting `oversample_detector=1`) and instead rebinning
            down the oversampled data.
        fov_pixels : float
            Field of view in pixels. Overrides fov_arcsec if both set.
        fov_arcsec : float
            Field of view, in arcseconds. Default is 2

        Returns
        -------
        osys : poppy.OpticalSystem
            an optical system instance representing the desired configuration.

        """

        _log.info("Creating optical system model:")

        if options is None:
            options = self.options
        if detector_oversample is None:
            detector_oversample = fft_oversample

        _log.debug("Oversample: %d  %d " % (fft_oversample, detector_oversample))
        optsys = poppy.OpticalSystem(
            name='{telescope}+{instrument}'.format(telescope=self.telescope, instrument=self.name),
            oversample=fft_oversample
        )
        if 'source_offset_r' in options:
            optsys.source_offset_r = options['source_offset_r']
        if 'source_offset_theta' in options:
            optsys.source_offset_theta = options['source_offset_theta']


        #---- set pupil OPD
        if isinstance(self.pupilopd, str):  # simple filename
            opd_map = self.pupilopd if os.path.exists( self.pupilopd) else os.path.join(self._datapath, "OPD",self.pupilopd)
        elif hasattr(self.pupilopd, '__getitem__') and isinstance(self.pupilopd[0], six.string_types): # tuple with filename and slice
            opd_map =  (self.pupilopd[0] if os.path.exists( self.pupilopd[0]) else os.path.join(self._datapath, "OPD",self.pupilopd[0]), self.pupilopd[1])
        elif isinstance(self.pupilopd, (fits.HDUList, poppy.OpticalElement)):
            opd_map = self.pupilopd # not a path per se but this works correctly to pass it to poppy
        elif self.pupilopd is None:
            opd_map = None
        else:
            raise TypeError("Not sure what to do with a pupilopd of that type:"+str(type(self.pupilopd)))

        #---- set pupil intensity
        if self.pupil is None:
            raise RuntimeError("The pupil shape must be specified in the "
                               "instrument class or by setting self.pupil")
        if isinstance(self.pupil, poppy.OpticalElement):
            # supply to POPPY as-is
            pupil_optic = optsys.addPupil(self.pupil)
        else:
            # wrap in an optic and supply to POPPY
            if isinstance(self.pupil, str): # simple filename
                if os.path.exists(self.pupil):
                    pupil_transmission = self.pupil
                else:
                    pupil_transmission = os.path.join(
                        self._WebbPSF_basepath,
                        self.pupil
                    )
            elif isinstance(self.pupil, fits.HDUList):
                # POPPY can use self.pupil as-is
                pupil_transmission = self.pupil
            else:
                raise TypeError("Not sure what to do with a pupil of "
                                "that type: {}".format(type(self.pupil)))
            #---- apply pupil intensity and OPD to the optical model
            pupil_optic = optsys.addPupil(
                name='{} Pupil'.format(self.telescope),
                transmission=pupil_transmission,
                opd=opd_map,
                opdunits='micron',
                rotation=self._rotation
            )
        self.pupil_radius = pupil_optic.pupil_diam / 2.0

        #---- Add defocus if requested
        if 'defocus_waves' in options:
           defocus_waves = options['defocus_waves']
           defocus_wavelength = float(options['defocus_wavelength']) if 'defocus_wavelength' in options else 2.0e-6
           _log.info("Adding defocus of %d waves at %.2f microns" % (defocus_waves, defocus_wavelength *1e6))
           lens = poppy.ThinLens(
               name='Defocus',
               nwaves=defocus_waves,
               reference_wavelength=defocus_wavelength,
               radius=self.pupil_radius
           )
           optsys.addPupil(optic=lens)

        #---- Create and insert virtual optic for field dependent aberration
        if 'pixel_position' in options:
            x_pixel, y_pixel = options['pixel_position']
            aperture_interpolator = self._apertures[self._aperture_name]
            aperture_interpolator.set_field_position(x_pixel, y_pixel)
            optsys.addPupil(aperture_interpolator)

        #--- add the detector element.
        if fov_pixels is None:
            if not np.isscalar(fov_arcsec): fov_arcsec = np.asarray(fov_arcsec) # cast to ndarray if 2D
            fov_pixels = np.round(fov_arcsec/self.pixelscale)
            if 'parity' in options:
                if options['parity'].lower() == 'odd'  and np.remainder(fov_pixels,2)==0: fov_pixels +=1
                if options['parity'].lower() == 'even' and np.remainder(fov_pixels,2)==1: fov_pixels +=1
        else:
            pass

        optsys.addDetector(self.pixelscale, fov_pixels = fov_pixels, oversample = detector_oversample, name=self.name+" detector")

        return optsys