import os
import itertools

import numpy as np
import astropy.convolution
from astropy.io import fits


class CreatePSFLibrary:

    # Class variables for NIRCam short vs long wave information:
    nrca_short_filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M', 'F187N',
                          'F200W', 'F210M', 'F212N']
    nrca_long_filters = ['F250M', 'F277W', 'F300M', 'F322W2', 'F323N', 'F335M', 'F356W', 'F360M', 'F405N', 'F410M',
                         'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']

    nrca_short_detectors = ['NRCA1', 'NRCA2', 'NRCA3', 'NRCA4', 'NRCB1', 'NRCB2', 'NRCB3', 'NRCB4']
    nrca_long_detectors = ['NRCA5', 'NRCB5']

    def __init__(self, instrument, filters="all", detectors="all", num_psfs=16, psf_location=None,
                 use_detsampled_psf=False, save=True, filename=None, overwrite=True, **kwargs):
        """
        Description:
        ------------
        Create a PSF library in the following format:
            For a given instrument, filter, and detector 1 file is produced in the form
            [i, y, x] where i is the PSF position on the detector grid and (y,x)
            is the 2D PSF.

        Parameters:
        -----------
        instrument: instance
            The instance of WebbPSF that is calling this class inside the psf_grid
            method.

        filters: str or list
            Which filter(s) you want to create a library for.

            Can be a string of 1 filter name, a list of filter names (as strings), or
            the default "all" will run through all the filters in the filter_list
            attribute of webbpsf.INSTR(). Spelling/capitalization must match what
            WebbPSF expects. See also special case for NIRCam.

        detectors: str or list
            Which detector(s) you want to create a library for.

            Can be a string of 1 detector name, a list of detector names (as strings), or
            the default "all" will run through all the detectors in the detector_list
            attribute of webbpsf.INSTR(). Spelling/capitalization must match what
            WebbPSF expects. See also special case for NIRCam.

        num_psfs: int
            The total number of fiducial PSFs to be created and saved in the files. This
            number must be a square number. Default is 16.
            E.g. num_psfs = 16 will create a 4x4 grid of fiducial PSFs.

        psf_location: tuple
            If num_psfs = 1, then this is used to set the (y,x) location of that PSF.
            Default is the center point for the detector.

        use_detsampled_psf: bool
            If True, the grid of PSFs returned will be detector sampled (made by binning down the
            oversampled PSF). If False, the PSFs will be oversampled by the factor defined by the
            oversample/detector_oversample/fft_oversample keywords. Default is False.

        save: bool
            True/False boolean if you want to save your file

        filename: str
            If "save" keyword is set to True, your current file will be saved under
            "{filename}_det_filt.fits". Default of None will save it in the current
            directory as: instr_det_filt_fovp#_samp#_npsf#.fits

        overwrite: bool
            True/False boolean to overwrite the output file if it already exists. Default
            is True.

        **kwargs
            This can be used to add any extra arguments to the WebbPSF calc_psf() method
            call.

        Special Case for NIRCam:
        For NIRCam, you can set detectors and filters with multiple options.
        You may set both filters and detectors = "all" just like the other instruments,
        and the short and long wave filter/detectors will be separated and run in the
        correct pairings.
        If you choose only certain filters (either by name or with "shortwave" or
        "longwave" to run all the shortwave/longwave filters), you may set detectors
        to "shortwave" or "longwave" or you can set it to be "all" and the script will
        pull the all possible detectors (ie either all the shortwave or all the longwave
        detectors).
        If you choose individual filters and detectors, they must match in
        shortwave or longwave. Mismatched lists of short and long wave filters and
        detectors will result in an error.

        Returns:
        -------
        Returns or saves 3D fits HDUlist object - 1 per instrument, detector, and filter

        Use:
        ----
        c = CreatePSFLibrary(instrument, filters, detectors, num_psfs, add_distortion,
                             fov_pixels, oversample, save, fileloc, filename, overwrite)
        c.create_files()

        """

        # Pull WebbPSF instance
        self.webb = instrument
        self.instr = instrument.name

        # Set psf_location if not already set
        if psf_location is None:
            psf_location = (int(self.webb._detector_npixels / 2), int(self.webb._detector_npixels / 2))

        # Set the filters and detectors based on the inputs
        self.filter_input = filters
        self.detector_input = detectors

        # A list of filters and a list of list of detectors (1 sublist per filter)
        self.filter_list = self._set_filters()
        self.detector_list = [self._set_detectors(filter) for filter in self.filter_list]

        # Set the locations on the detector of the fiducial PSFs
        self.location_list = self._set_psf_locations(num_psfs, psf_location)

        # For NIRCam: Check if filters/detectors match in terms of if they are longwave/shortwave
        if self.instr == "NIRCam":
            for filt, det_list in zip(self.filter_list, self.detector_list):
                for det in det_list:
                    if "5" in det and filt in CreatePSFLibrary.nrca_short_filters:
                        self._raise_value_error("short filter", det, filt)
                    elif "5" not in det and filt in CreatePSFLibrary.nrca_long_filters:
                        self._raise_value_error("long filter", det, filt)

        # Set PSF attributes for the 3 that will be used before the calc_psf call
        if "add_distortion" in kwargs:
            self.add_distortion = kwargs["add_distortion"]
        else:
            self.add_distortion = True
            kwargs["add_distortion"] = self.add_distortion

        if "oversample" in kwargs:
            self.oversample = kwargs["oversample"]
        elif len({"detector_oversample", "fft_oversample"}.intersection(kwargs.keys())) == 2:
            self.oversample = kwargs["detector_oversample"]
        elif len({"detector_oversample", "fft_oversample"}.intersection(kwargs.keys())) == 1:
            raise ValueError("Must pass either oversample keyword or detector_sample and fft_oversample keywords")
        else:
            self.oversample = 5
            kwargs["oversample"] = self.oversample

        if "fov_pixels" in kwargs:  # fov_pixels overrides fov_arcsec if both set -> same as in calc_psf
            self.fov_pixels = kwargs["fov_pixels"]
        elif "fov_arcsec" in kwargs:
            # From poppy.instrument.get_optical_system() -> line 581
            self.fov_pixels = int(np.round(kwargs["fov_arcsec"] / self.webb.pixelscale))
            if 'parity' in self.webb.options:
                if self.webb.options['parity'].lower() == 'odd' and np.remainder(self.fov_pixels, 2) == 0:
                    self.fov_pixels += 1
                if self.webb.options['parity'].lower() == 'even' and np.remainder(self.fov_pixels, 2) == 1:
                    self.fov_pixels += 1
        else:
            self.fov_pixels = 101
            kwargs["fov_pixels"] = self.fov_pixels

        self.use_detsampled_psf = use_detsampled_psf
        self._kwargs = kwargs

        # Set saving attributes
        self.save = save
        self.overwrite = overwrite
        self.filename = filename

    def _set_filters(self):
        """ Get the list of filters to create PSF library files for """

        # Set filter list to loop over
        if self.filter_input == "all":
            filter_list = self.webb.filter_list
        elif self.filter_input == "shortwave":
            filter_list = CreatePSFLibrary.nrca_short_filters
        elif self.filter_input == "longwave":
            filter_list = CreatePSFLibrary.nrca_long_filters
        elif type(self.filter_input) is str:
            filter_list = self.filter_input.split()
        elif type(self.filter_input) is list:
            filter_list = self.filter_input
        else:
            raise TypeError("Method of setting filters is not valid - see docstring for options")

        # If the user hand chose a filter list, check it's a valid list for the chosen instrument
        if self.filter_input not in ["all", "shortwave", "longwave"]:
            filt = set(filter_list).difference(set(self.webb.filter_list))
            if filt != set(): raise ValueError("Instrument {} doesn't have the filter(s) {}.".format(self.instr, filt))

        return filter_list

    def _set_detectors(self, filter):
        """ Get the list of detectors to include in the PSF library files """

        # Set detector list to loop over
        if self.detector_input == "all":
            if self.instr != "NIRCam":
                detector_list = self.webb.detector_list
            elif self.instr == "NIRCam" and filter in CreatePSFLibrary.nrca_short_filters:
                detector_list = CreatePSFLibrary.nrca_short_detectors
            elif self.instr == "NIRCam" and filter in CreatePSFLibrary.nrca_long_filters:
                detector_list = CreatePSFLibrary.nrca_long_detectors
        elif self.detector_input == "shortwave":
            detector_list = CreatePSFLibrary.nrca_short_detectors
        elif self.detector_input == "longwave":
            detector_list = CreatePSFLibrary.nrca_long_detectors
        elif type(self.detector_input) is str:
            detector_list = self.detector_input.split()
        elif type(self.detector_input) is list:
            detector_list = self.detector_input
        else:
            raise TypeError("Method of setting detectors is not valid - see docstring for options")

        # If the user hand chose a detector list, check it's a valid list for the chosen instrument
        if self.detector_input not in ["all", "shortwave", "longwave"]:
            det = set(detector_list).difference(set(self.webb.detector_list))
            if det != set(): raise ValueError("Instrument {} doesn't have the detector(s) {}.".format(self.instr, det))

        return detector_list

    @staticmethod
    def _raise_value_error(msg_type, det, filt):
        """Raise a specific ValueError based on mis-matched short/long wave detectors/filters"""

        if "short filter" in msg_type.lower():
            message = "You are trying to apply a shortwave filter ({}) to a longwave detector ({}). ".format(filt, det)
        if "long filter" in msg_type.lower():
            message = "You are trying to apply a longwave filter ({}) to a shortwave detector ({}). ".format(filt, det)

        raise ValueError(message + "Please change these entries so the filter falls within the detector band.")

    def _set_psf_locations(self, num_psfs, psf_location):
        """ Set the locations on the detector of the fiducial PSFs"""

        self.num_psfs = num_psfs
        max_size = self.webb._detector_npixels - 1

        if np.sqrt(self.num_psfs).is_integer():
            self.length = int(np.sqrt(self.num_psfs))
        else:
            raise ValueError("You must choose a square number of fiducial PSFs to create (E.g. 9, 16, etc.)")

        # Set the values
        if num_psfs == 1:
            # Want this case to be at the specified location
            location_list = [(psf_location[1], psf_location[0])]  # tuple of (x,y)
        else:
            loc_list = [int(round(num * max_size)) for num in np.linspace(0, 1, self.length, endpoint=True)]
            location_list = list(itertools.product(loc_list, loc_list))  # list of tuples (x,y) (for WebbPSF)

        return location_list

    def create_files(self):
        """
        This method is called in the create_files() method

        For a given instrument, filter, and detector 1 file is produced in the form
        [i, y, x] where i is the PSF position on the detector grid and (y,x)
        is the 2D PSF.

        Returns:
        -------
        This saves out the library files if requested and then returns a list of all the
        hdulist objects created (each in the form of [i, y, x], 1 per instrument/filter/
        detector).

        """

        # Set output mode and extension to use
        if self.use_detsampled_psf is True:
            self.webb.options['output_mode'] = 'Detector sampled image'
            self.oversample = 1
            if self.add_distortion:
                ext = "DET_DIST"
            else:
                ext = "DET_SAMP"
        elif self.use_detsampled_psf is False:
            self.webb.options['output_mode'] = 'Oversampled image'
            if self.add_distortion:
                ext = "OVERDIST"
            else:
                ext = "OVERSAMP"

        print(self.fov_pixels, self.oversample, self.add_distortion, self.use_detsampled_psf)
        print(self._kwargs)
        print(self)

        # Create kernel to smooth pixel based on oversample
        kernel = astropy.convolution.Box2DKernel(width=self.oversample)

        # For every filter
        final_list = []
        for filt, det_list in zip(self.filter_list, self.detector_list):
            print("\nStarting filter: {}".format(filt))

            # Set filter
            self.webb.filter = filt

            # For every detector
            for k, det in enumerate(det_list):
                print("  Running detector: {}".format(det))

                # Create an array to fill ([i, y, x])
                psf_size = self.fov_pixels * self.oversample
                psf_arr = np.empty((self.length**2, psf_size, psf_size))

                self.webb.detector = det

                # For each of the locations on the detector (loc = tuple = (x,y))
                for i, loc in enumerate(self.location_list):
                    self.webb.detector_position = loc  # (X,Y) - line 286 in webbpsf_core.py

                    # Create PSF
                    psf = self.webb.calc_psf(**self._kwargs)

                    # Convolve PSF with a square kernel
                    psf_conv = astropy.convolution.convolve(psf[ext].data, kernel)

                    # Add PSF to 5D array
                    psf_arr[i, :, :] = psf_conv

                # Write header
                header = fits.Header()

                header["INSTRUME"] = (self.instr, "Instrument name")
                header["DETECTOR"] = (det, "Detector name")
                header["FILTER"] = (filt, "Filter name")
                header["PUPILOPD"] = (self.webb.pupilopd, "Pupil OPD source name")

                header["FOVPIXEL"] = (self.fov_pixels, "Field of view in pixels (full array)")
                header["FOV"] = (psf[ext].header["FOV"], "Field of view in arcsec (full array) ")
                header["OVERSAMP"] = (psf[ext].header["OVERSAMP"], "Oversampling factor for FFTs in computation")
                header["DET_SAMP"] = (psf[ext].header["DET_SAMP"], "Oversampling factor for MFT to detector plane")
                header["NWAVES"] = (psf[ext].header["NWAVES"], "Number of wavelengths used in calculation")

                for h, loc in enumerate(self.location_list):  # these were originally written out in (x,y)
                    header["DET_YX{}".format(h)] = (str((loc[1], loc[0])),
                                                    "The #{} PSF's (y,x) detector pixel position".format(h))

                header["NUM_PSFS"] = (self.num_psfs, "The total number of fiducial PSFs")

                # Distortion information
                if self.add_distortion:
                    header["DISTORT"] = (psf[ext].header["DISTORT"], "SIAF distortion coefficients applied")
                    header["SIAF_VER"] = (psf[ext].header["SIAF_VER"], "SIAF PRD version used")

                    for key in list(psf[ext].header.keys()):
                        if "COEF_" in key:
                            header[key] = (psf[ext].header[key], "SIAF distortion coefficient for {}".format(key))

                    if self.instr in ["NIRCam", "NIRISS", "FGS"]:
                        header["ROTATION"] = (psf[ext].header["ROTATION"], "PSF rotated to match detector rotation")

                    if self.instr is "MIRI":
                        header["MIR_DIST"] = (psf[ext].header["MIR_DIST"], "MIRI detector scattering applied")
                        header["KERN_AMP"] = (psf[ext].header["KERN_AMP"],
                                              "Amplitude(A) in kernel function A * exp(-x / B)")
                        header["KERNFOLD"] = (psf[ext].header["KERNFOLD"],
                                              "e - folding length(B) in kernel func A * exp(-x / B)")

                # Pull values from the last made psf
                header["WAVELEN"] = (psf[ext].header["WAVELEN"], "Weighted mean wavelength in meters")
                header["DIFFLMT"] = (psf[ext].header["DIFFLMT"], "Diffraction limit lambda/D in arcsec")
                header["FFTTYPE"] = (psf[ext].header["FFTTYPE"], "Algorithm for FFTs: numpy or fftw")
                header["NORMALIZ"] = (psf[ext].header["NORMALIZ"], "PSF normalization method")
                header["JITRTYPE"] = (psf[ext].header["JITRTYPE"], "Type of jitter applied")
                header["JITRSIGM"] = (psf[ext].header["JITRSIGM"], "Gaussian sigma for jitter [arcsec]")
                header["TEL_WFE"] = (psf[ext].header["TEL_WFE"], "[nm] Telescope pupil RMS wavefront error")

                header["DATE"] = (psf[ext].header["DATE"], "Date of calculation")
                header["AUTHOR"] = (psf[ext].header["AUTHOR"], "username@host for calculation")
                header["VERSION"] = (psf[ext].header["VERSION"], "WebbPSF software version")
                header["DATAVERS"] = (psf[ext].header["DATAVERS"], "WebbPSF reference data files version ")

                # Add descriptor for how the file was made
                header["COMMENT"] = "For a given instrument, filter, and detector 1 file is produced in "
                header["COMMENT"] = "the form [i, y, x] where i is the PSF position on the detector grid "
                header["COMMENT"] = "and (y,x) is the 2D PSF. The order of PSFs can be found under the "
                header["COMMENT"] = "header DET_YX* keywords"

                # Add header labels
                header.insert("INSTRUME", ('', ''))
                header.insert("INSTRUME", ('COMMENT', '/ PSF Library Information'))

                header.insert("NORMALIZ", ('', ''))
                header.insert("NORMALIZ", ('COMMENT', '/ WebbPSF Creation Information'))

                header.insert("DATAVERS", ('COMMENT', '/ File Description'), after=True)
                header.insert("DATAVERS", ('', ''), after=True)

                # Combine the header and data
                hdu = fits.HDUList([fits.PrimaryHDU(psf_arr, header=header)])

                # Write file out
                if self.save:

                    # Set file information
                    if self.filename is None:
                        path = ""

                        # E.g. filename: nircam_nrca1_f090w_fovp1000_samp5_npsf16.fits
                        name = "{}_{}_{}_fovp{}_samp{}_npsf{}.fits".format(self.instr.lower(), det.lower(),
                                                                           filt.lower(), self.fov_pixels,
                                                                           self.oversample, self.num_psfs)
                        file = os.path.join(path, name)
                    else:
                        file = self.filename.split(".")[0]+"_{}_{}.fits".format(det.lower(),filt.lower())

                    print("  Saving file: {}".format(file))

                    hdu.writeto(file, overwrite=self.overwrite)

                final_list.append(hdu)

        # If only 1 hdulist object created, only return that. Else, return list of objects
        if len(self.filter_list) == 1 and len(self.detector_list[0]) == 1:
            return final_list[0]
        else:
            return final_list
