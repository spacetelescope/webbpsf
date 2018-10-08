import os
import itertools

import numpy as np
import astropy.convolution
from astropy.io import fits


class CreatePSFLibrary:
    """
    Class Description:
    -----------------
    Class to create a PSF library in the following format:
        For a given instrument, 1 file per filter in the form [SCA, j, i, y, x] where
        (j,i) is the PSF position on the detector grid (integer positions) and (y,x)
        is the 2D PSF.

    This class must be run separately for each instrument.

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

    Parameters:
    -----------
    instrument: str
        The name of the instrument you want to run. Can be any capitalization. Can
        only run 1 instrument at a time. Right now this class is only set up for NIRCam,
        NIRISS, and FGS (they are 2048x2048)

    filters: str or list
        Which filter(s) you want to create a library for.

        Can be a string of 1 filter name, a list of filter names (as strings), or
        the default "all" will run through all the filters in the filter_list
        attribute of webbpsf.INSTR(). Spelling/capitalization must match what
        webbpsf expects. See also special case for NIRCam. Default is "all"

    detectors: str or list
        Which detector(s) you want to create a library for.

        Can be a string of 1 detector name, a list of detector names (as strings), or
        the default "all" will run through all the detectors in the detector_list
        attribute of webbpsf.INSTR(). Spelling/capitalization must match what
        webbpsf expects. See also special case for NIRCam. Default is "all"

    num_psfs: int
        The total number of fiducial PSFs to be created and saved in the files. This
        number must be a square number. Default is 16.
        E.g. num_psfs = 16 will have the class create a 4x4 grid of fiducial PSFs.

    psf_location: tuple
        If num_psfs = 1, then this is used to set the (y,x) location of that PSF.
        Default is (1024,1024).

    add_distortion: bool
        If True, the PSF will have distortions applied: the geometric distortion from
        the detectors (using data from SIAF) and the rotation of the detectors with
        respect to the focal plane. Default is True.

    fov_pixels: int
        The field of view in undersampled detector pixels used by WebbPSF when
        creating the PSFs. Default is 101 pixels.

    oversample: int
        The oversampling factor used by WebbPSF when creating the PSFs. Default is 5.

    opd_type: str
        The type of OPD map you would like to use to create the PSFs. Options are
        "predicted" or "requirements" where the predicted map is of the expected
        WFE and the requirements map is slightly more conservative (has slightly
        larger WFE). Default is "requirements"

    opd_number: int
        The realization of the OPD map pulled from the OPD file. Options are an
        integer from 0 to 9, one for each of the 10 Monte Carlo realizations of
        the telescope included in the OPD map file. Default is 0.

    save: bool
        True/False boolean if you want to save your file

    fileloc: str
        Where to save your file if "save" keyword is set to True. Default of None
        will save in the current directory

    filename: str
        The name to save your current file under if "save" keyword is set to True.
        Default of None will save it in the form: "INSTR_FILT_fovp####_samp#_npsf##.fits"

    overwrite: bool
        True/False boolean to overwrite the output file if it already exists.

    **kwargs
        This can be used to add any extra arguments to the webbpsf calc_psf() method
        call.

    Use:
    ----
    c = CreatePSFLibrary(instrument, filters, detectors, num_psfs, add_distortion,
                         fov_pixels, oversample, save, fileloc, filename, overwrite)
    c.create_files()

    nis = CreatePSFLibrary("NIRISS") # will run all filters/detectors
    nis.create_files()

    """

    # Class variables for NIRCam short vs long wave information:
    nrca_short_filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M', 'F187N',
                          'F200W', 'F210M', 'F212N']
    nrca_long_filters = ['F250M', 'F277W', 'F300M', 'F322W2', 'F323N', 'F335M', 'F356W', 'F360M', 'F405N', 'F410M',
                         'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']

    nrca_short_detectors = ['NRCA1', 'NRCA2', 'NRCA3', 'NRCA4', 'NRCB1', 'NRCB2', 'NRCB3', 'NRCB4']
    nrca_long_detectors = ['NRCA5', 'NRCB5']

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
        """ Set the locations on the detector of the fiducial PSFs. Assumes a 2048x2048 detector"""

        # The locations these PSFs should be centered on for a 2048x2048 detector
        self.num_psfs = num_psfs

        if np.sqrt(self.num_psfs).is_integer():
            self.length = int(np.sqrt(self.num_psfs))
        else:
            raise ValueError("You must choose a square number of fiducial PSFs to create (E.g. 9, 16, etc.)")

        # Set the values
        if num_psfs == 1:
            # Want this case to be at the specified location
            ij_list = [(0, 0)]
            loc_list = [psf_location[1], psf_location[0]]  # list of x,y location
            location_list = [(psf_location[1], psf_location[0])]  # tuple of (x,y)
        else:
            ij_list = list(itertools.product(range(self.length), range(self.length)))
            loc_list = [int(round(num * 2047)) for num in np.linspace(0, 1, self.length, endpoint=True)]
            location_list = list(itertools.product(loc_list, loc_list))  # list of tuples (x,y) (for webbpsf)

        return ij_list, loc_list, location_list

    def __init__(self, webbinst, filters="all", detectors="all", num_psfs=16, psf_location=(1024, 1024),
                 add_distortion=True, fov_pixels=101, oversample=5, opd_type="requirements", opd_number=0,
                 save=True, fileloc=None, filename=None, overwrite=True,
                 **kwargs):

        # Pull webbpsf instance and instrument name
        self.webb = webbinst
        self.instr = webbinst.name

        # Set the filters and detectors based on the inputs
        self.filter_input = filters
        self.detector_input = detectors

        # A list of filters and a list of list of detectors (1 sublist per filter)
        self.filter_list = self._set_filters()
        self.detector_list = [self._set_detectors(filter) for filter in self.filter_list]

        # Set the locations on the detector of the fiducial PSFs
        self.ij_list, self.loc_list, self.location_list = self._set_psf_locations(num_psfs, psf_location)

        # For NIRCam: Check if filters/detectors match in terms of if they are longwave/shortwave
        if self.instr == "NIRCam":
            for filt, det_list in zip(self.filter_list, self.detector_list):
                for det in det_list:
                    if "5" in det and filt in CreatePSFLibrary.nrca_short_filters:
                        self._raise_value_error("short filter", det, filt)
                    elif "5" not in det and filt in CreatePSFLibrary.nrca_long_filters:
                        self._raise_value_error("long filter", det, filt)

        # Set PSF attributes
        self.add_distortion = add_distortion
        self.fov_pixels = fov_pixels
        self.oversample = oversample
        self.opd_type = opd_type
        self.opd_number = opd_number
        self._kwargs = kwargs

        # Set saving attributes
        self.save = save
        self.overwrite = overwrite
        self.fileloc = fileloc
        self.filename = filename

    def create_files(self):
        """
        This method is called in the create_files() method

        For a given instrument, 1 file per filter in the form [SCA, j, i, y, x] where
        (j,i) is the PSF position on the detector grid (integer positions) and (y,x)
        is the 2D PSF.

        Returns:
        -------
        This saves out the library files if requested and then returns a list of all the
        hdulist objects created (each in the form of [SCA, j, i, y, x], 1 per filter
        requested).

        """

        # Set extension to read based on distortion choice
        if self.add_distortion: ext = "OVERDIST"
        else: ext = "OVERSAMP"

        # Create kernel to smooth pixel based on oversample
        kernel = astropy.convolution.Box2DKernel(width=self.oversample)

        # Set output mode
        self.webb.options['output_mode'] = 'Oversampled Image'

        # Set OPD Map (pull most recent version with self.webb.opd_list call) - always predicted then requirements
        if self.opd_type.lower() == "requirements":
            opd = self.webb.opd_list[1]
        elif self.opd_type.lower() == "predicted":
            opd = self.webb.opd_list[0]
        self.webb.pupilopd = (opd, self.opd_number)

        # For every filter
        final_list = []
        for filt, det_list in zip(self.filter_list, self.detector_list):
            print("\nStarting filter: {}".format(filt))

            # Set filter
            self.webb.filter = filt

            # For every detector
            for k, det in enumerate(det_list):
                print("  Running detector: {}".format(det))

                # Create an array to fill ([SCA, j, i, y, x])
                psf_size = self.fov_pixels * self.oversample
                psf_arr = np.empty((self.length, self.length, psf_size, psf_size))

                self.webb.detector = det

                # For each of the 9 locations on the detector (loc = tuple = (x,y))
                for i, loc in enumerate(self.location_list):
                    print(loc)
                    self.webb.detector_position = loc  # (X,Y) - line 286 in webbpsf_core.py

                    # Create PSF
                    print(self.webb.filter, self.webb.detector, self.webb.detector_position)
                    psf = self.webb.calc_psf(add_distortion=self.add_distortion,
                                             fov_pixels=self.fov_pixels, oversample=self.oversample, **self._kwargs)

                    # Convolve PSF with a square kernel
                    psf_conv = astropy.convolution.convolve(psf[ext].data, kernel)

                    # Add PSF to 5D array
                    psf_arr[i, :, :] = psf_conv

                # Write header
                header = fits.Header()

                header["INSTRUME"] = (self.instr, "Instrument name")
                header["DETECTOR"] = (det, "Detector name")
                header["FILTER"] = (filt, "Filter name")
                header["PUPILOPD"] = (self.webb.pupilopd[0], "Pupil OPD source name")
                header["OPD_REAL"] = (self.webb.pupilopd[1], "Pupil OPD source realization from file")

                header["FOVPIXEL"] = (self.fov_pixels, "Field of view in pixels (full array)")
                header["FOV"] = (psf[ext].header["FOV"], "Field of view in arcsec (full array) ")
                header["OVERSAMP"] = (self.oversample, "Oversampling factor for FFTs in computation")
                header["NWAVES"] = (psf[ext].header["NWAVES"], "Number of wavelengths used in calculation")

                for h, loc in enumerate(self.location_list):  # these were originally written out in (i,j) and (x,y)
                    header["DET_YX{}".format(h)] = (str((loc[1], loc[0])),
                                                    "The #{} PSF's (y,x) detector pixel position".format(h))

                header["NUM_PSFS"] = (self.num_psfs, "The total number of fiducial PSFs")

                # The range of location values
                if self.num_psfs == 1:
                    # In this case, loc_list is the single x and y value (x may not equal y)
                    header["I0_X"] = (self.loc_list[0], "The x pixel value for i=0 (AXIS4)")
                    header["J0_Y"] = (self.loc_list[1], "The y pixel value for j=0 (AXIS3)")
                else:
                    last = len(self.loc_list) - 1
                    header["I0_X"] = (self.loc_list[0], "The x pixel value for i=0 (AXIS4)")
                    header["I{}_X".format(last)] = (self.loc_list[-1],
                                                    "The x pixel value for i={} (final value; AXIS4)".format(last))
                    header["J0_Y"] = (self.loc_list[0], "The y pixel value for j=0 (AXIS3)")
                    header["J{}_Y".format(last)] = (self.loc_list[-1],
                                                    "The y pixel value for j={} (final value; AXIS3)".format(last))

                # Distortion information
                if self.add_distortion:
                    header["ROTATION"] = (psf[ext].header["ROTATION"], "PSF rotated to match detector rotation")
                    header["DISTORT"] = (psf[ext].header["DISTORT"], "SIAF distortion coefficients applied")
                    header["SIAF_VER"] = (psf[ext].header["SIAF_VER"], "SIAF PRD version used")

                    for key in list(psf[ext].header.keys()):
                        if "COEF_" in key:
                            header[key] = (psf[ext].header[key], "SIAF distortion coefficient for {}".format(key))

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
                header["COMMENT"] = "For a given instrument, 1 file per filter in the form [SCA, j, i, y, x]"
                header["COMMENT"] = "where (j,i) is the PSF position on the detector grid (integer "
                header["COMMENT"] = "positions) and (y,x) is the 2D PSF. The order of the detectors can be "
                header["COMMENT"] = "found under the  header DETNAME* keywords and the order of the fiducial "
                header["COMMENT"] = "PSFs ((j,i) and (y,x)) under the header DET_JI*/DET_YX* keywords"

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
                    if self.fileloc is None:
                        #self.fileloc = os.path.expandvars('$MIRAGE_DATA/{}/test_webbpsf_library'.format(self.instr.lower()))
                        self.fileloc = "/Users/sosborne/Desktop/"

                    if self.filename is None:
                        # E.g. filename: nircam_nrca1_f090w_fovp1000_samp5_npsf16.fits
                        name = "{}_{}_{}_fovp{}_samp{}_npsf{}.fits".format(self.instr.lower(), det.lower(),
                                                                           filt.lower(), self.fov_pixels,
                                                                           self.oversample, self.num_psfs)
                        self.filepath = os.path.join(self.fileloc, name)
                    else:
                        self.filepath = os.path.join(self.fileloc, self.filename)

                    print("  Saving file: {}".format(self.filepath))

                    hdu.writeto(self.filepath, overwrite=self.overwrite)

                # Create something to return
                final_list.append(hdu)

        return final_list