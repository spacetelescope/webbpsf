from collections import OrderedDict
import os, sys
import astropy.io.fits as fits
from astropy.nddata import NDData
import numpy as np
import matplotlib.pyplot as plt

import scipy.interpolate as sciint

import logging

_log = logging.getLogger('webbpsf')

from . import conf

_DISABLE_FILE_LOGGING_VALUE = 'none'

_Strehl_perfect_cache = {}  # dict for caching perfect images used in Strehl calcs.


### Helper routines for logging: ###

class FilterLevelRange(object):
    def __init__(self, min_level, max_level):
        self.min_level = min_level
        self.max_level = max_level

    def filter(self, record):
        if record.levelno >= self.min_level and record.levelno <= self.max_level:
            return 1
        else:
            return 0


def restart_logging(verbose=True):
    """Restart logging using the same settings as the last WebbPSF
    session, as stored in the configuration system.

    Parameters
    ----------
    verbose : boolean
        Should this function print the new logging targets to
        standard output?
    """

    level = str(conf.logging_level).upper()
    lognames = ['webbpsf', 'poppy']

    root_logger = logging.getLogger()
    root_logger.handlers = []

    if level in ['DEBUG', 'INFO', 'WARN', 'ERROR', 'CRITICAL']:
        level_id = getattr(logging, level)  # obtain one of the DEBUG, INFO, WARN,
        # or ERROR constants
        if verbose:
            print("WebbPSF log messages of level {0} and above will be shown.".format(level))
    elif level == 'NONE':
        root_logger.handlers = []  # n.b. this will clear any handlers other libs/users configured
        return
    else:
        raise ValueError("Invalid logging level: {}".format(level))

    for name in lognames:
        logger = logging.getLogger(name)
        logger.setLevel(level_id)

    # set up screen logging
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.addFilter(FilterLevelRange(
        min_level=logging.DEBUG,
        max_level=logging.INFO
    ))

    stderr_handler = logging.StreamHandler(stream=sys.stderr)
    stderr_handler.addFilter(FilterLevelRange(
        min_level=logging.WARNING,
        max_level=logging.CRITICAL
    ))
    formatter = logging.Formatter(conf.logging_format_screen)
    stderr_handler.setFormatter(formatter)
    stdout_handler.setFormatter(formatter)

    root_logger.addHandler(stdout_handler)
    root_logger.addHandler(stderr_handler)

    if verbose:
        print("WebbPSF log outputs will be directed to the screen.")

    # set up file logging
    filename = conf.logging_filename
    if filename is None or filename.strip().lower() != _DISABLE_FILE_LOGGING_VALUE:
        hdlr = logging.FileHandler(filename)

        formatter = logging.Formatter(conf.logging_format_file)
        hdlr.setFormatter(formatter)

        root_logger.addHandler(hdlr)

        if verbose:
            print("WebbPSF log outputs will also be saved to file {}".format(filename))


def setup_logging(level='INFO', filename=None):
    """Allows selection of logging detail and output locations
    (screen and/or file)

    This is a convenience wrapper to Python's built-in logging package,
    as used by `webbpsf` and `poppy`. By default, this sets up log
    messages to be written to the screen, but the user can also
    request logging to a file.

    Editing the WebbPSF config file to set `autoconfigure_logging = True`
    (and any of the logging settings you wish to persist) instructs
    WebbPSF to apply your settings on import. (This is not
    done by default in case you have configured `logging` yourself
    and don't wish to overwrite your configuration.)

    For more advanced log handling, see the Python logging module's
    own documentation.

    Parameters
    -------------
    level : str
        Name of log output to show. Defaults to 'INFO', set to 'DEBUG'
        for more extensive messages, or to 'WARN' or 'ERROR' for fewer.
    filename : str, optional
        Filename to write the log output to. If not set, output will
        just be displayed on screen. (Default: None)

    Examples
    -----------

    >>> webbpsf.setup_logging(filename='webbpsflog.txt')

    This will save all log messages to 'webbpsflog.txt' in the current
    directory. If you later start another copy of webbpsf in a
    different directory, that session will also write to
    'webbpsflog.txt' in *that* directory. Alternatively you can specify
    a fully qualified absolute path to save all your logs to one
    specific file.

    >>> webbpsf.setup_logging(level='WARN')

    This will show only WARNING or ERROR messages on screen, and not
    save any logs to files at all (since the filename argument is None)
    """

    # implementation note: All this function actually does is apply the
    # defaults into the configuration system, then calls restart_logging to
    # do the actual work.
    level = str(level).upper()

    # The astropy config system will enforce the limited set of values for the logging_level
    # parameter by raising a TypeError on this next line if we feed in an invalid string.
    conf.logging_level = level

    if filename is None:
        # Use the string 'none' as a sentinel value for astropy.config
        filename = _DISABLE_FILE_LOGGING_VALUE

    conf.logging_filename = filename
    restart_logging(verbose=True)


### Helper routines for data handling and system setup: ###

MISSING_WEBBPSF_DATA_MESSAGE = """
 ***********  ERROR  ******  ERROR  ******  ERROR  ******  ERROR  ***********
 *                                                                          *
 *  WebbPSF requires several data files to operate.                         *
 *  These files could not be located automatically at this time, or this    *
 *  version of the software requires a newer set of reference files than    *
 *  you have installed.  For more details see:                              *
 *                                                                          *
 *        https://webbpsf.readthedocs.io/en/stable/installation.html        *
 *                                                                          *
 *  under "Installing the Required Data Files".                             *
 *  WebbPSF will not be able to function properly until the appropriate     *
 *  reference files have been downloaded to your machine and installed.     *
 *                                                                          *
 ****************************************************************************
"""


def get_webbpsf_data_path(data_version_min=None, return_version=False):
    """Get the WebbPSF data path

    Simply checking an environment variable is not always enough, since
    for packaging this code as a Mac .app bundle, environment variables are
    not available since .apps run outside the Terminal or X11 environments.

    Therefore, check first the environment variable WEBBPSF_PATH, and secondly
    check the configuration file in the user's home directory.

    If data_version_min is supplied (as a 3-tuple of integers), it will be
    compared with the version number from version.txt in the WebbPSF data
    package.
    """
    import os
    path_from_config = conf.WEBBPSF_PATH  # read from astropy configuration
    if path_from_config == 'from_environment_variable':
        path = os.getenv('WEBBPSF_PATH')
        if path is None:
            sys.stderr.write(MISSING_WEBBPSF_DATA_MESSAGE)
            raise EnvironmentError("Environment variable $WEBBPSF_PATH is not set!")
    else:
        path = path_from_config

    # at minimum, the path must be a valid directory
    if not os.path.isdir(path):
        sys.stderr.write(MISSING_WEBBPSF_DATA_MESSAGE)
        raise IOError("WEBBPSF_PATH ({}) is not a valid directory path!".format(path))

    if data_version_min is not None:
        # Check if the data in WEBBPSF_PATH meet the minimum data version
        version_file_path = os.path.join(path, 'version.txt')
        try:
            with open(version_file_path) as f:
                version_contents = f.read().strip()
                # keep only first 3 elements for comparison (allows "0.3.4.dev" or similar)
                parts = version_contents.split('.')[:3]
            version_tuple = tuple(map(int, parts))
        except (IOError, ValueError):
            sys.stderr.write(MISSING_WEBBPSF_DATA_MESSAGE)
            raise EnvironmentError(
                "Couldn't read the version number from {}. (Do you need to update the WebbPSF "
                "data? See https://webbpsf.readthedocs.io/en/stable/installation.html#data-install "
                "for a link to the latest version.)".format(version_file_path)
            )

        if not version_tuple >= data_version_min:
            sys.stderr.write(MISSING_WEBBPSF_DATA_MESSAGE)
            raise EnvironmentError(
                "WebbPSF data package has version {cur}, but {min} is needed. "
                "See https://webbpsf.readthedocs.io/en/stable/installation.html#data-install "
                "for a link to the latest version.".format(
                    cur=version_contents,
                    min='{}.{}.{}'.format(*data_version_min)
                )
            )

        if return_version:
            return (path, version_contents)

    return path


DIAGNOSTIC_REPORT = """
OS: {os}
CPU: {cpu}
Python version: {python}
numpy version: {numpy}
scipy version: {scipy}
astropy version: {astropy}
pysynphot version: {pysyn}

numexpr version: {numexpr}
pyFFTW version: {pyfftw}
Anaconda Accelerate version: {accelerate}

poppy version: {poppy}
webbpsf version: {webbpsf}

tkinter version: {tkinter}
wxpython version: {wxpython}


***************************************************************
Floating point type information for numpy.float:
{finfo_float}
Floating point type information for numpy.complex:
{finfo_complex}

***************************************************************
Numpy compilation and linking:

{numpyconfig}

"""


def system_diagnostic():
    """ return various helpful/informative information about the
    current system. For instance versions of python & available packages.

    Mostly undocumented function...
    """

    # There is probably a more clever way to do the following via introspection?

    import platform
    import os
    import poppy
    import numpy
    import scipy

    try:
        from .version import version
    except ImportError:
        version = ''

    try:
        import ttk
        ttk_version = ttk.__version__
    except ImportError:
        ttk_version = 'not found'

    try:
        import wx
        wx_version = wx.__version__
    except ImportError:
        wx_version = 'not found'

    try:
        import pyfftw
        try:
            pyfftw_version = pyfftw.__version__
        except AttributeError:
            # Back compatibility: Handle older versions with nonstandard version attribute name
            pyfftw_version = pyfftw.version

    except ImportError:
        pyfftw_version = 'not found'

    try:
        import pysynphot
        pysynphot_version = pysynphot.__version__
    except ImportError:
        pysynphot_version = 'not found'

    try:
        import astropy
        astropy_version = astropy.__version__
    except ImportError:
        astropy_version = 'not found'

    try:
        import numexpr
        numexpr_version = numexpr.__version__
    except ImportError:
        numexpr_version = 'not found'

    try:
        import accelerate
        accelerate_version = accelerate.__version__
    except ImportError:
        accelerate_version = 'not found'

    try:
        import psutil
        cpu_info = """
  Hardware cores: {hw}
  Logical core: {logical}
  Frequency: {freq} GHz
  Currently {percent}% utilized.
""".format(hw=psutil.cpu_count(logical=False),
           logical=psutil.cpu_count(logical=True),
           freq=psutil.cpu_freq()[0] / 1000,
           percent=psutil.cpu_percent())
    except:
        try:
            import multiprocessing
            cpu_info = "  Cores: {}".format(multiprocessing.cpu_count())
        except:
            cpu_info = "No CPU info available"

    # Get numpy config - the following is a modified version of
    # numpy.__config__.show()

    numpyconfig = ""
    for name, info_dict in numpy.__config__.__dict__.items():
        if name[0] == "_" or type(info_dict) is not type({}): continue
        numpyconfig += name + ":\n"
        if not info_dict:
            numpyconfig += "  NOT AVAILABLE\n"
        for k, v in info_dict.items():
            v = str(v)
            if k == "sources" and len(v) > 200:
                v = v[:60] + " ...\n... " + v[-60:]
            numpyconfig += "    %s = %s\n" % (k, v)

    result = DIAGNOSTIC_REPORT.format(
        os=platform.platform(),
        numpy=numpy.__version__,
        python=sys.version.replace("\n", " "),
        poppy=poppy.__version__,
        webbpsf=version,
        tkinter=ttk_version,
        wxpython=wx_version,
        pyfftw=pyfftw_version,
        pysyn=pysynphot_version,
        astropy=astropy_version,
        finfo_float=numpy.finfo(numpy.float),
        finfo_complex=numpy.finfo(numpy.complex),
        numexpr=numexpr_version,
        scipy=scipy.__version__,
        accelerate=accelerate_version,
        numpyconfig=numpyconfig,
        cpu=cpu_info
    )
    return result


### Helper routines for image manipulation: ###

def measure_strehl(HDUlist_or_filename=None, ext=0, slice=0, center=None, display=True, verbose=True, cache_perfect=False):
    """ Estimate the Strehl ratio for a PSF.

    This requires computing a simulated PSF with the same
    properties as the one under analysis.

    Note that this calculation will not be very accurate unless both PSFs are well sampled,
    preferably several times better than Nyquist. See
    `Roberts et al. 2004 SPIE 5490 <http://adsabs.harvard.edu/abs/2004SPIE.5490..504R>`_
    for a discussion of the various possible pitfalls when calculating Strehl ratios.

    WARNING: This routine attempts to infer how to calculate a perfect reference
    PSF based on FITS header contents. It will likely work for simple direct imaging
    cases with WebbPSF but will not work (yet) for more complicated cases such as
    coronagraphy, anything with image or pupil masks, etc. Code contributions to add
    such cases are welcomed.


    Parameters
    ----------
    HDUlist_or_filename : string
        Either a fits.HDUList object or a filename of a FITS file on disk
    ext : int
        Extension in that FITS file
    slice : int, optional
        If that extension is a 3D datacube, which slice (plane) of that datacube to use
    center : tuple
        center to compute around.  Default is image center. If the center is on the
        crosshairs between four pixels, then the mean of those four pixels is used.
        Otherwise, if the center is in a single pixel, then that pixel is used.
    verbose, display : bool
        control whether to print the results or display plots on screen.

    cache_perfect : bool
        use caching for perfect images? greatly speeds up multiple calcs w/ same config

    Returns
    ---------
    strehl : float
        Strehl ratio as a floating point number between 0.0 - 1.0

    """

    from .webbpsf_core import Instrument
    from poppy import display_psf

    if isinstance(HDUlist_or_filename, str):
        HDUlist = fits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, fits.HDUList):
        HDUlist = HDUlist_or_filename
    else:
        raise ValueError("input must be a filename or HDUlist")

    image = HDUlist[ext].data
    header = HDUlist[ext].header

    if image.ndim >= 3:  # handle datacubes gracefully
        image = image[slice, :, :]

    if center is None:
        # get exact center of image
        # center = (image.shape[1]/2, image.shape[0]/2)
        center = tuple((a - 1) / 2.0 for a in image.shape[::-1])

    # Compute a comparison image
    _log.info("Now computing image with zero OPD for comparison...")
    inst = Instrument(header['INSTRUME'])
    inst.filter = header['FILTER']
    inst.pupilopd = None  # perfect image
    inst.include_si_wfe = False  # perfect image
    inst.pixelscale = header['PIXELSCL'] * header['OVERSAMP']  # same pixel scale pre-oversampling
    cache_key = (header['INSTRUME'], header['FILTER'], header['PIXELSCL'], header['OVERSAMP'], header['FOV'], header['NWAVES'])
    try:
        comparison_psf = _Strehl_perfect_cache[cache_key]
    except KeyError:
        comparison_psf = inst.calc_psf(fov_arcsec=header['FOV'], oversample=header['OVERSAMP'], nlambda=header['NWAVES'])
        if cache_perfect: _Strehl_perfect_cache[cache_key] = comparison_psf

    comparison_image = comparison_psf[0].data

    if (int(center[1]) == center[1]) and (int(center[0]) == center[0]):
        # individual pixel
        meas_peak = image[center[1], center[0]]
        ref_peak = comparison_image[center[1], center[0]]
    else:
        # average across a group of 4
        bot = [int(np.floor(f)) for f in center]
        top = [int(np.ceil(f) + 1) for f in center]
        meas_peak = image[bot[1]:top[1], bot[0]:top[0]].mean()
        ref_peak = comparison_image[bot[1]:top[1], bot[0]:top[0]].mean()
    strehl = (meas_peak / ref_peak)

    if display:
        plt.clf()
        plt.subplot(121)
        display_psf(HDUlist, title="Observed PSF")
        plt.subplot(122)
        display_psf(comparison_psf, title="Perfect PSF")
        plt.gcf().suptitle("Strehl ratio = %.3f" % strehl)

    if verbose:
        print("Measured peak:  {0:.3g}".format(meas_peak))
        print("Reference peak: {0:.3g}".format(ref_peak))
        print("  Strehl ratio: {0:.3f}".format(strehl))

    return strehl


### Helper routines for display customization: ###
# use via poppy's display_annotate feature by assigning these to
# the display_annotate attribute of an OpticalElement class

def annotate_ote_entrance_coords(self, ax):
    """ Draw OTE V frame axes on first optical plane """
    color = 'yellow'
    loc = 3
    ax.arrow(-loc, -loc, .2, 0, color=color, width=0.005)
    ax.arrow(-loc, -loc, 0, .2, color=color, width=0.005)
    ax.text(-loc, -loc + 0.4, '+V3', color=color, size='small',
            horizontalalignment='center', verticalalignment='bottom')
    ax.text(-loc + 0.4, -loc, '+V2', color=color, size='small',
            horizontalalignment='left', verticalalignment='center')


def annotate_sky_pupil_coords(self, ax, show_NE=False, north_angle=45.):
    """ Draw OTE V frame axes projected onto the sky
    Optionally also draw a compass for north and east at some given
    position angle
    """
    color = 'yellow'
    loc = 2.9
    ax.arrow(-loc + 0.5, -loc, -.2, 0, color=color, width=0.005)
    ax.arrow(-loc + 0.5, -loc, 0, .2, color=color, width=0.005)
    ax.text(-loc + 0.5, -loc + 0.3, '+V3 on sky', color=color, size='small',
            horizontalalignment='center', verticalalignment='bottom')
    ax.text(-loc + 0.5 + 0.3, -loc, '+V2 on sky', color=color, size='small',
            horizontalalignment='left', verticalalignment='center')

    if show_NE:
        color2 = 'cyan'
        angle = np.deg2rad(north_angle)  # arbitrary
        dl = 0.3
        dx = np.sin(angle) * dl
        dy = np.cos(angle) * dl
        ax.arrow(-loc + 0.5, -loc, -dx, dy, color=color2, width=0.005)
        ax.arrow(-loc + 0.5, -loc, -dy, -dx, color=color2, width=0.005)
        ax.text(-loc + 0.5 - 2.3 * dx, -loc + 2.3 * dy, 'N', color=color2, size='small',
                horizontalalignment='center', verticalalignment='center')
        ax.text(-loc + 0.5 - 1.3 * dy, -loc - 1.3 * dx, 'E', color=color2, size='small',
                horizontalalignment='center', verticalalignment='center')



def _run_benchmark(timer, iterations=1):
    """ Common benchmarking core. Called from benchmark_imaging and benchmark_coronagraphy
    """
    import poppy
    defaults = (poppy.conf.use_fftw, poppy.conf.use_numexpr, poppy.conf.use_cuda, poppy.conf.use_opencl)

    # Time baseline performance in numpy
    print("Timing performance in plain numpy:")
    poppy.conf.use_fftw, poppy.conf.use_numexpr, poppy.conf.use_cuda, poppy.conf.use_opencl = (False, False, False, False)
    time_numpy = timer.timeit(number=iterations) / iterations
    print("  {:.2f} s".format(time_numpy))

    if poppy.accel_math._FFTW_AVAILABLE:
        print("Timing performance with FFTW:")
        poppy.conf.use_fftw = True
        time_fftw = timer.timeit(number=iterations) / iterations
        print("  {:.2f} s".format(time_fftw))
    else:
        time_fftw = np.NaN

    if poppy.accel_math._NUMEXPR_AVAILABLE:
        print("Timing performance with Numexpr:")
        poppy.conf.use_fftw = False
        poppy.conf.use_numexpr = True
        time_numexpr = timer.timeit(number=iterations) / iterations
        print("  {:.2f} s".format(time_numexpr))
    else:
        time_numexpr = np.NaN

    if poppy.accel_math._CUDA_AVAILABLE:
        print("Timing performance with CUDA + Numexpr:")
        poppy.conf.use_cuda = True
        poppy.conf.use_opencl = False
        time_cuda = timer.timeit(number=iterations) / iterations
        print("  {:.2f} s".format(time_cuda))
    else:
        time_cuda = np.NaN

    if poppy.accel_math._OPENCL_AVAILABLE:
        print("Timing performance with OpenCL + Numexpr:")
        poppy.conf.use_opencl = True
        poppy.conf.use_cuda = False
        time_opencl = timer.timeit(number=iterations) / iterations
        print("  {:.2f} s".format(time_opencl))
    else:
        time_opencl = np.NaN

    poppy.conf.use_fftw, poppy.conf.use_numexpr, poppy.conf.use_cuda, poppy.conf.use_opencl = defaults

    return {'numpy': time_numpy,
            'fftw': time_fftw,
            'numexpr': time_numexpr,
            'cuda': time_cuda,
            'opencl': time_opencl}


def benchmark_imaging(iterations=1, nlambda=1, add_distortion=True):
    """ Performance benchmark function for standard imaging """
    import poppy
    import timeit

    timer = timeit.Timer("psf = nc.calc_psf(nlambda=nlambda, add_distortion={})".format(add_distortion),
                         setup="""
import webbpsf
nc = webbpsf.NIRCam()
nc.filter='F360M'
nlambda={nlambda:d}""".format(nlambda=nlambda))
    print("Timing performance of NIRCam F360M with {} wavelengths, {} iterations".format(nlambda, iterations))


    return _run_benchmark(timer, iterations=iterations)


def benchmark_nircam_coronagraphy(iterations=1, nlambda=1, add_distortion=True):
    """ Performance benchmark function for standard imaging """
    import poppy
    import timeit

    timer = timeit.Timer("psf = nc.calc_psf(nlambda=nlambda, add_distortion={})".format(add_distortion),
                         setup="""
import webbpsf
nc = webbpsf.NIRCam()
nc.filter='F335M'
nc.image_mask='MASK335R'
nc.pupil_mask='MASKRND'
nlambda={nlambda:d}""".format(nlambda=nlambda))
    print("Timing performance of NIRCam MASK335R with {} wavelengths, {} iterations".format(nlambda, iterations))


    return _run_benchmark(timer, iterations=iterations)


def benchmark_miri_coronagraphy(iterations=1, nlambda=1):
    """ Performance benchmark function for standard imaging """
    import poppy
    import timeit

    timer = timeit.Timer("psf = miri.calc_psf(nlambda=nlambda)",
                         setup="""
import webbpsf
miri = webbpsf.MIRI()
miri.filter='F1065C'
miri.image_mask='FQPM1065'
miri.pupil_mask='MASKFQPM'
nlambda={nlambda:d}""".format(nlambda=nlambda))
    print("Timing performance of MIRI F1065C with {} wavelengths, {} iterations".format(nlambda, iterations))


    return _run_benchmark(timer, iterations=iterations)


def combine_docstrings(cls):
    """ Combine the docstrings of a method and earlier implementations of the same method in parent classes """
    for name, func in cls.__dict__.items():

        # Allow users to see the Poppy calc_psf docstring along with the JWInstrument version
        if name == 'calc_psf':
            jwinstrument_class = cls
            spacetelescope_class = cls.__base__

            ind0 = getattr(jwinstrument_class, 'calc_psf').__doc__.index("add_distortion")  # pull the new parameters
            ind1 = getattr(spacetelescope_class, 'calc_psf').__doc__.index("Returns")  # end of parameters

            func.__doc__ = getattr(spacetelescope_class, 'calc_psf').__doc__[0:ind1] + \
                           getattr(jwinstrument_class, 'calc_psf').__doc__[ind0:] + \
                           getattr(spacetelescope_class, 'calc_psf').__doc__[ind1:]

    return cls


def to_griddedpsfmodel(HDUlist_or_filename=None, ext_data=0, ext_header=0):
    """
    Create a photutils GriddedPSFModel object from either a FITS file or
    an HDUlist object. The input must have header keywords "DET_YX{}" and
    "OVERSAMP" (will already be present if psf_grid() is used to create
    the file).

    Parameters
    ----------
    HDUlist_or_filename : HDUList or str
        Either a fits.HDUList object or a filename of a FITS file on disk
    ext_data : int
        Extension of the data in the FITS file
    ext_header : int
        Extension of the header in the FITS file

    Returns
    -------
    model : GriddedPSFModel
        Photutils object with 3D data array and metadata with specified
        grid_xypos and oversampling keys
    """
    try:
        from photutils import GriddedPSFModel
    except ImportError:
        raise ImportError("This method requires photutils >= 0.6")

    if isinstance(HDUlist_or_filename, str):
        HDUlist = fits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, fits.HDUList):
        HDUlist = HDUlist_or_filename
    else:
        raise ValueError('Input must be a filename or HDUlist')

    data = HDUlist[ext_data].data
    header = HDUlist[ext_header].header

    # If there's only 1 PSF and the data is 2D, make the data 3D for photutils can use it
    if len(data.shape) == 2 and len(header['DET_YX*']) == 1:
        data = np.array([data])

    # Check necessary keys are there
    if not any("DET_YX" in key for key in header.keys()):
        raise KeyError("You are missing 'DET_YX{}' keys: which are the detector locations of the PSFs")
    if 'OVERSAMP' not in header.keys():
        raise KeyError("You are missing 'OVERSAMP' key: which is the oversampling factor of the PSFs")

    # Convert header to meta dict
    header = header.copy(strip=True)
    header.pop('COMMENT', None)
    header.pop('', None)
    header.pop('HISTORY', None)
    meta = OrderedDict((a, (b, c)) for (a, b, c) in header.cards)

    ndd = NDData(data, meta=meta, copy=True)

    # Edit meta dictionary for GriddedPSFLibrary specifics
    ndd.meta['grid_xypos'] = [((float(ndd.meta[key][0].split(',')[1].split(')')[0])),
                              (float(ndd.meta[key][0].split(',')[0].split('(')[1])))
                              for key in ndd.meta.keys() if "DET_YX" in key]  # from (y,x) to (x,y)

    if 'oversampling' not in ndd.meta:
        ndd.meta['oversampling'] = ndd.meta['OVERSAMP'][0]  # pull the value

    # Turn all metadata keys into lowercase
    ndd.meta = {key.lower(): ndd.meta[key] for key in ndd.meta}

    # Create model
    model = GriddedPSFModel(ndd)

    return model
