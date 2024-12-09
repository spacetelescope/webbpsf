import logging
import os
import sys
import warnings
from collections import OrderedDict

import astropy.io.fits as fits
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import scipy
from astropy.nddata import NDData

from . import conf

_log = logging.getLogger('webbpsf')


_DISABLE_FILE_LOGGING_VALUE = 'none'

_Strehl_perfect_cache = {}  # dict for caching perfect images used in Strehl calcs.


# Helper routines for logging: ###


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
            print('WebbPSF log messages of level {0} and above will be shown.'.format(level))
    elif level == 'NONE':
        root_logger.handlers = []  # n.b. this will clear any handlers other libs/users configured
        return
    else:
        raise ValueError('Invalid logging level: {}'.format(level))

    for name in lognames:
        logger = logging.getLogger(name)
        logger.setLevel(level_id)

    # set up screen logging
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.addFilter(FilterLevelRange(min_level=logging.DEBUG, max_level=logging.INFO))

    stderr_handler = logging.StreamHandler(stream=sys.stderr)
    stderr_handler.addFilter(FilterLevelRange(min_level=logging.WARNING, max_level=logging.CRITICAL))
    formatter = logging.Formatter(conf.logging_format_screen)
    stderr_handler.setFormatter(formatter)
    stdout_handler.setFormatter(formatter)

    root_logger.addHandler(stdout_handler)
    root_logger.addHandler(stderr_handler)

    if verbose:
        print('WebbPSF log outputs will be directed to the screen.')

    # set up file logging
    filename = conf.logging_filename
    if filename is None or filename.strip().lower() != _DISABLE_FILE_LOGGING_VALUE:
        hdlr = logging.FileHandler(filename)

        formatter = logging.Formatter(conf.logging_format_file)
        hdlr.setFormatter(formatter)

        root_logger.addHandler(hdlr)

        if verbose:
            print('WebbPSF log outputs will also be saved to file {}'.format(filename))


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
    ----------
    level : str
        Name of log output to show. Defaults to 'INFO', set to 'DEBUG'
        for more extensive messages, or to 'WARN' or 'ERROR' for fewer.
    filename : str, optional
        Filename to write the log output to. If not set, output will
        just be displayed on screen. (Default: None)

    Examples
    --------

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


# Helper routines for data handling and system setup: ###

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


def auto_download_webbpsf_data():
    import os
    import tarfile
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from urllib.request import urlretrieve

    # Create a default directory for the data files
    default_path = Path.home() / "data" / "webbpsf-data"
    default_path.mkdir(parents=True, exist_ok=True)

    os.environ["WEBBPSF_PATH"] = str(default_path)

    # Download the data files if the directory is empty
    if not any(default_path.iterdir()):
        warnings.warn("WebbPSF data files not found in default location, attempting to download them now...")

        with TemporaryDirectory() as tmpdir:
            # Download the data files to a temporary directory
            url = "https://stsci.box.com/shared/static/qxpiaxsjwo15ml6m4pkhtk36c9jgj70k.gz"
            filename = Path(tmpdir) / "webbpsf-data-LATEST.tar.gz"
            urlretrieve(url, filename)

            # Extract the tarball
            with tarfile.open(filename, "r:gz") as tar:
                tar.extractall(default_path.parent, filter="fully_trusted")

        if not any(default_path.iterdir()):
            raise IOError(f"Failed to get and extract WebbPSF data files to {default_path}")

    return default_path


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
    from pathlib import Path

    path_from_config = conf.WEBBPSF_PATH  # read from astropy configuration
    if path_from_config == 'from_environment_variable':
        path = os.getenv('WEBBPSF_PATH')
        if path is None:
            message = (
                'Environment variable $WEBBPSF_PATH is not set!\n'
                f'{MISSING_WEBBPSF_DATA_MESSAGE} searching default location..s'
            )
            warnings.warn(message)
            path = auto_download_webbpsf_data()
    else:
        path = path_from_config

    path = Path(path)

    # at minimum, the path must be a valid directory
    if not path.is_dir():
        raise IOError(f'WEBBPSF_PATH ({path}) is not a valid directory path!\n{MISSING_WEBBPSF_DATA_MESSAGE}')

    if data_version_min is not None:
        # Check if the data in WEBBPSF_PATH meet the minimum data version
        version_file_path = path / 'version.txt'
        try:
            with open(version_file_path) as f:
                version_contents = f.read().strip()
                # keep only first 3 elements for comparison (allows "0.3.4.dev" or similar)
                parts = version_contents.split('.')[:3]
            version_tuple = tuple(map(int, parts))
        except (IOError, ValueError):
            raise EnvironmentError(
                f"Couldn't read the version number from {version_file_path}. (Do you need to update the WebbPSF "
                'data? See https://webbpsf.readthedocs.io/en/stable/installation.html#data-install '
                'for a link to the latest version.)'
                f'\n{MISSING_WEBBPSF_DATA_MESSAGE}'
            )

        if not version_tuple >= data_version_min:
            min_ver = '{}.{}.{}'.format(*data_version_min)
            raise EnvironmentError(
                f'WebbPSF data package has version {version_contents}, but {min_ver} is needed. '
                'See https://webbpsf.readthedocs.io/en/stable/installation.html#data-install '
                'for a link to the latest version.'
                f'\n{MISSING_WEBBPSF_DATA_MESSAGE}'
            )

        if return_version:
            return str(path), version_contents

    return str(path)


DIAGNOSTIC_REPORT = """
OS: {os}
CPU: {cpu}
Python version: {python}
numpy version: {numpy}
scipy version: {scipy}
astropy version: {astropy}
stsynphot version: {stsyn}
pysynphot version: {pysyn}

numexpr version: {numexpr}
pyFFTW version: {pyfftw}
Anaconda Accelerate version: {accelerate}

poppy version: {poppy}
webbpsf version: {webbpsf}

tkinter version: {tkinter}
wxpython version: {wxpython}


***************************************************************
Floating point type information for numpy.double:
{finfo_float}
Floating point type information for numpy.cdouble:
{finfo_complex}

***************************************************************
Numpy compilation and linking:

{numpyconfig}

"""


def get_pupil_mask(npix=1024, label_segments=False):
    """Utility function to easily retrieve the pupil mask file for a given size.

    Parameters
    ----------
    npix : float
       Number of pixels desired in the pupil array. Must be 256, 512, 1024, 2048, etc.
    label_segments: bool
       Return map with segment IDs 1 through 18, rather than just 1 or 0 for in pupil or not.

    """
    basename = f'JWpupil_segments_RevW_npix{npix}.fits.gz' if label_segments else f'jwst_pupil_RevW_npix{npix}.fits.gz'

    fullname = os.path.join(get_webbpsf_data_path(), basename)

    if not os.path.exists(fullname):
        fullname = fullname.replace('.fits.gz', '.fits')
    try:
        data = fits.getdata(fullname)
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f'No pupil file found for npix={npix}. Check your WebbPSF data directory, '
            + "and make sure you're using a value which is a power of 2, at least 256"
        ) from e
    return data


def system_diagnostic():
    """return various helpful/informative information about the
    current system. For instance versions of python & available packages.

    Mostly undocumented function...
    """

    # There is probably a more clever way to do the following via introspection?

    import platform

    import numpy
    import poppy
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
        import stsynphot

        stsynphot_version = stsynphot.__version__
    except ImportError:
        stsynphot_version = 'not found'

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
""".format(
            hw=psutil.cpu_count(logical=False),
            logical=psutil.cpu_count(logical=True),
            freq=psutil.cpu_freq()[0] / 1000,
            percent=psutil.cpu_percent(),
        )
    except ImportError:
        try:
            import multiprocessing

            cpu_info = '  Cores: {}'.format(multiprocessing.cpu_count())
        except ImportError:
            cpu_info = 'No CPU info available'

    # Get numpy config - the following is a modified version of
    # numpy.__config__.show()

    numpyconfig = ''
    for name, info_dict in numpy.__config__.__dict__.items():
        if name[0] == '_' or not isinstance(info_dict, dict):
            continue
        numpyconfig += name + ':\n'
        if not info_dict:
            numpyconfig += '  NOT AVAILABLE\n'
        for k, v in info_dict.items():
            v = str(v)
            if k == 'sources' and len(v) > 200:
                v = v[:60] + ' ...\n... ' + v[-60:]
            numpyconfig += '    %s = %s\n' % (k, v)

    result = DIAGNOSTIC_REPORT.format(
        os=platform.platform(),
        numpy=numpy.__version__,
        python=sys.version.replace('\n', ' '),
        poppy=poppy.__version__,
        webbpsf=version,
        tkinter=ttk_version,
        wxpython=wx_version,
        pyfftw=pyfftw_version,
        stsyn=stsynphot_version,
        pysyn=pysynphot_version,
        astropy=astropy_version,
        finfo_float=numpy.finfo(numpy.double),
        finfo_complex=numpy.finfo(numpy.cdouble),
        numexpr=numexpr_version,
        scipy=scipy.__version__,
        accelerate=accelerate_version,
        numpyconfig=numpyconfig,
        cpu=cpu_info,
    )
    return result


# Helper routines for image manipulation: ###


def rms(opd, mask):
    """Compute RMS of an OPD over some given masked area"""
    return np.sqrt((opd[(mask != 0) & np.isfinite(opd)] ** 2).mean())


def measure_strehl(HDUlist_or_filename=None, ext=0, slice=0, center=None, display=True, verbose=True, cache_perfect=False):
    """Estimate the Strehl ratio for a PSF.

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
    -------
    strehl : float
        Strehl ratio as a floating point number between 0.0 - 1.0

    """

    from poppy import display_psf

    from .webbpsf_core import instrument

    if isinstance(HDUlist_or_filename, str):
        HDUlist = fits.open(HDUlist_or_filename)
    elif isinstance(HDUlist_or_filename, fits.HDUList):
        HDUlist = HDUlist_or_filename
    else:
        raise ValueError('input must be a filename or HDUlist')

    image = HDUlist[ext].data
    header = HDUlist[ext].header

    if image.ndim >= 3:  # handle datacubes gracefully
        image = image[slice, :, :]

    if center is None:
        # get exact center of image
        # center = (image.shape[1]/2, image.shape[0]/2)
        center = tuple((a - 1) / 2.0 for a in image.shape[::-1])

    # Compute a comparison image
    _log.info('Now computing image with zero OPD for comparison...')
    inst = instrument(header['INSTRUME'])
    inst.filter = header['FILTER']
    inst.pupilopd = None  # perfect image
    inst.include_si_wfe = False  # perfect image
    inst.pixelscale = header['PIXELSCL'] * header['OVERSAMP']  # same pixel scale pre-oversampling
    cache_key = (
        header['INSTRUME'],
        header['FILTER'],
        header['PIXELSCL'],
        header['OVERSAMP'],
        header['FOV'],
        header['NWAVES'],
    )
    try:
        comparison_psf = _Strehl_perfect_cache[cache_key]
    except KeyError:
        comparison_psf = inst.calc_psf(fov_arcsec=header['FOV'], oversample=header['OVERSAMP'], nlambda=header['NWAVES'])
        if cache_perfect:
            _Strehl_perfect_cache[cache_key] = comparison_psf

    comparison_image = comparison_psf[0].data

    if (int(center[1]) == center[1]) and (int(center[0]) == center[0]):
        # individual pixel
        meas_peak = image[center[1], center[0]]
        ref_peak = comparison_image[center[1], center[0]]
    else:
        # average across a group of 4
        bot = [int(np.floor(f)) for f in center]
        top = [int(np.ceil(f) + 1) for f in center]
        meas_peak = image[bot[1]: top[1], bot[0]: top[0]].mean()
        ref_peak = comparison_image[bot[1]: top[1], bot[0]: top[0]].mean()
    strehl = meas_peak / ref_peak

    if display:
        plt.clf()
        plt.subplot(121)
        display_psf(HDUlist, title='Observed PSF')
        plt.subplot(122)
        display_psf(comparison_psf, title='Perfect PSF')
        plt.gcf().suptitle('Strehl ratio = %.3f' % strehl)

    if verbose:
        print('Measured peak:  {0:.3g}'.format(meas_peak))
        print('Reference peak: {0:.3g}'.format(ref_peak))
        print('  Strehl ratio: {0:.3f}'.format(strehl))

    return strehl


def rescale_interpolate_opd(array, newdim):
    """Interpolates & rescales an input 2D array to any given size.
    Uses scipy.interpolate.RectBivariateSpline

    Parameters
    ----------
    array: float
         input array to interpolate
    newdim: int
         new size of the 2D square array (newdim x newdim)

    Returns
    -------
    newopd: new array interpolated to (newdim x newdim)

    """

    dim = array.shape[0]

    xmax, ymax = dim / 2, dim / 2
    x = np.arange(-xmax, xmax, 1)
    y = np.arange(-ymax, ymax, 1)
    X, Y = np.meshgrid(x, y)

    interp_spline = scipy.interpolate.RectBivariateSpline(y, x, array)

    dx, dy = float(dim) / float(newdim), float(dim) / float(newdim)

    x2 = np.arange(-xmax, xmax, dx)
    y2 = np.arange(-ymax, ymax, dy)
    # X2, Y2 = np.meshgrid(x2, y2)
    newopd = interp_spline(y2, x2)
    newopd = np.reshape(newopd, (newdim, newdim))

    return newopd


def border_extrapolate_pad(image, mask):
    """Extrapolate phases on an irregular aperture. Sort of an inelegant hack, but useful
    in some contexts, in particular to fill in phase values in segment gaps prior to rescaling
    or interpolation.

    Each invalid pixel adjacent to the aperture is replaced with the mean of the value
    of the valid pixels it is adjacent to, using square connectivity.

    Parameters
    ----------
    image: float ndarray
        image
    mask: int ndarray
        0 for invalid pixels
    """
    masked_im = np.copy(image)
    masked_im[mask == 0] = np.nan

    # 1 pixel wide border region
    # border = scipy.ndimage.morphology.binary_dilation(mask)^mask

    shifted = np.zeros((8, image.shape[0], image.shape[1]))
    i = 0
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            shifted[i] = np.roll(np.roll(masked_im, dy, axis=0), dx, axis=1)
            i += 1

    # Compress into a new 2D image
    with warnings.catch_warnings():
        # this will always produce a warning about taking a mean of an empty slice, which is safe to
        # ignore, so we do so.
        warnings.simplefilter('ignore')
        image_extrapolated = np.nanmean(shifted, axis=0)
    image_extrapolated = np.nan_to_num(image_extrapolated)

    # retain valid values in input image exactly
    image_extrapolated[mask] = image[mask]

    return image_extrapolated


# Helper routines for display customization: ###
# use via poppy's display_annotate feature by assigning these to
# the display_annotate attribute of an OpticalElement class


def annotate_ote_pupil_coords(self, ax, orientation='entrance_pupil'):
    """Draw OTE V frame axes on first optical plane"""
    color = 'yellow'

    xloc = 3
    if orientation == 'entrance_pupil':
        yloc = 3
        v3sign = +1
        v3verticalalignment = 'bottom'
    elif orientation == 'exit_pupil':
        yloc = 2.5
        v3sign = -1
        v3verticalalignment = 'top'
    else:
        raise ValueError(f"Unknown orientation {orientation}. Must be either 'entrance_pupil' or 'exit_pupil'. ")

    ax.arrow(-xloc, -yloc, 0.2, 0, color=color, width=0.005)
    ax.arrow(-xloc, -yloc, 0, 0.2 * v3sign, color=color, width=0.005)
    ax.text(
        -xloc,
        -yloc + 0.4 * v3sign,
        '+V3',
        color=color,
        size='small',
        horizontalalignment='center',
        verticalalignment=v3verticalalignment,
    )
    ax.text(-xloc + 0.4, -yloc, '+V2', color=color, size='small', horizontalalignment='left', verticalalignment='center')


def annotate_sky_pupil_coords(self, ax, show_NE=False, north_angle=45.0):
    """Draw OTE V frame axes projected onto the sky
    Optionally also draw a compass for north and east at some given
    position angle
    """
    color = 'yellow'
    loc = 2.9
    ax.arrow(-loc + 0.5, -loc, -0.2, 0, color=color, width=0.005)
    ax.arrow(-loc + 0.5, -loc, 0, 0.2, color=color, width=0.005)
    ax.text(
        -loc + 0.5,
        -loc + 0.3,
        '+V3 on sky',
        color=color,
        size='small',
        horizontalalignment='center',
        verticalalignment='bottom',
    )
    ax.text(
        -loc + 0.5 + 0.3,
        -loc,
        '+V2 on sky',
        color=color,
        size='small',
        horizontalalignment='left',
        verticalalignment='center',
    )

    if show_NE:
        color2 = 'cyan'
        angle = np.deg2rad(north_angle)  # arbitrary
        dl = 0.3
        dx = np.sin(angle) * dl
        dy = np.cos(angle) * dl
        ax.arrow(-loc + 0.5, -loc, -dx, dy, color=color2, width=0.005)
        ax.arrow(-loc + 0.5, -loc, -dy, -dx, color=color2, width=0.005)
        ax.text(
            -loc + 0.5 - 2.3 * dx,
            -loc + 2.3 * dy,
            'N',
            color=color2,
            size='small',
            horizontalalignment='center',
            verticalalignment='center',
        )
        ax.text(
            -loc + 0.5 - 1.3 * dy,
            -loc - 1.3 * dx,
            'E',
            color=color2,
            size='small',
            horizontalalignment='center',
            verticalalignment='center',
        )


def _run_benchmark(timer, iterations=1):
    """Common benchmarking core. Called from benchmark_imaging and benchmark_coronagraphy"""
    import poppy

    defaults = (poppy.conf.use_fftw, poppy.conf.use_numexpr, poppy.conf.use_cuda, poppy.conf.use_opencl)

    # Time baseline performance in numpy
    print('Timing performance in plain numpy:')
    poppy.conf.use_fftw, poppy.conf.use_numexpr, poppy.conf.use_cuda, poppy.conf.use_opencl = (False, False, False, False)
    time_numpy = timer.timeit(number=iterations) / iterations
    print('  {:.2f} s'.format(time_numpy))

    if poppy.accel_math._FFTW_AVAILABLE:
        print('Timing performance with FFTW:')
        poppy.conf.use_fftw = True
        time_fftw = timer.timeit(number=iterations) / iterations
        print('  {:.2f} s'.format(time_fftw))
    else:
        time_fftw = np.nan

    if poppy.accel_math._NUMEXPR_AVAILABLE:
        print('Timing performance with Numexpr:')
        poppy.conf.use_fftw = False
        poppy.conf.use_numexpr = True
        time_numexpr = timer.timeit(number=iterations) / iterations
        print('  {:.2f} s'.format(time_numexpr))
    else:
        time_numexpr = np.nan

    if poppy.accel_math._CUDA_AVAILABLE:
        print('Timing performance with CUDA + Numexpr:')
        poppy.conf.use_cuda = True
        poppy.conf.use_opencl = False
        time_cuda = timer.timeit(number=iterations) / iterations
        print('  {:.2f} s'.format(time_cuda))
    else:
        time_cuda = np.nan

    if poppy.accel_math._OPENCL_AVAILABLE:
        print('Timing performance with OpenCL + Numexpr:')
        poppy.conf.use_opencl = True
        poppy.conf.use_cuda = False
        time_opencl = timer.timeit(number=iterations) / iterations
        print('  {:.2f} s'.format(time_opencl))
    else:
        time_opencl = np.nan

    poppy.conf.use_fftw, poppy.conf.use_numexpr, poppy.conf.use_cuda, poppy.conf.use_opencl = defaults

    return {'numpy': time_numpy, 'fftw': time_fftw, 'numexpr': time_numexpr, 'cuda': time_cuda, 'opencl': time_opencl}


def benchmark_imaging(iterations=1, nlambda=1, add_distortion=True):
    """Performance benchmark function for standard imaging"""
    import timeit

    timer = timeit.Timer(
        'psf = nc.calc_psf(nlambda=nlambda, add_distortion={})'.format(add_distortion),
        setup="""
import webbpsf
nc = webbpsf.NIRCam()
nc.filter='F360M'
nlambda={nlambda:d}""".format(nlambda=nlambda),
    )
    print('Timing performance of NIRCam F360M with {} wavelengths, {} iterations'.format(nlambda, iterations))

    return _run_benchmark(timer, iterations=iterations)


def benchmark_nircam_coronagraphy(iterations=1, nlambda=1, add_distortion=True):
    """Performance benchmark function for standard imaging"""
    import timeit

    timer = timeit.Timer(
        'psf = nc.calc_psf(nlambda=nlambda, add_distortion={})'.format(add_distortion),
        setup="""
import webbpsf
nc = webbpsf.NIRCam()
nc.filter='F335M'
nc.image_mask='MASK335R'
nc.pupil_mask='MASKRND'
nlambda={nlambda:d}""".format(nlambda=nlambda),
    )
    print('Timing performance of NIRCam MASK335R with {} wavelengths, {} iterations'.format(nlambda, iterations))

    return _run_benchmark(timer, iterations=iterations)


def benchmark_miri_coronagraphy(iterations=1, nlambda=1):
    """Performance benchmark function for standard imaging"""
    import timeit

    timer = timeit.Timer(
        'psf = miri.calc_psf(nlambda=nlambda)',
        setup="""
import webbpsf
miri = webbpsf.MIRI()
miri.filter='F1065C'
miri.image_mask='FQPM1065'
miri.pupil_mask='MASKFQPM'
nlambda={nlambda:d}""".format(nlambda=nlambda),
    )
    print('Timing performance of MIRI F1065C with {} wavelengths, {} iterations'.format(nlambda, iterations))

    return _run_benchmark(timer, iterations=iterations)


def combine_docstrings(cls):
    """Combine the docstrings of a method and earlier implementations of the same method in parent classes"""
    for name, func in cls.__dict__.items():
        # Allow users to see the Poppy calc_psf docstring along with the JWInstrument version
        if name == 'calc_psf':
            jwinstrument_class = cls
            spacetelescope_class = cls.__base__

            ind0 = getattr(jwinstrument_class, 'calc_psf').__doc__.index('add_distortion')  # pull the new parameters
            ind1 = getattr(spacetelescope_class, 'calc_psf').__doc__.index('Returns')  # end of parameters

            func.__doc__ = (
                getattr(spacetelescope_class, 'calc_psf').__doc__[0:ind1]
                + getattr(jwinstrument_class, 'calc_psf').__doc__[ind0:]
                + getattr(spacetelescope_class, 'calc_psf').__doc__[ind1:]
            )

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
        from photutils.psf import GriddedPSFModel
    except ImportError:
        try:
            from photutils import GriddedPSFModel
        except ImportError:
            raise ImportError('This method requires photutils >= 0.6')

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
    if not any('DET_YX' in key for key in header.keys()):
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
    ndd.meta['grid_xypos'] = [
        ((float(ndd.meta[key][0].split(',')[1].split(')')[0])), (float(ndd.meta[key][0].split(',')[0].split('(')[1])))
        for key in ndd.meta.keys()
        if 'DET_YX' in key
    ]  # from (y,x) to (x,y)

    if 'oversampling' not in ndd.meta:
        ndd.meta['oversampling'] = ndd.meta['OVERSAMP'][0]  # pull the value

    # Turn all metadata keys into lowercase
    ndd.meta = {key.lower(): ndd.meta[key] for key in ndd.meta}

    # Create model
    model = GriddedPSFModel(ndd)

    return model


def determine_inst_name_from_v2v3(v2v3):
    """Given a (V2,V3) coordinate, look up which JWST SI FOV that point is within.

    In particular this is used as part of lookup for the OTE field dependence model.
    """
    # Figure out what instrument the field coordinate correspond to
    if (
        (v2v3[0] <= 4.7 * u.arcmin)
        and (v2v3[0] >= -0.9 * u.arcmin)
        and (v2v3[1] <= -10.4 * u.arcmin)
        and (v2v3[1] >= -12.9 * u.arcmin)
    ):
        instrument = 'FGS'
        _log.debug('Field coordinates determined to be in FGS field')
    elif (
        (v2v3[0] <= 2.6 * u.arcmin)
        and (v2v3[0] >= -2.6 * u.arcmin)
        and (v2v3[1] <= -6.2 * u.arcmin)
        and (v2v3[1] >= -9.4 * u.arcmin)
    ):
        instrument = 'NIRCam'
        _log.debug('Field coordinates determined to be in NIRCam field')
    elif (
        (v2v3[0] <= 8.95 * u.arcmin)
        and (v2v3[0] >= 3.7 * u.arcmin)
        and (v2v3[1] <= -4.55 * u.arcmin)
        and (v2v3[1] >= -9.75 * u.arcmin)
    ):
        instrument = 'NIRSpec'
        _log.debug('Field coordinates determined to be in NIRSpec field')
    elif (
        (v2v3[0] <= -6.2 * u.arcmin)
        and (v2v3[0] >= -8.5 * u.arcmin)
        and (v2v3[1] <= -5.2 * u.arcmin)
        and (v2v3[1] >= -7.3 * u.arcmin)
    ):
        instrument = 'MIRI'
        _log.debug('Field coordinates determined to be in MIRI field')
    elif (
        (v2v3[0] <= -3.70 * u.arcmin)
        and (v2v3[0] >= -6.0 * u.arcmin)
        and (v2v3[1] <= -10.5 * u.arcmin)
        and (v2v3[1] >= -12.8 * u.arcmin)
    ):
        instrument = 'NIRISS'
        _log.debug('Field coordinates determined to be in NIRISS field')
    else:
        raise ValueError(f'Given V2V3 coordinates {v2v3} do not fall within an instrument FOV region')

    return instrument


def label_wavelength(nwavelengths, wavelength_slices):
    # Allow up to 10,000 wavelength slices. The number matters because FITS
    # header keys can only have up to 8 characters. Backward-compatible.
    if nwavelengths < 100:
        label = 'WAVELN{:02d}'.format(wavelength_slices)
    elif nwavelengths < 10000:
        label = 'WVLN{:04d}'.format(wavelength_slices)
    else:
        raise ValueError('Maximum number of wavelengths exceeded. ' 'Cannot be more than 10,000.')
    return label


def get_target_phase_map_filename(apername):
    """Get WSS Target Phase Map for the specified aperture
    Note that the sensing maintenance program changed field point from NRC A3 to A1 around Dec 2024.
    """
    path = os.getenv('WEBBPSF_PATH')
    if apername == 'NRCA3_FP1':
        fn = 'wss_target_phase_fp1.fits'
    elif apername == 'NRCA1_FP6':
        fn = 'wss_target_phase_fp6.fits'
    else:
        raise ValueError(f"Target phase map not available for aperture = {apername}")

    was_targ_file = os.path.join(
        get_webbpsf_data_path(), 'NIRCam', 'OPD', fn)

    if not os.path.exists(was_targ_file):
        raise ValueError("File wss_target_phase_{}.fits, \
        not found under {}.".format(apername.split('_')[1].lower(),path))

    return was_targ_file