# See http://packages.python.org/distribute/setuptools.html
from setuptools import setup, find_packages

setupargs = {
    'name'          :       'webbpsf',
    'app'           :       'WebbPSF',
    'version'       :      	"0.0.0",  # will be replaced below
    'description'   :       'Create simulated point spread functions for the James Webb Space Telescope',
    'fullname'      :       'WebbPSF',
    'author'        :     	"Marshall Perrin",
    'author_email'  :      	"mperrin@stsci.edu",
    'url'           :  		"http://www.stsci.edu/~mperrin/software/webbpsf",
    'download_url'           :  		"http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-0.0.0.tar.gz",  # will be replaced below
    'platforms'     :      	["Linux","Mac OS X", "Win"],
    'requires'      :       ['pyfits','numpy', 'matplotlib', 'scipy', 'asciitable', 'poppy'],
    'packages'      :       ['webbpsf'],
    'entry_points'  :       {'gui_scripts': ['webbpsfgui = webbpsf.gui',]}, # should create exe file on Windows?
    'classifiers'   :   [
        "Programming Language :: Python",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Development Status :: 4 - Beta"
        ],
    'long_description': """

WebbPSF: Simulated Point Spread Functions for the James Webb Space Telescope
-------------------------------------------------------------------------------

WebbPSf produces simulated PSFs for the James Webb Space Telescope, NASA's next flagship
infrared space telescope. WebbPSF can simulate images for any of the four science instruments plus the
fine guidance sensor, including both direct imaging and coronagraphic modes. 

Developed by Marshall Perrin at STScI, 2010-2012. 

Documentation can be found online at http://www.stsci.edu/jwst/software/webbpsf/

WebbPSF requires a large amount of input data for its simulations, including optical path difference (OPD) maps,
filter transmission curves, and coronagraph Lyot mask shapes. These data files are not included in this
source distribution available from PYPI. Please see the main WebbPSF web page, linked above, to download
the required data tar file.

"""
    }

# read in the version number from the code itself, following 
# http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package

# don't actually import the _version.py, for the reasons described on that web page. 
import re
VERSIONFILE="webbpsf/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    setupargs['version'] = mo.group(1)
    setupargs['download_url'] = "http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-"+mo.group(1)+".tar.gz"
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))






PY2APP_OPTIONS = {'argv_emulation': True, 
                'iconfile': 'webbpsf_icon.icns'}




# Now actually call setup

setup(  options={'py2app':PY2APP_OPTIONS},
        setup_requires=['py2app'],
        **setupargs)
