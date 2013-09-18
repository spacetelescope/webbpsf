#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# --based on setup.py from astropy--

import glob
import os
import sys

import setuptools_bootstrap
from setuptools import setup, find_packages

import astropy
from astropy.setup_helpers import (register_commands, adjust_compiler,
                                   filter_packages, update_package_files,
                                   get_debug_option)
from astropy.version_helpers import get_git_devstr, generate_version_py

PACKAGENAME = 'webbpsf'
DESCRIPTION = 'Create simulated point spread functions for the James Webb Space Telescope'
LONG_DESCRIPTION = """
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
AUTHOR = 'Marshall Perrin'
AUTHOR_EMAIL = 'mperrin@stsci.edu'
LICENSE = 'BSD'
URL = 'http://www.stsci.edu/~mperrin/software/webbpsf'



# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.3.dev'

# Indicates if this version is a release version
RELEASE = 'dev' not in VERSION

if not RELEASE:
    VERSION += get_git_devstr(False)

#DOWNLOAD_BASE_URL = 'http://pypi.python.org/packages/source/w/webbpsf'

# Populate the dict of setup command overrides; this should be done before
# invoking any other functionality from distutils since it can potentially
# modify distutils' behavior.
cmdclassd = register_commands(PACKAGENAME, VERSION, RELEASE)

# Adjust the compiler in case the default on this platform is to use a
# broken one.
adjust_compiler(PACKAGENAME)

# Freeze build information in version.py
generate_version_py(PACKAGENAME, VERSION, RELEASE, get_debug_option())

## 
## 
## setupargs = {
##     'name'          :       PACKAGENAME,
##     'app'           :       'WebbPSF',
##     'version'       :      	VERSION,
##     'description'   :       '',
##     'fullname'      :       'WebbPSF',
##     'author'        :     	"Marshall Perrin",
##     'author_email'  :      	"mperrin@stsci.edu",
##     'url'           :  		"http://www.stsci.edu/~mperrin/software/webbpsf",
##     'download_url'           :  		"http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-0.0.0.tar.gz",  # will be replaced below
##     'platforms'     :      	["Linux","Mac OS X", "Win"],
## #    'requires'      :       ['pyfits','numpy', 'matplotlib', 'scipy', 'asciitable', 'poppy'],
## #    'packages'      :       ['webbpsf'],
## #    'entry_points'  :       {'gui_scripts': ['webbpsfgui = webbpsf.gui',]}, # should create exe file on Windows?
##     'classifiers'   :   [
##         "Programming Language :: Python",
##         "License :: OSI Approved :: BSD License",
##         "Operating System :: OS Independent",
##         "Intended Audience :: Science/Research",
##         "Topic :: Scientific/Engineering :: Astronomy",
##         'Topic :: Scientific/Engineering :: Physics',
##         "Development Status :: 4 - Beta"
##         ],
##     'long_description': """
## 
## WebbPSF: Simulated Point Spread Functions for the James Webb Space Telescope
## -------------------------------------------------------------------------------
## 
## WebbPSf produces simulated PSFs for the James Webb Space Telescope, NASA's next flagship
## infrared space telescope. WebbPSF can simulate images for any of the four science instruments plus the
## fine guidance sensor, including both direct imaging and coronagraphic modes. 
## 
## Developed by Marshall Perrin at STScI, 2010-2012. 
## 
## Documentation can be found online at http://www.stsci.edu/jwst/software/webbpsf/
## 
## WebbPSF requires a large amount of input data for its simulations, including optical path difference (OPD) maps,
## filter transmission curves, and coronagraph Lyot mask shapes. These data files are not included in this
## source distribution available from PYPI. Please see the main WebbPSF web page, linked above, to download
## the required data tar file.
## 
## """
##     }
## 
## read in the version number from the code itself, following 
## http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
#
## don't actually import the _version.py, for the reasons described on that web page. 
#import re
#VERSIONFILE="webbpsf/_version.py"
#verstrline = open(VERSIONFILE, "rt").read()
#VSRE = r"__version__ = ['\"]([^'\"]*)['\"]"
#mo = re.search(VSRE, verstrline, re.M)
#if mo:
#    setupargs['version'] = mo.group(1)
#    setupargs['download_url'] = "http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-"+mo.group(1)+".tar.gz"
#else:
#    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))
#


## # Use the find_packages tool to locate all packages and modules
packagenames = filter_packages(find_packages())

# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
           if os.path.basename(fname) != 'README.rst']

# Additional C extensions that are not Cython-based should be added here.
extensions = []

# A dictionary to keep track of all package data to install
package_data = {PACKAGENAME: ['data/*']}

# A dictionary to keep track of extra packagedir mappings
package_dirs = {}

# Update extensions, package_data, packagenames and package_dirs from
# any sub-packages that define their own extension modules and package
# data.  See the docstring for setup_helpers.update_package_files for
# more details.
update_package_files(PACKAGENAME, extensions, package_data, packagenames,
                     package_dirs)


print sys.argv

if 'py2app' in sys.argv:
    # standalone Mac .app bundle builds

    PY2APP_OPTIONS = {'argv_emulation': True, 
                    'iconfile': 'webbpsf_icon.icns'}

    setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      app = 'standalone_app.py',
      options={'py2app':PY2APP_OPTIONS},
      setup_requires=['py2app'],
      packages=packagenames,
      package_data=package_data,
      package_dir=package_dirs,
      ext_modules=extensions,
      scripts=scripts,
      requires=['astropy'],
      install_requires=['astropy'],
      provides=[PACKAGENAME],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url=URL,
      long_description=LONG_DESCRIPTION,
      cmdclass=cmdclassd,
      zip_safe=False,
      use_2to3=False
      )
    os.rename('dist/standalone_app.py', 'dist/WebbPSF.py')

else:

    setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      packages=packagenames,
      package_data=package_data,
      package_dir=package_dirs,
      ext_modules=extensions,
      scripts=scripts,
      requires=['astropy'],
      install_requires=['astropy'],
      provides=[PACKAGENAME],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url=URL,
      long_description=LONG_DESCRIPTION,
      cmdclass=cmdclassd,
      zip_safe=False,
      use_2to3=False
      )
