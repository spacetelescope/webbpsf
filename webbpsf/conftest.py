# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
except ImportError:  # In case this plugin is not installed
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}

# This really depends on how you set up your package version,
# modify as needed.
try:
    from webbpsf import __version__ as version
except ImportError:
    version = ''

def pytest_configure():
    PYTEST_HEADER_MODULES.pop('Pandas', None)
    PYTEST_HEADER_MODULES.pop('h5py', None)
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['Photutils'] = 'photutils'
    PYTEST_HEADER_MODULES['Poppy'] = 'poppy'
    TESTED_VERSIONS['Webbpsf'] = version
