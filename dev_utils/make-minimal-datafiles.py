#!/usr/bin/env python
# Make minimal data files
# This is used to make a stripped-down version of the data files for use on GitHub Actions

import os, sys
import astropy.io.fits as fits
import subprocess
import glob

try:
    inputfile = sys.argv[1]
except IndexError:
    print("""ERROR - no input data file provided.\n\nUsage: make-minimal-datafiles.py  <path_to_full_data.tar.gz> <version\n""")
    sys.exit(1)

try:
    version = sys.argv[2]
except IndexError:
    print("""ERROR - no version number provided.\n\nUsage: make-minimal-datafiles.py  <path_to_full_data.tar.gz> <version>\n""")
    sys.exit(1)

insts = ['FGS', 'NIRCam', 'NIRSpec','NIRISS','MIRI']

WORKING_DIR =  os.path.expanduser(f"~/tmp/minimal-webbpsf-data-{version}")
subprocess.call("mkdir -p "+WORKING_DIR, shell=True)


print("#### Expanding full tar file into temp directory ####")
os.chdir(WORKING_DIR)
subprocess.call("tar xvzf "+inputfile, shell=True)

# Remove this for >= 1.1.0:
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','JWST_OTE_OPD_RevAA_prelaunch_predicted.fits.gz'))


print("#### Removing extra optional pupil files ####")
# keep just the 1024 and 2048 ones for tests; don't need the rest
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','jwst_pupil_RevW_npix4096.fits.gz'))
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','jwst_pupil_RevW_npix16384.fits.gz'))
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','JWpupil_segments_RevW_npix4096.fits.gz'))
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','tricontagon.fits.gz'))

print("#### Creating tar file ####")
os.chdir(WORKING_DIR)
subprocess.call(f'tar cvzf minimal-webbpsf-data-{version}.tar.gz webbpsf-data', shell=True)
print(f"===>  {WORKING_DIR}/minimal-webbpsf-data-{version}.tar.gz ")


