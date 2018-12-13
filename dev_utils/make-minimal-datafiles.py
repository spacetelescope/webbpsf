#!/usr/bin/env python
# Make minimal data files
# This is used to make a stripped-down version of the data files for use on Travis CI

import os, sys
import astropy.io.fits as fits
import subprocess
import glob

try:
    inputfile = sys.argv[1]
except IndexError:
    print("""ERROR - no input data file provided.\n\nUsage: make-minimal-datafiles.py  <path_to_full_data.tar.gz>\n""")
    sys.exit(1)


insts = ['FGS', 'NIRCam', 'NIRSpec','NIRISS','MIRI']

WORKING_DIR =  os.path.expanduser("~/tmp/minimal-webbpsf-data")
subprocess.call("mkdir -p "+WORKING_DIR, shell=True)


print("#### Expanding full tar file into temp directory ####")
os.chdir(WORKING_DIR)
subprocess.call("tar xvzf "+inputfile, shell=True)

print("#### Trimming to only one OPD file per instrument ####")
for instr in insts:
    files = glob.glob(os.path.join(WORKING_DIR, 'webbpsf-data', instr, "OPD", "*.fits.gz"))
    files.sort()
    print(instr, files)

    # just save the lowest alphabetically of each of them
    for file_to_delete in files[1:]:
        print("Deleting "+file_to_delete)
        os.remove(file_to_delete)

    print("Trimming to only 1 datacube slice: "+files[0])

    f0 = fits.open(files[0])
    f0[0].data = f0[0].data[0]
    f0.writeto(files[0], overwrite=True)
    f0.close()

print("#### Removing extra optional pupil files ####")
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','jwst_pupil_revW_npix2048.fits.gz'))
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','jwst_pupil_revW_npix16384.fits.gz'))
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','tricontagon.fits.gz'))

print("#### Creating tar file ####")
os.chdir(WORKING_DIR)
subprocess.call('tar cvzf minimal-webbpsf-data.tar.gz webbpsf-data', shell=True)
print("===>  {0}/minimal-webbpsf-data.tar.gz ".format(WORKING_DIR))


