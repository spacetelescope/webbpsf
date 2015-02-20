# Make minimal data files
# This is used to make a stripped-down version of the data files for use on Travis CI

import os, sys
import astropy.io.fits as fits
import subprocess
import glob

insts = ['FGS', 'NIRCam', 'NIRSpec','NIRISS','MIRI']

WORKING_DIR =  os.path.expanduser("~/tmp/minimal-webbpsf-data")
subprocess.call("mkdir "+WORKING_DIR, shell=True)


print "#### Expanding full tar file into temp directory ####"
os.chdir(WORKING_DIR)
subprocess.call("tar xvzf ~/web/software/webbpsf/webbpsf-data-0.3.1.tar.gz", shell=True)

print "#### Trimming to only one OPD file per instrument ####"
for instr in insts:
    files = glob.glob(os.path.join(WORKING_DIR, 'webbpsf-data', instr, "OPD", "*.fits"))
    files.sort()
    print instr, files
    
    # just save the lowest alphabetically of each of them
    for file_to_delete in files[1:]:
        print "Deleting "+file_to_delete
        os.remove(file_to_delete)

    print "Trimming to only 1 datacube slice: "+files[0]

    f0 = fits.open(files[0], mode='update')
    f0[0].data = f0[0].data[0]
    f0.flush()
    f0.close()

print "#### Removing extra optional pupil files ####"
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','pupil_RevT.fits'))
os.remove(os.path.join(WORKING_DIR, 'webbpsf-data','tricontagon.fits'))

print "#### Creating tar file ####"
os.chdir(WORKING_DIR)
subprocess.call('tar cvzf minimal-webbpsf-data.tar.gz webbpsf-data', shell=True)
print "===>  ~/tmp/minimal-webbpsf-data.tar.gz "


