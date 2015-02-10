# Make minimal data files
# This is used to make a stripped-down version of the data files for use on Travis CI

import os, sys
import astropy.io.fits as fits
import subprocess
import glob

insts = ['FGS', 'NIRCam', 'NIRSpec','NIRISS','MIRI']


subprocess.call("mkdir ~/tmp/minimal-webbpsf-data", shell=True)
os.chdir(os.path.expanduser("~/tmp/minimal-webbpsf-data"))
subprocess.call("tar xvzf ~/web/software/webbpsf/webbpsf-data-0.3.0.tar.gz", shell=True)
os.chdir("webbpsf-data")


for instr in insts:
    files = glob.glob(os.path.join(instr, "OPD", "*.fits"))
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

os.remove('pupil_RevT.fits')
os.remove('tricontagon.fits')


subprocess.call('tar cvzf minimal-webbpsf-data.tar.gz', shell=True)
print "===>  minimal-webbpsf-data.tar.gz "


