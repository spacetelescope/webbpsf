import pyfits
import numpy as N
import glob

def convert_pupil(filename, diam=6.5):
    """ Various manipulations to bring the OPD files or pupil images produced from
    ITM or other tools into better compliance with the FITS standard, and with expected
    usage format for the JWST PSF modeling tools. This includes, for instance, adding
    the expected PUPLDIAM and PUPLSCAL keywords
    """
    f = pyfits.open(filename, mode='update')

    print "--- updating file: "+filename

    if 'PUPLSCAL' in f[0].header.keys():
        print "That file "+filename+" appears already to have a scale specified. No change needed."
        return

    # check for blank padding rows. 
    if (f[0].data[0:4,:].sum() ==0 ) and (f[0].data[:,0:4].sum() ==0) :
        newdata = f[0].data.copy()
        newdata = newdata[4:-4, 4:-4]
        f[0].header.add_history('Unnecessary blank padding rows trimmed off by trim_pupil.py')
        f[0].data = newdata


    # check for invalid/illegal chars in keyword names
    # for some reason the Ball software produces these
    # Replace the bogus ZCSG(9) with ZCSG_9 etc.
    h=f[0].header
    if '(' in " ".join(h.keys()):
        badkeys  = [k for k in h.keys() if '(' in k]
        for k in badkeys:
            f[0].header.rename_key(k, k.replace('(',"_").replace(')','' ) )
            print "\tReplaced keyword %s to %s" % (k, k.replace('(',"_").replace(')',''))

    # and make sure all are upper case
    for k in f[0].header.keys():

        if k != k.upper():
            print  "\tReplaced keyword %s to %s" % (k, k.upper() )
            f[0].header.rename_key(k, k.upper())

    f[0].header.update('PUPLDIAM', diam, 'Pupil diameter in meters')
    f[0].header.update('PUPLSCAL', diam*1.0 / f[0].data.shape[0], 'Pupil pixel scale in meters/pixel')
    f.close()

    print("Updated FITS file "+filename)


def convert_all(dir):
    fits = glob.glob(dir+"/*.fits")
    for f in fits:
        convert_pupil(f)
