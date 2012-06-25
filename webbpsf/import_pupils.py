import pyfits
import glob

def convert_pupil(filename, diam=6.5):
    """ Various manipulations to bring the OPD files or pupil images produced from
    ITM or other tools into better compliance with the FITS standard, and with expected
    usage format for the JWST PSF modeling tools. This includes, for instance, adding
    the expected PUPLDIAM and PUPLSCAL keywords
    """

    print "Converting pupil file %s" % filename

    f = pyfits.open(filename , mode='update')

    print "--- updating file: "+filename

    if 'PUPLSCAL' in f[0].header.keys():
        print "That file "+filename+" appears already to have a scale specified. No change needed."
        return

    # check for blank padding rows. 
    # CHANGE - no longer delete these - leave them as 1024x1024 for better FFTs
    if (f[0].data[0:4,:].sum() ==0 ) and (f[0].data[:,0:4].sum() ==0) :
        #newdata = f[0].data.copy()
        #newdata = newdata[4:-4, 4:-4]
        #f[0].header.add_history('Unnecessary blank padding rows trimmed off by trim_pupil.py')
        #f[0].data = newdata
        true_pupil_diam = 1016
    else:
        true_pupil_diam = 1024


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

    pupilscale = diam*1.0/true_pupil_diam
    file_diam = f[0].data.shape[0] * pupilscale

    #f[0].header.update('PUPLDIAM', file_diam, 'Pupil *file* diameter in meters')
    f[0].header.update('DIAM', diam, 'True Pupil diameter in meters (ignoring padding)')
    f[0].header.update('PUPLSCAL', pupilscale, 'Pupil pixel scale in meters/pixel')
    f.close(output_verify='ignore')
    #f.writeto(filename,output_verify='ignore')

    print("Updated FITS file "+filename)


def convert_all(dir):
    fits = glob.glob(dir+"/*.fits")
    for f in fits:
        convert_pupil(f)


if __name__ == "__main__":
    convert_all("/Users/mperrin/software/newJWPSF/data/MIRI/OPD")
    convert_all("/Users/mperrin/software/newJWPSF/data/MIRI/coronagraph")
    convert_all("/Users/mperrin/software/newJWPSF/data/NIRCam/OPD")
    convert_all("/Users/mperrin/software/newJWPSF/data/NIRSpec/OPD")
    convert_all("/Users/mperrin/software/newJWPSF/data/TFI/OPD")
    convert_all("/Users/mperrin/software/newJWPSF/data")
