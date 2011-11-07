#!/usr/bin/env python
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pywcs, pyfits
#import aplpy, atpy
#import RO.DS9
#from IPython.Debugger import Tracer; stop = Tracer()

import poppy


JWPSF_basepath = os.path.dirname(os.path.dirname(os.path.abspath(poppy.__file__))) +os.sep+"data"

def convert_pupil(filename, instname, wfe, outdir=JWPSF_basepath, includes="OPD", summary="OPD"):
    """ This is a copy of convert_one that just converts a single pupil file """
    if not os.path.exists(filename):
        print "File does not exist: %s" % filename

    diam = 6.5

    imstack = pyfits.getdata(filename)
    print filename, imstack.shape

    imstack.shape = (imstack.shape[0]/1024,1024,1024)

    imstack2 = imstack[-2::-1, :, :].copy() # get all the OPDs, and reverse order 
        #(so the one previously at the top becomes first in the new arrangement)
        # this reverse indexing syntax is confusing, but this does indeed return the
        # slices in order [9, 8, ..., 2, 1, 0] as desired.

    imstack3 = imstack2.transpose((0,2,1))[:,::-1,:] #transpose and flip Y to get +y = +V3

    plt.imshow(imstack3[0,:,:])
    pupil = imstack[-1, :, :].transpose()[::-1,:]





    hdu = pyfits.PrimaryHDU(pupil)
    hdu.header.update("EXTNAME", "PUPIL")



    # check for blank padding rows. 
    if (hdu.data[0:4,:].sum() ==0 ) and (hdu.data[:,0:4].sum() ==0) :
        true_pupil_diam = 1016
    else:
        true_pupil_diam = 1024

    hdu.header.update("SOURCE", os.path.basename(filename), "File this was created from")
    hdu.header.update("ORIENT", "Y axis is +V3, X axis is + or - V2 (arbitrary)")

    pupilscale = diam*1.0/true_pupil_diam
    file_diam = hdu.data.shape[1] * pupilscale

    hdu.header.update('PUPLDIAM', file_diam, 'Pupil *file* diameter in meters')
    hdu.header.update("DIAM", 6.5, "True Pupil diameter in meters (ignoring padding)") # I am not sure this is 100% correct
    hdu.header.update('PUPLSCAL', pupilscale, 'Pupil pixel scale in meters/pixel')


    outname = "%s/pupil_RevV.fits" % (outdir)

    hdu.writeto(outname, clobber=True)
    print "==>> %s" % outname





def convert_one(filename, instname, wfe, outdir=JWPSF_basepath, includes="OPD", summary="OPD"):
    if not os.path.exists(filename):
        print "File does not exist: %s" % filename


    diam = 6.5

    imstack = pyfits.getdata(filename)
    print filename, imstack.shape

    imstack.shape = (imstack.shape[0]/1024,1024,1024)

    #np.transpose(imstack, axes=[0,2,1]) # swap X and Y to align Y = +V3

    imstack2 = imstack[-2::-1, :, :].copy() # get all the OPDs, and reverse order 
        #(so the one previously at the top becomes first in the new arrangement)
        # this reverse indexing syntax is confusing, but this does indeed return the
        # slices in order [9, 8, ..., 2, 1, 0] as desired.

    imstack3 = imstack2.transpose((0,2,1))[:,::-1,:] #transpose and flip Y to get +y = +V3

    plt.imshow(imstack3[0,:,:])
    pupil = imstack[-1, :, :].transpose()[::-1,:]





    hdu = pyfits.PrimaryHDU(imstack3)
    pupilext = pyfits.ImageHDU(pupil)
    pupilext.header.update("EXTNAME", "PUPIL")



    # check for blank padding rows. 
    if (pupilext.data[0:4,:].sum() ==0 ) and (pupilext.data[:,0:4].sum() ==0) :
        true_pupil_diam = 1016
    else:
        true_pupil_diam = 1024

    hdu.header.update("SOURCE", os.path.basename(filename), "File this was created from")
    hdu.header.update("ORIENT", "Y axis is +V3, X axis is + or - V2 (arbitrary)")
    hdu.header.update("INSTRUME", instname, "Which JWST instrument is this for?")
    hdu.header.update("WFE_RMS", float(wfe), "RMS WFE [nm]")
    hdu.header.update("WFE_INCL", includes)
    hdu.header.update("SUMMARY", summary)
    hdu.header.update("EXTNAME", "OPDs")

    pupilscale = diam*1.0/true_pupil_diam
    file_diam = hdu.data.shape[1] * pupilscale

    hdu.header.update('PUPLDIAM', file_diam, 'Pupil *file* diameter in meters')
    hdu.header.update("DIAM", 6.5, "True Pupil diameter in meters (ignoring padding)") # I am not sure this is 100% correct
    hdu.header.update('PUPLSCAL', pupilscale, 'Pupil pixel scale in meters/pixel')


    outname = "%s/OPD_RevV_%s_%d.fits" % (outdir, instname, wfe)

    HL = pyfits.HDUList([hdu, pupilext])

    HL.writeto(outname, clobber=True)
    print "==>> %s" % outname






if __name__ == "__main__":
        revTdir = '/itar/jwst/tel/share/WFS&C/Simulated OPDs/OPD rev V'

        instr = [ ('MIRI', 'miri', [204,220,421]), 
                 ('FGS', 'fgs', [150, 163, 186]), 
                 ('NIRCam', 'nircam', [115, 132, 150]),
                 ('NIRCam', 'nircam', [123, 136, 155]),
                 ('NIRSpec', 'nirspec',[125, 145, 238]), 
                 ('TFI', 'tf', [144, 162, 180])]


        i,fn,wfes = instr[0]
        convert_pupil(revTdir+"/ipam_opd_%s%d_wo_im.fits"%(fn, wfes[0]), fn, wfes[0], outdir=JWPSF_basepath) 



        if 0:
            for i, fn, wfes in instr:
                #fn = instr[i][0]
                #wfes = instr[i][1]
                convert_one(revTdir+"/ipam_opd_%s%d_wo_im.fits"%(fn, wfes[0]), fn, wfes[0], outdir=JWPSF_basepath+os.sep+i+"/OPD", 
                        includes='OPD for OTE+ISIM with all reserves and stability; NO image motion',
                        summary = '%d nm RMS for OTE+ISIM'% (wfes[0]) )

                for wf in wfes[1:]:
                    if wf == wfes[-1]:
                        includes = 'Observatory OPD for OTE+ISIM+SI with all reserves and stability, and image motion as defocus'
                        summary = '%d nm RMS for OTE+ISIM + %s, plus image motion as defocus' % (wf, i)
                    else:
                        includes = 'Observatory OPD for OTE+ISIM with all reserves and stability, and image motion as defocus'
                        summary = '%d nm RMS for OTE+ISIM, plus image motion as defocus' % wf
                    convert_one(revTdir+"/ipam_opd_%s%d.fits"%(fn, wf), fn, wf, outdir=JWPSF_basepath+os.sep+i+"/OPD", includes=includes, summary=summary)



