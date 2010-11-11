#!/usr/bin/env  python

import utils 
from arrayinfo import  *
from numarray import  *
from OpticalTrainMW import *
from MakeDisk import *

# /Users/anand/Desktop/People/Remi/APLCtol/Gemini_APLC_4.7.fits

def do_ot(dict):

    ot = OpticalTrain(name=dict["trainname"],  \
                       pupil=dict["pup"],\
                       phase=zeros(dict["pup"].shape, type=Float64), \
                       oversample=dict["OVERSAMPLE"], verbose=0x42)
    # CONSTRUCT THE LYOT STOP AND IMAGE STOP (transparent or occulted)

    rn = dict["rn"]

    ot.setsavebox_pixels(saveboxpixels=dict["saveN"]) 
    ot.setpsavebox(psavebox=dict["saveN"])

    ot.saveIntensity(planename='entrance pupil', fn=rn+'Xap.fits', 
                     clobber='y')

    # FIRST IMAGE
    ot.toConjugate(operation='first image')
    ot.saveIntensity(planename='first image', fn=rn+'Xnc.fits', 
                     clobber='y')
    inpower = ot.totalPower('first image') 
    print "incoming on-axis power = %g"  %  inpower

    ot.deletePlane(planename='entrance pupil')
    

    # DIRECT WITHOUT FPM
    ot.amplitudeMask(mask=dict["imagestop"], operation='image stop')
    print  ot.planeinfo(planename='image stop')
    ot.deletePlane(planename='first image')
    ot.saveTransmissionMask(planename='image stop', fn=rn+'Xstop.fits', 
                     clobber='y')
    ot.deletePlane(planename='image stop')
    ot.saveIntensity(planename='after image stop', fn=rn+'Xims.fits', 
                     clobber='y')


    # To LYOT PUPIl
    ot.toConjugate(operation='Lyot pupil')
    ot.deletePlane(planename='after image stop')
    ot.saveIntensity(planename='Lyot pupil', fn=rn+'XLy.fits', 
                     clobber='y')

    # LYOT STOP 
    ot.amplitudeMask(mask=dict["lyotstop"], operation='lyot stop')
    ot.saveIntensity(planename='Lyot pupil', fn=rn+'Xly.fits', 
                     clobber='y')
    ot.saveIntensity(planename='after lyot stop', fn=rn+'Xply.fits', 
                     clobber='y')
    ot.deletePlane(planename='Lyot pupil')


    # FINAL BUT NO STOP IMAGE
    ot.toConjugate(operation='final image')
    ycpower = ot.totalPower('final image') 
    print "on-axis direct power = %g (not normalized)"  %  ycpower

    ot.displaytrain()

    ot.saveIntensity(planename='final image', fn=rn+'Xyc.fits', 
                     clobber='y')

    lyotthru = add.reduce(dict["lyotstop"].flat) / add.reduce(dict["clearaperture"].flat)
    print "lyotthru = %g"  %  lyotthru

    return inpower, ycpower, lyotthru




def doit(coro="generate_first"):

    """ Remi's fits file is RFN x RFN with his floating point oversampling  
        NFFTW is the fft array size we must use for eg fftw'ing the pupil.
        We fake out OpticalTrain by using an oversample of unity... do out own
        padding from RFN to NFFTW
    """


    remisfile = "/Users/anand/Desktop/People/Remi/APLCtol/Gemini_APLC_4.7.fits"
    print remisfile 
    fitscube = SimpleFitsRead(fn=remisfile, report=0x42)#, verbose=0x42)
    hdrList = fitscube.hdr.items()
    for h in hdrList:
        k,v = h
        print type(v), k, v
        
    # Remi to Anand Conversion: 
    roundoff = 0.001
    tdiapix = fitscube.hdr['tdiapix']
    NFFTW = fitscube.hdr['ftarray']
    oversample = int(roundoff + NFFTW/tdiapix)
    print "RS's oversample =", oversample, "\nRS's NFFTW =", NFFTW
    OVERSAMPLE = 1 # to fake out OT...

    # NAXIS3  = 5 / Pupil, Apodizer, FPM, Lyot Stop, PSF           
    remiap = fitscube.data[0,:,:].copy() 
    print "RS's saved file size =", remiap.shape
    RFN = remiap.shape[0]
    remiT = fitscube.data[1,:,:]
    remiLyot = fitscube.data[3,:,:]
    remiPSF = fitscube.data[4,:,:]
    remiQn = fitscube.data[5,:,:]

    remi_r = tdiapix/2


    # handoff ap vars to my traditional var names
    fftwshape = (NFFTW, NFFTW)
    s = fftwshape[0]
    ctr = (NFFTW/2 + 0.5, NFFTW/2 + 0.5)
    c = ctr[0]

    # old... remiFPM = fitscube.data[2,:,:]
    #kmbk_s = fitscube.hdr['FPMSIZE']
    skmbk_s = 4.7 # header says FPMSIZE 4.87067905815
    print "skmbk_s", skmbk_s, "  in pixels this is ", skmbk_s*oversample
    remiFPM = MakeDisk( center=(RFN/2, RFN/2), \
                        s=remiap.shape,  \
                        inside = 0, outside=1, \
                        radius=oversample * skmbk_s * 0.5,  \
                        t=Float64,
                        grey=1).data
    showfpmctr = remiFPM.copy()
    showfpmctr[256,256] = 1.0
    SimpleFitsWrite(fn="fpm.fits", data=showfpmctr, clobber='y')


    print "Remi's file arrays (%dx%d) are snippets of the full FFTW arrays..." % remiap.shape
    print "s/2-RFN/2 = %d - %d = %d"  % \
          (s/2,RFN/2, s/2-RFN/2)
    print "s/2+RFN/2 = %d + %d = %d"  % \
          (s/2,RFN/2, s/2+RFN/2)

    clearaperture = zeros(fftwshape, type=Float64)
    clearaperture[s/2-RFN/2:s/2+RFN/2, s/2-RFN/2:s/2+RFN/2] = remiap.copy()

    T = zeros(fftwshape, type=Float64)  # Transmission acts on Field Strength
    T[s/2-RFN/2:s/2+RFN/2, s/2-RFN/2:s/2+RFN/2] = remiT.copy()

    FPM = ones(fftwshape, type=Float64)
    FPM[s/2-RFN/2:s/2+RFN/2, s/2-RFN/2:s/2+RFN/2] = remiFPM.copy()

    Lyot = zeros(fftwshape, type=Float64)
    Lyot[s/2-RFN/2:s/2+RFN/2, s/2-RFN/2:s/2+RFN/2] = remiLyot.copy()

    PSF = zeros(fftwshape, type=Float64)
    PSF[s/2-RFN/2:s/2+RFN/2, s/2-RFN/2:s/2+RFN/2] = remiPSF.copy()

    apodaperture = T * clearaperture
    T2 = apodaperture * apodaperture
    print "Apodizer max %.2f" % T2.max(),
    print "Apodizer thuput %.4f" % (T2.sum() / clearaperture.sum())

    pup = apodaperture


    phase = zeros(fftwshape, type=Float64)


    trn = "/Users/anand/data/EXAOC/aplc_gpi_pupilwander/"
    rn_NOFPM = trn + "unocculted/"
    rn_FPM = trn + "occulted/"

    if coro == "generate_first":
        print " Coronagraph simulation needed - new subdirs will be created"

        utils.newdir(trn)
        utils.newdir(rn_NOFPM)
        utils.newdir(rn_FPM)


        dict = {}
        dict["rn"] = rn_NOFPM
        dict["trainname"] = "APLC GPI full corotrain sans FPM"
        dict["rn"] = rn_NOFPM
        dict["saveN"] = RFN
        dict["pup"] = pup
        dict["OVERSAMPLE"] = OVERSAMPLE
        dict["imagestop"] = ones(fftwshape, type=Float64)
        dict["lyotstop"] = Lyot
        dict["clearaperture"] = clearaperture

        inpower, ncpower, lyotthru = do_ot(dict)
        print "do_ot returns   inpower", inpower, "  ncpower", ncpower, "  lyotthru", lyotthru


        print "----------------------------------"
        print "----------------------------------"
        print "----------------------------------"

        dict = {}
        dict["rn"] = rn_FPM
        dict["trainname"] = "APLC GPI full corotrain avec FPM"
        dict["saveN"] = RFN
        dict["pup"] = pup
        dict["OVERSAMPLE"] = OVERSAMPLE
        dict["imagestop"] = FPM
        dict["lyotstop"] = Lyot
        dict["clearaperture"] = clearaperture

        inpower, ycpower, lyotthru = do_ot(dict)
        print "do_ot returns   inpower", inpower, "  ycpower", ycpower, "  lyotthru", lyotthru


    # No coronagraph simulation needed - all is on disk

    # for everything with an image stop, normalize to peak of
    # thru-Lyot unstopped PSF, in Macintosh & Krist style

    
    # yes, it IS yc... FPM was transparent
    nc = SimpleFitsRead(fn=rn_NOFPM+'Xyc.fits').data  
    yc = SimpleFitsRead(fn=rn_FPM+'Xyc.fits').data


    ncmax = nc.max()
    print "My files on disk...", 
    print yc.shape
    print "Max of Xyc = %.7g" % yc.max()
    print "Sum of Xyc = %.7g" % yc.sum()
    ycN = yc / nc.max()
    print "Max of ycN = %.7g afte normalization by NO_FPM" % ycN.max()
    print "Sum of ycN = %.7g afte normalization by NO_FPM" % ycN.sum()
    SimpleFitsWrite(fn=trn+"ycN.fits", data=ycN, clobber="y")

    Nas = nc.shape[0]
    Nrs = remiap.shape[0]
    print "Nas = ", nc.shape[0], " Nrs = ", remiap.shape[0]
    
    # normalize to Remi's PSF max 
    myPSF = yc
    print myPSF.shape
    print remiPSF.shape

    NormPeak = True
    anRS_PSF = remiPSF
    comment_remiPSF = """ 
    NormPeak = False =>equate sums
    Tot of AS PSF = 5.960419
    Tot of RS PSF = 0.01723918
    Multiply AS PSF by a factor 0.00289227598905
    RemiPSF.max = 0.0003182  RemiPSF.sum = 0.01724
    myPSFN.max = 0.0009038  myPSFN.sum = 0.01724
    diff(RemiPSF - AnandPSF) = -9.96e-18
    diff(RemiPSF - AnandPSF)stdev = 2.81e-06
    
    NormPeak = True =>equate peaks
    Multiply AS PSF by a factor 0.00101837831744
    RemiPSF.max = 0.0003182  RemiPSF.sum = 0.01724
    myPSFN.max = 0.0003182  myPSFN.sum = 0.00607
    diff(RemiPSF - AnandPSF) = 0.0112
    diff(RemiPSF - AnandPSF)stdev = 2.38e-06

    """
    anRS_PSF = remiQn  

    comment_remiQn = """ 
    NormPeak = False =>equate sums
    Tot of AS PSF = 5.960419
    Tot of RS PSF = 0.08555262
    Multiply AS PSF by a factor 0.0143534585402
    RemiPSF.max = 0.0004349  RemiPSF.sum = 0.08555
    myPSFN.max = 0.004485  myPSFN.sum = 0.08555
    diff(RemiPSF - AnandPSF) = -5.81e-17
    diff(RemiPSF - AnandPSF)stdev = 1.65e-05

    NormPeak = True =>equate peaks
    Max of AS PSF = 0.3124749
    Max of RS PSF = 0.0004348653
    Multiply AS PSF by a factor 0.00139168047996
    RemiPSF.max = 0.0004349  RemiPSF.sum = 0.08555
    myPSFN.max = 0.0004349  myPSFN.sum = 0.008295
    diff(RemiPSF - AnandPSF) = 0.0773
    diff(RemiPSF - AnandPSF)stdev = 4.52e-06

    """


    if NormPeak:
        print "Max of AS PSF = %.7g" % myPSF.max()
        print "Max of RS PSF = %.7g" % anRS_PSF.max()
        fac = anRS_PSF.max() / myPSF.max()
        print "Multiply AS PSF by a factor", fac
        myPSFN = myPSF * fac
    else:
        print "Tot of AS PSF = %.7g" % myPSF.sum()
        print "Tot of RS PSF = %.7g" % anRS_PSF.sum()
        fac = anRS_PSF.sum() / myPSF.sum()
        print "Multiply AS PSF by a factor", fac
        myPSFN = myPSF * fac

    #
    dPSF = anRS_PSF - myPSFN
    print "RemiPSF.max = %.4g  RemiPSF.sum = %.4g"  % \
    (anRS_PSF.max(), anRS_PSF.sum())
    print "myPSFN.max = %.4g  myPSFN.sum = %.4g"  % \
    (myPSFN.max(), myPSFN.sum())
    print "diff(RemiPSF - AnandPSF) = %.3g"  %  dPSF.sum()
    print "diff(RemiPSF - AnandPSF)stdev = %.3g"  %  dPSF.stddev()

    SimpleFitsWrite(fn=trn+"dyc.fits", data=dPSF.astype(Float32), clobber="y")

    rs_ap  = remiap
    rs_T   = remiT
    rs_fpm = remiFPM
    rs_ly  = remiLyot
    rs_6   = remiQn
    SimpleFitsWrite(fn=trn+"rs_ap.fits", data=rs_ap, clobber='y')
    SimpleFitsWrite(fn=trn+"rs_T.fits", data=rs_T, clobber='y')
    SimpleFitsWrite(fn=trn+"rs_fpm.fits", data=rs_fpm, clobber='y')
    SimpleFitsWrite(fn=trn+"rs_ly.fits", data=rs_ly, clobber='y')
    SimpleFitsWrite(fn=trn+"rs_yc.fits", data=remiPSF, clobber='y')
    SimpleFitsWrite(fn=trn+"rs_6.fits", data=rs_6, clobber='y')

    SimpleFitsWrite(fn=trn+"as_yc.fits", data=myPSFN.astype(Float32), clobber='y')

if __name__ == '__main__':

    #doit(coro="generate_first")
    doit(coro="None")

