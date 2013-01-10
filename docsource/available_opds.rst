.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Appendix: Available Optical Path Difference (OPD) files
================================================================


For each of the five instruments (four SIs plus FGS) there are three provided OPD files. These represent wavefronts as follows:

1. The OTE and ISIM intrinsic WFE
2. The above, plus a slight defocus to blur the image slightly to approximate image motion. 
3. The above #2, plus additional WFE due to SI internal optics. 

The latter is the largest WFE, and is the default file used in simulations unless another is explicitly chosen. For NIRCam only there is a second, duplicate 
set of these files with slightly improved WFE based on an optimistic case scenario for instrument and telescope alignment. 

Note that the trick of adding some (nonphysical) defocus to blur out the PSF is computationally easy and rapid, but does not give a high fidelity
representation of the true impact of image jitter. This is particularly true for coronagraphic observations. Future versions of WebbPSF will likely 
provide higher fidelity jitter models.

The units of the supplied OPD files are wavefront error in microns.


.. table:: Rev V OPDs

    =========================       ======= ================================================================================
                         File       RMS WFE                                                                         Contents
    =========================       ======= ================================================================================
        OPD_RevV_fgs_150.fits         150.0                OPD for OTE+ISIM with all reserves and stability; NO image motion
        OPD_RevV_fgs_163.fits         163.0 Observatory OPD for OTE+ISIM with all reserves and stability, and image motion as defocus
        OPD_RevV_fgs_186.fits         186.0 Observatory OPD for OTE+ISIM+SI with all reserves and stability, and image motion as defocus
       OPD_RevV_miri_204.fits         204.0                OPD for OTE+ISIM with all reserves and stability; NO image motion
       OPD_RevV_miri_220.fits         220.0 Observatory OPD for OTE+ISIM with all reserves and stability, and image motion as defocus
       OPD_RevV_miri_421.fits         421.0 Observatory OPD for OTE+ISIM+SI with all reserves and stability, and image motion as defocus
     OPD_RevV_nircam_115.fits         115.0                OPD for OTE+ISIM with all reserves and stability; NO image motion
     OPD_RevV_nircam_123.fits         123.0                OPD for OTE+ISIM with all reserves and stability; NO image motion
     OPD_RevV_nircam_132.fits         132.0 Observatory OPD for OTE+ISIM with all reserves and stability, and image motion as defocus
     OPD_RevV_nircam_136.fits         136.0 Observatory OPD for OTE+ISIM with all reserves and stability, and image motion as defocus
     OPD_RevV_nircam_150.fits         150.0 Observatory OPD for OTE+ISIM+SI with all reserves and stability, and image motion as defocus
     OPD_RevV_nircam_155.fits         155.0 Observatory OPD for OTE+ISIM+SI with all reserves and stability, and image motion as defocus
    OPD_RevV_nirspec_125.fits         125.0                OPD for OTE+ISIM with all reserves and stability; NO image motion
    OPD_RevV_nirspec_145.fits         145.0 Observatory OPD for OTE+ISIM with all reserves and stability, and image motion as defocus
    OPD_RevV_nirspec_238.fits         238.0 Observatory OPD for OTE+ISIM+SI with all reserves and stability, and image motion as defocus
         OPD_RevV_tf_144.fits         144.0                OPD for OTE+ISIM with all reserves and stability; NO image motion
         OPD_RevV_tf_162.fits         162.0 Observatory OPD for OTE+ISIM with all reserves and stability, and image motion as defocus
         OPD_RevV_tf_180.fits         180.0 Observatory OPD for OTE+ISIM+SI with all reserves and stability, and image motion as defocus
    =========================       ======= ================================================================================




--------------

Documentation last updated on |today|

