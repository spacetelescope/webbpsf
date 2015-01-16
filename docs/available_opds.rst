Appendix: Available Optical Path Difference (OPD) files
================================================================

For each of the five instruments (four SIs plus FGS) there are three provided OPD files. These represent wavefronts as follows:

1. The OTE and ISIM intrinsic WFE
2. The above, plus a slight defocus to blur the image slightly to approximate image motion. 
3. The above #2, plus additional WFE due to SI internal optics. 

The latter is the largest WFE, and is the default file used in simulations unless another is explicitly chosen. For NIRCam only there is a second, duplicate set of these files with slightly improved WFE based on an optimistic case scenario for instrument and telescope alignment.

The provided OPDs are based on the observatory design requirements, and were developed for the Mission Critical Design Review. The represent the nominal case of performance for JWST, and have not yet been updated with as-built details of mirror surface figures, etc. We intend to make updated OPD files available once suitable reference data have been provided to STScI. For now, see `Lightsey et al. 2014 <http://adsabs.harvard.edu/abs/2014SPIE.9143E..04L>`_ for recent predictions of JWST's likely performance

Note that the trick of adding some (nonphysical) defocus to blur out the PSF is computationally easy and rapid, but does not give a high fidelity
representation of the true impact of image jitter. This is particularly true for coronagraphic observations. Future versions of WebbPSF will likely 
provide higher fidelity jitter models.

The units of the supplied OPD files are wavefront error in microns.

.. table:: Rev V OPDs

    =========================  ==========     =======  ========================  ==========================  =======
                         File  Instrument     RMS WFE  Includes OTE + ISIM OPD?  Image motion (as defocus)?  SI OPD?
    =========================  ==========     =======  ========================  ==========================  =======
        OPD_RevV_fgs_150.fits  FGS              150.0  Yes                       No                          No
        OPD_RevV_fgs_163.fits  FGS              163.0  Yes                       Yes                         No
        OPD_RevV_fgs_186.fits  FGS              186.0  Yes                       Yes                         Yes
       OPD_RevV_miri_204.fits  MIRI             204.0  Yes                       No                          No
       OPD_RevV_miri_220.fits  MIRI             220.0  Yes                       Yes                         No
       OPD_RevV_miri_421.fits  MIRI             421.0  Yes                       Yes                         Yes
     OPD_RevV_nircam_115.fits  NIRCam           115.0  Yes, optimistic case      No                          No
     OPD_RevV_nircam_123.fits  NIRCam           123.0  Yes                       No                          No
     OPD_RevV_nircam_132.fits  NIRCam           132.0  Yes, optimistic case      Yes                         No
     OPD_RevV_nircam_136.fits  NIRCam           136.0  Yes                       Yes                         No
     OPD_RevV_nircam_150.fits  NIRCam           150.0  Yes, optimistic case      Yes                         Yes
     OPD_RevV_nircam_155.fits  NIRCam           155.0  Yes                       Yes                         Yes
    OPD_RevV_nirspec_125.fits  NIRSpec          125.0  Yes                       No                          No
    OPD_RevV_nirspec_145.fits  NIRSpec          145.0  Yes                       Yes                         No
    OPD_RevV_nirspec_238.fits  NIRSpec          238.0  Yes                       Yes                         Yes
     OPD_RevV_niriss_144.fits  NIRISS           144.0  Yes                       No                          No
     OPD_RevV_niriss_162.fits  NIRISS           162.0  Yes                       Yes                         No
     OPD_RevV_niriss_180.fits  NIRISS           180.0  Yes                       Yes                         Yes
    =========================  ==========     =======  ========================  ==========================  =======
