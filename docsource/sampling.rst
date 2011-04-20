
Appendix: Sampling Requirements for Numerical Accuracy
============================================================

The purpose of this appendix is to help you decide how many wavelengths and how much oversampling is required for your
particular science application. 

**Key Concepts:**


Obtaining high precision in PSF calculations requires treating both the multiwavelength 
nature of the selected bandpass and also the details of subpixel sampling and integration onto the detector pixels. 

*Note:* The current version of this code makes no attempt to incorporate detector effects such as pixel MTF and interpixel capacitance.
If you care about such effects, you should add them with another code. 

Multiwavelength effects scale the PSF linearly with respect to wavelength. Thus the absolute scale of this effect increases
linearly with distance from the PSF center. The larger a field of view you care about, the  more wavelengths you will need to include. 

Pixel sampling matters most near the core of the PSF, where the flux is changing very rapidly on small spatial scales. The closer 
to the core of the PSF you care about fine structure, the more finely sampled your PSF will need to be. 


**Some Useful Guidance:**


To evaluate what levels of sampling are needed in practice, for each NIRCam and MIRI filter we first computed a very highly oversampled image (nlambda=100, oversampling=16, field of view 5 arcsec), which we used as a "truth" image. 
(For practical purposes, we consider this level of sampling likely to be sufficiently fine that it's a good stand-in for an infinitely sampled PSF, but this is an assumption we have not quantitatively validated. )

We consider two types of measurement one might wish to make on the PSF: 
1) measuring the encircled energy curve to a given precision
2) measuring individual pixel flux levels to a given precision
The latter is substantially more challenging a measurement. 



Required sampling for NIRCam::

        Filter Name    SN=20 @r=1   SN=100 @r=1        SN=20 @r=3   SN=100 @r=3
               F070W    (2, 21)      (4, 50)            (2, 30)      (4, 75)
               F090W    (2, 21)      (4, 50)            (2, 21)      (4, 50)
               F115W    (2, 21)      (4, 50)            (2, 21)      (4, 75)
               F140M    (2, 13)      (4, 40)            (2, 21)      (4, 75)
              F150W2    (2, 21)      (4, 50)            (1, 75)      (2, 75)
               F150W    (2, 21)      (4, 50)            (2, 21)      (4, 75)
               F162M    (2, 13)      (4, 40)            (2, 21)      (4, 75)
               F164N    (4, 3)       (4, 3)             (2, 5)       (4, 13)
               F182M    (2, 13)      (4, 40)            (2, 21)      (4, 75)
               F187N    (2, 3)       (4, 3)             (2, 5)       (4, 13)
               F200W    (2, 21)      (4, 50)            (1, 40)      (4, 75)
               F210M    (2, 13)      (4, 40)            (2, 21)      (4, 75)
               F212N    (2, 3)       (4, 3)             (2, 3)       (4, 13)
               F225N    (2, 3)       (4, 3)             (2, 3)       (4, 9)
               F250M    (4, 5)       higher!            (2, 21)      (4, 50)
               F277W    (2, 13)      (4, 50)            (2, 21)      (4, 50)
               F300M    (2, 9)       (4, 30)            (2, 21)      (4, 75)
              F322W2    (2, 21)      (4, 50)            (1, 40)      (4, 50)
               F323N    (4, 3)       higher!            (2, 3)       (4, 5)
               F335M    (2, 9)       (4, 30)            (2, 13)      (4, 40)
               F356W    (2, 21)      (4, 50)            (1, 30)      (4, 50)
               F360M    (2, 9)       (4, 30)            (2, 21)      (4, 50)
               F405N    (4, 3)       higher!            (2, 3)       (4, 5)
               F410M    (2, 9)       (4, 21)            (2, 21)      (4, 75)
               F418N    (4, 3)       higher!            (2, 3)       (4, 5)
               F430M    (2, 3)       (4, 9)             (2, 13)      (4, 40)
               F444W    (2, 21)      (4, 50)            (1, 21)      (2, 75)
               F460M    (2, 3)       (4, 9)             (2, 9)       (4, 30)
               F466N    (2, 3)       higher!            (2, 3)       (4, 5)
               F470N    (2, 3)       (4, 3)             (2, 3)       (4, 5)
               F480M    (2, 5)       (4, 21)            (2, 21)      (4, 50)

And for MIRI::

         Filter Name    SN=20 @r=1   SN=100 @r=1        SN=20 @r=3   SN=100 @r=3
               F560W    (2, 9)       (4, 30)            (2, 13)      (4, 40)
               F770W    (2, 13)      (4, 40)            (1, 40)      (4, 50)
              F1000W    (2, 9)       (4, 21)            (1, 21)      (2, 50)
              F1065C    (2, 3)       (4, 3)             (1, 3)       (2, 9)
              F1130W    (2, 3)       (4, 5)             (1, 5)       (2, 13)
              F1140C    (2, 3)       (4, 3)             (1, 3)       (2, 9)
              F1280W    (1, 5)       (2, 13)            (1, 21)      (2, 50)
              F1500W    (1, 5)       (4, 21)            (1, 21)      (2, 75)
              F1550C    (2, 3)       (4, 3)             (1, 3)       (2, 9)
              F1800W    (1, 3)       (2, 13)            (1, 13)      (2, 40)
              F2100W    (1, 5)       (2, 21)            (1, 13)      (1, 40)
              F2300C    (1, 3)       (2, 13)            (1, 9)       (1, 30)
              F2550W    (1, 3)       (2, 9)             (1, 9)       (1, 21)


More later.
