Sampling Requirements for Numerical Accuracy
============================================================

The purpose of this appendix is to help you decide how many wavelengths and how much oversampling is required for your
particular science application.

Key Concepts
-----------------------------------------


Obtaining high accuracy and precision in PSF calculations requires treating both the multiwavelength
nature of the selected bandpass and also the details of subpixel sampling and integration onto the detector pixels.

*Note:* The current version of this code makes no attempt to incorporate detector effects such as pixel MTF and interpixel capacitance.
If you care about such effects, you should add them with another code.

Multiwavelength effects scale the PSF linearly with respect to wavelength. Thus the absolute scale of this effect increases
linearly with distance from the PSF center. The larger a field of view you care about, the more wavelengths you will need to include.

Pixel sampling matters most near the core of the PSF, where the flux is changing very rapidly on small spatial scales. The closer
to the core of the PSF you care about fine structure, the more finely sampled your PSF will need to be.


Some Useful Guidance
-------------------------------



We consider two types of measurement one might wish to make on the PSF:

1. measuring the encircled energy curve to a given precision
2. measuring individual pixel flux levels to a given precision

The latter is substantially more challenging a measurement. The below tables
present the number of (oversamplings, wavelengths) needed to achieve SNR=100 in
a single pixel at a given radius (where SNR in this context is calculated as
``(image-truth)/truth`` on a per-detector-pixel basis and then averaged in an annulus as
a function of radius).  This calculation is motivated by modeling coronagraphic
PSF subtraction, where we might hope to achieve 1-2 orders of magnitude
reduction in the PSF wings through PSF subtraction. Accurately simulating that
process demands a comparable level of fidelity in our PSF models. We also present tables giving the
requirements for SNR=20 in a given pixel for less demanding modeling tasks.


Note that we do not consider here the case of trying to model the PSF core at SNR=100/pixel. If you are
interested in doing so, I believe very fine subsampling would be needed. This might be most efficiently computed using a highly oversampled PSF for just the core, glued in to a larger image computed at lower angular resolution for the rest of the field of view. Investigating this is left as an exercise for another day.



Because NIRSpec, NIRISS, and FGS sample the PSF relatively coarsely, they will require a higher degree of oversampling in simulations than NIRCam to reach a given SNR level.
MIRI is fairly well-sampled.



Per-Instrument Sampling Requirements
-----------------------------------------

To evaluate what levels of sampling are needed in practice, for each NIRCam and MIRI filter we first computed a very highly oversampled image (nlambda=200, oversampling=16; field of view 5 arcsec for NIRCam and 12 arcsec for MIRI), which we used as a "truth" image.
(For practical purposes, we consider this level of sampling likely to be sufficiently fine that it's a good stand-in for an infinitely sampled PSF, but this is an assumption we have not quantitatively validated. However, since there are >200 subsamples in both pixel and wavelength space, the residuals ought to be <1/200 and thus these are sufficient for our purposes of testing SNR=100.)


These tables list the (oversampling, wavelengths) required to achieve the
specified SNR, in comparison with a 'truth' image based on simulations using
``oversampling = 16`` (i.e. 256 subpixels per detector pixel) and
``nlambda=200``.

Required sampling for NIRCam::

        NIRCam, SNR=100
                      r=0.5"       1.0"       2.0"       3.0"
             F070W    higher!    (4, 13)    (4, 21)    (4, 30)
             F090W    higher!    (4, 13)    (4, 21)    (4, 30)
             F115W    higher!     (4, 9)    (4, 21)    (4, 30)
             F140M    higher!     (4, 9)     (4, 9)    (4, 13)
            F150W2    higher!    (4, 30)    (2, 75)    (2, 75)
             F150W    higher!     (4, 9)    (4, 21)    (4, 21)
             F162M    higher!     (4, 9)     (4, 9)    (4, 13)
             F164N    higher!     (8, 3)     (8, 3)     (8, 3)
             F182M    higher!     (4, 9)     (4, 9)    (4, 13)
             F187N    higher!     (8, 1)     (4, 5)     (8, 3)
             F200W     (8, 5)     (4, 9)    (2, 21)    (2, 30)
             F210M     (8, 3)     (4, 5)     (4, 9)     (4, 9)
             F212N     (8, 1)     (4, 3)     (4, 3)    (4, 13)
             F225N     (8, 1)     (4, 3)     (4, 3)     (4, 5)
             F250M    higher!     (8, 5)    (4, 13)     (4, 9)
             F277W     (8, 5)     (4, 9)    (4, 13)    (4, 21)
             F300M     (8, 3)     (8, 5)     (4, 9)     (4, 9)
            F322W2     (8, 9)    (4, 21)    (4, 21)    (4, 30)
             F323N     (8, 1)     (8, 1)     (8, 3)     (8, 3)
             F335M     (8, 3)     (8, 5)     (4, 9)     (4, 9)
             F356W     (8, 5)     (4, 9)     (4, 9)    (4, 13)
             F360M     (8, 3)     (8, 5)     (4, 5)     (4, 9)
             F405N     (8, 1)     (8, 1)     (4, 9)     (8, 3)
             F410M     (8, 3)     (8, 5)     (4, 5)     (4, 9)
             F418N     (8, 1)     (8, 1)     (4, 5)     (8, 3)
             F430M     (8, 1)     (8, 3)     (4, 9)     (4, 9)
             F444W     (8, 5)     (4, 9)    (4, 13)    (2, 21)
             F460M     (8, 3)     (8, 5)     (4, 9)     (4, 9)
             F466N     (8, 1)     (8, 1)     (4, 3)     (4, 9)
             F470N     (8, 1)     (8, 1)     (4, 3)     (4, 3)
             F480M     (8, 3)    (4, 21)     (4, 5)     (4, 9)

        NIRCam, SNR=20
                      r=0.5"       1.0"       2.0"       3.0"
             F070W     (8, 3)     (2, 9)    (2, 21)    (2, 21)
             F090W     (8, 3)     (2, 9)    (2, 13)    (2, 21)
             F115W     (8, 3)     (2, 9)    (2, 13)    (2, 21)
             F140M     (8, 3)     (2, 5)     (2, 5)     (2, 9)
            F150W2     (8, 9)    (2, 21)    (1, 50)    (1, 75)
             F150W     (8, 3)     (2, 9)    (2, 13)    (2, 21)
             F162M     (8, 3)     (2, 5)     (2, 5)     (2, 9)
             F164N     (8, 1)     (4, 1)     (2, 3)     (4, 3)
             F182M     (8, 3)     (2, 3)     (2, 5)     (2, 9)
             F187N     (8, 1)     (4, 1)     (2, 3)     (2, 5)
             F200W     (4, 3)     (2, 5)    (1, 13)    (1, 21)
             F210M     (4, 3)     (2, 3)     (2, 5)     (2, 5)
             F212N     (4, 1)     (2, 1)     (2, 3)     (2, 3)
             F225N     (4, 1)     (2, 1)     (2, 3)     (2, 3)
             F250M     (8, 1)     (4, 3)     (2, 5)     (2, 5)
             F277W     (4, 3)     (2, 5)     (2, 9)    (2, 13)
             F300M     (4, 3)     (4, 3)     (2, 5)     (2, 5)
            F322W2     (4, 5)     (2, 9)    (2, 21)    (2, 21)
             F323N     (4, 1)     (4, 1)     (4, 1)     (4, 1)
             F335M     (4, 3)     (4, 3)     (2, 5)     (2, 5)
             F356W     (4, 3)     (2, 5)     (2, 9)     (2, 9)
             F360M     (4, 3)     (4, 3)     (2, 3)     (2, 5)
             F405N     (4, 1)     (4, 1)     (2, 3)     (2, 3)
             F410M     (4, 1)     (4, 3)     (2, 3)     (2, 5)
             F418N     (4, 1)     (4, 1)     (2, 1)     (4, 1)
             F430M     (4, 1)     (4, 3)     (2, 3)     (2, 5)
             F444W     (4, 3)     (2, 5)     (1, 9)    (1, 13)
             F460M     (4, 1)     (2, 5)     (2, 3)     (2, 5)
             F466N     (4, 1)     (4, 1)     (2, 1)     (2, 3)
             F470N     (4, 1)     (2, 3)     (2, 1)     (2, 1)
             F480M     (4, 1)     (2, 5)     (2, 3)     (2, 3)


We have not yet performed simulations for the case of NIRISS. The number of wavelengths used for each filter is set equal to
that used for NIRCam. This should certainly be adequate for the long-wavelength filters (given the NIRISS detector and NIRCam LW are
identical) but users may wish to investigate using finer sampling for the shorter wavelength filters that are very undersampled on NIRISS.


And for MIRI::


    MIRI, SNR = 100
                       r=1.0"      2.0"       3.5"       5.0"
             F560W     (4, 5)     (4, 9)    (4, 13)    (4, 13)
             F770W     (4, 5)     (2, 9)    (2, 13)    (2, 21)
            F1000W     (4, 3)     (4, 5)     (2, 9)     (2, 9)
            F1065C     (4, 3)     (4, 5)     (4, 5)     (2, 5)
            F1130W     (4, 3)     (4, 5)     (2, 5)     (2, 5)
            F1140C     (4, 3)     (4, 3)     (4, 5)     (2, 5)
            F1280W     (4, 3)     (2, 5)     (2, 9)     (2, 9)
            F1500W     (4, 3)     (2, 5)     (2, 9)     (2, 9)
            F1550C     (4, 3)     (2, 3)     (2, 3)     (2, 5)
            F1800W     (2, 3)     (2, 3)     (2, 9)     (2, 9)
            F2100W     (2, 3)     (2, 5)     (2, 9)     (1, 9)
            F2300C     (2, 3)     (2, 5)     (1, 9)     (1, 9)
            F2550W     (2, 3)     (1, 5)     (1, 9)     (1, 9)
               FND    (2, 30)    (2, 40)    (2, 50)    (2, 75)


    MIRI, SNR=20
                       r=1.0"      2.0"       3.5"       5.0"
             F560W     (2, 3)     (2, 5)     (2, 9)     (2, 9)
             F770W     (2, 3)     (1, 9)     (1, 9)     (1, 9)
            F1000W     (2, 3)     (1, 5)     (1, 5)     (1, 5)
            F1065C     (2, 1)     (2, 3)     (2, 3)     (1, 3)
            F1130W     (2, 1)     (2, 3)     (1, 3)     (1, 3)
            F1140C     (2, 1)     (2, 3)     (1, 3)     (1, 3)
            F1280W     (2, 3)     (1, 3)     (1, 5)     (1, 5)
            F1500W     (2, 3)     (1, 3)     (1, 5)     (1, 5)
            F1550C     (2, 1)     (1, 3)     (1, 3)     (1, 3)
            F1800W     (1, 3)     (1, 3)     (1, 5)     (1, 5)
            F2100W     (1, 3)     (1, 3)     (1, 5)     (1, 5)
            F2300C     (1, 3)     (1, 3)     (1, 5)     (1, 5)
            F2550W     (1, 3)     (1, 3)     (1, 3)     (1, 5)
               FND    (1, 13)    (1, 21)    (1, 40)    (1, 50)


The defaults for MIRI are set to 9 wavelengths for all filters, except for F560W and F770W which use 13 and FND which uses 40.

More later.
