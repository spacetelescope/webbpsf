.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Release Notes
================


Known Issues
--------------
* You may see various warning messages while running computations, like thus::

    No handlers could be found for logger "webbpsf"

    ((<pysynphot.spectrum.Box object at 0x1047132d0> * nircam,im,f200w)) does not have a defined 
    binset in the wavecat table. The waveset of the spectrum will be used instead.

    Warning: invalid value encountered in absolute

  These can safely be ignored. 
* The current development version of pysynphot prints text like the following when it starts up. This, too, can just be ignored. ::

    GRAPHTABLE:  /Users/mperrin/data/CDBS/mtab/z999999zz_tmg.fits
    COMPTABLE:  /Users/mperrin/data/CDBS/mtab/z999999zz_tmc.fits
    NOT DEFAULT -- Primary mirror area:  253260


Versions 0.2.1 - 0.2.3
-----------------------

* Smoother installation process (thanks to Anand Sivaramakrishan for initial testing)
* Semi-analytic coronagraphic algorithm added for TFI and NIRCam (Soummer et al. 2007)
* Advanced settings dialog box added to GUI
* NIRCam pixel scale auto-switching will no longer override custom user pixelscales.
* slight fix to pupil file pixel scales to reflect JWST flat-to-flat diameter=6.559 m rather than just "6.5m"
* Corrected NIRCam 430R occulter profile to exactly match flight design; other occulters still need to be tuned. Corrected for use of amplitude rather than intensity profiles (thanks to John Krist for comparison models). 
* added TFI NRM mode (thanks to Anand Sivaramakrishnan)




Version 0.2
------------

Initial Release, spring 2011. Questions, comments, criticism all welcome!

* Much improved pysynphot support.
* Reworked calling conventions for calcPSF() routine source parameters.
* poppy.calcPSFmultiprocessor merged in to regular poppy.calcPSF
* Minor bug fixes to selection of which wavelengths to compute for more even sampling
* Default OPDs are now the ones including SI WFE as well as OTE+ISIM.
* Improved fidelity for NIRCam coronagraphic occulter models including ND squares and substrate border.




Version 0.1
------------

Development, fall 2010.


**Included:**
 
* Revision V OPD files for OTE and SIs. Produced by Kong Ha at NASA GSFC, provided by Mark Clampin.
* Support for imaging mode in all SIs and FGS
* Basic support for coronagraphy with MIRI, NIRCam, and TFI. Further enhancements in fidelity to come later.  Coronagraphic calculations are done using the direct FFT method, not Soummer's semi-analytic method (though that may be implemented in the future?).
* Up-to-date science frame axes convention, including detector rotations for MIRI and NIRSpec.
* Tunable wavelengths and appropriate bandwidths for TFI.
* Partial support for modeling IFU PSFs through use of the 'monochromatic' parameter.


**The following factors are NOT included in these simulations:**

* PSF variations across the field of view of any instrument (though each one has its own distinct OPDs for the center of its FOV).
* Optical distortions.
* Any and all detector effects, including intrapixel sensitivity variations. There is no plan to include these at any point. Generate a subsampled PSF and use a separate detector model code instead. 
* Instrumental wavefront errors are not modeled separately, though they are included in some of the supplied RevV OPDs. 
* Coronagraphic masks are assumed to be perfect (i.e. the masks exactly match their design parameters.)
* TFI NRM mode.
* Edge effects near the center of the FQPMs.


Plans for Future Releases
--------------------------
* Full support for the NIRSpec and MIRI IFUs will be added in a future release
* Realistic (but time consuming) jitter models
* Possibly: separate handling of pre- and post- coronagraphic WFE in instruments, if it appears likely to be significant. 

