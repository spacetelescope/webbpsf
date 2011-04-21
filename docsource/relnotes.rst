.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Release Notes
================


Version 0.2
------------

* Much improved pysynphot support.
* Reworked calling conventions for calcPSF() routine source parameters.
* poppy.calcPSFmultiprocessor merged in to regular poppy.calcPSF
* Minor bug fixes to selection of which wavelengths to compute for more even sampling




Version 0.1
------------


Initial Release. Questions, comments, criticism all welcome!

**Included:**
 
* Revision V OPD files for OTE and SIs
* Support for imaging mode in all SIs and FGS
* Basic support for coronagraphy with MIRI, NIRCam, and TFI. Further enhancements in fidelity to come later.  Coronagraphic calculations are done using the direct FFT method, not Remi's semi-analytic method.
* Up-to-date science frame axes convention, including detector rotations for MIRI and NIRSpec.


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
* The NIRSpec and MIRI IFUs will be added in a future release
* Adjustible wavelengths and wavelength-dependent bandpass will be added to the TFI model.
* Realistic (but time consuming) jitter models
* Possibly: separate handling of pre- and post- coronagraphic WFE in instruments, if it appears likely to be significant. 

