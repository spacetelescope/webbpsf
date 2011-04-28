.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Introduction
============


Conceptually, the new JWST PSF simulation code has three layers of abstraction: 
 * A base package implements wavefront propagation through generic optical systems (provided by the Python module :py:mod:`POPPY <poppy>`).
 * An implementation of the specific details of JWST instruments using that base system (provided by :py:mod:`WebbPSF <webbpsf>`)
 * And a graphical user interface (provided by  :py:mod:`WebbPSFgui <webbpsfgui>`).

It is entirely possible (and indeed recommended for scripting) to just use the :py:mod:`WebbPSF <webbpsf>` interface without the GUI, but the
GUI will provide a quicker method for simple interactive exploratory calculations.


Why a new JWST PSF Simulator?
-----------------------------

Given that the JWPSF code has been available for several years now, one might ask why do we need a new PSF simulator? 
From a user's perspective this new code provides the following enhancements:

* Updated to the most recent JWST pupil and OPD models, Revision V.
* Added TFI and FGS.
* Updated lists of available filters.
* Ability to simulate coronagraphic observations with MIRI, NIRCam, and TFI. (Note that MIRI coronagraphy models were
  already available using the JWcorPSF code split from JWPSF, but with substantial limitations on computation such as
  a fixed oversampling factor. NIRCam and TFI coronagraphy were not supported)
* Includes the detector rotations, particularly for MIRI and NIRSpec
* Adds ability to set output image FOV size and pixel sampling, separate from the oversampling factor used for the optical propagation.
* Improved graphical user interface


Perhaps even more importantly, the underlying codebase has been entirely replaced and revamped. The most 
significant additions from a programmer's perspective include:

* Much cleaner object-oriented interface. Better abstraction of details across layers.
* Support for optics defined by analytic functions
* Support for coordinate rotations and rotated optics.
* Arbitrary oversampling for coronagraphic models.
* Matrix Fourier Transform algorithm from Soummer et al. implemented for arbitrary detector sampling
* Uses FFTW3 library for improved speed and efficient use of multiple processor cores. 

--------------



