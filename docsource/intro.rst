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

Software Requirements
-----------------------

Beyond the usual numpy/scipy/matplotlib core modules, the following are required:

* `pyfits`
* `ATPy`

These are optional but highly recommended:

* `pyFFTW3`. The code will work fine without it, but will be significantly slower.
* `pysynphot` (The updated, prerelease version for JWST!) enabled the 
   simulation of PSFs with proper spectral response to realistic source spectra.  

**Additional Requirements for the GUI**: 

* The graphical user interface requires the `ttk` enhanced version of the `Tkinter` widget library. 

`ttk` is not included by default on some installations of Python, for instance the default Mac OS Python 2.6 install. 
You may wish to either upgrade to a more current Python, or else compile and install `ttk` for your platform. This code
was developed using Python 2.7, which includes `ttk` by default, but it ought to work fine on any installations of
Python 2.5 or 2.6 provided `ttk` is available. Alternatively, you can just skip using the GUI; the optical modeling classes
themselves have no dependency on these widgets.

--------------



