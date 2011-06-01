.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Introduction
============


Conceptually, the new JWST PSF simulation code has three layers of abstraction: 
 * A base package implements wavefront propagation through generic optical systems (provided by the Python module :py:mod:`POPPY <poppy>`).
 * The specific details of JWST instruments are then implemented using that base system (provided by :py:mod:`WebbPSF <webbpsf>`)
 * And there is a graphical user interface (provided by  :py:mod:`WebbPSFgui <webbpsfgui>`).

It is entirely possible (and indeed recommended for scripting) to just use the :py:mod:`WebbPSF <webbpsf>` interface without the GUI, but the
GUI will provide a quicker method for simple interactive calculations.


Why a new JWST PSF Simulator?
-----------------------------

Given that the ``JWPSF`` package has been available for several years now, one might ask why do we need a new PSF simulator? 
From a user's perspective this new code provides the following enhancements:

* Updated to the most recent JWST pupil and OPD models, Revision V.
* Added TFI and FGS models.
* Updated lists of available filters.
* Added support for coronagraphic observations with MIRI, NIRCam, and TFI. (Note that MIRI coronagraphy models were
  already available using the ``JWcorPSF`` code split from ``JWPSF``, but with substantial limitations on computation such as
  a fixed oversampling factor. NIRCam and TFI coronagraphy were not supported)
* Includes the detector rotations, particularly for MIRI and NIRSpec
* Adds ability to set output image FOV size and pixel sampling, separate from the oversampling factor used for the optical propagation.
* New & improved graphical user interface.


Perhaps even more importantly, the underlying codebase has been entirely replaced and revamped. The most 
significant additions from a programmer's perspective include:

* Much cleaner object-oriented interface. Better abstraction of details across layers.
* Support for optics defined by analytic functions
* Support for coordinate rotations and rotated optics.
* Arbitrary oversampling for coronagraphic models.
* Matrix Fourier Transform algorithm from Soummer et al. implemented for arbitrary detector sampling.
* Uses ``FFTW3`` library for improved speed and efficient use of multiple processor cores. 
* Uses ``pysynphot`` library (same as the HST & Webb exposure time calculators) for consistent treatment of filter bandpasses and source spectra.


Algorithm Overview
---------------------

Most users may skip this section; read on only if you are interested in details of how the computations are performed. Otherwise, jump to :ref:`Quick Start <quickstart>`

The problem at hand is to transform supplied, precomputed OPDs (derived from a detailed optomechanical model
of the telescope)
into observed PSFs as seen with one or more of JWST's various detectors. This requires knowledge of the 
location and orientation of the detector planes, the properties of relevant optics such as bandpass filters and/or
coronagraphic image and pupil plane masks, and a model of light propagation between them.

Instrumental properties are taken from project documentation and the published
literature as appropriate; see the :ref:`References <references>` for detailed
provenance information. Optics may be described either numerically (for
instance, a FITS file containing a mask image for a Lyot plane or a FITS
bintable giving a spectral bandpass) or analytically (for instance, a
coronagraph occulter described as a circle of a given radius or a band-limited
mask function with given free parameters). 


WebbPSF computes PSFs under the assumption that JWST's instruments are well
described by Fraunhofer diffraction, as implemented using the usual Fourier
relationship between optical pupil and image planes (e.g. `Goodman et al. 1996
<http://books.google.com/books?id=ow5xs_Rtt9AC&printsec=frontcover#v=onepage&q&f=false>`_).
Two specific types of 2D Fourier transform are implemented: a Fast Fourier Transform and a discrete Matrix Fourier Transform.

The familiar Fast Fourier Transform (FFT) algorithm achieves its speed at the cost of imposing a specific fixed relationship between pixel
sampling in the pupil and image planes. As a result, obtaining finely sampled PSFs requires transforming very large arrays consisting
mostly of zero-padding. A more computationally attractive method is to use a discrete matrix Fourier transform, which
provides flexibility to compute PSFs on any desired output sampling without requiring any excess padding of the input arrays.
While this algorithm's computational cost grows as `O(N^3)` versus `O(N log N)` for the FFT, the FFT's apparent advantage is immediately lost
due to the need to resample the output onto the real pixel grid, which is an `O(N^2)` operation. By performing a matrix fourier transform 
directly to the desired output pixel scale, we can achieve arbitrarily fine sampling without the use of memory-intensive large padded arrays, and 
with lower overall computation time.

Further optimizations are available in coronagraphic mode using the semi-analytic coronagraphic propagation algorithm of Soummer et al. 2007. In this approach, rather than
propagating the entire wavefront from pupil to image and back to pupil in order to account for the coronagraphic masks, we can propagate only the subset of the wavefront that
is actually blocked by the image occulter and then subtract it from the rest of the wavefront at the Lyot plane. This relies on Babinet's principle to achieve the same final PSF
with more computational efficiency, particularly for the case of highly oversampled image planes (as is necessary to account for fine structure in image plane occulter masks). See Soummer et al. 2007 for a detailed description of this algorithm.






.. _quickstart:

Quick Start
------------
First, download and install the software (as described in the next page of this document).  Then just

>>> import webbpsf
>>> webbpsf.gui()

and you should be able to test drive things using the GUI: 

.. image:: ./fig_webbpsfgui_main.png
   :scale: 75%
   :align: right
   :alt: WebbPSFGui main window







