Introduction
============

**What this software does:**

* Uses OPD maps precomputed by detailed optical simulations of JWST and the Nancy Grace Roman Space Telescope (formerly WFIRST), and in the case of JWST
  based on instrument and telescope flight hardware cryo-vacuum test results.
* For JWST, computes PSF images with requested properties for any of JWST's instruments. Supports imaging, coronagraphy, and most spectrographic modes with all of JWST's instruments. IFUs are yet to come.
* For Roman, computes PSFs with the Wide Field Imager, based on recent GSFC optical models, including field- and wavelength-dependent aberrations.
  A preliminary version of the Coronagraph Instrument is also available.
* Provides a suite of tools for quantifying PSF properties such as FWHM, Strehl ratio, etc.

**What this software does NOT do:**

* Contain in itself any detailed thermal or optical model of JWST or Roman. For the results of end-to-end integrated simulations of JWST, see for instance `Predicted JWST imaging performance (Knight, Lightsey, & Barto; Proc. SPIE 2012) <http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1362264>`_. For Roman modeling, see `the Roman Reference Info page <http://wfirst.gsfc.nasa.gov/science/Instrument_Reference_Information.html>`_
* Model spectrally dispersed PSFs produced by any of the spectrograph gratings. It does, however, let you produce monochromatic PSFs in these modes, suitable for stitching together into spectra using some other software.
* Model most detector effects such as pixel MTF, intrapixel sensitivity variations, interpixel capacitance, or any noise sources. Add those separately with your favorite detector model code. (\*Note, one particularly significant
  detector scattering for MIRI imaging has now been added.)


Conceptually, this simulation code has two layers of abstraction:
 * A base package for wavefront propagation through generic optical systems (provided by :py:mod:`POPPY <poppy>`)
 * Models of the JWST and Roman instruments implemented on top of that base system (provided by :py:mod:`WebbPSF <webbpsf>`)

.. _intro_why_webbpsf:

Why WebbPSF?
------------

For any space telescope, an ability to predict the properties of
point spread functions (PSFs) is needed before launch for a wide
range of preparatory science studies and tool development.
Tools for producing high
quality model PSFs must be easily accessible to the entire astronomical
community.
WebbPSF provides an easy-to-use tool for PSF simulations of JWST and Roman, in
the style of the highly successful "Tiny Tim" PSF simulator for Hubble.

WebbPSF
simulations are based on a mixture of observatory design parameters and
as-built properties. The software provides a highly flexible and scriptable toolkit in
Python for simulating a very wide range of observing modes and science scenarios, using
efficient computational methods (including optional parallelization and use of GPUs). WebbPSF
is a key building block in higher-level observatory simulators, including the
JWST `Exposure Time Calculator <https://jwst.etc.stsci.edu>`_.


.. _intro_algorithms:

Algorithms Overview
-------------------

Read on if you're interested in details of how the computations are performed. Otherwise, jump to :ref:`Quick Start <quickstart>`.

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

Types of Fourier Transform Calculation in WebbPSF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  * Any direct imaging calculation, any instrument: Matrix DFT
  * NIRCam coronagraphy with circular occulters: Semi-Analytic Fast Coronagraphy and Matrix DFT
  * NIRCam coronagraphy with wedge occulters: FFT and Matrix DFT
  * MIRI Coronagraphy: FFT and Matrix DFT
  * NIRISS NRM, GR799XD: Matrix DFT
  * NIRSpec and NIRISS slit spectroscopy: FFT and Matrix DFT

See :ref:`Optimizing Performance and Parallelization <performance_and_parallelization>` in the POPPY documentation for more details on calculation performance.

Getting WebbPSF
---------------

See :ref:`installation`.

.. _quickstart:

Quick Start
------------

Once you have installed the software and data files, we recommend you begin with the 
`Jupyter Notebook quickstart tutorial <http://nbviewer.jupyter.org/github/spacetelescope/webbpsf/blob/develop/notebooks/WebbPSF_tutorial.ipynb>`_. Downloading and running that notebook is a great way to get started using WebbPSF.


