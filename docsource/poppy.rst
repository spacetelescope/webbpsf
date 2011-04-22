.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: poppy

=============================================
Physical Optics Propagation in PYthon (POPPY)
=============================================


**Introduction:**

This module implements an object-oriented system for modeling physical optics
propagation with diffraction, particularly for telescopic and coronagraphic
imaging. Right now only image and pupil planes are supported; intermediate
planes are a future goal.

This code makes use of the python standard module 'logging' for output information. Top-level details of the
calculation are output at level logging.INFO, while details of the propagation through each optical plane are 
printed at level logging.DEBUG. See the Python logging documentation for an explanation of how to redirect the 'poppy' 
logger to the screen, a textfile, or any other log destination of your choice.


**Key Concepts:**

To model optical propagation, Poppy implements an object-oriented system for
representing an optical train. There are a variety of :py:class:`OpticalElement` classes
representing both physical elements as apertures, mirrors, and apodizers, and
also implicit operations on wavefronts, such as rotations or tilts. Each
:py:class:`OpticalElement`  may be defined either via analytic functions (e.g. a simple
circular aperture) or by numerical input FITS files (e.g. the complex JWST
aperture with appropriate per-segment WFE). A series of such :py:class:`OpticalElements <OpticalElement>` is
chained together in order in an :py:class:`OpticalSystem` class. That class is capable of generating
:py:class:`Wavefronts <Wavefront>`  (another class) suitable for propagation through the desired elemeents 
(with correct array size and sampling), and then performing the optical propagation onto
the final image plane. 

Poppy presently assumes that these propagation steps can be modeled using Fraunhofer diffraction (far-field), such that
the relationship between pupil and image plane optics is given by two-dimensional Fourier transforms. Fresnel propagation is
not currently supported. 

Two different algorithmic flavors of Fourier transforms are used in Poppy. The
familiar FFT algorithm is used for transformations between pupil and image
planes in the general case. This algorithm is relatively fast (*O(N log(N))*)
but imposes strict constraints on the relative sizes and samplings of pupil and
image plane arrays. Obtaining fine sampling in the image plane requires very
large oversized pupil plane arrays and vice versa, and image plane pixel
sampling becomes wavelength dependent. To avoid these constraints, for
transforms onto the final :py:class:`Detector` plane, instead a Matrix Fourier Transform
(MFT) algorithm is used (See `Soummer et al. 2007 Optics Express <http://adsabs.harvard.edu/abs/2007OExpr..1515935S>`_).  This allows
computation of the PSF directly on the desired detector pixel scale or an
arbitrarily finely subsampled version therof. For equivalent array sizes *N*,
the MFT is slower than the FFT(*O(N^3)*), but in practice the ability to freely
choose a more appropriate *N* (and to avoid the need for post-FFT interpolation
onto a common pixel scale) more than makes up for this and the MFT is faster.




List of Classes
--------

.. inheritance-diagram:: poppy.Detector poppy.Wavefront poppy.OpticalSystem poppy.Rotation poppy.CircularAperture poppy.HexagonAperture poppy.SquareAperture poppy.IdealFieldStop poppy.IdealCircularOcculter poppy.IdealBarOcculter poppy.BandLimitedCoron poppy.IdealFQPM poppy.FQPM_FFT_aligner poppy.CompoundAnalyticOptic



.. _Wavefront:

Wavefront
---------

.. autoclass:: poppy.Wavefront
    :members:

.. OpticalSystem:

Optical System 
--------------

.. autoclass:: poppy.OpticalSystem
    :members:

.. OpticalElement:

Optical Elements
----------------

.. autoclass:: poppy.OpticalElement
   :members:

------

.. autoclass:: poppy.Rotation
   :show-inheritance:

.. autoclass:: poppy.AnalyticOpticalElement
   :show-inheritance:
.. autoclass:: poppy.CircularAperture
   :show-inheritance:
.. autoclass:: poppy.HexagonAperture
   :show-inheritance:
.. autoclass:: poppy.SquareAperture
   :show-inheritance:
.. autoclass:: poppy.IdealFieldStop
   :show-inheritance:
.. autoclass:: poppy.IdealCircularOcculter
   :show-inheritance:
.. autoclass:: poppy.IdealBarOcculter
   :show-inheritance:
.. autoclass:: poppy.BandLimitedCoron
   :show-inheritance:
.. autoclass:: poppy.IdealFQPM
   :show-inheritance:
.. autoclass:: poppy.FQPM_FFT_aligner
   :show-inheritance:
.. autoclass:: poppy.CompoundAnalyticOptic
   :show-inheritance:


------

.. autoclass:: poppy.Detector
   :show-inheritance:

--------------

Documentation last updated on |today|

