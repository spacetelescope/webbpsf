.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=============================================
Physical Optics Propagation in PYthon (POPPY)
=============================================

.. module:: poppy

This module implements an object-oriented system for modeling physical optics
propagation with diffraction, particularly for telescopic and coronagraphic
imaging. Right now only image and pupil planes are supported; intermediate
planes are a future goal.

This code makes use of the python standard module 'logging' for output information. Top-level details of the
calculation are output at level logging.INFO, while details of the propagation through each optical plane are 
printed at level logging.DEBUG. See the Python logging documentation for an explanation of how to redirect the 'poppy' 
logger to the screen, a textfile, or any other log destination of your choice.


List of Classes
--------
 * :ref:`Wavefront`
 * OpticalElement

   * AnalyticOpticalElement

     * BandLimitedCoron
     * IdealMonoFQPM
     * IdealFieldStop
     * IdealCircularOcculter
   * Detector
 * OpticalSystem


.. Wavefront:

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

.. autoclass:: poppy.AnalyticOpticalElement
   :show-inheritance:
.. autoclass:: poppy.CircularAperture
   :show-inheritance:
.. autoclass:: poppy.IdealFieldStop
   :show-inheritance:
.. autoclass:: poppy.IdealCircularOcculter
   :show-inheritance:
.. autoclass:: poppy.IdealBarOcculter
   :show-inheritance:
.. autoclass:: poppy.BandLimitedCoron
   :show-inheritance:
.. autoclass:: poppy.IdealMonoFQPM
   :show-inheritance:


------

.. autoclass:: poppy.Detector
   :show-inheritance:

--------------

Documentation last updated on |today|

