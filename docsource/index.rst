.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for JWST PSF Simulation
=====================================

This software computes simulated PSFs for the JWST instruments, taking into account detector pixel scales, rotations, filter profiles, and point source spectra. 
It is *not* a full optical model of JWST, but rather a tool for transforming optical path difference (OPD) maps, created with some other tool, into the resulting
PSFs as observed with JWST's instruments.


Contents
--------

.. toctree::
   :maxdepth: 2

   intro.rst
   relnotes.rst
   jwopt.rst
   gui.rst  
   poppy.rst
   fft_optimization.rst


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Documentation last updated on |today|

