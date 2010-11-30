.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for JWST PSF Simulation
=====================================



Conceptually, the new JWST PSF simulation code relies on several layers of abstraction: 
 * A base layer implementing wavefront propagation through generic optical systems (provided by the Python module `poppy`).
 * An implementation of the specific details of JWST instruments using that base system (provided by `jwopt`)
 * And a graphical user interface (provided by `newgui`).

It is entirely possible (and indeed recommended for scripting) to just use the `jwopt` interface without the GUI, but the
GUI will provide a quicker method for simple interactive or exploratory calculations.


.. toctree::
   :maxdepth: 2

   poppy.rst
   jwopt.rst
   gui.rst  


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Documentation last updated on |today|

