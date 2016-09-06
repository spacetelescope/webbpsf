Documentation for WebbPSF 
===============================

WebbPSF is a Python package that computes simulated PSFs for the JWST instruments (and now for WFIRST too!), taking into account detector pixel scales, rotations, filter profiles, and point source spectra. It is *not* a full optical model of JWST, but rather a tool for transforming optical path difference (OPD) maps, created with some other tool, into the resulting PSFs as observed with JWST's instruments. 

.. figure:: ./fig_instrument_comparison.png
   :scale: 45 %
   :align: center
   :alt: Sample PSFs for JWST's instruments. 

   Sample PSFs for JWST's instrument suite, all on the same angular scale and display stretch.

.. figure:: ./wfirst_figures/webbpsf-wfirst_page_header.png
   :scale: 70 %
   :align: center
   :alt: Sample PSFs for the filters in the WFIRST WFI.

   Sample PSFs for the filters in the WFIRST WFI. 



**What this software does:**

* Uses OPD map(s) precomputed by detailed optical simulations of these telescopes.
* Computes from those PSF images with requested properties for any of JWST's instruments
* Supports imaging, coronagraphy, and most spectrographic modes with all of JWST's instruments. IFUs are yet to come.
* For WFIRST, computes PSFs with the Wide Field Imager, based on recent GSFC optical models, including field- and wavelength-dependent aberrations.
  CGI models may be included in a future release.
* Provides a suite of tools for quantifying PSF properties such as FWHM, Strehl ratio, etc.

**What this software does NOT do:**

* Contain in itself any detailed thermal or optical model of JWST or WFIRST. For the results of end-to-end integrated simulations of JWST, see for instance `Predicted JWST imaging performance (Knight, Lightsey, & Barto; Proc. SPIE 2012) <http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1362264>`_. For WFIRST modeling, see `the WFIRST Reference Info page <http://wfirst.gsfc.nasa.gov/science/Instrument_Reference_Information.html>`_
* Model spectrally dispersed PSFs produced by any of the spectrograph gratings. It does, however, let you produce monochromatic PSFs in these modes, suitable for stitching together into spectra using some other software.
* Model any detector effects such as pixel MTF, intrapixel sensitivity variations, interpixel capacitance, or any noise sources. Add those separately with your favorite detector model code.


Getting Started with WebbPSF
----------------------------

The WebbPSF software system is composed of two Python packages: a lower-level optical propagation library (:py:mod:`POPPY <poppy>`) plus an implementation of the JWST instruments using that library (:py:mod:`WebbPSF <webbpsf>`).  This documentation explains the programming interfaces and graphical user interface of WebbPSF, as well as providing a :ref:`quick overview <poppy_overiew>` of POPPY.

.. admonition:: Quickstart IPython Notebook

   This documentation is complemented by an `Jupyter (IPython) Notebook format quickstart tutorial <http://nbviewer.jupyter.org/github/mperrin/webbpsf/blob/master/notebooks/WebbPSF_tutorial.ipynb>`_. Downloading and running that notebook is a great way to get started using WebbPSF.

:ref:`What's new in the latest release? <whatsnew>`


.. toctree::
   :maxdepth: 1

   intro.rst
   installation.rst
   webbpsf.rst
   jwst.rst
   wfirst.rst
   gui.rst
   more_examples.rst
   poppy.rst

.. admonition:: Getting Help

   For help using or installing webbpsf, you can contact the STScI Help Desk, help@stsci.edu. Note that WebbPSF is included in the `Ureka <http://ssb.stsci.edu/ureka>`_ python distribution, as well as being installable via :ref:`standard Python packaging tools <installation>`. For detailed aspects of the JWST models, contact Marshall Perrin at STScI; for WFIRST, contact Joseph Long at STScI. 

Advanced Usage
--------------

.. toctree::
   :maxdepth: 1

   api_reference.rst
   help.rst
   performance.rst
   sampling.rst
   fft_optimization.rst

Appendices and Reference
------------------------

.. toctree::
   :maxdepth: 1

   available_opds.rst
   references.rst
   relnotes.rst
   release.rst


.. admonition:: How to cite WebbPSF

    In addition to this documentation, WebbPSF is described in the following references.  Users of WebbPSF are encouraged to cite one of these.
   
    * Perrin et al. 2014, `"Updated point spread function simulations for JWST with WebbPSF" <http://adsabs.harvard.edu/abs/2014SPIE.9143E..3XP>`_,  Proc. SPIE. 9143, 
    * Perrin et al. 2012, `"Simulating point spread functions for the James Webb Space Telescope with WebbPSF", <http://adsabs.harvard.edu/abs/2012SPIE.8442E..3DP>`_ Proc SPIE 8842, and 
    * Perrin 2011, :download:`Improved PSF Simulations for JWST: Methods, Algorithms, and Validation <Improved_PSFs_for_Webb.pdf>`, JWST Technical report JWST-STScI-002469.

    In particular, the 2012 SPIE paper gives a broad overview, the 2014 SPIE paper presents comparisons to instrument cryotest data, and the Technical Report document describes in more detail the relevant optical physics, explains design decisions and motivation for WebbPSF's architecture, and presents extensive validation tests demonstrating consistency between WebbPSF and other PSF simulation packages used throughout the JWST project.


* :ref:`genindex`
* :ref:`search`

**Mailing List**

If you would like to receive email announcements of future versions, please contact Marshall Perrin, or visit `maillist.stsci.edu <https://maillist.stsci.edu/scripts/wa.exe?SUBED1=webbpsf-users&A=1>` to subscribe yourself to the "webbpsf-users@maillist.stsci.edu" list.
