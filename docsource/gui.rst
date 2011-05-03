.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: webbpsfgui
========================
Graphical User Interface
========================


.. autoclass:: newgui.JWPSF_GUI
   :members:


The WebbPSF GUI provides an easy interface to most, but not quite all, of the functionality of WebbPSF. 
(Not all of the advanced settings in the ``options`` structure for :py:class:`JWInstrument` are exposed in the GUI.)


To start the GUI:

>>> import webbpsf
>>> webbpsf.gui()

You can also just run directly from the commandline the ``webbpsfgui.py`` file contained within the webbpsf module. 


.. image:: ./fig_webbpsfgui_main.png
   :scale: 75%
   :align: right
   :alt: WebbPSFGui main window


The main window is divided into three regions:

* The top region allows control of the source spectral type and position. (Spectral type option requires optional install of `pysynphot` )
* The central, main region allows selection of instrument and configuration of instrument options. The options available here largely correspond to 
  attributes of the :py:class:`JWInstrument` classes.
* The lower region contains options for the PSF calculation itself such as pixel sampling and wavelengths. These correspond to parameters of the 
  :py:class:`JWInstrument.calcPSF` function call.


GUI Controls
--------------

The GUI buttons invoke actions as follows:


Compute PSF
^^^^^^^^^^^^

This invokes a PSF calculation with the given options. Each wavelength will be displayed in turn as it is computed, and finally the summed broadband PSF.
This resulting PSF is stored in memory for use by the next three buttons. 


Display PSF
^^^^^^^^^^^^
This button will redisplay the PSF if the window has closed or something else has been displayed.

.. image:: ./fig_display_psf.png
   :scale: 75%
   :align: right
   :alt: PSF display



Display Profiles
^^^^^^^^^^^^^^^^
This will display the PSF's radial profile and encircled energy profile.

.. image:: ./fig_display_profiles.png
   :scale: 75%
   :align: right
   :alt: PSF radial profiles display



Save PSF As...
^^^^^^^^^^^^^^

This will invoke a standard File Save dialog box allowing you to save your new PSF. 


Display Optics
^^^^^^^^^^^^^^


This will display a graphical representation of the optical train for the current instrument configuration.


.. image:: ./fig_nircam_coron_optics.png
   :scale: 75%
   :align: center
   :alt: Sample "Display Optics" results showing NIRCam coronagraphic optics.





--------------

Documentation last updated on |today|


