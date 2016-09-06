#############
Release Notes
#############

.. _known_issues:

Known Issues
--------------

See https://github.com/mperrin/webbpsf/issues for currently open issues and enhancement suggestions.

* Calculations at large radii (> 500 lambda/D ~ 30 arcsec for 2 microns) will show numerical artifacts from Fourier aliasing and the implicit repetition of
  the pupil entrance aperture in the discrete Fourier transform. If you need accurate PSF information at such large radii, please contact Marshall Perrin for
  higher resolution pupil data.


**The following factors are NOT included in these simulations:**

* Wavefront error (OPD) variations across the field of view of any JWST instrument (though each one has its own distinct OPDs for the center of its FOV).
  Field-dependent aberrations *are* included for WFIRST WFI.
* Optical distortions.
* Instrumental wavefront errors are not modeled separately, though they are included in some of the supplied RevV OPDs.
* Coronagraphic masks are assumed to be perfect (i.e. the masks exactly match their design parameters.)
* No edge effects near the center of the FQPMs. (However, these are believed to be negligible in practice based on detailed simulations by Remi Soummer.)
* Any and all detector effects, including intrapixel sensitivity variations. There is no plan to include these at any point in WebbPSF itself.  Generate a subsampled PSF and use a separate detector model code instead.

Road Map for Future Releases
--------------------------------
* Updates to JWST instrument and telecope models based on as-built cryotest data (expected mid 2016).
* Field dependence of PSFs over JWST instrument fields of view based on ISIM CV test data.
* Web interface based on Jupyter Notebook servers (first version available already for WFIRST)
* Improved spectroscopic simulations including prism/grating dispersions.
* Support for the NIRSpec and MIRI IFUs may be added in a future release, level of detail is still TBD.
* Improved models for pointing jitter.
* Possibly: separate handling of pre- and post- coronagraphic WFE in instruments, if this appears likely to be significant.

.. _whatsnew:

Version History and Change Log
-------------------------------

Version 0.5.1
=============

.. _rel0.5.1:

*Unreleased*

Version 0.5.0
=============

.. _rel0.5.0:

Released 2016 June 10. Various updates to instrument properties, improved
documentation, and overhaul of internals in preparation for measured WFE data on
JWST SIs. 

JWST updates: 

 * New documentation on :ref:`jwst_instruments`
 * Updated all JWST SI pixel scales to latest measured values from ISIM CV3 and
   STScI Science Instruments Aperture File. 
 * Add coordinate inversion to get the correct (inverted) orientation of the OTE
   exit pupil relative to the ISIM focal plane. This will show up as an extra
   intermediate optical plane in all PSF calculations from this point, with the
   OTE pupil obscuration flipped upside down in orientation relative to the
   entrance pupil. 

   * As a consequence of this, many optical planes displayed will now look
     "upside down" relative to prior versions of WebbPSF. This affects all
     coronagraphic Lyot masks for instance, the NIRISS CLEARP and NRM pupils, etc.
     This is as intended, and reflects the actual orientation of those optics in the
     internal pupil planes relative to a detector image that has been oriented to have
     +V3 up and +V2 left (e.g. 'SCI' frame orientation on the sky, with north up and east left 
     if the position angle is zero).

 * Added software infrastructure for using measured instrument WFE from ISIM
   cryo-tests - however the data files are not yet ready and approved. This
   functionality will be fully activated in a near-future release (later this summer).
 * Added attributes for detector selection and pixel positions to all SIs, backed with 
   latest science instrument aperture file mapping between detector pixels and angular positions
   on the JWST focal plane.
 * Improved automatic toggling based on selected filter of instrument properties such as
   NIRCam short/long channel and pixel scales, and NIRISS and MIRI pupil masks. 
 * *Thanks to Kyle van Gorkom, Anand Sivaramakrishnan, John Stansberry, Colin Cox, 
   Randal Telfer, and George Hartig for assisting with information and data to
   support these updates.*

WFIRST updates:
 
 * Updated to `GSFC Cycle 6 modeling results
   <http://wfirst.gsfc.nasa.gov/science/Inst_Ref_Info_Cycle6.html>`_ for WFI.
 * Some behind-the-scenes refactoring to implementation details for field dependent
   WFE to support code sharing between the JWST and WFIRST classes.
 * *Thanks to Alden Jurling for assisting with information and clarifications on the Cycle 6 models.*


General:
 
 * New `Python PEP8 style guide <https://www.python.org/dev/peps/pep-0008/>`_ compliant names have been added
   for most function calls, e.g. ``calc_psf`` instead of ``calcPSF``, ``display_psf`` instead of 
   ``display_PSF`` and so forth. For now these are synonymous and both forms will work. The new styling is
   preferred and at some future point (but not soon!) the older syntax may be removed.

Version 0.4.1
=============

.. _rel0.4.1:

Released 2016 April 04. Mostly minor bug fixes, plus some updates to better match orientations of output files.

 * Fix an bug that ignored the rotation of the MIRI coronagraph occulters, introduced by changes in ``poppy`` 0.4.0; (`#91 <https://github.com/mperrin/webbpsf/issue/91>`__; @kvangorkom, @josephoenix, @mperrin)
   and also flip the sign of that rotation from 4.5 degrees counterclockwise to 4.5 clockwise, to match the actual hardware (`#90 <https://github.com/mperrin/webbpsf/issue/90>`__; @kvangorkom, @josephoenix, @mperrin)
 * Also flip orientations of some NIRCam coronagraphic masks and improve modeling of NIRCam coronagraph ND squares and occulter bar mounting hardware (`#85 <https://github.com/mperrin/webbpsf/issue/85>`__; @mperrin);
   and remove two obsolete filter data files that don't correspond to any actual filters in NIRCam.
 * Relocate ``measure_strehl`` function code into ``webbpsf`` (`#88 <https://github.com/mperrin/webbpsf/issue/88>`__; Kathryn St.Laurent, @josephoenix, @mperrin)
 * Other minor bug fixes and improved error catching 
   (`#87 <https://github.com/mperrin/webbpsf/issue/87>`__; @mperrin)
   (`#95 <https://github.com/mperrin/webbpsf/issue/95>`__; @mperrin)
   (`#98 <https://github.com/mperrin/webbpsf/pull/98>`__; @josephoenix)
   (`#99 <https://github.com/mperrin/webbpsf/issue/99>`__; @mperrin)
 * Better document how to make monochromatic PSFs (`#92
   <https://github.com/mperrin/webbpsf/issue/92>`__; @mperrin) and fix broken
   link in docs (`#96 <https://github.com/mperrin/webbpsf/pull/96>`__;
   @josephoenix).





Version 0.4.0
=============

.. _rel0.4.0:

Released 2015 November 20

* **WFIRST WFI support added**:

  * including all WFI filters and filter-dependent pupil masks.
  * including field dependence based on GSFC Cycle 5 modeling (`#75 <https://github.com/mperrin/webbpsf/pull/75>`__, @josephoenix)
  * including initial/prototype GUI interface based on Jupyter/IPython notebook widgets (`#79 <https://github.com/mperrin/webbpsf/pull/79>`__, @josephoenix)

* Updated filter transmission files for MIRI (based on Glasse et al. 2015 PASP) and NIRISS (based on flight filter measurement data provided by Loic Albert).
  (`#66 <https://github.com/mperrin/webbpsf/issues/66>`_, `#78 <https://github.com/mperrin/webbpsf/issues/78>`_; @mperrin)
* Added utility to check for appropriate version of the data files and request an update if necessary  (`#76 <https://github.com/mperrin/webbpsf/pull/76>`__, @josephoenix)
* Some documentation updates, including new documentation for the WFIRST functionality (@josephoenix, @mperrin)
* Bug fixes for minor issues involving OPD file units (`#74 <https://github.com/mperrin/webbpsf/pull/74>`__, @josephoenix), cleaner logging output, and some Python 3 compatibility issues.

.. note::

    When updating to version 0.4 you will need to also update your WebbPSF data files
    to the latest version as well.



.. _rel0.3.3:

Version 0.3.3
=================

Released July 1, 2015

* **Python 3 compatibility added.** All tests pass on Python 3.4. (`#2 <https://github.com/mperrin/webbpsf/issues/2>`_)
* Fixed an issue that would prevent users from adding defocus to PSF calculations
* WebbPSF no longer attempts to display a welcome message on new installs; that idea proved to be less helpful than originally expected.
* Added a ``CLEAR`` filter option for NIRISS, since the corresponding clear position is actually in the filter wheel rather than the pupil mask wheel. Rather than an actual filter, the profile for ``CLEAR`` is 1.0 between 0.6 microns and 5.0 microns per the stated limits of the detector, and 0.0 everywhere else. (`#64 <https://github.com/mperrin/webbpsf/issues/64>`_)
* Multi-wavelength calculations across a filter were not choosing a sensible number of wavelengths from the tables included in ``webbpsf-data``. (`#68 <https://github.com/mperrin/webbpsf/issues/68>`_)

.. _rel0.3.2:

Version 0.3.2
=================

Released February 23, 2015

This is a bug-fix release to address an issue that rendered the GUI unusable.
(See `#55 <https://github.com/mperrin/webbpsf/pull/55>`_.) API usage was unaffected.

(Ask not what happened to 0.3.1.)

.. _rel0.3.0:

Version 0.3.0
=================

Released 2015 February

This is a major release of WebbPSF, with several additions to the optical
models (particularly for slit and slitless spectroscopy), and extensive software
improvements and under-the-hood infrastructure code updates. Many
default settings can now be customized by a text configuration file in your home
directory.


**Updates to the optical models**:


 * Initial support for spectroscopy: *NIRSpec fixed slit and some MSA spectroscopy*, *MIRI
   LRS spectroscopy* (for both slit and slitless modes), and *NIRISS
   single-object slitless spectroscopy*.   To model one of these modes,
   select the desired image plane stop (if any) plus the pupil plane stop for the
   grating. WebbPSF does not yet include any model for the spectral dispersion
   of the prisms, so you will want to perform monochromatic calculations for
   the desired wavelengths, and coadd the results together yourself into a
   spectrum appropriately. For example::

    >> nirspec.image_mask = 'S200A1'
    >> nirspec.pupil_mask = 'NIRSpec grating'
    >> monopsf = nirspec.calcPSF(monochromatic=3e-6, fov_arcsec=3)

    >> miri.image_mask = 'LRS slit'
    >> miri.pupil_mask = 'LRS grating'
    >> miripsf = miri.calcPSF(monochromatic=10e-6)

    >> niriss.pupil_mask = 'GR700XD'
    >> monopsf = niriss.calcPSF(monochromatic=1.5e-6, oversample=4)


   In fact the NIRSpec class now automatically defaults to having the NIRSpec
   grating pupil stop as the selected pupil mask, since that's always in the beam. For
   MIRI you must explicitly select the 'LRS grating' pupil mask, and may select
   the 'LRS slit' image stop.  For NIRISS you must select the 'GR700XD' grating
   as the pupil mask, though of course there is no slit for this one.

   *Please note* This is new/experimental code and these models have not been validated
   in detail against instrument hardware performance yet. Use with appropriate caution, and
   we encourage users and members of the instrument teams to provide input on how this
   functionality can be further improved.
   Note also that MIRI MRS and NIRSpec IFU are still unsupported.

   Thanks to Loic Albert (U de Montreal) and Anand Sivaramakrishnan for data
   and many useful discussions on NIRISS SOSS.
   Thanks to Klaus Pontoppidan for proposing the NIRSpec and MIRI support and
   useful discussions. Thanks to Erin Elliott for researching the NIRSpec
   grating wheel pupil stop geometry, and Charles Lajoie for information on the
   MIRI LRS pupil stop.

 * Added NIRISS CLEARP pupil mask; this includes the obscuration from the pupil alignment reference.
   Given the pupil wheel layout, this unavoidably must be in the beam for any NIRISS
   long-wave PSFs, and WebbPSF will automatically configure it in the necessary cases. Thanks to Anand Sivaramakrishnan.

 * Minor bug fix to weak lens code for NIRCam, which previously had an incorrect scaling factor.
   Weak lens defocus values updated to the as-built rather than ideal values (which differ by 3%, but the as built values are very well calibrated).

 * Added defocus option to all instruments, which can be used to simulate
   either internal focus mechanism moves or telescope defocus during MIMF. For
   example, set ::

    >> nircam.options['defocus_waves']=3
    >> nircam.options['defocus_wavelength']=2.0e-6

   to simulate 3 waves of defocus at 2 microns, equivalently 6 microns phase delay peak-to-valley in the wavefront.

 * Added new option to offset intermediate pupils (e.g. coronagraphic Lyot
   stops, spectrograph prisms/grisms, etc) in rotation as well as in
   centering::

    >> niriss.options['pupil_rotation'] = 2  # degrees counterclockwise

 * Added support for rectangular subarray calculations. You can invoke these by
   setting fov_pixels or fov_arcsec with a 2-element iterable::

    >> nc = webbpsf.NIRCam()
    >> nc.calcPSF('F212N', fov_arcsec=[3,6])
    >> nc.calcPSF('F187N', fov_pixels=(300,100) )

   Those two elements give the desired field size as (Y,X) following the usual
   Python axis order convention. This is motivated in particular by the rectangular
   subarrays used in some spectroscopic modes.



**Other Software Updates & Enhancements**:


* Required Python modules updated, now with dependency on `astropy <http::/www.astropy.org>`_:

    * ``astropy.io.fits`` replaces ``pyfits`` for FITS I/O.
    * ``astropy.io.ascii`` replaces ``asciitable`` for ASCII table I/O.
    * ``atpy`` is no longer required.
    * New ``astropy.config`` configuration system is used for persistent
      settings.  This includes saving accumulated FFTW 'wisdom' so that future
      FFT-based calculations will begin more rapidly.
    * ``lxml`` now required for XML parsing of certain config files
    * ``psutil`` strongly recommended for cross-platform detection of
      available free RAM to enable better parallelization.

* Improved packaging infrastructure. Thanks to Christine Slocum, Erik Bray, Mark Sienkiewicz, Michael Droetboom,
  and the developers of the `Astropy affiliated package template <https://github.com/astropy/package-template>`_.
  Thanks in particular to Christine Slocum for integration into the STScI SSB software distribution.

* Improvements to parallelization code. Better :ref:`documentation for parallelization <performance_and_parallelization>`.  PyFFTW3 replaced with pyFFTW for optimized
  FFTs (yes, those are two entirely different packages).

* Alternate GUI using the wxpython widget toolkit in place of the older/less
  functional Tkinter tool kit. Thanks to Klaus Pontoppidan for useful advice in
  wxpython. This should offer better cross-platform support and improved long
  term extensibility. The existing Tkinter GUI remains in place as well.

    * The calculation options dialog box now has an option to toggle between monochromatic and broadband calculations. In monochromatic mode, the "# of wavelengths" field is
      replaced by a "wavelength in microns" field.
    * There is also an option to toggle the field of view size between units of arcseconds and pixels.
    * Log messages giving details of calculations are now displayed in a window as part of the GUI as well.
    * The wx gui supports rectangular fields of view. Simply enter 2 elements separated by a comma in the 'Field of view' text box. As a convenience, these
      are interpreted as (X,Y) sizes. (Note that this is opposite of the convention used in the programming interface noted above; this is potentially confusing but
      seems a reasonable compromise for users of the webbpsf GUI who do not care to think about Python conventions in axis ordering. Comments on this topic are welcome.)

* Improved configuration settings system. Many settings such as default
  oversampling, default field of view size, and output file format can now be
  set in a configuration file for persistence between sessions. So if you
  always want e.g. 8x oversampling, you can now make that the default. An
  example configuration file with default values will be created automatically the first
  time you run webbpsf now, including informative comments describing possible settings.
  This file will be in your astropy config directory, typically something like "~/.astropy/config".

    * New 'Preferences' dialog allows changing these persistent defaults through the GUI.

* New function webbpsf.setup_logging() adds some more user-friendliness to the
  underlying python logging system. This includes persistent log settings
  between sessions. See updated documentation in the :py:mod:`webbpsf` page.

* The first time it is invoked on a computer, WebbPSF will display a welcome
  message providing some information of use to new users. This includes checking
  whether the requisite data files have been installed properly, and alerting users
  to the location of the configuration file, among other things.

* Refactoring of instrument class and rebalancing where the lines between WebbPSF and POPPY had been blurry.

* Some bugfixes in the example code. Thanks to Diane Karakla, Anand Sivaramakrishnan, Schuyler Wolff.

* Various updates & enhancements to this documentation. More extensive documentation for POPPY now available as well. Doc theme derived from astropy.

* Improved unit test suite and test coverage. Integration with Travis CI for continuous testing: https://travis-ci.org/mperrin/webbpsf

* Updated to astropy package helpers framework 0.4.4


Version 0.2.8
=================

Released May 18, 2012

* Repaired functionality for saving intermediate opticals planes
* Coronagraph pupil shear shifts now use scipy.ndimage.shift instead of numpy.roll to avoid wrapping pixels around the edge of the array.
* Significant internal code reorganizations and cleanup:

        * switched package building to use `setuptools` instead of `distutils`/`stsci_distutils_hack`
        * `poppy` now installed as a separate package to more easily allow direct use.
        * new `Instrument` class in poppy provides much of the functionality previously in JWInstrument, to make it
          easier to model generic non-JWST instruments using this code.
        * Better packaging in general, with more attention to public/private API consistency
        * Built-in test suite available via `python setup.py test`

* Minor fix to MIRI ND filter transmission curve (Note: MIRI ND data is available on internal STScI data ditribution only)
* Binset now specified when integrating across bandpasses in pysynphoteliminating a previous warning message for that calculation.
* Stellar spectra are now by default drawn from the PHOENIX models catalog rather than the Castelli & Kurucz 2004 models. This is because the PHOENIX models have better spectral sampling at mid-infrared wavelengths.
* Default centroid box sizes are now consistent for measure_centroid() and the markcenter option to display_PSF(). (Thanks to Charles Lajoie for noting the discrepancy)
* TFI class (deprecated in version 0.2.6) now removed.

Version 0.2.7
=================

Released December 6, 2011

* Bug fix for installation problems in previous release 0.2.6 (thanks to Anand Sivaramakrishnan and Kevin Flaherty for bringing the problem to my attention).

* Updated FITS keywords for consistency with JWST Data Management System (DMS) based on DMS Software Design Review 1.

  * "PUPIL" keyword now is used for pupil mechanisms instead of OTE pupil intensity filename; the filename is available in "PUPILINT" now, for consistency with the OPD filename in "PUPILOPD" now.
  * "CORONMSK" instead of CORON
  * Some minor instrument-specific FITS keywords added via new _instrument_fits_header() functions for each instrument object.
  * For instance, NIRCam PSFs now have "MODULE" and "CHANNEL" keywords (eg. "MODULE = A", "CHANNEL = Short"). Note that there is no optical difference between modules A and B in this version of webbpsf.

* Added support for weak lenses in NIRCam. Note that the +4 lens is in the filter wheel and is coated with a narrowband interference filter similar to but wider than F212N.
  WebbPSF currently does not model this, and will let you simulate weak lens observations with any filter you want. As always, it's up to the user to determine whether
  a given webbpsf configuration corresponds to an actual physically realizable instrument mode.



Version 0.2.6
=================

Released November 7, 2011

* Updated & renamed TFI -> NIRISS.

  * Removed etalon code.
  * Added in filters transmissions copied from NIRCam
  * Removed coronagraphic Lyot pupils. Note: the coronagraphic occulting spots are machined into the pickoff mirror so will still fly, and thus are retained in the NIRISS model.
  * Slitless spectroscopy not yet supported; check back in a future version.
  * Fix to FITS header comments for NIRISS NRM mask file for correct provenance information.

  * TFI class still exists for back compatibility but will no longer be maintained, and may be removed in a future version of webbpsf.

* Strehl measurement code caches computed perfect PSFs for improved speed when measuring many files.
* Added GUI options for flat spectra in F_nu and F_lambda. (Thanks to Christopher Willmer at Steward Observatory for this suggestion)
* "display_psf" function renamed to "display_PSF" for consistency with all-uppercase use of PSF in all function names.
* numpy and pylab imports changed to 'np' and 'plt' for consistency with astropy guidelines (http://astropy.wikispaces.com/Astropy+Coding+Guidelines)
* poppy.py library updates (thanks to Anand Sivaramakrishnan for useful discussions leading to several of these improvements):

  * :py:class:`Rotation` angles can be specified in either degrees or radians. Added units parameters to Rotations.__init__
  * :py:class:`OpticalElement` objects created from FITS files use the filename as a default optic name instead of "unnamed optic".
  * :py:class:`FITSOpticalElement` class created, to separate FITS file reading functionality from the base OpticalElement class.
    This class also adds a 'pixelscale' keyword to directly specify the pixel scale for such a file, if not present in the FITS header.
  * Removed redundant 'pupil_scale' attribute: 'pixelscale' is now used for both image and pupil plane pixel scales.
  * unit test code updates & improvements.

* Miscellaneous minor documentation improvements.




Version 0.2.5
==============

Initial public release, June 1 2011. Questions, comments, criticism all welcome!

* Improved spectrum display
* Improved display of intermediate results during calculations.

Versions 0.2.1 - 0.2.3
=======================

* Smoother installation process (thanks to Anand Sivaramakrishan for initial testing)
* Semi-analytic coronagraphic algorithm added for TFI and NIRCam circular occulters (Soummer et al. 2007)
* Advanced settings dialog box added to GUI
* NIRCam pixel scale auto-switching will no longer override custom user pixelscales.
* slight fix to pupil file pixel scales to reflect JWST flat-to-flat diameter=6.559 m rather than just "6.5 m"
* Corrected NIRCam 430R occulter profile to exactly match flight design; other occulters still need to be tuned. Corrected all for use of amplitude rather than intensity profiles (thanks to John Krist for comparison models).
* added TFI NRM mode (thanks to Anand Sivaramakrishnan)


Version 0.2
============

Initial STScI internal release, spring 2011. Questions, comments, criticism all welcome!

* Much improved pysynphot support.
* Reworked calling conventions for calcPSF() routine source parameters.
* poppy.calcPSFmultiprocessor merged in to regular poppy.calcPSF
* Minor bug fixes to selection of which wavelengths to compute for more even sampling
* Default OPDs are now the ones including SI WFE as well as OTE+ISIM.
* Improved fidelity for NIRCam coronagraphic occulter models including ND squares and substrate border.




Version 0.1
============

Development, fall 2010.

* Support for imaging mode in all SIs and FGS
* Support for coronagraphy with MIRI, NIRCam, and TFI. Further enhancements in fidelity to come later.  Coronagraphic calculations are done using the direct FFT method, not Soummer's semi-analytic method (though that may be implemented in the future?).
* Up-to-date science frame axes convention, including detector rotations for MIRI and NIRSpec.
* Tunable wavelengths and appropriate bandwidths for TFI.
* Partial support for modeling IFU PSFs through use of the 'monochromatic' parameter.
* Revision V OPD files for OTE and SIs. Produced by Ball Aerospace for Mission CDR, provided by Mark Clampin.




