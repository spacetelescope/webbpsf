#############
Release Notes
#############

.. _known_issues:

Known Issues
--------------

See https://github.com/spacetelescope/webbpsf/issues for currently open issues and enhancement suggestions.

* Calculations at large radii (> 500 lambda/D ~ 30 arcsec for 2 microns) will
  show numerical artifacts from Fourier aliasing and the implicit repetition of
  the pupil entrance aperture in the discrete Fourier transform. If you need
  accurate PSF information at such large radii, please contact Marshall Perrin
  or Marcio Melendez for higher resolution pupil data.

**The following factors are NOT included in these simulations:**

* Coronagraphic masks are assumed to be perfect (i.e. the masks exactly match their design parameters.)
* Most detector effects, such as intrapixel sensitivity variations or interpixel capacitance. There are currently no plans to include these WebbPSF itself.  Generate a subsampled PSF and use a separate detector model code instead. The one exception is a scattering artifact in the MIRI imager detector substrate.

Road Map for Future Releases
--------------------------------
* Continued validation and updates as needed based on further analyses of instrument and telescope hardware test data.
* Support for the NIRSpec and MIRI IFUs may be added in a future release; level of detail is still TBD.
* Improved models for telescope WFE evolution over time.
* Possibly: separate handling of pre- and post- coronagraphic WFE in instruments, or pre- and post- NIRSpec MSA plane WFE; pending receipt of test data and models from the instrument teams.

.. _whatsnew:

Version History and Change Log
-------------------------------


Version 0.8.0
=============

*2018 Dec 15*

This release focused on software engineering improvements, rather than changes in any of the optical models or reference data. (In particular, there are NO changes in the reference data files; the contents of the WebbPSF version 0.8 data zip file are identical to the reference data as distributed for version 0.7.  This version of WebbPSF will work with either of those interchangably.).

.. admonition:: Python version support: Python 3 required

        This version drops support for Python 2.7. The minimum supported version of Python is now 3.5.

**New functionality:**

- *Added new capability to create grids of fiducial, distorted PSFs* spanning a chosen instrument/detector. This new ``psf_grid`` method is meant to be used as the first step of using the ``photutils`` package to do PSF-fitting photometry on simulated JWST PSFs. This method will output a list of or single ``photutils`` ``GriddedPSFModel`` object(s) which can then be read into ``photutils`` to apply interpolation to the grid and simulate a spatially dependent PSF anywhere on the instrument. See this `Jupyter notebook <https://github.com/spacetelescope/webbpsf/blob/master/notebooks/Gridded_PSF_Library.ipynb>`_ for examples. This method requires ``photutils`` version 0.6 or higher. [`#241, <https://github.com/spacetelescope/webbpsf/pull/241>` _, @shanosborne with inputs from @mperrin, @larrybradley, @hcferguson, and @eteq]

**Bug fixes and small changes:**

- *Improved the application of distortion to PSFs* to allow distorted PSFs to be created when the output mode is set to only “oversampled” or only “detector-sampled.”  When either of these modes is set in the options dictionary, the output will be an HDUList object with two extensions, where the 1st extension is the same PSF as in the 0th extension but with distortion applied. [`#229, <https://github.com/spacetelescope/webbpsf/pull/229>` _, @shanosborne]
- Also fixed distorted PSFs which were shifted off-center compared to their undistorted counterparts. These distorted PSFs had always been created in the correct detector location, but the values in the array returned by ``calc_psf`` were shifted off from the center. This bug was particularly apparent when the PSFs were set with a location near the edge of the detector. [`#219, <https://github.com/spacetelescope/webbpsf/pull/219>` _, @shanosborne]
- Fix FITS output from JWST OTE linear model, plus typo fixes and PEP8 improvements [#232, @laurenmarietta]
- Display code added for the PSF grid functionality mentioned above [#247, @mperrin]

**Software and Package Infrastructure Updates:**

- Removed Python 2.7 compatibility code, use of six and 2to3 packages, and Python 2 test cases on Travis (#236, #239, @mperrin, @kjbrooks]
- Packaging re-organized for consistency with current STScI package template (#240, @robelgeda)
- Documentation template updated for consistency with current STScI docs template (#250, @robelgeda)
- Documentation added or updated for a variety of features [#248, @mperrin]
- Various smaller code cleanup and doc improvements, including code cleanup for better Python PEP8 style guide compliance [#227, #255, @shanosborne]
- Updated to newer syntax for specifying pupil shifts of optical elements [#257, @mperrin]
- Unit tests added for defocused instruments, including the NIRCam weak lenses [#256, @mperrin]
- Updated astropy-helpers submodule to 3.0.2 [#249, @mperrin]
- Software development repo on Github shifted to within the `spacetelescope organization <https://github.com/spacetelescope/poppy>`_.


--------




Version 0.7.0
=============

*2018 May 30*


Note, when upgrading to this version you will need to update to the latest data files as well. This is
handled automatically if you use `conda`, otherwise you will need to download and install the data from:
`webbpsf-data-0.7.0.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-data-0.7.0.tar.gz>`_.

.. admonition:: Python version support: Future releases will require Python 3.

    Please note, this is the *final* release of WebbPSF to support Python 2.7. All
    future releases will require Python 3.5+. See `here <https://python3statement.org>`_ for more information on migrating to Python 3.

.. admonition:: Deprecated function names will go away in next release.

    This is also the *final* release of WebbPSF to support the older, deprecated
    function names with mixed case that are not compatible with the Python PEP8
    style guide (e.g. ``calcPSF`` instead of ``calc_psf``, etc). Future versions will
    require the use of the newer syntax.


**General:**

- Improved numerical performance in calculations  using new accelerated
  math functions in ``poppy``. It is highly recommended that users install the
  ``numexpr`` package, which enables significant speed boosts in typical
  propagations. ``numexpr`` is easily installable via Anaconda. Some use cases,
  particularly for coronagraphy or slit spectroscopy, can also benefit from GPU
  acceleration. See the latest ``poppy`` release notes for more.

**JWST optical model improvements:**


- *Models of field-dependent wavefront error are now included for all the SIs.*
  The OPD information is derived from the ISIM CV3 test campaign at Goddard, as
  described extensively in David Aronstein et al. "Science Instrument Wavefront
  Error and Focus: Results Summary from the ISIM Cryogenic Vacuum Tests:",
  JWST-RPT-032131. (See also `the SPIE paper version
  <http://adsabs.harvard.edu/abs/2016SPIE.9904E..09A>`_.) The measured SI
  wavefront errors are small, some tens of nanometers, and are in general less
  than the telescope WFE at given location. This information on SI WFE is
  provided to help inform modeling for what potential variations in PSFs
  across the field of view might look like, in broad trends. However it should
  _not_ be taken as precise guarantee of the exact amplitudes or functional form of
  those variations. The WFE was measured at a small handful of particular field
  points during CV3, and the resulting Zernike coefficients are interpolated to
  produce _estimated_ wavefront maps at all other field points across the focal
  planes.  Density and precision of the available measurements vary
  substantially between instruments.  [@mperrin, with contributions from
  @josephoenix in prior releases, and from @robelgeda and @JarronL for the
  interpolation between field points. [`#121
  <https://github.com/mperrin/webbpsf/pull/121>`_, `#187
  <https://github.com/mperrin/webbpsf/pull/187>`_]
- *Added new capabilities for modeling distortions of the image planes*, which
  cause slight deflections in the angles of diffractive features.  The result
  of geometric distortion is that detector pixels are not ideal square sections
  of the sky; they're slightly skewed parallelograms.  (See `the ACS handbook
  <http://www.stsci.edu/hst/acs/documents/handbooks/current/c05_imaging7.html#357374>`_
  for examples of what this looks like for Hubble PSFs) For the JWST
  instruments, this effect is largest for FGS, and fairly small but noticeable
  for the other SIs. See `this Jupyter notebook <https://github.com/mperrin/webbpsf/blob/master/notebooks/Distortion_examples.ipynb>`_ for 
  examples of the effect on JWST PSFs. Note that the distorted PSFs are added as *additional extensions*
  in the output FITS file, so you will need to read from extension 2 or 3 if you want the
  PSF with the distortion included; extensions 0 and 1 remain consistent with prior versions.  The distortion information is taken from the Science
  Instrument Aperture file (SIAF) reference data maintained at STScI. As a
  result the ``pysiaf`` package is a new dependency required for using
  ``webbpsf``.  The distortion calculations can add 1-3 seconds to each PSF calculation, and double the size of the output FITS files;
  if modeling distortion is not needed for your use case, you can deactivate this by setting ``add_distortion=False`` in calls to ``calc_psf``.  [ `#209 <https://github.com/mperrin/webbpsf/pull/209>`_,
  @shanosborne]
- *Added small nonzero pupil shears* for most instruments, based on measurements
  from the ISIM CV3 and OTIS cryo tests, adjusted for gravity release to produce
  predicted on-orbit pupil shears. See JWST-RPT-028027 and JWST-RPT-037134. For most
  imaging mode PSFs, this has _no_ practical effect because the SI internal pupils are
  oversized to provide tolerance, and the measured shears are well below that amount. 
  It has a small but nonzero effect for long-wave NIRISS filters with the CLEARP pupil
  obscuration.  The greatest effect is for MIRI coronagraphy since MIRI's Lyot stops were
  not undersized to allow for pupil shear, but even so the impact is small for the < 1% 
  expected shift.  Note that for NIRCam, the expected pupil shear is set to precisely 
  zero, given the expectation that NIRCam's steerable pickoff mirror will be used in flight 
  to achieve precise pupil alignment. 
  [`#212, <https://github.com/mperrin/webbpsf/pull/212>`_, @shanosborne, with inputs from
  Melendez, Telfer, and Hartig]
- *For MIRI only*, added new capability for modeling blurring due to
  *scattering of light within the MIRI imager detector substrate itself*. This
  acts as a cross-shaped convolution kernel, strongest at the shortest
  wavelengths. See MIRI document MIRI-TN-00076-ATC for details on the relevant
  physics and detector calibration.   This is implemented as part of the distortion framework, though
  it is different physics. See `this Jupyter notebook <https://github.com/mperrin/webbpsf/blob/master/notebooks/Distortion_examples.ipynb>`_ for
  example output. For F560W through F1000W this is a much more obvious effect than the subtle distortions. [`#209,
  <https://github.com/mperrin/webbpsf/pull/209>`_, @shanosborne]
- *Added new capabilities for modeling mirror moves of the JWST primary
  segments and secondary mirror*, using a linear optical model to adjust OPDs.
  Added a new `notebook demonstrating these capabilities
  <https://github.com/mperrin/webbpsf/blob/master/notebooks/Simulated%20OTE%20Mirror%20Move%20Demo.ipynb>`_.
  Note this code allows simulation of arbitrary mirror motions within a
  simplified linear range, and relies on user judgement what those mirror
  motions should be; it is not a detailed rigorous optomechanical model of the
  observatory.  [Code by @mperrin, with some fixes by Geda in <`#185
  <https://github.com/mperrin/webbpsf/pull/185>`_] 
- All the instrument+filter relative spectral response functions have been
  updated to values derived from the official validated JWST ETC reference
  data, using the Pandeia ETC release version 1.2.2. [@mperrin]


**WFIRST optical model improvements:**

- *The WFI optical model has been updated to use optical data from the Cycle 7
  design revision for WFI*. This includes a change in the instrument field of
  view layout relative to the axes, as shown `here
  <https://github.com/mperrin/webbpsf/pull/184>`_. [`#184
  <https://github.com/mperrin/webbpsf/pull/184>`_, @robelgeda]
- Added R062 filter. 
- Updated ``pupil_mask`` attribute for toggling between the masked and
  non-masked pupils now works the same way as that attribute does for the JWST
  instrument classes. Note, most users will not need to deal with this manually
  as the WFI class will by default automatically select the correct pupil based
  on the selected filter. [`#203
  <https://github.com/mperrin/webbpsf/issue/203>`_, @robelgeda]


**Bug fixes and minor changes:**

- All JWST instruments: Added new feature for importing OPD files produced with the JWST Wavefront Analysis System software [`#208 <https://github.com/mperrin/webbpsf/pull/208>`_, @skyhawk172] 
- All JWST instruments: Fix to generalize OPD loading code to handle either compressed or uncompressed OPDs [`#173 <https://github.com/mperrin/webbpsf/pull/173>`_, @JarronL]
- All JWST instruments: Fix to properly load the default number of wavelengths per calculation from the filters.tsv file, rather than defaulting to 10 wavelengths regardless. [@shanosborne])
- All JWST instrument: Fix to more correctly handle non-integer-pixel positions of the PSF when writing DET_X and DET_Y header keywords (`#205 <https://github.com/mperrin/webbpsf/pull/205>`_, @shanosborne]
- NIRCam and MIRI coronagraphy: Automatically set the detector coordinates and SI WFE maps based on the location of a selected coronagraph occulter. [`#181 <https://github.com/mperrin/webbpsf/pull/181>`_, @mperrin]
- NIRCam coronagraphy: Fix a sign error in offsets for the NIRCam coronagraph SWB occulters [`#172 <https://github.com/mperrin/webbpsf/issue/172>`_, @mperrin]. 
- NIRCam coronagraphy: Fix a half-percent throughput error in the round occulter masks [`#206  <https://github.com/mperrin/webbpsf/issue/206>`_, @mperrin]
- NIRCam coronagraphy: Fix an issue with transmission of the coronagraph bars precisely along the y axis, due to a typo [`#190  <https://github.com/mperrin/webbpsf/issue/190>`_, @JarronL]
- NIRCam coronagraphy: New option for shifting the coronagraph masks relative to the source, rather than vice versa. This is mostly of use for edge cases such as PSF library generation for the ETC, and is probably not of widespread utility. [`#191 <https://github.com/mperrin/webbpsf/issue/191>`_, @mperrin]
- NIRISS: Fix the `pupil_rotation` option so it works for NIRISS too, in particular for NRM/AMI. [`#118  <https://github.com/mperrin/webbpsf/issue/118>`_, @mperrin]
- NIRSpec: Very incomplete initial rudimentary support for the NIRSpec IFU, specifically just implementing the field stop for the IFU aperture. [@mperrin]
- Updated to newer version of the astropy_helpers package infrastructure [@sosey]
- Various smaller code cleanup and doc improvements, including code cleanup for better Python PEP8 style guide compliance [@mperrin, @shanosborne, @robelgeda, @douglase]
- The ``utils.system_diagnostic`` function now checks and reports on a few more things that might be useful in diagnosing performance issues. 


--------



.. _rel0.6.0:

Version 0.6.0
=============

*2017 August 11*

**JWST optical models:**

- Substantial update to the optical models for the telescope, to incorporate
  measurements of the as-built optics plus the latest expectations for
  alignments in flight.  The reference data layout has changed: each instrument
  now includes only two OPD files, a ``predicted`` and a ``requirements`` OPD.
  Ex: ``OPD_RevW_ote_for_NIRCam_predicted.fits.gz``. The OPD files are now
  derived from measured flight mirror surfaces (for high spatial frequencies),
  plus statistical models for their alignment in flight following wavefront
  sensing and control (for mid and lower spatial frequencies), as described in
  :doc:`jwst`.  Each OPD file still contains 10 different realizations of the
  statistical part.
- The NIRISS ``auto_pupil`` feature now recognizes that the ``CLEAR`` filter is used with the ``GR700XD`` pupil mask  [#151]
- Correctly convert wavelengths to microns when computing NIRISS ZnS index of refraction [#149]
- Aperture definitions now come from a copy of the SIAF bundled in ``jwxml`` rather than in the WebbPSF reference data.
- An alpha version of a linear optical model for adjusting OPDs is now provided for power-users, but currently unsupported and not documented.

**WFIRST optical models:**

- Addition of a model for the WFIRST CGI (Coronagraph Instrument) shaped pupil coronagraph by @neilzim [#154]

**General:**

- Jitter is now enabled by default (approximated by convolution with 0.007 arcsec FWHM Gaussian)
- Source offsets can now be specified as ``source_offset_x`` and ``source_offset_y`` in ``instrument.options`` (in addition to the existing ``instrument.options[‘source_offset_r’]`` and ``instrument.options[‘source_offset_theta’]``)
- The Astropy Helpers have been updated to v2.0.1 to fix various install-time issues.

.. _rel0.5.1:

Version 0.5.1
=============

Released 2016 November 2. Bug fix release to solve some issues that manifested
for AstroConda users.

 - Fixed a few missed version number->0.5.0 edits in install docs
 - Updated install instructions for Ureka->Astroconda change
 - Clarified release instructions for data packages
 - Fixed ConfigParser import in setup.py
 - Documented PSF normalization options better. (#112)
 - Updated Travis-CI config, consistent with poppy#187
 - Made a display tweak for the primary V2V3 annotation
 - Removed redundant ``calcPSF`` in favor of just using the superclass ``calc_psf`` (#132)
 - Updated ``measure_strehl`` to turn off SI WFE for perfect PSF calcs
 - Enforced Python 3.0+ compliance on code with ``__future__`` imports
 - Used ``six.string_types`` for Python 3.x compliance
 - Add version specs to dependencies in ``setup.py``
 - Made ``jwxml`` a dependency in ``setup.py``

.. _rel0.5.0:

Version 0.5.0
=============

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

.. _rel0.4.1:

Version 0.4.1
=============

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

.. _rel0.4.0:

Version 0.4.0
=============

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




