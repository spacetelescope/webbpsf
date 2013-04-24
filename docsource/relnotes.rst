.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Release Notes
######################


Known Issues
--------------
* You may see various warning messages while running computations, like thus::

    No handlers could be found for logger "webbpsf"

    ((<pysynphot.spectrum.Box object at 0x1047132d0> * nircam,im,f200w)) does not have a defined 
    binset in the wavecat table. The waveset of the spectrum will be used instead.

    Warning: invalid value encountered in absolute

  These can safely be ignored. 

* Calculations at large radii (> 500 lambda/D ~ 30 arcsec for 2 microns) will show numerical artifacts from Fourier aliasing and the implicit repetition of 
  the pupil entrance aperture in the discrete Fourier transform. If you need accurate PSF information at such large radii, please contact Marshall Perrin for
  higher resolution pupil data. 


**The following factors are NOT included in these simulations:**

* PSF variations across the field of view of any instrument (though each one has its own distinct OPDs for the center of its FOV).
* Optical distortions.
* Instrumental wavefront errors are not modeled separately, though they are included in some of the supplied RevV OPDs. 
* Coronagraphic masks are assumed to be perfect (i.e. the masks exactly match their design parameters.)
* No edge effects near the center of the FQPMs. (However, these are believed to be negligible in practice based on detailed simulations by Remi Soummer.)
* Any and all detector effects, including intrapixel sensitivity variations. There is no plan to include these at any point in WebbPSF itself.  Generate a subsampled PSF and use a separate detector model code instead. 

Plans for Future Releases
--------------------------
* Full support for the NIRSpec and MIRI IFUs may be added in a future release. Likewise for grisms.
* Realistic (but time consuming) jitter models (This code now available in beta form if you need it; contact Marshall.)
* Integration with OPD generation software and detector noise models.
* Possibly: separate handling of pre- and post- coronagraphic WFE in instruments, if this appears likely to be significant. 
* Python 3 support will be added as soon as it is needed, but is not an immediate priority. Any users who would like to run webbpsf under python 3, please let me know.


Version History and Change Log
-------------------------------


Version 0.3.0
=================

Released ?????


**Changes and Updates to the optical models**:


 * Bug fix to weak lens code for NIRCam, which previously had an incorrect scaling factor.  
 * Added defocus option to all instruments, which can be used to simulate either internal focus mechanism moves or telescope defocus during MIMF. For example, set ::
 
    >> nircam.options['defocus_waves']=3
    >> nircam.options['defocus_wavelength']=2.0e-6
    
   to simulate 3 waves of defocus at 2 microns, equivalently 6 microns phase delay peak-to-valley in the wavefront.

 * Added new option to offset intermediate pupils (e.g. coronagraphic Lyot stops, spectrograph prisms/grisms, etc) in rotation as well as in centering::

    >> niriss.options['pupil_rotation'] = 2  # degrees counterclockwise  

 * Added support for rectangular subarray calculations. You can invoke these by setting fov_pixels or fov_arcsec with a 2-element iterable::

    >> nc = webbpsf.NIRCam()
    >> nc.calcPSF('F212N', fov_arcsec=[3,6])
    >> nc.calcPSF('F187N', fov_pixels=(300,100) )

   Those two elements give the desired field size as (Y,X) following the usual Python axis order convention.


.. comment
  * Added a new model for NIRISS single-object slitless spectroscopy (SOSS).  Wide field slitless 


**Other Software Updates & Enhancements**: 


* Required Python modules updated, now with dependency on `astropy <http::/www.astropy.org>`_:

    * ``astropy.io.fits`` replaces ``pyfits`` for FITS I/O. 
    * ``astropy.io.ascii`` replaces ``asciitable`` for ASCII table I/O.
    * ``atpy`` is no longer required.
    * New ``astropy.config`` configuration system is used for persistent settings.


* New GUI using the wxpython widget toolkit in place of the older/less function Tkinter tool kit. Thanks to Klaus Pontoppidan for useful advice in wxpython. This should offer 
  better cross-platform support and improved long term extensibility. (For now, the existing Tkinter GUI remains in place but is deprecated and further development is not planned.) 

    * The advanced options dialog box now has an option to toggle between monochromatic and broadband calculations. In monochromatic mode, the "# of wavelengths" field is 
      replaced by a "wavelength in microns" field. 
    * There is also an option to toggle the field of view size between arcseconds and pixels. 
    * Log messages giving details of calculations are now displayed in a window as part of the GUI as well. 
    * The wx gui supports rectangular fields of view. Simply enter 2 elements separated by a comma in the 'Field of view' text box. As a convenience, these 
      are interpreted as (X,Y) sizes. (Note that this is opposite of the convention used in the programming interface noted above; this is potentially confusing but 
      seems a reasonable compromise for users of the webbpsf GUI who do not care to think about Python conventions in axis ordering.)


* New function webbpsf.setup_logging() adds some more user-friendliness to the
  underlying python logging system. This includes persistent log settings
  between sessions. See updated documentation in the :py:mod:`webbpsf` page. 

* Many settings such as default oversampling, default field of view size, and 
  output file format can now be set in a configuration file for persistence
  between sessions. So if you always want e.g. 8x oversampling, you can now
  make that the default.

* Some bugfixes in the example code. Thanks to Diane Karakla, Anand Sivaramakrishnan.

* Various minor updates & enhancements to this documentation.


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




