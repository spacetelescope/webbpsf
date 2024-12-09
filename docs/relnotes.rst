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
* Coronagraphic masks are assumed to be perfect (i.e. the masks exactly match their design parameters.)

------------------

.. _whatsnew:

Version History and Change Log
-------------------------------


Version 1.5.0
=============
*2024 December*
WebbPSF IS BEING MIGRATED TO A NEW REPOSITORY: STPSF (Space Telescope PSF)
To reflect its broader support for Roman as well as James Webb, WebbPSF is being migrated to a new repository: STPSF (Space Telescope PSF).
This transition is being done in such a way as to maintain back-compatibility for existing code, and existing installations will continue to run as-is.
This release will be among the final releases to WebbPSF and the initial release to STPSF.

This release requires new data webbpsf-data-1.5.0 https://stsci.box.com/s/oyoih9ibuaowksxzaia02ua4t6g1bbp1
For more information on this new data please view https://github.com/spacetelescope/webbpsf/pull/937

**What's Changed**
* Add IFU+datacubes docs page. by @mperrin in https://github.com/spacetelescope/webbpsf/pull/917

* Update for Photutils 2.0 by @mperrin in https://github.com/spacetelescope/webbpsf/pull/925

* MNT: Use hash for Action workflow versions and update if needed by @pllim in https://github.com/spacetelescope/webbpsf/pull/916

* Dependabot updates PLUS pyproject upper limit for photutils.  by @dependabot in https://github.com/spacetelescope/webbpsf/pull/926

* update ifu docs to describe about coord_system and pixel scale by @mperrin in https://github.com/spacetelescope/webbpsf/pull/923

* update photutils upper bound by @BradleySappington in https://github.com/spacetelescope/webbpsf/pull/935

* Add automatic data download by @WilliamJamieson in https://github.com/spacetelescope/webbpsf/pull/932

* switch a test case float comparison to np.allclose for robustness by @mperrin in https://github.com/spacetelescope/webbpsf/pull/934

* Add support for new sensing point and target phase map by @obi-wan76 in https://github.com/spacetelescope/webbpsf/pull/937

* Minor update to docs about the WFI optical model by @Skyhawk172 in https://github.com/spacetelescope/webbpsf/pull/936

* update monthly trending plot to allow shifted date ranges by @mperrin in https://github.com/spacetelescope/webbpsf/pull/943

* Readme advise by @BradleySappington in https://github.com/spacetelescope/webbpsf/pull/944

* update MAST query for OPDs to also retrieve OPDs from FP6 by @mperrin in https://github.com/spacetelescope/webbpsf/pull/945

## New Contributors
* @WilliamJamieson made their first contribution in https://github.com/spacetelescope/webbpsf/pull/932

**Full Changelog**: https://github.com/spacetelescope/webbpsf/compare/v1.4.0...1.5.0


Version 1.4.0
=============

*2024 September*

We are pleased to announce the release of the latest version of WebbPSF version 1.4.0, now available on PyPi and GitHub. This release comes with new features and improvements including but not limited to:

**What's Changed**
* Improve trending plots, particularly wfe_histogram arrows display by @mperrin in https://github.com/spacetelescope/webbpsf/pull/870

* build(deps): bump photutils from 1.12.0 to 1.13.0 by @dependabot in https://github.com/spacetelescope/webbpsf/pull/872

* more strict units handling; fixes some issues for astropy 6.0.0 compatibility by @mperrin in https://github.com/spacetelescope/webbpsf/pull/879

* add dedicated cache workflow that can be called by other projects by @zacharyburnett in https://github.com/spacetelescope/webbpsf/pull/877

* fix syntax error in `.github/workflows/download_data.yml` by @zacharyburnett in https://github.com/spacetelescope/webbpsf/pull/882

* use `pathlib` to handle a `WEBBPSF_PATH` with double slashes by @zacharyburnett in https://github.com/spacetelescope/webbpsf/pull/886

* Fix #888, typo in NIRSpec disperser list by @mperrin in https://github.com/spacetelescope/webbpsf/pull/889

* Use MASK_NRM to match AMI data by @rcooper295 in https://github.com/spacetelescope/webbpsf/pull/893

* minor: trending function add cycle 3 routine WFS PIDs by @mperrin in https://github.com/spacetelescope/webbpsf/pull/878

* Update NIRSpec aperture PAs by @mperrin in https://github.com/spacetelescope/webbpsf/pull/897

* build(deps): bump scipy from 1.13.0 to 1.14.0 by @dependabot in https://github.com/spacetelescope/webbpsf/pull/874

* Update citations by @BradleySappington in https://github.com/spacetelescope/webbpsf/pull/898

* using z range in requirements.txt by @BradleySappington in https://github.com/spacetelescope/webbpsf/pull/909

* Docs improvements, including for matching PSFs to data by @mperrin in https://github.com/spacetelescope/webbpsf/pull/899

* Enhance capabilities to simulate PSFs with larger JWST pupil  by @obi-wan76 in https://github.com/spacetelescope/webbpsf/pull/908

* IFU mode improvements, continued by @mperrin in https://github.com/spacetelescope/webbpsf/pull/890

* Roman docs & figures updates by @Skyhawk172 in https://github.com/spacetelescope/webbpsf/pull/910

* add CI tests for mast_wss, in test_mast_wss.py by @mperrin in https://github.com/spacetelescope/webbpsf/pull/911


**Full Changelog**: https://github.com/spacetelescope/webbpsf/compare/v1.3.0...v1.4.0

Version 1.3.0
=============

*2024 May*

This release comes with new features and improvements for JWST PSF models, including but not limited to:

1. Improved support for NIRSpec and MIRI IFU PSF calculations, including addition of a ``mode`` attribute for toggling between imaging mode and IFU mode simulations; an option for much faster (but slightly simplified) IFU datacube PSF calculations. Spectral bandpass information was added for the IFU bands for both NIRSpec and MIRI. For NIRSpec, IFU mode PSF outputs are rotated by an additional 90 degrees to match the convention used in pipeline-output s3d datacubes made using the IFUAlign orientation. This extra rotation can be optionally disabled if so desired; see the NIRSpec class docstring.
2. Improved PSF models for MIRI imager, in particular including an empirical model for the field-dependent shifts of the cruciform artifact seen at short wavelengths, and improved fidelity for modeling MIRI LRS slit PSFs.
3. For all JWST instruments, simulation of interpixel capacitance effects is now included for oversampled outputs as well as for the detector-sampled outputs.
4. Various improvements to the `setup_sim_to_match_data` function to automatically configure PSF simulations correctly for a more complete set of JWST observing modes.
5. Additional trending functions for assessing JWST wavefront error stability over time around specified science observations, and improvements to existing trending plots.

Please note, the minimum supported version of Python is increaed to Python 3.10, consistent with the minimumfor astropy.

**What's Changed**

* Add interpixel capacitance effects (IPC) for both distortion extension in the simulated PSF and adding some more per-instrument specializations in setup_sim_to_match_file by @obi-wan76 in https://github.com/spacetelescope/webbpsf/pull#768

* setup_sim_to_match_file fix setting aperture name for SW/LW parallel coronagraphy  by @mperrin in https://github.com/spacetelescope/webbpsf/pull/752

* Trending: add option to plot OTE-only WFE in wfe_histogram_plot by @Skyhawk172 in https://github.com/spacetelescope/webbpsf/pull/750

* #727_RMS_Label fixed the issue that rms_label is not initialized. by @bchen2 in https://github.com/spacetelescope/webbpsf/pull/730

* use PyPI upload workflow from OpenAstronomy by @zacharyburnett in https://github.com/spacetelescope/webbpsf/pull/749

* import GriddedPSFModel from photutils.psf by @braingram in https://github.com/spacetelescope/webbpsf/pull/755

* update NRC coron dispersion coeffs by @JarronL in https://github.com/spacetelescope/webbpsf/pull/766

* Update trending plot to subtract piston before displaying proposed correction by @mperrin in https://github.com/spacetelescope/webbpsf/pull/756

* add aperture name support for NIRSpec slit and IFU apertures; also, 100x faster data cube calculations. by @mperrin in https://github.com/spacetelescope/webbpsf/pull/767

* Support MIRI LRS slit in aperturename and setup_sim_to_match_file functions by @mperrin in https://github.com/spacetelescope/webbpsf/pull/781

* Coronagraph calcs: add 'coron_include_pre_lyot_plane' option for extra output plane by @mperrin in https://github.com/spacetelescope/webbpsf/pull/778

* Add `nrc_ta_image_comparison` function to trending tools by @mperrin in https://github.com/spacetelescope/webbpsf/pull/789

* Support the slightly non-square dimensions of the MIRI detector (1024, 1032) pixels in size. Fixes #676 by @mperrin in https://github.com/spacetelescope/webbpsf/pull/803

* filter custom opbtable for month by @BradleySappington in https://github.com/spacetelescope/webbpsf/pull/802

* Allow setup_sim_to_match_data to select choice of WFS before or after by @mperrin in https://github.com/spacetelescope/webbpsf/pull/819

* [SCSB-145] require Python 3.10 by @zacharyburnett in https://github.com/spacetelescope/webbpsf/pull/817

* add show_wfs_during_program function by @mperrin in https://github.com/spacetelescope/webbpsf/pull/798

* add delta_wfe_around_time function by @mperrin in https://github.com/spacetelescope/webbpsf/pull/826

* Performance enhancement: avoid repeated slow loads of the SIAF by @mperrin in https://github.com/spacetelescope/webbpsf/pull/825

* NRC TA plot enhancements by @mperrin in https://github.com/spacetelescope/webbpsf/pull/794

* add functions to download WFSC WL image data by @mperrin in https://github.com/spacetelescope/webbpsf/pull/827

* improve MIRI LRS model details by @mperrin in https://github.com/spacetelescope/webbpsf/pull/787

* #839 Visit Id naming smarts in get_visit_nrc_ta_image by @kulpster85 in https://github.com/spacetelescope/webbpsf/pull/840

* #835 - Improvements to trending.plot_wfs_obs_delta & wfe_histogram_pl… by @kulpster85 in https://github.com/spacetelescope/webbpsf/pull/836

* avoid linear interpolation between WFE values by @mperrin in https://github.com/spacetelescope/webbpsf/pull/834

* Update for WFI to use pysiaf instead of soc_roman_tools by @Skyhawk172 in https://github.com/spacetelescope/webbpsf/pull/848

* Infrastructure improvements for improved IFU sims by @mperrin in https://github.com/spacetelescope/webbpsf/pull/770

* Allow specifying NIRCam LW detectors equivalently like 'NRCA5' or 'NRCALONG' by @mperrin in https://github.com/spacetelescope/webbpsf/pull/849

* Improved model for MIRI cruciform artifact, part 1 by @mperrin in https://github.com/spacetelescope/webbpsf/pull/837

* Manual lint by @BradleySappington in https://github.com/spacetelescope/webbpsf/pull/847

* Implement support for NIRCam DHS sub apertures by @mperrin in https://github.com/spacetelescope/webbpsf/pull/845

**New Contributors**
* @braingram made their first contribution in https://github.com/spacetelescope/webbpsf/pull/740

* @bchen2 made their first contribution in https://github.com/spacetelescope/webbpsf/pull/731

* @eteq made their first contribution in https://github.com/spacetelescope/webbpsf/pull/807

* @bryce-wedig made their first contribution in https://github.com/spacetelescope/webbpsf/pull/815

**Full Changelog**: https://github.com/spacetelescope/webbpsf/compare/v1.2.1...v1.3.0.rc1


Version 1.2.1
=============
Minor documentation updates

Version 1.2.0
=============

*2023 August*

We are pleased to announce the release of the latest version of WebbPSF version 1.2.0, now available on PyPi and GitHub. This release comes with new features and improvements including but not limited to:

1. The addition of detector effects for JWST simulations. H2RG detector effects are included in two flavors, a simple ad hoc Gaussian convolution to model charge diffusion effects and a set of convolution kernels to model interpixel capacitance (IPC) and post-pixel coupling effects. We have found that these effects greatly improve the agreement between observations and simulations. See `JWST Detector Effects for more details. <https://webbpsf.readthedocs.io/en/latest/jwst_detector_effects.html>`_

2. A new utility function for simulating matching PSFs to science data. See `Matching PSF sims to in-flight JWST data <https://webbpsf.readthedocs.io/en/latest/jwst_matching_psfs_to_data.html>`_.

3. Implement geometric distortion for Roman using the Roman SIAF.

4. Various improvements for OTE trending.

**What's Changed**

* Fixed trending histogram binning so that bars add up to 1.0 by @Skyhawk172 in https://github.com/spacetelescope/webbpsf/pull/634

* Add phase retrieval crosscheck plot and wfs obs delta plot by @mperrin in https://github.com/spacetelescope/webbpsf/pull/650

* Add opdtable as positional param to monthly_trending_plot by @kulpster85 in https://github.com/spacetelescope/webbpsf/pull/600

* Update to read SI pixelscales directly from pysiaf by @mperrin in https://github.com/spacetelescope/webbpsf/pull/626

* Update/enhance trending plot to show WSS proposed corrections by @mperrin in https://github.com/spacetelescope/webbpsf/pull/642

* Add notebooks for plotting JWST SI WFE, and JWST SI MIMF field points by @mperrin in https://github.com/spacetelescope/webbpsf/pull/652

* Add H2RG detector effects sim framework by @mperrin and @obi-wan76 in https://github.com/spacetelescope/webbpsf/pull/671

* Tune detector effects model parameters to better match measured ePSFs by @mperrin in https://github.com/spacetelescope/webbpsf/pull/693

* Non-standard pixel sizes for distortion by @JarronL in https://github.com/spacetelescope/webbpsf/pull/669

* Add setup_sim_to_match_data function by @mperrin in https://github.com/spacetelescope/webbpsf/pull/706

* Add trending plot function "show_wfs_around_obs" by @mperrin in https://github.com/spacetelescope/webbpsf/pull/705

* Additional fixes to trending.py by @Skyhawk172 in https://github.com/spacetelescope/webbpsf/pull/688

* Implement distortion for Roman by @Skyhawk172 in https://github.com/spacetelescope/webbpsf/pull/668

**Full Changelog**: https://github.com/spacetelescope/webbpsf/compare/v1.1.1...v1.2.0

Note, this release requires updating your WebbPSF data files to version 1.2.0, `webbpsf-data-1.2.0.tar.gz <https://stsci.box.com/shared/static/34g3slaq4jidgccqj25qqo80tlk6tubl.gz>`_


Version 1.1.1
=============
*2022 December 14*

Minor bug fix release and improvements in JWST wavefront trending plots.

**James Webb Space Telescope improvements**:

 * Fix a units issue and filename inconsistency in one of the data files for NIRCam wavefront error at the wavefront sensing field point. (:issue:`612`, :pr:`613:` by :user:`mperrin`, :user:`obi-wan76`)
 * Improvements in OTE wavefront trending plots and  phase decomposition tools (:pr:`598` by :user:`kulpster85`, :pr:`599`, :pr:`601` by :user:`mperrin`, :pr:`603` by :user:`Skyhawk172:`,
   :pr:`621` by :user:`obi-wan76`)
 * Bug fixes for OTE field dependence flag (:pr:`595` by :user:`mperrin`)
 * Updates various package dependencies to upstream latest versions.


Version 1.1.0
=============
*2022 September 23*

*First release with JWST in flight optical performance!*  Updates and tools added after completion of commissioning.

Note, this release requires updating your WebbPSF data files to version 1.1.0. See :ref:`here <data_install>` .

This release's upgraded requirements drop support for Python 3.7, meaning conda installation is temporarily unavailable since the AstroConda channel is not equipped for newer Python versions. Installation with pip works as normal.

**James Webb Space Telescope OTE model improvements**:

 * Add feature to use measured OPDs from wavefront sensing in flight, including retrieval from MAST. See :doc:`jwst_measured_opds`. (:pr:`556`, :pr:`559`, :pr:`560`, :pr:`571` by :user:`mperrin; :pr:`563` by :user:`rcooper295`; :pr:`579` by :user:`obi-wan76`)
 * Add functions to trend and display wavefronts over time. See :doc:`jwst_measured_opds`.
 * Updated default line-of-sight jitter for JWST observations to 1 milliarcsecond instead of 6 (1 sigma per axis).
 * Updated default OPD to be an actual measured on-orbit OPD from early in cycle 1 science operations.

**Software and Package Infrastructure Updates:**

 * Add support for Python 3.10; drop support for Python 3.7 (:pr:`549` by :user:`shanosborne`)
 * Fixes to a few minor plotting bugs (:pr:`537` by :user:`shanosborne`; :pr:`581`, :pr:`582` by :user:`mperrin`)
 * Thanks to :user:`jsoref` for contributing :pr:`520` with spelling corrections, and :user:`NaincyKumariKnoldus` for fixing a bad link in the docs.
 * Add unit test for the coronagraph mask shift option (:pr:`578` by :user:`mperrin`)


Version 1.0.0
=============
*2021 December 10*

For JWST, this release includes updates to WebbPSF just prior to the launch. For Roman, it includes updates to use the Cycle 9 optical model results.

**James Webb Space Telescope OTE model improvements**:

* Updates in sign conventions for representing WFE, for strict consistency with the JWST WSS and other tools. Much of this was implemented by upstream changes in ``poppy``; see `this page in the poppy docs <https://poppy-optics.readthedocs.io/en/latest/sign_conventions_for_coordinates_and_phase.html>`_ for details.  (:pr:`397`, :pr:`419` by :user:`mperrin`, :pr:`418` by :user:`Skyhawk172`)
* Significant update to JWST OTE optical models, to reflect more recent 2020 optical modeling of the as-built observatory (the "PSR2020" integrated modeling cycle). These have noticeably lower WFE than the prior models (which were intentionally conservative, but ended up being more conservative than intended); typically the WFE is lower by some tens of nanometers in the new "prelaunch_predicted" OPDs. See details in :ref:`jwst_ote_details`. We will all learn together in 2022 how well these models predict the observatory's performance in flight. (:pr:`512`, :user:`mperrin`).
* Add models of OTE field dependence from the nominal OTE design and as-built optics (:pr:`389` by :user:`grbrady`, :pr:`505` by :user:`mperrin`) and from any misalignment of the secondary mirror, such as would be measured and corrected in MIMF (:pr:`392` by :user:`Skyhawk172`). These additions were enabled by more consistent use of JWST Linear Optical Model framework behind the scenes (:pr:`378` by :user:`mperrin`). This model of field dependence plus the updated OTE OPD files should yield a more comprehensive and precise model of PSF variations across the observatory.
* Add an option to use a lookup table of field dependent OPDs from Ball's ITM tool (for JWST team internal use in
  pre-launch wavefront team practices and rehearsals). (:pr:`425` by :user:`Skyhawk172`, :pr:`474` by :user:`mperrin`)
* Update the JWST OTE Linear Model to allow more flexible pupil sampling, in particular using higher sampling to reduce Fourier aliasing in certain FGS calculations (:pr:`440` by :user:`kjbrooks`)
* New capability for visualizing the JWST optical budget terms as represented in WebbPSF. See :doc:`jwst_optical_budgets`.


**James Webb Space Telescope instrument model improvements**:

* MIRI: Minor updates to pixel scale and rotation (:pr:`456` by :user:`mperrin`),
  an improved model of the MIRI imager detector cross artifact (:pr:`417` by :user:`mperrin`)
  and correctly label MIRI's P750L prism for the LRS mode as a prism, not a grating (:pr:`477` by :user:`mperrin` and :user:`skendrew`)
* MIRI: Add capability for shifting MIRI coronagraph masks, consistent with NIRCam sim capabilities (:pr:`428` by :user:`JarronL`)
* NIRCam: Higher fidelity model of NIRCam weak lenses, including field dependence, non-linear interactions between lenses,
  and as-built measured performances. (:pr:`496` by :user:`mperrin`, using results of calibration work by Randal Telfer)
* All SIs: Substantial performance improvements speeding up the calculation of optical distortion (:pr:`429`, :user:`jarronL`)

**Nancy Grace Roman Space Telescope and instrument model improvements**:

* Use of Cycle 9 optical and integrated modeling results, including updated Zernike coefficients, pupil images, and filter throughputs.
* Updated :py:obj:`~webbpsf.RomanInstrument` pointing stability to 12 milliarcseconds per axis, following new predictions [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* :py:obj:`WFI` wavelength range now covers 0.48 - 2.3 µm [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Added ``WFI``'s new F213 filter [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Renamed ``WFI``'s ``'P120'`` filter to ``'PRISM'`` [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Split ``WFI``'s ``'G150'`` filter into ``'GRISM0'`` and ``'GRISM1'`` components to represent the transmission for the grism's  undispersed zeroth order and dispersed first order, respectively [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Renamed WFI pupil masks to ``'SKINNY'`` (formerly ``'RIM_MASK'`` in version 0.9.2), ``'WIDE'`` (formerly ``'FULL_MASK'``), ``'GRISM'``, and ``'PRISM'`` (also formerly captured in ``'RIM_MASK'``) [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Created new :py:meth:`~webbpsf.WFI.lock_pupil()` and :py:meth:`~webbpsf.WFI.lock_pupil_mask()` methods for advanced users who prefer to disable automated selections and instead stick with a specific pupil file or mask, respectively. The corresponding ``WFI.unlock_pupil()`` and ``WFI.unlock_pupil_mask()`` methods return the class to its normal behavior [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Locked ``WFI.pupil`` and ``WFI.pupil_mask`` attributes from direct assignment given the new lock/unlock schema [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Renamed ``WFI.override_aberrations()`` to :py:meth:`~webbpsf.WFI.lock_aberrations()` and ``WFI.reset_override_aberrations()`` to :py:meth:`~webbpsf.WFI.unlock_aberrations()` to reinforce the new lock/unlock schema [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Condensed and refactored existing tests [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* New algorithm for field point nearest approximation/extrapolation [:pr:`466` by :user:`ojustino` with :user:`robelgeda`]
* Renamed ``CGI`` class to :py:obj:`RomanCoronagraph` [:pr:`516`, :pr:`517`, :user:`ojustino` with :user:`mperrin`]

**Software and Package Infrastructure Updates:**

* Software engineering improvements to meet STScI INS-JWST Software Standards (:pr:`404` by :user:`shanosborne`)
* Migrate optional dependency for synthetic photometry from pysynphot to synphot (:pr:`424` by :user:`shanosborne`)
* Deprecated the ``jwxml`` package, and moved the SUR (Segment Update Request) parsing code from that package into WebbPSF (:pr:`390` by :user:`shanosborne`)
* Various minor bug fixes (:pr:`410`, :pr:`422`, :pr:`427`, :pr:`497` by :user:`mperrin`, :pr:`423` by :user:`kjbrooks`, :pr:`493` by :user:`JarronL`)
* Updates to recommended (not minimum) dependency versions. Drop support for Python 3.6. (various PRs by :user:`shanosborne`)
* Remove deprecated older code including the GUIs (:pr:`439` by :user:`mperrin`)
* Streamline test suite to keep CI runtimes manageable (:pr:`459` by :user:`mperrin`)

------------------


Version 0.9.2
=============
*2021 July 23*

This release only improves a subset of WFIRST functionality; additional improvements to both WFIRST (including renaming to Roman) and JWST models will be at the upcoming 1.0.0 major release.

**WFIRST Improvements**

- New Grism and Prism filters: [:pr:`416`, :pr:`471`, :user:`robelgeda`]

    - `GRISM_FILTER = 'G150'`
    - `PRISM_FILTER = 'P120'`

- Changing filters to `G150` or  `P120` changes the mode of the WFI and the aberrations files (unless there is a user aberrations override) [:pr:`416`, :pr:`471`, :user:`robelgeda`]
- New `WFI.mode`: Class property that returns the current mode of the WFI instance by passing the current filter to `WFI. _get_filter_mode`. WFI modes are: [:pr:`416`, :pr:`471`, :user:`robelgeda`]

    -  Imaging
    -  Grism
    -  Prism

- New `WFI.override_aberrations(aberrations_path)`: Overrides and locks the current aberrations with aberrations at `aberrations_path`. Lock means changing the filter/mode has no effect on the aberrations. [:pr:`416`, :pr:`471`, :user:`robelgeda`]
- New `WFI.reset_override_aberrations()`: Releases `WFI.override_aberrations` lock and start using the default aberrations. [:pr:`416`, :pr:`471`, :user:`robelgeda`]
- New Tests for mode and filter switching. [:pr:`416`, :pr:`471`, :user:`robelgeda`]
- New Field point nearest point approximation (extrapolation). [:pr:`416`, :pr:`471`, :user:`robelgeda`]

**Software and Package Infrastructure Updates:**

- This release uses Github Actions CI and removes TravisCI. [:pr:`455`, :user:`shanosborne`, :pr:`471`, :user:`robelgeda`]

--------

Version 0.9.1
=============
*2020 June 22*

This minor release resolves several bugs and occasional installation issues and updates behind-the-scenes package infrastructure for consistency with current astropy and numpy releases. There are small improvements to a few aspects of JWST models as detailed below (in particular for wavelength dispersion in NIRCam LW coronagraphy and in tools for modeling time-dependent WFE) but the vast majority of JWST PSF calculations are not changed in any way.

There are no changes in reference data, so the WebbPSF reference data files for 0.9.0 should continue to be used with this release.

.. admonition:: Python version support: Python 3.6+ required

        This version drops support for Python 3.5. The minimum supported version of Python is now 3.6.


**JWST Improvements**

- *Apply wavelength dependent offsets for NIRCam coronagraphic PSFs* due to the dispersion from the optical wedge in the coronagraphic pupil masks. This primarily affects the LW channel with approximately 0.015 mm/um dispersion. The SW channel is almost a factor of 10 smaller and mostly negligible, but has been included for completeness. [:pr:`347`, :user:`JarronL`]
- *Improved models for OTE wavefront variations over time* by adding utility functions for decomposing WFE models into piston, tip, tilt motions in the JWST control coordinate system, adding a model for frill-induced WFE drift, adding a model for IEC-heater-induced WFE drift, and adding an option to adjust amplitude of OTE backplane thermal drift model for B.O.L. vs E.O.L. expected amplitudes. [:pr:`340`, :user:`mperrin`]
- *Add new* ``aperturename`` *attribute* for JWST instruments which returns the SIAF aperture name used for transforming between the detector position and instrument field of view on the sky. [:pr:`360`, :user:`mperrin`]. Relatedly, improves setting of detector geometry for NIRCam to automatically set the SIAF aperture name based on detector, filter, and coronagraph image mask and pupil mask settings. This can be turned off by setting ``auto_apname=False``. [:pr:`351`, :user:`JarronL`]
- Add model for image jitter with JWST in coarse point mode under two different assumptions about LOS stability. This is relevant only for commissioning simulations. [:pr:`345`, :pr:`346`, :user:`mperrin`]
- Documentation updates, in particular adding :ref:`figures of JWST instrument internal wavefront error models <jwst_instruments>`. [:pr:`369`, :user:`mperrin`]

**General bug fixes and small changes:**

- Allow FGS detector to be set to ``GUIDER1`` and ``GUIDER2``, while still supporting old method of setting the detector (using ``FGS1`` and ``FGS2``) [:pr:`361`, :user:`mperrin`]
- Add ``allow_huge=True`` option to ``astropy.convolution.convolve_fft`` call when applying MIRI distortion so it can handle large arrays when calculating PSFs in very large FOV by using a higher resolution pupil and OPD. [:pr:`354`, :user:`obi-wan76`]
- Fixed bug that caused an error when plotting OPDs using the ``display_opd`` function [:pr:`362`, :user:`shanosborne`]
- Update default NIRSpec detector coordinates to be the S1600A1 square aperture coordinates in imaging mode, rather than an implausible location outside of the MSA field of view. [:pr:`348`, :user:`mperrin`]
- Updated Simulated OTE Mirror Move Demo notebook. [:pr:`343`, :user:`kjbrooks`]
- Improved the reproducibility of the thermal slew model with small updates to the ``update_opd`` and ``move_jsc_acf`` functions. [:pr:`339`, :user:`mperrin`]

**Software and Package Infrastructure Updates:**

- *The minimum Python version is now 3.6.* [:pr:`353`, :user:`mperrin`]
- Removed dependency on ``astropy-helpers`` sub-package [:pr:`337`, :user:`shanosborne`]
- Fixed problem that resulted in the ``otelm/`` and ``tests/surs/`` sub-directories not installing correctly. [:pr:`356`, :user:`shanosborne`]
- Removed python 3.5 testing and added python 3.8 testing in Travis continuous integration. [:pr:`352`, :user:`mperrin`]
- Documentation added and/or updated for a variety of features, including referencing the newly renamed Nancy Grace Roman Space Telescope (formerly WFIRST). [:pr:`364`, :pr:`360`, :pr:`330`,  :user:`shanosborne, mperrin`]

--------




Version 0.9.0
=============
*2019 November 25*

Note, when upgrading to this version you will need to update to the latest data files as well. This is handled automatically if you use `conda`, otherwise you will need to download and install the data from: `webbpsf-data-0.9.0.tar.gz <https://stsci.box.com/shared/static/qcptcokkbx7fgi3c00w2732yezkxzb99.gz>`_.


**JWST Improvements**

- *Added a new capability to model the impact of thermal variations*, from telescope slews relative to the sun, onto mirror alignments and therefore onto PSFs. This new ``thermal_slew`` method  can be used to create a delta OPD for some elapsed time after the slew at either the maximum slew angle, some specified angle, or with a scaling factor applied to maximum case. Once combined with an input OPD (requirements or predicted), the new shape of the mirrors can be used to simulate predicted PSFs some time after a slew. See this `Jupyter notebook (ex1) <https://github.com/spacetelescope/webbpsf/blob/stable/notebooks/Example%20Construction%20of%20OPDs%20from%20Delta%20Time%20After%20Slew.ipynb>`_ for examples. [:pr:`269`, :user:`kjbrooks`]
- *Improved wavefront error extrapolation method for field points near FOV corners* that are outside the bounds of Zernike reference table data, in order to provide more seamless extrapolation.  [:pr:`283`, :user:`JarronL`]
- *Improvements in NIRCam optical model*: Updated polynomial model for NIRCam defocus versus wavelength. Adds Zernike coefficients for the wavefront error at NIRCam coronagraphy field points. [:pr:`283`, :user:`JarronL`]
- NIRISS NRM mask was flipped along the X axis to match the as-built instrument and measured PSFs [:pr:`275`, :user:`KevinVolkSTScI`, :user:`anand0xff`, :user:`mperrin`]
- Updated FGS throughput values to use data from the instrument sub-level testing that was done by Comdev/Honeywell, detector quantum efficiency as measured by Teledyne, and the OTE throughput from Lightsey 2012. The throughput file was also updated to include the WAVEUNIT keyword, which removes a warning. [:user:`shanosborne`]]

**WFIRST Improvements**

- *The WFI optical model has been updated to use optical data from the Cycle 8 design revision.* These include updated Zernike coefficients for field-dependent wavefront error, and masked and unmasked pupil images for each SCA, and updated filter throughputs (consistent with values used in Pandeia 1.4.2). The correct pupil file will automatically be selected for each calculation based on the chosen detector position and filter.   The pupil files are consistent with those provided in the WFI cycle 8 reference information, but have been resampled onto a common pixel scale.  See `WFIRST instrument model details <roman.html>`_ for more.  [:pr:`309` :user:`robelgeda`]
- Note, WFI's filters have been renamed so they all begin with “F”; `see the table here <https://github.com/spacetelescope/webbpsf/pull/309>`_ .
- *The WFI wavelength range has now been extended to cover the 0.48 - 2.0 µm range.* [:pr:`309` :user:`robelgeda`]
- *Expanded ``psf_grid`` method’s functionality so it can also be used to make grids of WFIRST PSFs.* Note that focal plane distortion is not yet implemented for WFIRST PSFs and so ``add_distortion`` keyword should not be used for this case. [:pr:`294`, :user:`shanosborne`]
- *The WFIRST F062 filter bandpass red edge was corrected* from 8000A to 7600A, and associated unit tests were updated to include F062  [:pr:`288`, :user:`robelgeda`]
- The WFI simulations now include the pointing jitter model, using the predicted WFI pointing stability of 14 milliarcseconds per axis. [:pr:`322`, :user:`mperrin`]

**General bug fixes and small changes:**

- *Many improvements in the PSF Grid functionality for generating photutils.GriddedPSFModels*:

  - New options in ``psf_grid`` to specify both/either the output filename and output directory location. See this `Jupyter notebook (ex2) <https://github.com/spacetelescope/webbpsf/blob/stable/notebooks/Gridded_PSF_Library.ipynb>`_ for examples. [:pr:`294`, :user:`shanosborne`]
  - sFfilenames when saving out a ``psf_grid`` FITS object which has it’s ``filename`` parameter set will now end with ``_det.fits`` instead of the previous ``_det_filt.fits`` [:pr:`294`, :user:`shanosborne`]
  - Update added to ``utils.to_griddedpsfmodel`` where a 2-dimensional array input with a header containing only 1 ``DET_YX`` keyword can be turned into ``GriddedPSFModel`` object without error as it  implies the case of a PSF grid with num_psfs = 1. [:pr:`294`, :user:`shanosborne`]
  - Remove deletion of ``det_yx`` and ``oversamp`` keywords from ``psf_grid`` output to allow for easier implementation in certain cases. Normal case users will have extra keywords but will not change functionality [:pr:`291`, :user:`shanosborne`]
  - Updated normalization of PSFs from ``psf_grid`` to be in surface brightness units, independent of oversampling in order to match the expectation of ``photutils.GriddedPSFModel``. This is different than webbpsf's default in which PSFs usually sum to 1 so the counts/pixel varies based on sampling. [:pr:`311`, :user:`mperrin`]
  - Fix bug in how ``pupilopd`` keyword is saved and include extra keywords ``opd_file``, ``opdslice``, ``coronmsk``, and ``pupil`` in the ``psf_grid`` output, both the GriddedPSFModel meta data and FITS object's header [:pr:`284`, :pr:`293`, :pr:`299`, :user:`shanosborne`]

- The ``set_position_from_aperture_name`` method now correctly sets the detector position parameter in the science frame [:pr:`281`, :user:`shanosborne`, :user:`JarronL`, :user:`mperrin`]
- Fix OPD HDUList output from the ``as_fits`` method inside the OPD class to include the previously existing header information [:pr:`270` :user:`laurenmarietta`]
- Added support for secondary mirror moves to the move_sur() method through the move_sm_local method [:pr:`295`, :user:`AldenJurling`]
- Remove ``units`` keyword from ``get_opd`` method, now the wave input needs to be a Wavefront object [:pr:`304`, :user:`shanosborne`]

**Software and Package Infrastructure Updates:**

- Added ``environment.yml`` file [:pr:`321`, :user:`shanosborne`, :user:`mperrin`]
- Remove leftover deprecated syntax ``_getOpticalSystem`` for ``_get_optical_system`` and ``display_PSF`` for ``display_psf`` [:pr:`280`, :pr:`294`, :user:`mperrin`, :user:`shanosborne`]
- Various smaller code cleanup and doc improvements, including code cleanup for better Python PEP8 style guide compliance [:user:`mperrin`, :user:`shanosborne`, :user:`robelgeda`]
- Documentation added and/or updated for a variety of features [:pr:`277`, :pr:`280`, :pr:`318`, :user:`mperrin, @shanosborne`]


--------




Version 0.8.0
=============

*2018 Dec 15*

This release focused on software engineering improvements, rather than changes in any of the optical models or reference data. (In particular, there are NO changes in the reference data files; the contents of the WebbPSF version 0.8 data zip file are identical to the reference data as distributed for version 0.7.  This version of WebbPSF will work with either of those interchangeably.).

.. admonition:: Python version support: Python 3 required

        This version drops support for Python 2.7. The minimum supported version of Python is now 3.5.

**New functionality:**

- *Added new capability to create grids of fiducial, distorted PSFs* spanning a chosen instrument/detector. This new ``psf_grid`` method is meant to be used as the first step of using the ``photutils`` package to do PSF-fitting photometry on simulated JWST PSFs. This method will output a list of or single ``photutils`` ``GriddedPSFModel`` object(s) which can then be read into ``photutils`` to apply interpolation to the grid and simulate a spatially dependent PSF anywhere on the instrument. See this `Jupyter notebook (ex3) <https://github.com/spacetelescope/webbpsf/blob/stable/notebooks/Gridded_PSF_Library.ipynb>`_ for examples. This method requires ``photutils`` version 0.6 or higher. [`#241, <https://github.com/spacetelescope/webbpsf/pull/241>` _, @shanosborne with inputs from @mperrin, @larrybradley, @hcferguson, and @eteq]

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
    future releases will require Python 3.5+. `See here <https://python3statement.org>`_ for more information on migrating to Python 3.

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
  for the other SIs. See this `Jupyter notebook (ex4) <https://github.com/mperrin/webbpsf/blob/stable/notebooks/Distortion_examples.ipynb>`_ for
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
  it is different physics. See this `Jupyter notebook (ex5) <https://github.com/mperrin/webbpsf/blob/stable/notebooks/Distortion_examples.ipynb>`_ for
  example output. For F560W through F1000W this is a much more obvious effect than the subtle distortions. [`#209,
  <https://github.com/mperrin/webbpsf/pull/209>`_, @shanosborne]
- *Added new capabilities for modeling mirror moves of the JWST primary
  segments and secondary mirror*, using a linear optical model to adjust OPDs.
  Added a new `notebook demonstrating these capabilities
  <https://github.com/mperrin/webbpsf/blob/stable/notebooks/Simulated%20OTE%20Mirror%20Move%20Demo.ipynb>`_.
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

* Improvements to parallelization code. Better documentation for parallelization.  PyFFTW3 replaced with pyFFTW for optimized
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

* Minor fix to MIRI ND filter transmission curve (Note: MIRI ND data is available on internal STScI data distribution only)
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




