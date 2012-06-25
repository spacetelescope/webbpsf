

=============================
More Examples
=============================


Any user of Webbpsf is invited to submit snippets of example code for sharing here. 


Typical Usage Cases
^^^^^^^^^^^^^^^^^^^^^^^


Iterating over multiple OPDs and filters
----------------------------------------

Perhaps you want to calculate PSFs for all filters of a given instrument, using all 10 available simulated OPDs:

    >>> def niriss_psfs():
    >>>     niriss = webbpsf.NIRISS()
    >>> 
    >>>     opdname = niriss.pupilopd
    >>> 
    >>>     for i in range(10):
    >>>         niriss.pupilopd = (opdname,i)
    >>>         for filtname in niriss.filter_list:
    >>>             niriss.filter=filtname
    >>>             fov=18
    >>>             outname = "PSF_NIRISS_%scen_wfe%d.fits" % (filtname, i)
    >>>             psf = webbpsf.calc_or_load_PSF(outname, niriss, nlambda=1, oversample=4, fov_arcsec=fov, rebin=True, display=True)
    >>> 



NIRCam coronagraphy with an offset source
-----------------------------------------

>>> from webbpsf import *
>>> nc = NIRCam()
>>> nc.image_mask='MASK430R'
>>> nc.pupil_mask='CIRCLYOT'
>>> nc.options['source_offset_r'] = 0.020       # source is 20 mas from center of coronagraph     
>>> nc.options['source_offset_theta'] = 45      # at a position angle of 45 deg
>>> nc.calcPSF('coronagraphic.fits', oversample=8)   # create highly oversampled output image


Create monochromatic PSFs across an instrument's entire wavelength range
-----------------------------------------------------------------------------
Monochromatic PSFs with steps of 0.1 micron from 5-28.3 micron.

>>> m = webbpsf.MIRI()
>>> m.pupilopd = 'OPD_RevV_miri_421.fits'       # select an OPD
>>>                                             # looks inside $WEBBPSF_DATA/MIRI/OPD by default
>>>                                             # or you can specify a full path name. 
>>> m.options['parity'] = 'odd'                 # please make an output PSF with its center
>>>                                             # aligned to the center of a single pixel
>>>
>>> waves = N.linspace(5.0, 28.3, 234)*1e-6     # iterate over wavelengths in meters
>>> for w in waves:
>>>     m.calcPSF(fov_arcsec=30, oversample=4, rebin=True, monochromatic=wavelength, display=True,
>>>                outfile='psf_MIRI_mono_%.1fum_revV_opd1.fits' % (wavelength*1e6))


Simulate NIRCam coronagraphic acquisition images
--------------------------------------------------



    >>> def compute_psfs():
    >>>     nc = webbpsf.NIRCam()
    >>> 
    >>>     # acq filter, occulting mask, lyot, coords of acq ND square
    >>>     sets = [('F182M', 'MASKSWB', 'WEDGELYOT', -10,  7.5),
    >>>             ('F182M', 'MASK210R', 'CIRCLYOT', -7.5, 7.5),
    >>>             ('F335M', 'MASKLWB', 'WEDGELYOT',  7.5, 7.5),
    >>>             ('F335M', 'MASK335R', 'CIRCLYOT', -10,  7.5)]
    >>> 
    >>>     nlambda = 9     
    >>>     oversample = 2  
    >>> 
    >>>     calc_oversample=4
    >>> 
    >>> 
    >>>     fov_arcsec = 25
    >>> 
    >>> 
    >>> 
    >>>     for param in sets:
    >>>         nc.filter = param[0]
    >>>         nc.image_mask = param[1]
    >>>         nc.pupil_mask = param[2]
    >>>         source_offset_x = param[3]
    >>>         source_offset_y = param[4]
    >>> 
    >>> 
    >>>         source_offset_r = N.sqrt(source_offset_x**2+ source_offset_y**2)
    >>>         source_offset_theta = N.arctan2(source_offset_x, source_offset_y)*180/N.pi
    >>>         nc.options['source_offset_r'] = source_offset_r
    >>>         nc.options['source_offset_theta'] = source_offset_theta
    >>> 
    >>> 
    >>>         filename = "PSF_NIRCam_%s_%s_%s_offset.fits" % (param[0], param[1], param[2])
    >>>         result = nc.calcPSF(nlambda=nlambda, oversample=oversample, calc_oversample=calc_oversample, fov_arcsec=fov_arcsec, outfile=filename, display=False)
    >>> 
    >>> 

Iterate a calculation over all MIRI corongraphic modes
-------------------------------------------------------

    >>> def miri_psfs_coron():
    >>>     miri = webbpsf.MIRI()
    >>> 
    >>>     for filtwave in [1065, 1140, 1550, 2300]:
    >>> 
    >>>         miri.filter='F%4dC' % filtwave
    >>>         if filtwave<2000:
    >>>             miri.image_mask='FQPM%4d' % filtwave
    >>>             miri.pupil_mask='MASKFQPM'
    >>>             fov=24
    >>>         else:
    >>>             miri.image_mask='LYOT2300'
    >>>             miri.pupil_mask='MASKLYOT'
    >>>             fov=30
    >>> 
    >>> 
    >>>         offset_x = 0.007 # arcsec
    >>>         offset_y = 0.007 # arcsec
    >>> 
    >>>         miri.options['source_offset_r'] = N.sqrt(offset_x**2+offset_y**2) # offset in arcsec
    >>>         miri.options['source_offset_theta'] = N.arctan2(-offset_x, offset_y)*180/N.pi # PA in deg
    >>> 
    >>> 
    >>>         outname = "PSF_MIRI_%s_x%+05.3f_y%+05.3f.fits" % (miri.image_mask, offset_x, offset_y)
    >>>         psf = webbpsf.calc_or_load_psf(outname, miri, oversample=4, fov_arcsec=fov, display=True)
    >>> 
    >>> 

Make plots of encircled energy in PSFs at various wavelengths
----------------------------------------------------------------

    >>> def miri_psfs_for_ee():
    >>>     miri = webbpsf.MIRI()
    >>> 
    >>>     opdname = miri.pupilopd
    >>> 
    >>>     for i in range(10):
    >>>         miri.pupilopd = (opdname,i)
    >>>         for wave in [5.0, 7.5, 10, 14]:
    >>> 
    >>>             fov=18
    >>> 
    >>>             outname = "PSF_MIRI_%.1fum_wfe%d.fits" % (wave, i)
    >>>             psf = webbpsf.calc_or_load_psf(outname, miri, monochromatic=wave*1e-6, oversample=4, fov_arcsec=fov, rebin=True, display=True)
    >>> 
    >>> 
    >>> 
    >>> def plot_ee_curves():
    >>>     pl.clf()
    >>>     for iw, wave in enumerate([5.0, 7.5, 10, 14]):
    >>> 
    >>>         ees60 = []
    >>>         ees51 = []
    >>>         ax = pl.subplot(2,2,iw+1)
    >>>         for i in range(10):
    >>>             name = "PSF_MIRI_%.1fum_wfe%d.fits" % (wave, i)
    >>>             webbpsf.display_EE(name, ax=ax, mark_levels=False)
    >>> 
    >>>             eefn = webbpsf.measure_EE(name)
    >>>             ees60.append(eefn(0.60))
    >>>             ees51.append(eefn(0.51))
    >>> 
    >>>         ax.text(1, 0.6, 'Mean EE inside 0.60": %.3f' % np.asarray(ees60).mean())
    >>>         ax.text(1, 0.5, 'Mean EE inside 0.51": %.3f' % np.asarray(ees51).mean())
    >>> 
    >>>         ax.set_title("Wavelength = %.1f $\mu$m" % wave)
    >>> 
    >>>         ax.axvline(0.6, ls=":", color='k')
    >>>         ax.axvline(0.51, ls=":", color='k')
    >>> 
    >>> 
    >>>     pl.tight_layout()
    >>> 

Advanced Usage Tricks
^^^^^^^^^^^^^^^^^^^^^^^

This section serves as a catch-all for other example codes, possibly more esoteric in application. 

Writing out only downsampled images
-----------------------------------

Perhaps you may want to calculate the PSF using oversampling, but to save disk space you only want to write out the PSF downsampled to detector resolution.

   >>> result =  inst.calcPSF(args, ...)
   >>> result['DET_SAMP'].writeto(outputfilename)

Or if you really care about writing it as a primary HDU rather than an extension, replace the 2nd line with

   >>> pyfits.PrimaryHDU(data=result['DET_SAMP'].data, header=result['DET_SAMP'].header).writeto(outputfilename)


Providing your own OPDs from some other source
-----------------------------------------------






Modifying existing OPDs to add defocus
----------------------------------------

Perhaps you want to modify the OPD used for a given instrument, for instance to
add a defocus. You can do this by subclassing one of the existing instrument
classes to patch over the _getOpticalSystem function. An OpticalSystem is
basically a list so it's straightforward to just add another optic there. In
this example it's a lens for defocus but you could just as easily add another
FITSOpticalElement instead to read in a disk file.


    >>> class TF_with_defocus(webbpsf.TFI):
    >>>         def __init__(self, \*args, \*\*kwargs):
    >>>                 webbpsf.TFI.__init__(self, \*args, \*\*kwargs)
    >>>                 # modify the following as needed to get your desired defocus
    >>>                 self.defocus_waves = 0
    >>>                 self.defocus_lambda = 4e-6
    >>>         def _getOpticalSystem(self, \*args, \*\*kwargs):
    >>>                 osys = webbpsf.TFI._getOpticalSystem(self, \*args, \*\*kwargs)
    >>>                 lens = poppy.ThinLens(name='my lens', nwaves=self.defocus_waves, reference_wavelength=self.defocus_lambda)  
    >>>                 lens.planetype=poppy.PUPIL # needed to flag plane location for the propagation algorithms
    >>>                 osys.planes.insert(1, lens)
    >>>                 return osys
    >>> 
    >>> tf2 = TF_with_defocus()
    >>> tf2.defocus= 4  # means 4 waves of defocus at the wavelength defined by tf2.defocus_lambda
    >>> psf = tf2.calcPSF()
    >>> 


Simulate coronagraphy with pupil shear, saving the wavefront intensity in the Lyot pupil plane
------------------------------------------------------------------------------------------------


This is an example of a much more complicated calculation, including code to generate publication-quality plots. 



There are two functions here, one that creates a simulated PSF for a given amount of shear, and one that makes some nice plots of it.

    >>> def miri_psf_sheared(shearx=0, sheary=0, nopds = 1, display=True, overwrite=False, \*\*kwargs):
    >>>     """ Compute MIRI coronagraphic PSFs assuming pupil shear between the MIRI lyot mask and the OTE
    >>> 
    >>>     Parameters
    >>>     ------------
    >>>     shearx, sheary: float
    >>>         Shear across the pupil expressed in percent, i.e. shearx=3 means the coronagraph pupil is sheared by 3% of the primary.
    >>> 
    >>>     """
    >>>     miri = webbpsf.MIRI()
    >>> 
    >>>     miri.options['pupil_shift_x'] = shearx/100 # convert shear amount to float between 0-1
    >>>     miri.options['pupil_shift_y'] = sheary/100
    >>> 
    >>>     opdname = miri.pupilopd         # save default OPD name for use in iterating over slices
    >>> 
    >>>     filtsets = [('F1065C', 'FQPM1065', 'MASKFQPM'), ('F2300C','LYOT2300','MASKLYOT')]
    >>> 
    >>>     fov=10
    >>> 
    >>>     for i in range(nopds):
    >>>         miri.pupilopd = (opdname,i)
    >>>         for filt, im_mask, pup_mask in filtsets:
    >>>             print("Now computing OPD %d for %s, %s, %s" % (i, filt, im_mask, pup_mask))
    >>>             miri.filter=filt
    >>>             miri.image_mask = im_mask
    >>>             miri.pupil_mask = pup_mask
    >>> 
    >>> 
    >>>             outname = "PSF_MIRI_%s_wfe%d_shx%.1f_shy%.1f.fits" % (filt, i, shearx, sheary)
    >>>             outname_lyot = outname.replace("PSF_", 'LYOTPLANE_')
    >>> 
    >>> 
    >>>             if os.path.exists(outname) and not overwrite:
    >>>                 print ("File %s already exists. Skipping and continuing for now... set overwrite=True to recalculate" % outname)
    >>>                 return
    >>> 
    >>>             psf, intermediates = miri.calcPSF(oversample=4, fov_arcsec=fov, rebin=True, display=display, return_intermediates=True, \*\*kwargs)
    >>> 
    >>>             lyot_intensity = intermediates[4]
    >>> 
    >>>             psf.writeto(outname, clobber=True)
    >>>             lyot_intensity.writeto(outname_lyot, clobber=True, includepadding=False)
    >>> 
    >>> 
    >>> def plot_sheared_psf(shearx=1.0, sheary=0, lyotmax=1e-5, psfmax = 1e-3, diffmax=10):
    >>>     i = 0
    >>>     filtsets = [('F1065C', 'FQPM1065', 'MASKFQPM')]#, ('F2300C','LYOT2300','MASKLYOT')]
    >>> 
    >>>     pl.clf()
    >>>     pl.subplots_adjust(left=0.02, right=0.98, wspace=0.3)
    >>>     for filt, im_mask, pup_mask in filtsets:
    >>>         perfectname = "PSF_MIRI_%s_wfe%d_shx%.1f_shy%.1f.fits" % (filt, i, 0,0)
    >>>         perfectname_lyot = perfectname.replace("PSF_", 'LYOTPLANE_')
    >>> 
    >>> 
    >>>         outname = "PSF_MIRI_%s_wfe%d_shx%.1f_shy%.1f.fits" % (filt, i, shearx, sheary)
    >>>         outname_lyot = outname.replace("PSF_", 'LYOTPLANE_')
    >>> 
    >>>         if not os.path.exists(outname):
    >>>             print "File %s does not exist, skipping" % outname
    >>>             return False
    >>> 
    >>> 
    >>>         #psf = pyfits.open(outname)
    >>>         #perfpsf = pyfits.open(perfectname)
    >>>         lyot = pyfits.open(outname_lyot)
    >>>         perflyot = pyfits.open(perfectname_lyot)
    >>> 
    >>>         wzero = np.where(lyot[0].data == 0)
    >>>         wzero = np.where(lyot[0].data < 1e-15)
    >>>         lyot[0].data[wzero] = np.nan
    >>>         wzero = np.where(perflyot[0].data == 0)
    >>>         perflyot[0].data[wzero] = np.nan
    >>> 
    >>>         cmap = matplotlib.cm.jet
    >>>         cmap.set_bad('gray')
    >>> 
    >>> 
    >>> 
    >>>         # plot comparison perfect case Lyot Intensity
    >>>         ax = pl.subplot(231)
    >>>         pl.imshow(perflyot[0].data, vmin=0, vmax=lyotmax, cmap=cmap)
    >>>         pl.title("Lyot plane, no shear")
    >>>         ax.yaxis.set_ticklabels("")
    >>>         ax.xaxis.set_ticklabels("")
    >>> 
    >>>         wg = np.where(np.isfinite(perflyot[0].data))
    >>>         ax.set_xlabel("Residual flux = %.1f%%" % (perflyot[0].data[wg].sum()*100))
    >>> 
    >>>         # plot shifted pupil Lyot intensity
    >>>         ax = pl.subplot(234)
    >>>         pl.imshow(lyot[0].data, vmin=0, vmax=lyotmax, cmap=cmap)
    >>>         pl.title("Lyot plane, shear (%.1f, %.1f)" % (shearx, sheary))
    >>>         ax.yaxis.set_ticklabels("")
    >>>         ax.xaxis.set_ticklabels("")
    >>>         wg = np.where(np.isfinite(lyot[0].data))
    >>>         ax.set_xlabel("Residual flux = %.1f%%" % (lyot[0].data[wg].sum()*100))
    >>> 
    >>> 
    >>> 
    >>>         # Radial profile plot
    >>>         pl.subplot(233)
    >>> 
    >>>         radius, profperf = webbpsf.radial_profile(perfectname, ext=1)
    >>>         radius2, profshear = webbpsf.radial_profile(outname, ext=1)
    >>> 
    >>>         # normalize all radial profiles to peak=1 for an unocculted source
    >>>         radiusu, profunocc = webbpsf.radial_profile('PSF_MIRI_F1065C_wfe0_noshear_unocculted.fits', ext=1, center=(43.3, 68.6)) # center is in pixel coords
    >>> 
    >>>         peakunocc = profunocc.max()
    >>>         profperf /= peakunocc
    >>>         profshear/= peakunocc
    >>>         profunocc/= peakunocc
    >>> 
    >>> 
    >>>         pl.semilogy(radius, profperf, label="No shear")
    >>>         pl.semilogy(radius2, profshear, label="shear (%.1f, %.1f)" % (shearx, sheary))
    >>>         pl.semilogy(radiusu, profunocc, label="Unocculted", ls=":" )
    >>> 
    >>> 
    >>>         pl.xlabel("Separation [arcsec]")
    >>>         pl.ylabel("Relative Intensity")
    >>>         pl.legend(loc='upper right')
    >>>         pl.gca().set_xlim(0,6)
    >>> 
    >>> 
    >>>         # plot comparison perfect case PSF - detector sampled
    >>>         pl.subplot(232)
    >>>         webbpsf.display_PSF(perfectname, ext=1, vmax=psfmax)
    >>>         pl.title("PSF, no shear")
    >>> 
    >>>         # plot shifted pupil PSF - detector sampled
    >>>         pl.subplot(235)
    >>>         webbpsf.display_PSF(outname, ext=1, vmax=psfmax)
    >>>         pl.title("PSF, shear (%.1f, %1.f)" % (shearx, sheary))
    >>>         pl.xlabel("Separation [arcsec]")
    >>>         # difference PSf
    >>>         pl.subplot(236)
    >>>         webbpsf.display_PSF_difference(outname, perfectname, ext1=1, ext2=1, vmax=diffmax, vmin=-0.1, normalize_to_second=True)
    >>>         pl.title('Relative PSF increase')
    >>>         pl.xlabel("Separation [arcsec]")
    >>> 
    >>> 
    >>>         #pl.tight_layout()
    >>>         return True
    >>> 
    >>> 
    >>> 





..
  Copy in some examples here from test_webbpsf and validate_webbpsf ? 


