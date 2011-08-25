

=============================
More Examples
=============================



NIRCam coronagraphy with an offset source

>>> from webbpsf import *
>>> nc = NIRCam()
>>> nc.image_mask='MASK430R'
>>> nc.pupil_mask='CIRCLYOT'
>>> nc.options['source_offset_r'] = 0.020       # source is 20 mas from center of coronagraph     
>>> nc.options['source_offset_theta'] = 45      # at a position angle of 45 deg
>>> nc.calcPSF('coronagraphic.fits', oversample=8)   # create highly oversampled output image


Create monochromatic MIRI PSFs across its entire wavelength range, with steps of 0.1 micron from 5-28.3 micron.

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



Copy in some examples here from test_webbpsf and validate_webbpsf ? 


