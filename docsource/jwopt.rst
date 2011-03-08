.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. module:: jwopt

==========================
The JWInstrument interface
==========================


This module provides the primary interface for interactive non-GUI use and programming. It provides 
five classes corresponding to the JWST instruments with consistent interfaces for each.  Each instrument knows the details of its own
detector sampling, available filter complement (and other mechanisms), and orientation with respect to the OTE.

This module also contains a handful of utility functions for measuring and plotting PSF properties.

 * :ref:`Dive into some examples <examples>`
 
   * :ref:`First steps <examples_1>`
   * :ref:`OPD files <examples_2>`
   * :ref:`Input source spectra <examples_3>`
   * :ref:`Field of view and centering <examples_4>`
   * :ref:`Pixel scales and sampling <examples_5>`
   * :ref:`Offset source positions <examples_6>`

 * :ref:`Reference: common interface <jwinst>`
 * :ref:`Reference: individual instruments <specific_instruments>`
 * :ref:`Utility functions <helpers>`




.. _examples :

Examples
--------

.. _examples_1 :

**First steps:**

Simple PSFs are easily obtained. 

>>> nc = NIRCam()
>>> nc.filter =  'F200W'
>>> psf = nc.calcPSF(oversample=4)      # returns a pyfits.HDUlist containing PSF and header
>>> pylab.imshow(psf[0].data)           # display it on screen yourself, or
>>> display_PSF(psf)                    # use this convenient function to make a nice log plot with labeled axes
>>>
>>> psf = nc.calcPSF(filter='F470N', oversample=4)    # this is just a shortcut for setting the filter, then computing a PSF
>>>
>>> psf.writeto("myPSF.fits", clobber=True)   # the returned PSF is just a regular pyfits HDUlist object

For interactive use, you can have the PSF displayed as it is computed:

>>> nc.calcPSF(display=True)                          # will make nice plots with matplotlib.

More complicated instrumental configurations are available by setting the instrument's attributes. For instance,
one can create an instance of MIRI and configure it for coronagraphic observations, thus:

>>> miri = MIRI()
>>> miri.filter = 'F1065C'
>>> miri.image_mask = 'FQPM1065'
>>> miri.pupil_mask = 'MASKFQPM'
>>> miri.calcPSF('outfile.fits')

.. _examples_2 :

**Choosing Wavefront Error (OPD) files**

This data library included with this code provides several Rev V OPDs for each instrument. You can also easily select another file
by changing the pupilopd attribute. 

>>> print tfi.pupilopd
'TFI_OPDisim1.fits'
>>> tfi.pupilopd = "/path/to/my/OPDfile.fits"

OPD files may be formatted as datacubes, with the third dimension being (for instance) time or different statistical realizations.
For such a file, specify the OPD as a tuple containing filename and slice:

>>> tfi.pupilopd = ("MyOPDCube.fits", 5)  # use 5th plane in the cube

If you select a cube OPD without explicitly naming the slice, the zeroth slice will be used.


.. _examples_3 :



**Choosing an input source spectrum:**

The spectrum of the observed source can be specified by giving a `source` parameter in the call to :py:meth:`calcPSF() <JWInstrument.calcPSF>`. The following are valid sources:

1. A dictionary with elements `source["wavelengths"]` and `source["weights"]` giving the wavelengths in meters and the relative weights for each. These should be numpy arrays or lists.

   >>> src = {'wavelengths': [2.0e-6, 2.1e-6, 2.2e-6], 'weights': [0.3, 0.5, 0.2]}
   >>> nc.calcPSF(source=src, outfile='psf_for_src.fits')

2. A tuple or list containing the numpy arrays `(wavelength, weights)` instead.
3. A `pysynphot.Spectrum` object. This is the best option.  **but not yet fully implemented - will be in version 0.2**


.. _examples_4 :


**Adjusting field of view and pixel centering:**

The field of view may be specified in units of either arcseconds or pixels.  If specified in arcseconds, 
the resulting FoV will be whatever the closest integer number of pixels is. For instance, 

>>> mynircam = NIRCam()
>>> result   = mynircam.calcPSF(fov_npixels = 512)
>>> result2  = mynircam.calcPSF(fov_arcsec=7, oversample=2, filter='F250M')


By default, the PSF will be located at the exact center of the output array. This means that if the PSF is computed on an array with an odd number of pixels, the
PSF will be centered exactly on the central pixel. If the PSF is computed on an array with even size, it will be centered on the "crosshairs" at the intersection of the central four pixels.
If one of these is particularly desirable to you, set the parity option appropriately:

>>>  mynircam.options['parity'] = 'even'
>>>  mynircam.options['parity'] = 'odd'

Setting one of these options will ensure that a field of view specified in arcseconds 
results in a number of pixels with the desired parity (i.e. the closest integer number of pixels with that parity is selected).
Alternatively, 
you may of course just set the desired number of pixels explicitly in the call to calcPSF().  Note that the odd/even parity options are applied with respect to
the *actual detector pixels*, so if you use oversampling then the subsampled array's parity will also depend on the oversampling factor.

.. _examples_5 :



**Pixel scales, sampling, and oversampling:**

The derived instrument classes all known their own instrumental pixel scales and compute PSFs at the proper scale by default. However, you can change the output 
pixel scale in a variety of ways, if desired. See the :ref:`JWInstrument.calcPSF` documentation for more details.

1. Set the `oversample` parameter to calcPSF. This will produce a PSF with a pixel grid more finely sampled by this factor. 
   `oversample=1` is the native detector scale, `oversample=2` means divide each pixel into 2x2 finer pixels, and so forth.
   You can automatically obtain both the oversampled PSF and a version rebinned onto the real detector pixel scale by setting `rebin=True` 
   in the call to calcPSF:

   >>> hdulist = instrument.calcPSF(oversample=2, rebin=True)   # hdulist will contain a primary HDU with the 
   >>>                                                          # oversampled data, plus an image extension 
   >>>                                                          # with the PSF rebinned down to regular sampling.

   

2. For coronagraphic calculations, it is possible to set different oversampling factors for different parts of the calculation. See the `calc_oversample` and `detector_oversample` parameters documented for `JWInstrument.calcPSF`. This
   is of no applicability to regular imaging calculations (for which `oversample` is effectively a synonym for `detector_oversample`).

   >>> tfi.calcPSF(calc_oversample=8, detector_oversample= 2)    # model the occulter with very fine pixels, then save the 
   >>>                                                           # data on a coarser (but still oversampled) scale

3. Or, if you need even more flexibility, just change the `instrument.pixelscale` attribute to be whatever arbitrary scale you require. 

   >>> instrument.pixelscale = 0.0314159


.. _examples_6 :



**Offset source positions**

The PSF may be shifted off-center by adjusting the offset of the stellar source. This is done in polar coordinates:

>>> instrument.options['source_offset_r'] = 0.3         # offset in arcseconds
>>> instrument.options['source_offset_theta'] = 45.     # degrees counterclockwise from detector +Y

If these options are set, the offset is applied relative to the center of the output array (as determined by the array size and parity as discussed above).  The angle convention is the usual astronomical one going counterclockwise from north: theta=0 is +Y, theta=90 is -X, and so forth. 

For coronagraphic modes, the coronagraph occulter is always assumed to be at the center of the output array. Therefore, these options let you offset the source away from the coronagraph.
 


.. _jwinst:

The JWInstrument generic class
--------------------------------

.. autoclass:: jwopt.JWInstrument
   :members:


.. _specific_instruments:

Specific Instruments
--------------------
.. autoclass:: jwopt.NIRCam
.. autoclass:: jwopt.NIRSpec
.. autoclass:: jwopt.MIRI
.. autoclass:: jwopt.TFI
.. autoclass:: jwopt.FGS

--------------

.. autofunction:: jwopt.Instrument


.. _helpers:

Helper Functions for Display & Measurement
------------------------------------------
.. autofunction:: jwopt.display_PSF
.. autofunction:: jwopt.radial_profile
.. autofunction:: jwopt.measure_radial
.. autofunction:: jwopt.measure_EE
.. autofunction:: jwopt.measure_fwhm
.. autofunction:: jwopt.measure_sharpness
.. autofunction:: jwopt.measure_anisotropy



--------------

Documentation last updated on |today|

