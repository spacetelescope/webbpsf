.. _jwst_instruments:


*****************************
JWST Instrument Model Details
*****************************


The following describes specific details for the various JWST instrument classes. See also :ref:`the references page <references>` for information on data sources. 


One general note is that the ``webbpsf`` class interfaces do not attempt to exactly
model the implementation details of all instrument mechanisms, particularly for
NIRCam and NIRISS that each have multiple wheels. The
``filter`` attribute of a given class is used to select any and all filters,
even if as a practical matter a given filter is physically installed in a
"pupil" wheel instead of a "filter" wheel. Likewise any masks that affect the
aperture shape are selected via the ``pupil_mask`` attribute even if physically
an optic is in a so-called "filter" wheel.

All classes share some common attributes:

 * ``filter`` which is the string name of the currently selected filter. 
   The list of available filters is provided as the ``filter_list`` attribute.
 * ``image_mask`` which is the name of a selected image plane element such as a
   coronagraph mask or spectrograph slit, or ``None`` to indicate no 
   such element is present.  
   The list of available options is provided as the ``image_mask_list`` attribute.
 * ``pupil_mask`` likewise allows specification of any optic that modifies the pupil plane
   subsequent to the image mask, such as a coronagraphic Lyot stop or spectrograph grating stop.
   The list of available options is provided as the ``pupil_mask_list`` attribute.
 * Each SI has a ``detector`` attribute that can be used to select among its
   multiple detectors (if more than one are present in that SI), and a
   ``detector_position`` attribute which is a 2-tuple giving the pixel coordinates
   on that detector for the center location in any calculated output PSF.
   Note that the ``detector_position`` value should be
   specified using the Python (Y,X) axes order convention.

.. admonition:: Instrument measured WFE coming in a future release

    This version of WebbPSF does not yet contain data on the measured SI
    wavefront errors. However it does contain software infrastructure to
    support that functionality in the near future, once the data files have
    been fully reviewed by the project and SI teams. So for the time being,
    some options related to detector field positions and SI wavefront error
    are present in WebbPSF, but do not have any effect on the output PSFs.
    Stay tuned!


.. warning::
    WebbPSF provides some sanity checking on user inputs, but does not strive to
    strictly forbid users from trying to simulate instrument configurations that 
    may not be achievable in practice.  Users are responsible for knowing the
    available modes well enough to be aware if they are trying to
    simulate an inconsistent or physically impossible configuration. 


NIRCam
======

Imaging
--------

NIRCam is one of the more complicated classes in ``webbpsf``, and has several unique selectable options to model the two copies of NIRCam each with two channels..  

The ``detector`` attribute can be used to select between any of the ten detectors,
A1-A5 and B1-B5.  Additional attributes are then automatically set for ``channel``
("short" or "long") and module ("A" or "B") but these cannot be set directly;
just set the desired detector and the channel and module are inferred
automatically. 

The choice of ``filter`` also impacts the channel selection: If you choose a
long-wavelength filter such as F460M, then the detector will automatically
switch to the long-wave detector for the current channel. For example, if the
detector was previously set to A2, and the user enters ``nircam.filter = "F460M"``
then the detector will automatically change to A5.  If the user later selects
``nircam.filter = "F212N"`` then the detector will switch to A1 (and the user will
need to manually select if a different short wave detector is desired).  This
behavior on filter selection can be disabled by setting ``nircam.auto_channel = False``. 


Coronagraph Masks
------------------

The coronagraph image-plane masks and pupil-plane Lyot masks are all included as options.
These are based on the nominal design properties as provided by the NIRCam team,
not on any specific measurements of the as-built masks. The simulations of the occulting mask
fields also include the nearby neutral density squares for target acquisitions.

WebbPSF won't prevent users from simulating configuration using a coronagraph
image mask without the Lyot stop, but that's not something that can be done for
real with NIRCam. 


Weak Lenses for Wavefront Sensing
---------------------------------

WebbPSF includes models for the three weak lenses used for wavefront sensing, including the 
pairs of lenses that can be used together simultaneously.

The convention is such that the "negative" 8 waves lens is concave, the
"positive" two lenses are convex. Thus positive weak lenses move best focus
in front of the detector, or equivalently the electric field imaged on the detector
becomes behind or beyond best focus. Negative weak lenses move best focus behind the detector, 
or equivalently the image on the detector is moved closer to the OTE exit pupil
than best focus. 

Note that the weak lenses are in the short wave channel only; webbpsf won't stop
you from simulating a LW image with a weak lens, but that's not a real
configuration that can be acheived with NIRCam.


SI WFE
------

(Not yet available)

The SI internal WFE measurements are distinct for each of the modules and
channels. When enabled, these are added to the final pupil of the optical
train, i.e. after the coronagraphic image planes. 


NIRSpec
=======

Imaging and spectroscopy
------------------------

webbpsf models the optics of NIRSpec, mostly in **imaging** mode or for monochromatic PSFs that can be assembled into spectra using other tools.

This is not a substitute for a spectrograph model, but rather a way of
simulating a PSF as it would appear with NIRSpec in imaging mode (e.g. for
target acquisition).  It can also be used to produce monochromatic PSFs
appropriate for spectroscopic modes, but other software must be used for
assembling those monochromatic PSFs into a spectrum.

Slits: webbpsf includes models of each of the fixed slits in NIRSpec (S200A1, S1600A1, and so forth), plus a
few patterns with the MSA: (1) a single open shutter, (2) three adjacent
open shutters to make a mini-slit, and (3) all shutters open at once.
Other MSA patterns could be added if requested by users. 

By default the ``pupil_mask`` is set to the "NIRSpec grating" pupil mask.  In
this case a rectangular pupil mask 8.41x7.91 m as projected onto the primary is
added to the optical system at the pupil plane after the MSA. This is an
estimate of the pupil stop imposed by the outer edge of the grating clear
aperture, estimated based on optical modeling by Erin Elliot and Marshall
Perrin.


SI WFE
------

(Not yet available)

SI WFE will most likely be added to the entrance pupil, prior to the MSA image plane. This model is still under development. 

NIRISS
======


Imaging and AMI
----------------

WebbPSF models the direct imaging and nonredundant aperture masking interferometry modes of NIRISS in the usual manner.

Note that long wavelength filters (>2.5 microns) are used with a pupil
obscuration which includes the pupil alignment reference fixture. This is called
the "CLEARP" pupil.

Based on the selected filter, webbpsf will automatically toggle the
``pupil_mask`` between "CLEARP" and the regular clear pupil (i.e.
``pupil_mask = None``).


Slitless Spectroscopy
---------------------

webbpsf provides preliminary support for 
the single-object slitless
spectroscopy ("SOSS") mode using the GR700XD cross-dispersed grating. Currently
this includes the clipping of the pupil due to the undersized grating and its
mounting hardware, and the cylindrical lens that partially defocuses the light
in one direction.

.. warning::

    Prototype implementation - Not yet fully tested or verified.

Note that WebbPSF does not model the spectral dispersion in any of NIRISS'
slitless spectroscopy modes.  For wide-field slitless spectroscopy, this
can best be simulated by using webbpsf output PSFs as input to the aXe
spectroscopy code. Contact Van Dixon at STScI for further information.
For SOSS mode, contact Loic Albert at Universite de Montreal.

The other two slitless spectroscopy grisms use the regular pupil and do not require any special
support in WebbPSF; just calculate monochromatic PSFs at the desired wavelengths
and assemble them into spectra using tools such as aXe. 

Coronagraph Masks
------------------

NIRISS includes four coronagraphic occulters, machined as features on its
pick-off mirror. These were part of its prior incarnation as TFI, and are not
expected to see much use in NIRISS. However they remain a part of the physical
instrument and we retain in webbpsf the capability to simulate them.

SI WFE
-------

(Not yet available)

The SI internal WFE measurements are distinct for each of the modules and
channels. When enabled, these are added to the final pupil of the optical
train, i.e. after the coronagraphic image planes. 


MIRI
====

Imaging
-------

webbpsf models the MIRI imager; currently there is no specific support for MRS,
however monochromatic PSFS computed for the imager may be used as a reasonable
proxy for PSF properties at the entrance to the MRS slicers.


Coronagraphy
-------------

webbpsf includes models for all three FQPM coronagraphs and the Lyot
coronagraph. In practice, the wavelength selection filters and the Lyot stop are
co-mounted. webbpsf models this by automatically setting the ``pupil_mask``
element to one of the coronagraph masks or the regular pupil when the ``filter``
is changed. If you want to disable this behavior, set ``miri.auto_pupil = False``.


LRS Spectroscopy
----------------

webbpsf includes models for the LRS slit and the subsequent pupil stop on the
grism in the wheels. Users should select ``miri.image_mask = "LRS slit"`` and ``miri.pupil_mask = 'P750L LRS grating'``.
That said, the LRS simulations have not been extensively tested yet; 
feedback is appreciated about any issues encountered.


SI WFE
------

(Not yet available)

The SI internal WFE measurements, when enabled, are added to the final pupil of the optical
train, i.e. after the coronagraphic image planes. 


FGS
===

The FGS class does not have any selectable optical elements (no filters or
image or pupil masks of any kind). Only the detector is selectable, between
either 'FGS1' or 'FGS2'. 

SI WFE
------

(Not yet available)

