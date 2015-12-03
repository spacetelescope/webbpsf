.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. _references:

Appendix: Instrument Property References
================================================================

We give here references for the instrumental properties assumed in PSF
computations, with particular attention to coronagraphic optics. It also notes
several places where the current models or available files are limited in some
manner that might be improved in a future release. 


*Note: The WebbPSF software and all of its associated data files are entirely ITAR-free.*

OTE
----

The supplied OPDs are the Mission CDR OPD simulation set, produced in March
2010 by Ball Aerospace staff (Paul Lightsey et al.) via the IPAM optical model
using Zernike WFE coefficients consistent with Revision V of the JWST
optical error budget.

**Note:** The provided files included no header metadata, and in particular no
pixel scale, so one was assumed based on the apparent pupil diameter in the
files. The estimated uncertainty in this scale is 1 part in 1000, so users concerned with measurements of PSF FWHMs etc at that level should be cautious. 

The current model pixel scale, roughly 6 mm/pixel, is too coarse to resolve well the edge roll-off around the border of each segment. We make no
attempt to include such effects here at this time. An independent study using much more finely sampled pupils has shown that the effect of segment edge roll-off is to scatter ~2% of the light from the PSF core out to large radii, primarily in the form of increased intensity along the diffraction spikes (Soummer et al. 2009, Technical Report JWST-STScI-001755)


NIRCam
------

NIRCam focal plane scale: 0.317 arcsec/pixel (short wave) 0.0648 arcsec/pixel (long wave).  From `STScI NIRCam web page <http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/#channels>`_. 

The coronagraph optics models are based on the NIRCam instrument team's series of SPIE papers describing the coronagraph designs and flight hardware. 
(Krist et al. 2007, 2009, 2010 Proc. SPIE), as clarified through cross checks with information provided by the NIRCam instrument team (Krist, private communication 2011).  Currently, the models include only the 5 arcsec square ND acquisition boxes and not the second set of 2 arcsec squares. 

.. comment
    Note that the NIRCam wedge BLCs both have 'flat' regions with constant FWHM at the extreme left and right
    sides of the wedge, as well as the region in the middle with varying FWHM. Though the widths of these flat 
    regions are not explicitly stated in either of Krist's papers, by inspection of the figures they appear to be
    ~ 2.5 arcsec wide, so the actual wedge is 15 arcsec in length.  **Note:** This should be double-checked with John Krist.
    **John says "Do not reference or distribute my memo. " so don't say the following **
    in the file "JWST NIRCam Lyot Stop Definitions" dated January 22, 2007. The
    provided mask data were in the form of pupil plane coordinates normalized
    by the telescope radius. A Python script was used to convert these
    coordinates into pixel mask files 1024x1024 pixels in size. This
    transformation included a bit of anti-aliasing such that greyscale values
    are used for pixels right along the border of curved or diagonal edges.
    However, this algorithm could probably be improved further.


Weak lenses: The lenses are nominally +- 8 and +4 waves at 2.14 microns. The as built defocus values are as follows based on component testing:  7.76198,
-7.74260, 3.90240. 


NIRSpec
--------
NIRspec field of view rotation: 41.5 degrees. Matt Lallo, draft SIAF information; and Ball SI Fields for WFS&C document, J. Scott Knight


NIRISS
-------

Occulting spots: Assumed to be perfect circles with diameters 0.58, 0.75, 1.5,
and 2.0 arcsec. Doyon et al. 2010 SPIE 7731. While these are not likely to see
extensive use with NIRISS, they are indeed still present in the hardware, so we
retain the ability to simulate them. 

NRM occulter mask: Provided by Anand Sivaramakrishnan. 



MIRI
------

MIRIM focal plane scale, 0.11 arcsec/pix:                 MIRI Optical Bench Assembly (OBA) Design Description, MIRI-DD-00001-AEU, 2.2.1

MIRIM field of view rotation, 4.561 degrees:              MIRI Optical Bench Assembly (OBA) Design Description, MIRI-DD-00001-AEU

Coronagraph pupils rotated to match,  4.56 degrees:  MIRI-DD-00001-AEU  5.7.8.2.1

Coronagraphic FOVs,  30.0 arcsec for Lyot, 24.0x23.8 arcsec for FQPMs: MIRI-DD-00001-AEU 2.2.1

Lyot coronagraph occulting spot diameter,               4.25 arcsec:      

Lyot coronagraph support bar width, 0.46 mm = 0.722 arcsec:              Anthony Boccaletti private communication December 2010 to Perrin and Hines

Lyot mask files:                                         Anthony Boccaletti private communication to Remi Soummer

LRS slit size (4.7 x 0.51 arcsec):     MIRI-TR-00001-CEA. And LRS Overview presentation by Silvia Scheithaur to MIRI team meeting May 2013. 

LRS P750L grating aperture mask (3.8% oversized tricontagon): MIRI OBA Design Description, MIRI-DD-00001-AEU




Instrument + Filter Throughputs
---------------------------------

Where possible, instrumental relative spectral responses were derived from the
Pysynphot CDBS files used for the development version of the JWST Exposure Time Calculators (ETCs),
normalized to peak transmission = 1.0 (because absolute throughput is not
relevant for PSF calculations). Not all filters are yet supported in Pysynphot,
however.  


For the following filters we take information from alternate sources other than the CDBS::

   Instrument    Filter         Source
   -----------  -------------   ----------------------------------------------------------------------------------------------------------
   NIRCam       F150W2          Top-hat function based on filter properties list at http://ircamera.as.arizona.edu/nircam/features.html
   NIRCam       F322W2          Top-hat function based on filter properties list at http://ircamera.as.arizona.edu/nircam/features.html
   NIRSpec      F115W           Assumed to be identical to the NIRCam one
   NIRSpec      F140X           NIRSpec "BBA" transmission curve traced from NIRSpec GWA FWA Assembly Report, NIRS-ZEO-RO-0051, section 6.3.2
   MIRI         F*W filters     Data published in Glasse et al. 2015 PASP Vol 127 No. 953, p. 688 Fig 2
   MIRI         F*C filters     Data published in Bouchet et al. 2015 PASP Vol 127 No. 953, p. 612 Fig 3
   NIRISS       all filters     Measurement data provided by Loic Albert of the NIRISS team
   FGS          none            Assumed top-hat function based on detector cut-on and cut-off wavelengths. 


The MIRI wide filters (F*W) are total system photon conversion efficiencies
including filter, telescope, instrument, and detector throughputs, normalized
to unity.  The MIRI coronagraphic filters are just the filters themselves, but
the detector and optics throughputs are relatively flat with wavelength
compared to the narrow coronagraphic filters. These are sufficiently accurate for
typical coronagraphic modeling but be aware of that caveat if attempting precise photometric 
calculations.

For the NIRCam and NIRSpec filters called out in the table above, the provided throughputs do not include the detector QE or OTE/SI optics throughputs versus wavelength. 

All other filters do include these effects, to the extent that they are accurately 
captured in the Calibration Database in support of the ETCs. 
