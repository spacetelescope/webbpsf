.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Appendix: Instrument Property References
================================================================

We give here references for the instrumental properties assumed in PSF
computations, with particular attention to coronagraphic optics. It also notes
several places where the current models or available files are limited in some
manner that might be improved in a future release. 


OTE
----

The supplied OPDs are Revision V OPDs provided by Kong Ha at NASA GSFC computed using IPAM. The telescope pupil function was derived as a mask of all pixels which are non-zero in the OPD files. 

**Note:** The provided files included no header metadata, and in particular no pixel scale, so one was assumed based on the apparent pupil diameter in the files.

NIRCam
------

Coronagraph models based on Krist et al. 2007, 2009 SPIE papers, which gives the equations for sombrero and sinc band-limited masks.
These masks are specified in terms of a width parameter :math:`\sigma`. The equations are nonlinear and cannot be straightforwardly inverted
to solve for :math:`\sigma` for a given desired half-width at half-max (HWHM). ::


        Three circular occulters: HWHM = 6 lambda/D at 2.1, 3.35, 4.3
                                       = 0.4, 0.64, 0.8 arcsec (avg)
                                       assuming D_tel=6.5m exactly:
                                        = 0.3998, 0.6378, 0.8187 arcsec
        These values of HWHM are achieved using the following values for :math:`\sigma` in the sombrero BLC equation:
                                :math:`\sigma` = 5.2530, 3.2927, 2.5652



        Two linear bar occulters: Wedges vary from HWHM = 2 lam/D to 6 lam/D at 2.1 and 4.6 micron
                    2.1e-6:    HWHM = 0.13327 to 0.3998
                    4.6e-6:    HWHM = 0.27290 to 0.8187

        For these, I numerically solved the sinc BLC equation in inverse to find :math:`\sigma` as a function of HWHM. 
        I fit a fourth-degree polynomial to the results, and used this to smoothly interpolate to find :math:`\sigma` as 
        a function of position along the wedge. 


Note that the NIRCam wedge BLCs both have 'flat' regions with constant FWHM at the extreme left and right
sides of the wedge, as well as the region in the middle with varying FWHM. Though the widths of these flat 
regions are not explicitly stated in either of Krist's papers, by inspection of the figures they appear to be
~ 2.5 arcsec wide, so the actual wedge is 15 arcsec in length.  **Note:** This should be double-checked with John Krist.


NIRCam Lyot stop information was provided by John Krist in the file "JWST NIRCam Lyot Stop Definitions" dated January 22, 2007. The provided mask
data were in the form of pupil plane coordinates normalized by the telescope radius. A Python script was used to convert these coordinates into
pixel mask files 1024x1024 pixels in size. This transformation included a bit of anti-aliasing such that greyscale values are used for pixels right along the 
border of curved or diagonal edges.  However, this algorithm could probably be improved further.



MIRI
------

MIRIM focal plane scale, 0.11 arcsec/pix:                 MIRI Optical Bench Assembly (OBA) Design Description, MIRI-DD-00001-AEU, 2.2.1

MIRIM field of view rotation, 4.561 degrees:              MIRI Optical Bench Assembly (OBA) Design Description, MIRI-DD-00001-AEU

Pupil rotation,  4.56 degrees:  MIRI-DD-00001-AEU  5.7.8.2.1

Coronagraphic FOVs,  30.0 arcsec for Lyot, 24.0x23.8 arcsec for FQPMs: MIRI-DD-00001-AEU 2.2.1

Lyot coronagraph occulting spot diameter,               4.25 arcsec:      

Lyot coronagraph support bar width, 0.46 mm = 0.722 arcsec:              Anthony Boccaletti private communication December 2010


Lyot mask files:                                         Anthony Boccaletti private communication to Remi Soummer



TFI
----

TFI Etalon spectral resolution model:            From Craig Haley at ComDev, provided by Alex Fullerton

The transmission of TFI is modeled as a Gaussian with peak 1.0 and FWHM corresponding to the spectral resolution at the given wavelength. **Note:** In a future version of this software this should be improved to match the Airy function for an Etalon as given in "An Introduction to the TFI Etalon", JWST-STScI-002059.


TFI occulting spots: Assumed to be perfect circles with diameters 0.58, 0.75, 1.5, and 2.0 arcsec. 

Lyot occulter masks were provided by David Lafreniere. **Note:** The stated pixel scale is 6.44716 mm/pixel, which is slightly discrepant from the assumed pixel 
scale for the IPAM OPDs (differing by ~1 part in 1000). This discrepancy should be resolved in future versions of this software.


Instrument + Filter Throughputs
---------------------------------

Where possible, these were derived from the Pysynphot CDBS files used for the
JWST Exposure Time Calculators (ETCs), normalized to peak transmission = 1.0
(because absolute throughput is not relevant for PSF calculations). Not all
filters are yet supported in Pysynphot, however.::

   Instrument    Filter         Source
   -----------  --------        ----------------------------------------------------------------------------------------------------------
   NIRCam       F150W2          Top-hat function based on filter properties list at http://ircamera.as.arizona.edu/nircam/features.html
   NIRCam       F322W2          Top-hat function based on filter properties list at http://ircamera.as.arizona.edu/nircam/features.html
   MIRI         F1065C          MIRI test team spreadsheet provided to Christine Chen, obtained from STScI Coron WG site
   MIRI         F1140C          MIRI test team spreadsheet provided to Christine Chen, obtained from STScI Coron WG site
   MIRI         F1550C          MIRI test team spreadsheet provided to Christine Chen, obtained from STScI Coron WG site
   MIRI         F2300C          MIRI test team spreadsheet provided to Christine Chen, obtained from STScI Coron WG site
   MIRI         FND             MIRI test team spreadsheet provided to Christine Chen, obtained from STScI Coron WG site
   NIRSpec      F115W           Assumed to be identical to the NIRCam one
   NIRSpec      F140X           Top-hat function based on stated filter bandpass in NIRSpec Docs
   FGS          none            Assumed top-hat function based on detector cut-on and cut-off wavelengths. 


The above filters' throughputs do not include the detector QE or OTE/SI optics throughputs versus wavelength (or the throughput of the 
Germanium FQPM substrates for the MIRI coronagraphic filters). All other filters do include these effects to the extent that they are accurately 
captured in the Calibration Database in support of the ETCs. 


--------------

Documentation last updated on |today|

