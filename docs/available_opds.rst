Appendix: Available Optical Path Difference (OPD) files
================================================================


On-Orbit OPD Measurements
-------------------------

The default OPD in WebbPSF is now an on-orbit measured OPD from early in cycle 1 science operations: 

 * `JWST_OTE_OPD_cycle1_example_2022-07-30.fits`: The specific OPD measurement selected, from 2022 July 30, intentionally represents a slightly conservative alignment state (after a couple weeks of relative stability, shortly before a next wavefront correction). This is about a 75th percentile  performance level, based on operational experience during the first few months of cycle 1. Most observations will have OPDs similar to this, or slightly better; a smaller number will have not quite as good. 

This intentionally-slightly-conservative OPD was chosen for use in the cycle 2 ETC and is set as the default here for consistency. 
For modeling PSFs of science observations, it is generally best to retrieve and use the 
`measured in-flight wavefronts <https://webbpsf.readthedocs.io/en/latest/jwst_measured_opds.html>`_ based on the date of 
the science observation. 




Pre-Launch Model OPDs (Obsolete)
---------------------------------

Prior versions of this page documented the _pre-flight_ models for JWST wavefronts. These models are now happily obsolete, and can be replaced by 
`actual measured in-flight wavefronts <https://webbpsf.readthedocs.io/en/latest/jwst_measured_opds.html>`_.
    

Most pre-launch predicted OPDs are no longer included in the WebbPSF data files distribution. The only such files still
included are

 * `JWST_OTE_OPD_RevAA_prelaunch_predicted.fits.gz`:  This contains predicted OPDs based on the best available optical models at the time of launch, including the calibrated, known surface figures of each mirror as measured in ground testing. For most purposes this is superceded by the on-orbit measurements. It is retained for back-compatibility mostly. 
    
    
For any cases in which you may need further information about the other pre-launch predicted OPD files, 
consult the `older versions of this page hosted on readthedocs<https://webbpsf.readthedocs.io/en/v0.9.0/available_opds.html>`_.
    
