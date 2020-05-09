Data files for OTE Linear Models

JWST_influence_functions_control_with_sm.fits:
	Coefficients from the WSS linear control model, derived by D. S. Acton et al. 

seg_sens.txt:
	Coefficients from the independent linear optical model by Erin Elliott. 


coarse_track_sim_pointing.fits    coarse_track2_sim_pointing.fits:
	Two models of JWST line of sight jitter in coarse pointing mode.
	Not strictly a part of the OTE linear model, but used with it to simulate
	unstacked single segment PSFs in some stages of commissioning.

See comment text at top of opds.py for more information.
