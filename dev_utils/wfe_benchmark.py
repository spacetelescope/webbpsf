# Note: meant to just copy/paste into ipython command line

# Import WebbPSF and setup multiprocessing
import webbpsf
webbpsf.webbpsf_core.poppy.conf.use_fftw = False
webbpsf.webbpsf_core.poppy.conf.use_multiprocessing = True
ncores = 8
webbpsf.webbpsf_core.poppy.conf.n_processes = ncores 

inst = webbpsf.NIRCam()
inst.filter = "F430M"
inst.detector_position = (1024,1024)

# Baseline test: No SI WFE, no distortion
inst.include_si_wfe = False
%timeit psf = inst.calc_psf(add_distortion=False, monochromatic=4.3e-6)                                                                                
# Result: 911 ms ± 5.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
%timeit psf = inst.calc_psf(add_distortion=False, nlambda=ncores)                                                                                
# Result: 5.62 s ± 177 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

# Turn on SI WFE at center of detector (will use interpolation)
inst.include_si_wfe = True
%timeit psf = inst.calc_psf(add_distortion=False, monochromatic=4.3e-6)                                                                                
# Result: 1.41 s ± 12.6 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
%timeit psf = inst.calc_psf(add_distortion=False, nlambda=ncores)                                                                                
# Result: 6.1 s ± 96.9 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

# Use pixel (0,0) to force extrapolation algorithm
inst.detector_position = (0,0)
%timeit psf = inst.calc_psf(add_distortion=False, monochromatic=4.3e-6)                                                                                
# Result: 1.8 s ± 12.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
%timeit psf = inst.calc_psf(add_distortion=False, nlambda=ncores)                                                                                
# Result: 6.53 s ± 85.1 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)