from astropy.io import fits
import numpy as N
import pylab as P
import logging

from .. import webbpsf_core


def validate_vs_jwpsf_nircam():

    models = [ ('NIRCam','F200W', 'f200w_perfect_offset', '/Users/mperrin/software/jwpsf_v3.0/data/NIRCam/OPD/perfect_opd.fits', 0.034,True),
            ('NIRCam','F200W', 'f200w_perfect', '/Users/mperrin/software/jwpsf_v3.0/data/NIRCam/OPD/perfect_opd.fits', 0.034,False),
            ('NIRCam','F200W', 'f200w', '/Users/mperrin/software/jwpsf_v3.0/data/NIRCam/OPD/nircam_obs_w_rsrv1.fits', 0.034,True),
                ('MIRI','F1000W', 'f1000w', '/Users/mperrin/software/jwpsf_v3.0/data/MIRI/OPD/MIRI_OPDisim1.fits', 0.11,True)]


    fig = P.figure(1, figsize=(13,8.5), dpi=80)
    oversamp=4
    for params in models:

        nc = webbpsf_core.Instrument(params[0])
        nc.filter = params[1]
        nc.pupilopd = params[3] #'/Users/mperrin/software/jwpsf_v3.0/data/NIRCam/OPD/nircam_obs_w_rsrv1.fits'
        nc.pixelscale = params[4] #0.034 # this is wrong, but compute this way to match JWPSF exactly
        if params[5]:
            # offset by half a pixel to match the JWPSF convention
            nc.options['source_offset_r'] = params[4]/2 * N.sqrt(2)/oversamp  # offset half a pixel each in X and Y
            nc.options['source_offset_theta'] = -45


        jw_fn = 'jwpsf_%s_%s.fits' % (params[0].lower(), params[2].lower())
        my_fn = 'test_vs_' + jw_fn

        if not os.path.exists( my_fn):
            my_psf = nc.calcPSF(my_fn, oversample=oversamp, fov_pixels=512./oversamp)
        else:
            my_psf = fits.open(my_fn)

        jw_psf = fits.open(jw_fn)
        jw_psf[0].header.update('PIXELSCL', jw_psf[0].header['CDELT1']*3600)


        P.clf()
        #P.subplots_adjust(top=0.95, bottom=0.05, left=0.01, right=0.99)
        P.subplot(231)
        titlestr = "%s %s, \n"%  (params[0], params[2])
        poppy.display_PSF(my_psf, title=titlestr+"computed with WebbPSF" , colorbar=False)
        P.subplot(232)
        poppy.display_PSF(jw_psf, title=titlestr+"computed with JWPSF" , colorbar=False)
        P.subplot(233)
        poppy.display_PSF_difference(my_psf,jw_psf, title=titlestr+'Difference Image', colorbar=False)

        imagecrop = 30*params[4]

        P.subplot(234)
        poppy.display_PSF(my_psf, title=titlestr+"computed with WebbPSF", colorbar=False, imagecrop=imagecrop)
        centroid = poppy.measure_centroid(my_psf)
        P.gca().set_xlabel("centroid = (%.3f,%.3f)" % centroid)

        P.subplot(235)
        poppy.display_PSF(jw_psf, title=titlestr+"computed with JWPSF", colorbar=False, imagecrop=imagecrop)
        centroid = poppy.measure_centroid(jw_psf)
        P.gca().set_xlabel("centroid = (%.3f,%.3f)" % centroid)

        P.subplot(236)
        poppy.display_PSF_difference(my_psf,jw_psf, title='Difference Image', colorbar=False, imagecrop=imagecrop)

        P.savefig("results_vs_jwpsf_%s_%s.pdf" % (params[0], params[2]))



        # oversampling = 4
        # wavelength range = 5 for JWPSF

        #stop()



def compare_revs_tv(nlambda=20, oversample=16, fov_arcsec=4):
    iname = 'NIRCam'
    inst = webbpsf_core.Instrument(iname)

    def ee_find_radius(ee_fn, value):
        "find the radius at which a given EE occurs. bruce force & crude "
        radius = N.arange(1000.)/1000.
        ee = ee_fn(radius)
        wmin = N.argmin(N.abs(ee-value))
        return radius[wmin]


    for filt in ['F070W', 'F200W']:
        inst.filter = filt
        base_opd = inst.pupilopd
        nopd = 3 # for testing
        fwhm_v = N.zeros(nopd)
        fwhm_t = N.zeros(nopd)
        eer_v = N.zeros((nopd,3))
        eer_t = N.zeros((nopd,3))

        for i in range(nopd):
            inst.pupilopd = (base_opd, i) # iterate through slices
            PV = webbpsf.calc_or_load_psf('test_%s_%s_revV_%02d.fits' % (iname, filt, i+1), inst, nlambda=nlambda, oversample=oversample, fov_arcsec=fov_arcsec)
            inst.pupilopd = '/Users/mperrin/software/webbpsf/data/OPD_RevT/nircam_obs_w_rsrv%d.fits' % (i+1)
            PT = webbpsf.calc_or_load_psf('test_%s_%s_revT_%02d.fits' % (iname, filt, i+1), inst, nlambda=nlambda, oversample=oversample, fov_arcsec=fov_arcsec)

            fwhm_v[i] = webbpsf.measure_fwhm(PV)
            fwhm_t[i] = webbpsf.measure_fwhm(PT)
            ee_fn_v = webbpsf.measure_EE(PV)
            ee_fn_t = webbpsf.measure_EE(PT)
            for j, val in enumerate([0.4, 0.5, 0.6]):
                eer_v[i,j] = ee_find_radius(ee_fn_v, val)
                eer_t[i,j] = ee_find_radius(ee_fn_t, val)

        mean_fwhm_v = fwhm_v.mean()
        mean_fwhm_t = fwhm_t.mean()
        mean_eer_v = eer_v.mean(axis=0)
        mean_eer_t = eer_t.mean(axis=0)
        print("Filter: "+filt)
        print("   Rev T:  FWHM=%.4f    EE(0.4, 0.5, 0.6) = %.3f, %0.3f, %.3f" % (mean_fwhm_t, mean_eer_t[0], mean_eer_t[1], mean_eer_t[2]))
        print("   Rev V:  FWHM=%.4f    EE(0.4, 0.5, 0.6) = %.3f, %0.3f, %.3f" % (mean_fwhm_v, mean_eer_v[0], mean_eer_v[1], mean_eer_v[2]))
    stop()


def compare_pupils_tv( oversample=8, vmax=1e-5, skipone=True):
    P.clf()
    inst = webbpsf.NIRCam()
    inst.pupilopd=None

    fov_arcsec = 10
    nlambda=30

    inst.pupil = 'tricontagon.fits'
    psf_tri = webbpsf.calc_or_load_psf('test_NIRCam_perfect_tricontagon_o%d.fits' % oversample, inst, nlambda=nlambda, oversample=oversample, fov_arcsec=fov_arcsec)


    if not skipone:
        ax = P.subplot(1,3,1)
        webbpsf.display_PSF(psf_tri, normalize='peak', colorbar=False, title='Tricontagon')


        for i, rev in enumerate(['T','V']):
            inst.pupil = "pupil_Rev%s.fits" % rev
            psf = webbpsf.calc_or_load_psf('test_NIRCam_perfect_rev%s_o%d.fits'  % (rev, oversample),inst,  nlambda=nlambda, oversample=oversample, fov_arcsec=fov_arcsec)
            P.subplot(1,3,i+2)
            webbpsf.display_PSF(psf, normalize='peak', colorbar = (rev =='V'), title='OTE Rev '+rev)


        stop()
        P.clf()

    psf_V = fits.open('test_NIRCam_perfect_rev%s_o%d.fits'  % ('V', oversample))
    psf_T = fits.open('test_NIRCam_perfect_rev%s_o%d.fits'  % ('T', oversample))
    P.subplot(221)
    webbpsf.display_PSF_difference(psf_V, psf_tri, vmax=vmax, title="Rev V - tricontagon")
    P.subplot(222)
    webbpsf.display_PSF_difference(psf_V, psf_T,vmax=vmax, title="Rev V - Rev T")


    ax3 = P.subplot(223)
    ax3.set_ylabel('Azimuthally averaged profile')
    ax3.set_xlabel('Separation (arcsec)')

    ax4 = P.subplot(224)
    ax4.set_ylabel('Fractional Encircled Energy')
    ax4.set_xlabel('Separation (arcsec)')


    for p, label in zip([psf_tri, psf_T, psf_V], ['Tri', "Rev T", 'Rev V']):
        rad, rp = webbpsf.radial_profile(p)
        ee_fn = webbpsf.measure_EE(p)

        ax3.plot(rad,rp/rp.max())
        ax4.plot(rad, ee_fn(rad), label=label)

        print webbpsf.measure_fwhm(p)
    ax4.legend(loc='lower right')
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    #ax3.axhline([psf_V[0].data.max()*0.5], ls=":")

    #ax3.set_xbound(0,4)
    ax3.set_xbound(0.01,4)
    ax4.set_xbound(0,4)
    P.draw()


    stop()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(name)-12s: %(levelname)-8s %(message)s',)

