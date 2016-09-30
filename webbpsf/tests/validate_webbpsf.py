from __future__ import division, print_function, absolute_import, unicode_literals
import os
from astropy.io import fits
import numpy as N
import matplotlib
import pylab as P
import logging
import poppy


from .. import webbpsf_core


__doc__="""

Validation Tests for Webb PSF. These functions perform simulations using WebbPSF and compare 
their results with the output of other simulations - e.g. earlier simulations by JWPSF, or ones
from the SI teams, etc.


"""


def validate_vs_russ_plot7(base_opd = 'OPD_RevV_nircam_155.fits'):
    """ Validate against plots from Makidon et al. 2007 JWST-STScI-001157 

    """

    waves_flux = N.array([ 0.6672437,  1.0602914,  1.458226 ,  1.9633503,  3.5440383, 4.7919989])
    flux_in_150_mas = N.array([ 0.7121211,  0.7727273,  0.7878786,  0.8484846,  0.6931819, 0.6893938])

    waves_fwhm = N.array([ 0.6672437,  1.0602914,  1.458226 ,  1.9550256,  3.5440383, 4.7919989]) 
    fwhm_arcsec = N.array([ 0.0225   ,  0.0347727,  0.0477273,  0.0627273,  0.1125   , 0.1520454])


    #fig, axarr = P.subplots(fignum=1, 
    P.clf()
    P.subplots_adjust(hspace=0.04)
    ax1 = P.subplot(211)
    ax1.plot(waves_flux, flux_in_150_mas, 'bx-', markersize=6, markeredgewidth=2)
    ax1.set_ybound(0, 1.0)
    ax1.set_xscale('log')

    ax2 = P.subplot(212, sharex=ax1)
    ax2.plot(waves_fwhm, fwhm_arcsec, 'bx-', markersize=6, markeredgewidth=2)
    ax2.set_ybound(0, 0.18)
    ax2.set_yticks([0.0, 0.05, 0.10, 0.15])
    ax2.set_xbound(0.6, 6)
    ax2.set_xlabel('Wavelength ($\mu$m)')
    ax2.set_ylabel('FWHM (arcsec)')
    ax1.set_ylabel('Flux within 0.15 arcsec')


    filt_list = ['F070W', 'F090W', 'F115W', 'F150W', 'F200W', 'F250M', 'F300M', 'F360M', 'F410M', 'F430M', 'F460M', 'F480M']

    nircam = webbpsf_core.NIRCam()
    waves = []
    fwhms = []
    ee15s = []
    nlambda=5   # empirically this gives essentially indistinguishable results as using the full thing.
    oversample=4 # empirically going higher than this to 8 makes no appreciable difference.
    fov = 8.0    # make it relatively large to get the denominator right for the EE.
    #base_opd = nircam.pupilopd

    for filt in filt_list:
        nircam.filter = filt
        psf = nircam.calcPSF( nlambda=nlambda, fov_arcsec=fov, oversample=oversample) # we only need a tiny PSF here
        waves.append(psf[0].header['WAVELEN'] * 1e6)
        #fwhms.append( webbpsf.measure_fwhm(psf) )
        #ee = webbpsf.measure_EE(psf)
        #ee15s.append(  ee(0.15) )

        # try other OPDs to check scatter
        my_ees = []
        my_fwhms = []
        for i in range(10): #[1,2,3,4,5]:
            nircam.pupilopd = (base_opd,i)
            psf = nircam.calcPSF( nlambda=nlambda, fov_arcsec=fov, oversample=oversample) # we only need a tiny PSF here
            ee = poppy.measure_EE(psf)
            my_ees.append(ee(0.15))
            my_fwhms.append(poppy.measure_fwhm(psf))
            ax1.plot([waves[-1]], [ee(0.15)], 'r+')
            ax2.plot([waves[-1]], [  my_fwhms[-1]], 'r+')
            P.draw()
        ee15s.append( N.array(my_ees).mean())
        fwhms.append( N.array(my_fwhms).mean())

        ax1.plot(waves, ee15s, 'ro-')
        ax2.plot(waves, fwhms, 'ro-')



        P.draw()
    l1 =ax2.plot(waves_fwhm, fwhm_arcsec, 'bx-', markersize=6, markeredgewidth=2, label='Makidon et al. 2007')
    l2 = ax2.plot(waves, fwhms, 'ro-',label='This work')

    ax1.set_ybound(0.6, 0.9)
    #ax1.set_ybound(0, 1.0)
    ax2.set_ybound(0, 0.18)
    ax2.set_yticks([0.0, 0.05, 0.10, 0.15])
    ax2.set_xbound(0.6, 6)
    ax2.set_xticks([0.7, 0.8, 0.9, 1.0, 2,3, 4, 5])
    ax2.set_xticklabels([0.7, 0.8, 0.9, 1.0, 2,3, 4, 5])
    ax1.set_xticklabels([])

    ax2.legend([l1, l2], ['Makidon et al. 2007', 'This work'], loc='lower right', frameon=False)
    P.draw()
 
    P.savefig('results_makidon_2007_fig7.pdf')
    stop()

    
def validate_vs_russ_plot6(nlambda=5, ax=None):
    """ Compare against the radial profiles in Plot 6. 
    
    """
    P.clf()
    P.subplots_adjust(wspace=0.01, hspace=0.04)
        #validate_vs_russ_plot6('F070W', nlambda=nlambda, ax=ax) #20)
        #ax2 = P.subplot(122)
        #validate_vs_russ_plot6('F200W', nlambda=nlambda, ax=ax2)
        #ax2.set_xlabel("Radial separation (arcsec)")
        #ax2.set_yticklabels([])

        #return
    for i, filter_ in enumerate(['F070W', 'F200W']):
        ax = P.subplot(2,2,i+1)
        ax.set_xlabel("Radial separation (arcsec)")
        if i == 0: ax.set_ylabel("Azimuthally averaged profile")

        if filter_ == 'F070W':
            # first a whole mess of data traced from those plots via Graphclick. Stuck here just to avoid having tons of random extra data files for this.
            rad_prof= N.array([ 0.0039118,  0.0061468,  0.0083818,  0.0106169,  0.0128519, 0.0162044,  0.0173219,  0.0206745,  0.021792 ,  0.0251445, 0.026262 ,  0.0329671,  0.0363196,  0.0396722,  0.0396722,
                0.0396722,  0.0396722,  0.0486123,  0.0486123,  0.0486123, 0.0575524,  0.0575524,  0.0575524,  0.0664925,  0.0664925, 0.0664925,  0.0754326,  0.0754326,  0.0754326,  0.0843727,
                0.0843727,  0.0843727,  0.0933128,  0.0933128,  0.0933128, 0.1022529,  0.1022529,  0.1022529,  0.111193 ,  0.111193 , 0.111193 ,  0.1201331,  0.1201331,  0.1201331,  0.1201331,
                0.1290732,  0.1290732,  0.1290732,  0.1380133,  0.1380133, 0.1380133,  0.1469534,  0.1469534,  0.1469534,  0.1558935, 0.1558935,  0.1558935,  0.1648336,  0.1648336,  0.1648336,
                0.1737737,  0.1737737,  0.1737737,  0.1827138,  0.1827138, 0.1827138,  0.1916539,  0.1916539,  0.1916539,  0.200594 , 0.200594 ,  0.200594 ,  0.200594 ,  0.2095341,  0.2095341,
                0.2095341,  0.2184742,  0.2184742,  0.2184742,  0.2274143, 0.2274143,  0.2274143,  0.2363545,  0.2363545,  0.2363545, 0.2452946,  0.2452946,  0.2452946,  0.2542346,  0.2542346,
                0.2542346,  0.2631748,  0.2631748,  0.2631748,  0.2721149, 0.2721149,  0.2721149,  0.281055 ,  0.281055 ,  0.281055 , 0.281055 ,  0.2899951,  0.2899951,  0.2899951,  0.2989352,
                0.2989352,  0.2989352,  0.3078753,  0.3078753,  0.3078753, 0.3168154,  0.3168154,  0.3168154,  0.3257555,  0.3257555, 0.3257555,  0.3346956,  0.3346956,  0.3346956,  0.3436357,
                0.3436357,  0.3436357,  0.3525758,  0.3525758,  0.3525758, 0.3615159,  0.3615159,  0.3615159,  0.3615159,  0.370456 , 0.370456 ,  0.370456 ,  0.3793961,  0.3793961,  0.3793961,
                0.3883362,  0.3883362,  0.3883362,  0.3972763,  0.3972763, 0.3972763,  0.4062164,  0.4062164,  0.4062164,  0.4151565, 0.4151565,  0.4151565,  0.4240966,  0.4240966,  0.4240966,
                0.4330367,  0.4330367,  0.4330367,  0.4419768,  0.4419768, 0.4419768,  0.4419768,  0.4509169,  0.4509169,  0.4509169, 0.459857 ,  0.459857 ,  0.459857 ,  0.4687971,  0.4687971,
                0.4687971,  0.4777372,  0.4777372,  0.4777372,  0.4866773, 0.4866773,  0.4866773,  0.4956174,  0.4956174,  0.4956174, 0.5045575,  0.5045575,  0.5045575,  0.5134977,  0.5134977,
                0.5134977,  0.5224378,  0.5224378,  0.5224378,  0.5224378, 0.5313779,  0.5313779,  0.5313779,  0.540318 ,  0.540318 , 0.540318 ,  0.5492581,  0.5492581,  0.5492581])
            prof = N.array([  5.07173700e-01,   4.60630100e-01,   3.96958800e-01, 3.21763100e-01,   2.41059000e-01,   1.88674300e-01, 1.47673200e-01,   1.28379400e-01,   1.04975200e-01,
                 8.73531000e-02,   7.46228000e-02,   7.01890000e-02, 5.99600000e-02,   5.39830000e-02,   5.39830000e-02, 5.39830000e-02,   5.39830000e-02,   3.67316000e-02,
                 3.67316000e-02,   3.67316000e-02,   2.41336000e-02, 2.41336000e-02,   2.41336000e-02,   1.53109000e-02, 1.53109000e-02,   1.53109000e-02,   1.04180000e-02,
                 1.04180000e-02,   1.04180000e-02,   7.87360000e-03, 7.87360000e-03,   7.87360000e-03,   6.16260000e-03, 6.16260000e-03,   6.16260000e-03,   5.35740000e-03,
                 5.35740000e-03,   5.35740000e-03,   4.65750000e-03, 4.65750000e-03,   4.65750000e-03,   3.90970000e-03, 3.90970000e-03,   3.90970000e-03,   3.90970000e-03,
                 3.28200000e-03,   3.28200000e-03,   3.28200000e-03, 2.75500000e-03,   2.75500000e-03,   2.75500000e-03, 2.31270000e-03,   2.31270000e-03,   2.31270000e-03,
                 2.23310000e-03,   2.23310000e-03,   2.23310000e-03, 2.15630000e-03,   2.15630000e-03,   2.15630000e-03, 2.23310000e-03,   2.23310000e-03,   2.23310000e-03,
                 2.08220000e-03,   2.08220000e-03,   2.08220000e-03, 2.01050000e-03,   2.01050000e-03,   2.01050000e-03, 1.87460000e-03,   1.87460000e-03,   1.87460000e-03,
                 1.87460000e-03,   1.68770000e-03,   1.68770000e-03, 1.68770000e-03,   1.51950000e-03,   1.51950000e-03, 1.51950000e-03,   1.41680000e-03,   1.41680000e-03,
                 1.41680000e-03,   1.23170000e-03,   1.23170000e-03, 1.23170000e-03,   1.10890000e-03,   1.10890000e-03, 1.10890000e-03,   9.64000000e-04,   9.64000000e-04,
                 9.64000000e-04,   8.38100000e-04,   8.38100000e-04, 8.38100000e-04,   7.54500000e-04,   7.54500000e-04, 7.54500000e-04,   6.55900000e-04,   6.55900000e-04,
                 6.55900000e-04,   6.55900000e-04,   6.11600000e-04, 6.11600000e-04,   6.11600000e-04,   5.31700000e-04, 5.31700000e-04,   5.31700000e-04,   4.78700000e-04,
                 4.78700000e-04,   4.78700000e-04,   4.46300000e-04, 4.46300000e-04,   4.46300000e-04,   4.16100000e-04, 4.16100000e-04,   4.16100000e-04,   3.88000000e-04,
                 3.88000000e-04,   3.88000000e-04,   3.61800000e-04, 3.61800000e-04,   3.61800000e-04,   3.37300000e-04, 3.37300000e-04,   3.37300000e-04,   3.14500000e-04,
                 3.14500000e-04,   3.14500000e-04,   3.14500000e-04, 2.93200000e-04,   2.93200000e-04,   2.93200000e-04, 2.73400000e-04,   2.73400000e-04,   2.73400000e-04,
                 2.54900000e-04,   2.54900000e-04,   2.54900000e-04, 2.54900000e-04,   2.54900000e-04,   2.54900000e-04, 2.37700000e-04,   2.37700000e-04,   2.37700000e-04,
                 2.21600000e-04,   2.21600000e-04,   2.21600000e-04, 2.14000000e-04,   2.14000000e-04,   2.14000000e-04, 2.06600000e-04,   2.06600000e-04,   2.06600000e-04,
                 1.92700000e-04,   1.92700000e-04,   1.92700000e-04, 1.92700000e-04,   1.86000000e-04,   1.86000000e-04, 1.86000000e-04,   1.79600000e-04,   1.79600000e-04,
                 1.79600000e-04,   1.73500000e-04,   1.73500000e-04, 1.73500000e-04,   1.67500000e-04,   1.67500000e-04, 1.67500000e-04,   1.56200000e-04,   1.56200000e-04,
                 1.56200000e-04,   1.50800000e-04,   1.50800000e-04, 1.50800000e-04,   1.45600000e-04,   1.45600000e-04, 1.45600000e-04,   1.45600000e-04,   1.45600000e-04,
                 1.45600000e-04,   1.40600000e-04,   1.40600000e-04, 1.40600000e-04,   1.40600000e-04,   1.40600000e-04, 1.40600000e-04,   1.40600000e-04,   1.35800000e-04,
                 1.35800000e-04,   1.35800000e-04,   1.26600000e-04, 1.26600000e-04,   1.26600000e-04])
            # fix the erroneous tracing of these points, which had the top Y axis marked at the 1.58 point instead of 1.0 due to how the plot axes are drawn.
            # in linear coordinates, the erroneous top of the axes was 1.04 higher than the real top.
            # So un-log scale and then re-log scale. 
            prof = 10**((N.log10(prof/1e-5))/5*1.04*5)*1e-5



            rad_perf = N.array([ 0.0184394,  0.0184394,  0.0206745,  0.021792 ,  0.0229095, 0.0251445,  0.026262 ,  0.0273795,  0.0307321,  0.0329671, 0.0340846,  0.0363196,  0.0407897,  0.0419072,  0.0452597,
                0.0452597,  0.0474948,  0.0486123,  0.0497298,  0.0508473, 0.0530823,  0.0541998,  0.0553174,  0.0586699,  0.0609049, 0.0642575,  0.06761  ,  0.069845 ,  0.0720801,  0.0743151,
                0.0765501,  0.0799026,  0.0810202,  0.0832552,  0.0843727, 0.0854902,  0.0877252,  0.0910778,  0.0944303,  0.0955478, 0.0989004,  0.1011354,  0.1044879,  0.106723 ,  0.108958 ,
                0.1123105,  0.1167806,  0.1178981,  0.1201331,  0.1212506, 0.1246032,  0.1257207,  0.1301907,  0.1346608,  0.1368958, 0.1391308,  0.1469534,  0.150306 ,  0.1536585,  0.1558935,
                0.1614811,  0.1670687,  0.1693037,  0.1737737,  0.1894189, 0.193889 ,  0.196124 ,  0.2061816,  0.2117692,  0.2128867, 0.2195918,  0.2240618,  0.2251793,  0.2274143,  0.2307669,
                0.2318844,  0.2341194,  0.237472 ,  0.2408245,  0.244177 , 0.2531171,  0.2709973,  0.2765849,  0.281055 ,  0.2944651, 0.2978176,  0.3112278,  0.3123453,  0.3156979,  0.3324606,
                0.3358131,  0.3436357,  0.3469882,  0.3492233,  0.3525758, 0.3581634,  0.3592809,  0.3603984,  0.3637509,  0.3738085, 0.3771611,  0.3793961,  0.3827486,  0.3939238,  0.3983938,
                0.4017463,  0.4050989,  0.416274 ,  0.4207441,  0.4274491, 0.4464469,  0.4520344,  0.455387 ,  0.4676796,  0.4699146, 0.4732672,  0.4822073,  0.4866773,  0.4922649,  0.4967349,
                0.4967349,  0.50344  ,  0.5079101,  0.5112626,  0.5146151, 0.5190852,  0.5224378,  0.5291429,  0.5347304,  0.5358479, 0.538083 ])
            prof_perf = N.array([  9.87379000e-02,   7.79603000e-02,   4.90290000e-02, 3.87117000e-02,   2.45597000e-02,   1.55812000e-02, 1.95620000e-02,   1.27407000e-02,   1.61363000e-02,
                 2.02588000e-02,   2.41336000e-02,   2.23058000e-02, 1.82393000e-02,   1.49142000e-02,   1.40281000e-02, 1.14707000e-02,   9.05690000e-03,   7.21390000e-03,
                 5.74590000e-03,   4.57670000e-03,   4.38070000e-03, 3.03340000e-03,   2.33300000e-03,   1.82600000e-03, 1.42920000e-03,   1.22090000e-03,   1.12850000e-03,
                 1.17890000e-03,   1.23170000e-03,   1.53280000e-03, 1.90770000e-03,   1.99300000e-03,   2.17530000e-03, 2.31270000e-03,   2.43740000e-03,   2.43740000e-03,
                 2.33300000e-03,   2.04600000e-03,   1.92450000e-03, 1.68770000e-03,   1.32100000e-03,   1.21030000e-03, 1.04300000e-03,   9.72500000e-04,   8.09200000e-04,
                 6.67500000e-04,   6.50200000e-04,   6.11600000e-04, 5.36400000e-04,   5.00100000e-04,   4.87100000e-04, 4.31000000e-04,   3.40300000e-04,   2.68700000e-04,
                 2.52700000e-04,   2.66300000e-04,   2.95800000e-04, 3.01000000e-04,   3.01000000e-04,   2.88200000e-04, 2.52700000e-04,   2.50500000e-04,   2.44000000e-04,
                 2.44000000e-04,   2.33600000e-04,   2.25500000e-04, 2.15900000e-04,   1.91000000e-04,   1.84400000e-04, 1.76500000e-04,   1.57500000e-04,   1.54800000e-04,
                 1.48200000e-04,   1.43100000e-04,   1.27700000e-04, 1.21200000e-04,   1.17000000e-04,   1.00800000e-04, 9.23792550e-05,   8.03097440e-05,   7.68716620e-05,
                 7.68716620e-05,   7.48798850e-05,   7.29397060e-05, 7.23042150e-05,   6.92088410e-05,   5.91227580e-05, 6.17670150e-05,   5.91227580e-05,   5.23057220e-05,
                 4.87692550e-05,   4.27701000e-05,   4.16619060e-05, 4.05824230e-05,   3.95309070e-05,   3.55903360e-05, 3.43661240e-05,   3.26083060e-05,   3.17634030e-05,
                 3.14866660e-05,   3.28949020e-05,   3.43661240e-05, 3.49728750e-05,   3.65370320e-05,   3.62187030e-05, 3.55903360e-05,   3.49728750e-05,   3.12123370e-05,
                 3.09403990e-05,   3.09403990e-05,   3.06708290e-05, 2.96158310e-05,   2.85971280e-05,   2.62010380e-05, 2.52997920e-05,   2.44295450e-05,   2.14244460e-05,
                 2.21876410e-05,   2.25793740e-05,   1.98018660e-05, 1.84630300e-05,   1.56349180e-05,   1.59109580e-05, 1.67686740e-05,   1.86253040e-05,   1.86253040e-05,
                 1.75186500e-05,   1.54986980e-05,   1.52298080e-05, 1.56349180e-05,   1.63341880e-05]) 
            prof_perf = 10**((N.log10(prof_perf/1e-5))/5*1.04*5)*1e-5


            rad_ee = N.array([ 0.       ,  0.0022541,  0.0045082,  0.0067623,  0.0090164, 0.0112705,  0.0135246,  0.0157787,  0.0180328,  0.0202869, 0.022541 ,  0.0270492,  0.0293033,  0.0338115,  0.0405738,
                    0.045082 ,  0.0495902,  0.0563525,  0.0608607,  0.067623 , 0.0766393,  0.0856557,  0.1014344,  0.1127049,  0.1172131, 0.1239754,  0.1307377,  0.1487705,  0.1622951,  0.1803279,
                    0.1893443,  0.2028689,  0.2209016,  0.2434426,  0.2502049, 0.2704918,  0.2862705,  0.2997951,  0.3155738,  0.3336066, 0.3426229,  0.3606558,  0.3786885,  0.3877049,  0.4057377,
                    0.4237705,  0.4327869,  0.4508197,  0.4688525,  0.4778689, 0.4959016,  0.5139344,  0.5297131,  0.5454918])
            ee = N.array([ 0.0936175,  0.1361595,  0.1787016,  0.2212436,  0.2485921, 0.2698631,  0.2926535,  0.3124052,  0.3321568,  0.3519085, 0.3640634,  0.3898925,  0.4142022,  0.4537056,  0.4947283,
                    0.5311929,  0.5524639,  0.5767736,  0.5919672,  0.6041221, 0.6208351,  0.637548 ,  0.6664158,  0.6831288,  0.6876869, 0.6998417,  0.7074385,  0.7256708,  0.7378257,  0.756058 ,
                    0.7697322,  0.7864452,  0.8016388,  0.8198711,  0.8259485, 0.8381034,  0.8472196,  0.8548164,  0.8593745,  0.8654519, 0.8684906,  0.874568 ,  0.8806455,  0.8806455,  0.8867229,
                    0.8897616,  0.8928003,  0.8958391,  0.9019165,  0.9019165, 0.907994 ,  0.9110327,  0.9140714,  0.9186295])

#            rad_ee = N.array([0.00225,  0.0045, 0.0124304,  0.0213705,  0.0303106,  0.0392507,  0.0392857, 0.0481908,  0.052381 ,  0.066071 ,  0.0750111,  0.0839513, 0.0928914,  0.1018315,  0.1107716,  0.1197117,  0.1286518,
#                    0.1375919,  0.146532 ,  0.1554721,  0.1644122,  0.1702381, 0.1822924,  0.1912325,  0.2001726,  0.2091127,  0.2180528, 0.2269929,  0.235933 ,  0.2448731,  0.2538132,  0.2627533,
#                    0.2716934,  0.2806336,  0.2985137,  0.3074538,  0.3163939, 0.3253341,  0.3342742,  0.3432143,  0.3521544,  0.3610945, 0.3700346,  0.3789747,  0.3879148,  0.3968549,  0.4147351,
#                    0.4236752,  0.4326153,  0.4415554,  0.4504955,  0.4594356, 0.4683757,  0.4773158,  0.4862559,  0.495196 ,  0.5041361, 0.5130762,  0.5238097,  0.5369048,  0.545635 ])
#            ee = N.array([ 0.1346, 0.1984, 0.2303598,  0.2941729,  0.3762183,  0.4400313,  0.4747775, 0.5008057,  0.545994 ,  0.5798124,  0.6041222,  0.6253932, 0.6405867,  0.6557804,  0.670974 ,  0.6831288,  0.6952837,
#                    0.7043998,  0.7165547,  0.7287096,  0.734787 ,  0.7448071, 0.7530193,  0.7651742,  0.7773291,  0.7864453,  0.7955614, 0.8046775,  0.810755 ,  0.8168324,  0.8259486,  0.832026 ,
#                    0.8381035,  0.8411421,  0.8502584,  0.8563358,  0.8563358, 0.8624132,  0.8654519,  0.8684906,  0.8715293,  0.874568 , 0.874568 ,  0.8806455,  0.8806455,  0.8867229,  0.8867229,
#                    0.8928003,  0.8928003,  0.8928003,  0.8988778,  0.8988778, 0.8988778,  0.9049552,  0.9049552,  0.9049552,  0.9110327, 0.9110327,  0.9139466,  0.9169139,  0.9169139])

            rad_ee_perf = N.array([ 0.0034903,  0.0068429,  0.0090779,  0.015783 ,  0.0191355, 0.022488 ,  0.0314281,  0.0358982,  0.0403682,  0.0470733,
                    0.0560134,  0.063836 ,  0.0772462,  0.0817162,  0.0895388, 0.0906563,  0.1040665,  0.1130066,  0.1286518,  0.1398269,
                    0.1554721,  0.1666472,  0.1834099,  0.2135828,  0.2415206, 0.2705759,  0.3007488,  0.3309216,  0.359977 ,  0.3879148,
                    0.4169701,  0.444908 ,  0.4761983,  0.5052537,  0.5331915, 0.5477192])
            ee_perf = N.array([ 0.3169633,  0.457504 ,  0.5577817,  0.6360288,  0.6785708, 0.6983225,  0.737066 ,  0.7773292,  0.8191115,  0.8373438,
                    0.8502584,  0.853297 ,  0.8616535,  0.87001  ,  0.8814053, 0.8859632,  0.9003971,  0.9026762,  0.9095133,  0.912552 ,
                    0.9178698,  0.9216682,  0.9262264,  0.9353424,  0.9406602, 0.9452182,  0.9482571,  0.9520553,  0.9550941,  0.9566135,
                    0.9596521,  0.9619312,  0.9657295,  0.967249 ,  0.9687683, 0.9687683])

        else:

            rad_prof = N.array([ 0.0022541,  0.0067623,  0.0112705,  0.0157787,  0.0202869, 0.0247951,  0.0293033,  0.0315574,  0.0360656,  0.0371926, 0.0439549,  0.047336 ,  0.0495901,  0.0540983,  0.0563524,
                0.0619877,  0.0631147,  0.0631147,  0.0653688,  0.0664959, 0.06875  ,  0.0721311,  0.0743852,  0.0766393,  0.0788934, 0.0811475,  0.0856557,  0.0901639,  0.092418 ,  0.0980532,
                0.1003073,  0.1048155,  0.1093237,  0.1127049,  0.1205942, 0.1307376,  0.139754 ,  0.1442622,  0.1532786,  0.162295 , 0.1713114,  0.1803278,  0.1882171,  0.1893442,  0.1949794,
                0.1983605,  0.2017417,  0.2107581,  0.2163933,  0.2197745, 0.2254097,  0.2344261,  0.2434425,  0.2524589,  0.2592212, 0.2648564,  0.2727458,  0.277254 ,  0.2828892,  0.2885244,
                0.2919056,  0.300922 ,  0.3043031,  0.3099384,  0.3167007, 0.3212089,  0.3279711,  0.3358605,  0.3392417,  0.3460039, 0.347131 ,  0.3561474,  0.3640367,  0.3753072,  0.3820695,
                0.3899588,  0.3967211,  0.4046105,  0.4147539,  0.4237703, 0.4282785,  0.4418031,  0.4463112,  0.4519465,  0.4587088, 0.4643441,  0.4733604,  0.4778686,  0.4947743,  0.5082989,
                0.5150612,  0.5229506,  0.5285859,  0.5353481,  0.5409833, 0.5466186])
            prof = N.array([  6.14843700e-01,   5.83394500e-01,   5.43950300e-01, 4.89727500e-01,   4.40909000e-01,   3.90071900e-01, 3.30321800e-01,   2.94803100e-01,   2.34813000e-01,
                 2.04134300e-01,   1.62594100e-01,   1.29507700e-01, 9.96057000e-02,   6.37477000e-02,   4.98946000e-02, 2.68056000e-02,   2.17278000e-02,   2.17278000e-02,
                 1.93915000e-02,   1.70061000e-02,   1.45278000e-02, 1.33105000e-02,   1.10761000e-02,   1.00596000e-02, 1.15715000e-02,   1.33105000e-02,   1.55812000e-02,
                 1.82393000e-02,   2.06165000e-02,   2.33034000e-02, 2.47755000e-02,   2.52129000e-02,   2.47755000e-02, 2.39233000e-02,   2.02588000e-02,   1.42757000e-02,
                 9.71360000e-03,   6.38210000e-03,   4.34260000e-03, 3.39890000e-03,   2.75500000e-03,   2.15630000e-03, 1.92450000e-03,   1.74790000e-03,   1.57360000e-03,
                 1.41680000e-03,   1.35610000e-03,   1.34430000e-03, 1.46720000e-03,   1.62970000e-03,   1.74790000e-03, 2.08220000e-03,   2.15630000e-03,   2.15630000e-03,
                 2.08220000e-03,   1.94140000e-03,   1.74790000e-03, 1.64400000e-03,   1.51950000e-03,   1.39220000e-03, 1.28670000e-03,   1.10890000e-03,   1.04300000e-03,
                 9.72500000e-04,   8.75500000e-04,   8.09200000e-04, 7.54500000e-04,   7.03500000e-04,   6.79300000e-04, 6.55900000e-04,   6.38900000e-04,   5.90600000e-04,
                 5.13400000e-04,   4.46300000e-04,   3.88000000e-04, 3.58600000e-04,   3.43300000e-04,   3.46300000e-04, 3.68200000e-04,   3.88000000e-04,   3.98300000e-04,
                 3.98300000e-04,   3.91400000e-04,   3.74700000e-04, 3.68200000e-04,   3.55500000e-04,   3.46300000e-04, 3.40300000e-04,   3.37300000e-04,   3.28600000e-04,
                 3.22900000e-04,   3.14500000e-04,   3.06400000e-04, 2.98400000e-04,   2.90700000e-04,   2.88200000e-04]) 
            # fix the erroneous tracing of these points, which had the top Y axis marked at the 1.58 point instead of 1.0 due to how the plot axes are drawn.
            # in linear coordinates, the erroneous top of the axes was 1.04 higher than the real top.
            # So un-log scale and then re-log scale. 
            prof = 10**((N.log10(prof/1e-5))/5*1.04*5)*1e-5


            prof_perf = None

            rad_ee = N.array([ 0.       ,  0.       ,  0.0056352,  0.0123975,  0.0202869, 0.0270492,  0.0304303,  0.0360656,  0.045082 ,  0.045082 ,
                    0.0495902,  0.0518443,  0.0597336,  0.0653689,  0.0732582, 0.0788934,  0.0856557,  0.0935451,  0.1104508,  0.1183402,
                    0.136373 ,  0.1442623,  0.1566598,  0.1713115,  0.1803279, 0.2006148,  0.2163935,  0.2434426,  0.2558402,  0.2761271,
                    0.2919057,  0.3065574,  0.3347336,  0.3550205,  0.3696721, 0.3933402,  0.4418033,  0.4632172,  0.480123 ,  0.504918 ,
                    0.527459 ,  0.5466189,  0.5680328])
            ee = N.array([ 0.0601916,  0.0601916,  0.0913384,  0.163508 ,  0.2911341, 0.3784973,  0.3974893,  0.5038444,  0.5372703,  0.5737349,
                    0.5866495,  0.6056415,  0.6269125,  0.6383077,  0.6436254, 0.6497028,  0.658819 ,  0.6793303,  0.7317483,  0.7568177,
                    0.7887242,  0.800879 ,  0.8099952,  0.8191115,  0.8229098, 0.829747 ,  0.8358244,  0.8502583,  0.8548164,  0.8654519,
                    0.872289 ,  0.874568 ,  0.8844439,  0.8859633,  0.8897616, 0.8950794,  0.9019165,  0.9034359,  0.9057149,  0.9110327,
                    0.914831 ,  0.9163504,  0.7499806]) 

            rad_ee_perf = N.array([ 0.0030056,  0.0180328,  0.0293033,  0.0338115,  0.0360656, 0.0383197,  0.0428279,  0.0484631,  0.0518443,  0.0586066,
                    0.0732582,  0.0822746,  0.0856557,  0.091291 ,  0.0957992, 0.1025615,  0.1093238,  0.114959 ,  0.1228484,  0.1318648,
                    0.1375   ,  0.1498975,  0.1566598,  0.1724385,  0.1848361, 0.1994877,  0.2130123,  0.2276639,  0.2366803,  0.251332 ,
                    0.2535861,  0.2580943,  0.2761271,  0.2817623,  0.298668 , 0.3088115,  0.3245902,  0.3369877,  0.3561476,  0.367418 ,
                    0.3843238,  0.3978484,  0.4215164,  0.4485656,  0.480123 , 0.5094262,  0.5364754,  0.5477459])
            ee_perf = N.array([ 0.0612045,  0.2961987,  0.436233 ,  0.4962476,  0.535751 , 0.5767736,  0.6170366,  0.65654  ,  0.6770514,  0.699082 ,
                    0.7081982,  0.7127563,  0.7241515,  0.7302289,  0.7446628, 0.7682129,  0.7841662,  0.8061969,  0.8251889,  0.8396228,
                    0.8464599,  0.8570954,  0.8608938,  0.8646922,  0.8677309, 0.8715293,  0.8738084,  0.8791261,  0.8844439,  0.8928003,
                    0.8965988,  0.8981181,  0.9064746,  0.910273 ,  0.914831 , 0.9171101,  0.9216682,  0.9254666,  0.9277456,  0.9277456,
                    0.931544 ,  0.931544 ,  0.9345827,  0.9361021,  0.9391408, 0.9429392,  0.9459779,  0.9459779])

            #wave_perf = wave
            #prof_perf = prof*0 # didn't bother tracing this one...
             
        ax = P.gca() 
        P.semilogy(rad_prof, prof, 'b-', label='Makidon et al. 2007')
        if prof_perf is not None:
            P.semilogy(rad_perf, prof_perf, 'b--',  label='Perfect PSF')
        P.title(filter_)

        ax.set_ybound(1e-5, 1)
        ax.set_xbound(0, 0.55)
        P.draw()

        oversample=10 # by testing, this needs to be really high. This is because it's normalized to the very peak central pixel...
        fov = 4.0    # make it relatively large to get the denominator right for the EE.
     
        nc = webbpsf_core.NIRCam()
        nc.pupilopd='OPD_RevV_nircam_123.fits'
        nc.filter = filter_
        psf = nc.calcPSF(oversample=oversample, nlambda=nlambda, fov_arcsec=fov)
        my_rad, my_prof = poppy.radial_profile(psf)
        ax.plot(my_rad, my_prof/my_prof.max(),'r-', label = 'This work')
        ax.set_xbound(0, 0.55)
        P.draw()

        nc.pupilopd=None
        psf_perf = nc.calcPSF(oversample=oversample, nlambda=nlambda, fov_arcsec=fov)
        my_rad_perf, my_prof_perf = poppy.radial_profile(psf_perf)
        ax.plot(my_rad_perf, my_prof_perf/my_prof_perf.max(),'r--', label = 'Perfect PSF')

        ax.set_ybound(1e-5, 1)
        ax.set_xbound(0, 0.55)
        ax.legend( loc='upper right', frameon=False)

        if i == 1: ax.set_yticklabels([])
        ax.set_xticklabels([])

        ########## bottom plot
        # we can just compute the encircled energy from the radial profiles...

        ax2 = P.subplot(2,2,i+3)
        if i == 0: ax2.set_ylabel("Fractional Encircled Energy")

#        def prof_to_ee(radius, profile):
#            dr = radius[1:] - radius[:-1]
#            flux = profile[:-1] * 2 * N.pi * radius[:-1]  * dr
#            ee = N.cumsum(flux)
#            ee /= ee.max() # normalize to the max of whatever we've got..
#            return (radius[:-1], ee)
#
#        rad_ee, ee = prof_to_ee(rad_prof, prof)
        ax2.plot(rad_ee, ee, 'b-', label='Makidon et al. 2007')
        ax2.plot(rad_ee_perf, ee_perf, 'b--', label='Perfect PSF')

        ee_fn = poppy.measure_EE(psf)
        my_rad_ee = N.linspace(0.0, 0.55, 100)
        ax2.plot(my_rad_ee, ee_fn(my_rad_ee), 'r-', label='This work')

        ee_fn_perf = poppy.measure_EE(psf_perf)
        ax2.plot(my_rad_ee, ee_fn_perf(my_rad_ee), 'r--', label='This work')

        ax2.set_xlabel("Radial separation (arcsec)")
        if i == 1: ax2.set_yticklabels([])
        ax2.set_xbound(0, 0.55)
        
        
        #stop()






        P.draw()
        #stop()



################################################################################
#
#       Comparisons with simulations by SI Teams
#

#---- NIRCam ------------------

def validate_vs_krist_blc(which='spot'):
    """ Compare the BLC occulter structure (not the PSFs, the actual occulter profile)
    with the data files provided by John Krist
    """
    if which=='spot':
        image_mask =  'MASK430R'
        pupil_mask = 'CIRCLYOT'
    else:
        image_mask =  'MASKLWB'
        pupil_mask = 'WEDGELYOT'
 
    nc = webbpsf_core.NIRCam()
    nc.pupilopd=None
    nc.filter = 'F460M'
    nc.image_mask = image_mask
    nc.pupil_mask = pupil_mask
    nc.options['no_sam'] = True
    nc.pixelscale = 0.065 # match the Krist sims exacly. vs 0.648 official

    cor_vmin = 1e-12
    cor_vmax=1e-5


    P.clf()
    mask1 = 'nircam_4600nm_%s_occ.fits' % wedge
    mask1f = fits.open(mask1)

    #P.subplot(332)
    os = nc._getOpticalSystem()
    #os.planes[1].display(ax=P.gca(),  what='intensity', colorbar_orientation='vertical')
    #P.gca().set_title('')
    #P.gca().set_xbound(-8,8)
    #P.gca().set_ybound(-8,8)

    wf = poppy.Wavefront(wavelength=4.6e-6,  npix = mask1f[0].data.shape[0], pixelscale = mask1f[0].header['PIXSIZE'])
    trans = os.planes[1].getPhasor(wf)**2

    npix = mask1f[0].data.shape[0]
    pixels = (N.arange(npix) - npix/2.) *  mask1f[0].header['PIXSIZE']

    ax = P.subplot(211)
    P.plot(pixels, mask1f[0].data[:, npix/2], label='Krist')
    P.plot(pixels, trans[:, npix/2], ls="--", color='r', label="Perrin")
    ax.set_ylabel("Transmission")

    ax.legend(loc='lower right')

    ax2 = P.subplot(212,sharex=ax)

    P.semilogy(pixels, 1-mask1f[0].data[:, npix/2], label='Krist')
    P.semilogy(pixels, 1-trans[:, npix/2],  ls="--", color='r', label="WebbPSF")
    P.semilogy(pixels, trans[:, npix/2] - mask1f[0].data[:, npix/2] ,  ls=":", color='r', label="Difference")
    ax2.set_ylabel("1 - Transmission")
    ax2.set_xlabel("Radius [arcsec]")
    ax.set_xbound(-8,8)
    #ax2.set_xbound(-7.3,-6.6)
    #ax2.set_ybound(1e-4, 1e-3)
    P.draw()

    P.savefig('results_nircam_blc430r_profile_comparison.pdf')


    diff = trans[npix/2, :] - mask1f[0].data[npix/2, :] 
    print("Max diff: %.3g" % diff.max())



def validate_vs_krist_sims(clobber=False, normalize=False, which='spot', no_sam=False):
    """ Compare with PSFs provided by John Krist 
    """


    if which=='spot':
        image_mask =  'MASK430R'
        pupil_mask = 'CIRCLYOT'
    else:
        image_mask =  'MASKLWB'
        pupil_mask = 'WEDGELYOT'
 
    P.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.05)
    
    nc = webbpsf_core.NIRCam()
    nc.pupilopd=None
    nc.filter = 'F460M'
    nc.image_mask = image_mask 
    nc.pupil_mask = pupil_mask
    nc.options['no_sam'] = no_sam
    nc.pixelscale = 0.065 # match the Krist sims exacly. vs 0.648 official

    cor_vmin = 1e-12
    cor_vmax=1e-5


    P.clf()

    fig = P.gcf()
    fig.text(0.2, 0.95, 'Krist', horizontalalignment='center', size=18)
    fig.text(0.50, 0.95, 'Perrin', horizontalalignment='center', size=18)
    fig.text(0.80, 0.95, 'Difference P-K', horizontalalignment='center', size=18)

    fig.text(0.05, 1./6, 'off-axis 4.6$\mu$m', verticalalignment='center', rotation='vertical' , size=18)
    fig.text(0.05, 0.48, 'occulted 4.6$\mu$m', verticalalignment='center', rotation='vertical' , size=18)
    fig.text(0.05, 5./6-0.05, image_mask + ' occulter', verticalalignment='center', rotation='vertical' , size=18)


    P.subplot(331)
    mask1 = 'nircam_4600nm_%s_occ.fits' % which
    mask1f = fits.open(mask1)
    poppy.display_PSF(mask1f, title="", pixelscale='PIXSIZE', vmin=0, vmax=1, scale='linear', cmap=matplotlib.cm.gray)

    P.subplot(332)
    os = nc._getOpticalSystem()
    os.planes[1].display(ax=P.gca(),  what='intensity', colorbar_orientation='vertical')
    P.gca().set_title('')
    P.gca().set_xbound(-8,8)
    P.gca().set_ybound(-8,8)

    wf = poppy.Wavefront(wavelength=4.6e-6,  npix = mask1f[0].data.shape[0], pixelscale = mask1f[0].header['PIXSIZE'])
    trans = os.planes[1].getPhasor(wf)**2

    P.subplot(333)

    if normalize: 
        to_plot = (trans-mask1f[0].data) / (trans+mask1f[0].data)/2
        vmin, vmax = -1, 1
    else: 
        to_plot = (trans-mask1f[0].data) 
        vmin, vmax = -1e-3, 1e-3
 
    poppy.imshow_with_mouseover(to_plot, cmap=matplotlib.cm.gray, vmin=vmin, vmax=vmax, extent=[-8, 8, -8, 8])
    P.colorbar(P.gca().images[0], orientation='vertical')
    try:
        fits.PrimaryHDU(trans).writeto('test_nircam_4600nm_%s_occ.fits' % which, clobber=clobber)
    except:
        pass





    #---- occulted --
    P.subplot(334)
    k1 = 'nircam_4600nm_%s.fits' % which
    k1f = fits.open(k1)

    print("Total of %s is %f" % (k1, k1f[0].data.sum()))
    poppy.display_PSF(k1f,  title="", pixelscale='SAMPLING', vmin=cor_vmin, vmax=cor_vmax)

    P.subplot(335)
    my1 = 'test_'+k1
    mypsf1 = webbpsf_core.calc_or_load_psf('test_'+k1, nc, nlambda=1,monochromatic=4.6e-6, oversample=4, fov_pixels=247, clobber=clobber)
    print("Total of %s is %f" % (my1, mypsf1[0].data.sum()))
    #nc.calcPSF(nlambda=1)
    poppy.display_PSF(mypsf1, ext=1, title="", adjust_for_oversampling=True, vmin=cor_vmin, vmax=cor_vmax)

    P.subplot(336)
    poppy.display_PSF_difference( mypsf1, k1f,  ext2=0, ext1=1, title="", vmax=1e-7, normalize=normalize)


    #---- unocculted --
    P.subplot(337)
    k2 = 'nircam_4600nm_%s_fieldpsf.fits' % which
    k2f = fits.open(k2)
    poppy.display_PSF(k2f,  title="", pixelscale='SAMPLING')
    print("Total of %s is %f" % (k2, k2f[0].data.sum()))


    nc.image_mask = None # make a coronagraphic-off-axis type PSF but still on-axis in the array
    my2 = 'test_'+k2
    mypsf2 = webbpsf_core.calc_or_load_psf('test_'+k2, nc, nlambda=1, monochromatic=4.6e-6, oversample=4, fov_pixels=247, clobber=clobber)
    P.subplot(338)
    poppy.display_PSF(mypsf2,  title="", ext=1, adjust_for_oversampling=True)
    print("Total of %s is %f" % (my2, mypsf2[0].data.sum()))
 
    P.subplot(339)
    poppy.display_PSF_difference( mypsf2, k2f,  ext2=0, ext1=1, title="", vmax=1e-5, normalize=normalize)



    print("shape of %s is %s" % (k1, k1f[0].data.shape))
    print("shape of %s is %s" % (my1, mypsf1[1].data.shape))

    P.savefig('results_nircam_coron_comparison_%s.pdf' % which)
    stop()
     

#---- MIRI ---------------------

def validate_vs_cavarroc_2008():
    inst = 'MIRI'
    filter = 'F1065C'

    miri = webbpsf_core.MIRI()
    miri.filter = filter
    miri.image_mask = 'FQPM1065'
    miri.pupil_mask = 'MASKFQPM'

    miri.options['source_offset_r'] = 0.005
    miri.options['source_offset_theta'] = 45

    miri.calcPSF('test_%s_%s_offset%.3f.fits'  % (inst, filter,  miri.options['source_offset_r']))



    #--- Acq image 


def validate_miri_coron():
    """Validate MIRI coronagraphic performance against the simulations in Cavarroc et al. 2008 SPIE


    """


    miri = webbpsf_core.MIRI()

    miri.filter = 'FND'
    miri.image_mask = 'FQPM1140'


    im_nd_onaxis = miri.calcPSF(fov_arcsec=10)

    miri.options['source_offset_x'] = 2
    miri.options['source_offset_y'] = 2

    im_nd_offaxis = miri.calcPSF(fov_arcsec=10)

    P.subplot(211)
    poppy.display_PSF(im_nd_offaxis, colorbar=False)
    P.subplot(212)
    poppyt.display_PSF(im_nd_onaxis, colorbar=False)

    stop()


#---- TFI

def validate_tfi_coron():
    """ Verify that simulated performance of the TFI coronagraph is consistent with that shown in 
    Figure 10 of Doyon et al. 2010 SPIE 7731

    """

    return NotImplementedError('TFI has been deprecated.')
	#     tfi = webbpsf_core.TFI()
	# 
	#     # create test files
	#     tfi.etalon_wavelength = 4.6
	#     tfi.pupil_shift_x = 0.025 # assumep 2.5% pupil shear for consistency with Doyon et al. 
	#     image_masks = [None, 'CORON058', 'CORON075', 'CORON150', 'CORON150', 'CORON200']
	#     pupil_masks = [None, 'MASKC71N', 'MASKC71N', 'MASKC66N', 'MASKC21N', 'MASKC21N']
	#     throughputs = [1.0, 0.71, 0.71, 0.66, 0.21, 0.21]
	#     fns = ['test_tfi_4.6um_unocculted.fits', 'test_tfi_4.6um_058_c71n.fits', 'test_tfi_4.6um_075_c71n.fits', 
	#             'test_tfi_4.6um_150_c66n.fits', 'test_tfi_4.6um_150_c21n.fits', 'test_tfi_4.6um_200_c21n.fits']
	#     colors = ['black','red','purple','yellow','orange', 'cyan']
	# 
	#     for i in range(len(fns)):
	#         if not os.path.exists(fns[i]):
	#             tfi.image_mask = image_masks[i]
	#             tfi.pupil_mask = pupil_masks[i]
	#             tfi.calcPSF(fns[i], fov_arcsec=16, oversample=2, calc_oversample=4)
	# 
	#     # plot radial profiles
	#     seps = N.array([0.5, 1.0, 1.5, 2.0, 2.5, 5.0])  #arcsecs
	#     contrasts = N.array([6.9, 7.9, 9.2, 10.0, 10.9, 12.3]) # mags
	# 
	#     P.clf()
	#     peak_cts = None
	#     for i in range(len(fns)):
	#         radius, profile = webbpsf.radial_profile(fns[i])
	#         radius2, stds= webbpsf.radial_profile(fns[i], stddev=True)
	#         if peak_cts is None:
	#             peak_cts = profile.max()
	#         profile /= (peak_cts*throughputs[i])
	#         stds /= (peak_cts*throughputs[i])
	#         #P.semilogy(radius, profile, color=colors[i])
	#         P.semilogy(radius, 5*stds, color=colors[i], lw=5 if i >0 else 1)
	#         #stop()
	# 
	#     P.xlabel("Separation (arcsec)")
	#     P.ylabel("Contrast ratio (5$\sigma$)")
	#     ax = P.gca()
	#     ax.set_xlim(0,8)
	#     ax.set_ylim(1e-6, 1e0)



################################################################################
#
#       Comparisons vs JWPSF


def validate_vs_jwpsf_nircam():
    """ Compare results from WebbPSF with earlier simulations produced with JWPSF
    """

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

# Revision V and T pupils
def compare_revs_tv(nlambda=20, oversample=16, fov_arcsec=4):
    """ Compare PSF properties like FWHM and EE, for Rev T and Rev V pupil shapes"""
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
            PV = webbpsf_core.calc_or_load_psf('test_%s_%s_revV_%02d.fits' % (iname, filt, i+1), inst, nlambda=nlambda, oversample=oversample, fov_arcsec=fov_arcsec)
            inst.pupilopd = '/Users/mperrin/software/webbpsf/data/OPD_RevT/nircam_obs_w_rsrv%d.fits' % (i+1)
            PT = webbpsf_core.calc_or_load_psf('test_%s_%s_revT_%02d.fits' % (iname, filt, i+1), inst, nlambda=nlambda, oversample=oversample, fov_arcsec=fov_arcsec)

            fwhm_v[i] = poppy.measure_fwhm(PV)
            fwhm_t[i] = poppy.measure_fwhm(PT)
            ee_fn_v = poppy.measure_EE(PV)
            ee_fn_t = poppy.measure_EE(PT)
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
    """
    Compare PSFs with the Rev T and Rev V pupil shapes
    """
    P.clf()
    inst = webbpsf_core.NIRCam()
    inst.pupilopd=None

    fov_arcsec = 10
    nlambda=30

    inst.pupil = 'tricontagon.fits'
    psf_tri = webbpsf_core.calc_or_load_psf('test_NIRCam_perfect_tricontagon_o%d.fits' % oversample, inst, nlambda=nlambda, oversample=oversample, fov_arcsec=fov_arcsec)


    if not skipone:
        ax = P.subplot(1,3,1)
        poppy.display_PSF(psf_tri, normalize='peak', colorbar=False, title='Tricontagon')


        for i, rev in enumerate(['T','V']):
            inst.pupil = "pupil_Rev%s.fits" % rev
            psf = webbpsf_core.calc_or_load_psf('test_NIRCam_perfect_rev%s_o%d.fits'  % (rev, oversample),inst,  nlambda=nlambda, oversample=oversample, fov_arcsec=fov_arcsec)
            P.subplot(1,3,i+2)
            poppy.display_PSF(psf, normalize='peak', colorbar = (rev =='V'), title='OTE Rev '+rev)


        stop()
        P.clf()

    psf_V = fits.open('test_NIRCam_perfect_rev%s_o%d.fits'  % ('V', oversample))
    psf_T = fits.open('test_NIRCam_perfect_rev%s_o%d.fits'  % ('T', oversample))
    P.subplot(221)
    poppy.display_PSF_difference(psf_V, psf_tri, vmax=vmax, title="Rev V - tricontagon")
    P.subplot(222)
    poppy.display_PSF_difference(psf_V, psf_T,vmax=vmax, title="Rev V - Rev T")


    ax3 = P.subplot(223)
    ax3.set_ylabel('Azimuthally averaged profile')
    ax3.set_xlabel('Separation (arcsec)')

    ax4 = P.subplot(224)
    ax4.set_ylabel('Fractional Encircled Energy')
    ax4.set_xlabel('Separation (arcsec)')


    for p, label in zip([psf_tri, psf_T, psf_V], ['Tri', "Rev T", 'Rev V']):
        rad, rp = poppy.radial_profile(p)
        ee_fn = poppy.measure_EE(p)

        ax3.plot(rad,rp/rp.max())
        ax4.plot(rad, ee_fn(rad), label=label)

        print(poppy.measure_fwhm(p))
    ax4.legend(loc='lower right')
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    #ax3.axhline([psf_V[0].data.max()*0.5], ls=":")

    #ax3.set_xbound(0,4)
    ax3.set_xbound(0.01,4)
    ax4.set_xbound(0,4)
    P.draw()


    stop()

def validate_nircam_ee():

    # from BALL-JWST-SYST-05-003, Systems Eng Rpt on OTE Geometric Optical Model Config

    # Encircled energies at 1 micron, field position (0,0), 
    # page 18
    # from Code V model I think
    EE_fraction = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    EE_radius = [.004592, 0.009561, 0.015019, 0.021147, 0.028475, 0.037626, 0.052406, 0.081993, 0.204422] 


    # same thing, page 24, from OSLO model
    EE_radius = [0.006815, 0.013630, 0.020411, 0.026182, 0.031953, 0.037725, 0.058640, 0.090230, 0.210803]

    # TODO write this function sometime...


################################################################################

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(name)-12s: %(levelname)-8s %(message)s',)

