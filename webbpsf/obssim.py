#!/usr/bin/env python
from __future__ import division, print_function, absolute_import, unicode_literals
"""
obssim.py

    Observation Simulator wrapper for webbPSF. 

    This package lets you easily script simulations of things more complicated than just a single point source. 


"""
import os
import numbers
import numpy as np
import scipy.interpolate, scipy.ndimage
import matplotlib.pyplot as plt
import matplotlib
import pysynphot
import logging
import poppy

import webbpsf_core


_log = logging.getLogger('webbpsf')
#
###########################################################################
#
#

class TargetScene(object):
    """ This class allows the user to specify some scene consisting of a central star
    plus one or more companions at specified separation, spectral type, etc. It automates the
    necessary calculations to perform a simulated JWST observation of that target. 

    pysynphot is required for this.



    """


    def __init__(self):
        self.sources = []

    def addPointSource(self, sptype_or_spectrum, name="unnamed source", separation=0.0, PA=0.0, normalization=None):
        """ Add a point source to the list for a given scene

        Parameters
        -----------
        sptype_or_spectrum : string or pysynphot.Spectrum
            spectrum of the source
        name : str
            descriptive string
        separation : float
            arcsec
        PA : float
            deg from N
        normalization : scalar or tuple TBD
            Simple version: this is a float to multiply the PSF by.
            Complex version: Probably tuple of arguments to spectrum.renorm(). 



        How normalization works:  
            First the PSF for that source is calculated, using calcPSF(norm='first')
            i.e. the input intensity through the telescope pupil is set to 1. 
            The resulting output PSF total counts will be proportional to the 
            throughput through the OTE+SI (including filters, coronagraphs etc)

            Then we apply the normalization:
                1) if it's just a number, we just multiply by it.
                2) if it's something else: Then we use a separate bandpass object and parameters 
                   passed in here to figure out the overall normalization, and apply that as a 
                   multiplicative factor to the resulting PSF itself?
        """
        if type(sptype_or_spectrum) is str:
            spectrum = poppy.specFromSpectralType(sptype_or_spectrum)
        else:
            spectrum = sptype_or_spectrum

        self.sources.append(   {'spectrum': sptype_or_spectrum, 'separation': separation, 'PA': PA, 
            'normalization': normalization, 'name': name})

    def calcImage(self, instrument, outfile=None, noise=False, rebin=True, clobber=True, 
            PA=0, offset_r=None, offset_PA=0.0, **kwargs):
        """ Calculate an image of a scene through some instrument


        Parameters
        -----------
        instrument : webbpsf.jwinstrument instance
            A configured instance of an instrument class
        outfile : str
            filename to save to
        rebin : bool
            passed to calcPSF
        PA : float
            postion angle for +Y direction in the output image
        offset_r, offset_PA : float
            Distance and angle to offset the target center from the FOV center.
            This is to simulate imperfect acquisition + alignment. 
        noise : bool
            add read noise? TBD
        clobber : bool
            overwrite existing files? default True


        It may also be useful to pass arguments to the calcPSF() call, which is supported through the **kwargs 
        mechanism. Such arguments might include fov_arcsec, fov_pixels, oversample, etc.
        """

        sum_image = None
        image_PA = PA

        for obj in self.sources:
            _log.info('Now propagating for '+obj['name'])
            # set  companion spectrum and position
            src_spectrum = obj['spectrum']

            if offset_r is None:
                instrument.options['source_offset_r'] = obj['separation']
                instrument.options['source_offset_theta'] = obj['PA'] - image_PA
            else:
                # combine the actual source position with the image offset position.
                obj_x = obj['separation'] * np.cos(obj['PA'] * np.pi/180)
                obj_y = obj['separation'] * np.sin(obj['PA'] * np.pi/180)
                offset_x = offset_r * np.cos(offset_PA * np.pi/180)
                offset_y = offset_r * np.sin(offset_PA * np.pi/180)

                src_x = obj_x + offset_x
                src_y = obj_y + offset_y
                src_r = np.sqrt(src_x**2+src_y**2)
                src_pa = np.arctan2(src_y, src_x) * 180/np.pi
                instrument.options['source_offset_r'] = src_r
                instrument.options['source_offset_theta'] = src_pa - image_PA
                #stop()

            _log.info('  post-offset & rot pos: %.3f  at %.1f deg' % (instrument.options['source_offset_r'], instrument.options['source_offset_theta']))


            src_psf =  instrument.calcPSF(source = src_spectrum, outfile=None, save_intermediates=False, rebin=rebin, 
                **kwargs)

            # figure out the flux ratio
            if obj['normalization'] is not None:
                # use the explicitly-provided normalization:
                if isinstance(obj['normalization'], numbers.Number):
                    src_psf[0].data *= obj['normalization']
                else:
                    raise NotImplemented("Not Yet")
            else:
                # use the flux level already implicitly set by the source spectrum.
                # i.e. figure out what the flux of the source is, inside the selected bandpass
                bp = instrument._getSynphotBandpass()
                effstim_Jy = pysynphot.Observation(src_spectrum, bp).effstim('Jy')
                src_psf[0].data *= effstim_Jy
 
            # add the scaled companion PSF to the stellar PSF:
            if sum_image is None:
                sum_image = src_psf
                sum_image[0].header.add_history("obssim : Creating an image simulation with multiple PSFs")
                sum_image[0].header['IMAGE_PA'] = ( image_PA,'PA of scene in simulated image')
                sum_image[0].header['OFFSET_R'] = (0 if offset_r is None else offset_r ,'[arcsec] Offset of target center from FOV center')
                sum_image[0].header['OFFSETPA'] = (0 if offset_PA is None else offset_PA ,'[deg] Position angle of target offset from FOV center')

                if offset_r is None:
                    sum_image[0].header.add_history("Image is centered on target (perfect acquisition)")
                else:
                    sum_image[0].header.add_history("Image is offset %.2f arcsec at PA=%.1f from target" % (offset_r, offset_PA))

            else:
                sum_image[0].data += src_psf[0].data
            #update FITS header history
            sum_image[0].header.add_history("Added source %s at r=%.3f, theta=%.2f" % (obj['name'], obj['separation'], obj['PA']))
            sum_image[0].header.add_history("                with effstim = %.3g Jy" % effstim_Jy)
            sum_image[0].header.add_history("                counts in image: %.3g" % src_psf[0].data.sum())
            sum_image[0].header.add_history("                pos in image: %.3g'' at %.1f deg" % (instrument.options['source_offset_r'],  instrument.options['source_offset_theta'])  )


        if noise:
            raise NotImplemented("Not Yet")

        sum_image[0].header['NSOURCES'] = ( len(self.sources), "Number of point sources in sim")
            #add noise in image - photon and read noise, mainly.
       
        # downsample? 
        if rebin and sum_image[0].header['DET_SAMP'] > 1:
            # throw away the existing rebinned extension
            sum_image.pop() 
            # and generate a new one from the summed image
            _log.info(" Downsampling summed image to detector pixel scale.")
            rebinned_sum_image = sum_image[0].copy()
            detector_oversample = sum_image[0].header['DET_SAMP']
            rebinned_sum_image.data = poppy.rebin_array(rebinned_sum_image.data, rc=(detector_oversample, detector_oversample))
            rebinned_sum_image.header['OVERSAMP'] = ( 1, 'These data are rebinned to detector pixels')
            rebinned_sum_image.header['CALCSAMP'] = ( detector_oversample, 'This much oversampling used in calculation')
            rebinned_sum_image.header['EXTNAME'] = ( 'DET_SAMP')
            rebinned_sum_image.header['PIXELSCL'] *= detector_oversample
            sum_image.append(rebinned_sum_image)



        if outfile is not None:
            sum_image[0].header["FILENAME"] = ( os.path.basename (outfile), "Name of this file")
            sum_image.writeto(outfile, clobber=clobber)
            _log.info("Saved image to "+outfile)
        return sum_image

    def display(self):
        plt.clf()
        for obj in self.sources:
            X = obj['separation'] * -np.sin(obj['PA'] * np.pi/180)
            Y = obj['separation'] * np.cos(obj['PA'] * np.pi/180)

            plt.plot([X],[Y],'*')
            plt.text(X,Y, obj['name'])




def test_obssim(nlambda=3, clobber=False):
    s = TargetScene()

    s.addPointSource('G0V', name='G0V star', separation = 0.1, normalization=1.)
    s.addPointSource('K0V', name='K0V star', separation = 1.0, PA=45,  normalization=0.4)
    s.addPointSource('M0V', name='M0V star', separation = 1.5, PA=245,  normalization=0.3)

    inst = webbpsf_core.NIRCam()

    for filt in ['F115W', 'F210M', 'F360M']:
        inst.filter = filt
        outname = "test_scene_%s.fits"% filt
        if not os.path.exists(outname) or clobber:
            s.calcImage(inst, outfile=outname, fov_arcsec=5, nlambda=nlambda)






if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')



