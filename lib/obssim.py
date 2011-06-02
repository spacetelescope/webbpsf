#!/usr/bin/env python
"""
obssim.py

    Observation Simulator wrapper for webbPSF. 

    This package lets you easily script simulations of things more complicated than just a single point source. 


"""
import os
import numbers
import numpy as N
import scipy.interpolate, scipy.ndimage
import pylab as P
import matplotlib
import pysynphot
import atpy
import pyfits
import logging
import webbpsf

_log = logging.getLogger('obssim')
_log.setLevel(logging.DEBUG)
_log.setLevel(logging.INFO)
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
            spectrum = specFromSpectralType(sptype_or_spectrum)
        else:
            spectrum = sptype_or_spectrum

        self.sources.append(   {'spectrum': sptype_or_spectrum, 'separation': separation, 'PA': PA, 
            'normalization': normalization, 'name': name})

    def calcImage(self, instrument, outfile=None, noise=False, rebin=True, PA=0, clobber=True, **kwargs):
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
            # set  companion spectrum and position
            src_spectrum = obj['spectrum']
            instrument.options['source_offset_r'] = obj['separation']
            instrument.options['source_offset_theta'] = obj['PA'] - image_PA

            src_psf =  instrument.calcPSF(source = src_spectrum, outfile=None, save_intermediates=False, rebin=rebin, 
                **kwargs)

            # figure out the flux ratio
            if obj['normalization'] is not None:
                if isinstance(obj['normalization'], numbers.Number):
                    src_psf[0].data *= obj['normalization']
                else:
                    raise NotImplemented("Not Yet")

            # add the scaled companion PSF to the stellar PSF:
            if sum_image is None:
                sum_image = src_psf
                sum_image[0].header.add_history("obssim : Creating an image simulation with multiple PSFs")
            else:
                sum_image[0].data += src_psf[0].data
            #update FITS header history
            sum_image[0].header.add_history("Added source %s at r=%.3f, theta=%.2f" % (obj['name'], obj['separation'], obj['PA']))
        if noise:
            raise NotImplemented("Not Yet")

        sum_image[0].header.update('NSOURCES', len(self.sources), "Number of point sources in sim")
            #add noise in image - photon and read noise, mainly.
       
        # downsample? 
        if rebin and sum_image[0].header['DET_SAMP'] > 1:
            # throw away the existing rebinned extension
            sum_image.pop() 
            # and generate a new one from the summed image
            _log.info(" Downsampling summed image to detector pixel scale.")
            rebinned_sum_image = sum_image[0].copy()
            detector_oversample = sum_image[0].header['DET_SAMP']
            rebinned_sum_image.data = webbpsf.rebin_array(rebinned_sum_image.data, rc=(detector_oversample, detector_oversample))
            rebinned_sum_image.header.update('OVERSAMP', 1, 'These data are rebinned to detector pixels')
            rebinned_sum_image.header.update('CALCSAMP', detector_oversample, 'This much oversampling used in calculation')
            rebinned_sum_image.header.update('EXTNAME', 'DET_SAMP')
            rebinned_sum_image.header['PIXELSCL'] *= detector_oversample
            sum_image.append(rebinned_sum_image)



        if outfile is not None:
            sum_image[0].header.update ("FILENAME", os.path.basename (outfile),
                           comment="Name of this file")
            sum_image.writeto(outfile, clobber=clobber)
            _log.info("Saved image to "+outfile)
        return sum_image


def test_obssim(nlambda=3, clobber=False):
    s = TargetScene()

    s.addPointSource('G0V', name='G0V star', separation = 0.1, normalization=1.)
    s.addPointSource('K0V', name='K0V star', separation = 1.0, PA=45,  normalization=0.4)
    s.addPointSource('M0V', name='M0V star', separation = 1.5, PA=245,  normalization=0.3)

    inst = webbpsf.NIRCam()

    for filt in ['F115W', 'F210M', 'F360M']:
        inst.filter = filt
        outname = "test_scene_%s.fits"% filt
        if not os.path.exists(outname) or clobber:
            s.calcImage(inst, outfile=outname, fov_arcsec=5, nlambda=nlambda)




if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')




def specFromSpectralType(sptype, return_list=False):
    """Get Pysynphot Spectrum object from a spectral type string.

    """
    lookuptable = {
        "O3V":   (50000, 0.0, 5.0),
        "O5V":   (45000, 0.0, 5.0),
        "O6V":   (40000, 0.0, 4.5),
        "O8V":   (35000, 0.0, 4.0),
        "O5I":   (40000, 0.0, 4.5),
        "O6I":   (40000, 0.0, 4.5),
        "O8I":   (34000, 0.0, 4.0),
        "B0V":   (30000, 0.0, 4.0),
        "B3V":   (19000, 0.0, 4.0),
        "B5V":   (15000, 0.0, 4.0),
        "B8V":   (12000, 0.0, 4.0),
        "B0III": (29000, 0.0, 3.5),
        "B5III": (15000, 0.0, 3.5),
        "B0I":   (26000, 0.0, 3.0),
        "B5I":   (14000, 0.0, 2.5),
        "A0V":   (9500, 0.0, 4.0),
        "A5V":   (8250, 0.0, 4.5),
        "A0I":   (9750, 0.0, 2.0),
        "A5I":   (8500, 0.0, 2.0),
        "F0V":   (7250, 0.0, 4.5),
        "F5V":   (6500, 0.0, 4.5),
        "F0I":   (7750, 0.0, 2.0),
        "F5I":   (7000, 0.0, 1.5),
        "G0V":   (6000, 0.0, 4.5),
        "G5V":   (5750, 0.0, 4.5),
        "G0III": (5750, 0.0, 3.0),
        "G5III": (5250, 0.0, 2.5),
        "G0I":   (5500, 0.0, 1.5),
        "G5I":   (4750, 0.0, 1.0),
        "K0V":   (5250, 0.0, 4.5),
        "K5V":   (4250, 0.0, 4.5),
        "K0III": (4750, 0.0, 2.0),
        "K5III": (4000, 0.0, 1.5),
        "K0I":   (4500, 0.0, 1.0),
        "K5I":   (3750, 0.0, 0.5),
        "M0V":   (3750, 0.0, 4.5),
        "M2V":   (3500, 0.0, 4.5),
        "M5V":   (3500, 0.0, 5.0),
        "M0III": (3750, 0.0, 1.5),
        "M0I":   (3750, 0.0, 0.0),
        "M2I":   (3500, 0.0, 0.0)}


    if return_list:
        sptype_list = lookuptable.keys()
        def sort_sptype(typestr):
            letter = typestr[0]
            lettervals = {'O':0, 'B': 10, 'A': 20,'F': 30, 'G':40, 'K': 50, 'M':60}
            value = lettervals[letter]*1.0
            value += int(typestr[1])
            if "III" in typestr: value += .3
            elif "I" in typestr: value += .1
            elif "V" in typestr: value += .5
            return value
        sptype_list.sort(key=sort_sptype)
        return sptype_list

    try:
        keys = lookuptable[sptype]
    except:
        raise LookupError("Lookup table does not include spectral type %s" % sptype)

    return pysynphot.Icat('ck04models',keys[0], keys[1], keys[2])




