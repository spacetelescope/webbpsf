#!/usr/bin/env python
"""
obssim.py

    Observation Simulator wrapper for webbPSF. 

    This package lets you easily script simulations of things more complicated than just a single point source. 


"""
import numpy as N
import scipy.interpolate, scipy.ndimage
import pylab as P
import matplotlib
import atpy
import pyfits
import webbpsf


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

    def addPointSource(self, sptype_or_spectrum, separation=0.0, PA=0.0, normalization=None):
        if type(sptype_or_spectrum) is str:
            spectrum = specFromSpectralType(sptype_or_spectrum)
        else:
            spectrum = sptype_or_spectrum

        self.companions.append(   {'spectrum': sptype_or_spectrum, 'separation': separation, 'PA': PA, normalization)

    def calcImage(self, instrument, outfile=None, noise=False, rebin=True, **kwargs):
        """ Calculate an image of a scene through some instrument


        Parameters
        -----------
        instrument : webbpsf.jwinstrument instance
            A configured instance of an instrument class

        It may also be useful to pass arguments to the calcPSF() call, which is supported through the **kwargs 
        mechanism. Such arguments might include fov_arcsec, fov_pixels, oversample, etc.
        """

        star_psf = instrument.calcPSF(source = self.sourcespectrum, outfile=None, save_intermediates=False, rebin=rebin, **kwargs)

        for comp in companion:
            # set  companion spectrum and position
            comp_spectrum = None
            comp_psf =  instrument.calcPSF(source = comp_spectrum, outfile=None, save_intermediates=False, rebin=rebin, **kwargs)

            # figure out the flux ratio

            # add the scaled companion PSF to the stellar PSF:
            star_psf[0].data += comp_psf[0].data * comp_flux_ratio
            #update FITS header history
        if noise:
            pass
            #add noise in image - photon and read noise, mainly.
       
        # downsample? 
        if rebin and detector_oversample > 1:
            # throw away the existing rebinned extension

            # and generate a new one from the summed image
            _log.info(" Downsampling to detector pixel scale.")
            rebinned_result = result[0].copy()
            rebinned_result.data = rebin_array(rebinned_result.data, rc=(detector_oversample, detector_oversample))
            rebinned_result.header.update('OVERSAMP', 1, 'These data are rebinned to detector pixels')
            rebinned_result.header.update('CALCSAMP', detector_oversample, 'This much oversampling used in calculation')
            rebinned_result.header.update('EXTNAME', 'DET_SAMP')
            rebinned_result.header['PIXELSCL'] *= detector_oversample
            result.append(rebinned_result)



        if outfile is not None:
            result[0].header.update ("FILENAME", os.path.basename (outfile),
                           comment="Name of this file")
            result.writeto(outfile, clobber=clobber)
            _log.info("Saved result to "+outfile)
        return result

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




