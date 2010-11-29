#!/usr/bin/env python
"""
Physical Optics Propagation in PYthon (POPPY)


This package implements an object-oriented system for modeling physical optics
propagation with diffraction, particularly for telescopic and coronagraphic
imaging. Right now only image and pupil planes are supported; intermediate
planes are a future goal.

Classes:
 * Wavefront
 * OpticalElement
   * AnalyticOpticalElement
     * BandLimitedCoron
     * IdealMonoFQPM
     * IdealFieldStop
     * IdealCircularOcculter
   * Detector
 * OpticalSystem


    Code by Marshall Perrin <mperrin@stsci.edu>

"""

#import os, sys
import numpy as N
import pylab as p
import pyfits
import scipy.special
import copy
import matplotlib
from matplotlib.colors import LogNorm  # for log scaling of images, with automatic colorbar support

import SFT

from IPython.Debugger import Tracer; stop = Tracer()


# constants for types of plane
PUPIL = 1
IMAGE = 2
DETECTOR = 3 # specialized type of image plane.


#conversions
ARCSECtoDEGREES = 1. / 3600.
RADIANStoDEGREES = 180. / N.pi
RADIANStoARCSEC = 180.*60*60 / N.pi
MICRONStoMETERS = 1.e-6
MICRONStoNANOMETERS = 1000.


def padToOversample(array, oversample):
    """ Add zeros around the edge of an array.

    Parameters
    ----------
    array :  ndarray
        A 2D array representing some image
    oversample : int
        Padding factor for expanding the array

    Returns
    -------
    padded_array : ndarray
        A larger array containing mostly zeros but with the input array in the center.
    """
    npix = array.shape[0]
    padded = N.zeros(shape=(npix*oversample, npix*oversample))
    n0 = float(npix)*(oversample - 1)/2
    n1 = n0+npix
    padded[n0:n1, n0:n1] = array
    return padded

def removePadding(array,oversample):
    " Remove zeros around the edge of an array "
    npix = array.shape[0] / oversample
    n0 = float(npix)*(oversample - 1)/2
    n1 = n0+npix
    return array[n0:n1,n0:n1].copy()

#------
class Wavefront():
    """ A class representing a monochromatic wavefront that can be transformed between
    pupil and image planes (but not to intermediate planes, yet).

    In a pupil plane, a wavefront object 'wf' has
        wf.diam         diameter in meters
        wf.pixelscale   scale in meters/pixel
    In an image plane, it has
        wf.fov          field of view in arcseconds
        wf.pixelscale    scale in arcsec/pixel


    Use the wf.propagateTo() method to transform a wavefront between conjugate planes.


    """
    def __init__(self,wavelength=2e-6, npix=1024, dtype=N.complex128, diam=8.0, oversample=2, pixelscale=None):
        """
        Create a wavefront object. By default, is in a pupil plane. Set pixelscale= to make an image plane.

        Parameters
        ----------
        wavelength : float
            Wavelength of light in meters
        npix : int
            Size parameter for wavefront array to create, per side.
        diam : float, optional
            For PUPIL wavefronts, sets physical size corresponding to npix. Units are meters.
            At most one of diam or pixelscale should be set when creating a wavefront.
        pixelscale : float, optional
            For IMAGE PLANE wavefronts, use this pixel scale.
        oversample : int, optional
            how much to oversample by in FFTs. Default is 2.


        dtype : numpy.dtype, optional
            default is double complex.

        """

        if wavelength > 1e-4:
            raise ValueError("The specified wavelength is implausibly large. Remember to specify the desired wavelength in *meters*.")

        self.oversample = oversample

        self.wavelength = float(wavelength)                 # wavelen in meters, obviously
        self.diam= float(diam)                              # pupil plane size in meters
        self.fov = None                                     # image plane size in arcsec

        if pixelscale is None:
            self.pixelscale = self.diam / npix                  # scale in meters/pix or arcsec/pix, as appropriate
            self.planetype = PUPIL                              # are we at image or pupil?
        else:
            self.pixelscale = pixelscale
            self.planetype = IMAGE
        self.wavefront = N.ones((npix,npix), dtype=dtype)   # the actual complex wavefront array
        self.ispadded = False                               # is the wavefront padded for oversampling?
        self.history=[]
        self.history.append("Created wavefront: wavelen=%g m, diam=%f m" %(self.wavelength, self.diam))
        self.history.append(" using array size %s" % (self.wavefront.shape,) )
        self.location='Entrance'

    def __str__(self):
        # TODO add switches for image/pupil planes
        return """Wavefront:
        wavelength = %f microns
        shape = (%d,%d)
        sampling = %f meters/pixel
        """ % (self.wavelength/1e-6, self.wavefront.shape[0], self.wavefront.shape[1], self.pixelscale )

    def copy(self):
        return copy.deepcopy(self)

    def normalize(self):
        "Set this wavefront's total intensity to 1 "
        #print "Wavefront normalized"
        self.wavefront /= N.sqrt(self.totalIntensity)

    def __imul__(self, optic):
        "Multiply a Wavefront by an OpticalElement or scalar"
        if (isinstance(optic,float)) or isinstance(optic,int):
            self.wavefront *= optic # it's just a scalar
            self.history.append("Multiplied WF by scalar value "+str(optic))
            return self


        if (not isinstance(optic, OpticalElement)) :
            raise ValueError('Wavefronts can only be *= multiplied by OpticalElements or scalar values')

        if isinstance(optic,Detector):
            # detectors don't modify a wavefront.
            return self

        phasor = optic.getPhasor(self)

        if len(phasor) > 1:
            assert self.wavefront.shape == phasor.shape

        self.wavefront *= phasor
        if optic.verbose: print "  Multiplied WF by phasor for "+str(optic)
        self.history.append("Multiplied WF by phasor for "+str(optic))
        self.location='after '+optic.name
        return self

    def __mul__(self, optic):
        """ Multiply a wavefront by an OpticalElement or scalar """
        new = self.copy()
        new *= optic
        return new

    def __iadd__(self,wave):
        "Add another wavefront to this one"
        if not isinstance(wave,Wavefront):
            raise ValueError('Wavefronts can only be summed with other Wavefronts')

        if not self.wavefront.shape[0] == wave.wavefront.shape[0]:
            raise ValueError('Wavefronts can only be added if they have the same size and shape')

        self.wavefront += wave.wavefront
        self.history.append("Summed with another wavefront!")
        return self

    def asFITS(self, what='intensity'):
        """ Return a wavefront as a pyFITS HDUList object

        Parameters
        -----------
        what : string
            what kind of data to write. Must be one of 'parts', 'intensity', 'complex'.
            The default is to write a file containing intensity.

        """

        if what.lower() =='parts':
            outarr = N.zeros((2,self.shape[0], self.shape[1]))
            outarr[0,:,:] = self.amplitude
            outarr[1,:,:] = self.phase
            outFITS = pyfits.HDUList(pyfits.PrimaryHDU(outarr))
            outFITS[0].header.update('PLANE1', 'Wavefront Amplitude')
            outFITS[0].header.update('PLANE2', 'Wavefront Phase')
        elif what.lower() =='intensity':
            outFITS = pyfits.HDUList(pyfits.PrimaryHDU(self.intensity))
            outFITS[0].header.update('PLANE1', 'Wavefront Intensity')
        elif what.lower()  == 'complex':
            outFITS = pyfits.HDUList(pyfits.PrimaryHDU(self.wavefront))
            outFITS[0].header.update('PLANE1', 'Wavefront Complex Phasor ')

        outFITS[0].header.update('WAVELEN', self.wavelength, 'Wavelength in meters')
        outFITS[0].header.update('DIFFLMT', self.wavelength/self.diam*206265., 'Diffraction limit lambda/D in arcsec')
        outFITS[0].header.update('OVERSAMP', self.oversample, 'Oversampling factor for computation')
        if self.planetype ==IMAGE:
            outFITS[0].header.update('PIXELSCL', self.pixelscale, 'Pixel scale in arcsec/pixel')
            outFITS[0].header.update('FOV', self.fov, 'Field of view in arcsec (full array)')
        else:
            outFITS[0].header.update('PIXELSCL', self.pixelscale, 'Pixel scale in meters/pixel')
            outFITS[0].header.update('DIAM', self.diam, 'Pupil diameter in meters (not incl padding)')

        for h in self.history: outFITS[0].header.add_history(h)

        return outFITS

    def writeto(self,filename, clobber=True, **kwargs):
        """Write a wavefront to a FITS file.

        Parameters
        -----------
        filename : string
            filename to use
        what: string
            what to write. Must be one of 'parts', 'intensity', 'complex'
        clobber: bool, optional
            overwhat existing? default is True

        """
        self.asFITS(**kwargs).writeto(filename, clobber=clobber)
        print("  Wavefront saved to %s" % filename)

    def display(self,what='intensity', nrows=1,row=1,showpadding=False,imagezoom=5.0):
        """Display wavefront on screen

        Parameters
        ----------
        what : string
           What to display. Must be one of {intensity, phase, best}.
           'Best' implies 'display the phase if there is OPD, or else
           display the intensity for a perfect pupil.


        """

        intens = N.ma.masked_array(self.intensity, mask=(self.intensity==0))
        phase = N.ma.masked_array(self.phase, mask=(intens==0))
        amp = self.amplitude

        if not showpadding and self.ispadded:
            if self.planetype == PUPIL:
                intens = removePadding(intens,self.oversample)
                phase = removePadding(phase,self.oversample)
                amp = removePadding(amp,self.oversample)


        if self.planetype == PUPIL:
            extent = [0,self.pixelscale*intens.shape[0], 0,self.pixelscale*intens.shape[1]]
            unit = "m"
            norm=matplotlib.colors.Normalize(vmin=0)
        else:
            halffov = self.pixelscale*intens.shape[0]/2
            extent = [-halffov, halffov, -halffov, halffov]
            unit="arcsec"
            norm=LogNorm(vmin=1e-8,vmax=1e-1)

        cmap = matplotlib.cm.jet
        cmap.set_bad('k', 0.8)


        if what =='best':
            if self.planetype ==IMAGE: what = 'intensity' # always show intensity for image planes
            elif phase.sum() == 0: what = 'intensity' # for perfect pupils
            else: what='phase' # for aberrated pupils

        if what == 'intensity':
            nc = int(N.ceil(N.sqrt(nrows)))
            nr = int(N.ceil(float(nrows)/nc))
            p.subplot(nr,nc,int(row))
            p.imshow(intens,extent=extent, norm=norm, cmap=cmap)
            p.title("Intensity "+self.location)
            p.xlabel(unit)
            p.colorbar(orientation='vertical', shrink=0.8)

            if self.planetype ==IMAGE:
                p.axhline(0,ls="k:")
                p.axvline(0,ls="k:")
                ax = p.gca()
                imsize = min( (imagezoom, halffov))
                ax.set_xbound(-imsize, imsize)
                ax.set_ybound(-imsize, imsize)
        elif what =='phase':
            nc = int(N.ceil(N.sqrt(nrows)))
            nr = int(N.ceil(float(nrows)/nc))
            p.subplot(nr,nc,int(row))
            p.imshow(phase,extent=extent, norm=norm, cmap=cmap)
            p.title("Phase "+self.location)
            p.xlabel(unit)
            p.colorbar(orientation='vertical', shrink=0.8)



        else:
            p.subplot(nrows,2,(row*2)-1)
            p.imshow(amp,extent=extent,cmap=cmap)
            p.title("Wavefront amplitude")
            p.ylabel(unit)
            p.colorbar(orientation='vertical',shrink=0.8)

            p.subplot(nrows,2,row*2)
            p.imshow(self.phase,extent=extent, cmap=cmap)
            p.ylabel(unit)
            p.title("Wavefront phase")

        p.draw()

    # add convenient properties for intensity, phase, amplitude, total_flux
    @property
    def amplitude(self):
        return N.abs(self.wavefront)

    @property
    def intensity(self):
        return N.abs(self.wavefront)**2

    @property
    def phase(self):
        return N.angle(self.wavefront)

    @property
    def shape(self):
        return self.wavefront.shape

    @property
    def totalIntensity(self):
        return self.intensity.sum()

    # methods for wavefront propagation:
    def propagateTo(self, optic):
        """Propagates a wavefront object to the next optic in the list.
        Modifies this wavefront object itself.
        """
        if self.planetype == optic.planetype:
            print "  Wavefront and optic %s already at same plane type, no propagation needed." % optic.name
            return
        else:
            msg = "  Propagating wavefront to %s. " % str(optic)
            print(msg)
            self.history.append(msg)


        if optic.planetype == DETECTOR:
            self._propagateMFT(optic)
        else:
            self._propagateFFT(optic)

        self.location='before '+optic.name

    def _propagateFFT(self, optic):
        """ Propagate from pupil to image or vice versa using a padded FFT """


        if self.oversample > 1 and not self.ispadded: #add padding for oversampling
            assert self.oversample == optic.oversample
            self.wavefront = padToOversample(self.wavefront, self.oversample)
            self.ispadded = True
            if optic.verbose: print "    Padded WF array for oversampling by %dx" % self.oversample
            self.history.append("    Padded WF array for oversampling by %dx" % self.oversample)

        if self.planetype == PUPIL and optic.planetype == IMAGE:
            FFT_direction = 'FORWARD'
            # do FFT
            #print "Pre-FFT total intensity: "+str(self.totalIntensity)
                # due to annoying normalization convention in numpy fft, we have to divide by the number of pixels
                # when doing this fft step:
            self.wavefront = N.fft.fft2(self.wavefront) / self.wavefront.shape[0]
            self.wavefront = N.fft.fftshift(self.wavefront)
            #print "Post-FFT total intensity: "+str(self.totalIntensity)
            # update keywords
            self.planetype=IMAGE
            self.pixelscale = self.wavelength/ self.diam / self.oversample * RADIANStoARCSEC
            self.fov = self.wavefront.shape[0] * self.pixelscale

            self.history.append('   FFT %s,  to IMAGE  scale=%f' %(self.wavefront.shape, self.pixelscale))
            #self.fov = det.fov_arcsec
            #self.pixelscale = det.fov_arcsec / det_calc_size_pixels


        elif self.planetype == IMAGE and optic.planetype ==PUPIL:
            FFT_direction= 'BACKWARD'
            # do FFT
            #print "Pre-FFT total intensity: "+str(self.totalIntensity)
                # due to annoying normalization convention in numpy fft, we have to divide by the number of pixels
                # when doing this fft step:
            self.wavefront = N.fft.ifftshift(self.wavefront) * self.wavefront.shape[0]
            self.wavefront = N.fft.ifft2(self.wavefront)
            #print "Post-FFT total intensity: "+str(self.totalIntensity)
            # update keywords
            self.planetype=PUPIL
            self.pixelscale = self.diam *self.oversample / self.wavefront.shape[0]
            self.history.append('   FFT %s,  to PUPIL scale=%f' %(self.wavefront.shape, self.pixelscale))
            #self.diam = det.fov_arcsec
            #self.pixelscale = det.fov_arcsec / det_calc_size_pixels

    def _propagateMFT(self, det):
        """ Compute from pupil to an image using the Soummer et al. 2007 MFT algorithm"""

        assert self.planetype == PUPIL
        assert det.planetype == DETECTOR

        if self.ispadded:
            #pupil plane is padded - trim that out since it's not needed
            self.wavefront = removePadding(self.wavefront, self.oversample)
            self.ispadded = False


        # the arguments for the SFT are
        # - wavefront (assumed to fill the input array)
        # - focal plane size in lambda/D units
        # - number of pixels on a side in focal plane array.

        lamD = self.wavelength / self.diam * RADIANStoARCSEC
        #print "lam/D = %f arcsec" % lamD

        det_fov_lamD = det.fov_arcsec / lamD
        det_calc_size_pixels = det.fov_npix * det.oversample

        mfter = SFT.SlowFourierTransform()
        msg= '    Propagating w/ MFT: %.4f"/pix     fov=%.3f lam/D    npix=%d' % \
            (det.pixelscale/det.oversample, det_fov_lamD, det_calc_size_pixels)
        print(msg)
        self.history.append(msg)
        self.wavefront = mfter.perform(self.wavefront, det_fov_lamD, det_calc_size_pixels)

        self.planetype=IMAGE
        self.fov = det.fov_arcsec
        self.pixelscale = det.fov_arcsec / det_calc_size_pixels

    def tilt(self, Xangle=0.0, Yangle=0.0):
        """ Tilt a wavefront in X and Y. """
        if self.planetype==IMAGE:
            raise NotImplementedError("Are you sure you want to tilt a wavefront in an IMAGE plane?")

        #Compute the tilt of the wavefront required to shift it by some amount in the image plane.

        tiltphasor = 1.

        self.wavefront *= tiltphasor
        self.history.append("Tilted wavefront")


#------
class OpticalElement():
    """ Defines an arbitrary optic, based on amplitude transmission and/or OPD files.

    This optic could be a pupil or field stop, an aberrated mirror, a phase mask, etc.
    The OpticalElement class follows the behavoior of the Wavefront class, using units
    of meters/pixel in pupil space and arcsec/pixel in image space.


    NOTE: All mask files must be **squares**.
    """
    def __init__(self, name="unnamed optic", transmission=None, opd= None, verbose=True, planetype=None, oversample=1,opdunits="meters", shift=None):
        """
        Parameters
        ----------
        name : string
            descriptive name for optic
        tranmission, opd : string
            FITS filenames for the transmission (from 0-1) and opd (in meters)
        opdunits: string
            units for the OPD file. Default is 'meters'.
        verbose : bool
            whether to print stuff while computing
        planetype : int
            either IMAGE or PUPIL
        oversample : int
            how much to oversample beyond Nyquist.

        shift : tuple of floats
            2-tuple containing X and Y fractional shifts for the pupil.


        """


        self.name = name
        self.verbose=verbose

        self.amplitude_header = None
        self.opd_header = None

        self.planetype = planetype      # pupil or image
        self.oversample = oversample    # oversampling factor, none by default
        self.ispadded = False           # are we padded w/ zeros for oversampling the FFT?

        # Read amplitude and/or OPD from disk
        if opd is None and transmission is None:   # no input files, so just make a scalar
            print "No input files specified. "
            print "Creating a null optical element. Are you sure that's what you want to do?"
            self.amplitude = N.asarray([1.])
            self.opd = N.asarray([0.])
            self.pixelscale = 0
            self.name = "-empty-"
        else:
            if transmission is not None:        # load transmission file.
                self.amplitude_file = transmission
                self.amplitude, self.amplitude_header = pyfits.getdata(self.amplitude_file, header=True)
                if len (self.amplitude.shape) != 2 or self.amplitude.shape[0] != self.amplitude.shape[1]:
                    raise ValueError, "OPD image must be 2-D and square"
                if self.verbose:
                    print(self.name+": Loaded amplitude from "+self.amplitude_file)
            else:                               # else if only OPD set, create an array of 1s with same size.
                opd_shape = pyfits.getdata(opd).shape
                self.amplitude = N.ones(opd_shape)

            if opd is not None:         # Load OPD file.
                # if OPD is specified as a tuple, treat the first element as the filename and 2nd as the slice of a cube.
                if isinstance(opd, basestring):
                    self.opd_file = opd
                    self.opd_slice = 0
                    datacube = False
                else:
                    self.opd_file = opd[0]
                    self.opd_slice = opd[1]
                    datacube = True

                self.opd, self.opd_header = pyfits.getdata(self.opd_file, header=True)

                if datacube:
                    self.opd = self.opd[self.opd_slice, :,:]

                if len (self.opd.shape) != 2 or self.opd.shape[0] != self.opd.shape[1]:
                    raise ValueError, "OPD image must be 2-D and square"
                if not self.verbose:
                    print(self.name+": Loaded opd from "+self.opd_file)


                # convert OPD into meters
                if opdunits.lower() == 'meters' or opdunits.lower() == 'meter':
                    pass # no need to rescale
                elif opdunits.lower() == 'micron' or opdunits.lower() == 'um':
                    self.opd *= 1e-6
            else:                   #else if only amplitude set, create an array of 0s with same size.
                self.opd = N.zeros(self.amplitude.shape)

            assert self.amplitude.shape == self.opd.shape
            assert self.amplitude.shape[0] == self.amplitude.shape[1]

            # if a shift is specified and we're NOT a null (scalar) optic, then do the shift:
            if shift is not None and len(self.amplitude.shape) ==2:
                if abs(shift[0]) > 0.5 or abs(shift[1])> 0.5:
                    raise ValueError("""You have asked for an implausibly large shift. Remember, shifts should be specified as
                      decimal values between 0.0 and 1.0, a fraction of the total optic diameter. """)
                rolly = int(N.round(self.amplitude.shape[0] * shift[1])) #remember Y,X order for shape, but X,Y order for shift
                rollx = int(N.round(self.amplitude.shape[1] * shift[0]))
                print "Requested optic shift of (%6.3f, %6.3f) %%" % (shift)
                print "Actual shift applied   = (%6.3f, %6.3f) %%" % (rollx*1.0/self.amplitude.shape[1], rolly *1.0/ self.amplitude.shape[0])

                self.amplitude = N.roll(self.amplitude, rolly, axis=0)
                self.amplitude = N.roll(self.amplitude, rollx, axis=1)
                self.opd       = N.roll(self.opd,       rolly, axis=0)
                self.opd       = N.roll(self.opd,       rollx, axis=1)




            if self.planetype == PUPIL:
                try:
                    self.pupil_scale = self.amplitude_header['PUPLSCAL']
                    self.pupil_diam = self.amplitude_header['PUPLDIAM']
                    self.pixelscale = self.pupil_scale # synonyms
                except:
                    raise ValueError('That pupil appears to be missing the required PUPLDIAM keyword.')
            elif self.planetype == IMAGE:
                 try:
                    self.pixelscale = self.amplitude_header['PIXSCALE']
                 except:
                    raise ValueError('That image appears to be missing the required PIXSCALE keyword.')

    def getPhasor(self,wave):
        """ Compute a complex phasor from an OPD, given a wavelength.

        Parameters
        ----------
        wave : float or obj
            either a scalar wavelength or a Wavefront object

        """

        if isinstance(wave, Wavefront):
            wavelength=wave.wavelength
        else:
            wavelength=wave
        scale = 2. * N.pi / wavelength

        # set the self.phasor attribute:
        # first check whether we need to interpolate to do this.
        float_tolerance = 0.0001  #how big of a relative scale mismatch before resampling?
        if hasattr(wave,'pixelscale') and abs(wave.pixelscale -self.pixelscale)/self.pixelscale >= float_tolerance:
            print wave.pixelscale, self.pixelscale
            raise ValueError("Non-matching pixel scale for wavefront and optic! Need to add interpolation / rescaling ")
            if self.has_attr('_resampled_scale') and abs(self._resampled_scale-wave.pixelscale)/self._resampled_scale >= float_tolerance:
                # we already did this same resampling, so just re-use it!
                self.phasor = self._resampled_amplitude * N.exp (1.j * self._resampled_opd * scale)
            else:
                raise NotImplementedError("Need to implement resampling.")

        else:
            # compute the phasor directly, without any need to rescale.
            self.phasor = self.amplitude * N.exp (1.j * self.opd * scale)



        # check whether we need to pad before returning or not.
        if self.planetype == PUPIL and wave.ispadded:
            return padToOversample(self.phasor, wave.oversample)
        else:
            return self.phasor

    def display(self, nrows=1, row=1, phase=False, wavelength=None):
        "Display plots showing an optic's transmission and OPD"
        if self.planetype == PUPIL:
            pixelscale = self.pupil_scale
            units = "meters"
        else:
            pixelscale = self.pixelscale
            units = "arcsec"

        extent = [0,pixelscale*self.amplitude.shape[0], 0,pixelscale*self.amplitude.shape[1]]

        if nrows == 1: orient = "horizontal"
        else: orient = "vertical"

        cmap = matplotlib.cm.jet
        cmap.set_bad('k', 0.8)
        norm_amp=matplotlib.colors.Normalize(vmin=0, vmax=1)
        norm_opd=matplotlib.colors.Normalize(vmin=-0.5e-6, vmax=0.5e-6)

        ampl = N.ma.masked_equal(self.amplitude, 0)
        opd= N.ma.masked_array(self.opd, mask=(self.amplitude ==0))

        p.subplot(nrows, 2, row*2-1)
        p.imshow(ampl, extent=extent, cmap=cmap, norm=norm_amp)
        p.title("Transmissivity for "+self.name)
        p.colorbar(orientation=orient, ticks=[0,0.25, 0.5, 0.75, 1.0])

        p.subplot(nrows, 2, row*2)
        p.imshow(opd, extent=extent, cmap=cmap, norm=norm_opd)
        p.title("OPD for "+self.name)
        cb = p.colorbar(orientation=orient, ticks=N.array([-0.5, -0.25, 0, 0.25, 0.5])*1e-6)
        cb.set_label('meters')

    def __str__(self):
        typestrs = ['', 'Pupil plane', 'Image plane', 'Detector']
        if self.planetype is PUPIL:
            return "Pupil plane: %s (%dx%d pixels, diam=%f m)" % (self.name, self.shape[0], self.shape[0], self.pupil_diam)
        elif self.planetype is IMAGE:
            return "Image plane: %s (%dx%d pixels, scale=%f arcsec/pixel'')" % (self.name, self.shape[0], self.shape[0], self.pixelscale)
        else:
            return "Optic: "+self.name

    @property
    def shape(self):
        if hasattr(self, 'amplitude'):
            return self.amplitude.shape
        else: return [1,1]

    def apply(self, wavefront):
        phasor = self.getPhasor(wavefront.wavelength)
        if not N.isscalar(phasor):
            assert wavefront.shape == phasor.shape

class AnalyticOpticalElement(OpticalElement):
    """ Defines an abstract analytic optical element.
        This class is useless on its own and must have a getPhasor routine
        provided by a derived subclass.
        It exists mostly to provide some default keyword handling in __init__
    """
    def __init__(self, name="unnamed", verbose=True, oversample=1, transmission=None, opd=None, planetype=None):
        self.name = name
        self.verbose=verbose
        self.planetype=planetype
        self.oversample = oversample # this is irrelevant for analytic functions, but is included for compatibility with
                                     # the base OpticalElement class.
        if transmission is not None:
            raise ValueError("AnalyticOpticalElements should not have a transmission= argument specified.")
        if opd is not None:
            raise ValueError("AnalyticOpticalElements should not have a opd= argument specified.")

    def __str__(self):
        typestrs = ['', 'Pupil plane', 'Image plane', 'Detector']
        if self.planetype is PUPIL:
            return "Pupil plane: %s (Analytic)" % (self.name)
        elif self.planetype is IMAGE:
            return "Image plane: %s (Analytic)" % (self.name)
        else:
            return "Optic: "+self.name

    def getPhasor(self,wave):
        raise NotImplementedError("getPhasor must be supplied by a derived subclass")

    def display(self, nrows=1, row=1, phase=False, wavelength=2e-6):
        "Display an Analytic optic by first computing it onto a grid..."
        if self.planetype is PUPIL:
            unit="meters"
            halffov = 4.0
        else:
            unit="arcsec"
            halffov = 2.0

        w = Wavefront(wavelength=wavelength, npix=512,  pixelscale = 2.*halffov/512)

        phasor = self.getPhasor(w)
        self.amplitude = N.abs(phasor)
        phase = N.angle(phasor) * 2*N.pi
        self.opd = phase

        stop()  # rewrite this to set properties appropriately then call parent class display

        extent = [-halffov, halffov, -halffov, halffov]

        p.subplot(nrows,2,(row*2)-1)
        p.imshow(trans,extent=extent)
        p.title("Transmission for "+self.name)
        p.ylabel(unit)
        p.colorbar(orientation='vertical')

        p.subplot(nrows,2,row*2)
        p.imshow(phase,extent=extent)
        p.ylabel(unit)
        p.title("Phase for "+self.name)

class BandLimitedCoron(AnalyticOpticalElement):
    """ Defines an ideal band limited coronagraph occulting mask

    """
    def __init__(self, name="unnamed BLC", kind='circular', sigma=1, **kwargs):
        AnalyticOpticalElement.__init__(self,**kwargs)
        self.name = name
        self.verbose=verbose
        self.planetype=IMAGE

        self.kind = kind.lower()        # either circular or linear
        if self.kind not in ['circular', 'linear']:
            raise ValueError("Invalid kind of BLC: "+self.kind)
        self.sigma = sigma              # size parameter. See section 2.1 of Krist et al. SPIE 2009

    def getPhasor(self,wave):
        """ Compute the amplitude transmission appropriate for a BLC for some given pixel spacing
        corresponding to the supplied Wavefront
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("BLC getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == IMAGE)

        #phasor = N.zeros(wave.shape)
        y, x = N.indices(wf.shape)
        y-= wf.shape[0]/2
        if self.kind == 'circular':
            x-= wf.shape[1]/2
            r = N.sqrt(x**2+y**2) * wave.pixelscale
            sigmar = sigma*r
            self.transmission = (1-  (2*scipy.special.jn(1,sigmar)/sigmar)**2)**2

        elif self.kind == 'linear':
            r = abs(y)*pixelscale # FIXME TODO change scale linearly across mask
            sigmar = sigma*r
            self.transmission = (1-  (2*           N.sin(1,sigmar)/sigmar)**2)**2

        return self.transmission

class IdealMonoFQPM(AnalyticOpticalElement):
    """ Defines an ideal monochromatic 4-quadrant phase mask coronagraph

    """
    def __init__(self, name="unnamed FQPM ", wavelength=10.65e-6, **kwargs):
        AnalyticOpticalElement.__init__(self,**kwargs)
        self.name = name

        self.central_wavelength =wavelength

    def getPhasor(self,wave):
        """ Compute the amplitude transmission appropriate for a 4QPM for some given pixel spacing
        corresponding to the supplied Wavefront
        """

        if not isinstance(wave, Wavefront):
            raise ValueError("4QPM getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == IMAGE)

        # TODO this computation could be sped up a lot w/ optimzations
        phase = N.empty(wave.shape)
        n0 = wave.shape[0]/2
        phase[:n0,:n0]=0.5
        phase[n0:,n0:]=0.5
        phase[n0:,:n0]=0
        phase[:n0,n0:]=0

        retardance = phase*self.central_wavelength/wave.wavelength
        FQPM_phasor = N.exp(1.j * 2* N.pi * retardance)
        return FQPM_phasor

class IdealFieldStop(AnalyticOpticalElement):
    """ Defines an ideal square field stop
    """

    def __init__(self, name="unnamed field stop",  size=20., **kwargs):
        AnalyticOpticalElement.__init__(self,**kwargs)
        self.name = name
        self.size = size            # size of square stop in arcseconds.
        self.pixelscale=0

    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the field stop.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("IdealFieldStop getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == IMAGE)

        #phasor = N.zeros(wave.shape)
        y, x = N.indices(wave.shape)
        y -= wave.shape[0]/2
        x -= wave.shape[1]/2
        self.transmission = N.ones(wave.shape)

        halfsize_pixels = self.size  / wave.pixelscale / 2
        w_outside = N.where( (abs(y) > halfsize_pixels)  | (abs(x) > halfsize_pixels))
        self.transmission[w_outside] = 0

        return self.transmission

class IdealCircularOcculter(AnalyticOpticalElement):
    """ Defines an ideal circular occulter
    """

    def __init__(self, name="unnamed occulter",  radius=1.0, **kwargs):
        AnalyticOpticalElement.__init__(self,**kwargs)
        self.name = name
        self.radius = radius    # radius of circular occulter in arcseconds.
        self.pixelscale=0

    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the occulter.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("IdealCircleOcculter getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == IMAGE)

        y, x = N.indices(wave.shape)
        y -= wave.shape[0]/2
        x -= wave.shape[1]/2
        r = N.sqrt(x**2+y**2) * wave.pixelscale
        self.transmission = N.zeros(wave.shape)

        w_inside = N.where( r <= self.radius)
        self.transmission[w_inside] = 0

        return self.transmission

class Detector(OpticalElement):
    """ A Detector is a specialized type of OpticalElement that forces a wavefront
    onto a specific fixed pixelization.
    """
    def __init__(self, pixelscale, fov_npix=None, fov_arcsec=None, oversample=1, name="Unnamed"):
        self.name=name
        self.planetype = DETECTOR
        self.pixelscale = float(pixelscale)
        self.oversample = oversample

        if fov_npix is None and fov_arcsec is None:
            raise ValueError("Either fov_npix or fov_arcsec must be specified!")
        elif fov_npix is not None:
            self.fov_npix = int(fov_npix)
            self.fov_arcsec = self.fov_npix * self.pixelscale
        else:
            # set field of view to closest value possible to requested,
            # consistent with having an integer number of pixels
            self.fov_npix = round(fov_arcsec / self.pixelscale)
            self.fov_arcsec = self.fov_npix * self.pixelscale

        self.shape = (self.fov_npix, self.fov_npix)

        self.amplitude = 1
        self.opd = 0

    def __str__(self):
        return "Detector plane: %s (%dx%d, %f arcsec/pixel)" % (self.name, self.fov_npix, self.fov_npix, self.pixelscale)

#------
class OpticalSystem():
    """ A class representing a series of optical elements,
    either Pupil, Image, or Detector planes, through which light
    can be propagated.

    The difference between
    Image and Detector planes is that Detectors have fixed pixels
    in terms of arcsec/pixel regardless of wavelength (computed via
    MFT) while Image planes have variable pixels scaled in terms of
    lambda/D. Pupil planes are some fixed size in meters, of course.

    Parameters
    ----------
    name :
        string, descriptive name of optical system
    oversample :
        int. Either how many times *above* Nyquist we should be
        (for pupil or image planes), or how many times a fixed
        detector pixel will be sampled. E.g. oversample=2 means
        image plane sampling lambda/4*D (twice Nyquist) and
        detector plane sampling 2x2 computed pixels per real detector
        pixel.  Default is 2.
    verbose :
        bool, whether to print stuff while computing




    """
    def __init__(self, name="unnamed system", verbose=True, oversample=2):
        self.name = name
        self.verbose=verbose
        self.planes = []                    # List of OpticalElements
        self.oversample = oversample

        self.source_tilt = N.zeros((2))     # off-axis tilt of the source, in ANGULAR units.

        self.intermediate_wfs = None        #
        if self.verbose:
            print "Created OpticalSystem: "+self.name

    # Methods for adding or manipulating optical planes:

    def addPupil(self, **kwargs):
        " Add a pupil plane optic from file(s) giving transmission or OPD"
        self.planes.append(OpticalElement(planetype=PUPIL, oversample=self.oversample, **kwargs))
        if self.verbose: print "Added pupil plane: "+self.planes[-1].name

    def addImage(self, optic=None, function=None, **kwargs):
        """ Add an image plane optic, either
          1) from file(s) giving transmission or OPD
                [set arguments transmission=filename and/or opd=filename]
          2) from an analytic function
                [set function='circle, fieldstop, bandlimitedcoron, or FQPM'
                and set additional kwargs to define shape etc.
          3) from an already-created OpticalElement object
                [set optic=that object]

        """
        if optic is None:
            if function == 'CircularOcculter':
                fn = IdealCircularOcculter
            elif function == 'fieldstop':
                fn = IdealFieldStop
            elif function == 'BandLimitedCoron':
                fn = BandLimitedCoron
            elif function == 'FQPM':
                fn = IdealMonoFQPM
            elif function is not None:
                raise ValueError("Analytic mask type %s is unknown." % function)
            else: # create image from files specified in kwargs
                fn = OpticalElement

            optic = fn(planetype=IMAGE, oversample=self.oversample, **kwargs)
        else:
            optic.planetype = IMAGE

        self.planes.append(optic)
        if self.verbose: print "Added image plane: "+self.planes[-1].name

    def addDetector(self, *args, **kwargs):
        self.planes.append(Detector(oversample=self.oversample, *args, **kwargs))
        if self.verbose: print "Added detector: "+self.planes[-1].name


        return "Optical system '%s' containing %d optics" % (self.name, len(self.planes))

    def list(self):
        print str(self)+"\n\t"+ "\n\t".join([str(p) for p in self.planes])

    def __getitem__(self, num):
        return self.planes[num]

    # methods for dealing with wavefronts:
    def inputWavefront(self, wavelength=2e-6):
        """Create a Wavefront object suitable for sending through a given optical system, based on
        the size of the first optical plane, assumed to be a pupil.

        Uses self.source_tilt to assign an off-axis tilt, if requested.

        """
        inwave = Wavefront(wavelength=wavelength,
                npix = self.planes[0].shape[0],
                diam = self.planes[0].pupil_diam,
                oversample=self.oversample)

        if abs(self.source_tilt).sum() >= 0:
            raise NotImplementedError("Need to implement shift of target")

        return inwave

    def propagate_mono(self, wavelength=2e-6, normalize='first', save_intermediates=False, display_intermediates=False, intermediate_fn='wave_step_%03d.fits', poly_weight=None):
        """ Propagate a wavefront through some number of optics.
        Returns a pyfits.HDUList object.

        Parameters
        ----------
        * wavelength    Wavelength in meters
        * normalize     how to normalize the wavefront?
                        'first' = set total flux = 1 after the first optic, presumably a pupil
                        'last' = set total flux = 1 after the entire optical system.

        * poly_weight   is this being called as part of a polychromatic calculation?
                        if not, set this to None. if so, set this to the weight for
                        that wavelength.
        * display_intermediates
        * save_intermediates

        poly_weight controls whether intermediate optical planes are actually saved to disk by this routine
        (for the monochromatic case) or are handled in calcPSF (for the polychromatic case).


        """
        #if not isinstance(wavefront, Wavefront):
            #raise TypeError("propagate must be called with a valid Wavefront object.")


        if self.verbose: print "\n*** Propagating wavelength = %g meters" % wavelength
        wavefront = self.inputWavefront(wavelength)

        if save_intermediates and poly_weight is None:
            self.intermediate_wfs=[]
            print "reset intermediates"
        #if display_intermediates and poly_weight is None: p.clf()
        # need to CLF due to obnoxious color bar re-creation otherwise
        if display_intermediates: p.clf()

        # do the propagation:
        count = 0
        for optic in self.planes:
            wavefront.propagateTo(optic)
            wavefront *= optic
            count += 1


            if normalize.lower()=='first' and count == 1:
                wavefront.normalize()


            print "  Flux === "+str(wavefront.totalIntensity)

            if save_intermediates:
                if len(self.intermediate_wfs) < count:
                    self.intermediate_wfs.append(wavefront.copy())
                else:
                    self.intermediate_wfs[count-1] += wavefront.copy()*poly_weight
                if poly_weight is None: self.intermediate_wfs[count-1].writeto(intermediate_fn % count, what='parts')
            if display_intermediates:
               if save_intermediates: self.intermediate_wfs[count-1].display(what='best',nrows=len(self.planes),row=count)
               else: wavefront.display(what='best',nrows=len(self.planes),row=count)

        # prepare output arrays
        if normalize.lower()=='last':
                wavefront.normalize()

        return wavefront.asFITS()

    def calcPSF(self, source, save_intermediates=False, **kwargs):
        """Calculate a multi-wavelength PSF over some weighted
        sum of wavelengths.

        Parameters
        ---------
        source :
            a dict containing 'wavelengths' and 'weights' list.
            TBD - replace w/ pysynphot observation object
        s # always show intensity for image planesave_intermediates :
            Bool, whether to output intermediate optical planes


        Returns
        -------
        outfits :
            a pyfits.HDUList
        """


        if save_intermediates:
            self.intermediate_wfs = []
            print 'reset intermediates in calcPSF'
            #raise ValueError("Saving intermediates for multi-wavelen not yet implemented!!")

        # loop over wavelengths
        if self.verbose: print "** Calculating PSF with %d wavelengths" % (len(source['wavelengths']))
        outFITS = None

        normwts =  N.asarray(source['weights'])
        normwts /= normwts.sum()

        for wavelen, weight in zip(source['wavelengths'], normwts):
            mono_psf = self.propagate_mono(wavelen, poly_weight=weight, save_intermediates=save_intermediates, **kwargs)
            # add mono_psf into the output array

            if outFITS is None:
                outFITS = mono_psf
                outFITS[0].data = mono_psf[0].data*weight
            else:
                outFITS[0].data += mono_psf[0].data *weight

        if save_intermediates:
            for i in range(len(self.intermediate_wfs)):
                self.intermediate_wfs[i].writeto('wave_step_%03d.fits' % i )

        waves = N.asarray(source['wavelengths'])
        wts = N.asarray(source['weights'])
        mnwave = (waves*wts).sum() / wts.sum()
        outFITS[0].header.update('WAVELEN', mnwave, 'Weighted mean wavelength in meters')
        outFITS[0].header.update('NWAVES',waves.size, 'Number of wavelengths used in calculation')
        for i in range(waves.size):
            outFITS[0].header.update('WAVE'+str(i), waves[i], "Wavelength "+str(i))
            outFITS[0].header.update('WGHT'+str(i), wts[i], "Wavelength weight "+str(i))

        if self.verbose: print "** PSF Calculation completed."
        return outFITS

    def display(self, **kwargs):

        planes_to_display = [p for p in self.planes if not isinstance(p, Detector)]
        nplanes = len(planes_to_display)
        for i in range(nplanes):
            self.planes[i].display(nrows=nplanes, row=1, **kwargs)


def test_MFT():
    osys = OpticalSystem("Perfect JW", oversample=2)
    osys.addPupil(transmission="/Users/mperrin/software/newJWPSF/data/pupil.fits", name='JW Pupil')
    osys.addDetector(0.032, fov_npix=128)


    out = osys.propagate_mono(wavelength=2e-6)
    out.writeto('test_mono.fits',clobber=True)


    source = dict()
    source['wavelengths'] = [2.0e-6, 2.1e-6, 2.2e-6]
    source['weights'] = [0.3, 0.5, 0.2]
    out2 = osys.calcPSF(source)
    out2.writeto('test_MW.fits', clobber=True)



def test_poppy():
    p.clf()

    osys = OpticalSystem("Perfect JW", oversample=4)
    osys.addPupil(transmission="/Users/mperrin/software/newJWPSF/data/pupil.fits", name='JW Pupil')
    osys.addImage(function='fieldstop', name='20 arcsec stop', size=20)
    osys.addImage(function='FQPM',wavelength=10.65e-6)
    osys.addPupil(transmission="/Users/mperrin/software/newJWPSF/data/MIRI/coronagraph/MIRI_FQPMLyotStop.fits", name='MIRI FQPM Lyot')
    osys.addDetector(0.032, name="Detector", fov_npix=128)

    out = osys.propagate_mono(wavelength=10.65e-6, display_intermediates=True, save_intermediates=True)
    out.writeto('test_fft.fits',clobber=True)





if __name__ == "__main__":

    #test_MFT()
    test_poppy()

