#!/usr/bin/env python
"""
=============================================
Physical Optics Propagation in PYthon (POPPY)
=============================================


This package implements an object-oriented system for modeling physical optics
propagation with diffraction, particularly for telescopic and coronagraphic
imaging. Right now only Fraunhoffer diffraction is supported, providing image and pupil planes; 
intermediate planes are a future goal.

Classes:
--------
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


Module-level configuration constants
------------------------------------

_USE_FFTW3 : bool
    Should the FFTW3 library be used? Set automatically to True if fftw3 is importable, else False.
_USE_MULTIPROC : bool
    Should we use python multiprocessing to spawn multiple simultaneous processes to speed calculations?
    Note that FFTW3 automatically makes use of multiple cores so this should not be used at the same time as
    _USE_FFTW3. Default is False.
_MULTIPROC_NPROCESS : int
    If the above is set, how many processes should be used? Default is 4. Note that you probably don't want to 
    set this too large due to memory limitations; i.e. don't just blindly set this to 16 on a 16-core machine unless
    you have many, many GB of RAM...

_TIMETESTS : bool
    Print out simple benchmarking of elapsed time to screen. Default False
_FLUXCHECK : bool
    Print out total flux after each step of a propagation. Useful for debugging, mostly.

_IMAGECROP : float
    Default zoom in region for image displays. Default is 5.0


"""

from __future__ import division
#import os, sys
import multiprocessing
import copy
import numpy as N
import pylab as p
import pylab as P
import pyfits
import scipy.special
import scipy.ndimage.interpolation
import matplotlib
import time
from matplotlib.colors import LogNorm  # for log scaling of images, with automatic colorbar support
import SFT

__version__ = '0.2.2'

try:
    from IPython.Debugger import Tracer; stop = Tracer()
except:
    def stop():
        pass

import logging
_log = logging.getLogger('poppy')
_log.addHandler(logging.NullHandler())


try:
    import fftw3
    _USE_FFTW3 = True
    _FFTW3_INIT = {}
except:
    _USE_FFTW3 = False

_USE_MULTIPROC = False # auto set this cleverly?
_MULTIPROC_NPROCESS = min(4,multiprocessing.cpu_count())  # Caution: Do not make this too large on high-CPU-count machines
                                                          # because this is a memory-intensive calculation and you willg
                                                          # just end up thrashing IO and swapping out a ton, so everything
                                                          # becomes super slow.

_TIMETESTS=False #set to true for benchmarking
_FLUXCHECK= False
_IMAGECROP = 5.0 # default image display is 5 arcsec

# constants for types of plane
PUPIL = 1
IMAGE = 2
DETECTOR = 3 # specialized type of image plane.
ROTATION = 4 # not a real optic, just a coordinate transform
typestrs = ['', 'Pupil plane', 'Image plane', 'Detector']


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
    padded = N.zeros(shape=(npix*oversample, npix*oversample), dtype=array.dtype)
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

def _wrap_propagate_for_multiprocessing(args):
    """ This is an internal helper routine for parallelizing computations across multiple processors.
   g
    Python's multiprocessing module allows easy execution of tasks acrossg
    many CPUs or even distinct machines. It relies on Python's pickle mechanism to
    serialize and pass objects between processes. One annoying side effect of this is
    that object instance methods cannot easily be pickled, and thus cannot be easilyg
    invoked in other processes.g

    Here, we work around that by pickling the entire object and argument list, packed
    as a tuple, transmitting that to the new process, and then unpickling that,g
    unpacking the results, and *then* at last making our instance method call.g
    """
    self, wavelength, weight, kwargs, usefftw3flag = args
    _USE_FFTW3 = usefftw3flag

    return args[0].propagate_mono(wavelength, poly_weight=weight, save_intermediates=False, **kwargs)

def imshow_with_mouseover(image, ax=None,  *args, **kwargs):
    """ Wrapper for pyplot.imshow that sets up a custom mouseover display formatter
    so that mouse motions over the image are labeled in the status bar with
    pixel numerical value as well as X and Y coords.

    Why this behavior isn't the matplotlib default, I have no idea...
    """
    if ax is None: ax = p.gca()
    myax = ax.imshow(image, *args, **kwargs)
    aximage = ax.images[0].properties()['array']
    # need to account for half pixel offset of array coordinates for mouseover relative to pixel center,
    # so that the whole pixel from e.g. ( 1.5, 1.5) to (2.5, 2.5) is labeled with the coordinates of pixel (2,2)


    # We use the extent and implementation to map back from the data coord to pixel coord
    # There is probably an easier way to do this...
    imext = ax.images[0].get_extent()  # returns [-X, X, -Y, Y]
    imsize = ax.images[0].get_size()   # returns [sY, sX]g
    # map data coords back to pixel coords:
    #pixx = (x - imext[0])/(imext[1]-imext[0])*imsize[1]
    #pixy = (y - imext[2])/(imext[3]-imext[2])*imsize[0]
    # and be sure to clip appropriatedly to avoid array bounds errors
    report_pixel = lambda x, y : "(%6.3f, %6.3f)     %g" % \
        (x,y,   aximage[N.floor( (y - imext[2])/(imext[3]-imext[2])*imsize[0]  ).clip(0,imsize[0]-1),\
                        N.floor( (x - imext[0])/(imext[1]-imext[0])*imsize[1]  ).clip(0,imsize[1]-1)])

        #(x,y, aximage[N.floor(y+0.5),N.floor(x+0.5)])   # this works for regular pixels w/out an explicit extent= call
    ax.format_coord = report_pixel

    return ax


#------
class Wavefront(object):
    """ A class representing a monochromatic wavefront that can be transformed between
    pupil and image planes (but not to intermediate planes, yet).

    In a pupil plane, a wavefront object `wf` has
        * `wf.diam`,         a diameter in meters
        * `wf.pixelscale`,   a scale in meters/pixel
    In an image plane, it has
        * `wf.fov`,          a field of view in arcseconds
        * `wf.pixelscale`,   a  scale in arcsec/pixel


    Use the `wf.propagateTo()` method to transform a wavefront between conjugate planes. This will update those properties as appropriate.

    By default, `Wavefronts` are created in a pupil plane. Set `pixelscale=#` to make an image plane instead.

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
        Note that final propagations to Detectors use a different algorithmg
        and, optionally, a separate oversampling factor.
    dtype : numpy.dtype, optional
        default is double complex.

    """

    #planetype = IMAGE
    "Is this a PUPIL or IMAGE plane? Uses constants from poppy to define."
    #intensity = 0
    "numpy.ndarray.  The intensity of the wavefront as a function of position."
    #amplitude = 0
    "numpy.ndarray.  Amplitude of the electric field"
    #phase = 0
    "numpy.ndarray.  The phase of the wavefront, in units of radians."

    def __init__(self,wavelength=2e-6, npix=1024, dtype=N.complex128, diam=8.0, oversample=2, pixelscale=None):

        if wavelength > 1e-4:
            raise ValueError("The specified wavelength %f is implausibly large. Remember to specify the desired wavelength in *meters*." % wavelength)

        self._last_transform_type=None # later used to track MFT vs FFT pixel coord centering in coordinates()
        self.oversample = oversample

        self.wavelength = float(wavelength)                 # wavelen in meters, obviously
        """Wavelength in meters """
        self.diam= float(diam)                              # pupil plane size in meters
        """Diameter in meters. Applies to a pupil plane only."""
        self.fov = None                                     # image plane size in arcsec
        """Field of view in arcsec. Applies to an image plane only."""
        self.pixelscale = None
        "Pixel scale, in arcsec/pixel or meters/pixel depending on plane type"

        if pixelscale is None:
            self.pixelscale = self.diam / npix                  # scale in meters/pix or arcsec/pix, as appropriate
            self.planetype = PUPIL                              # are we at image or pupil?
        else:
            self.pixelscale = pixelscale
            self.planetype = IMAGE
        self._image_centered='array_center'                     # one of 'array_center', 'pixel', 'corner'
                                                                # This records where the coordinate origin is
                                                                # in image planes, and depends on how the imageg
                                                                # plane was produced (e.g. FFT implies pixel)
        "Are FFT'ed image planes centered on a pixel or on a corner between pixels? "
        self.wavefront = N.ones((npix,npix), dtype=dtype)   # the actual complex wavefront array
        self.ispadded = False                               # is the wavefront padded for oversampling?
        self.history=[]
        "List of strings giving a descriptive history of actions performed on the wavefront. Saved to FITS headers."
        self.history.append("Created wavefront: wavelen=%g m, diam=%f m" %(self.wavelength, self.diam))
        self.history.append(" using array size %s" % (self.wavefront.shape,) )
        self.location='Entrance'
        "A descriptive string for where a wavefront is instantaneously located (e.g. 'before occulter'). Used mostly for titling displayed plots."

    def __str__(self):
        # TODO add switches for image/pupil planes
        return """Wavefront:
        wavelength = %f microns
        shape = (%d,%d)
        sampling = %f meters/pixel
        """ % (self.wavelength/1e-6, self.wavefront.shape[0], self.wavefront.shape[1], self.pixelscale )

    def copy(self):
        "Return a copy of the wavefront as a different object."
        return copy.deepcopy(self)

    def normalize(self):
        "Set this wavefront's total intensity to 1 "
        #print "Wavefront normalized"
        self.wavefront /= N.sqrt(self.totalIntensity)

    def __imul__(self, optic):
        "Multiply a Wavefront by an OpticalElement or scalar"
        if isinstance(optic,Rotation):
            return self # a rotation doesn't actually affect the wavefront via multiplication,
                        # but instead via forcing a call to rotate() elsewhere...
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

        if not N.isscalar(phasor) and phasor.size>1:  # actually isscalar() does not handle the case of a 1-element array properly
            #print self.wavefront.shape, phasor.shape
            #stop()
            assert self.wavefront.shape == phasor.shape

        self.wavefront *= phasor
        msg =  "  Multiplied WF by phasor for "+str(optic)
        _log.debug(msg)
        self.history.append(msg)
        self.location='after '+optic.name
        return self

    def __mul__(self, optic):
        """ Multiply a wavefront by an OpticalElement or scalar """
        new = self.copy()
        new *= optic
        return new
    __rmul__ = __mul__  # either way works.


    def __iadd__(self,wave):
        "Add another wavefront to this one"
        if not isinstance(wave,Wavefront):
            raise ValueError('Wavefronts can only be summed with other Wavefronts')

        if not self.wavefront.shape[0] == wave.wavefront.shape[0]:
            raise ValueError('Wavefronts can only be added if they have the same size and shape')

        self.wavefront += wave.wavefront
        self.history.append("Summed with another wavefront!")
        return self

    def __add__(self,wave):
        new = self.copy()
        new += wave
        return new

    def asFITS(self, what='intensity'):
        """ Return a wavefront as a pyFITS HDUList object

        Parameters
        -----------
        what : string
            what kind of data to write. Must be one of 'parts', 'intensity', 'complex'.
            The default is to write a file containing intensity.

        """
        if what.lower() =='all':
            outarr = N.zeros((3,self.shape[0], self.shape[1]))
            outarr[0,:,:] = self.intensity
            outarr[1,:,:] = self.amplitude
            outarr[2,:,:] = self.phase
            outFITS = pyfits.HDUList(pyfits.PrimaryHDU(outarr))
            outFITS[0].header.update('PLANE1', 'Wavefront Intensity')
            outFITS[0].header.update('PLANE2', 'Wavefront Amplitude')
            outFITS[0].header.update('PLANE3', 'Wavefront Phase')
        elif what.lower() =='parts':
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
        outFITS[0].header.update('OVERSAMP', self.oversample, 'Oversampling factor for FFTs in computation')
        outFITS[0].header.update('DET_SAMP', self.oversample, 'Oversampling factor for MFT to detector plane')
        if self.planetype ==IMAGE:
            outFITS[0].header.update('PIXELSCL', self.pixelscale, 'Scale in arcsec/pix (after oversampling)')
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
        what : string
            what to write. Must be one of 'parts', 'intensity', 'complex'
        clobber : bool, optional
            overwhat existing? default is True

        Returns
        -------
        outfile: file on disk
            The output is written to disk.

        """
        self.asFITS(**kwargs).writeto(filename, clobber=clobber)
        _log.info("  Wavefront saved to %s" % filename)

    def display(self,what='intensity', nrows=1,row=1,showpadding=False,imagecrop=None, colorbar=False, crosshairs=True, ax=None):
        """Display wavefront on screen

        Parameters
        ----------
        what : string
           What to display. Must be one of {intensity, phase, best}.
           'Best' implies 'display the phase if there is OPD, or else
           display the intensity for a perfect pupil.

        nrows : int
            Number of rows to display in current figure (used for showing steps in a calculation)
        row : int
            Which row to display this one in?
        imagecrop : float, optional
            For image planes, set the maximum # of arcseconds to display. Default is 5, so
            only the innermost 5x5 arcsecond region will be shown.
        showpadding : bool, optional
            Show the entire padded arrays, or just the good parts? Default is False
        colorbar : bool
            Display colorbar
        ax : matplotlib Axes
            axes to display into

        Returns
        -------
        figure : matplotlib figure
            The current figure is modified.


        """
        if imagecrop is None: imagecrop = _IMAGECROP

        #intens = N.ma.masked_array(self.intensity, mask=(self.intensity==0))
        intens = self.intensity.copy()
        phase = self.phase.copy()
        phase[N.where(intens ==0)] = N.nan
        #phase = N.ma.masked_array(self.phase, mask=(intens==0))
        amp = self.amplitude

        if not showpadding and self.ispadded:
            if self.planetype == PUPIL:
                intens = removePadding(intens,self.oversample)
                phase = removePadding(phase,self.oversample)
                amp = removePadding(amp,self.oversample)


        # extent specifications need to include the *full* data region, including the half pixel on either
        # side outside of the pixel center coordinates.  And remember to swap Y and X
        extent = N.array([-0.5,intens.shape[1]-1+0.5, -0.5,intens.shape[0]-1+0.5]) * self.pixelscale
        if self.planetype == PUPIL:
            unit = "m"
        else:
            # make coordinates relative to center.
            # image plane coordinates depend on whether theg
            if self._image_centered == 'array_center' or self._image_centered=='corner':
                cenx = (intens.shape[0]-1)/2.
            elif self._image_centered == 'pixel':
                cenx = (intens.shape[0])/2.

            extent -= cenx*self.pixelscale
            halffov = intens.shape[0]/2.*self.pixelscale #for use later
            #(N.array([-0.5,intens.shape[1]-1+0.5, -0.5,intens.shape[0]-1+0.5])-cenx) * self.pixelscale
            #halffov = self.pixelscale*intens.shape[0]/2
            #extent = [-halffov, halffov, -halffov, halffov]
            unit="arcsec"

        # implement semi-intellegent selection of what to display
        if what =='best':
            if self.planetype ==IMAGE:
                what = 'intensity' # always show intensity for image planes
            elif phase[N.where(N.isfinite(phase))].sum() == 0:
                what = 'intensity' # for perfect pupils
            elif int(row) > 1: what='intensity'  # show intensity for coronagraphic downstream propagation.
            else: what='phase' # for aberrated pupils

        # compute plot parameters
        nc = int(N.ceil(N.sqrt(nrows)))
        nr = int(N.ceil(float(nrows)/nc))
        if nrows - nc*(nc-1) == 1: # avoid just one alone on a row by itself...
            nr -= 1
            nc += 1

        # now display the chosen selection..
        if what == 'intensity':
            if self.planetype == PUPIL:
                norm=matplotlib.colors.Normalize(vmin=0)
                cmap = matplotlib.cm.gray
                cmap.set_bad('0.0')
            else:
                norm=LogNorm(vmin=1e-8,vmax=1e-1)
                cmap = matplotlib.cm.jet
                cmap.set_bad(cmap(0))

            if ax is None:
                ax = p.subplot(nr,nc,int(row))
            imshow_with_mouseover(intens, ax=ax, extent=extent, norm=norm, cmap=cmap)
            title = "Intensity "+self.location
            title = title.replace('after', 'after\n')
            title = title.replace('before', 'before\n')
            p.title(title)
            p.xlabel(unit)
            if colorbar: p.colorbar(ax.images[0], orientation='vertical', shrink=0.8)

            if self.planetype ==IMAGE:
                if crosshairs:
                    p.axhline(0,ls=":", color='k')
                    p.axvline(0,ls=":", color='k')
                imsize = min( (imagecrop, halffov))
                ax.set_xbound(-imsize, imsize)
                ax.set_ybound(-imsize, imsize)
        elif what =='phase':
            # Display phase in waves.
            cmap = matplotlib.cm.jet
            cmap.set_bad('0.3')
            norm=matplotlib.colors.Normalize(vmin=-0.25,vmax=0.25)
            if ax is None:
                ax = p.subplot(nr,nc,int(row))
            imshow_with_mouseover(phase/(N.pi*2), ax=ax, extent=extent, norm=norm, cmap=cmap)
            p.title("Phase "+self.location)
            p.xlabel(unit)
            if colorbar: p.colorbar(ax2.images[0], orientation='vertical', shrink=0.8)


        else:
            p.subplot(nrows,2,(row*2)-1)
            p.imshow(amp,extent=extent,cmap=cmap)
            p.title("Wavefront amplitude")
            p.ylabel(unit)
            if colorbar: p.colorbar(orientation='vertical',shrink=0.8)

            p.subplot(nrows,2,row*2)
            p.imshow(self.phase,extent=extent, cmap=cmap)
            p.ylabel(unit)
            p.title("Wavefront phase")

        p.draw()

    # add convenient properties for intensity, phase, amplitude, total_flux
    @property
    def amplitude(self):
        return N.abs(self.wavefront)
    "Amplitude of the wavefront "

    @property
    def intensity(self):
        return N.abs(self.wavefront)**2
    "Intensity of the wavefront"

    @property
    def phase(self):
        return N.angle(self.wavefront)
    "Phase in radians"

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
            _log.debug("  Wavefront and optic %s already at same plane type, no propagation needed." % optic.name)
            return
        else:
            msg = "  Propagating wavefront to %s. " % str(optic)
            _log.debug(msg)
            self.history.append(msg)

        #_log.debug('==== at %s, to %s' % (typestrs[self.planetype], typestrs[optic.planetype]))
        if optic.planetype == ROTATION:     # rotate
            self.rotate(optic.angle)
            self.location='after '+optic.name
        elif optic.planetype == DETECTOR and self.planetype ==PUPIL:    # MFT pupil to detector
            self._propagateMFT(optic)
            self.location='before '+optic.name
        elif optic.planetype == PUPIL and self.planetype ==IMAGE and self._last_transform_type =='MFT': # inverse MFT detector to pupil
            self._propagateMFTinverse(optic)
            self.location='before '+optic.name
        elif self.planetype==IMAGE and optic.planetype == DETECTOR:
            raise NotImplemented('image plane directly to detector propagation (resampling!) not implemented yet')
        else:
            self._propagateFFT(optic)           # FFT pupil to image or image to pupil
            self.location='before '+optic.name

    def _propagateFFT(self, optic):
        """ Propagate from pupil to image or vice versa using a padded FFT """


        if self.oversample > 1 and not self.ispadded: #add padding for oversampling, if necessary
            assert self.oversample == optic.oversample
            self.wavefront = padToOversample(self.wavefront, self.oversample)
            self.ispadded = True
            if optic.verbose: _log.debug("    Padded WF array for oversampling by %dx" % self.oversample)
            self.history.append("    Padded WF array for oversampling by %dx" % self.oversample)

        # Set up for computation - figure out direction & normalization
        if self.planetype == PUPIL and optic.planetype == IMAGE:
            FFT_direction = 'forward'
            normalization_factor = 1./ self.wavefront.shape[0] # correct for numpy fft
            numpy_fft = N.fft.fft2
            numpy_fftshift = N.fft.fftshift
            #(pre-)update state:
            self.planetype=IMAGE
            self.pixelscale = self.wavelength/ self.diam / self.oversample * RADIANStoARCSEC
            self.fov = self.wavefront.shape[0] * self.pixelscale
            self.history.append('   FFT %s,  to IMAGE  scale=%f' %(self.wavefront.shape, self.pixelscale))

        elif self.planetype == IMAGE and optic.planetype ==PUPIL:
            FFT_direction = 'backward'
            if _USE_FFTW3:
                normalization_factor =  1./self.wavefront.shape[0] # correct for FFTW3 FFT
            else:
                normalization_factor =  self.wavefront.shape[0] # correct for numpy fft
            numpy_fft = N.fft.ifft2
            numpy_ifftshift = N.fft.ifftshift
            #(pre-)update state:
            self.planetype=PUPIL
            self.pixelscale = self.diam *self.oversample / self.wavefront.shape[0]
            self.history.append('   FFT %s,  to PUPIL scale=%f' %(self.wavefront.shape, self.pixelscale))


        # do FFT
        if _FLUXCHECK: _log.debug("\tPre-FFT total intensity: "+str(self.totalIntensity))
        if _TIMETESTS: t0 = time.time()

        if FFT_direction =='backward': self.wavefront = numpy_ifftshift(self.wavefront)

        if _USE_FFTW3:
            _log.debug("using FFTW3 FFT")
            # Benchmarking on a Mac Pro (8 cores) indicated that the fastest performance comes from
            # in-place FFTs, and that it is safe to ignore byte alignment issues for these arrays
            # (indeed, even beneficial in many cases) contrary to the suggestion of the FFTW docs
            # which say that aligning arrays helps. Not sure why, but it's true!
            # See the discussion of FFTs in the documentation.
            #wfold = self.copy()
            if (self.wavefront.shape, FFT_direction) not in _FFTW3_INIT.keys():
                # The first time you run FFTW3 to transform a given size, it does a speed test to determine optimal algorithm
                # that is destructive to your chosen array. So only do that test on a copy, not the real array:
                _log.info("Evaluating FFT optimal algorithm for %s, direction=%s" % (str(self.wavefront.shape), FFT_direction))
                fftplan = fftw3.Plan(self.wavefront.copy(), None, nthreads = multiprocessing.cpu_count(),direction=FFT_direction, flags=['measure'])
                _FFTW3_INIT[(self.wavefront.shape, FFT_direction)] = True

            fftplan = fftw3.Plan(self.wavefront, None, nthreads = multiprocessing.cpu_count(),direction=FFT_direction, flags=['measure'])
            fftplan.execute() # execute the plan
                #print "After  FFTW Flux 2: %f" % (abs(outarr)**2).sum()
            # due to FFTW normalization convention, must divide by number of pixels per side.
                #print "After  FFTW Flux 1: %f" % (self.totalIntensity)
        else:
            _log.debug("using numpy FFT")
            self.wavefront = numpy_fft(self.wavefront)

        #wave1 = self.copy()
        #stop()
        if FFT_direction == 'forward':
            self.wavefront = numpy_fftshift(self.wavefront)
            # FFT produces pixel-centered images by default, unless the _image_centered param has already been set by an FQPM_FFT_aligner class
            if self._image_centered != 'corner': self._image_centered = 'pixel'
        self.wavefront = self.wavefront *normalization_factor
        self._last_transform_type = 'FFT'

        if _TIMETESTS:
            t1 = time.time()
            _log.debug("\tTIME %f s\t for the FFT" % (t1-t0))

        if _FLUXCHECK: _log.debug("\tPost-FFT total intensity: "+str(self.totalIntensity))

        #stop()


    def _propagateMFT(self, det):
        """ Compute from pupil to an image using the Soummer et al. 2007 MFT algorithm"""

        assert self.planetype == PUPIL
        assert det.planetype == DETECTOR

        if self.ispadded:
            #pupil plane is padded - trim that out since it's not needed
            self.wavefront = removePadding(self.wavefront, self.oversample)
            self.ispadded = False
        self._preMFT_pupil_shape =self.wavefront.shape  #save for possible inverseMFT
        self._preMFT_pupil_pixelscale = self.pixelscale #save for possible inverseMFT


        # the arguments for the SFT are
        # - wavefront (assumed to fill the input array)
        # - focal plane size in lambda/D units
        # - number of pixels on a side in focal plane array.

        lamD = self.wavelength / self.diam * RADIANStoARCSEC
        #print "lam/D = %f arcsec" % lamD

        det_fov_lamD = det.fov_arcsec / lamD
        det_calc_size_pixels = det.fov_pixels * det.oversample

        mft = SFT.SlowFourierTransform(choice='ADJUSTIBLE', verbose=False)
        msg= '    Propagating w/ MFT: %.4f"/pix     fov=%.3f lam/D    npix=%d' % \
            (det.pixelscale/det.oversample, det_fov_lamD, det_calc_size_pixels)
        _log.debug(msg)
        self.history.append(msg)
        det_offset = det.det_offset if hasattr(det, 'det_offset') else (0,0)

        # det_offset controls how to shift the PSF.
        # it gives the coordinates (X, Y) relative to the exact center of the array
        # for the location of the phase center of a converging perfect spherical wavefront.
        # This is where a perfect PSF would be centered. Of course any tilts, comas, etc, from the OPD
        # will probably shift it off elsewhere for an entirely different reason, too.
        self.wavefront = mft.perform(self.wavefront, det_fov_lamD, det_calc_size_pixels, offset=det_offset)
        self._last_transform_type = 'MFT'

        self.planetype=IMAGE
        self.fov = det.fov_arcsec
        self.pixelscale = det.fov_arcsec / det_calc_size_pixels

    def _propagateMFTinverse(self, pupil, pupil_npix=None):
        """ Compute from an image to a pupil using the Soummer et al. 2007 MFT algorithm
        This allows transformation back from an arbitrarily-sampled 'detector' plane to a pupil. 

        This is only used inside the semi-analytic coronagraphy algorithm.

        """

        assert self.planetype == IMAGE
        assert pupil.planetype == PUPIL

        # the arguments for the SFT are
        # - wavefront (assumed to fill the input array)
        # - focal plane size in lambda/D units
        # - number of pixels on a side in focal plane array.

        lamD = self.wavelength / self.diam * RADIANStoARCSEC
        #print "lam/D = %f arcsec" % lamD

        det_fov_lamD = self.fov / lamD
        #det_calc_size_pixels = det.fov_pixels * det.oversample

        # try to transform to whatever the intrinsic scale of the next pupil is. 
        # but if this ends up being a scalar (meaning it is an AnalyticOptic) then
        # just go back to our own prior shape and pixel scale.
        if pupil_npix == None:
            if pupil.shape is not None and pupil.shape[0] != 1:
                pupil_npix = pupil.shape[0]
            else:
                pupil_npix = self._preMFT_pupil_shape[0]

        mft = SFT.SlowFourierTransform(choice='ADJUSTIBLE', verbose=False)
        msg= '    Propagating w/ InvMFT: %.4f"/pix     fov=%.3f lam/D    pupil npix=%d' % \
            (self.pixelscale, det_fov_lamD, pupil_npix)
        _log.debug(msg)
        self.history.append(msg)
        det_offset = (0,0)  # det_offset not supported for InvMFT

        self.wavefront = mft.inverse(self.wavefront, det_fov_lamD, pupil_npix)
        self._last_transform_type = 'InvMFT'

        self.planetype=PUPIL
        self.pixelscale = self.diam / self.wavefront.shape[0]


    def tilt(self, Xangle=0.0, Yangle=0.0):
        """ Tilt a wavefront in X and Y.

        Recall from Fourier optics (although this is straightforwardly rederivable by drawing triangles...)
        that for a wavefront tilted by some angle theta in radians, for a point r meters from the center of
        the pupil has
            extra_pathlength = sin(theta) * r
            extra_waves = extra_pathlength/ wavelength = r * sin(theta) / wavelength

        So we calculate the U and V arrays (corresponding to r for the pupil, in meters from the center)
        and then multiply by the appropriate trig factors for the angle.

        The sign convention is chosen such that positive Yangle tilts move the star upwards in theg
        array at the focal plane. (This is sort of an inverse of what physically happens in the propagation
        to or through focus, but we're ignoring that here and trying to just work in sky coords)

        Parameters
        ----------
        Xangle, Yangle : float
            tilt angles, specified in arcseconds

        """
        if self.planetype==IMAGE:
            raise NotImplementedError("Are you sure you want to tilt a wavefront in an IMAGE plane?")

        if N.abs(Xangle) > 0 or N.abs(Yangle)>0:
            xangle_rad = Xangle * (N.pi/180/60/60)
            yangle_rad = Yangle * (N.pi/180/60/60)

            npix = self.wavefront.shape[0]
            V, U = N.indices(self.wavefront.shape, dtype=float)
            V -= (npix-1)/2.0
            V *= self.pixelscale
            U -= (npix-1)/2.0
            U *= self.pixelscale

            tiltphasor = N.exp( 2j*N.pi * (U * xangle_rad + V * yangle_rad)/self.wavelength)

        else:
            _log.warn("Wavefront.tilt() called, but requested tilt was zero. No change.")
            tiltphasor = 1.

        #stop()
        #Compute the tilt of the wavefront required to shift it by some amount in the image plane.




        self.wavefront *= tiltphasor
        self.history.append("Tilted wavefront")

    def rotate(self, angle=0.0):
        """Rotate a wavefront by some amount

        Parameters
        ----------
        angle : float
            Angle to rotate, in degrees counterclockwise.

        """
        #self.wavefront = scipy.ndimage.interpolation.rotate(self.wavefront, angle, reshape=False)
        # Huh, the ndimage rotate function does not work for complex numbers. That's weird.
        # so let's treat the real and imaginary parts individually
        # FIXME TODO or would it be better to do this on the amplitude and phase?
        rot_real = scipy.ndimage.interpolation.rotate(self.wavefront.real, angle, reshape=False)
        rot_imag = scipy.ndimage.interpolation.rotate(self.wavefront.imag, angle, reshape=False)
        self.wavefront = rot_real + 1.j*rot_imag

        self.history.append('Rotated by %f degrees, CCW' %(angle))


    def coordinates(self):
        """ Return Y, X coordinates for this wavefront, in the manner of numpy.indices()

        This function knows about the offset resulting from FFTs. Use it whenever computing anything
        measures in wavefront coordinates.

        Returns
        -------
        Y, X :  array_like
            Wavefront coordinates in either meters or arcseconds for pupil and image, respectively

        """
        y, x = N.indices(self.shape, dtype=float)

        # in most cases, the x and y values are centered around the exact center of the array.
        # This is not true in general for FFT-produced image planes where the center is in the
        # middle of one single pixel (the 0th-order term of the FFT), even though that means that the
        # PSF center is slightly offset from the array center.
        # On the other hand, if we used the FQPM FFT Aligner optic, then that forces the PSF center to
        # the exact center of an array.
        if self.planetype == PUPIL:
            y-= (self.shape[0]-1)/2.
            x-= (self.shape[1]-1)/2.
        elif self.planetype == IMAGE:
            # The following are just relevant for the FFT-created images, not for the Detector MFT image at the end.
            if self._last_transform_type == 'FFT':
                # FFT array sizes will always be even, right?
                if self._image_centered=='pixel':  # so this goes to an integer pixel
                    y-= (self.shape[0])/2.
                    x-= (self.shape[1])/2.
                elif self._image_centered=='array_center' or self._image_centered=='corner':  # and this goes to a pixel center
                    y-= (self.shape[0]-1)/2.
                    x-= (self.shape[1]-1)/2.
            else:
                # MFT produced images are always exactly centered.
                y-= (self.shape[0]-1)/2.
                x-= (self.shape[1]-1)/2.


        x *= self.pixelscale
        y *= self.pixelscale
        return y, x

#------
class OpticalElement():
    """ Defines an arbitrary optic, based on amplitude transmission and/or OPD files.

    This optic could be a pupil or field stop, an aberrated mirror, a phase mask, etc.
    The OpticalElement class follows the behavoior of the Wavefront class, using units
    of meters/pixel in pupil space and arcsec/pixel in image space.

    Parameters
    ----------
    name : string
        descriptive name for optic
    tranmission, opd : string
        FITS filenames for the transmission (from 0-1) and opd (in meters)
    opdunits : string
        units for the OPD file. Default is 'meters'.
    verbose : bool
        whether to print stuff while computing
    planetype : int
        either IMAGE or PUPIL
    oversample : int
        how much to oversample beyond Nyquist.
    shift : tuple of floats, optional
        2-tuple containing X and Y fractional shifts for the pupil.
    rotation : float
        Rotation for that optic



    *NOTE:* All mask files must be *squares*.
    """
    #pixelscale = None
    #"float attribute. Pixelscale in arcsec or meters per pixel. Will be 'None' for null or analytic optics."


    def __init__(self, name="unnamed optic", transmission=None, opd= None, verbose=True, planetype=None, oversample=1,opdunits="meters", shift=None, rotation=None):

        self.name = name
        """ string. Descriptive Name of this optic"""
        self.verbose=verbose

        self.opd_file = None
        self.amplitude_file = None
        self.amplitude_header = None
        self.opd_header = None

        self.planetype = planetype      # pupil or image
        self.oversample = oversample    # oversampling factor, none by default
        self.ispadded = False           # are we padded w/ zeros for oversampling the FFT?

        # Read amplitude and/or OPD from disk
        if opd is None and transmission is None:   # no input files, so just make a scalar
            _log.warn("No input files specified. ")
            _log.warn("Creating a null optical element. Are you sure that's what you want to do?")
            self.amplitude = N.asarray([1.])
            self.opd = N.asarray([0.])
            self.pixelscale = None
            self.name = "-empty-"
        else:
            if transmission is not None:        # load transmission file.
                self.amplitude_file = transmission
                self.amplitude, self.amplitude_header = pyfits.getdata(self.amplitude_file, header=True)
                if len (self.amplitude.shape) != 2 or self.amplitude.shape[0] != self.amplitude.shape[1]:
                    raise ValueError, "OPD image must be 2-D and square"
                _log.info(self.name+": Loaded amplitude from "+self.amplitude_file)
            else:                               # else if only OPD set, create an array of 1s with same size.
                opd_shape = pyfits.getdata(opd).shape
                self.amplitude = N.ones(opd_shape)

            if opd is not None:         # Load OPD file.
                # if OPD is specified as a tuple, treat the first element as the filename and 2nd as the slice of a cube.
                if isinstance(opd, basestring):
                    self.opd_file = opd
                    self.opd_slice = 0
                    extra_msg = ""
                else:
                    self.opd_file = opd[0]
                    self.opd_slice = opd[1]
                    extra_msg = ", plane %d " % self.opd_slice

                self.opd, self.opd_header = pyfits.getdata(self.opd_file, header=True)

                if len(self.opd.shape) ==3: # we were given a datacube, so grab just one slice.
                    self.opd = self.opd[self.opd_slice, :,:]

                if len (self.opd.shape) != 2 or self.opd.shape[0] != self.opd.shape[1]:
                    raise ValueError, "OPD image must be 2-D and square"
                msg =  self.name+": Loaded OPD from "+self.opd_file+extra_msg
                _log.info(msg)


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
                _log.info("Requested optic shift of (%6.3f, %6.3f) %%" % (shift))
                _log.info("Actual shift applied   = (%6.3f, %6.3f) %%" % (rollx*1.0/self.amplitude.shape[1], rolly *1.0/ self.amplitude.shape[0]))
                self._shift = (rollx*1.0/self.amplitude.shape[1], rolly *1.0/ self.amplitude.shape[0])

                self.amplitude = N.roll(self.amplitude, rolly, axis=0)
                self.amplitude = N.roll(self.amplitude, rollx, axis=1)
                self.opd       = N.roll(self.opd,       rolly, axis=0)
                self.opd       = N.roll(self.opd,       rollx, axis=1)

            # Likewise, if a rotation is specified and we're NOT a null (scalar) optic, then do the rotation:
            if rotation is not None and len(self.amplitude.shape) ==2:

                # do rotation with interpolation, but try to clean up some of the artifacts afterwards.
                # this is imperfect at best, of course...

                self.amplitude = scipy.ndimage.interpolation.rotate(self.amplitude, rotation, reshape=False).clip(min=0,max=1.0)
                wnoise = N.where(( self.amplitude < 1e-3) & (self.amplitude > 0))
                self.amplitude[wnoise] = 0
                self.opd       = scipy.ndimage.interpolation.rotate(self.opd,       rotation, reshape=False)
                _log.info("  Rotated optic by %f degrees counter clockwise." % rotation)
                #pyfits.PrimaryHDU(self.amplitude).writeto("test_rotated_amp.fits", clobber=True)
                #pyfits.PrimaryHDU(self.opd).writeto("test_rotated_opt.fits", clobber=True)
                self._rotation = rotation


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
        if self.pixelscale is not None and hasattr(wave,'pixelscale') and abs(wave.pixelscale -self.pixelscale)/self.pixelscale >= float_tolerance:
            _log.debug("Pixelscales: %f, %f" % (wave.pixelscale, self.pixelscale))
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
        # note: do not pad the phasor if it's just a scalar!
        if self.planetype == PUPIL and wave.ispadded and self.phasor.size !=1:
            return padToOversample(self.phasor, wave.oversample)
        else:
            return self.phasor

    def display(self, nrows=1, row=1, what='both', phase=False, wavelength=None, crosshairs=True, ax=None, colorbar_orientation=None):
        """Display plots showing an optic's transmission and OPD

        Parameters
        ----------
        what : str
            What do display: 'intensity', 'phase', or 'both'
        ax : matplotlib.Axes
            Axes to display into

        """
        if colorbar_orientation is None:
            colorbar_orientation= "horizontal" if nrows == 1 else 'vertical'
        cmap = matplotlib.cm.gray
        cmap.set_bad('0.0')
        cmap_opd = matplotlib.cm.jet
        cmap_opd.set_bad('0.3')
        norm_amp=matplotlib.colors.Normalize(vmin=0, vmax=1)
        norm_opd=matplotlib.colors.Normalize(vmin=-0.5e-6, vmax=0.5e-6)


        if self.amplitude.shape is () or self.amplitude.size == 1:
            # can't really display a null or scalar optical element?
            # really the best way to do this would be to create a subclass for a
            # scalar optical element...

            #--this code is probably now obsoleted by ScalarElement?? --
            tmp = N.ones((10,10))
            tmp2 = N.ones((10,10))
            ax = p.subplot(nrows, 2, row*2-1)
            tmp[:,:] = self.amplitude

            imshow_with_mouseover(tmp, ax=ax, cmap=cmap, norm=norm_amp )
            ax.set_xticklabels([""]*10)
            ax.set_yticklabels([""]*10)
            p.ylabel(self.name+"\n")
            cb = p.colorbar(ax.images[0], orientation=colorbar_orientation, ticks=[0,0.25, 0.5, 0.75, 1.0])

            ax2 = p.subplot(nrows, 2, row*2)
            tmp2[:,:] = self.opd
            imshow_with_mouseover(tmp2, ax=ax2, cmap=cmap_opd, norm=norm_opd )
            ax2.set_xticklabels([""]*10)
            ax2.set_yticklabels([""]*10)
            cb = p.colorbar(ax2.images[0], orientation=colorbar_orientation, ticks=N.array([-0.5, -0.25, 0, 0.25, 0.5])*1e-6)
            if crosshairs:
                for a in [ax, ax2]:
                    a.axhline(0,ls=":", color='k')
                    a.axvline(0,ls=":", color='k')

            return


        if self.planetype == PUPIL:
            pixelscale = self.pupil_scale
            units = "[meters]"
        else:
            pixelscale = self.pixelscale
            units = "[arcsec]"
        if nrows > 1: units = self.name+"\n"+units


        halfsize = pixelscale*self.amplitude.shape[0]/2
        #extent = [0,pixelscale*self.amplitude.shape[0], 0,pixelscale*self.amplitude.shape[1]]
        extent = [-halfsize, halfsize, -halfsize, halfsize]

        _log.debug("Display pixel scale = %.3f " % pixelscale)

        #ampl = N.ma.masked_equal(self.amplitude, 0)
        ampl = self.amplitude
        #opd= N.ma.masked_array(self.opd, mask=(self.amplitude ==0))
        opd = self.opd.copy()
        opd[N.where(self.amplitude ==0)] = N.nan

        if what=='intensity' or what=='both':
            if ax is None:
                ax = p.subplot(nrows, 2, row*2-1)
            imshow_with_mouseover(ampl, ax=ax, extent=extent, cmap=cmap, norm=norm_amp)
            if nrows == 1:
                p.title("Transmissivity for "+self.name)
            p.ylabel(units)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))
            cb = p.colorbar(ax.images[0], orientation=colorbar_orientation, ticks=[0,0.25, 0.5, 0.75, 1.0])
            cb.set_label('transmission')
            if crosshairs:
                ax.axhline(0,ls=":", color='k')
                ax.axvline(0,ls=":", color='k')


        if what=='phase' or what=='both':
            if ax is None:
                ax2 = p.subplot(nrows, 2, row*2-1)
            else:
                ax2 = ax
    
            ax2 = p.subplot(nrows, 2, row*2)
            imshow_with_mouseover(opd, ax=ax2, extent=extent, cmap=cmap_opd, norm=norm_opd)
            p.ylabel(units)
            ax2.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))
            ax2.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))
            if nrows == 1:
                p.title("OPD for "+self.name)
            cb = p.colorbar(ax2.images[0], orientation=colorbar_orientation, ticks=N.array([-0.5, -0.25, 0, 0.25, 0.5])*1e-6)
            cb.set_label('meters')
            if crosshairs:
                ax2.axhline(0,ls=":", color='k')
                ax2.axvline(0,ls=":", color='k')

    def __str__(self):
        if self.planetype is PUPIL:
            #return "Pupil plane: %s (%dx%d pixels, diam=%f m)" % (self.name, self.shape[0], self.shape[0], self.pupil_diam)
            return "Pupil plane: %s " % (self.name)
        elif self.planetype is IMAGE:
            desc = "(%dx%d pixels, scale=%f arcsec/pixel)" % (self.shape[0], self.shape[0], self.pixelscale) if self.pixelscale is not None else "(Analytic)"
            return "Image plane: %s %s" % (self.name, desc)
        else:
            return "Optic: "+self.name

    @property
    def shape(self):
        if hasattr(self, 'amplitude'):
            return self.amplitude.shape
        else: return None

    def apply(self, wavefront):
        phasor = self.getPhasor(wavefront.wavelength)
        if not N.isscalar(phasor):
            assert wavefront.shape == phasor.shape

class Rotation(OpticalElement):
    """ Performs a rotation of the axes in the optical train.

    This is not an actual optic itself, of course, but can be used to model
    a rotated optic by appling a Rotation before and/or after light is incident
    on that optic.


    This is basically a placeholder to indicate the need for a rotation at a
    given part of the optical train. The actual rotation computation is performed
    in the Wavefront object's propagation routines.


    Parameters
    ----------
    angle : float
        Rotation angle, in degrees counterclockwise.

    """
    def __init__(self, angle=0.0, verbose=True):
        self.name = "Rotation by %.2f degrees" % angle
        self.angle = angle
        self.planetype = ROTATION
        self.verbose =verbose

    def __str__(self):
        return "Rotation by %f degrees counter clockwise" % self.angle

    def getPhasor(self,wave):
        return 1.0  #no change in wavefront (apart from the rotation)
        # returning this is necessary to allow the multiplication in propagate_mono to be OK

    def display(self, nrows=1, row=1, **kwargs):
        p.subplot(nrows, 2, row*2-1)
        p.text(0.3,0.3,self.name)

        #raise NotImplementedError("display is not applicable for a Rotation.")


class AnalyticOpticalElement(OpticalElement):
    """ Defines an abstract analytic optical element, i.e. one definable by some formula rather than
        by an input OPD or pupil file.

        This class is useless on its own; instead use its various subclasses that implement appropriate getPhasor functions.
        It exists mostly to provide some behaviors & initialization common to all analytic optical elements.

        Parameters
        ----------
        name, verbose, oversample, planetype : various
            Same as for OpticalElement
        transmission, opd : string
            These are *not allowed* for Analytic optical elements, and this class will raise an error if you try to set one.


    """
    def __init__(self, name="unnamed", verbose=True, oversample=1, transmission=None, opd=None, planetype=None):
        self.name = name
        self.verbose=verbose
        self.planetype=planetype
        self.oversample = oversample # this is irrelevant for analytic functions, but is included for compatibility with
                                     # the base OpticalElement class.
        self.shape = None # no explicit shape required
        self.pixelscale = None
        if transmission is not None:
            raise ValueError("AnalyticOpticalElements should not have a transmission= argument specified.")
        if opd is not None:
            raise ValueError("AnalyticOpticalElements should not have a opd= argument specified.")

    def __str__(self):
        if self.planetype is PUPIL:
            return "Pupil plane: %s (Analytic)" % (self.name)
        elif self.planetype is IMAGE:
            return "Image plane: %s (Analytic)" % (self.name)
        else:
            return "Optic: "+self.name

    def getPhasor(self,wave):
        raise NotImplementedError("getPhasor must be supplied by a derived subclass")

    def display(self, nrows=1, row=1,  wavelength=2e-6, npix=512, **kwargs):
        "Display an Analytic optic by first computing it onto a grid..."

        if self.planetype is PUPIL:
            unit="meters"
            if hasattr(self, 'pupil_diam'): self.diam = self.pupil_diam
            else : self.diam = 6.5
            halffov = self.diam/2
            w = Wavefront(wavelength=wavelength, npix=npix,  diam = self.diam)
            self.pupil_scale = self.diam/npix
            self.pixelscale = self.pupil_scale

        else:
            unit="arcsec"
            if hasattr(self, '_default_display_size'):
                halffov = self._default_display_size/2
            else:
                halffov = 2.0
            self.pixelscale = 2.*halffov/npix
            w = Wavefront(wavelength=wavelength, npix=npix,  pixelscale = self.pixelscale)


        # set attributes appropriately as if this were a regular OPticalElement
        phasor = self.getPhasor(w)
        self.amplitude = N.abs(phasor)
        phase = N.angle(phasor) * 2*N.pi
        self.opd = phase *wavelength

        #then call parent class display
        OpticalElement.display(self,nrows=nrows, row=row, **kwargs)

        # now un-set everything back cause this is analytic and these are unneeded
        self.pixelscale = None
        self.pupil_scale = None
        self.diam = None
        self.opd = None
        self.amplitude = None

class ScalarTransmission(AnalyticOpticalElement):
    """ Either a null optic (empty plane) or some perfect ND filter...
    But most commonly this is just used as a null optic placeholder """
    def __init__(self, name="-empty-", transmission=1.0, **kwargs):
        AnalyticOpticalElement.__init__(self,name=name, **kwargs)
        self.transmission = transmission
    def getPhasor(self, wave):
        return self.transmission

class InverseTransmission(OpticalElement):
    """ Given any arbitrary OpticalElement with transmission T(x,y)
    return the inverse transmission 1 - T(x,y)

    This is a useful ingredient in the SemiAnalyticCoronagraph algorithm.
    """
    def __init__(self, optic=None):
        if optic is None or not hasattr(optic, 'getPhasor'):
            raise ValueError("Need to supply an valid optic to invert!")
        self.optic = optic
        self.name = "1 - "+optic.name
        self.planetype = optic.planetype
        self.shape = optic.shape
        self.pixelscale = optic.pixelscale
        self.oversample = optic.oversample
    def getPhasor(self, wave):
        return 1- self.optic.getPhasor(wave)


class ThinLens(AnalyticOpticalElement):
    """ An idealized thin lens, implemented as a Zernike defocus term.

    """
    def __init__(self, name='Thin lens', nwaves=4.0, reference_wavelength=2e-6, **kwargs):
        AnalyticOpticalElement.__init__(self,name=name, **kwargs)
        self.planetype=PUPIL

        self.reference_wavelength = reference_wavelength
        self.nwaves = nwaves
        self.max_phase_delay = reference_wavelength * nwaves


    def getPhasor(self, wave):
        y, x = wave.coordinates()
        r = N.sqrt(x**2+y**2)
        # get the normalized radius, assuming the input wave
        # is a square
        #if wave.shape[0] != wave.shape[1]:
            #raise NotImplemented("ThinLens only implemented for square wavefronts.")
        max_r = r[N.where(wave.intensity > 0)].max()
        r_norm = r / max_r 
        #r_norm = r / (wave.shape[0]/2.*wave.pixelscale)

        defocus_zernike = N.sqrt(3)* (2* r_norm**2 - 1)  *  (self.nwaves *self.reference_wavelength/wave.wavelength)
        lens_phasor = N.exp(1.j * 2* N.pi * defocus_zernike)
        #stop()
        return lens_phasor




class BandLimitedCoron(AnalyticOpticalElement):
    """ Defines an ideal band limited coronagraph occulting mask.

        Parameters
        ----------
        name : string
            Descriptive name
        kind : string
            Either 'circular' or 'linear'. The linear ones are custom shaped to NIRCAM's design
            with flat bits on either side of the linear tapered bit.
        sigma : float
            The numerical size parameter, as specified in Krist et al. 2009 SPIE
        wavelength : float
            Wavelength this BLC is optimized for, only for the linear ones.

    """
    def __init__(self, name="unnamed BLC", kind='circular', sigma=1, wavelength=None, **kwargs):
        AnalyticOpticalElement.__init__(self,name=name, **kwargs)
        self.planetype=IMAGE

        self.kind = kind.lower()        # either circular or linear
        if self.kind not in ['circular', 'linear', 'nircamwedge', 'nircamcircular']:
            raise ValueError("Invalid kind of BLC: "+self.kind)
        self.sigma = float(sigma)              # size parameter. See section 2.1 of Krist et al. SPIE 2007, 2009
        if wavelength is not None:
            self.wavelength = float(wavelength)    # wavelength, for selecting the linear wedge option only
        self._default_display_size= 20.        # default size for onscreen display, sized for NIRCam

    def getPhasor(self,wave):
        """ Compute the amplitude transmission appropriate for a BLC for some given pixel spacing
        corresponding to the supplied Wavefront
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("BLC getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == IMAGE)

        y, x = wave.coordinates()
        if self.kind == 'circular':
            # larger sigma implies narrower peak? TBD verify if this is correct
            #
            r = N.sqrt(x**2+y**2)
            sigmar = self.sigma*r
            self.transmission = (1-  (2*scipy.special.jn(1,sigmar)/sigmar)**2)**2
        if self.kind == 'nircamcircular':
            # larger sigma implies narrower peak? TBD verify if this is correct
            #
            r = N.sqrt(x**2+y**2)
            sigmar = self.sigma*r
            sigmar.clip(N.finfo(sigmar.dtype).tiny)  # avoid divide by zero -> NaNs
            self.transmission = (1-  (2*scipy.special.jn(1,sigmar)/sigmar)**2)**2

            # add in the ND squares. Note the positions are not exactly the same in the two wedges.
            # See the figures  in Krist et al. of how the 6 ND squares are spaced among the 5 corongraph regions
            # Also add in the opaque border of the coronagraph mask holder.
            if self.sigma > 4:
                wnd = N.where((y > 5) & (((x < -5)&(x>-10)) | ((x > 7.5)&(x<12.5))) ) # MASK210R has one in the corner and one half in the other corner
                wborder  = N.where((N.abs(y) > 10) | (x < -10) )                        #left end of mask holder
            else:
                wnd = N.where((y > 5) & (N.abs(x) > 7.5) & (N.abs(x) < 12.5))        # the others have two halves on in each corner.
                wborder  = N.where(N.abs(y) > 10)

            self.transmission[wnd] = 1e-3
            self.transmission[wborder] = 0


        elif self.kind == 'linear':
            #raise(NotImplemented("Generic linear not implemented"))
            sigmar = self.sigma * N.abs(y)
            sigmar.clip( N.finfo(sigmar.dtype).tiny)  # avoid divide by zero -> NaNs
            self.transmission = (1-  (N.sin(sigmar)/sigmar)**2)**2


        elif self.kind == 'nircamwedge':
            # This is hard-coded to the wedge-plus-flat-regions shape for NIRCAM

            # we want a scale factor that goes from  2 to 6 with 1/5th of it as a fixed part on either end
            #scalefact = N.linspace(1,7, x.shape[1]).clip(2,6)

            # the scale fact should depent on X coord in arcsec, scaling across a 20 arcsec FOV.
            # map flat regions to 2.5 arcsec each?
            # map -7.5 to 2, +7.5 to 6. slope is 4/15, offset is +9.5
            scalefact = (2+(-x+7.5)*4/15).clip(2,6)

            #scalefact *= self.sigma / 2 #;2.2513
            #scalefact *= 2.2513
            #scalefact.shape = (1, x.shape[1])
            # This does not work - shape appears to be curved not linear.
            # This is NOT a linear relationship. See calc_blc_wedge in test_poppy.

            if N.abs(self.wavelength - 2.1e-6) < 0.1e-6:
                polyfitcoeffs = N.array([  2.01210737e-04,  -7.18758337e-03,   1.12381516e-01,
                                          -1.00877701e+00,   5.72538509e+00,  -2.12943497e+01,
                                           5.18745152e+01,  -7.97815606e+01,   7.02728734e+01])
            elif N.abs(self.wavelength - 4.6e-6) < 0.1e-6:
                polyfitcoeffs = N.array([  9.16195583e-05,  -3.27354831e-03,   5.11960734e-02,
                                          -4.59674047e-01,   2.60963397e+00,  -9.70881273e+00,
                                           2.36585911e+01,  -3.63978587e+01,   3.20703511e+01])
            else:
                raise NotImplemented("No defined NIRCam wedge BLC mask for that wavelength?  ")

            sigmas = scipy.poly1d(polyfitcoeffs)(scalefact)

            sigmar = sigmas * N.abs(y)
            sigmar.clip( N.finfo(sigmar.dtype).tiny)  # avoid divide by zero -> NaNs
            self.transmission = (1-  (N.sin(sigmar)/sigmar)**2)**2
            # the bar should truncate at +- 10 arcsec:
            woutside = N.where(N.abs(x) > 10)
            self.transmission[woutside] = 1.0


            # add in the ND squares. Note the positions are not exactly the same in the two wedges.
            # See the figures  in Krist et al. of how the 6 ND squares are spaced among the 5 corongraph regions
            # Also add in the opaque border of the coronagraph mask holder.
            if N.abs(self.wavelength - 2.1e-6) < 0.1e-6:
                wnd = N.where((y > 5) & (((x < -5)&(x>-10)) | ((x > 7.5)&(x<12.5))) ) # half ND square on each side
                wborder  = N.where(N.abs(y) > 10)
            elif N.abs(self.wavelength - 4.6e-6) < 0.1e-6:
                wnd = N.where((y > 5) & (((x < -7.5)&(x>-12.5)) | (x > 5)) )
                wborder  = N.where((N.abs(y) > 10) | (x > 10) )     # right end of mask holder

            self.transmission[wnd] = 1e-3
            self.transmission[wborder] = 0



        if not N.isfinite(self.transmission.sum()):
            #print "There are NaNs in the BLC mask - error!"
            #stop()
            _log.debug("There are NaNs in the BLC mask - correcting to zero. (DEBUG LATER?)")
            self.transmission[N.where(N.isfinite(self.transmission) == False)] = 0
        return self.transmission


class FQPM_FFT_aligner(AnalyticOpticalElement):
    """  Helper class for modeling FQPMs accurately

    Adds (or removes) a slight wavelength- and pixel-scale-dependent tilt
    to a pupil wavefront, to ensure the correct alignment of the image plane
    FFT'ed PSF with the desired quad pixel alignment for the FQPM.

    This is purely a computational convenience tool to work around the
    pixel coordinate restrictions imposed by the FFT algorithm,
    not a representation of any physical optic.

    Parameters
    ----------
    direction : string
        'forward' or 'backward'

    """
    def __init__(self, name="FQPM FFT aligner", direction='forward', **kwargs):
        AnalyticOpticalElement.__init__(self, name=name, **kwargs)
        direction = direction.lower()
        if direction != 'forward' and direction !='backward': raise ValueError("Invalid direction %s, must be either forward or backward." % direction)
        self.direction = direction
        #self.displayable = False

    def getPhasor(self,wave):
        """ Compute the required tilt needed to get the PSF centered on the corner between
        the 4 central pixels, not on the central pixel itself.
        """

        if not isinstance(wave, Wavefront):
            raise ValueError("FQPM getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == PUPIL )


        fft_im_pixelscale = wave.wavelength/ wave.diam / wave.oversample * RADIANStoARCSEC
        required_offset = -fft_im_pixelscale *0.5
        if self.direction == 'backward':
            required_offset *= -1
            wave._image_centered='pixel'
        else:
            wave._image_centered='corner'
        wave.tilt(required_offset, required_offset)

        # gotta return something... so return a value that will not affect the wave any more.
        align_phasor = 1.0
        return align_phasor

class IdealFQPM(AnalyticOpticalElement):
    """ Defines an ideal 4-quadrant phase mask coronagraph, with its retardance
    set perfectly to 0.5 waves at one specific wavelength and varying linearly on
    either side of that.

    Parameters
    ----------
    name : string
        Descriptive name
    wavelength : float
        Wavelength in meters for which the FQPM was designed


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

        #outFITS = pyfits.HDUList(pyfits.PrimaryHDU(retardance))
        #outFITS.writeto('retardance_fqpm.fits', clobber=True)
        #_log.info("Retardance is %f waves" % retardance.max())
        FQPM_phasor = N.exp(1.j * 2* N.pi * retardance)
        return FQPM_phasor

class IdealFieldStop(AnalyticOpticalElement):
    """ Defines an ideal square field stop

    Parameters
    ----------
    name : string
        Descriptive name
    size : float
        Size of the field stop, in arcseconds. Default 20.
    angle : float
        Position angle of the field stop sides relative to the detector +Y direction, in degrees.

    """

    def __init__(self, name="unnamed field stop",  size=20., angle=0, **kwargs):
        AnalyticOpticalElement.__init__(self,**kwargs)
        self.name = name
        self.size = size            # size of square stop in arcseconds.
        #self.pixelscale=0
        self.angle=angle
        self._default_display_size = size*1.2

    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the field stop.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("IdealFieldStop getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == IMAGE)

        #phasor = N.zeros(wave.shape)
        y, x= wave.coordinates()
        #y, x = N.indices(wave.shape)
        #y -= wave.shape[0]/2
        #x -= wave.shape[1]/2
        xnew =  x*N.cos(N.deg2rad(self.angle)) + y*N.sin(N.deg2rad(self.angle))
        ynew = -x*N.sin(N.deg2rad(self.angle)) + y*N.cos(N.deg2rad(self.angle))
        x,y = xnew, ynew
        #del xnew
        #del ynew


        halfsize = self.size   / 2
        w_outside = N.where( (abs(y) > halfsize)  | (abs(x) > halfsize))
        del x # for large arrays, cleanup very promptly, before allocating self.transmission
        del y
        self.transmission = N.ones(wave.shape)
        self.transmission[w_outside] = 0

        return self.transmission

class IdealCircularOcculter(AnalyticOpticalElement):
    """ Defines an ideal circular occulter (opaque circle)

    Parameters
    ----------
    name : string
        Descriptive name
    radius : float
        Radius of the occulting spot, in arcseconds. Default is 1.0

    """

    def __init__(self, name="unnamed occulter",  radius=1.0, **kwargs):
        AnalyticOpticalElement.__init__(self,**kwargs)
        self.name = name
        self.radius = radius    # radius of circular occulter in arcseconds.
        self._default_display_size = 10
        #self.pixelscale=0

    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the occulter.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == IMAGE)

        y, x= wave.coordinates()
        #y, x = N.indices(wave.shape)
        #y -= wave.shape[0]/2
        #x -= wave.shape[1]/2
        r = N.sqrt(x**2+y**2) #* wave.pixelscale
        w_inside = N.where( r <= self.radius)

        del x
        del y
        del r
        self.transmission = N.ones(wave.shape)
        self.transmission[w_inside] = 0

        return self.transmission

class IdealBarOcculter(AnalyticOpticalElement):
    """ Defines an ideal bar occulter (like in MIRI's Lyot coronagraph)

    Parameters
    ----------
    name : string
        Descriptive name
    width : float
        width of the bar stop, in arcseconds. Default is 1.0
    angle : float
        position angle of the bar, rotated relative to the normal +y direction.

    """

    def __init__(self, name="bar occulter",  width=1.0, angle= 0, **kwargs):
        AnalyticOpticalElement.__init__(self,**kwargs)
        self.name = name
        self.width = width
        self.angle = angle
        #self.pixelscale=0
        self._default_display_size = 10

    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the occulter.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == IMAGE)

        y, x= wave.coordinates()
        #y, x = N.indices(wave.shape)
        #y -= wave.shape[0]/2
        #x -= wave.shape[1]/2

        xnew = x*N.cos(N.deg2rad(self.angle)) + y*N.sin(N.deg2rad(self.angle))
        del x
        del y
        w_inside = N.where( N.abs(xnew) <= self.width/2)
        del xnew
        self.transmission = N.ones(wave.shape)
        self.transmission[w_inside] = 0

        return self.transmission

class ParityTestAperture(AnalyticOpticalElement):
    """ Defines a circular pupil aperture with boxes cut out.
    This is mostly a test aperture

    Parameters
    ----------
    name : string
        Descriptive name
    radius : float
        Radius of the pupil, in meters. Default is 1.0

    pad_factor : float, optional
        Amount to oversize the wavefront array relative to this pupil.
        This is in practice not very useful, but it provides a straightforward way
        of verifying during code testing that the amount of padding (or size of the circle)
        does not make any numerical difference in the final result.

    """

    def __init__(self, name=None,  radius=1.0, pad_factor = 1.5, **kwargs):
        if name is None: name = "Circle, radius=%.2f m" % radius
        AnalyticOpticalElement.__init__(self,name=name, **kwargs)
        self.radius = radius
        self.pupil_diam = pad_factor * 2* self.radius # for creating input wavefronts - let's pad a bit


    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the occulter.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("CircularAperture getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == PUPIL)

        y, x = wave.coordinates()
        r = N.sqrt(x**2+y**2) #* wave.pixelscale
        del x
        del y


        w_outside = N.where( r > self.radius)
        del r
        self.transmission = N.ones(wave.shape)
        self.transmission[w_outside] = 0

        w_box1 = N.where( r> (self.radius*0.5) & N.abs(x) < self.radius*0.1 & y < 0)
        w_box1 = N.where( r> (self.radius*0.75) & N.abs(y) < self.radius*0.2 & x< 0)
        self.transmission[w_box1] = 0
        self.transmission[w_box2] = 0

        return self.transmission


class CircularAperture(AnalyticOpticalElement):
    """ Defines an ideal circular pupil aperture

    Parameters
    ----------
    name : string
        Descriptive name
    radius : float
        Radius of the pupil, in meters. Default is 1.0

    pad_factor : float, optional
        Amount to oversize the wavefront array relative to this pupil.
        This is in practice not very useful, but it provides a straightforward way
        of verifying during code testing that the amount of padding (or size of the circle)
        does not make any numerical difference in the final result.

    """

    def __init__(self, name=None,  radius=1.0, pad_factor = 1.5, **kwargs):
        if name is None: name = "Circle, radius=%.2f m" % radius
        AnalyticOpticalElement.__init__(self,name=name, **kwargs)
        self.radius = radius
        self.pupil_diam = pad_factor * 2* self.radius # for creating input wavefronts - let's pad a bit


    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the occulter.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("CircularAperture getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == PUPIL)

        y, x = wave.coordinates()
        #y, x = N.indices(wave.shape)
        #y -= wave.shape[0]/2
        #x -= wave.shape[1]/2
        r = N.sqrt(x**2+y**2) #* wave.pixelscale
        del x
        del y


        w_outside = N.where( r > self.radius)
        del r
        self.transmission = N.ones(wave.shape)
        self.transmission[w_outside] = 0

        #_log.debug("Circle: radius = %.3f m = %.3f pixels " % (self.radius, self.radius/wave.pixelscale))
        #_log.debug("        center = %f" % (wave.shape[0]/2 ))
        #w_inside = N.where( r <= self.radius)
        #_log.debug("        minx: %d" % w_inside[0].min())
        #stop()

        return self.transmission

class HexagonAperture(AnalyticOpticalElement):
    """ Defines an ideal hexagonal pupil aperture

    Specify either the side length (= corner radius) or the
    flat-to-flat distance.

    Parameters
    ----------
    name : string
        Descriptive name
    flattoflat : float, optional
        Flat-to-flat distance of the pupil, in meters. Default is 1.0
    side : float, optional
        side length of hexagon, in meters.

    """

    def __init__(self, name=None,  flattoflat=None, side=None, **kwargs):
        if flattoflat is None and side is None:
            self.side = 1.0
        elif side is not None:
            self.side = float(side)
        else:
            self.side = float(flattoflat)/N.sqrt(3.)
        self.pupil_diam = 2* self.side # for creating input wavefronts
        if name is None: name = "Hexagon, side length= %.1f m" % self.side

        AnalyticOpticalElement.__init__(self,name=name, **kwargs)


    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the occulter.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("CircularAperture getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == PUPIL)

        y, x = wave.coordinates()
        absy = N.abs(y)

        self.transmission = N.zeros(wave.shape)

        w_rect = N.where((N.abs(x) <= 0.5*self.side) & (absy <= N.sqrt(3)/2*self.side))
        w_left_tri  = N.where(   (x<=-0.5*self.side) & (x >= -1*self.side) & (absy <= (x+1*self.side)*N.sqrt(3)))
        w_right_tri = N.where(   (x>= 0.5*self.side) & (x <=  1*self.side) & (absy <= (1*self.side-x)*N.sqrt(3)))
        self.transmission[w_rect] = 1
        self.transmission[w_left_tri] = 1
        self.transmission[w_right_tri] = 1

        return self.transmission


class SquareAperture(AnalyticOpticalElement):
    """ Defines an ideal square pupil aperture

    Parameters
    ----------
    name : string
        Descriptive name
    size: float
        half width of the square, in meters. Default is 1.0

    """

    def __init__(self, name=None,  size=1.0, **kwargs):
        if name is None: name = "Square, side length= %.1f m" % size*2
        AnalyticOpticalElement.__init__(self,name=name,**kwargs)
        self.size = size
        self.pupil_diam = 2* self.size # for creating input wavefronts

    def getPhasor(self,wave):
        """ Compute the transmission inside/outside of the occulter.
        """
        if not isinstance(wave, Wavefront):
            raise ValueError("getPhasor must be called with a Wavefront to define the spacing")
        assert (wave.planetype == PUPIL)

        y, x = wave.coordinates()

        w_outside = N.where( (abs(y) > self.size) | (abs(x) > self.size) )
        del y
        del x

        self.transmission = N.ones(wave.shape)
        self.transmission[w_outside] = 0
        return self.transmission


class CompoundAnalyticOptic(AnalyticOpticalElement):
    """ Define a compound analytic optical element made up of the combination
    of two or more individual optical elements.

    This is just a convenience routine for semantic organization of optics.
    It can be useful to keep the list of optical planes cleaner, but
    you can certainly just add a whole bunch of planes all in a row without
    using this class to group them.

    All optics should be of the same plane type (pupil or image); propagation between
    different optics contained inside one compound is not supported.

    Parameters
    ----------
    opticslist : list
        A list of AnalyticOpticalElements to be merged together.

    """
    def __init__(self, name="unnamed", verbose=True, oversample=1, transmission=None, opd=None, planetype=None,
        opticslist=None):
        AnalyticOpticalElement.__init__(self,name=name, verbose=verbose, oversample=oversample, transmission=transmission,opd=opd, planetype=planetype)

        #self.operation = operation
        self.opticslist = []
        self._default_display_size = 3
        for optic in opticslist:
            if not isinstance(optic, AnalyticOpticalElement):
                raise ValueError("Supplied optics list to CompoundAnalyticOptic can only contain AnalyticOptics")
            else:
                if hasattr(optic, '_default_display_size'):
                    self._default_display_size = max(self._default_display_size, optic._default_display_size)
                self.opticslist.append(optic)

    def getPhasor(self,wave):
        #phasor = self.opticslist[0].getPhasor(wave)
        for optic in self.opticslist[:-1]:
            #nextphasor = optic.getPhasor(wave)
            #phasor *= nextphasor #FIXME this does not work... for instance if you have an aperture mask (so all the phase is zero)
                                  # then you really want to multiply the amplitude transmissions and add the phases.
                                  # simpler to just multiply the wave instead here:
            wave *= optic
        # and just hand back the last one to the calling routine:
        return self.opticslist[-1].getPhasor(wave)


class Detector(OpticalElement):
    """ A Detector is a specialized type of OpticalElement that forces a wavefront
    onto a specific fixed pixelization.

    Parameters
    ----------
    name : string
        Descriptive name
    pixelscale : float
        Pixel scale in arcsec/pixel
    fov_pixels, fov_arcsec : float
        The field of view may be specified either in arcseconds or by a number of pixels. Either is acceptable
        and the pixel scale is used to convert as needed.
    oversample : int
        Oversampling factor beyond the detector pixel scale
    offset : tuple (X,Y)
        Offset for the detector center relative to a hypothetical off-axis PSF. Specifying this lets you
        pick a different sub-region for the detector to compute, if for some reason you are computing a small
        subarray around an off-axis source. (Has not been tested!)

    """
    def __init__(self, pixelscale, fov_pixels=None, fov_arcsec=None, oversample=1, name="Detector", offset=None):
        self.name=name
        self.planetype = DETECTOR
        self.pixelscale = float(pixelscale)
        self.oversample = oversample

        if fov_pixels is None and fov_arcsec is None:
            raise ValueError("Either fov_pixels or fov_arcsec must be specified!")
        elif fov_pixels is not None:
            self.fov_pixels = int(fov_pixels)
            self.fov_arcsec = self.fov_pixels * self.pixelscale
        else:
            # set field of view to closest value possible to requested,
            # consistent with having an integer number of pixels
            self.fov_pixels = round(fov_arcsec / self.pixelscale)
            self.fov_arcsec = self.fov_pixels * self.pixelscale
        if hasattr(self.fov_pixels, '__getitem__'):
            self.shape = self.fov_pixels[0:2]
        else:
            self.shape = (self.fov_pixels, self.fov_pixels)

        self.amplitude = 1
        self.opd = 0

    def __str__(self):
        return "Detector plane: %s (%dx%d, %f arcsec/pixel)" % (self.name, self.shape[1], self.shape[0], self.pixelscale)

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
    name : string
        descriptive name of optical system
    oversample : int
        Either how many times *above* Nyquist we should be
        (for pupil or image planes), or how many times a fixed
        detector pixel will be sampled. E.g. `oversample=2` means
        image plane sampling lambda/4*D (twice Nyquist) and
        detector plane sampling 2x2 computed pixels per real detector
        pixel.  Default is 2.
    verbose : bool
        whether to print stuff while computing




    """
    def __init__(self, name="unnamed system", verbose=True, oversample=2):
        self.name = name
        self.verbose=verbose
        self.planes = []                    # List of OpticalElements
        self.oversample = oversample

        self.source_offset_r = 0 # = N.zeros((2))     # off-axis tilt of the source, in ARCSEC
        self.source_offset_theta = 0 # in degrees CCW

        self.intermediate_wfs = None        #
        if self.verbose:
            _log.info("Initialized OpticalSystem: "+self.name)

    # Methods for adding or manipulating optical planes:

    def addPupil(self, optic=None, function=None, **kwargs):
        """ Add a pupil plane optic from file(s) giving transmission or OPD

          1) from file(s) giving transmission or OPD
                [set arguments `transmission=filename` and/or `opd=filename`]
          2) from an analytic function
                [set `function='Circle', 'Square'`
                and set additional kwargs to define shape etc.
          3) from an already-created OpticalElement object
                [set `optic=that object`]

        Parameters
        ----------
        optic : poppy.OpticalElement, optional
            An already-created OpticalElement you would like to add
        function: string, optional
            Name of some analytic function to add.
            Optional `kwargs` can be used to set the parameters of that function.
            Allowable function names are Circle, Square
        opd, transmission : string, optional
            Filenames of FITS files describing the desired optic.


        Note: Now you can use the optic arguement for either an OpticalElement or a string function name,
        and it will do the right thing depending on type.  Both existing arguments are left for back compatibility for now.



        Any provided parameters are passed to :ref:`OpticalElement`.


        """
        if optic is None and function is not None: # support legacy 'function' input
            optic = function

        if isinstance(optic, OpticalElement):
            optic.planetype = PUPIL
        elif isinstance(optic, str):
            if optic.lower() == 'circle':
                fn = CircularAperture
            elif optic.lower() == 'square':
                fn = SquareAperture
            elif optic.lower() == 'hexagon' or optic.lower()=='hex':
                fn = HexagonAperture
            elif optic.lower()[0:4] == 'fqpm':
                fn = FQPM_FFT_aligner

            else:
                raise ValueError("Unknown pupil function type: %s. Perhaps you meant to set transmission= or opd= instead?" % optic)
            optic = fn(planetype=PUPIL, oversample=self.oversample, **kwargs)
        elif optic is None: # create image from files specified in kwargs
            optic = OpticalElement(planetype=PUPIL, oversample=self.oversample, **kwargs)
        else:
            raise TypeError("Not sure how to handle an Optic input of that type...")

        self.planes.append(optic)
        if self.verbose: _log.info("Added pupil plane: "+self.planes[-1].name)

    def addImage(self, optic=None, function=None, **kwargs):
        """ Add an image plane optic, either

          1) from file(s) giving transmission or OPD
                [set arguments `transmission=filename` and/or `opd=filename`]
          2) from an analytic function
                [set `function='circle, fieldstop, bandlimitedcoron, or FQPM'`
                and set additional kwargs to define shape etc.
          3) from an already-created OpticalElement object
                [set `optic=that object`]

        Parameters
        ----------
        optic : poppy.OpticalElement
            An already-created OpticalElement you would like to add
        function: string
            Name of some analytic function to add.
            Optional `kwargs` can be used to set the parameters of that function.
            Allowable function names are CircularOcculter, fieldstop, BandLimitedCoron, FQPM
        opd, transmission : string
            Filenames of FITS files describing the desired optic.


        Note: Now you can use the optic arguement for either an OpticalElement or a string function name,
        and it will do the right thing depending on type.  Both existing arguments are left for back compatibility for now.



        """

        if isinstance(optic, str):
            function = optic
            optic = None

        if optic is None:
            if function == 'CircularOcculter':
                fn = IdealCircularOcculter
            elif function == 'BarOcculter':
                fn = IdealBarOcculter
            elif function == 'fieldstop':
                fn = IdealFieldStop
            elif function == 'BandLimitedCoron':
                fn = BandLimitedCoron
            elif function == 'FQPM':
                fn = IdealFQPM
            elif function is not None:
                raise ValueError("Analytic mask type '%s' is unknown." % function)
            else: # create image from files specified in kwargs
                fn = OpticalElement

            optic = fn(planetype=IMAGE, oversample=self.oversample, **kwargs)
        else:
            optic.planetype = IMAGE
            optic.oversample = self.oversample # these need to match...

        self.planes.append(optic)
        if self.verbose: _log.info("Added image plane: "+self.planes[-1].name)

    def addRotation(self, angle=0):
        self.planes.append(Rotation(angle=angle))
        if self.verbose: _log.info("Added rotation plane: "+self.planes[-1].name)


    def addDetector(self, pixelscale, oversample=None, **kwargs):
        """ Add a Detector object to an optical system.
        By default, use the same oversampling as the rest of the optical system,
        but the user can override to a different value if desired by setting `oversample`.


        Other arguments are passed to the init method for Detector().

        Parameters
        ----------
        pixelscale : float
            Pixel scale in arcsec/pixel
        oversample : int, optional
            Oversampling factor for *this detector*, relative to hardware pixel size.
            Optionally distinct from the default oversampling parameter of the OpticalSystem.

        """

        if oversample is None:
            oversample = self.oversample
        self.planes.append(Detector(pixelscale, oversample=oversample, **kwargs))
        if self.verbose: _log.info("Added detector: "+self.planes[-1].name)


        #return "Optical system '%s' containing %d optics" % (self.name, len(self.planes))

    def list(self):
        print str(self)+"\n\t"+ "\n\t".join([str(p) for p in self.planes])

    def __getitem__(self, num):
        return self.planes[num]

    # methods for dealing with wavefronts:
    def inputWavefront(self, wavelength=2e-6):
        """Create a Wavefront object suitable for sending through a given optical system, based on
        the size of the first optical plane, assumed to be a pupil.

        If the first optical element is an Analytic pupil (i.e. has no pixel scale) then
        an array of 1024x1024 will be created (not including oversampling).

        Uses self.source_offset to assign an off-axis tilt, if requested.

        Parameters
        ----------
        wavelength : float
            Wavelength in meters

        Returns
        -------
        wavefront : poppy.Wavefront instance
            A wavefront appropriate for passing through this optical system.

        """

        npix = self.planes[0].shape[0] if self.planes[0].shape is not None else 1024
        diam = self.planes[0].pupil_diam if hasattr(self.planes[0], 'pupil_diam') else 8

        inwave = Wavefront(wavelength=wavelength,
                npix = npix,
                diam = diam,
                oversample=self.oversample)

        if N.abs(self.source_offset_r) > 0:
            offset_x = self.source_offset_r *-N.sin(self.source_offset_theta*N.pi/180)  # convert to offset X,Y in arcsec
            offset_y = self.source_offset_r * N.cos(self.source_offset_theta*N.pi/180)  # using the usual astronomical angle convention
            inwave.tilt(Xangle=offset_x, Yangle=offset_y)
        return inwave

    def propagate_mono(self, wavelength=2e-6, normalize='first', save_intermediates=False, display_intermediates=False, intermediate_fn='wave_step_%03d.fits', poly_weight=None):
        """ Propagate a wavefront through some number of optics.
        Returns a pyfits.HDUList object.

        Parameters
        ----------
        wavelength : float
            Wavelength in meters
        normalize : string, {'first', 'last'}
            how to normalize the wavefront?
            * 'first' = set total flux = 1 after the first optic, presumably a pupil
            * 'last' = set total flux = 1 after the entire optical system.
            * 'first=2' = set total flux = 2 after the first optic (used for debugging only)
        display_intermediates : bool
            Should intermediate steps in the calculation be displayed on screen? Default False
        save_intermediates : bool
            Should intermediate steps in the calculation be saved to disk? Default False.
            If this is True, then setting `poly_weight` controls whether intermediate optical planes are actually saved to *disk* by this routine
            (for the monochromatic case) or are passed back up via memory and handled in `calcPSF` (for the polychromatic case).
        poly_weight : float
            is this being called as part of a polychromatic calculation?
            if not, set this to None. if so, set this to the weight for
            that wavelength. (This is used only for properly normalizing the
            multiwavelength FITs file written to disk if save_intermediates is set.)


        """

        if _TIMETESTS:
            t_start = time.time()
        if self.verbose: _log.info(" Propagating wavelength = %g meters" % wavelength)
        wavefront = self.inputWavefront(wavelength)

        if save_intermediates and poly_weight is None:
            self.intermediate_wfs=[]
            _log.debug("reset intermediates")

        # do the propagation:

        count = 0
        for optic in self.planes:
            wavefront.propagateTo(optic)
            wavefront *= optic
            count += 1


            if normalize.lower()=='first' and count == 1:
                wavefront.normalize()

            if normalize.lower()=='first=2' and count == 1:
                wavefront.normalize()
                wavefront *= N.sqrt(2) 



            if _FLUXCHECK: _log.debug("  Flux === "+str(wavefront.totalIntensity))

            if save_intermediates:
                if len(self.intermediate_wfs) < count:
                    self.intermediate_wfs.append(wavefront.copy())
                else:
                    self.intermediate_wfs[count-1] += wavefront.copy()*poly_weight
                if poly_weight is None: self.intermediate_wfs[count-1].writeto(intermediate_fn % count, what='parts')
            if display_intermediates:
                if _TIMETESTS:
                    t0 = time.time()
                #if save_intermediates: self.intermediate_wfs[count-1].display(what='best',nrows=len(self.planes),row=count)
                wavefront.display(what='best',nrows=len(self.planes),row=count, colorbar=False)
                if _TIMETESTS:
                    t1 = time.time()
                    _log.debug("\tTIME %f s\t for displaying the wavefront." % (t1-t0))


        # prepare output arrays
        if normalize.lower()=='last':
                wavefront.normalize()

        if _TIMETESTS:
            t_stop = time.time()
            _log.debug("\tTIME %f s\tfor propagating one wavelength" % (t_stop-t_start))

        return wavefront.asFITS()

    def calcPSFmultiproc(self, source, nprocesses=4, save_intermediates=False, **kwargs):
        """Calculate a multi-wavelength PSF over some weighted
        sum of wavelengths, across multiple processors

        This version uses Python's `multiprocessing` package to span tasks across
        available processor cores.

        Any additional `kwargs` will be passed on to `propagate_mono()`


        Notes
        -----

        Don't set the number of processes too high, or your machine will exhaust its memory and grind to a halt.
        Figure about 600 MB per process, due to the large arrays used in the FFTing, and multiple copies etc.
        Yes, this is surprisingly high overhead - but the raw wavefront array is typically 2048^2 * 16
        (since it's type dcomplex) = 67 MB, so this is just ~ 8-10 copies of the array floating arround.
        TODO: be more clever and efficient with all this somehow?
        TODO: write an auto tool to optimize the number of processes automatically?

        Parameters
        ----------
        source : dict
            a dict containing 'wavelengths' and 'weights' list.
            *TBD - replace w/ pysynphot observation object*
        nprocesses : int
            Number of processes. Don't make this too large, or you will overfill your memory and
            things will grind painfully to a half as everything swaps to disk and back.
        save_intermediates : bool
            whether to output intermediate optical planes to disk. Default is False


        Returns
        -------
        outfits :
            a pyfits.HDUList
        """

        if _USE_FFTW3:
            _log.warn('IMPORTANT WARNING: Python multiprocessing and fftw3 do not appear to play well together. This is likely to crash intermittently')
            _log.warn('   We suggest you set   poppy._USE_FFTW3 = False   if you want to use calcPSFmultiproc().')


        if save_intermediates:
            raise NotImplementedError("Can't save intermediate steps if using parallelized code")
        self.intermediate_wfs = []
            #print 'reset intermediates in calcPSF'

        # loop over wavelengths
        if self.verbose: _log.info("Calculating PSF with %d wavelengths, using multiprocessing across %d processes" % (len(source['wavelengths']), _MULTIPROC_NPROCESS  ))
        outFITS = None

        normwts =  N.asarray(source['weights'], dtype=float)
        normwts /= normwts.sum()

        if len(tuple(source['wavelengths'])) != len(tuple(source['weights'])):
            raise ValueError("Input source has different number of weights and wavelengths...")
        # do *NOT* just blindly try to create as many processes as one has CPUs, or one per wavelength either
        # This is a memory-intensive task so that can end up swapping to disk and thrashing IO
        pool = multiprocessing.Pool(nprocesses )

        # build a single iterable containing the required function arguments
        _log.info("Beginning multiprocessor job")
        iterable = [(self, wavelen, weight, kwargs, _USE_FFTW3) for wavelen, weight in zip(source['wavelengths'], normwts)]
        results = pool.map(_wrap_propagate_for_multiprocessing, iterable)
        _log.info("Finished multiprocessor job")
        pool.close()

        # Sum all the results up into one array, using the weights
        outFITS = results[0]
        outFITS[0].data *= normwts[0]
        _log.info("got results for wavelength channel %d / %d" % (0, len(tuple(source['wavelengths']))) )
        for i in range(1, len(normwts)):
            _log.info("got results for wavelength channel %d / %d" % (i, len(tuple(source['wavelengths']))) )
            outFITS[0].data += results[i][0].data * normwts[i]

        # build output FITS header
        waves = N.asarray(source['wavelengths'])
        wts = N.asarray(source['weights'])
        mnwave = (waves*wts).sum() / wts.sum()
        outFITS[0].header.update('WAVELEN', mnwave, 'Weighted mean wavelength in meters')
        outFITS[0].header.update('NWAVES',waves.size, 'Number of wavelengths used in calculation')
        for i in range(waves.size):
            outFITS[0].header.update('WAVE'+str(i), waves[i], "Wavelength "+str(i))
            outFITS[0].header.update('WGHT'+str(i), wts[i], "Wavelength weight "+str(i))
        outFITS[0].header.add_history("Multiwavelength PSF calc on %d processors completed." % nprocesses)

        if self.verbose: _log.info(" PSF Calculation completed.")
        return outFITS

    def calcPSF(self, wavelength=1e-6, weight=None,
        save_intermediates=False, save_intermediates_what='all', display= False, return_intermediates=False, **kwargs):
        """Calculate a PSF, either
        - multi-wavelength PSF over some weighted sum of wavelengths (if you provide a source parameter)
        - monochromatic (if you provide just a wavelen= parameter)

        Any additional `kwargs` will be passed on to `propagate_mono()`

        Parameters
        ----------
        source : dict
            a dict containing 'wavelengths' and 'weights' list.
            *TBD - replace w/ pysynphot observation object*
        wavelen : float, optional
            wavelength in meters for monochromatic calculation.
        save_intermediates : bool
            whether to output intermediate optical planes to disk. Default is False
        save_intermediate_what : string
            What to save - phase, intensity, amplitude, complex, parts, all. default is all.
        return_intermediates: bool
            return intermediate wavefronts as well as PSF?
        display : bool
            whether to display when finished or not.


        Returns
        -------
        outfits :
            a pyfits.HDUList
        """

        #if source is None and wavelength is None:
            #wavelength = 1e-6
        #if source is None and wavelength is not None:
            #source = {'weights': [1], 'wavelengths': [wavelength]}
        if not hasattr(wavelength,'__iter__'):
            wavelength = [wavelength]
        if weight is None:
            weight = [1.0] * len(wavelength)

        if len(tuple(wavelength)) != len(tuple(weight)):
            raise ValueError("Input source has different number of weights and wavelengths...")

        if display: p.clf()

        if save_intermediates or return_intermediates:
            self.intermediate_wfs = []
            _log.debug('reset intermediates in calcPSF')
            #raise ValueError("Saving intermediates for multi-wavelen not yet implemented!!")

        # loop over wavelengths
        if self.verbose: _log.info("Calculating PSF with %d wavelengths" % (len(wavelength)))
        outFITS = None

        normwts =  N.asarray(weight, dtype=float)
        normwts /= normwts.sum()

        if _USE_MULTIPROC and len(wavelength) > 1 : ######### Parallellized computation ############
            if _USE_FFTW3:
                _log.warn('IMPORTANT WARNING: Python multiprocessing and fftw3 do not appear to play well together. This is likely to crash intermittently')
                _log.warn('   We suggest you set   poppy._USE_FFTW3 = False   if you want to use calcPSFmultiproc().')

            if save_intermediates:
                raise NotImplementedError("Can't save intermediate steps if using parallelized code")
            self.intermediate_wfs = []

            # do *NOT* just blindly try to create as many processes as one has CPUs, or one per wavelength either
            # This is a memory-intensive task so that can end up swapping to disk and thrashing IO
            pool = multiprocessing.Pool(_MULTIPROC_NPROCESS)

            # build a single iterable containing the required function arguments
            _log.info("Beginning multiprocessor job")
            iterable = [(self, wavelen, wt, kwargs, _USE_FFTW3) for wavelen, wt in zip(wavelength, normwts)]
            results = pool.map(_wrap_propagate_for_multiprocessing, iterable)
            _log.info("Finished multiprocessor job")
            pool.close()

            # Sum all the results up into one array, using the weights
            outFITS = results[0]
            outFITS[0].data *= normwts[0]
            _log.info("got results for wavelength channel %d / %d" % (0, len(tuple(wavelength))) )
            for i in range(1, len(normwts)):
                _log.info("got results for wavelength channel %d / %d" % (i, len(tuple(wavelength))) )
                outFITS[0].data += results[i][0].data * normwts[i]
            outFITS[0].header.add_history("Multiwavelength PSF calc on %d processors completed." % nprocesses)





        else:  ########## single-threaded computations (may still use multi cores if FFTW3 enabled ######


            for wavelen, wave_weight in zip(wavelength, normwts):
            #for wavelen, weight in zip(source['wavelengths'], normwts):
                mono_psf = self.propagate_mono(wavelen, poly_weight=wave_weight, save_intermediates=save_intermediates or return_intermediates, **kwargs)
                # add mono_psf into the output array:

                if outFITS is None:
                    outFITS = mono_psf
                    outFITS[0].data = mono_psf[0].data*wave_weight
                else:
                    outFITS[0].data += mono_psf[0].data *wave_weight

            if save_intermediates:
                for i in range(len(self.intermediate_wfs)):
                    self.intermediate_wfs[i].writeto('wavefront_plane_%03d.fits' % i, what=save_intermediates_what )

        if display:
            cmap = matplotlib.cm.jet
            cmap.set_bad('0.3')
            #cmap.set_bad('k', 0.8)
            halffov =outFITS[0].header['PIXELSCL']*outFITS[0].data.shape[0]/2
            extent = [-halffov, halffov, -halffov, halffov]
            unit="arcsec"
            norm=LogNorm(vmin=1e-8,vmax=1e-1)
            p.xlabel(unit)

            imshow_with_mouseover(outFITS[0].data, extent=extent, norm=norm, cmap=cmap)


        # TODO update FITS header for oversampling here if detector is different from regular?
        waves = N.asarray(wavelength)
        wts = N.asarray(weight)
        mnwave = (waves*wts).sum() / wts.sum()
        outFITS[0].header.update('WAVELEN', mnwave, 'Weighted mean wavelength in meters')
        outFITS[0].header.update('NWAVES',waves.size, 'Number of wavelengths used in calculation')
        for i in range(waves.size):
            outFITS[0].header.update('WAVE'+str(i), waves[i], "Wavelength "+str(i))
            outFITS[0].header.update('WGHT'+str(i), wts[i], "Wavelength weight "+str(i))
        if _USE_FFTW3:
            ffttype = "pyFFTW3"
        else:
            ffttype = "numpy.fft"
        outFITS[0].header.update('FFTTYPE',ffttype, 'Algorithm for FFTs: numpy or fftw')

        if self.verbose: _log.info("PSF Calculation completed.")
        if return_intermediates:
            return outFITS, self.intermediate_wfs
        else:
            return outFITS


    def display(self, **kwargs):
        """ Display all elements in an optical system on screen.

        Any extra arguments are passed to the `optic.display()` methods of each element.

        """

        planes_to_display = [p for p in self.planes if not isinstance(p, Detector)]
        nplanes = len(planes_to_display)
        for i in range(nplanes):
            _log.info("Displaying plane %s in row %d" % (self.planes[i].name, i))
            self.planes[i].display(nrows=nplanes, row=i+1, **kwargs)



class SemiAnalyticCoronagraph(OpticalSystem):
    """ A subclass of OpticalSystem that implements a specialized propagation
    algorithm for coronagraphs whose occulting mask has limited and small support in
    the image plane. Algorithm from Soummer et al. (2007)

    The way to use this class is to build an OpticalSystem class the usual way, and then
    cast it to a SemiAnalyticCoronagraph, and then you can just call calcPSF on that in the
    usual fashion.

    Parameters
    -----------
    ExistingOpticalSystem : OpticalSystem
        An optical system which can be converted into a SemiAnalyticCoronagraph. This
        means it must have exactly 4 planes, in order Pupil, Image, Pupil, Detector.
    oversample : int
        Oversampling factor in intermediate image plane. Default is 8
    occulter_box : float
        half size of field of view region entirely including the occulter, in arcseconds. Default 1.0
        This can be a tuple or list to specify a rectangular region [deltaY,deltaX] if desired.


    Notes
    ------

    Note that this algorithm is only appropriate for certain types of Fourier transform,
    namely those using occulters limited to a sub-region of the image plane.
    It is certainly appropriate for TFI, and probably the right choice for NIRCam as well, but
    is of no use for MIRI's FQPMs.



    """

    def __init__(self, ExistingOpticalSystem, oversample=8, occulter_box = 1.0):

        if len(ExistingOpticalSystem.planes) != 4:
            raise ValueError("Input optical system must have exactly 4 planes to be convertible into a SemiAnalyticCoronagraph")
        self.name = "SemiAnalyticCoronagraph for "+ExistingOpticalSystem.name
        self.verbose = ExistingOpticalSystem.verbose
        self.source_offset_r = ExistingOpticalSystem.source_offset_r
        self.source_offset_theta = ExistingOpticalSystem.source_offset_theta
        self.planes = ExistingOpticalSystem.planes

        # SemiAnalyticCoronagraphs have some fixed planes, so give them reasonable names.
        self.inputpupil = self.planes[0]
        self.occulter = self.planes[1]
        self.lyotplane = self.planes[2]
        self.detector = self.planes[3]

        self.mask_function = InverseTransmission(self.occulter)

        for i, typecode in enumerate([PUPIL, IMAGE, PUPIL, DETECTOR]):
            if not self.planes[i].planetype == typecode:
                raise ValueError("Plane %i is not of the right type for a semianalytic coronagraph calculation: should be %s but is %s." % (i, typestrs[typecode], typestrs[self.planes[i].planetype]))


        self.oversample = oversample

        if hasattr(occulter_box, '__getitem__'):
            occulter_box = N.array(occulterbox) # cast to numpy array so the multiplication by 2 just below will work
        self.occulter_box = occulter_box

        self.occulter_det = Detector(self.detector.pixelscale/self.oversample, fov_arcsec = self.occulter_box*2, name='Oversampled Occulter Plane')

    def propagate_mono(self, wavelength=2e-6, normalize='first', save_intermediates=False, display_intermediates=False, intermediate_fn=None, poly_weight=None):
        """

        Parameters
        -------------
        wavelength : float
            Wavelength in meters
        normalize : string, {'first', 'last'}
            how to normalize the wavefront?
            * 'first' = set total flux = 1 after the first optic, presumably a pupil
            * 'last' = set total flux = 1 after the entire optical system.


        save_intermediates, display_intermediates, intermediate_fn, poly_weight : bools
            Ignored in current version of code?

        """
        if _TIMETESTS:
            t_start = time.time()
        if self.verbose: _log.info(" Propagating wavelength = %g meters" % wavelength)
        wavefront = self.inputWavefront(wavelength)

        if save_intermediates:
            raise NotImplemented("not yet")

        #------- differences from regular propagation begin here --------------
        wavefront *= self.inputpupil

        if normalize.lower()=='first':
            wavefront.normalize()


        if display_intermediates:
            nrows = 6
            P.clf()
            wavefront.display(what='best',nrows=nrows,row=1, colorbar=False)


        # determine FOV region bounding the image plane occulting stop.
        # determine number of pixels across that to use ("N_B")
        # calculate the MFT to the N_B x N_B occulting region.
        wavefront_cor = wavefront.copy()
        wavefront_cor.propagateTo(self.occulter_det)

        if display_intermediates: wavefront_cor.display(what='best',nrows=nrows,row=2, colorbar=False)


        # Multiply that by M(r) =  1 - the occulting plane mask function
        wavefront_cor *= self.mask_function

        if display_intermediates: wavefront_cor.display(what='best',nrows=nrows,row=3, colorbar=False)

        # calculate the MFT from that small region back to the full Lyot plane

        wavefront_lyot = wavefront_cor.copy()
        wavefront_lyot.propagateTo(self.lyotplane)

        if display_intermediates: wavefront_lyot.display(what='best',nrows=nrows,row=4, colorbar=False)
        # combine that with the original pupil function
        wavefront_combined = wavefront + (-1)*wavefront_lyot

        wavefront_combined *= self.lyotplane
        wavefront_combined.location = 'after combined Lyot pupil'

        if display_intermediates: wavefront_combined.display(what='best',nrows=nrows,row=5, colorbar=False)

        # propagate to the real detector in the final image plane.
        wavefront_combined.propagateTo(self.detector)

        if display_intermediates: wavefront_combined.display(what='best',nrows=nrows,row=6, colorbar=False)

        #------- differences from regular propagation end here --------------

        # prepare output arrays
        if normalize.lower()=='last':
                wavefront_combined.normalize()

        if _TIMETESTS:
            t_stop = time.time()
            _log.debug("\tTIME %f s\tfor propagating one wavelength" % (t_stop-t_start))

        return wavefront_combined.asFITS()


## TODO: Move all these test fns into the test suite .py file
def test_MFT():
    osys = OpticalSystem("Perfect JW", oversample=2)
    osys.addPupil(transmission="/Users/mperrin/software/newJWPSF/data/pupil.fits", name='JW Pupil')
    osys.addDetector(0.032, fov_pixels=128)


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
    osys.addDetector(0.032, name="Detector", fov_pixels=128)

    out = osys.propagate_mono(wavelength=10.65e-6, display_intermediates=True, save_intermediates=True)
    out.writeto('test_fft.fits',clobber=True)



def test_fftw3():
   #test_MFT()
    osys = OpticalSystem("Perfect JW", oversample=4)
    osys.addPupil(transmission="/Users/mperrin/software/newJWPSF/data/pupil.fits", name='JW Pupil')
    osys.addImage(function='fieldstop', name='20 arcsec stop', size=20)
    osys.addImage(function='FQPM',wavelength=10.65e-6)
    osys.addPupil(transmission="/Users/mperrin/software/newJWPSF/data/MIRI/coronagraph/MIRI_FQPMLyotStop.fits", name='MIRI FQPM Lyot')
    osys.addDetector(0.032, name="Detector", fov_pixels=128)



    nlam = 20
    nlam = 3
    source = {'wavelengths': N.linspace(10,15,nlam)*1e-6, 'weights': nlam* [1]}
    #source = {'wavelengths': N.arange(10,11,0.25)*1e-6, 'weights': 4* [1]}
    mono = {'wavelengths': [source['wavelengths'].mean()], 'weights': [1]}

    import time

    #print "-- mono, regular --"
    #t1 = time.time()
    #osys.calcPSF(mono).writeto('test_mono.fits', clobber=True)
    #print "-- poly, single process --"
    #t2 = time.time()
    #osys.calcPSF(source).writeto('test_single.fits', clobber=True)
    #print "-- poly, multi process --"
    #t3 = time.time()
    #osys.calcPSFmultiproc(source).writeto('test_multi.fits', clobber=True)
    #t4 = time.time()
    #print "for %d wavelengths: " % nlam
    #for t, v in zip(['mono:', 'poly single:', 'poly multi:'], [t2-t1, t3-t2, t4-t3]):
        #print "  Executed %s in %f seconds." % (t,v)

    #print "poly single relative computation time:\t%f" % ( (t3-t2)/(t2-t1)/nlam )
    #print "multi/single relative computation time:\t%f" % ( (t4-t3)/(t3-t2) )


    print "-- poly, singleprocess, numpy fft--"
    t2 = time.time()
    _USE_FFTW3 = False
    osys.calcPSF(source).writeto('test_numpyfft.fits', clobber=True)
    print "-- poly, single process, fftw3 --"
    t3 = time.time()
    _USE_FFTW3 = True
    osys.calcPSF(source).writeto('test_fftw3.fits', clobber=True)
    t4 = time.time()


    print "for %d wavelengths: " % nlam
    for t, v in zip(['Numpy FFT', 'FFTW3'], [t3-t2, t4-t3]):
        print "  Executed %s in %f seconds." % (t,v)



def test_speed():

    npix = 111
    offset = (0.0, 0)
    osys = OpticalSystem("Circle JW", oversample=1)
    #osys.addPupil(function='Square', name = "6.5m square", size=6.5)
    osys.addPupil(function='Circle', name = "6.5m circle", radius=6.5/2)
    #osys.addImage()
    #osys.addDetector(0.032, name="Detector", fov_pixels=128)
    osys.addDetector(0.032/2, name="Detector", fov_pixels=npix)

    osys.planes[-1].det_offset = offset

    #_USE_FFTW3 = True
    #_USE_FFTW3 = False
    #_TIMETESTS= True

    #res = osys.propagate_mono(2e-6,display_intermediates=False)

    src = {'wavelengths': [2e-6], 'weights': [1.0]}
    res = osys.calcPSF(src, display=True)


    #P.clf()
    #P.imshow(res[0].data)
    res.writeto('test_ci_np%d_off%s.fits' %(npix, offset), clobber=True)




if __name__ == "__main__":
    import pylab as P
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')



    if 0:
        for npix in (111,112):
        #for offset in [ (0,0)]:
            for offset in [ (0,0), (0,0.5), (1,0), (1,1.5)]:

                osys = OpticalSystem("Circle JW", oversample=1)
                #osys.addPupil(function='Square', name = "6.5m square", size=6.5)
                osys.addPupil(function='Circle', name = "6.5m circle", radius=6.5/2)
                #osys.addImage()
                #osys.addDetector(0.032, name="Detector", fov_pixels=128)
                osys.addDetector(0.032/2, name="Detector", fov_pixels=npix)

                osys.planes[-1].det_offset = offset

                #_USE_FFTW3 = True
                #_USE_FFTW3 = False
                #_TIMETESTS= True

                #res = osys.propagate_mono(2e-6,display_intermediates=False)

                src = {'wavelengths': [2e-6], 'weights': [1.0]}
                res = osys.calcPSF(src, display=True)


                #P.clf()
                #P.imshow(res[0].data)
                res.writeto('test_ci_np%d_off%s.fits' %(npix, offset), clobber=True)

    #test_speed()
