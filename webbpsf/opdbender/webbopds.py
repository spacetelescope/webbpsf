#!/usr/bin/env python
#import os, sys
from __future__ import division, print_function

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy
import astropy.table
import astropy.io.fits as fits
try:
    from IPython.core.debugger import Tracer; stop = Tracer()

except:
    pass
import logging
_log = logging.getLogger('webbopds')


from poppy import zernike

__doc__ = """


Linear Model for JWST based on influence functions provided by Erin Elliott

"""


################################################################################
class OPD(object):
    """ Base class for JWST Optical Path Difference (OPD) files

    Provides methods for analyzing and displaying OPDS. 

    If you want to manipulate one using the linear optical model, see the OPDbender class.


    Key class attributes include
    .OPD            the OPD data (as a fits.PrimaryHDU)
    .pupilmask      a mask for the pupil (as a fits.ImageHDU)
    .pixelscale     Scale in meters/pixel
    .header         reference to the OPD data's header


    """
    def __init__(self, opdfile=None, ext=0, slice=0, pupilfile=None, segment_mask_file='JWpupil_segments.fits', name=None):
        """
        Parameters
        ----------
        opdfile : str or fits.HDUList
            FITS file to load an OPD from. The OPD must be specified in microns. 
        ext : int, optional
            FITS extension to load OPD from
        slice : int, optional
            slice of a datacube to load OPD from, if the selected extension contains a datacube. 
        pupilfile : str
            FITS file for pupil mask, with throughput from 0-1. If not explicitly provided, will be inferred from 
            wherever is nonzero in the OPD file.
        """

        segnames = ['A'+str(i+1) for i in range(6)]         # names of each segment in order
        for i in range(6): segnames.append('B'+str(i+1))
        for i in range(6): segnames.append('C'+str(i+1))
        self.segnames = np.array(segnames)
        mypath = os.path.dirname(os.path.abspath( __file__ ))+os.sep
        self._segment_masks = fits.getdata(mypath + segment_mask_file)

        seg_rotangles= np.concatenate([np.arange(6)*60, np.arange(6)*60, np.arange(6)*60])  # rotation angle of each segment relative to the A1/B1/C1 set
        self._rotations =  {}
        for i in range(18): self._rotations[segnames[i]] = seg_rotangles[i]

        if pupilfile is None and opdfile is None:
            _log.warn('Neither a pupil mask nor OPD were specified. Using the default JWST pupil.')
            pupilfile = os.path.join(mypath,"JWpupil_1024.fits")
        if opdfile is not None:
            # load the OPD
            if isinstance(opdfile, fits.HDUList):
                self._opdHDU = opdfile[ext]
                if name is None: name='OPD from fits object'
            if isinstance(opdfile, fits.PrimaryHDU) or isinstance(opdfile, fits.ImageHDU):
                self._opdHDU = opdfile
                if name is None: name='OPD from fits object'
            elif isinstance(opdfile, str):
                self._opdHDU = fits.open(opdfile)[ext]
                if name is None: name='OPD from '+opdfile
            else:
                raise TypeError("not sure what to do with OPD of type "+str(type(opdfile)))


            _log.info("Loaded OPD from %s, extension %s " % (opdfile, ext))
            self._opdHDU.header.add_history("Loaded OPD from %s, extension %s " % (opdfile, ext))
            try:
                self.pixelscale = self._opdHDU.header['PUPLSCAL']
            except:
                raise RuntimeError("No PIXELSCL value from header")
                # self.pixelscale = 6.559/ 1016. # FIXME where was this from?

            if self._opdHDU.data.ndim > 2:
                self._opdHDU.data = self._opdHDU.data[slice, :,:]
                _log.info("OPD file contains a datacube, selecting OPD # %d" % slice)

            if (self._segment_masks.shape != self._opdHDU.data.shape):
                _log.warn('Segment mask file is different size from OPD; per-segment math will not work properly.')
            #assert(self._segment_masks.shape == self._opdHDU.data.shape)


            if pupilfile is not None:
                # load the pupil
                if isinstance(pupilfile, fits.HDUList):
                    self._pupilHDU = pupilfile[ext]
                    if name is None: name='pupil from fits object'
                if isinstance(pupilfile, fits.PrimaryHDU) or isinstance(pupilfile, fits.ImageHDU):
                    self._pupilHDU = pupilfile
                    if name is None: name='pupil from fits object'
                elif isinstance(pupilfile, str):
                    self._pupilHDU = fits.open(pupilfile)[ext]
                    if name is None: name='pupil from '+pupilfile
                else:
                    raise TypeError("not sure what to do with pupil of type "+str(type(pupilfile)))


            else:
                self._pupilHDU = fits.ImageHDU(np.ones_like(self._opdHDU.data).astype(np.uint8))
                wz = np.where(self._opdHDU.data == 0)
                self._pupilHDU.data[wz] = 0

        else:
            # load the pupil, make the OPD all zeros.
            #pupil = fits.open(pupilfile)[ext]
            self._pupilHDU =  fits.open(pupilfile)[ext]
            self._opdHDU = self._pupilHDU.copy() #copy.deepcopy(pupil)
            self._opdHDU.data[:] = 0
            self._opdHDU.header.add_history('Pupil mask: '+str(pupilfile))
            self._opdHDU.header.add_history('No OPD supplied; created a zero OPD based on pupil mask')
            if name is None: name="Null OPD from pupil in "+pupilfile

        self.name=name
        self.header = self._opdHDU.header # convenience reference

    @property
    def data(self):
        return self._opdHDU.data

    @property
    def pupil(self):
        return self._pupilHDU.data


    def copy(self):
        """ Make a copy of a wavefront object """
        from copy import deepcopy
        return deepcopy(self)

    # We can add and subtract OPD objects using the following.
    def __add__(self, x):
        """ Operator add """

        output = self.copy()
        if isinstance(x, OPD):
            output._opdHDU.data += x.data
            output.name += " + "+x.name
        else:            
            output._opdHDU.data += x
            output.name += " + "+str(x)
        return output

    def __radd__(self,x):
        return self.__add__(x)

    def __sub__(self, x):
        """ Operator subtract """

        output = self.copy()
        if isinstance(x, OPD):
            output._opdHDU.data -= x.data
            output.name += " - "+x.name
        else:            
            output._opdHDU.data -= x
            output.name += " - "+str(x)
        return output


    def asFITS(self, include_pupil=True):
        """ Return an OPD as a fits.HDUList object

        Parameters
        -----------
        include_pupil : bool
            Include the pupil mask as a FITS extension?
        """

        output = fits.HDUList([self._opdHDU])
        output[0].header.update('EXTNAME','OPD')
        if include_pupil: 
            #puphdu= fits.ImageHDU(self._pupilHDU, name='PUPIL')
            self._pupilHDU.header.update('EXTNAME','PUPIL')
            output.append(self._pupilHDU)
            #hdus.append(puphdu)
            #hdus[0].header.update("EXTEND",'True', 'File may contain extensions')
        return output

    def writeto(self, outname, clobber=True, **kwargs):
        """ Write OPD to a FITS file on disk """
        self.asFITS(**kwargs).writeto(outname, clobber=clobber)

    #---- display and analysis
    def powerspectrum(self, max_cycles=50, sampling=5, vmax=100, iterate=False):
        """ 
        Compute the spatial power spectrum of an aberrated wavefront. 

        Produces nice plots on screen. 

        Returns an array [low, mid, high] giving the RMS spatial frequencies in the
        different JWST-defined spatial frequency bins:
            low:   <=5 cycles/aperture
            mid:   5 < cycles <= 30
            high:  30 < cycles 

        """
        import SFT

        def tvcircle(radius=1, xcen=0, ycen=0, center=None,**kwargs):
            """
                draw a circle on an image.

                    radius
                    xcen
                    ycen
                    center=     tuple in (Y,X) order.
            """
            if center is not None:
                xcen=center[1]
                ycen=center[0]
            t = np.arange(0, np.pi * 2.0, 0.01)
            t = t.reshape((len(t), 1))
            x = radius * np.cos(t) + xcen
            y = radius * np.sin(t) + ycen
            plt.plot(x,y, **kwargs)

        cmap = matplotlib.cm.jet
        cmap.set_bad('0.3')

        plt.clf()
        plt.subplot(231)
        self.draw(title='full wavefront', clear=False, colorbar=False, vmax=vmax)


        ps_pixel_size = 1./sampling # how many cycles per pixel
        trans = SFT.SFT3(self.data, max_cycles*2, max_cycles*2*sampling)

        #if iterate:
            #inverse1 = SFT.SFT3(trans, max_cycles*2, self._opdHDU.data.shape[0], inverse=True)
            #plt.subplot(232)
            #plt.imshow(inverse1.real , vmin=(-vmax)/1000., vmax=vmax/1000, cmap=cmap)  # vmax is in nm, but WFE is in microns, so convert
            #stop()


        abstrans = np.abs(trans)

        extent = [-max_cycles, max_cycles, -max_cycles, max_cycles]


        plt.subplot(233)
        plt.imshow(abstrans, extent=extent)
        plt.title("Power Spectrum of the phase")
        plt.ylabel("cycles/aperture")
        tvcircle(radius=5, color='k', linewidth=1) #, ls='--')
        tvcircle(radius=30, color='k', linewidth=1) #2, ls='--')
        plt.gca().set_xbound(-max_cycles, max_cycles)
        plt.gca().set_ybound(-max_cycles, max_cycles)


        y, x = np.indices(abstrans.shape)
        y -= abstrans.shape[0]/2.
        x -= abstrans.shape[1]/2.
        r = np.sqrt(x**2+y**2) * ps_pixel_size

        mask = np.ones_like(self.data)
        mask[np.where(self._pupilHDU.data == 0)] = np.nan
        wgood = np.where(self._pupilHDU.data != 0)


        components = []
        for i, label in enumerate(['low', 'mid', 'high']):
            plt.subplot(2,3,i+4)
            if label =='low':
                condition = r <= 5
            elif label =='mid':
                condition = (r > 5) &  (r <=30)
            else:
                condition = (r > 30)
            filtered = trans * condition

            inverse = SFT.SFT3(filtered, max_cycles*2, self._opdHDU.data.shape[0], inverse=True)
            inverse = inverse[::-1, ::-1] # I thought SFT did this but apparently this is necessary to get the high freqs right...

            plt.imshow(inverse.real * mask, vmin=(-vmax)/1000., vmax=vmax/1000, cmap=cmap)  # vmax is in nm, but WFE is in microns, so convert
            plt.title(label +" spatial frequencies")
            rms = (np.sqrt( (inverse.real[wgood]**2).mean())*1000)

            components.append(rms)
            plt.xlabel("%.3f nm RMS WFE" %  rms )

        return np.asarray(components)

    def draw(self, ax=None, labelsegs=True, vmax=150., colorbar=True, clear=True, title= None, unit='nm', cbpad=None, colorbar_orientation='vertical'):
        """ Draw on screen the perturbed OPD
        
        Parameters
        -----------
        ax : matplotlib.Axes
            axes instance to display into. 
        labelsegs : bool
            draw segment name labels on each segment? default True. 

        clear : bool
            Clear plot window before display? default true
        unit : str
            Unit for WFE. default is 'nm'

        """

        _log.info("in draw")
        if unit == 'nm':
            scalefact = 1000.
        elif unit =='micron':
            scalefact = 1.0
        else: 
            raise ValueError("unknown unit keyword")

        if clear: 
            plt.clf()
        cmap = matplotlib.cm.jet
        cmap.set_bad('0.3')

        mask = np.ones_like(self._opdHDU.data)
        mask[np.where(self._pupilHDU.data == 0)] = np.nan

        try: 
            pupilscale = self._opdHDU.header['PUPLSCAL']
            s = self._opdHDU.data.shape
            extent = [a*pupilscale for a in [-s[0]/2, s[0]/2, -s[1]/2,s[1]/2]]
        except:
            extent=None

        if ax is None:
            ax = plt.gca()


        plot = ax.imshow(self._opdHDU.data*mask*scalefact, vmin=-vmax, vmax=vmax, cmap=cmap) #, extent=extent)

        _log.debug("drawing OPD. Vmax is %f, data max is %f " % (vmax, self._opdHDU.data.max()))

        #if self.remove_piston_tip_tilt: title +", Piston/tip/tilt removed"
        if title is None:
            title=self.name
        ax.set_title(title)
        ax.set_xlabel("RMS WFE = %.1f nm" % self.rms())

        if labelsegs:
            Y, X = np.indices(self._opdHDU.data.shape)
            for iseg in range(18):
                wseg = np.where(self._segment_masks == iseg+1)
                cx = np.mean([X[wseg].min(), X[wseg].max()])
                cy = np.mean([Y[wseg].min(), Y[wseg].max()])
                plt.text(cx,cy, self.segnames[iseg], color='k', horizontalalignment='center',verticalalignment='center')
        if colorbar:
            if cbpad is None:
                cbpad = 0.05 if colorbar_orientation=='vertical' else 0.15
            #pts = plt.gca().get_position().get_points()
            #col_ax = plt.gcf().add_axes([0.92, pts[0,1]+0.08, 0.02, 0.62 ])
            cb= plt.colorbar(plot, ax=ax, pad=cbpad, orientation=colorbar_orientation)
            cb.set_label("WFE [%s]" % unit)
        plt.draw()

    def label_seg(self, segment, ax=None):
        Y, X = np.indices(self._opdHDU.data.shape)

        base = {'A':0, 'B':6,'C':12}
        try:
            iseg = base[segment.upper()[0]]+int(segment[1])
        except:
            return
        wseg = np.where(self._segment_masks == iseg)
        pupilscale = self._opdHDU.header['PUPLSCAL']
        cx = (np.mean([X[wseg].min(), X[wseg].max()])  -512) * pupilscale
        cy = (np.mean([Y[wseg].min(), Y[wseg].max()])  -512) * pupilscale

        if ax is None: ax = plt.gca()
        print(segment, cx, cy)
        label = ax.text(cx,cy, segment, color='k', horizontalalignment='center',verticalalignment='center')
        ax.get_figure().canvas.draw()
        return label

    def zern_seg(self, segment, vmax=150, unit='nm'):
        """ Show the Zernike terms applied to a given segment"""

        # the actual wavefront is always in units of microns.
        if unit == 'nm':
            scalefact = 1000.
        elif unit =='micron':
            scalefact = 1.0
        else: 
            raise ValueError("unknown unit keyword")




        nzerns = 11
        title = "Zernike components for "+segment
        zerns = np.zeros((nzerns))
        if segment+"-decenter" in self.state.keys():
            zerns += self._move(segment, type='decenter', vector=self.state[segment+"-decenter"], return_zernikes=True)
            title += ", displacement=[%.2e, %.2e, %.2e] um" % tuple(self.state[segment+"-decenter"])
        if segment+"-tilt" in self.state.keys():
            zerns += self._move(segment, type='tilt', vector=self.state[segment+"-tilt"], return_zernikes=True)
            title += ", tilts =[%.5f, %.5f, %.5f] urad" % tuple(self.state[segment+"-tilt"])

        _log.info("Zerns: "+str(zerns))
        fig = plt.gcf()
        fig.clf()

        npix=200
        hexap = zernike.hex_aperture(npix)
        hexap[np.where(hexap == 0)] = np.nan
        cmap = matplotlib.cm.jet
        cmap.set_bad('0.5', alpha=0.0)

        for j in np.arange(nzerns)+1:
            ax = fig.add_subplot(3, 4, j, frameon=False, xticks=[], yticks=[])
            # n, m = zernike.noll_indices(j)
            Z = zernike.zernike1(j, npix=npix)
            ax.imshow(Z * zerns[j-1] * hexap*scalefact, vmin=-1*vmax, vmax=vmax, cmap=cmap)
            ax.text(npix*0.95, npix*0.8, "$Z%d$" % (j), fontsize=20, horizontalalignment='right')
            ax.text(npix*0.95, npix*0.1, "%.2e" % (zerns[j-1]), fontsize=15, horizontalalignment='right')

        fig.text(0.5, 0.95, title, horizontalalignment='center', fontsize=15)
        fig.text(0.95, 0.15, 'Segment RMS WFE = %.2f %s ' % (self.rms(segment)*(scalefact/1000), unit), fontsize=10, horizontalalignment='right')
        #fig.text(0.95, 0.05, 'OPDs scaled from %.2e - %.2e um' % (-vmax, vmax), horizontalalignment='right', fontsize=10)
        fig.text(0.95, 0.05, 'OPDs scaled from %.2f - %.2f %s' % (-vmax, vmax, unit), horizontalalignment='right', fontsize=10)

        plt.draw()

    def rms(self, segment=None):
        """ Return RMS WFE in nanometers 
        
        Parameters
        ----------
        segment : string
            Segment name, to compute RMS for a single segment. Leave unspecified (None) to compute
            RMS WFE for the entire aperture. Segments are identified by name: A1-A6, B1-B6, C1-C6
        """
        if segment is None:
            # RMS for whole aperture
            wgood = np.where(self._pupilHDU.data != 0)
        else:
            assert(segment in self.segnames)
            iseg = np.where(self.segnames == segment)[0][0]+1  # segment index from 1 - 18
            wseg = np.where(self._segment_masks == iseg)
            wgood = wseg
        #return self._opdHDU.data[wgood].std() * 1000
        return np.sqrt(np.nanmean((self._opdHDU.data[wgood]**2))) * 1000

    def ptv(self, segment=None):
        """return peak to valley WFE in nanometers

        Parameters
        ----------
        segment : string
            Segment name, to compute RMS for a single segment. Leave unspecified (None) to compute
            for the entire aperture.
        """
        if segment is None:
            # RMS for whole aperture
            wgood = np.where(self._pupilHDU.data != 0)
        else:
            assert(segment in self.segnames)
            iseg = np.where(self.segnames == segment)[0][0]+1  # segment index from 1 - 18
            wseg = np.where(self._segment_masks == iseg)
            wgood = wseg
 
        peak = self._opdHDU.data[wgood].max()
        valley =  self._opdHDU.data[wgood].min()
        return (peak-valley) * 1000



    def estimated_Strehl(self, wavelength, verbose=True):
        """ Compute an estimated Strehl given a wavelength in meters
        
        Parameters
        -----------
        wavelength : float
            in meters
        verbose : bool
            should I print out an informative message?
            """
        # rms = sqrt((-1.0)*alog(strehl))*wavelength/(2*!pi)

        rms = self.rms() * 1e-9
        strehl  = np.exp(-(rms * 2 * np.pi / wavelength)**2)
        if verbose: 
            print("For wavelength = {0} meters, estimated Strehl = {1} for {2} nm rms WFE".format(wavelength, strehl, rms*1e9))
        return strehl


################################################################################

class OPDbender(OPD):
    """ Perturb an existing wavefront OPD file, by applying changes in WFE
    derived from Erin Elliott's linear optical model.
    
    Note that this model is strictly accurate only for a single NIRCam field 
    point, but we can apply it elsewhere and it should not be too far off for small motions.

    TBD to quantify the range of applicability and uncertainties.
    
    """

    def __init__(self, opdfile=None, ext=0, slice=0, pupilfile=None, segment_mask_file='JWpupil_segments.fits', rm_piston_tilt=False, zero=False):
        """
        Parameters
        ----------
        opdfile : str or fits.HDUList
            FITS file to load an OPD from. The OPD must be specified in microns. 
        ext : int, optional
            FITS extension to load OPD from
        slice : int, optional
            slice of a datacube to load OPD from, if the selected extension contains a datacube. 
        pupilfile : str
            FITS file for pupil mask, with throughput from 0-1. If not explicitly provided, will be inferred from 
            wherever is nonzero in the OPD file.


        """

        OPD.__init__(self, opdfile=opdfile, ext=ext, slice=slice, pupilfile=pupilfile, segment_mask_file=segment_mask_file)

        mypath = os.path.dirname(os.path.abspath( __file__ ))+os.sep
        self._sensitivities = astropy.table.Table.read(os.path.join(mypath,  'seg_sens.txt'), format='ascii', delimiter='\t')
        self.state = {}
        self.remove_piston_tip_tilt = rm_piston_tilt

#
#
#        if pupilfile is None and opdfile is None:
#            _log.warn('Neither a pupil mask nor OPD were specified. Using the default JWST pupil.')
#
#            import webbpsf
#            pupilfile = webbpsf.JWInstrument()._WebbPSF_basepath+os.sep+"pupil_RevV.fits"
#        if opdfile is not None:
#            # load the OPD
#            if isinstance(opdfile, fits.HDUList):
#                self._opdHDU = opdfile[ext]
#            else:
#                self._opdHDU = fits.open(opdfile)[ext]
#            _log.info("Loaded OPD from %s, extension %s " % (opdfile, ext))
#            self._opdHDU.header.add_history("Loaded OPD from %s, extension %s " % (opdfile, ext))
#            try:
#                self.pixelscale = self._opdHDU.header['PIXELSCL']
#            except:
#                self.pixelscale = 6.559/ 1016.
#
#            if self._opdHDU.data.ndim > 2:
#                self._opdHDU.data = self._opdHDU.data[slice, :,:]
#                _log.info("OPD file contains a datacube, selecting OPD # %d" % slice)
#            assert(self._segment_masks.shape == self._opdHDU.data.shape)
#
#
#            if pupilfile is None:
#                self._pupilHDU = np.ones_like(self._opdHDU.data).astype(np.uint8)
#                wz = np.where(self._opdHDU.data == 0)
#                self._pupilHDU[wz] = 0
#        else:
#            # load the pupil, make the OPD all zeros.
#            pupil = fits.open(pupilfile)[ext]
#            self._pupilHDU = pupil.data
#            self._opdHDU = pupil.copy() #copy.deepcopy(pupil)
#            self._opdHDU.data[:] = 0
#            self._opdHDU.header.add_history('Pupil mask: '+pupilfile)
#            self._opdHDU.header.add_history('No OPD supplied; created a zero OPD based on pupil mask')
#
        self._opdHDU_orig = self._opdHDU.copy()
        if zero: self.zero()


    #---- overall state manipulation

    def reset(self):
        """ Reset an OPD to the state it was loaded from disk.  

        i.e. undo all segment moves. 
        """
        self._opdHDU.data = self._opdHDU_orig.data.copy()
        _log.info("Reset to unperturbed OPD")

    def zero(self):
        self._opdHDU.data *= 0
        self._opdHDU_orig.data *= 0
        self.state = {}
        _log.info("Set OPD to zero WFE!")

    def load_state(self, newstate):
        _log.info("Loading new state.")
        self.zero()
        for k in newstate.keys():
            segname, motion =  k.split('-')
            self._move(segname, motion, newstate[k])

    def print_state_v1(self):
        keys = self.state.keys()
        keys.sort()
        if keys is None:
            print("No perturbations")
        else:
            print("state = {")
            for k in keys:
                print("    '%s' : %s," % (k, repr(self.state[k])))
            print("    }")

    def print_state(self, type='Report'):
        keys = self.state.keys()

        if type == 'Report':
            print("Segment positions in Report coordinates: (microns for decenter, microradians for alpha & beta, milliradians for gamma):")
            print("  \t %10s %10s %10s %10s %10s %10s" % ("X dec", "Y dec", "Z dec", "alpha", "beta", "gamma"))
            for segment in self.segnames:
                if segment+"-tilt" in keys:
                    tilts = self.state[segment+"-tilt"].tolist()
                else: 
                    tilts = [0,0,0]
                if segment+"-decenter" in keys:
                    decenters = self.state[segment+"-decenter"].tolist()
                else: 
                    decenters = [0,0,0]

                print("%2s\t %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" % tuple([segment]+decenters+tilts))
        elif type == 'Matlab':
            print("Segment positions in Matlab coordinates: (millimeters and degrees)")
            print("  \t %12s %12s %12s %12s %12s %12s" % ("X dec", "Y dec", "Z dec", "alpha", "beta", "gamma"))
            for segment in self.segnames:
                if segment+"-tilt" in keys:
                    tilts = self.state[segment+"-tilt"].tolist()
                else: 
                    tilts = [0,0,0]
                if segment+"-decenter" in keys:
                    decenters = self.state[segment+"-decenter"].tolist()
                else: 
                    decenters = [0,0,0]

                decenters = (np.array([decenters[1], decenters[0], -1 * decenters[2]])/1000).tolist()
                tilts[2] *= 1000 # convert millirad to microrad for consistency 
                tilts = (np.array([tilts[1], tilts[0], -1 * tilts[2]])*1e-6*180/np.pi).tolist()

                print("%2s\t %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e" % tuple([segment]+decenters+tilts))

    #---- segment manipulation via linear model

    def _get_seg_sensitivities(self, segment='A1', type='decenter'):
        assert(segment in self.segnames)
        refseg = 2 if segment[0] == 'C' else 1
        fields = ['%s_%s_%s%d' % (axis, type, segment[0], refseg) for axis in ['X', 'Y', 'Z']]

        return np.asarray([self._sensitivities[f] for f in fields])

    def _record(self, segment='A1', type='decenter', values=np.zeros((3))):
        """ update the state structure with the current segment positions """
        key = "%s-%s" % (segment, type)
        if key in self.state.keys():
            self.state[key] += values
        else:
            self.state[key] = values

    def print_(self):
        keys = self.state.keys()
        keys.sort()
        if keys is None:
            print("No perturbations")
        else:
            for k in keys:
                print("%s\t%s" % (k, str(self.state[k])))

    def _apply_zernikes_to_seg(self, segment, zernike_coeffs):
        """ Apply Zernike perturbations to a given segment """
        assert(segment in self.segnames)

        iseg = np.where(self.segnames == segment)[0][0]+1  # segment index from 1 - 18
        wseg = np.where(self._segment_masks == iseg)

        # determine the X and Y zernike coefficients for each segment
        # determine the center of each segment, as the mean of the min and max X and Y values
        Y, X = np.indices(self._opdHDU.data.shape)
        cx = np.mean([X[wseg].min(), X[wseg].max()])
        cy = np.mean([Y[wseg].min(), Y[wseg].max()])

        seg_radius = (X[wseg].max()-X[wseg].min())/2.0

        _log.debug("Segment %s is centered at pixel loc (%.1f, %.1f) with radius %.1f pix" % (segment, cx, cy, seg_radius))

        Xw = X[wseg].astype(np.float64)
        Yw = Y[wseg].astype(np.float64)
        Yw -= cy
        Xw -= cx
        ang = self._rotations[segment]*np.pi/180  # This definitely has to be a positive sign, 
                                                 # to get the right X, Y for locally-defined zernikes
        Xr = Xw * np.cos(ang) + Yw * np.sin(ang)
        Yr = Xw *-np.sin(ang) + Yw * np.cos(ang)

        theta = np.arctan2(Yr, Xr)
        Rw = np.sqrt(Xr**2 + Yr**2)/seg_radius

        #zernike_coeffs *= 0
        #zernike_coeffs[2]=1
        if self.remove_piston_tip_tilt:
            zernike_coeffs[0:3] = 0
        for i in range(len(zernike_coeffs)):
            zern = zernike.zernike1(i+1, rho=Rw, theta=theta, outside=0.0) * zernike_coeffs[i]
            #print "Z%d = %f" % (i+1, zernike_coeffs[i])
            self._opdHDU.data[wseg] += zern
            #self._opdHDU.data[wseg] = Yw

        outtxt="Zs=["+", ".join(['%.1e'%z for z in zernike_coeffs])+"]"
        #print outtxt
        _log.debug("     "+outtxt)

    def _move(self, segment, type='tilt', vector=None, draw=False, return_zernikes=False):
        """ 
        Internal function that actually modifies the OPD by applying the linear optical model.
        Don't use this directly, use tilt() or displace() instead. 

        """

        #print "Segment rotation for %s is %f" % (segment, self._rotations[segment])
        ang = self._rotations[segment] * np.pi/180
        local_coordX = vector[0] * np.cos(ang) + vector[1] * np.sin(ang)
        local_coordY = vector[0] *-np.sin(ang) + vector[1] * np.cos(ang)
        local_vector = np.array([local_coordX, local_coordY, vector[2]])
        if type=='tilt':
            local_vector[2] /= 1000 # convert Z tilt to milliradians instead of microradians because that is what the sensitivity tables use
            units = 'microradians for tip/tilt, milliradians for clocking'
        else:
            units = 'microns'

        sensitivities = self._get_seg_sensitivities(segment, type=type)

        zernike_coeffs_2d = sensitivities * local_vector[:, np.newaxis]
        zernike_coeffs = zernike_coeffs_2d.sum(axis=0)

        if return_zernikes: 
            return zernike_coeffs
        else:
            _log.info("Segment %s requested %s: %s in %s" % (segment, type, str(vector), units))
            _log.debug("    local %s: %s"  % (type, str(local_vector)))
            self._record(segment, type, vector)
            self._apply_zernikes_to_seg(segment, zernike_coeffs)

            if draw: self.draw()

    def tilt(self, segment, tiltX=0.0, tiltY=0.0, tiltZ=0.0, unit='urad', draw=False):
        """ Tilt/rotate a segment some angle around X, Y, or Z.

        Parameters
        -----------
        segment : str
            Segment name, e.g. 'A1'
        tiltX, tiltY, tiltZ : floats
            Tilt angle, in microradians by default.
        unit : str
            Unit for displacements. Can be 'urad', 'radian', 'arcsec', 'arcmin', 'milliarcsec'
        draw: bool
            Display after moving? 

        """
        tilts = np.array([tiltX, tiltY, tiltZ]).astype(float)

        self._opdHDU.header.add_history('Rotation: %s %s' % (str(tuple(tilts)), unit))


        # new sensitivity matrices are in urad for  alpha and beta, mrad for gamma.
        # first just convert all to urad.
        if unit.endswith('s'): unit = unit[:-1]
        unit = unit.lower()
        #print "UNIT IS "+unit
        if unit == 'urad': pass
        elif unit =='milliarcsec': tilts *= (1e6*np.pi/ (180.*60*60*1000))
        elif unit =='arcsec': tilts *= (1e6*np.pi/ (180.*60*60))
        elif unit =='arcmin': tilts *= (1e6*np.pi/ (180.*60))
        elif unit == 'radian' or unit=='rad': tilts*= 1e6
        else: raise NotImplemented('unknown unit')

        #tilts[2]  /= 1000 # convert Z tilt to milliradians instead because that is what the sensitivity tables use
        # now this happens in the _move routine.

        self._move(segment, 'tilt', tilts, draw=draw)

    def displace(self, segment, distX=0.0, distY=0.0, distZ=0.0, unit='micron', draw=False):
        """ Move a segment some distance in X, Y, and Z.

        Parameters
        -----------
        segment : str
            Segment name, e.g. 'A1'
        distX, distY, distZ : floats
            Displacement distance, in microns by default.
        unit : str
            Unit for displacements. Can be 'micron', 'millimeter','nanometer', 'mm', 'nm', 'um'
        draw : bool
            Display after moving? 

        """
        vector = np.array([distX, distY, distZ])

        self._opdHDU.header.add_history('Displacement: %s %s' % (str(tuple(vector)), unit))
        if unit.endswith('s'): unit = unit[:-1]
        unit = unit.lower()
        if unit == 'micron' or unit =='um': pass
        elif unit =='millimeters' or unit =='millimeter' or unit =='mm': vector *= 1000
        elif unit =='nm' or unit =='nanometer' or unit =='nanometers' : vector /= 1000
        else: raise ValueError("Unknown unit for length: %s" % unit) 


        self._move(segment, 'decenter', vector, draw=draw)

    def sinewave(self, cyclesX = 0.0, cyclesY = 0.0, amplitude = 100.,  draw=True, stop=False):
        """ Add a (nonphysical) sine-wave phase of a given spatial frequency in X and Y, specified
            in terms of cycles per aperture.


            Those cycles for spatial frequency are defined in accordance with the JWST OTE spec,
            i.e. cycles per 6.605 m circumscribing circle, rather than the actual aperture shape.

            cyclesX, Y : float
                cycles / aperture
            amplitude : float
                amplitude in nm

        """

        if cyclesX==0 and cyclesY==0:
            _log.info("Must specify either X or Y cycles - doing nothing since neither was given")
            return

        ref_ap_diam = 6.605
        Y, X = np.indices(self._opdHDU.data.shape, dtype=float)

        Y -= Y.mean() # center
        X -= X.mean() # center
        #Y *= 2*np.pi* cyclesY/(ref_ap_diam/self.pixelscale)
        #X *= 2*np.pi* cyclesX/(ref_ap_diam/self.pixelscale)

        wsegs = np.where(self._segment_masks != 0 )

        #self._opdHDU.data[wsegs] += (np.cos(X)[wsegs] + np.cos(Y)[wsegs]) * (amplitude/1000.)  # convert amplitude to microns

        self._opdHDU.data[wsegs] += np.cos(   2*np.pi*( X*cyclesX + Y*cyclesY ) /(ref_ap_diam/self.pixelscale) )[wsegs] * (amplitude/1000.)  # convert amplitude to microns

        #self._opdHDU.data[wsegs] = 0.2 # debug

        _log.info("added sine wave: (%.2f, %.2f, %.2f)" % (cyclesX, cyclesY, amplitude))
        #stop()

        if draw: self.draw()
        if stop: Tracer()()

    def perturb_all(self, draw=True, verbose=True, **kwargs):
        """ Randomly perturb all segments


        There is no good easy way to dial in a desired level of WFE right now. 
        Best/easiest is to try iteratively adjusting the multiplier keyword and 
        repeat until you get the desired level of RMS WFE.
 

        """
        for seg in self.segnames:
            self.perturb(seg, **kwargs)

        if draw:
            plt.clf()
            self.draw(colorbar=True, vmax=1)
        if verbose:
            print("")
            self.print_state(type='Report')
            print("")
            self.print_state(type='Matlab')

    def perturb(self, segment, max=False, multiplier=1.0):
        """ Randomly perturb a segment 


        There is no good easy way to dial in a desired level of WFE right now. 
        Best/easiest is to try iteratively adjusting the multiplier keyword and 
        repeat until you get the desired level of RMS WFE.
        
        Parameters
        ----------
        segment : str
            segment name
        multiplier : float
            Scale factor for making larger or smaller motions.
        max : bool
            force the motion to the largest possible motion consistent with the linear optical
            model bounds. 

        """
        bound = [13, 13, 0.2, 0.9, 0.9, 200] # bounds on linear regime, the region over which the
                                             # linear model can be considered trustworthy.
        stepsize = [1, 1, 0.1, 0.1, 0.1, 1]  # what relative size motion should we put in the various terms?

        
        if max:
            steps = bound

        else: 
            # stick in a random step between -1 and 1 times the desired step size
            steps = (np.random.random_sample(6)-0.5)*2 * np.array(stepsize) * multiplier

        self.displace( segment, steps[0], steps[1], steps[2])
        self.tilt( segment, steps[3], steps[4], steps[5])

        rms = self.rms(segment)
        print("After perturbation, segment %s has RMS WFE = %.1f nm"  % (segment, rms))

#--------------------------------------------------------------------------------
def segment_primary(infile='JWpupil.fits'):
    """ Given a FITS file containing the JWST aperture, create an array
    with the PMSAs numbered from 1 - 18. 

    This is used to create the 'JWpupil_segments.fits' file used by the OPDbender code. 
    You should not need to run this routine unless you want to change the model of the
    JWST pupil.

    The basic algorithm is general, but right now this routine contains
    a bunch of hard-coded values tuned for the 1024x1024 Rev V pupil, 
    so probably you will need to edit it for any other pupil. 

    M. P. 2011-02-15


    """

    def renumber_array(array):
        uniq, inds = np.unique1d(array, return_inverse=True)
        #res = np.zeros_like(res2.shape)
        newuniq = range(len(uniq))
        res = np.array([newuniq[i] for i in inds]).astype(np.int16)
        res.shape = array.shape
        #for i in range(len(uniq)):
            #print "Segment %s has %d pixels" % (i, np.where( res == i)[0].size)
        return res
     
    im = fits.getdata(infile)
    markers = np.zeros_like(im).astype(np.int16)
    xm, ym = np.ogrid[0:1024:102, 0:1024:100]
    markers[xm, ym]= np.arange(xm.size*ym.size).reshape((xm.size,ym.size))
    res2 = scipy.ndimage.watershed_ift(im.astype(np.uint8), markers)
    res2[xm, ym] = res2[xm-1, ym-1] # remove the isolate seeds


    res = renumber_array(res2)

    # split and merge segments as needed - hard coded for Rev V pupil

    res3 = res.copy() #np.zeros_like(res2.shape)
    merge_inds = [[0,3], [2,4], [11,15], [18,19]]
    for pair in merge_inds:
        res3[ np.where( res == pair[1]) ] = pair[0]

    Y, X = np.indices(res.shape)
    w10h = np.where( (res == 10) & ( Y >= 512))
    res3[w10h] = 31
    w10h = np.where( (res == 12) & ( Y >= 512))
    res3[w10h] = 32
    
    res3 = renumber_array(res3)



    plt.clf()
    #plt.subplot(121)
    #plt.imshow(res)
    plt.subplot(121)
    plt.imshow(res3)
 

    for i in range(19):
        w = np.where(res3 == i)
        mx = X[w].mean()
        my = Y[w].mean()
        plt.text(mx, my, str(i), color='k')

    segs = ['A'+str(i+1) for i in range(6)]
    for i in range(6): segs.append('B'+str(i+1))
    for i in range(6): segs.append('C'+str(i+1))
    seg_inds = [9, 18, 10, 4, 8, 17, 15, 14, 5, 1, 3, 11, 13, 7, 2, 0, 6, 12]


    result = np.zeros_like(res3)
    mxs = []
    mys = []
    for i in range(18):
        w = np.where(res3 == seg_inds[i])
        result[w] = i+1
        mxs.append(X[w].mean())
        mys.append(Y[w].mean())

    plt.subplot(122)
    plt.imshow(result)
    for i in range(18): plt.text(mxs[i], mys[i], segs[i], color='k', horizontalalignment='center',verticalalignment='center')
 
    hdu = fits.PrimaryHDU((result*im).astype(np.uint8))
    for i in range(18):
        hdu.header.update('PMSA_'+str(i+1), segs[i])

    # TODO copy relevant keywords and history from input FITS file header
    hdu.writeto("JWpupil_segments.fits", clobber=True)


#--------------------------------------------------------------------------------

def test_OPDbender():

    plt.figure(1)
    tel =OPDbender()
    tel.displace('A1', 1,0,0,draw=False)
    tel.displace('A2', 0,1,0,draw=False)
    tel.displace('A3', 0,0,.03,draw=False)
    tel.displace('A4', 0,-1,0,draw=False)
    tel.displace('A5', 1,-1,0,draw=False)

    tel.tilt('B1',.1,0,0)
    tel.tilt('B2',0,.1,0)
    tel.tilt('B3',0,0,100)

    tel.draw()

    plt.figure(2)
    tel.zern_seg('B3')


    print("\n\n\n")
    print("RMS WFE is ", tel.rms())

    tel.print_state()


def test2_OPDbender(filename='OPD_RevV_nircam_132.fits'):
    orig = OPDbender(filename)

    plot_kwargs = {'colorbar_orientation':'horizontal', 'clear': False}

    plt.clf()
    plt.subplot(131)
    orig.draw(title='Input OPD from \n'+filename, **plot_kwargs)

    perturbed = orig.copy()
    perturbed.perturb_all(multiplier=0.2, draw=False)

    plt.subplot(132)
    perturbed.draw(title='OPD after small random perturbation', **plot_kwargs)

    plt.subplot(133)
    diff = (perturbed-orig)
    diff.draw(title='Difference ({0:.1f} nm rms)'.format(diff.rms()), **plot_kwargs)




################################################################################
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(name)-12s: %(levelname)-8s %(message)s',)

    pass

    ob =OPDbender('OPD_RevV_nircam_132.fits', rm_piston_tilt=False)

