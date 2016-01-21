
""" jwxml: Various Python classes for parsing JWST-related information in XML files

* SUR: a segment update request file (mirror move command from the WAS to the MCS)
* Update: a single mirror update inside of a SUR
* SIAF: a SIAF file (Science Instrument Aperture File, listing the defined apertures for a given instrument)
* Aperture: a single aperture inside a SIAF


"""
import numpy as np
import matplotlib.pyplot as plt
try:
    from lxml import etree
    HAVE_LXML = True
except ImportError:
    import xml.etree.cElementTree as etree
    HAVE_LXML = False

import logging
import unittest
import os
_log = logging.getLogger('jwxml')
from . import conf


#---------------------------------------------------------------------------------
#   Mirror Move related classes

class Segment_Update(object):
    """ Class for representing one single mirror update (will be inside of groups in SURs)
    """
    def __init__(self, xmlnode):
        if xmlnode.attrib['type'] != 'pose': raise NotImplemented("Only Pose updates supported yet")

        self.id = int(xmlnode.attrib['id'])
        self.type = xmlnode.attrib['type']
        self.segment = xmlnode.attrib['seg_id'][0:2]
        self.absolute = xmlnode.attrib['absolute'] =='true'
        self.coord= xmlnode.attrib['coord'] #local or global
        self.stage_type= xmlnode.attrib['stage_type']  # recenter_fine, fine_only, none

        self.units = dict()
        self.moves = dict()
        for move in iterchildren(xmlnode):
            #print(move.tag, move.text )
            self.moves[move.tag] =float(move.text)
            self.units[move.tag] = move.attrib['units']
            #X_TRANS, Y_TRANS, PISTON, X_TILT, Y_TILT, CLOCK
        #allowable units:
		#units="id"
		#units="meters"
		#units="none"
		#units="radians"
		#units="sag"
		#units="steps"
		#
        # pose moves will only ever have meters/radians as units
    def __str__(self):
        return ("Update %d, %s, %s: "% (self.id, 'absolute' if self.absolute else 'relative', self.coord)) + str(self.moves)
    def shortstr(self):
        outstr = ("Update %d: %s, %s, %s {"% (self.id, self.segment, 'absolute' if self.absolute else 'relative', self.coord))

        outstr+= ", ".join([ coordname+"=%.3g" % self.moves[coordname] for coordname in ['PISTON','X_TRANS','Y_TRANS','CLOCK', 'X_TILT','Y_TILT']])
        #for coordname in ['PISTON','X_TRANS','Y_TRANS','CLOCK', 'X_TILT','Y_TILT']:
            #outstr+=coordname+"=%.3g" % self.moves[coordname]
        outstr+="}"
        return outstr

    @property
    def xmltext(self):
        """ The XML text representation of a given move """
        text= '        <UPDATE id="{0.id}" type="{0.type}" seg_id="{0.segment}" absolute="{absolute}" coord="{0.coord}" stage_type="{0.stage_type}">\n'.format( self, absolute = str(self.absolute).lower())
        for key in ['X_TRANS','Y_TRANS','PISTON','X_TILT', 'Y_TILT', 'CLOCK']:
            if key in self.moves:
                text+='            <{key}  units="{unit}">{val:E}</{key}>\n'.format(key=key, unit=self.units[key], val=self.moves[key])
        text+= '        </UPDATE>\n'
        return text

    def toGlobal(self):
        """ Return moves cast to global coordinates """
        if self.coord =='global':
            return self.moves
        else:
            raise NotImplemented("Error")


    def toLocal(self):
        """ Return moves cast to local coordinates """
        if self.coord =='local':
            return self.moves
        else:
            raise NotImplemented("Error")
            # TO implement based on Ball's 'pmglobal_to_seg' in ./wfsc_core_algs/was_core_pmglobal_to_seg.pro
            # or the code in ./segment_control/mcs_hexapod_obj__define.pro


class SUR(object):
    """ Class for parsing/manipulating Segment Update Request files

    """
    def __init__(self, filename):
        """ Read a SUR from disk """
        self.filename=filename

        self._tree = etree.parse(filename)

        for tag in ['creator','date','time','version', 'operational']:
            self.__dict__[tag] = self._tree.getroot().attrib[tag]
        for element in self._tree.getroot().iter():
            if element.tag =='CONFIGURATION_NAME':  self.configuration_name = element.text
            if element.tag =='CORRECTION_ID':  self.correction_id = element.text

        self.groups = []
        for grp in self._tree.getroot().iter('GROUP'):
            myupdates = []
            for update in grp.iter('UPDATE'):
                myupdates.append(Segment_Update(update))
            self.groups.append(myupdates)

    def __str__(self):
        outstr = "SUR %s\n" % self.filename #, type=%s, coords=%s\n" % (self.filename, 'absolute' if self.absolute else 'relative', self.coord)
        for igrp, grp in enumerate(self.groups):
            outstr+= "\tGroup %d\n" % (igrp+1)
            for update in grp:
                outstr+= "\t\t"+str(update)+"\n"
        return outstr

    @property
    def xmltext(self):
        """ The XML text representation of a given move """
        text = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<SEGMENT_UPDATE_REQUEST creator="?" date="{date}" time="{time}" version="0.0.1" operational="false" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="../../setup_files/schema/segment_update_request.xsd">
    <CONFIGURATION_NAME>{self.configuration_name}</CONFIGURATION_NAME>
    <CORRECTION_ID>{self.correction_id}</CORRECTION_ID>\n""".format(self=self, date='YYYY-MM-DD', time='HH:MM:SS')
    # FIXME add date and time keywords for real
        for igrp, grp in enumerate(self.groups):
            text+='    <GROUP id="{id}">\n'.format(id=igrp+1)
            for update in grp:
                text+=update.xmltext
            text+='    </GROUP>\n'
        text+= '</SEGMENT_UPDATE_REQUEST>'
        return text


    #@property
    #def name(self): return self._tree.getroot().attrib['name']


#---------------------------------------------------------------------------------
#  SIAF related classes

class Aperture(object):
    """ An Aperture, as parsed from the XML.
    All XML nodes are converted into object attributes.

    See JWST-STScI-001550 for the reference on which this implementation was based.

    4 Coordinate systems:
        * Detector:  pixels, in raw detector read out axes orientation ("Det")
        * Science:   pixels, in conventional DMS axes orientation ("Sci")
        * Ideal:     arcsecs relative to aperture reference location. ("Idl")
        * Telescope: arcsecs V2,V3 ("Tel")


    Example
    ========

    ap = some_siaf['desired_aperture_name']     # extract one aperture from a SIAF

    ap.Det2Tel(1024, 512)                       # convert pixel coordinates to sky Tel coords.
                                                # takes pixel coords, returns arcsec
    ap.Idl2Sci( 10, 3)                          # convert Idl coords to Sci pixels
                                                # takes arcsec, returns pixels

    # there exist functions for all of the possible {Tel,Idl,Sci,Det}2{Tel,Idl,Sci,Det} combinations.
    # you can also specify frames by string:
    ap.convert(1024, 512, frame_from='Det', frame_to='Tel')  # same as first example above

    ap.corners('Tel')                           # Get detector corners in Tel frame
    ap.center('Tel')                            # Get the reference point defined in the SIAF
                                                # this is typically the center of this region

    ap.plot('Idl', annotate=True, title=True)   # plot coordinates in Idl frame
    ap.plotDetectorChannels()                   # color in the readout channels


    """
    def __init__(self, xmlnode, instrument=None):

        self.instrument=instrument
        convfactors = {'RADIANS': 1, 'DEGREES': np.pi/180, 'ARCSECS': np.pi/180/60/60}

        for node in iterchildren(xmlnode):
            tag = node.tag.replace('{http://www.stsci.edu/SIAF}','')
            if len(node.getchildren()) ==0:
                # if doens't have children,
                try:
                    value = float(node.text) # do we care about ints vs floats?
                except (ValueError,TypeError):
                    value=node.text
                self.__dict__[tag] = value
            else:
                # if does have children:
                if '{http://www.stsci.edu/SIAF}units' in [c.tag for c in node.getchildren()]:
                    # this will be an angle/units pair. units are either in arcsec or degrees. Convert to radians in any case for internal use.
                    unit  = node.find('{http://www.stsci.edu/SIAF}units').text
                    value =float( node.find('{http://www.stsci.edu/SIAF}value').text) * convfactors[unit]
                    self.__dict__[tag] = value
                elif '{http://www.stsci.edu/SIAF}elt' in [c.tag for c in node.getchildren()]:
                    #  an array of values which should go to an NDarray
                    elts = [float(c.text) for c in iterchildren(node, '{http://www.stsci.edu/SIAF}elt')]
                    self.__dict__[tag] = np.asarray(elts)

                else:
                    raise NotImplemented("Not sure how to parse that node.")

        # pack things into NDarrays for convenient access
        # first the vertices
        self.XIdlVert = np.asarray((self.XIdlVert1, self.XIdlVert2,self.XIdlVert3,self.XIdlVert4))
        self.YIdlVert = np.asarray((self.YIdlVert1, self.YIdlVert2,self.YIdlVert3,self.YIdlVert4))

        # then the transformation coefficients

        if self.Sci2IdlDeg is not None:
            self.Sci2IdlDeg = int(self.Sci2IdlDeg)
            self.Sci2IdlCoeffs_X = np.zeros( (self.Sci2IdlDeg+1, self.Sci2IdlDeg+1))
            self.Sci2IdlCoeffs_Y = np.zeros( (self.Sci2IdlDeg+1, self.Sci2IdlDeg+1))
            self.Idl2SciCoeffs_X = np.zeros( (self.Sci2IdlDeg+1, self.Sci2IdlDeg+1))
            self.Idl2SciCoeffs_Y = np.zeros( (self.Sci2IdlDeg+1, self.Sci2IdlDeg+1))
            for i in range(1,self.Sci2IdlDeg+1):
                for j in range(0,i+1):
                    #if self.AperName == 'FGS2_FULL_CNTR':
                        #print('Sci2IdlX{0:1d}{1:1d}'.format(i,j), self.__dict__['Sci2IdlX{0:1d}{1:1d}'.format(i,j)])
                    self.Sci2IdlCoeffs_X[i,j] = self.__dict__['Sci2IdlX{0:1d}{1:1d}'.format(i,j)]
                    self.Sci2IdlCoeffs_Y[i,j] = self.__dict__['Sci2IdlY{0:1d}{1:1d}'.format(i,j)]
                    self.Idl2SciCoeffs_X[i,j] = self.__dict__['Idl2SciX{0:1d}{1:1d}'.format(i,j)]
                    self.Idl2SciCoeffs_Y[i,j] = self.__dict__['Idl2SciY{0:1d}{1:1d}'.format(i,j)]

    def __repr__(self):
        return "<jwxml.Aperture object AperName={0} >".format(self.AperName)


    #--- the actual fundamental transformation code follows in these next routines:
    def Det2Sci(self, XDet, YDet):
        """ Detector to Science, following Section 4.1 of JWST-STScI-001550"""
        XDet = np.asarray(XDet, dtype=float)
        YDet = np.asarray(YDet, dtype=float)
        ang = np.deg2rad(self.DetSciYAngle)
        XSci = self.XSciRef + self.DetSciParity* ((XDet - self.XDetRef)* np.cos(ang) + (YDet-self.YDetRef) * np.sin(ang))
        YSci = self.YSciRef -                     (XDet - self.XDetRef)* np.sin(ang) + (YDet-self.YDetRef) * np.cos(ang)
        return XSci, YSci

    def Sci2Det(self, XSci, YSci):
        """ Science to Detector, following Section 4.1 of JWST-STScI-001550"""
        XSci = np.asarray(XSci, dtype=float)
        YSci = np.asarray(YSci, dtype=float)

        ang = np.deg2rad(self.DetSciYAngle)
        XDet = self.XDetRef + self.DetSciParity * (XSci - self.XSciRef ) * np.cos(ang) - (YSci - self.YSciRef ) * np.sin(ang)
        YDet = self.YDetRef + self.DetSciParity * (XSci - self.XSciRef ) * np.sin(ang) + (YSci - self.YSciRef ) * np.cos(ang)
        return XDet, YDet

    def Sci2Idl(self, XSci, YSci):
        """ Convert Sci to Idl
        input in pixel, output in arcsec """
        dX = np.asarray(XSci, dtype=float) - self.XSciRef
        dY = np.asarray(YSci, dtype=float) - self.YSciRef

        degree = self.Sci2IdlDeg
        #CX = self.Sci2IdlCoefX
        #CY = self.Sci2IdlCoefY

        #XIdl = CX[0]*dX + CX[1]*dY + CX[2]*dX**2 + CX[3]*dX*dY + CX[4]*dY**2
        #YIdl = CY[0]*dY + CY[1]*dY + CY[2]*dY**2 + CY[3]*dY*dY + CY[4]*dY**2
        XIdl = np.zeros_like(np.asarray(XSci), dtype=float)
        YIdl = np.zeros_like(np.asarray(YSci), dtype=float)

        for i in range(1,degree+1):
            for j in range(0,i+1):
                XIdl += self.Sci2IdlCoeffs_X[i,j] * dX**(i-j) * dY**j
                YIdl += self.Sci2IdlCoeffs_Y[i,j] * dX**(i-j) * dY**j


        return XIdl, YIdl

    def Idl2Sci(self, XIdl, YIdl):
        """ Convert Idl to  Sci
        input in arcsec, output in pixels """
        XIdl = np.asarray(XIdl, dtype=float)
        YIdl = np.asarray(YIdl, dtype=float)

        degree = self.Sci2IdlDeg
        #dX = XIdl #Idl origin is by definition 0
        #dY = YIdl #Idl origin is by definition 0

        XSci = np.zeros_like(np.asarray(XIdl), dtype=float)
        YSci = np.zeros_like(np.asarray(YIdl), dtype=float)

        for i in range(1,degree+1):
            for j in range(0,i+1):
                XSci += self.Idl2SciCoeffs_X[i,j] * XIdl**(i-j) * YIdl**j
                YSci += self.Idl2SciCoeffs_Y[i,j] * XIdl**(i-j) * YIdl**j



        #CX = self.Idl2SciCoefX
        #CY = self.Idl2SciCoefY

        #XSci = CX[0]*dX + CX[1]*dY + CX[2]*dX**2 + CX[3]*dX*dY + CX[4]*dY**2
        #YSci = CY[0]*dY + CY[1]*dY + CY[2]*dY**2 + CY[3]*dY*dY + CY[4]*dY**2
        return XSci + self.XSciRef, YSci + self.YSciRef
        #return XSci, YSci

    def Idl2Tel(self, XIdl, YIdl):
        """ Convert Idl to  Tel

        input in arcsec, output in arcsec

        WARNING
        --------
        This is an implementation of the planar approximation, which is adequate for most
        purposes but may not be for all. Error is about 1.7 mas at 10 arcminutes from the tangent
        point. See JWST-STScI-1550 for more details.
        """
        XIdl = np.asarray(XIdl, dtype=float)
        YIdl = np.asarray(YIdl, dtype=float)

        #print(self.V2Ref, self.V3Ref)
        #rad2arcsec = 1./(np.pi/180/60/60)

        #V2Ref and V3Ref are now in arcseconds in the XML file
        ang = np.deg2rad(self.V3IdlYAngle)
        V2 = self.V2Ref + self.VIdlParity * XIdl * np.cos(ang) + YIdl * np.sin(ang)
        V3 = self.V3Ref - self.VIdlParity * XIdl * np.sin(ang) + YIdl * np.cos(ang)
        return V2, V3

    def Tel2Idl(self,V2, V3):
        """ Convert Tel to Idl

        input in arcsec, output in arcsec

        This transformation involves going from global V2,V3 to local angles with respect to some
        reference point, and possibly rotating the axes and/or flipping the parity of the X axis.


        WARNING
        --------
        This is an implementation of the planar approximation, which is adequate for most
        purposes but may not be for all. Error is about 1.7 mas at 10 arcminutes from the tangent
        point. See JWST-STScI-1550 for more details.
        """

        #rad2arcsec = 1./(np.pi/180/60/60)
        dV2 = np.asarray(V2, dtype=float)-self.V2Ref
        dV3 = np.asarray(V3, dtype=float)-self.V3Ref
        ang = np.deg2rad(self.V3IdlYAng)

        XIdl = self.VIdlParity * (dV2 * np.cos(ang) - dV2 * np.sin(ang))
        YIdl =                    dV2 * np.sin(ang) + dV3 * np.cos(ang)
        return XIdl, YIdl

    #--- and now some compound transformations that are less fundamental. This just nests calls to the above.

    def Det2Idl(self, *args):
        return self.Sci2Idl(*self.Det2Sci(*args))
    def Det2Tel(self, *args):
        return self.Idl2Tel(*self.Sci2Idl(*self.Det2Sci(*args)))
    def Sci2Tel(self, *args):
        return self.Idl2Tel(*self.Sci2Idl(*args))
    def Idl2Det(self, *args):
        return self.Sci2Det(*self.Idl2Sci(*args))
    def Tel2Sci(self, *args):
        return self.Idl2Sci(*self.Tel2Idl(*args))
    def Tel2Det(self, *args):
        return self.Sci2Det(*self.Idl2Sci(*self.Tel2Idl(*args)))

    #--- now, functions other than direct coordinate transformations
    def convert(self, X, Y, frame_from=None, frame_to=None):
        """ Generic conversion routine, that calls one of the
        specific conversion routines based on the provided frame names as strings. """

        if frame_from is None: raise ValueError("You must specify a frame_from value : Tel, Idl, Sci, Det")
        if frame_to is None: raise ValueError("You must specify a frame_to value : Tel, Idl, Sci, Det")


        if frame_from == frame_to: return X, Y  # null transformation

        #frames = ['Det','Sci', 'Idl','Tel']
        function = eval('self.%s2%s' % (frame_from, frame_to))
        return function(X,Y)

    def corners(self, frame='Idl'):
        " Return coordinates of the aperture outline"
        return self.convert(self.XIdlVert, self.YIdlVert, 'Idl', frame)

    def center(self, frame='Tel'):
        """ Return the defining center point of the aperture"""
        return self.convert(self.V2Ref, self.V3Ref, 'Tel', frame)


    def plot(self, frame='Idl', label=True, ax=None, title=True, units='arcsec', annotate=False, color=None):
        """ Plot this one aperture

        Parameters
        -----------
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
        label : bool
            Add text label stating aperture name
        units : str
            one of 'arcsec', 'arcmin', 'deg'
        annotate : bool
            Add annotations for detector (0,0) pixels
        title : str
            If set, add a label to the plot indicating which frame was plotted.

        """


        if units is None:
            units='arcsec'

        # should we flip the X axis direction at the end of this function?
        need_to_flip_axis = False # only flip if we created the axis
        if ax is None:
            ax = plt.gca()
            ax.set_aspect('equal')
            if frame=='Idl' or frame=='Tel':
                need_to_flip_axis = True # *and* we're displaying some coordinates in angles relative to V2.
                ax.set_xlabel('V2 [{0}]'.format(units))
                ax.set_ylabel('V3 [{0}]'.format(units))

            elif frame=='Sci' or frame=='Det':
                ax.set_xlabel('X pixels [{0}]'.format(frame))
                ax.set_ylabel('Y pixels [{0}]'.format(frame))



        x, y = self.corners(frame=frame)

        if units.lower() == 'arcsec':
            scale=1
        elif units.lower() =='arcmin':
            scale=01./60
        elif units.lower() =='deg':
            scale=01./60/60
        else:
            raise ValueError("Unknown units: "+units)


        x2 = np.concatenate([x, [x[0]]]) # close the box
        y2 = np.concatenate([y, [y[0]]])
        #print((x2,y2))
        ax.plot(x2*scale,y2*scale, color=color) # convert arcsec to arcmin


        if need_to_flip_axis:
            #print("flipped x axis")
            #ax.set_xlim(ax.get_xlim()[::-1])
            pass

        if label:
            rotation = 30 if self.AperName.startswith('NRC') else 0 # partially mitigate overlapping NIRCam labels
            ax.text(x.mean()*scale, y.mean()*scale, self.AperName,
                verticalalignment='center', horizontalalignment='center', rotation=rotation,
                color=ax.lines[-1].get_color())
        if title:
            ax.set_title("{0} frame".format(frame))
        if annotate:
            self.plotDetectorOrigin(frame=frame)

    def plotDetectorOrigin(self, frame='Idl', which='both'):
        """ Draw red and blue squares to indicate the raw detector
        readout and science frame readout, respectively

        Parameters
        -----------
        which : str
            Which detector origin to plot: 'both', 'Det', 'Sci'
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
        """

        # raw detector frame
        if which.lower() == 'det' or which.lower()=='both':
            c1, c2 = self.convert( 0, 0, 'Det', frame)
            plt.plot(c1, c2, color='red', marker='s', markersize=9)

        # science frame
        if which.lower() == 'sci' or which.lower()=='both':
            c1, c2 = self.convert( 0, 0, 'Sci', frame)
            plt.plot(c1, c2, color='blue', marker='s')

    def plotDetectorChannels(self, frame='Idl', color='0.5', alpha=0.3, evenoddratio=0.5, **kwargs):
        """ Mark on the plot the various detector readout channels

        These are depicted as alternating light/dark bars to show the
        regions read out by each of the output amps.

        Parameters
        ----------
        frame : str
            Optional if you have already called plot() to specify a
            coordinate frame.
       """

        import matplotlib
        if self.instrument == 'MIRI': npixels = 1024
        else: npixels = 2048
        ch = npixels/4

        ax = plt.gca()
        pts = ((0, 0), (ch,0), (ch,npixels), (0, npixels))
        for chan in range(4):
            plotpoints = np.zeros((4,2))
            for i,xy in enumerate(pts):
                plotpoints[i] = self.convert(xy[0]+chan*ch,xy[1],'Det',frame)
            rect = matplotlib.patches.Polygon(plotpoints, closed=True,
                    alpha=(alpha if np.mod(chan,2) else alpha*evenoddratio),
                    facecolor=color, edgecolor='none', lw=0)
            ax.add_patch(rect)


class SIAF(object):
    """ Science Instrument Aperture File

    This is a class interface to SIAF information stored in an XML file.
    It lets you read (only) the SIAF information, retrieve apertures,
    plot them, and transform coordinates accordingly.

    This class is basically just a container. See the Aperture class for
    the detailed implementation of the transformations.

    Briefly, this class acts like a dict containing Aperture objects, accessible
    using their names defined in the SIAF

    Examples
    ---------

    fgs_siaf = SIAF('FGS')
    fgs_siaf.apernames                # returns a list of aperture names
    ap = fgs_siaf['FGS1_FULL_CNTR']   # returns an aperture object
    ap.plot(frame='Tel')              # plot one aperture
    fgs_siaf.plot()                   # plot all apertures in this file

    """
    def __init__(self, instr='NIRISS', basepath=None, filename=None, **kwargs):
            #basepath="/Users/mperrin/Dropbox/JWST/Optics Documents/SIAF/"
        """ Read a SIAF from disk

        Parameters
        -----------
        instr : string
            one of 'NIRCam', 'NIRSpec', 'NIRISS', 'MIRI', 'FGS'; case sensitive.
        basepath : string
            Directory to look in for SIAF files
        filename : string, optional
            Alternative method to specify a specific SIAF XML file.
        """

        if instr not in ['NIRCam', 'NIRSpec', 'NIRISS', 'MIRI', 'FGS']:
            raise ValueError("Invalid instrument name: {0}. Note that this is case sensitive.".format(instr))

        self.instrument=instr

        if basepath is None:
            from utils import get_webbpsf_data_path
            basepath = os.path.join(get_webbpsf_data_path(), instr)

        if filename is None:
            self.filename=os.path.join(basepath, instr+'_SIAF.xml')
        else:
            self.filename = filename

        self.apertures = {}

        self._tree = etree.parse(self.filename)

        self._last_plot_frame=None



        #for entry in self._tree.getroot().iter('{http://www.stsci.edu/SIAF}SiafEntry'):
        for entry in self._tree.getroot().iter('SiafEntry'):
            aperture = Aperture(entry, instrument=self.instrument)
            self.apertures[aperture.AperName] = aperture

    def __getitem__(self, key):
        return self.apertures[key]

    def __len__(self):
        return len(self.apertures)

    @property
    def apernames(self):
        """ List of aperture names defined in this SIAF"""
        return self.apertures.keys()


    def _getFullApertures(self):
        """ Return whichever subset of apertures correspond to the entire detectors. This is a
        helper function for the various plotting routines following"""
        fullaps = []
        if self.instrument =='NIRCam':
            fullaps.append( self.apertures['NRCA5_FULL'])
            fullaps.append( self.apertures['NRCB5_FULL'])
            #for letter in ['A', 'B']:
                #for number in range(1,6):
                    #fullaps.append(self.apertures['NRC{letter}{number}_FULL_CNTR'.format(letter=letter, number=number)])
        elif self.instrument =='NIRSpec':
            #fullaps.append( self.apertures['NRS1_FULL'])
            #fullaps.append( self.apertures['NRS2_FULL'])
            fullaps.append( self.apertures['NRS_FULL_MSA1'])
            fullaps.append( self.apertures['NRS_FULL_MSA2'])
            fullaps.append( self.apertures['NRS_FULL_MSA3'])
            fullaps.append( self.apertures['NRS_FULL_MSA4'])
        elif self.instrument =='NIRISS':
            fullaps.append( self.apertures['NIS-CEN'])
        elif self.instrument =='MIRI':
            fullaps.append( self.apertures['MIRIM_FULL_CNTR'])
        elif self.instrument =='FGS':
            fullaps.append( self.apertures['FGS1_FULL'])
            fullaps.append( self.apertures['FGS2_FULL'])
        return fullaps


    def plot(self, frame='Tel', names=None, label=True, units=None, clear=True, annotate=False,
        subarrays=True):
        """ Plot all apertures in this SIAF

        Parameters
        -----------
        names : list of strings
            A subset of aperture names, if you wish to plot only a subset
        subarrays : bool
            Plot all the minor subarrays if True, else just plot the "main" apertures
        label : bool
            Add text labels stating aperture names
        units : str
            one of 'arcsec', 'arcmin', 'deg'
        clear : bool
            Clear plot before plotting (set to false to overplot)
        annotate : bool
            Add annotations for detector (0,0) pixels
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
        """
        if clear: plt.clf()
        ax = plt.subplot(111)
        ax.set_aspect('equal')

        # which list of apertures to iterate over?
        if subarrays:
            iterable = self.apertures.itervalues
        else:
            iterable = self._getFullApertures

        for ap in iterable():
            if names is not None:
                if ap.AperName not in names: continue

            ap.plot(frame=frame, label=label, ax=ax, units=None)
            if annotate:
                ap.plotDetectorOrigin(frame=frame)
        ax.set_xlabel('V2 [arcsec]')
        ax.set_ylabel('V3 [arcsec]')



        if frame =='Tel' or frame=='Idl':
            # enforce V2 increasing toward the left
            ax.autoscale_view(True,True,True)
            xlim = ax.get_xlim()
            if xlim[1] > xlim[0]: ax.set_xlim(xlim[::-1])
            ax.set_autoscalex_on(True)

        self._last_plot_frame = frame

    def plotDetectorOrigin(self, which='both', frame=None):
        """ Mark on the plot the detector's origin in Det and Sci coordinates

        Parameters
        -----------
        which : str
            Which detector origin to plot: 'both', 'Det', 'Sci'
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
            Optional if you have already called plot() to specify a
            coordinate frame.

        """
        if frame is None: frame = self._last_plot_frame
        for ap in self._getFullApertures():
            ap.plotDetectorOrigin(frame=frame, which=which)


    def plotDetectorChannels(self, frame=None):
        """ Mark on the plot the various detector readout channels

        These are depicted as alternating light/dark bars to show the
        regions read out by each of the output amps.

        Parameters
        ----------
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
            Optional if you have already called plot() to specify a
            coordinate frame.

        """
        if frame is None: frame = self._last_plot_frame
        for ap in self._getFullApertures():
            ap.plotDetectorChannels(frame=frame)


def plotAllSIAFs(subarrays = True, showorigin=True, showchannels=True, **kwargs):
    """ Plot All instrument """

    for instr in ['NIRCam','NIRISS','NIRSpec','FGS','MIRI']:
        aps =SIAF(instr, **kwargs)
        print("{0} has {1} apertures".format(aps.instrument, len(aps)))

        aps.plot(clear=False, subarrays=subarrays, **kwargs)
        if showorigin: aps.plotDetectorOrigin()
        if showchannels: aps.plotDetectorChannels()


def plotMainSIAFs(showorigin=False, showchannels=False, label=False, **kwargs):
    col_imaging = 'blue'
    col_coron = 'green'
    col_msa = 'magenta'

    nircam = SIAF('NIRCam')
    niriss= SIAF('NIRISS')
    fgs = SIAF('FGS')
    nirspec = SIAF('NIRSpec')
    miri = SIAF('MIRI')

    im_aps = [ nircam['NRCA5_FULL'], nircam['NRCB5_FULL'], niriss['NIS-CEN'], miri['MIRIM_FULL_ILLCNTR'],
            fgs['FGS1_FULL'], fgs['FGS2_FULL']]

    coron_aps = [nircam['NRCA2_MASK210R'], nircam['NRCA4_MASKSWB'],
            nircam['NRCA5_MASK335R'],
            nircam['NRCA5_MASK430R'],
            nircam['NRCA5_MASKLWB'],
            nircam['NRCB3_MASKSWB'],
            nircam['NRCB1_MASK210R'],
            nircam['NRCB5_MASK335R'],
            nircam['NRCB5_MASK430R'],
            nircam['NRCB5_MASKLWB'],
            miri['MIRIM_MASK1065_CNTR'],
            miri['MIRIM_MASK1140_CNTR'],
            miri['MIRIM_MASK1550_CNTR'],
            miri['MIRIM_MASKLYOT_CNTR']]
    msa_aps = [nirspec['NRS_FULL_MSA'+str(n+1)] for n in range(4)]

    for aplist, col in zip(  [im_aps, coron_aps, msa_aps], [col_imaging, col_coron, col_msa]):
        for ap in aplist:
            ap.plot(color=col, frame='Tel', label=label, **kwargs)


    # ensure V2 increases to the left
            #ax.set_xlim(
    ax = plt.gca()
    xlim = ax.get_xlim()
    if xlim[0] < xlim[1]:
        ax.set_xlim(xlim[::-1])
class Test_SIAF(unittest.TestCase):

    def assertAlmostEqualTwo(self, tuple1, tuple2):
        self.assertAlmostEqual(tuple1[0], tuple2[0], places=1)
        self.assertAlmostEqual(tuple1[1], tuple2[1], places=1)


    def _test_up(self):
        siaf = SIAF("JwstSiaf-2010-10-05.xml")
        startx = 1023
        starty = 1024

        nca = siaf['NIRCAM A']

        self.assertAlmostEqualTwo( nca.Det2Sci(startx,starty), (1020.,1020.))
        print("Det2Sci OK")

        self.assertAlmostEqualTwo( nca.Det2Idl(startx,starty), (0.0, 0.0))
        print("Det2Idl OK")
        self.assertAlmostEqualTwo( nca.Det2Tel(startx,starty), (87.50, -497.10))
        print("Det2Tel OK")

    def _test_down(self):
        siaf = SIAF("JwstSiaf-2010-10-05.xml")
        startV2 = 87.50
        startV3 = -497.10
        nca = siaf['NIRCAM A']

        self.assertAlmostEqualTwo( nca.Sci2Det(1020., 1020), (1023.,1024.))
        print("Sci2Det OK")

        self.assertAlmostEqualTwo( nca.Tel2Idl(startV2, startV3), (0.0, 0.0))
        print("Tel2Idl OK")
        self.assertAlmostEqualTwo( nca.Tel2Sci(startV2, startV3), (1020., 1020.))
        print("Tel2Sci OK")
        self.assertAlmostEqualTwo( nca.Tel2Det(startV2, startV3), (1023.,1024.))
        print("Tel2Det OK")

    def test_inverses(self):
        siaf = SIAF("JwstSiaf-2010-10-05.xml")
        nca = siaf['NIRCAM A']

        self.assertAlmostEqualTwo( nca.Det2Sci(*nca.Sci2Det(1020., 1020)), (1020., 1020) )
        self.assertAlmostEqualTwo( nca.Sci2Det(*nca.Det2Sci(1020., 1020)), (1020., 1020) )
        print("Det <-> Sci OK")

        self.assertAlmostEqualTwo( nca.Tel2Idl(*nca.Idl2Tel(10., 10)), (10., 10) )
        self.assertAlmostEqualTwo( nca.Idl2Tel(*nca.Tel2Idl(10., 10)), (10., 10) )
        print("Tel <-> Idl OK")

        self.assertAlmostEqualTwo( nca.Tel2Sci(*nca.Sci2Tel(10., 10)), (10., 10) )
        self.assertAlmostEqualTwo( nca.Sci2Tel(*nca.Tel2Sci(10., 10)), (10., 10) )
        print("Tel <-> Sci OK")


# The ElementTree implementation in xml.etree does not support
# Element.iterchildren, so provide this wrapper instead
# This wrapper does not currently provide full support for all the arguments as
# lxml's iterchildren
def iterchildren(element, tag=None):
    if HAVE_LXML:
        return element.iterchildren(tag)
    else:

        if tag is None:
            return iter(element)

        def _iterchildren():
            for child in element:
               if child.tag == tag:
                   yield child

        return _iterchildren()



if __name__== "__main__":
    logging.basicConfig(level=logging.DEBUG,format='%(name)-10s: %(levelname)-8s %(message)s')
    #unittest.main()
    #sur = SUR('/itar/jwst/wss/MMS_Delivery/09_MMS_Source_Code/wfsc_mcs~1.1.2/wfsc_mcs/fqt/tc7/sur_ok_rel_gl.xml')

    s = SIAF()



