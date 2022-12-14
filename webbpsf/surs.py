"""
Mirror Move related classes originally from jwxml/mirrors.py
"""

import logging
_log = logging.getLogger('webbpsf')

try:
    from lxml import etree

    def iterchildren(element, tag=None):
        return element.iterchildren(tag)
except ImportError:
    import xml.etree.cElementTree as etree

    # The ElementTree implementation in xml.etree does not support
    # Element.iterchildren, so provide this wrapper instead
    # This wrapper does not currently provide full support for all the arguments as
    # lxml's iterchildren
    def iterchildren(element, tag=None):
        if tag is None:
            return iter(element)

        def _iterchildren():
            for child in element:
                if child.tag == tag:
                    yield child

        return _iterchildren()


class SegmentUpdate(object):
    """
    Class for representing one single mirror update (will be inside of groups in SURs)
    Allowable units: "id", "meters", "none", "radians", "sag", "steps"
    Pose moves will only ever have meters/radians as units
    """
    def __init__(self, xmlnode=None):
        self.units = dict()
        self.moves = dict()

        if xmlnode is not None:
            # parse from some XML in a SUR file
            self.id = int(xmlnode.attrib['id'])
            self.type = xmlnode.attrib['type']
            self.segment = xmlnode.attrib['seg_id'][0:2]
            self.absolute = xmlnode.attrib['absolute'] == 'true'
            self.coord = xmlnode.attrib['coord']  # local or global
            self.stage_type = xmlnode.attrib['stage_type']  # recenter_fine, fine_only, none

            for move in iterchildren(xmlnode):
                self.moves[move.tag] = float(move.text)
                self.units[move.tag] = move.attrib['units']
        else:
            # create an empty SUR
            self.id = None
            self.type = None
            self.segment = None
            self.absolute = None
            self.coord = None
            self.stage_type = None

    def __str__(self):
        return ("Update %d, move %s, %s, %s: " % (self.id, self.segment, 'absolute' if self.absolute else 'relative', self.coord)) + \
               str(self.moves)

    def shortstr(self):
        outstr = ("Update %d: %s, %s, %s {" % (self.id, self.segment, 'absolute' if self.absolute else 'relative',
                                               self.coord))

        outstr += ", ".join([coordname+"=%.3g" % self.moves[coordname] for coordname in
                             ['PISTON', 'X_TRANS', 'Y_TRANS', 'CLOCK', 'X_TILT', 'Y_TILT']])
        outstr += "}"
        return outstr

    @property
    def xmltext(self):
        """ The XML text representation of a given move """
        text = '        <UPDATE id="{0.id}" type="{0.type}" seg_id="{0.segment}" absolute="{absolute}" ' \
               'coord="{0.coord}" stage_type="{0.stage_type}">\n'.format(self, absolute=str(self.absolute).lower())
        for key in ['X_TRANS', 'Y_TRANS', 'PISTON', 'X_TILT', 'Y_TILT', 'CLOCK']:
            if key in self.moves:
                text += '            <{key}  units="{unit}">{val:E}</{key}>\n'.format(key=key, unit=self.units[key],
                                                                                      val=self.moves[key])
        text += '        </UPDATE>\n'
        return text

    def to_global(self):
        """ Return moves cast to global coordinates """
        if self.coord == 'global':
            return self.moves
        else:
            raise NotImplemented("Error")

    def to_local(self):
        """ Return moves cast to local coordinates """
        if self.coord == 'local':
            return self.moves
        else:
            raise NotImplemented("Error")
            # TO implement based on Ball's 'pmglobal_to_seg' in ./wfsc_core_algs/was_core_pmglobal_to_seg.pro
            # or the code in ./segment_control/mcs_hexapod_obj__define.pro


class SUR(object):
    """ Class for parsing/manipulating Segment Update Request files
    """
    def __init__(self, filename=None):
        """ Read a SUR from disk """


        self.filename = filename
        self.groups = []

        if filename is not None:
            # Parse an existing file
            self._tree = etree.parse(filename)

            for tag in ['creator', 'date', 'time', 'version', 'operational']:
                self.__dict__[tag] = self._tree.getroot().attrib[tag]
            for element in self._tree.getroot().iter():
                if element.tag == 'CONFIGURATION_NAME':
                    self.configuration_name = element.text
                if element.tag == 'CORRECTION_ID':
                    self.correction_id = element.text

            for grp in self._tree.getroot().iter('GROUP'):
                myupdates = []
                for update in grp.iter('UPDATE'):
                    myupdates.append(SegmentUpdate(update))
                self.groups.append(myupdates)
        else:
            # Blank empty SUR
            self.correction_id = None
            self.configuration_name = None
            self.creator = 'WebbPSF'
            self.date = None
            self.time = None
            self.version = None
            self.operational = False
            self.groups.append([])

    def __str__(self):
        outstr = "SUR %s\n" % self.filename
        for igrp, grp in enumerate(self.groups):
            outstr += "\tGroup %d\n" % (igrp+1)
            for update in grp:
                outstr += "\t\t"+str(update)+"\n"
        return outstr

    @property
    def ngroups(self):
        return len(self.groups)
    @property
    def nmoves(self):
        return sum([len(g) for g in self.groups])

    def describe(self):
        return f"SUR with {self.ngroups} groups and {self.nmoves} moves"

    @property
    def xmltext(self):
        """ The XML text representation of a given move """
        text = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
                <SEGMENT_UPDATE_REQUEST creator="?" date="{date}" time="{time}" version="0.0.1" operational="false" 
                xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="../../setup_files/
                schema/segment_update_request.xsd">
                <CONFIGURATION_NAME>{self.configuration_name}</CONFIGURATION_NAME>
                <CORRECTION_ID>{self.correction_id}</CORRECTION_ID>\n""".format(self=self,
                                                                                date='YYYY-MM-DD', time='HH:MM:SS')
        # FIXME add date and time keywords for real
        for igrp, grp in enumerate(self.groups):
            text += '    <GROUP id="{id}">\n'.format(id=igrp+1)
            for update in grp:
                text += update.xmltext
            text += '    </GROUP>\n'
        text += '</SEGMENT_UPDATE_REQUEST>'
        return text
