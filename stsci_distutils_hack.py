#
# $HeadURL: https://www.stsci.edu/svn/ssb/stsci_python/stsci_python/trunk/pytools/lib/stsci_distutils_hack.py $
# $Rev: 7917 $
#
# Implements setup.py code common to many of our packages.
#
# The new standard stsci module setup.py is just
#
#   import pytools.stsci_distutils_hack
#   pytools.stsci_distutils_hack.run( pytools_version = "XX" )
#
# where XX is the version of pytools you expect for the install to work
#

from __future__ import division # confidence high

"""
Special handling for stsci_python package installation.

stsci_python is distributed as a single package, but it contains
packages that are also distributed separately.  When we use this
module to install our package, we can use the exact same definition
file to control the setup.py of the individual package _and_ the
setup.py of stsci_python.

This module also preserves revision control data in the installed
or distributed files.

If you are not a developer at STScI, this module is probably not of
much interest to you.

"""

__docformat__ = 'restructuredtext'


######## ######## ######## ######## ######## ######## ######## ########
#
# actually perform the install
#
# NOTE: This is not used to install pytools itself!

import sys

def run( pytools_version = None ) :
    """
    Perform a stsci_python install based on the information in defsetup.py

    * gather our subversion revision number and the install time

    * perform the install

    usage:

        import pytools.stsci_distutils_hack
        pytools.stsci_distutils_hack.run(pytools_version = "3.0")

    """

    if not hasattr(sys, 'version_info') or sys.version_info < (2,3,0,'alpha',0):
        raise SystemExit, "Python 2.3 or later required."

    if pytools_version :
        # Only try to import pytools if we are asked to check for a version.
        #
        # ( We may have been extracted from pytools and bundled with a package.
        # In that case, we do not want to risk finding some _other_ pytools
        # and comparing that version. )
        import pytools

        # bug: should use distutils version comparator to perform ">" comparisons
        if ( pytools.__version__ != pytools_version ) :
            print "wrong version of pytools!"
            print "have ",pytools.__version__
            print "want ",pytools_version
            sys.exit(1)


    from distutils.core import setup
    from defsetup import setupargs, pkg

    # collect our subversion information
    __set_svn_version__()

    # save the date when we last ran setup.py
    __set_setup_date__()

    if "version" in sys.argv :
        sys.exit(0)

    # If they have multiple packages, we have to allow them to give a list.
    # That is the unusual case, so we let them give a string if they have a single
    # package.
    if isinstance(pkg,str) :
        pkg = [ pkg ]

    # If they have multiple packages, they have to specify package_dir.  Otherwise,
    # we can create one for them.
    if not 'package_dir' in setupargs :
        setupargs['package_dir'] = { pkg[0] : 'lib' }

    setup(
        name =              pkg[0],
        packages =          pkg,

        **setupargs

        )



######## ######## ######## ######## ######## ######## ######## ########
#
# This part fixes install_data to put data files in the same directory
# with the python library files, which is where our packages want
# them.
#
# This is essentially "smart_install_data" as used in the old
# setup.py files, except that it also understands wildcards
# and os-specific paths.  This means the module author can
# ask for data files with
#       "data/generic/*"
# instead of
#       glob.glob(os.path.join('data', 'generic', '*'))


import os
import glob

import distutils.util

import distutils.command.install_data

o =  distutils.command.install_data.install_data

# same trick as smart_install_data used: save the old run() method and
# insert our own run method ahead of it

o.old_run = o.run

def new_run ( self ) :
        """
        Hack for distutils to cause install_data to be in the same directory
        as the python library files.  Our packages expect this.
        """

        # We want our data files in the directory with the library files
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')


        # self.data_files is a list of
        #       ( destination_directory, [ source_file, source_file, source_file ] )
        #
        # We want to do wildcard expansion on all the file names.
        #
        l = [ ]
        for f in self.data_files :
            ( dest_dir, files ) = f
            fl = [ ]
            for ff in files :
                ff = distutils.util.convert_path(ff)
                ff = glob.glob(ff)
                fl.extend(ff)
            dest_dir = distutils.util.convert_path(dest_dir)
            l.append( ( dest_dir, fl ) )
        self.data_files = l

        # now use the original run() function to finish
        return distutils.command.install_data.install_data.old_run(self)

o.run = new_run


######## ######## ######## ######## ######## ######## ######## ########
#
# Function to collect svn version information - used to be stsci_python/version.py
# with multiple copies in the system.
#
import os.path
import re

#
# This is the entry point.  All you need to do is call this function from your
# setup.py according to the example above.  It will create a file called
# lib/svn_version.py ;  After that, you can
#
#   # find out what subversion information applies to yourpackage
#   import yourpackage.svn_version
#   print yourpackage.svn_version.__svn_version__
#   print yourpackage.svn_version.__full_svn_info__
#

def __set_svn_version__(path="./", fname='svn_version.py' ) :
    #
    # path is the package where the version information will be stored.  Default
    # is "this package", but from a higher level package, you can specify a directory
    # of a package to process
    #
    # fname is the name of the file to store the version information in.  Never change
    # this.
    #

    info = None
    rev = __get_svn_rev__(path)
    version_file = os.path.join(path,'lib',fname)

    # if we are unable to determine the revision, we default to leaving the
    # revision file unchanged.  Otherwise, we fill it in with whatever
    # we have

    if rev is None:
        if os.path.exists(version_file) :
            return
        revision = 'Unable to determine SVN revision'
    else:
        if ( rev == 'exported' or rev == 'unknown' ) and os.path.exists(version_file) :
            return
        revision = str(rev)

    info = __get_full_info__(path)

    # now we can write the version information

    f = open(version_file,'w')
    f.write("__svn_version__ = %s\n" % repr(revision))

    # info will be a multi-line string.  We are not using repr(info)
    # for readability; the output of "svn info" can not contain '''
    # unless you are doing something bad.
    f.write("\n__full_svn_info__ = '''\n%s'''\n\n" % info)
    f.close()


def __get_svn_rev__(path):
    m = None
    try:
        # with popen3,  stderr goes into a pipe where we ignore it,
        # This means the user does not see errors.
        cmd = 'svnversion '+path
        (sin, sout, serr) = os.popen3(cmd)

        # pick up the first line of output
        m=sout.read().strip()

        # if it looks like valid svnversion output, return it
        if m == 'exported' :
            return m
        if re.match('^[0-9][0-9:]*[A-Z]*$',m) :
            return m

        # if we get here, it was not valid - that probably means
        # an error of some kind.
    except:
        pass

    return None

def __get_full_info__(path):
    info = None
    try:
        # with popen3,  stderr goes into a pipe where we ignore it,
        # This means the user does not see errors.
        (sin, sout, serr) = os.popen3('svn info %s' % path)

        # pick up all the lines of output
        info = [l.strip() for l in sout.readlines()]

        # if no output, there was an error and we don't know anything
        if len(info) == 0 :
            return "unknown"

        # there was output, so join it all together
        return '\n'.join(info)

    except:
        pass

    return "unknown"

######## ######## ######## ######## ######## ######## ######## ########
#
# note when we last ran setup.py -- what we really want is when the
# software was installed, but we can use the time we ran setup.py as
# a proxy for that.
#

def __set_setup_date__( path="./", fname='svn_version.py') :
    import datetime
    file = os.path.join(path,'lib',fname)
    d = datetime.datetime.now()
    l = [ ]
    try :
        # we don't expect this to fail ever, but it might
        f = open(file,"r")
        for line in f :
            if line.find("# setupdate") < 0 :
                l.append(line)
        f.close()
    except IOError :
        pass
    f=open(file,"w")
    for line in l :
        f.write(line)

    f.write("%s # setupdate\n" % "import datetime")
    f.write("%s # setupdate\n" % ("setupdate = "+repr(d)))
    f.close()

