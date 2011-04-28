#!/usr/bin/env python

# We use the local copy of stsci_distutils_hack, unless
# the user asks for the stpytools version

import os
if os.getenv("USE_STPYTOOLS") :
    import pytools.stsci_distutils_hack as H
    pytools_version = "3.0"
else :
    import stsci_distutils_hack as H
    pytools_version = None

H.run(pytools_version = pytools_version)


