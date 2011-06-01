.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Requirements & Installation
============================


Software Requirements
-----------------------

Beyond the usual numpy/scipy/matplotlib core modules, the following are required:

* `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits>`_
* `ATPy <http://atpy.github.com/>`_, which in turn requires `vo <https://trac6.assembla.com/astrolib>`_ and `asciitable <http://cxc.harvard.edu/contrib/asciitable/>`_
  
These are optional but highly recommended:

* `pyFFTW3 <http://pypi.python.org/pypi/PyFFTW3/0.2.1>`_. The code will work fine without it, but will be significantly slower.
* `pysynphot <https://trac6.assembla.com/astrolib>`_ enables the simulation of PSFs with proper spectral response to realistic source spectra.  Without this, PSF fidelity is reduced.

Additional requirement for the GUI: 

* The graphical user interface requires the ``ttk`` enhanced version of the ``Tkinter`` widget library. 

``ttk`` is not included by default on some installations of Python, for instance the default Mac OS Python 2.6 install. 
You may wish to either upgrade to a more current Python, or else compile and install ``ttk`` for your platform. ``WebbPSF``
was developed using Python 2.7, which includes ``ttk`` by default, but it ought to work fine on any installations of
Python 2.5 or 2.6 provided ``ttk`` is available. Alternatively, you can just skip using the GUI; the optical modeling classes
themselves have no dependency on these widgets.



Obtaining WebbPSF
-------------------------

Download the following files:

* `webbpsf-0.2.5.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-0.2.5.tar.gz>`_
* `webbpsf-data-0.2.5.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-data-0.2.5.tar.gz>`_  [417 MB]

Installing WebbPSF
--------------------

1. Untar ``webbpsf-0.2.tar.gz`` into a temporary working directory. 
2. Run ``python setup.py install`` in that directory. This will install ``webbpsf`` into your Python path. 
3. Untar ``webbpsf-data-0.2.tar.gz`` into a directory of your choosing.
4. Set the environment variable ``WEBBPSF_PATH`` to point to that directory. e.g. ``setenv WEBBPSF_PATH $HOME/data/webbpsf-data``.
5. You should now be able to do ``import webbpsf`` in a Python session. 

Installing or Updating pysynphot
-------------------------------
To install or update ``pysynphot``, do the following. (See also http://stsdas.stsci.edu/pysynphot/ and https://trac6.assembla.com/astrolib). If you already have ``pysynphot`` 
installed, it will probably work fine without this update, but computations may be slower than the current version. 

.. comment 
        work without this update but computations will be slower than the current version, so we recommend updating it. 

* Download `pysynphot-0.7jwst.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/pysynphot-0.7jwst.tar.gz>`_  

Note that the above-linked file is an unofficial, pre-release version, provided courtesy of STScI
Science Software Branch (Laidler, Greenfield, et al.) and not an "officially
supported" release. (This file tracks the internal development subversion
repository of pysynphot as of revision 2007 on April 21 2011). An official release of an updated pysynphot is expected in the near future.


.. comment
        you should still do these steps to update it to support all the JWST instruments transmission profiles. 

1. Untar ``pysynphot-0.7jwst.tar.gz`` into a temporary working directory. 
2. run ``python setup.py install`` in that directory.  You can delete the setup files there after you do this step. 
3. If this is your initial installation of ``pysynphot`` you need to install the CDBS files. See the `pysynphot installation guide <https://trac6.assembla.com/astrolib/wiki/PysynphotInstallationGuide>`_. The necessary files are available from https://trac6.assembla.com/astrolib; follow the download links for "throughput files" and "model spectra". If you already have CDBS installed, then you're all set and can skip this step.


WebbPSF includes its own normalized copies of the new JWST instrumental throughputs from the development CDBS at STScI.
If you have JWST throughput files available in your ``$PYSYN_CDBS`` directory (likely true only for internal users at STScI), those will be used in preference to the WebbPSF internal files, but
this is not required and so we do not distribute those files now.

.. comment
        3. Untar ``CDBS-for-webb.tar.gz`` in a directory of your choosing. (Typically replacing into your current CDBS directory if already present)
        4. Set the environment variable ``PYSYN_CDBS`` to point to that directory. e.g. ``setenv PYSYN_CDBS $HOME/data/CDBS``.


