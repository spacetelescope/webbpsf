.. JWST-PSFs documentation master file, created by
   sphinx-quickstart on Mon Nov 29 15:57:01 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Requirements & Installation
============================


Software Requirements
-----------------------

**Python**: Python 2.6 or higher is required. Python 2.7 is required for the GUI (see below) and strongly recommended overall. WebbPSF is not yet Python 3 compatible.


**Python modules**: Beyond the usual numpy/scipy/matplotlib core modules, the following are required. 

* **Either** `astropy <http://astropy.org>`_, in particular its ``astropy.io.fits`` and ``astropy.io.ascii`` modules, 
* **or** the following individual python modules to provide similar functionality:

        * `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits>`_
        * `asciitable <http://cxc.harvard.edu/contrib/asciitable/>`_

  
These are optional but recommended:

* `pysynphot <https://trac6.assembla.com/astrolib>`_ enables the simulation of PSFs with proper spectral response to realistic source spectra.  Without this, PSF fidelity is reduced. See below for :ref:`installation instructions for pysynphot <pysynphot_install>`. 
* `pyFFTW3 <http://pypi.python.org/pypi/PyFFTW3/0.2.1>`_. The FFTW library will significantly speed up the FFTs used in coronagraphic simulations. Since direct imaging simulations use a discrete matrix FT instead, direct imaging simulation speed is unchanged.  pyFFTW3 is highly recommended if you expect to perform many coronagraphic calculations.

**Additional requirement for the GUI:** The :ref:`graphical user interface<gui>` requires 

* **Either**  the `wxpython <http://www.wxpython.org>`_ interface to the ``wxwidgets`` widget library. 


* **or**  the `ttk <http://docs.python.org/2/library/ttk.html>`_ enhanced version of the ``Tkinter`` widget library. 
  (``ttk`` is not included by default on some installations of Python, for instance Mac OS 10.6's default of Python 2.6. 
  You may wish to either upgrade to a more current Python, or else compile and install ``ttk`` for your platform. ``WebbPSF``
  was developed using Python 2.7, which includes ``ttk`` by default, but it ought to work fine on any installations of
  Python 2.5 or 2.6 provided ``ttk`` is available.)

Similar GUIs are implemented in both widget tool kits, with the wxpython GUI
providing some additional functionality. Future development will likely
concentrate on the wxpython toolkit, but for now both are supported.

Alternatively, you can just skip using the GUI; the optical modeling classes
themselves have no dependency on these widgets.



Installing WebbPSF
----------------------

Installing WebbPSF via PYPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since version 0.2.8, WebbPSF and its underlying optical library ``poppy`` are both
installable from the `Python Package Index <http://pypi.python.org/pypi>`_ via
the standard toolchain using `pip
<http://www.pip-installer.org/en/latest/index.html>`_ or `easy_install <http://pypi.python.org/pypi/setuptools>`_.  This is the recommended installation
method if you already have a working copy of python, numpy, and matplotlib on your computer. 

1. Invoke pip in the usual manner::

   $ pip install webbpsf
   [... progress report ...]

   ``Successfully installed webbpsf``

2. You should now be able to do ``import webbpsf`` in a Python session. 

3. Future versions may be installed with ``pip install --upgrade webbpsf`` when they become available.

However, this installs only the program code. You still must download and install the data files, as :ref:`described below <data_install>`. 


Installing WebbPSF manually
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If for some reason you do not wish to use PYPI, you can just install the source file directly:


1. Download the following file: `webbpsf-0.2.9.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-0.2.9.tar.gz>`_
2. Untar ``webbpsf-0.2.9.tar.gz`` into a temporary working directory. 
3. Run ``python setup.py install`` in that directory. This will install ``webbpsf`` into your Python path. 

   If you lack the filesystem permissions to write into the system python directory 
   (for instance, on a machine you don't have root on), you can do ``python setup.py install --user`` to install locally
   in your home directory.
4. You should now be able to do ``import webbpsf`` in a Python session. 



Installing WebbPSF development versions, and/or contributing to its development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As of version 0.2.8, the `webbpsf source code repository <https://github.com/mperrin/webbpsf>`_ is hosted at GitHub, as is the repository for `poppy <https://github.com/mperrin/poppy>`_. Users may clone, fork, and pull diffs in the usual manner. Pull requests with code enhancements welcomed!  

.. _data_install:

Installing the Required Data Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Download the following file:  `webbpsf-data-0.2.6.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-data-0.2.6.tar.gz>`_  [417 MB]
2. Untar ``webbpsf-data-0.2.x.tar.gz`` into a directory of your choosing.
3. Set the environment variable ``WEBBPSF_PATH`` to point to that directory. e.g. ``setenv WEBBPSF_PATH $HOME/data/webbpsf-data`` for tcsh/csh, or ``WEBBPSF_PATH=$HOME/data/webbpsf-data; export WEBBPSF_PATH`` for bash.


.. _pysynphot_install:

Installing or updating pysynphot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pysynphot is an optional dependency, but is highly recommended. 

To install or update ``pysynphot``, do the following. (See also http://stsdas.stsci.edu/pysynphot/ and https://trac6.assembla.com/astrolib). If you already have ``pysynphot`` 
installed, it will probably work fine without this update, but computations may be slower if you have a version earlier than 0.8.  WebbPSF has most recently been tested using pysynphot v 0.8.3

.. comment 
        work without this update but computations will be slower than the current version, so we recommend updating it. 

1. Download the most recent version of pysynphot from https://trac6.assembla.com/astrolib. 
2. Untar that file into a temporary working directory. 
3. run ``python setup.py install`` in that directory.  You can delete the setup files there after you do this step. 
4. If this is your initial installation of ``pysynphot`` you need to install the CDBS files. See the `pysynphot installation guide <https://trac6.assembla.com/astrolib/wiki/PysynphotInstallationGuide>`_. The necessary files are available from https://trac6.assembla.com/astrolib; follow the download links for "throughput files" and "model spectra". If you already have CDBS installed, then you're all set and can skip this step.


WebbPSF includes its own normalized copies of the new JWST instrumental
throughputs from the development CDBS at STScI.  If you have JWST throughput
files available in your ``$PYSYN_CDBS`` directory (likely true only for
internal users at STScI), those will be used in preference to the WebbPSF
internal files, but this is not required.

.. comment
        3. Untar ``CDBS-for-webb.tar.gz`` in a directory of your choosing. (Typically replacing into your current CDBS directory if already present)
        4. Set the environment variable ``PYSYN_CDBS`` to point to that directory. e.g. ``setenv PYSYN_CDBS $HOME/data/CDBS``.


Note for STScI Internal Users
---------------------------------


Webbpsf is installed centrally on the WITServ computers for use by all members of the Webb instrument teams. 

The directory ``/witserv/data10/software`` contains shared software tools, currently a handful of Python modules, WebbPSF, and a copy of CDBS. 

The file ``/witserv/data10/software/README.txt`` gives a little bit of documentation, but briefly, it should be enough to add the line::

    source /witserv/data10/software/setup.tcsh

to your shell startup files on witserv* (assuming you're using tcsh), and then you should be able to run ``webbpsfgui`` from the command line, or start a python session and ``import webbpsf``.

Prerelease access to updated versions of the CDBS files may be available; contact Marshall if interested. 



