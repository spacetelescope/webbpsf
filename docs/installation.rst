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

* `astropy <http://astropy.org>`_, 0.2 or more recent.
* `lxml <https://pypi.python.org/pypi/lxml>`_ for XML parsing

  
These are optional but recommended:

* `pysynphot <https://trac6.assembla.com/astrolib>`_ enables the simulation of PSFs with proper spectral response to realistic source spectra.  Without this, PSF fidelity is reduced. See below for :ref:`installation instructions for pysynphot <pysynphot_install>`. 
* `psutil <https://pypi.python.org/pypi/psutil>`_ for cross-platform memory management enables better automatic selection of numbers of processes for multiprocess calculations.
* `pyFFTW <https://pypi.python.org/pypi/pyFFTW>`_. The FFTW library can speed up the FFTs used in coronagraphic simulations and slit spectroscopy. Since direct imaging simulations use a discrete matrix FFT instead, direct imaging simulation speed is unchanged.  pyFFTW is recommended if you expect to perform many coronagraphic calculations, particularly for MIRI.  (Note: WebbPSF previously made use of the PyFFTW3 package, which is *different* from pyFFTW. The latter is more actively maintained and supported today, hence the switch.) See the :ref:`performance_and_parallelization` page for more details.

Additional requirement for the GUI: The :ref:`graphical user interface<gui>` requires 

* **Either**  the `wxpython <http://www.wxpython.org>`_ interface to the ``wxwidgets`` widget library (recommended), 


* **or**  the `ttk <http://docs.python.org/2/library/ttk.html>`_ enhanced version of the ``Tkinter`` widget library. 
  (``ttk`` is not included by default on some installations of Python, for instance Mac OS 10.6's default of Python 2.6. 
  You may wish to either upgrade to a more current Python, or else compile and install ``ttk`` for your platform. ``WebbPSF``
  was developed using Python 2.7, which includes ``ttk`` by default, but it ought to work fine on any installations of
  Python 2.5 or 2.6 provided ``ttk`` is available.)

Similar GUIs are implemented in both widget tool kits, with the wxpython GUI
providing additional functionality. Future development will 
concentrate on the wxpython toolkit only, but for now both are supported.

Alternatively, you can just skip using the GUI; the optical modeling classes
themselves have no dependency on these widgets.


Installing WebbPSF
----------------------

.. warning::
  For a prerelease version, you can't install off of PyPI yet. Please use the download links below.

WebbPSF and its underlying optical library ``poppy`` are both
installable from the `Python Package Index <http://pypi.python.org/pypi>`_ via
the standard toolchain using `pip
<http://www.pip-installer.org/en/latest/index.html>`_.  This is the easiest installation
method if you already have a working copy of python, numpy, and matplotlib on your computer. 


Simply invoke pip in the usual manner::

   $ pip install webbpsf
   [... progress report ...]

   ``Successfully installed webbpsf``

You should now be able to do ``import webbpsf`` in a Python session to start WebbPSF. 

However, the above installs only the program code. You still must download and install the data files, as :ref:`described below <data_install>`. 

Future versions may be installed with ``pip install --upgrade webbpsf`` when they become available.

.. note::
  If you wish to install webbpsf on a machine for which you do not have administrative access, you can do so by using Python's
  built-in `"--user" mechanism  <http://docs.python.org/2/install/#alternate-installation-the-user-scheme>`_
  for installing packages into your home directory. ::

    $ pip install webbpsf --user


Installing WebbPSF manually
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If for some reason you don't wish to use PYPI, you can just install from the source directly:

1. Download the following files.

 * http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-0.3rc2.tar.gz
 * http://www.stsci.edu/~mperrin/software/webbpsf/poppy-0.3rc2.tar.gz

.. comment 
	#.. pypi-release:: webbpsf
	#   :prefix: Download
	#   :class: note
	#
	#
	#.. pypi-release:: poppy
	#   :prefix: Download
	#   :class: note


2. Untar each into a temporary working directory. 
3. Run ``python setup.py install`` in each of those directories to install first ``poppy`` and then ``webbpsf``. 


You should now be able to do ``import webbpsf`` in a Python session to start WebbPSF. 

However, the above installs only the program code. You still must download and install the data files, as :ref:`described below <data_install>`. 


.. note::
   If you lack the filesystem permissions to write into the system python directory 
   (for instance, on a machine you don't have root on), you can do ``python setup.py install --user`` to install locally
   in your home directory.


Installing WebbPSF development versions, and/or contributing to its development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The `webbpsf source code repository <https://github.com/mperrin/webbpsf>`_ is hosted at GitHub, as is the repository for `poppy <https://github.com/mperrin/poppy>`_. Users may clone or fork in the usual manner. Pull requests with code enhancements welcomed.  

.. _data_install:

Installing the Required Data Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Download the following file:  `webbpsf-data-0.3.0.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-data-0.3.0.tar.gz>`_  [417 MB]
2. Untar ``webbpsf-data-0.3.0.tar.gz`` into a directory of your choosing.
3. Set the environment variable ``WEBBPSF_PATH`` to point to that directory. e.g. ``setenv WEBBPSF_PATH $HOME/data/webbpsf-data`` for tcsh/csh, or ``WEBBPSF_PATH=$HOME/data/webbpsf-data; export WEBBPSF_PATH`` for bash.




Optional: sign up to receive announcement of updates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. comment 
    The first time you import WebbPSF, it will ask you whether you want to sign up for announcements of new versions::
            >>> import webbpsf
            *********************************************
            *           WebbPSF Initialization          *
            *********************************************
            This appears to be the first time you have used WebbPSF.
            Would you like to register your email address to
            stay informed of future versions, updates, etc?
            This will also register some basic information about
            your system (OS, Python version, WebbPSF version, etc.)
            to help us better support this software.
            Register? [Y/n]


This is entirely optional, but you may wish to sign up to the mailing list ``webbpsf-users@stsci.edu``.
This is a low-traffic moderated announce-only list, to which we will periodically post announcements of updates to this software.  

To subscribe, email `majordomo@stsci.edu` with the message body text ``"subscribe webbpsf-users"``. 


.. _pysynphot_install:

Installing or updating pysynphot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pysynphot is an optional dependency, but is highly recommended. 

To install or update ``pysynphot``, do the following. (See also http://stsdas.stsci.edu/pysynphot/ and https://trac6.assembla.com/astrolib). WebbPSF has most recently been tested using pysynphot 0.9.5 but is known to work well with earlier versions as well.


.. warning::
   You may have trouble installing pysynphot, as the zip file of the source on pypi is broken. This has been
   communicated upstream but not yet fixed. You may have more luck installing from an updated zip file 
   on testpypi: https://testpypi.python.org/pypi/pysynphot/0.9.5
   To install this, use this command::

     pip install -i https://testpypi.python.org/pypi pysynphot

.. comment 
        work without this update but computations will be slower than the current version, so we recommend updating it. 
    1. Download the most recent version of pysynphot from https://trac6.assembla.com/astrolib. 
    2. Untar that file into a temporary working directory. 
    3. run ``python setup.py install`` in that directory.  You can delete the setup files there after you do this step. 

If this is your initial installation of ``pysynphot`` you need to install the CDBS files. See the `pysynphot installation guide <https://trac6.assembla.com/astrolib/wiki/PysynphotInstallationGuide>`_. The necessary files are available from https://trac6.assembla.com/astrolib; follow the download links for "throughput files" and "model spectra". If you already have CDBS installed, then you're all set and can skip this step.


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



