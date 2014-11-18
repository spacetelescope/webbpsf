Requirements & Installation
============================

WebbPSF uses the Python Package Index (PyPI) to distribute new versions. If you have Python 2.7 and ``pip`` installed, you can easily install or upgrade to the latest stable version of WebbPSF with::

    $ pip install numpy && pip install -U webbpsf

*If you do not have access to install packages system-wide, replace* ``pip install`` *with* ``pip install --user``.

Once the WebbPSF code has been installed, you can then proceed to :ref:`installing the required data files <data_install>`. (Why install NumPy first? To work around a bug! See `numpy/numpy#2434 <https://github.com/numpy/numpy/issues/2434>`_ for details.)

For ease of installation, we recommend a scientific Python distribution like `Ureka <http://ssb.stsci.edu/ureka/>`_. Ureka includes NumPy, SciPy, matplotlib, and other packages that can be tricky to compile on your own machine.

.. admonition:: Optional: sign up to receive announcement of updates

    This is entirely optional, but you may wish to sign up to the mailing list ``webbpsf-users@stsci.edu``. This is a low-traffic moderated announce-only list, to which we will periodically post announcements of updates to this software.

    To subscribe, email majordomo@stsci.edu with the message body text ``subscribe webbpsf-users``.

Software Requirements
-----------------------

**Required Python version**: Python 2.6 or higher is required. Python 2.7 is required for the GUI (see below) and strongly recommended overall. WebbPSF is not yet Python 3 compatible.

**Required Python packages**:

* Recent versions of `NumPy, SciPy <http://www.scipy.org/scipylib/download.html>`_ and `matplotlib <http://matplotlib.org>`_, if not installed already.
* `Astropy <http://astropy.org>`_, 0.4 or more recent.
* `POPPY <https://pypi.python.org/pypi/poppy>`_, 0.3.1 or more recent.

**Recommended Python packages**:

* `pysynphot <https://pypi.python.org/pypi/pysynphot>`_ enables the simulation of PSFs with proper spectral response to realistic source spectra.  Without this, PSF fidelity is reduced. See below for :ref:`installation instructions for pysynphot <pysynphot_install>`.
* `ttk <http://docs.python.org/2/library/ttk.html>`_ (included with Python 2.7) or `wxPython <http://www.wxpython.org>`_ is required for the :ref:`graphical user interface<gui>`.

Pysynphot is recommended for most users. If you are only using WebbPSF through the API, you do not need a GUI toolkit, but ttk should be present by default on Python 2.7 and newer. The optional packages below are only worth adding for speed improvements if you are spending substantial time running calculations.

**Optional Python packages**:

* `psutil <https://pypi.python.org/pypi/psutil>`_ enables slightly better automatic selection of numbers of processes for multiprocess calculations.
* `pyFFTW <https://pypi.python.org/pypi/pyFFTW>`_. The FFTW library can speed up the FFTs used in coronagraphic simulations and slit spectroscopy. Since direct imaging simulations use a discrete matrix FFT instead, direct imaging simulation speed is unchanged.  pyFFTW is recommended if you expect to perform many coronagraphic calculations, particularly for MIRI.

(Note: WebbPSF previously made use of the PyFFTW3 package, which is *different* from pyFFTW. The latter is more actively maintained and supported today, hence the switch.) See the :ref:`performance_and_parallelization` page for more details.

Installing WebbPSF
----------------------

The latest stable version of WebbPSF and its underlying optical library POPPY are both installable from the `Python Package Index <http://pypi.python.org/pypi>`_ via the standard toolchain using `pip <https://pip.pypa.io/>`_.  This is the easiest installation method if you already have a working copy of Python, NumPy, and matplotlib on your computer. (Alternatively, see the section on :ref:`alternate_install`.)

Simply invoke pip in the usual manner::

    $ pip install webbpsf
    [... progress report ...]

    Successfully installed webbpsf

However, ``pip install webbpsf`` only installs the program code. You still must download and install the data files, as :ref:`described below <data_install>`. To obtain source spectra for calculations, you should also follow :ref:`installation instructions for pysynphot <pysynphot_install>`.

Future versions may be installed with ``pip install --upgrade webbpsf`` when they become available.

.. note::
  If you wish to install webbpsf on a machine for which you do not have administrative access, you can do so by using Python's
  built-in `"--user" mechanism  <http://docs.python.org/2/install/#alternate-installation-the-user-scheme>`_
  for installing packages into your home directory. ::

    $ pip install webbpsf --user

.. warning::
  If you get the message ``SystemError: Cannot compile 'Python.h'. Perhaps you need to install python-dev|python-devel.`` during install *even when Python.h is available*, this means ``setup.py`` was unable to install NumPy. This can sometimes be fixed by executing ``pip install numpy`` separately, before installing webbpsf. See the bug report at `numpy/numpy#2434 <https://github.com/numpy/numpy/issues/2434>`_ for details.

.. _pysynphot_install:

Installing or updating pysynphot
---------------------------------

Pysynphot is an optional dependency, but is highly recommended. 

To install or update ``pysynphot``, simply invoke ``pip install -U pysynphot``. WebbPSF has most recently been tested using pysynphot 0.9.5 but is known to work well with earlier versions as well.

If you already have the CDBS data package installed, or are using WebbPSF at STScI, then you're all set and can skip the rest of this section.

If this is your initial installation of ``pysynphot``, you will need to install the CDBS files. These are available from STScI in DMG form for Mac users, as well as in gzipped tar format.

**Installing CDBS on Mac:** To obtain the DMG, consult the "Installing CDBS locally on a Mac" section of http://ssb.stsci.edu/ssb_software.shtml. Download the DMG and open it to find ``cdbs.pkg``. Running this graphical installer will place the CDBS files in ``/usr/stsci/stdata``. Set the environment variable ``PYSYN_CDBS`` to point to that directory, e.g. ``setenv PYSYN_CDBS /usr/stsci/stdata`` for tcsh/csh or ``export PYSYN_CDBS="/usr/stsci/stdata"`` for bash.

**Installing CDBS from tar archives**: To obtain the tar files, consult http://www.stsci.edu/hst/observatory/crds/cdbs_throughput.html. Download the archives numbered ``synphot[1-6].tar.gz`` and extract them to a directory such as ``$HOME/data/CDBS``.
Set the environment variable ``PYSYN_CDBS`` to point to that directory. e.g. ``setenv PYSYN_CDBS $HOME/data/CDBS`` for tcsh/csh or ``export PYSYN_CDBS="$HOME/data/CDBS"``.

WebbPSF includes its own normalized copies of the new JWST instrumental throughputs from the development CDBS at STScI.  If you have JWST throughput files available in your ``$PYSYN_CDBS`` directory (likely true only for internal users at STScI), those will be used in preference to the WebbPSF internal files, but this is not required.

.. _data_install:

Installing the Required Data Files
---------------------------------------------

Files containing such information as the JWST pupil shape, instrument throughputs, and aperture positions are distributed separately from WebbPSF. To run WebbPSF, you must download these files and tell WebbPSF where to find them using the ``WEBBPSF_PATH`` environment variable.

1. Download the following file:  `webbpsf-data-0.3.0.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-data-0.3.0.tar.gz>`_  [417 MB]
2. Untar ``webbpsf-data-0.3.0.tar.gz`` into a directory of your choosing.
3. Set the environment variable ``WEBBPSF_PATH`` to point to that directory. e.g. ::

    setenv WEBBPSF_PATH $HOME/data/webbpsf-data

   for tcsh/csh, or ::

    WEBBPSF_PATH=$HOME/data/webbpsf-data; export WEBBPSF_PATH

   for bash. (You will probably want to add this to your ``.cshrc`` or ``.bashrc``.)

You should now be able to successfully ``import webbpsf`` in a Python session, or start the GUI with the command ``webbpsfgui``.

.. admonition:: For STScI Users Only

  Users at STScI may access WebbPSF through the standard `SSB software distributions <http://ssb.stsci.edu/ssb_software.shtml>`_. 
  In particular, webbpsf and its required dependencies are now included in SSBDEV and will soon be in SSBX.  To make use of this,
  it should be sufficient to:

    1. Install SSBDEV and select that version of Python (e.g. ``ur_setup common ssbdev``)
    2. ``setenv WEBBPSF_PATH /grp/jwst/ote/webbpsf-data``  
    3. ``setenv PYSYN_CDBS /grp/hst/cdbs`` 

.. _alternate_install:

Alternate Installation Methods
---------------------------------------

Installing with `conda <http://conda.pydata.org>`_ or `miniconda <http://conda.pydata.org/miniconda.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Many users have expressed a preference for the `Anaconda <https://store.continuum.io/cshop/anaconda/>`_ distribution of scientific computing tools. Fortunately, it is straightforward to install WebbPSF into a ``conda`` environment.

1. Create a new environment for WebbPSF to live in::

    conda create -n webbpsf numpy scipy matplotlib pip

2. Activate the environment so that the next command takes effect in the new environment::

    source activate webbpsf

3. Install WebbPSF with pip::

    pip install webbpsf

You must next download and install the data files, as described in :ref:`data_install`. To obtain source spectra for calculations, you should also follow :ref:`installation instructions for pysynphot <pysynphot_install>`.

Later, when you open a new terminal window, remember to run ``source activate webbpsf`` before running ``webbpsfgui`` or attempting to ``import webbpsf``. You may also install webbpsf in the default environment, if that is more convenient for you. Simply ensure the packages listed in step 1 are installed with ``conda install``, then ``pip install webbpsf``.

Installing a pre-release version or contributing to WebbPSF development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `WebbPSF source code repository <https://github.com/mperrin/webbpsf>`_ is hosted at GitHub, as is the repository for `POPPY <https://github.com/mperrin/poppy>`_. Users may clone or fork in the usual manner. Pull requests with code enhancements welcomed.

To install the current development version of WebbPSF, you can use ``pip`` to install directly from a ``git`` repository. To install WebbPSF and POPPY from ``git``, uninstall any existing copies of WebbPSF and POPPY, then invoke pip as follows::

    $ pip install -e git+https://github.com/mperrin/poppy.git#egg=poppy \
       -e git+https://github.com/mperrin/webbpsf.git#egg=webbpsf

This will create directories ``./src/poppy`` and ``./src/webbpsf`` in your current directory containing the cloned repository. If you have commit access to the repository, you may want to clone via ssh with a URL like ``git+ssh://git@github.com:mperrin/webbpsf.git``. Documentation of the available options for installing directly from Git can be found in the `pip documentation <http://pip.readthedocs.org/en/latest/reference/pip_install.html#git>`_.

Remember to :ref:`install the required data files <data_install>`, if you have not already installed them.

Installing WebbPSF manually
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If for some reason you don't wish to use PyPI, you can just install from the source directly:

1. Download the following files.

 * http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-0.3rc2.tar.gz
 * http://www.stsci.edu/~mperrin/software/webbpsf/poppy-0.3rc2.tar.gz

2. Untar each into a temporary working directory.
3. Run ``python setup.py install`` in each of those directories to install first ``poppy`` and then ``webbpsf``.

You should now be able to do ``import webbpsf`` in a Python session to start WebbPSF.

However, the above installs only the program code. You still must download and install the data files, as :ref:`described below <data_install>`.

.. note::
   If you lack the filesystem permissions to write into the system Python directory (for instance, on a machine you don't have root on), you can do ``python setup.py install --user`` to install locally in your home directory.
