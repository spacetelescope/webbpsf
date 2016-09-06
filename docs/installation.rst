.. _installation:

Requirements & Installation
===========================

The latest version of WebbPSF can be installed from either the `Python Package Index <https://pypi.python.org>`_ (PyPI) or the ``conda`` package management system.

.. admonition:: Optional: sign up to receive announcement of updates

   This is entirely optional, but you may wish to sign up to the mailing list ``webbpsf-users@stsci.edu``. This is a low-traffic moderated announce-only list, to which we will periodically post announcements of updates to this software.

   To subscribe, visit  the `maillist.stsci.edu server <https://maillist.stsci.edu/scripts/wa.exe?SUBED1=Webbpsf-users&A=1>`_

Installing with AstroConda
--------------------------

For ease of installation, we recommend using `AstroConda <http://astroconda.readthedocs.io/en/latest/>`_, which includes WebbPSF as well as NumPy, SciPy, and other packages that can be tricky to compile on your own.

Installing AstroConda according to `their instructions <http://astroconda.readthedocs.io/en/latest/installation.html>`_ will also install WebbPSF, its dependencies, and the reference data it needs.

Installing with conda (but no AstroConda)
-----------------------------------------

If you already use ``conda``, but do not want to install the full suite of STScI software, you can simply add the AstroConda *channel* and install as follows (for a conda environment named ``myenv``)::

   $ conda config --add channels http://ssb.stsci.edu/astroconda
   $ source activate myenv
   (myenv)$ conda install webbpsf

Upgrading to the latest version is done with ``conda update webbpsf``.

.. warning::

   You must install WebbPSF into a specific environment; our conda package will not work if installed into the "root" environment. Activate the environment where you want to install WebbPSF, and ``conda install webbpsf`` as above.

Installing with pip
-------------------

**If you have Python 2.7, 3.4, or 3.5 already installed another way**, WebbPSF and its underlying optical library POPPY may be installed from the `Python Package Index <http://pypi.python.org/pypi>`_ in the usual manner for Python packages. ::

    $ pip install webbpsf --upgrade
    [... progress report ...]

    Successfully installed webbpsf

However, ``pip install webbpsf`` only installs the program code. You still must download and install the data files, as :ref:`described below <data_install>`. To obtain source spectra for calculations, you should also follow :ref:`installation instructions for pysynphot <pysynphot_install>`.

If you wish to install webbpsf on a machine for which you do not have administrative access, you can do so by using Python's
built-in `"--user" mechanism  <http://docs.python.org/2/install/#alternate-installation-the-user-scheme>`_
for installing packages into your home directory. ::

    $ pip install webbpsf --user

.. warning::

   If you get the message ``SystemError: Cannot compile 'Python.h'. Perhaps you need to install python-dev|python-devel.`` during install *even when Python.h is available*, this means ``setup.py`` was unable to install NumPy. This can sometimes be fixed by executing ``pip install numpy`` separately, before installing WebbPSF. See the bug report at `numpy/numpy#2434 <https://github.com/numpy/numpy/issues/2434>`_ for details.

Software Requirements
---------------------

**Required Python version**: WebbPSF is supported on both Python 2.7 and 3.4+.

**Required Python packages**:

* Recent versions of `NumPy, SciPy <http://www.scipy.org/scipylib/download.html>`_ and `matplotlib <http://matplotlib.org>`_, if not installed already.
* `Astropy <http://astropy.org>`_, 1.0 or more recent.
* `POPPY <https://pypi.python.org/pypi/poppy>`_, 0.5.0 or more recent.

**Recommended Python packages**:

* `pysynphot <https://pypi.python.org/pypi/pysynphot>`_ enables the simulation
  of PSFs with proper spectral response to realistic source spectra.  Without
  this, PSF fidelity is reduced. See below for :ref:`installation instructions
  for pysynphot <pysynphot_install>`.  Pysynphot is recommended for most users.

**Optional Python packages**:

Some calculations with POPPY can benefit from the optional packages `psutil <https://pypi.python.org/pypi/psutil>`_ and `pyFFTW <https://pypi.python.org/pypi/pyFFTW>`_, but these are not needed in general. See `the POPPY installation docs <http://pythonhosted.org//poppy/installation.html>`_ for more details. 
These optional packages are only worth adding for speed improvements if you are spending substantial time running calculations.

.. _pysynphot_install:

Installing or updating pysynphot
--------------------------------

Pysynphot is an optional dependency, but is highly recommended.  Installation instructions can be found `here in the POPPY docs <http://pythonhosted.org//poppy/installation.html#installing-or-updating-pysynphot>`_.

.. _data_install:

Installing the Required Data Files
----------------------------------

Files containing such information as the JWST pupil shape, instrument throughputs, and aperture positions are distributed separately from WebbPSF. To run WebbPSF, you must download these files and tell WebbPSF where to find them using the ``WEBBPSF_PATH`` environment variable.

1. Download the following file:  `webbpsf-data-0.5.0.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-data-0.5.0.tar.gz>`_  [approx. 430 MB]
2. Untar ``webbpsf-data-0.5.0.tar.gz`` into a directory of your choosing.
3. Set the environment variable ``WEBBPSF_PATH`` to point to that directory. e.g. ::

       setenv WEBBPSF_PATH $HOME/data/webbpsf-data

   for tcsh/csh, or::

       WEBBPSF_PATH=$HOME/data/webbpsf-data; export WEBBPSF_PATH

   for bash. (You will probably want to add this to your ``.cshrc`` or ``.bashrc``.)

You should now be able to successfully ``import webbpsf`` in a Python session, or start the GUI with the command ``webbpsfgui``.

.. warning::

   If you have previously installed the data files for an earlier version of webbpsf, and then update to a newer version, the
   software may prompt you that you must download and install a new updated version of the data files.

.. admonition:: For STScI Users Only

   Users at STScI may access the required data files from the Central Storage network. If you use bash, export these environment variables::

      export WEBBPSF_PATH="/grp/jwst/ote/webbpsf-data"
      export PYSYN_CDBS="/grp/hst/cdbs"

   If you are using tcsh::

      setenv WEBBPSF_PATH /grp/jwst/ote/webbpsf-data
      setenv PYSYN_CDBS /grp/hst/cdbs

.. _install_dev_version:

Installing a pre-release version or contributing to WebbPSF development
-----------------------------------------------------------------------

The `WebbPSF source code repository <https://github.com/mperrin/webbpsf>`_ is hosted at GitHub, as is the repository for `POPPY <https://github.com/mperrin/poppy>`_. Users may clone or fork in the usual manner. Pull requests with code enhancements welcomed.

To install the current development version of WebbPSF, you can use ``pip`` to install directly from a ``git`` repository. To install WebbPSF and POPPY from ``git``, uninstall any existing copies of WebbPSF and POPPY, then invoke pip as follows::

    $ pip install -e git+https://github.com/mperrin/poppy.git#egg=poppy \
       -e git+https://github.com/mperrin/webbpsf.git#egg=webbpsf

This will create directories ``./src/poppy`` and ``./src/webbpsf`` in your current directory containing the cloned repository. If you have commit access to the repository, you may want to clone via ssh with a URL like ``git+ssh://git@github.com:mperrin/webbpsf.git``. Documentation of the available options for installing directly from Git can be found in the `pip documentation <http://pip.readthedocs.org/en/latest/reference/pip_install.html#git>`_.

Remember to :ref:`install the required data files <data_install>`, if you have not already installed them.
