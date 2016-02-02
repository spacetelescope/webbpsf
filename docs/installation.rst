.. _installation:

Requirements & Installation
============================

WebbPSF uses the `Python Package Index <https://pypi.python.org>`_ (PyPI) to distribute new versions. For ease of installation, we recommend a scientific Python distribution like `Ureka <http://ssb.stsci.edu/ureka/>`_ that includes NumPy, SciPy, and other packages that can be tricky to compile on your own.

If you have Python 2.7 or 3.4, pip, and NumPy already installed, you can easily install or upgrade WebbPSF with::

    $ pip install -U webbpsf

Once the WebbPSF code has been installed, you should proceed to :ref:`install the required data files <data_install>`. You may also wish to :ref:`set up Pysynphot <pysynphot_install>`, a recommended optional dependency that improves PSF fidelity.

For detailed instructions and software requirements, read on.

.. admonition:: Optional: sign up to receive announcement of updates

   This is entirely optional, but you may wish to sign up to the mailing list ``webbpsf-users@stsci.edu``. This is a low-traffic moderated announce-only list, to which we will periodically post announcements of updates to this software.

   To subscribe, visit  the `maillist.stsci.edu server <https://maillist.stsci.edu/scripts/wa.exe?SUBED1=Webbpsf-users&A=1>`_

Software Requirements
-----------------------

**Required Python version**: WebbPSF is supported on both Python 2.7 and 3.4.

**Required Python packages**:

* Recent versions of `NumPy, SciPy <http://www.scipy.org/scipylib/download.html>`_ and `matplotlib <http://matplotlib.org>`_, if not installed already.
* `Astropy <http://astropy.org>`_, 1.0 or more recent.
* `POPPY <https://pypi.python.org/pypi/poppy>`_, 0.4.1 or more recent.

**Recommended Python packages**:

* `pysynphot <https://pypi.python.org/pypi/pysynphot>`_ enables the simulation
  of PSFs with proper spectral response to realistic source spectra.  Without
  this, PSF fidelity is reduced. See below for :ref:`installation instructions
  for pysynphot <pysynphot_install>`.  Pysynphot is recommended for most users. 


**Optional Python packages**:

Some calculations with POPPY can benefit from the optional packages `psutil <https://pypi.python.org/pypi/psutil>`_ and `pyFFTW <https://pypi.python.org/pypi/pyFFTW>`_, but these are not needed in general. See `the POPPY installation docs <http://pythonhosted.org//poppy/installation.html>`_ for more details. 
These optional packages are only worth adding for speed improvements if you are spending substantial time running calculations.


Installing WebbPSF
----------------------

WebbPSF and its underlying optical library POPPY may be installed from the
`Python Package Index <http://pypi.python.org/pypi>`_ in the usual manner for
Python packages. :: 

    $ pip install webbpsf --upgrade
    [... progress report ...]

    Successfully installed webbpsf

However, ``pip install webbpsf`` only installs the program code. You still must download and install the data files, as :ref:`described below <data_install>`. To obtain source spectra for calculations, you should also follow :ref:`installation instructions for pysynphot <pysynphot_install>`.

If you wish to install webbpsf on a machine for which you do not have administrative access, you can do so by using Python's
built-in `"--user" mechanism  <http://docs.python.org/2/install/#alternate-installation-the-user-scheme>`_
for installing packages into your home directory. ::

    $ pip install webbpsf --user

.. admonition:: For STScI Users Only

   Users at STScI may also access WebbPSF through the standard `SSB software distributions <http://ssb.stsci.edu/ssb_software.shtml>`_. Set the environment variable ``WEBBPSF_PATH`` to ``/grp/jwst/ote/webbpsf-data``, select the desired version of SSB (``ssbx``, ``ssbdev``, or ``ssbrel``), and the webbpsf package will be available.

The source code is hosted in `this repository on GitHub <https://github.com/mperrin/webbpsf>`_. See :ref:`below <install_dev_version>` for more on how
to install from the development source.


.. warning::
  If you get the message ``SystemError: Cannot compile 'Python.h'. Perhaps you need to install python-dev|python-devel.`` during install *even when Python.h is available*, this means ``setup.py`` was unable to install NumPy. This can sometimes be fixed by executing ``pip install numpy`` separately, before installing webbpsf. See the bug report at `numpy/numpy#2434 <https://github.com/numpy/numpy/issues/2434>`_ for details.

.. _pysynphot_install:

Installing or updating pysynphot
---------------------------------

Pysynphot is an optional dependency, but is highly recommended.  Installation instructions can be found `here in the POPPY docs <http://pythonhosted.org//poppy/installation.html#installing-or-updating-pysynphot>`_.

.. _data_install:

Installing the Required Data Files
---------------------------------------------

Files containing such information as the JWST pupil shape, instrument throughputs, and aperture positions are distributed separately from WebbPSF. To run WebbPSF, you must download these files and tell WebbPSF where to find them using the ``WEBBPSF_PATH`` environment variable.

1. Download the following file:  `webbpsf-data-0.4.0.tar.gz <http://www.stsci.edu/~mperrin/software/webbpsf/webbpsf-data-0.4.0.tar.gz>`_  [approx. 430 MB]
2. Untar ``webbpsf-data-0.4.0.tar.gz`` into a directory of your choosing.
3. Set the environment variable ``WEBBPSF_PATH`` to point to that directory. e.g. ::

    setenv WEBBPSF_PATH $HOME/data/webbpsf-data

   for tcsh/csh, or ::

    WEBBPSF_PATH=$HOME/data/webbpsf-data; export WEBBPSF_PATH

   for bash. (You will probably want to add this to your ``.cshrc`` or ``.bashrc``.)

You should now be able to successfully ``import webbpsf`` in a Python session, or start the GUI with the command ``webbpsfgui``.


.. warning::
  If you have previously installed the data files for an earlier version of webbpsf, and then update to a newer version, the
  software may prompt you that you must download and install a new updated version of the data files. 

.. admonition:: For STScI Users Only

  Users at STScI may access the required data files from the Central Storage network. 

    1. ``setenv WEBBPSF_PATH /grp/jwst/ote/webbpsf-data``  
    2. ``setenv PYSYN_CDBS /grp/hst/cdbs`` 

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

.. _install_dev_version:

Installing a pre-release version or contributing to WebbPSF development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `WebbPSF source code repository <https://github.com/mperrin/webbpsf>`_ is hosted at GitHub, as is the repository for `POPPY <https://github.com/mperrin/poppy>`_. Users may clone or fork in the usual manner. Pull requests with code enhancements welcomed.

To install the current development version of WebbPSF, you can use ``pip`` to install directly from a ``git`` repository. To install WebbPSF and POPPY from ``git``, uninstall any existing copies of WebbPSF and POPPY, then invoke pip as follows::

    $ pip install -e git+https://github.com/mperrin/poppy.git#egg=poppy \
       -e git+https://github.com/mperrin/webbpsf.git#egg=webbpsf

This will create directories ``./src/poppy`` and ``./src/webbpsf`` in your current directory containing the cloned repository. If you have commit access to the repository, you may want to clone via ssh with a URL like ``git+ssh://git@github.com:mperrin/webbpsf.git``. Documentation of the available options for installing directly from Git can be found in the `pip documentation <http://pip.readthedocs.org/en/latest/reference/pip_install.html#git>`_.

Remember to :ref:`install the required data files <data_install>`, if you have not already installed them.

