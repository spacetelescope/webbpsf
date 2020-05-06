.. _installation:

Requirements & Installation
===========================

.. note:: 

   This is entirely optional, but you may wish to sign up to the mailing list ``webbpsf-users@maillist.stsci.edu``. This is a very low-traffic moderated announce-only list, to which we will periodically post announcements of updates to this software.

   To subscribe, visit  the `maillist.stsci.edu server <https://maillist.stsci.edu/scripts/wa.exe?SUBED1=Webbpsf-users&A=1>`_


Recommended: Installing via AstroConda
--------------------------------------------------

For ease of installation, we recommend using `AstroConda <http://astroconda.readthedocs.io/en/latest/>`_, an astronomy-optimized software distribution for scientific Python built on Anaconda. Install AstroConda according to `their instructions <http://astroconda.readthedocs.io/en/latest/installation.html>`_, then activate the environment with::

   $ source activate astroconda

(Note: if you named your environment something other than ``astroconda``, change the above command appropriately.)

Next, install WebbPSF (along with all its dependencies and required reference data) with::

   (astroconda)$ conda install webbpsf

Updates to the latest version can be done as for any conda package::

   (astroconda)$ conda update webbpsf


.. _install-with-conda:

Installing with conda (but not AstroConda)
-------------------------------------------

If you already use ``conda``, but do not want to install the full suite of STScI software, you can simply add the AstroConda *channel* and install WebbPSF as follows (creating a new environment named ``webbpsf-env``)::

   $ conda config --add channels http://ssb.stsci.edu/astroconda
   $ conda create -n webbpsf-env webbpsf
   $ source activate webbpsf-env

Upgrading to the latest version is done with ``conda update -n webbpsf-env --all``.

.. warning::

   You *must* install WebbPSF into a specific environment (e.g. ``webbpsf-env``); our conda package will not work if installed into the default "root" environment.

.. _install_pip:

Installing with pip
-------------------

WebbPSF and its underlying optical library POPPY may be installed from the `Python Package Index <http://pypi.python.org/pypi>`_ in the usual manner for Python packages. ::

    $ pip install --upgrade webbpsf
    [... progress report ...]

    Successfully installed webbpsf

Note that ``pip install webbpsf`` only installs the program code. **If you install via pip, you must manually download and install the data files, as** :ref:`described <data_install>` **below.**
To obtain source spectra for calculations, you should also follow :ref:`installation instructions for pysynphot <pysynphot_install>`.


.. _pysynphot_install:

Installing or updating pysynphot
--------------------------------

Pysynphot is an optional dependency, but is highly recommended.  Pysynphot is best installed via AstroConda. Further installation instructions can be found in `the pysynphot docs <https://pysynphot.readthedocs.io/en/latest/#installation-and-setup>`_ or `a discussion in the POPPY docs <http://poppy-optics.readthedocs.io/en/stable/installation.html#installing-or-updating-pysynphot>`_.

.. _data_install:

Installing the Required Data Files
----------------------------------

*If you install via pip or manually*, you must install the data files yourself. If you install via Conda, the data files are automatically installed, in which case you can skip this section. 

Files containing such information as the JWST pupil shape, instrument throughputs, and aperture positions are distributed separately from WebbPSF. To run WebbPSF, you must download these files and tell WebbPSF where to find them using the ``WEBBPSF_PATH`` environment variable.

1. Download the following file:  `webbpsf-data-0.9.0.tar.gz <https://stsci.box.com/shared/static/qcptcokkbx7fgi3c00w2732yezkxzb99.gz>`_  [approx. 240 MB]
2. Untar ``webbpsf-data-0.9.0.tar.gz`` into a directory of your choosing.
3. Set the environment variable ``WEBBPSF_PATH`` to point to that directory. e.g. ::

   export WEBBPSF_PATH=$HOME/data/webbpsf-data

for bash. (You will probably want to add this to your ``.bashrc``.)

You should now be able to successfully ``import webbpsf`` in a Python session, or start the GUI with the command ``webbpsfgui``.

.. warning::

   If you have previously installed the data files for an earlier version of WebbPSF, and then update to a newer version, the
   software may prompt you that you must download and install a new updated version of the data files.

.. Note:: 

   **For STScI Users Only:** Users at STScI may access the required data files from the Central Storage network. Set the following environment variables in your ``bash`` shell. (You will probably want to add this to your ``.bashrc``.) ::

      export WEBBPSF_PATH="/grp/jwst/ote/webbpsf-data"
      export PYSYN_CDBS="/grp/hst/cdbs"

Software Requirements
---------------------


See `the environment.yml specification file <https://github.com/spacetelescope/webbpsf/blob/master/environment.yml>`_ for the required package dependencies. 

**Required Python version**: WebbPSF 0.8 and above require Python 3.5 or higher.

The major dependencies are the standard `NumPy, SciPy <http://www.scipy.org/scipylib/download.html>`_, `matplotlib <http://matplotlib.org>`_ stack, and `Astropy <http://astropy.org>`_

**Recommended Python packages**:

* `pysynphot <https://pypi.python.org/pypi/pysynphot>`_ enables the simulation
  of PSFs with proper spectral response to realistic source spectra.  Without
  this, PSF fidelity is reduced. See above for :ref:`installation instructions
  for pysynphot <pysynphot_install>`.  Pysynphot is recommended for most users.

**Optional Python packages**:

Some calculations with POPPY can benefit from the optional packages `psutil <https://pypi.python.org/pypi/psutil>`_ and `pyFFTW <https://pypi.python.org/pypi/pyFFTW>`_, but these are not needed in general. See `the POPPY installation docs <http://poppy-optics.readthedocs.io/en/stable/installation.html>`_ for more details.
These optional packages are only worth adding for speed improvements if you are spending substantial time running calculations. 

Additional packages are needed for the optional use of GPUs to accelerate calculations. See the POPPY documentation.

.. _install_dev_version:

Installing a pre-release version or contributing to WebbPSF development
-----------------------------------------------------------------------

The `WebbPSF source code repository <https://github.com/spacetelescope/webbpsf>`_ is hosted at GitHub, as is the repository for `POPPY <https://github.com/spacetelescope/poppy>`_. Users may clone or fork in the usual manner. Pull requests with code enhancements welcomed.

To install the current development version of WebbPSF, you can use ``pip`` to install directly from a ``git`` repository. To install WebbPSF and POPPY from ``git``, uninstall any existing copies of WebbPSF and POPPY, then invoke pip as follows::

    $ pip install -e git+https://github.com/spacetelescope/poppy.git#egg=poppy \
       -e git+https://github.com/spacetelescope/webbpsf.git#egg=webbpsf

This will create directories ``./src/poppy`` and ``./src/webbpsf`` in your current directory containing the cloned repository. If you have commit access to the repository, you may want to clone via ssh with a URL like ``git+ssh://git@github.com:spacetelescope/webbpsf.git``. Documentation of the available options for installing directly from Git can be found in the `pip documentation <http://pip.readthedocs.org/en/latest/reference/pip_install.html#git>`_.

Remember to :ref:`install the required data files <data_install>`, if you have not already installed them.
