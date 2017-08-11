**********************************
Releasing a new version of WebbPSF
**********************************

Prerequisites
=============

 * Is the build `passing on Travis? <https://travis-ci.org/mperrin/webbpsf>`_
 * Are you up to date with ``master`` on the upstream branch (mperrin/webbpsf)?
 * Do you have `twine <https://pypi.python.org/pypi/twine>`_ installed?
 * Do you have access to WebbPSF on PyPI with the owner or maintainer role?
 * Do you have your ``~/.pypirc`` filled out? (`more info <https://python-packaging-user-guide.readthedocs.org/en/latest/distributing.html#register-your-project>`_)

Releasing new data packages
===========================

 #. Run ``dev_utils/make-data-sdist.sh`` (details below) to make a gzipped tarred archive of the WebbPSF data
 #. If the new data package is **required** (meaning you can't run WebbPSF without it, or you can run but may get incorrect results), you must bump ``DATA_VERSION_MIN`` in ``__init__.py`` to ``(0, X, Y)``
 #. Extract the resulting data archive and check that you can run the WebbPSF tests with ``WEBBPSF_PATH`` pointing to it
 #. Copy the data archive into public web space
 #. ``cd`` to ``/grp/jwst/ote`` and remove the ``webbpsf-data`` symlink
 #. Copy the archive into ``/grp/jwst/ote/`` and extract it to ``/grp/jwst/ote/webbpsf-data``
 #. Rename the folder to ``webbpsf-data-0.x.y``
 #. Create a symbolic link at ``/grp/jwst/ote/webbpsf-data`` to point to the new folder
 #. Update the URL in ``installation.rst`` under :ref:`data_install`

Invoke ``dev_utils/make-data-sdist.sh`` one of the following ways to make a gzipped tarred archive of the WebbPSF data suitable for distribution.

**If you are on the Institute network:** ::

   $ cd webbpsf/dev_utils/
   $ ./make-data-sdist.sh 0.X.Y
   $ cp ./webbpsf-data-0.X.Y.tar.gz /path/to/public/web/directory/

**If you're working from a local data root:** ::

   $ cd webbpsf/dev_utils/
   $ DATAROOT="/Users/you/webbpsf-data-sources/" ./make-data-sdist.sh 0.X.Y
   $ cp ./webbpsf-data-0.X.Y.tar.gz /where/ever/you/want/

Releasing new versions on PyPI
==============================

 #. Edit ``relnotes.rst`` to add a release date and reference anchor (e.g. ``.. _rel0.X.Y:``) to the section for this release
 #. Update the link to "What's new" in ``index.rst``
 #. Add any important notes to the appropriate section in the release notes
 #. Edit ``setup.py`` in this repository to remove ``.dev`` from the version number in the ``VERSION`` variable
 #. Build a source distribution with ``python setup.py build sdist``
 #. Copy the resulting file (``webbpsf-0.X.Y.tar.gz``) to a new folder, extract it, and ``cd`` there
 #. Run ``python setup.py test`` (preferably in a new ``virtualenv`` containing only the WebbPSF dependencies) and verify that the test suite passes with the code you're about to release
 #. If that runs as expected, ``cd`` back to your ``webbpsf`` repository and run ``twine upload dist/webbpsf-0.X.Y.tar.gz`` for your new version
 #. Verify that the latest version is visible and others are hidden on the `PyPI package editing page <https://pypi.python.org/pypi?%3Aaction=pkg_edit&name=webbpsf>`_

Finishing the release
^^^^^^^^^^^^^^^^^^^^^

 #. Commit your edits to ``relnotes.rst`` and ``setup.py``
 #. Tag that commit as the release with ``git tag v0.X.Y`` and push the tags to origin and upstream with ``git push --tags origin`` and ``git push --tags upstream``
 #. Edit ``setup.py`` to increment the version number in the ``VERSION`` variable and re-add the ``.dev`` suffix
 #. Edit ``relnotes.rst`` to add a new heading for the upcoming version
 #. Commit your edits with a message like "Back to development: version 0.X.Y+1"
 #. Email an announcement to ``webbpsf-users@stsci.edu``

Releasing a new version through AstroConda
==========================================

To test that an Astroconda package builds, you will need ``conda-build``::

   $ conda install conda-build

#. Fork (if needed) and clone https://github.com/astroconda/astroconda-contrib
#. If there is a new version of POPPY available to package, edit `poppy/meta.yaml <https://github.com/astroconda/astroconda-contrib/blob/master/poppy/meta.yaml>`_ to reflect the new ``version`` and ``git_tag``.
#. If the minimum needed version of the webbpsf-data package has changed in ``webbpsf/__init__.py``, edit `webbpsf-data/meta.yaml <https://github.com/astroconda/astroconda-contrib/blob/master/webbpsf-data/meta.yaml>`_ to reflect the new ``version`` and ``url``.
#. Edit `webbpsf/meta.yaml <https://github.com/astroconda/astroconda-contrib/blob/master/webbpsf/meta.yaml>`_ to reflect the new versions of POPPY and webbpsf-data, if necessary.
#. Edit in the ``git_tag`` name from ``git tag`` in the PyPI release instructions (``v0.X.Y``).
#. Verify that you can build the package from the astroconda-contrib directory: ``conda build -c http://ssb.stsci.edu/astroconda webbpsf``
#. Commit your changes to a new branch and push to GitHub.
#. Create a pull request against ``astroconda/astroconda-contrib``.
#. Wait for SSB to build the conda packages.
#. (optional) Create a new conda environment to test the package installation following :ref:`these instructions <install-with-conda>`.
