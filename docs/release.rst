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
 #. Update the link in ``installation.rst`` under :ref:`data_install`

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

Uploading updated docs
^^^^^^^^^^^^^^^^^^^^^^

 #. On the same `PyPI package editing page <https://pypi.python.org/pypi?%3Aaction=pkg_edit&name=webbpsf>`_, there is a field for uploading updated documentation, so keep that window open.
 #. In your terminal, ``cd`` to ``webbpsf/docs`` and run ``make html``
 #. ``cd _build/html``
 #. Build a zip file for upload: ``zip -r docs-to-upload.zip ./*``
 #. Choose the new zip file from the file field on the `PyPI package editing page <https://pypi.python.org/pypi?%3Aaction=pkg_edit&name=webbpsf>`_ and click "Upload Documentation"

Finishing the release
^^^^^^^^^^^^^^^^^^^^^

 #. Commit your edits to ``relnotes.rst`` and ``setup.py``
 #. Tag that commit as the release with ``git tag v0.X.Y`` and push the tags to origin and upstream with ``git push --tags origin`` and ``git push --tags upstream``
 #. Edit ``setup.py`` to increment the version number in the ``VERSION`` variable and re-add the ``.dev`` suffix
 #. Edit ``relnotes.rst`` to add a new heading for the upcoming version
 #. Commit your edits with a message like "Back to development: version 0.X.Y+1"