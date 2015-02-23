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

**TODO: explain how to use the scripts in dev_utils**

Releasing new versions on PyPI
==============================

 1. Edit ``relnotes.rst`` to add a release date and reference anchor (e.g. ``.. _rel0.X.Y:``) to the section for this release
 2. Add any important notes to the appropriate section in the release notes
 3. Edit ``setup.py`` in this repository to remove ``.dev`` from the version number in the ``VERSION`` variable
 4. Build a source distribution with ``python setup.py build sdist``
 5. Copy the resulting file (``webbpsf-0.X.Y.tar.gz``) to a new folder, extract it, and ``cd`` there
 6. Run ``python setup.py test`` (preferably in a new ``virtualenv`` containing only the WebbPSF dependencies) and verify that the test suite passes with the code you're about to release
 7. If that runs as expected, ``cd`` back to your ``webbpsf`` repository and run ``twine upload dist/webbpsf-0.X.Y.tar.gz`` for your new version
 8. Verify that the latest version is visible and others are hidden on the `PyPI package editing page <https://pypi.python.org/pypi?%3Aaction=pkg_edit&name=webbpsf>`_

Uploading updated docs
^^^^^^^^^^^^^^^^^^^^^^

 1. On the same `PyPI package editing page <https://pypi.python.org/pypi?%3Aaction=pkg_edit&name=webbpsf>`_, there is a field for uploading updated documentation, so keep that window open.
 2. In your terminal, ``cd`` to ``webbpsf/docs`` and run ``make html``
 3. ``cd _build/html``
 4. Build a zip file for upload: ``zip -r docs-to-upload.zip ./*``
 5. Choose the new zip file from the file field on the `PyPI package editing page <https://pypi.python.org/pypi?%3Aaction=pkg_edit&name=webbpsf>`_ and click "Upload Documentation"

Finishing the release
^^^^^^^^^^^^^^^^^^^^^

 1. Commit your edits to ``relnotes.rst`` and ``setup.py``
 2. Tag that commit as the release with ``git tag v0.X.Y`` and push the tags to origin and upstream with ``git push --tags origin`` and ``git push --tags upstream``
 3. Edit ``setup.py`` to increment the version number in the ``VERSION`` variable and re-add the ``.dev`` suffix
 4. Edit ``relnotes.rst`` to add a new heading for the upcoming version
 5. Commit your edits with a message like "Back to development: version 0.X.Y+1"