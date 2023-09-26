***************************************************
Developer Notes: Releasing a new version of WebbPSF
***************************************************

Prerequisites
=============

 * Is the `develop` build `passing on Github Actions? <https://github.com/spacetelescope/webbpsf/actions>`_ with all desired release items included?

Releasing new data packages
===========================

 #. Run ``dev_utils/master_data_release.sh`` (details below) to make a gzipped tarred archive of the WebbPSF data
 #. If the new data package is **required** (meaning you can't run WebbPSF without it, or you can run but may get incorrect results), you must bump ``DATA_VERSION_MIN`` in ``__init__.py`` to ``(0, X, Y)``
 #. Extract the resulting data archive and check that you can run the WebbPSF tests with ``WEBBPSF_PATH`` pointing to it
 #. Copy the data archive into public web space. This now means on Box. The following steps need to be performed in this sequence in order to preserve the naming conventions.
     #. Find webbpsf-data-LATEST.tar.gz, and click on "more options" and "Update Version".  Choose the newest version of webbpsf-data-#.#.#.tar.gz
     #. This will change the name of webbpsf-data-LATEST.tar.gz to be what you just uploaded, rename the file back to "webbpsf-data-LATEST.tar.gz"
     #. Upload to Box a separate version of webbpsf-data-#.#.#.tar.gz shared data folder for future storage.
     #. Upload to Box the minimal-webbpsf-data-#.#.#.tar.gz shared data folder.
     #. Verify the shared link of webbpsf-data-latest.tar.gz is the same that exists in ``docs/installation.rst`` ("copy shared link" then "link settings")

 #. A shared copy will be automatically configured in STScI Central Store with updated symlink ``/grp/jwst/ote/webbpsf-data``
 #. Update the URL in ``installation.rst`` under :ref:`data_install`

Details for using `master_data_release.sh`:
-------------------------------------

Invoke ``dev_utils/master_data_release.sh`` one of the following ways to make a gzipped tarred archive of the WebbPSF data suitable for distribution.

**If you are on the Institute network:** ::

   $ cd webbpsf/dev_utils/
   $ ./master_data_release.sh 0.X.Y

**If you're working from a local data root:** ::

   $ cd webbpsf/dev_utils/
   $ DATAROOT="/Users/you/webbpsf-data-sources/" ./make-data-sdist.sh 0.X.Y
   $ cp ./webbpsf-data-0.X.Y.tar.gz /where/ever/you/want/

Releasing new versions
======================

If you are making a release for `poppy` at the same time as a release in WebbPSF, do that first.
Update the dependency requirement to the new version of poppy, in ``webbpsf/pyproject.toml`.

When you are ready, proceed with the WebbPSF release as follows:

#. Get the `develop` branch into the state that you want, including all PRs merged, updated release notes. This includes all tests passing both locally and on GitHub Actions.
#. Tag the commit with `v<version>`, being sure to sign the tag with the `-s` option.
   * ``git tag -s v<version> -m "Release v<version>"``

#. Push tag to github, on `develop`
#. On github, make a PR from `develop` to `stable` (this can be done ahead of time and left open, until all individual PRs are merged into `develop`.).
#. After verifying that PR is complete and tests pass, merge it. (Once merged, both the `stable` and `develop` branches should match).
#. Release on Github:

   #. On Github, click on "[N] Releases".
   #. Select "Draft a new release".
   #. Specify the version number, title, and brief description of the release.
   #. Press "Publish Release".

#. Release to PyPI. This should now happen automatically on GitHub Actions. This will be triggered by a GitHub Actions build of a tagged commit on the `stable` branch, so it will happen automatically on the prior step for the PR into `stable`.

.. note::

  Once conda installation is working again, find this page in the documentation
  for version 1.0.0 and adapt the steps from the "Releasing a new version
  through AstroConda" section to the new process.

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


Finishing the release
=====================

 #. Email an announcement to ``webbpsf-users@maillist.stsci.edu``


