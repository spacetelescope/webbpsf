WebbPSF: Simulated Point Spread Functions for the James Webb Space Telescope
============================================================================

WebbPSF produces simulated PSFs for the James Webb Space Telescope, NASA's next
flagship infrared space telescope. WebbPSF can simulate images for any of the
four science instruments plus the fine guidance sensor, including both direct
imaging and coronagraphic modes.

Developed by Marshall Perrin, Joseph Long, and collaborators, 2010-2015.

Documentation can be found online at https://pythonhosted.org/webbpsf/

WebbPSF requires a large amount of input data for its simulations, including
optical path difference (OPD) maps, filter transmission curves, and coronagraph
Lyot mask shapes. These data files are not included in this source distribution.
Please see the documentation to download the required data files.

This is intended to be an `Astropy <http://astropy.org/>`_ affiliated package.


Status reports for developers
-----------------------------

.. image:: https://pypip.in/v/webbpsf/badge.png
    :target: https://pypi.python.org/pypi/webbpsf

.. image:: https://pypip.in/d/webbpsf/badge.png
    :target: https://pypi.python.org/pypi/webbpsf

.. image:: https://travis-ci.org/mperrin/webbpsf.png?branch=master
    :target: https://travis-ci.org/mperrin/webbpsf
    :alt: Test Status

.. image:: https://coveralls.io/repos/mperrin/webbpsf/badge.svg
    :target: https://coveralls.io/r/mperrin/webbpsf
    :alt: Test Coverage
