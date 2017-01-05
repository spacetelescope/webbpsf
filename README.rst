WebbPSF: Simulated Point Spread Functions for JWST and WFIRST
=============================================================

.. image:: docs/readme_fig.png

.. image:: https://img.shields.io/pypi/v/webbpsf.svg
   :target: https://pypi.python.org/pypi/webbpsf
   :alt: Badge showing current released PyPI version

.. image:: https://travis-ci.org/mperrin/webbpsf.svg?branch=master
   :target: https://travis-ci.org/mperrin/webbpsf
   :alt: Badge showing continuous integration test status

.. image:: https://coveralls.io/repos/github/mperrin/webbpsf/badge.svg?branch=master
   :target: https://coveralls.io/github/mperrin/webbpsf?branch=master

WebbPSF produces simulated PSFs for the James Webb Space Telescope, NASA's next
flagship infrared space telescope. WebbPSF can simulate images for any of the
four science instruments plus the fine guidance sensor, including both direct
imaging and coronagraphic modes.

WebbPSF also supports simulating PSFs for the upcoming Wide Field Infrared Survey Telescope (WFIRST),
including its Wide Field Instrument and a preliminary version of the Coronagraph Instrument.

Developed by Marshall Perrin, Joseph Long, Neil Zimmerman, and collaborators, 2010-2017.

Documentation can be found online at https://pythonhosted.org/webbpsf/

WebbPSF requires input data for its simulations, including optical path
difference (OPD) maps, filter transmission curves, and coronagraph Lyot mask
shapes. These data files are not included in this source distribution.
Please see the documentation to download the required data files.

*This is intended to be an* `Astropy <http://astropy.org/>`_ *affiliated package.*
