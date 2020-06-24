WebbPSF: Simulated Point Spread Functions for the James Webb and Nancy Grace Roman Space Telescopes
===================================================================================================

.. image:: docs/readme_fig.png

.. image:: https://img.shields.io/pypi/v/webbpsf.svg
   :target: https://pypi.python.org/pypi/webbpsf
   :alt: Badge showing current released PyPI version

.. image:: https://travis-ci.org/spacetelescope/webbpsf.svg?branch=master
   :target: https://travis-ci.org/spacetelescope/webbpsf
   :alt: Badge showing continuous integration test status

.. image:: https://codecov.io/gh/spacetelescope/webbpsf/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/spacetelescope/webbpsf

.. image:: https://img.shields.io/badge/ascl-1504.007-blue.svg?colorB=262255
   :target: http://ascl.net/1504.007

WebbPSF produces simulated PSFs for the James Webb Space Telescope, NASA's next
flagship infrared space telescope. WebbPSF can simulate images for any of the
four science instruments plus the fine guidance sensor, including both direct
imaging and coronagraphic modes.

WebbPSF also supports simulating PSFs for the upcoming Nancy Grace Roman Space Telescope (formerly WFIRST),
including its Wide Field Instrument and a preliminary version of the Coronagraph Instrument.

Developed by Marshall Perrin, Joseph Long, Neil Zimmerman, Robel Geda, Shannon
Osborne, Marcio Melendez Hernandez, Lauren Chambers, Keira Brooks, Charles-Phillipe Lajoie, Jarron Leisenring, Alden Jurling, and collaborators, 2010-2020.

Documentation can be found online at https://webbpsf.readthedocs.io

WebbPSF requires input data for its simulations, including optical path
difference (OPD) maps, filter transmission curves, and coronagraph Lyot mask
shapes. These data files are not included in this source distribution.
Please see the documentation to download the required data files.
