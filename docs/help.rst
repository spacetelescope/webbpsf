Diagnostics & Troubleshooting
=============================

If something does not work right, the first place to look is the :ref:`known_issues` section of the release notes. The next place to check is the `GitHub issues <https://github.com/spacetelescope/webbpsf/issues>`_ page, where another user may have reported the problem.

To report a new issue, you will need a free GitHub account. Alternatively, you may report the issue via email to the project maintainers. Include code that exhibits the issue to facilitate debugging.

WebbPSF includes a helper function that will return a report with information that may be useful for troubleshooting. An example of its usage is given below::

   In [1]: import webbpsf
   WebbPSF log messages of level INFO and above will be shown.
   WebbPSF log outputs will be directed to the screen.

   In [2]: print webbpsf.system_diagnostic()

   OS: Darwin-13.4.0-x86_64-i386-64bit
   Python version: 2.7.8 (default, Oct  2 2014, 13:50:25)  [GCC 4.2.1 Compatible Apple LLVM 6.0 (clang-600.0.51)]
   numpy version: 1.9.1
   poppy version: 0.3.3.dev335
   webbpsf version: 0.3rc4

   tkinter version: 0.3.1
   wxpython version: not found

   astropy version: 0.4.2
   pysynphot version: 0.9.6
   pyFFTW version: 0.9.2

   Floating point type information for numpy.float_:
   Machine parameters for float64
   ---------------------------------------------------------------------
   precision= 15   resolution= 1.0000000000000001e-15
   machep=   -52   eps=        2.2204460492503131e-16
   negep =   -53   epsneg=     1.1102230246251565e-16
   minexp= -1022   tiny=       2.2250738585072014e-308
   maxexp=  1024   max=        1.7976931348623157e+308
   nexp  =    11   min=        -max
   ---------------------------------------------------------------------

   Floating point type information for numpy.complex_:
   Machine parameters for float64
   ---------------------------------------------------------------------
   precision= 15   resolution= 1.0000000000000001e-15
   machep=   -52   eps=        2.2204460492503131e-16
   negep =   -53   epsneg=     1.1102230246251565e-16
   minexp= -1022   tiny=       2.2250738585072014e-308
   maxexp=  1024   max=        1.7976931348623157e+308
   nexp  =    11   min=        -max
   ---------------------------------------------------------------------

