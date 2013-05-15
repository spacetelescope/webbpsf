
Diagnostics & Troubleshooting
===============================


If something does not work right, webbpsf includes a few simple diagnostic tools. 

The first checks the installed versions of the required and optional packages. This can be helpful information to provide as part of a trouble report. ::

  >>> print webbpsf._system_diagnostic()
    OS: Darwin-10.8.0-i386-64bit
    Python version: 2.7.3 (default, Oct 23 2012, 16:00:23)  [GCC 4.2.1 (Apple Inc. build 5664)]
    poppy version: 0.3dev
    webbpsf version: 0.3dev

    tkinter version: 0.3.1
    wxpython version: 2.9.4.0

    astropy version: 0.3.dev3921
    pysynphot version: 0.8.3
    pyfits version: 3.2.dev-0.3.dev3921
    FFTW3 version: yes, present



