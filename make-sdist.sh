#!/bin/sh

python setup.py sdist
cd dist
tar xvzf webbpsf-0.2.tar.gz
\cp ../stsci_distutils_hack.py webbpsf-0.2
\cp ../defsetup.py webbpsf-0.2
tar -cvz --exclude .svn  -f webbpsf-0.2.tar.gz webbpsf-0.2
\cp webbpsf-0.2.tar.gz ~/web/software/webbpsf
cd ..
mv data webbpsf-data
tar -cvz  --exclude .svn --exclude OPD_RevT --exclude .DS_Store -f webbpsf-data-0.2.tar.gz webbpsf-data
mv  webbpsf-data data
\cp webbpsf-data-0.2.tar.gz ~/web/software/webbpsf

