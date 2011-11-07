#!/bin/sh
# Script to make a distributable version of WebbPSF, with various packaging tweaks

VER="0.2.6"
TAR=/usr/bin/tar  # make sure to use the BSD version, required for the -L option

# Create a source distribution
#   use python sdist, but then re-make the tar file so we can
#   include the stsci_distutils_hack and defsetup files
python setup.py sdist
cd dist
tar xvzf webbpsf-$VER.tar.gz
\cp ../stsci_distutils_hack.py webbpsf-$VER
\cp ../defsetup.py webbpsf-$VER
$TAR -cvz --exclude .svn --exclude '*.pyc' -f webbpsf-$VER.tar.gz webbpsf-$VER
\cp webbpsf-$VER.tar.gz ~/web/software/webbpsf
mv webbpsf-$VER.tar.gz ..
cd ..


