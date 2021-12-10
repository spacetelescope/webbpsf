#!/bin/bash
# Top-level script to make a distributable version of the data files

# See /itar/jwst/tel/share/webbpsf/webbpsf-data-source/README_DEVEL.txt

if ! [[ $1 ]]; then
  echo "Provide a version string, e.g.:"
  echo "    ./master_data_release.sh  0.3.3"
  exit 1
fi

VER="$1"
TMPDIR="/tmp/webbpsf-data"

./make-data-sdist.sh $VER

./make-minimal-datafiles.py  ${PWD}/webbpsf-data-${VER}.tar.gz $VER


echo
echo "================================================="
echo "OUTPUT FILES:"
echo
echo ${PWD}/webbpsf-data-${VER}.tar.gz
echo ~/tmp/minimal-webbpsf-data-${VER}/minimal-webbpsf-data-${VER}.tar.gz
echo
echo You probably want to test if those look as expected, and if so then copy into the Box folder 'webbpsf_data_public'
echo "================================================="
echo

