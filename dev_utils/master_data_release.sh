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

./make-minimal-datafiles.py  ${PWD}/webbpsf-data-${VER}.tar.gz


echo "Copying webbpsf-data-${VER}.tar.gz   ==>>  /grp/webpages/mperrin/software/webbpsf"
\cp ${PWD}/webbpsf-data-${VER}.tar.gz /grp/webpages/mperrin/software/webbpsf


echo "Copying minimal-webbpsf-data.tar.gz   ==>> /grp/webpages/mperrin/software/webbpsf"
\cp /Users/mperrin/tmp/minimal-webbpsf-data/minimal-webbpsf-data.tar.gz /grp/webpages/mperrin/software/webbpsf



