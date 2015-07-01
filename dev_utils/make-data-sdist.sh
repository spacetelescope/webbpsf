#!/bin/sh
# Script to make a distributable version of WebbPSF, with various packaging tweaks
set -e

if ! [[ $1 ]]; then
  echo "Provide a version string, e.g.:"
  echo "    ./make-data-sdist.sh 0.3.3"
  exit 1
fi

if ! [[ $DATAROOT ]]; then
  DATAROOT=/itar/jwst/tel/share/webbpsf/webbpsf-data-source/
fi
echo "Using data from $DATAROOT"

# If on Mac OS, tell tar to not include ._* files for
# HFS-specific extended attributes
export COPYFILE_DISABLE=1

TMPDIR="/tmp/webbpsf-data"

mkdir -p  $TMPDIR
rsync -avz --exclude '._*' --exclude '_Obsolete' "$DATAROOT" $TMPDIR

VER=$1
echo "$VER" > $TMPDIR/version.txt
echo "Saving version number $VER to version.txt"

# Create the data tarfile
# make a copy of the filter file excluding FND to appease MIRI PI requirements
# Also exclude various things we don't want to distribute, like .svn, the old OPDs, and the data source directories

# set MIRI filters to the square profiles for broad distribution
# link to tophat filters
echo "Setting up for simplified MIRI filter profiles"
rm -fv $TMPDIR/MIRI/filters
ln -s $TMPDIR/MIRI/tophat_filters $TMPDIR/MIRI/filters

# create public distributable tar file
tar -cvz -C $TMPDIR/..  \
    --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store \
    --exclude sources --exclude "*FND*" --exclude "*_filters" --exclude "*py" \
    --exclude "_Obsolete" --exclude README_DEVEL.txt \
    -f "webbpsf-data-$VER.tar.gz" webbpsf-data

# Make a copy with the complete MIRI filter profiles, for internal or CoroWG use

# link to measured filters
echo "Setting up for measured MIRI filter profiles"
rm -fv $TMPDIR/MIRI/filters
ln -s $TMPDIR/MIRI/measured_filters $TMPDIR/MIRI/filters

# create internal distributable tar file
tar -cvz -C $TMPDIR/..  \
    --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store \
    --exclude "*_filters" --exclude "*py" \
    --exclude "_Obsolete" --exclude README_DEVEL.txt \
    -f "webbpsf-data-internal-$VER.tar.gz" webbpsf-data-source

echo "Public file output to:    $(pwd)/webbpsf-data-$VER.tar.gz"
echo "Internal file output to:  $(pwd)/webbpsf-data-internal-$VER.tar.gz"

