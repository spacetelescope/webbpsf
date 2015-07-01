#!/bin/sh
# Script to make a distributable version of WebbPSF, with various packaging tweaks
set -e

if ! [[ $VER ]]; then
  echo "Provide a version, e.g.:"
  echo "    VER=\"0.3.3\" ./make-data-sdist.sh"
  exit 1
fi

if ! [[ $DATAROOT ]]; then
  DATAROOT=/itar/jwst/tel/share/webbpsf/webbpsf-data-source
fi

if ! [[ $DEST ]]; then
  DEST="$(pwd)"
fi

# If on Mac OS, tell tar to not include ._* files for
# HFS-specific extended attributes
export COPYFILE_DISABLE=1

TMPDIR="/tmp/webbpsf-data"

mkdir $TMPDIR
rsync -avz "$DATAROOT" $TMPDIR

# Remove existing ._* files
find $TMPDIR -name "._*" -exec rm -i {} \;

# Create the data tarfile
# make a copy of the filter file excluding FND to appease MIRI PI requirements
# Also exclude various things we don't want to distribute, like .svn, the old OPDs, and the data source directories

# set MIRI filters to the square profiles for broad distribution
# link to tophat filters
echo "Setting up for simplified MIRI filter profiles"
\rm $TMPDIR/MIRI/filters
ln -s $TMPDIR/MIRI/tophat_filters $TMPDIR/MIRI/filters


# create public distributable tar file
tar -cvz -L -C $TMPDIR/..  \
    --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store \
    --exclude sources --exclude "*FND*" --exclude "*_filters" --exclude "*py" \
    --exclude "_Obsolete" --exclude README_DEVEL.txt \
    -f "webbpsf-data-$VER.tar.gz" webbpsf-data

cp "$TMPDIR/webbpsf-data-$VER.tar.gz" "$DEST/"

# Make a copy with the complete MIRI filter profiles, for internal or CoroWG use

# link to measured filters
echo "Setting up for measured MIRI filter profiles"
\rm $TMPDIR/MIRI/filters
ln -s $TMPDIR/MIRI/measured_filters $TMPDIR/MIRI/filters
# creat internal distributable tar file
tar -cvz -L -C $TMPDIR/..  \
    --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store \
    --exclude "*_filters" --exclude "*py" \
    --exclude "_Obsolete" --exclude README_DEVEL.txt \
    -f "webbpsf-data-internal-$VER.tar.gz" webbpsf-data-source

cp "$TMPDIR/webbpsf-data-internal-$VER.tar.gz" "$DEST/"

echo "Public file output to:    $DEST/webbpsf-data-$VER.tar.gz"
echo "Internal file output to:  $DEST/webbpsf-data-internal-$VER.tar.gz"

