#!/bin/sh
# Script to make a distributable version of WebbPSF, with various packaging tweaks

VER="0.3.1"
TAR=/usr/bin/tar  # make sure to use the BSD version, required for the -L option
DATAROOT=/itar/jwst/tel/share/webbpsf/webbpsf-data-source

# Create the data tarfile
# make a copy of the filter file excluding FND to appease MIRI PI requirements
# Also exclude various things we don't want to distribute, like .svn, the old OPDs, and the data source directories

# set MIRI filters to the square profiles for broad distribution
# link to tophat filters
echo "Setting up for simplified MIRI filter profiles"
\rm -r $DATAROOT/MIRI/filters
ln -s $DATAROOT/MIRI/tophat_filters $DATAROOT/MIRI/filters


# create public distributable tar file
$TAR -cvz -L -C $DATAROOT/..  \
    --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store \
    --exclude sources --exclude "*FND*" --exclude "*_filters" --exclude "*py" \
    --exclude "_Obsolete" --exclude README_DEVEL.txt \
    -s "/webbpsf-data-source/webbpsf-data/" \
    -f webbpsf-data-$VER.tar.gz webbpsf-data-source


# Make a copy with the complete MIRI filter profiles, for internal or CoroWG use

# link to measured filters
echo "Setting up for measured MIRI filter profiles"
\rm $DATAROOT/MIRI/filters
ln -s $DATAROOT/MIRI/measured_filters $DATAROOT/MIRI/filters
# creat internal distributable tar file
$TAR -cvz -L -C $DATAROOT/..  \
    --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store \
    --exclude "*_filters" --exclude "*py" \
    --exclude "_Obsolete" --exclude README_DEVEL.txt \
    -s "/webbpsf-data-source/webbpsf-data/"  \
    -f webbpsf-data-internal-$VER.tar.gz webbpsf-data-source



#\cp webbpsf/dist/webbpsf-data-public-$VER.tar.gz ~/web/software/webbpsf/webbpsf-data-$VER.tar.gz

echo "Public file output to:    $PWD/webbpsf-data-$VER.tar.gz"
echo "Internal file output to:  $PWD/webbpsf-data-internal-$VER.tar.gz"

