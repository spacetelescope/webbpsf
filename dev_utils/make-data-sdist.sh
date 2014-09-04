#!/bin/sh
# Script to make a distributable version of WebbPSF, with various packaging tweaks

VER="0.3.0"
TAR=/usr/bin/tar  # make sure to use the BSD version, required for the -L option


# Create the data tarfile
# make a copy of the filter file excluding FND to appease Gillian Wright, and exclude the FND from the tar file
# Also exclude various things we don't want to distribute, like .svn, the old OPDs, and the data source directories

cd ~/software/
cp webbpsf-data/filters.txt webbpsf
mv webbpsf-data/filters.txt .

# remove FND, and set MIRI filters to the square profiles
cat filters.txt | grep -v FND >  webbpsf-data/filters.txt
# link to tophat filters
\rm -r webbpsf-data/MIRI/filters
ln -s tophat_filters webbpsf-data/MIRI/filters
# create public distributable tar file
$TAR -cvz -L  --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store --exclude sources --exclude "*FND*" --exclude "*_filters" --exclude "*py" --exclude "_Obsolete" -f webbpsf/dist/webbpsf-data-public-$VER.tar.gz webbpsf-data


# Make a copy with more real data, for internal or CoroWG use

rm webbpsf-data/filters.txt # put the real filters file back
mv filters.txt webbpsf-data/filters.txt
# link to measured filters
\rm webbpsf-data/MIRI/filters
ln -s measured_filters webbpsf-data/MIRI/filters
# create internal distributable tar file
$TAR -cvz -L  --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store -f webbpsf/dist/webbpsf-data-internal-$VER.tar.gz webbpsf-data



\cp webbpsf/dist/webbpsf-data-public-$VER.tar.gz ~/web/software/webbpsf/webbpsf-data-$VER.tar.gz

