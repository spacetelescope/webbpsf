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

echo
echo "Copying latest data to /grp/jwst/ote for internal stsci use..."
main_directory="/grp/jwst/ote"
new_directory="$main_directory/webbpsf-data-$VER"
symlink_directory="/grp/jwst/ote/webbpsf-data"

cp "$PWD/webbpsf-data-$VER.tar.gz" "$main_directory"
tar -xzf "$PWD/webbpsf-data-$VER.tar.gz" -C "$main_directory"
mv "$main_directory/webbpsf-data" "$new_directory"
ln -s "$new_directory" "$symlink_directory"

./make-minimal-datafiles.py  ${PWD}/webbpsf-data-${VER}.tar.gz $VER

echo
echo "================================================="
echo "Data extracted for internal use with updated symlink  $symlink_directory -> $new_directory"
echo
echo "OUTPUT FILES:"
echo
echo ${PWD}/webbpsf-data-${VER}.tar.gz
echo ~/tmp/minimal-webbpsf-data-${VER}/minimal-webbpsf-data-${VER}.tar.gz
echo
echo You probably want to test if those look as expected, and if so then copy into the Box folder 'webbpsf_data_public'
echo "================================================="
echo

