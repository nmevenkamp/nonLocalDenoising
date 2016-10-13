#!/bin/sh

# Export hg repository to EXPORTDIR and change to that directory
printf '* Running hg archive... '
if /usr/bin/hg archive -tfiles $EXPORTDIR
then
    printf 'ok\n'
else
    printf 'failed\n'
    exit 1
fi

printf "* Changing directory to $EXPORTDIR\n"
cd $EXPORTDIR

# Remove files related to internal (repository) structure
printf '* Removing HG remains\n'
/bin/rm .hgeol
/bin/rm .hgignore
# Use -f to ignore (potentially) non-existent file
/bin/rm -f .hg_archival
/bin/rm exclusion.opensource

# Add short license to files?
printf '* Adding license header to files\n'
./internal/util/addShortLicenseToFiles.py

# Remove internal code and remake folder
printf '* Removing internal code... '
/bin/rm -rf internal/
if [ -d internal/ ]
then
    printf 'failed to remove directory\n'
    exit 1
else
    printf 'ok\n'
fi

printf '* Removing makefile.local files... '
find -name makefile.local -print0 | xargs -0 -I {} /bin/rm {}
if find -name makefile.local | grep -q makefile.local
then
    printf 'failed\n'
    exit 1
else
    printf 'ok\n'
fi

# Fix exceptionTest
printf '* Modifying exceptionTest.txt\n'
sed -i 's/exceptionTest.cpp, line 8/exceptionTest.cpp, line 20/' selfTest/aol/runExceptionTest.sh

# Remove internal colorize script
printf '* Removing colorize script\n'
sed -i 's/ | perl ${QUOCPATH}\/internal\/util\/colorize.pl//' util/cmakeParseError.sh

# Turn on open source documentation flag
printf '* Enabling OSDOC\n'
sed -i 's/\(SET( DOXYGEN_ENABLED_SECTIONS ".*\)\(" )\)/\1 OSDOC\2/' CMakeLists.txt

printf '* Disabling cuda, grape, mercurial, metis\n'
# Disable USE_MERCURIAL, USE_CUDA, USE_GRAPE. TODO: remove CUDA option.
sed -i 's/OPTION ( USE_MERCURIAL "" ON )/OPTION ( USE_MERCURIAL "" OFF )/' CMakeLists.txt
sed -i 's/OPTION ( USE_CUDA "" ON )/OPTION ( USE_CUDA "" OFF )/' CMakeLists.txt
sed -i 's/OPTION ( USE_GRAPE "" ON )/OPTION ( USE_GRAPE "" OFF )/' CMakeLists.txt
sed -i 's/OPTION ( USE_METIS "" ON )/OPTION ( USE_METIS "" OFF )/' CMakeLists.txt

printf '* Removing internal files from doxygen exclude lists\n'
# doxygen.conf.cmake is not internal but may contain paths of internal files
awk 'm = /.*EXCLUDE.*/ { for(i = 1; i <= NF; i++) { if ( $i !~ /.*internal.*/ ) { if ( i != NF ) { printf "%s ",$i } else { printf "%s\n",$i } } } } !m { print }' doc/doxygen.conf.cmake > doc/doxygen.conf.cmake.edited
mv doc/doxygen.conf.cmake.edited doc/doxygen.conf.cmake

# Make cmake open source selection the default selection
printf '* Copying cmake.selection.opensource to cmake.selection.default\n'
cp cmake.selection.opensource cmake.selection.default

printf '* Generating package... '
if [ ! -d $PACKAGEDIR ]
then
    mkdir $PACKAGEDIR
    if [ $? != 0 ]
    then
	printf "Could not create $PACKAGEDIR\n"
        exit 1
    fi
fi
# Create tar.bz2 archive
DIRNAME=$(basename $PWD)
cd ..
/bin/tar -cjf "$PACKAGEDIR/quocmesh-opensource-nightly.tar.bz2" $DIRNAME
cd $PACKAGEDIR
# Create Checksums
/usr/bin/sha256sum quocmesh-opensource-nightly.tar.bz2 > SHA256SUM
/usr/bin/sha512sum quocmesh-opensource-nightly.tar.bz2 > SHA512SUM
printf 'ok\n'

printf 'QuocMesh open source export complete\n'
exit 0
