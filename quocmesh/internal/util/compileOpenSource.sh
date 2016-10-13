#!/bin/sh

if [ ! -d $BUILDDIR ]
then
    mkdir $BUILDDIR
    if [ $? != 0 ]
    then
	printf "Could not create $BUILDDIR\n"
	exit 1
    fi
fi
cd $BUILDDIR

printf '* Running CMake\n'
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_STLDEBUG=1 -DGENERATE_INCLUDE_TEST=1 $EXPORTDIR >build.log 2>&1

if [ $? != 0 ]
then
    printf 'CMake failed!\n'
    exit 1
fi

printf '* Starting compilation\n'
time make -j3  >>build.log 2>&1
if [ $? != 0 ]
then
    printf 'Build failed!\n'
    exit 1
fi

printf '* Running tests\n'
time make -j3 test >>build.log 2>&1
if [ $? != 0 ]
then
    printf 'Tests failed!\n'
    exit 1
fi

# doxygen is required to build the documentation
printf '* Running CMake (doxygen)\n'
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_STLDEBUG=1 -DGENERATE_INCLUDE_TEST=1 -DUSE_DOXYGEN=1 $EXPORTDIR >>build.log 2>&1

if [ $? != 0 ]
then
    printf 'CMake (doxygen) failed!\n'
    exit 1
fi

printf '* Building documentation\n'
make doclib >>build.log 2>&1
if [ $? != 0 ]
then
    printf 'make doclib failed!\n'
    exit 1
fi
make docall >>build.log 2>&1
if [ $? != 0 ]
then
    printf 'make docall failed!\n'
    exit 1
fi

# Copy the documentation to the package directory and compress.
printf '* Copying and compressing documentation\n'
/bin/cp -ar ./doc $PACKAGEDIR
cd $PACKAGEDIR
/bin/tar -cjf ./quocmesh-opensource-nightly-doc.tar.bz2 ./doc

# Create Checksums
/usr/bin/sha256sum quocmesh-opensource-nightly-doc.tar.bz2 >> SHA256SUM
/usr/bin/sha512sum quocmesh-opensource-nightly-doc.tar.bz2 >> SHA512SUM

printf 'QuocMesh open source build complete\n'
exit 0
