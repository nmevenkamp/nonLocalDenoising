#!/bin/sh

# ATTENTION: Before using, please give correct paths:
export SOURCEDIR='/path/to/quocmesh-source-dir'
export EXPORTDIR='/path/to/quocmesh-opensource'
export PACKAGEDIR='/path/to/quocmesh-nightly-package-dir'
export BUILDDIR='/path/to/quocmesh-opensource-build'
EXPORTSCRIPT='/path/to/generateOpenSource.sh'
BUILDSCRIPT='/path/to/compileOpenSource.sh'

# Clean up last export and build
if [ -d $EXPORTDIR ]
then
    rm -rf $EXPORTDIR
fi

if [ -d $BUILDDIR ]
then
    rm -rf $BUILDDIR
fi

# Create new export and build
cd $SOURCEDIR

printf '== QuocMesh open source export ==\n'
sh $EXPORTSCRIPT
if [ $? != 0 ]
then
    printf 'QuocMesh open source export failed!\n'
    exit 1
fi
printf '\n'

printf '== QuocMesh open source build ==\n'
sh $BUILDSCRIPT
if [ $? != 0 ]
then
     printf 'QuocMesh open source build failed!\n'
     exit 1
fi

unset SOURCEDIR
unset EXPORTDIR
unset PACKAGEDIR
unset BUILDDIR
