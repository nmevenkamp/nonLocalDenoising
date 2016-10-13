#!/bin/bash
export CC="/Applications/MacPorts/bin/gcc-mp-6"
export CXX="/Applications/MacPorts/bin/g++-mp-6"
cmake -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=1 -DDYNAMIC_LINKING=1 -DUSE_FFTW=1 -DUSE_PNG=1 -DUSE_TIFF=1 -DUSE_BLAS=1 -DUSE_LAPACK=1 -DUSE_BOOST=1 -DBUILD_AND_USE_WAVELET=1 -DUSE_CIMG=1 ../quocmesh