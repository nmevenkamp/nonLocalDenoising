# Non-local denoising
This repository contains implementations of the Non-local means and BM3D denoising algorithms.
The implementation is based on the [QuocMesh software library](http://numod.ins.uni-bonn.de/software/quocmesh/).

In order to COMPILE the code under UNIX, follow these steps:

1) In quocGCC/go.sh, specify path to compiler binaries (CC=... & CXX=...), e.g. GCC
2) execute go.sh
3) run make

In order to RUN the denoising executables, follow these steps:

1) (Optional) adjust the parameter files (nlm.par, bm3d.par, localFilter.par)
2) Execute quocGCC/internal/projects/mevenkamp/<Filter> <pathToParameterFile>
   where <Filter> can be either "BM3D", "NLM" or "LocalFilter" (without the "")

The results will be stored in a new folder inside the output directory (which is specified in the parameter file).
The name of the new folder is a combination of the base name of the input file, the filter name and its parameters.

Results are stored in the QuocMesh double binary image format (.q2bz) and can be converted to png using
    quocGCC/tools/image/converter/convertMultipleQuoc2DToMaxContrastPNG <pathToQ2bzFile>
