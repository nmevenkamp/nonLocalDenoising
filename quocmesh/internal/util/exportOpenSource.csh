#!/bin/tcsh

# call from main quocmesh directory of a clean checkout

# todo: cross-check changes of makefile.selection.default/opensource

setenv SOURCEDIR `pwd`
setenv EXPORTDIR ../quocmesh-opensource-01.3/

# first copy whole doc directory, cleaned up below
svn cp                 doc                 $EXPORTDIR/

# "create" empty directories
mkdir $EXPORTDIR/cmake                && svn add $EXPORTDIR/cmake
mkdir $EXPORTDIR/examples             && svn add $EXPORTDIR/examples
mkdir $EXPORTDIR/exampleProjects      && svn add $EXPORTDIR/exampleProjects
mkdir $EXPORTDIR/external             && svn add $EXPORTDIR/external
mkdir $EXPORTDIR/modules              && svn add $EXPORTDIR/modules
mkdir $EXPORTDIR/selfTest             && svn add $EXPORTDIR/selfTest
mkdir $EXPORTDIR/tools                && svn add $EXPORTDIR/tools
mkdir $EXPORTDIR/util                 && svn add $EXPORTDIR/util


# examples OK for release
svn cp examples/adaptiveGrids                       $EXPORTDIR/examples/
svn cp examples/AmbrosioTortorelli                  $EXPORTDIR/examples/
svn cp examples/anisotropicDiffusion                $EXPORTDIR/examples/
svn cp examples/anisotropyVisualization             $EXPORTDIR/examples/
svn cp examples/BartonNackman                       $EXPORTDIR/examples/
svn cp examples/cimg                                $EXPORTDIR/examples/
svn cp examples/debugging                           $EXPORTDIR/examples/
svn cp examples/DirichletBCs                        $EXPORTDIR/examples/
svn cp examples/foxToolkit                          $EXPORTDIR/examples/
svn cp examples/heatequation                        $EXPORTDIR/examples/
svn cp examples/imageManipulation                   $EXPORTDIR/examples/
svn cp examples/ipopt                               $EXPORTDIR/examples/
svn cp examples/multigrid                           $EXPORTDIR/examples/
svn cp examples/MumfordShahSegmentation             $EXPORTDIR/examples/
svn cp examples/mutualInformationRegistration       $EXPORTDIR/examples/
svn cp examples/Newton                              $EXPORTDIR/examples/
svn cp examples/openmesh                            $EXPORTDIR/examples/
svn cp examples/openmp                              $EXPORTDIR/examples/
svn cp examples/plotter                             $EXPORTDIR/examples/
svn cp examples/rowIterator                         $EXPORTDIR/examples/
svn cp examples/RudinOsherFatemiModel               $EXPORTDIR/examples/
svn cp examples/simultaneouslyAssembled             $EXPORTDIR/examples/
svn cp examples/sseTest                             $EXPORTDIR/examples/
svn cp examples/sweepingSignedDist                  $EXPORTDIR/examples/
svn cp examples/testdata                            $EXPORTDIR/examples/
svn cp examples/tpcfe                               $EXPORTDIR/examples/
svn cp examples/tvSmoothing                         $EXPORTDIR/examples/
svn cp examples/utils                               $EXPORTDIR/examples/
svn cp examples/vectorMatrixOps                     $EXPORTDIR/examples/
svn cp examples/virtualBaseClasses                  $EXPORTDIR/examples/


# externals not containing code
svn cp external/ahmed                               $EXPORTDIR/external
svn cp external/suitesparse                         $EXPORTDIR/external
svn cp external/cimg                                $EXPORTDIR/external
svn cp external/fox                                 $EXPORTDIR/external
svn cp external/grape                               $EXPORTDIR/external
svn cp external/hypre                               $EXPORTDIR/external
svn cp external/ipopt                               $EXPORTDIR/external
svn cp external/openmesh                            $EXPORTDIR/external
svn cp external/qt                                  $EXPORTDIR/external
svn cp external/vtk                                 $EXPORTDIR/external
svn cp external/xml2                                $EXPORTDIR/external
# see below for externals which will not reveive our license block

# modules OK for release
svn cp modules/aol                                  $EXPORTDIR/modules/
svn cp modules/eikonal                              $EXPORTDIR/modules/
svn cp modules/grape                                $EXPORTDIR/modules/
svn cp modules/multigrid                            $EXPORTDIR/modules/
svn cp modules/narrowband                           $EXPORTDIR/modules/
svn cp modules/openmesh                             $EXPORTDIR/modules/
svn cp modules/qcsm                                 $EXPORTDIR/modules/
svn cp modules/quoc                                 $EXPORTDIR/modules/
svn cp modules/tpcfe                                $EXPORTDIR/modules/


# finished projects OK for release
svn cp finishedProjects/AmbrosioTortorelli3D        $EXPORTDIR/exampleProjects/
# svn cp finishedProjects/AllenCahnShapeOptimization  $EXPORTDIR/exampleProjects/ # this currently does not work
svn cp finishedProjects/boneelast                   $EXPORTDIR/exampleProjects/
svn cp finishedProjects/mcm2d                       $EXPORTDIR/exampleProjects/
# svn cp finishedProjects/shapeGeodesics              $EXPORTDIR/exampleProjects/ # this currently does not work
# svn cp finishedProjects/shapeStatistics             $EXPORTDIR/exampleProjects/ # this currently does not work
svn cp finishedProjects/tpcfe                       $EXPORTDIR/exampleProjects/
svn cp finishedProjects/vecfieldInpainting          $EXPORTDIR/exampleProjects/
svn cp finishedProjects/vesseltree                  $EXPORTDIR/exampleProjects/
svn cp finishedProjects/WulffShapeLevelsetAniso     $EXPORTDIR/exampleProjects/

# current projects OK for release
svn cp projects/vtkfox                              $EXPORTDIR/exampleProjects/

# selfTests for modules released
svn cp selfTest/aol                                 $EXPORTDIR/selfTest/
svn cp selfTest/eikonal                             $EXPORTDIR/selfTest/
svn cp selfTest/grape                               $EXPORTDIR/selfTest/
svn cp selfTest/multigrid                           $EXPORTDIR/selfTest/
svn cp selfTest/narrowband                          $EXPORTDIR/selfTest/
svn cp selfTest/openmesh                            $EXPORTDIR/selfTest/
svn cp selfTest/qcsm                                $EXPORTDIR/selfTest/
svn cp selfTest/quoc                                $EXPORTDIR/selfTest/
svn cp selfTest/tpcfe                               $EXPORTDIR/selfTest/

# tools OK for release
svn cp tools/benchmark                              $EXPORTDIR/tools/
svn cp tools/grape                                  $EXPORTDIR/tools/
svn cp tools/imageManipulation                      $EXPORTDIR/tools/


# utils OK for release
svn cp util/build-clang-3.0.sh                      $EXPORTDIR/util/
svn cp util/build-gcc-4.6.0.sh                      $EXPORTDIR/util/
svn cp util/cmakeParseError.sh                      $EXPORTDIR/util/
svn cp util/configure                               $EXPORTDIR/util/
svn cp util/detailedTableOfProjects.sh              $EXPORTDIR/util/
svn cp util/findRawString.py                        $EXPORTDIR/util/
svn cp util/remove_templates.pl                     $EXPORTDIR/util/
svn cp util/removeTrailingWhitespace.py             $EXPORTDIR/util/
svn cp util/replaceExpressions.py                   $EXPORTDIR/util/
svn cp util/replaceTabStops.sh                      $EXPORTDIR/util/
svn cp util/runtest.sh                              $EXPORTDIR/util/
svn cp util/setSvnProperties.sh                     $EXPORTDIR/util/
svn cp util/standardizeSourceFormatEverywhere.sh    $EXPORTDIR/util/
svn cp util/standardizeSourceFormat.sh              $EXPORTDIR/util/
svn cp util/suppress.conf                           $EXPORTDIR/util/
svn cp util/suppress.pl                             $EXPORTDIR/util/
svn cp util/tableOfProjects.sh                      $EXPORTDIR/util/
svn cp util/testIndividualInclusion.sh              $EXPORTDIR/util/


# files in main directory OK for release
svn cp makefile                                     $EXPORTDIR/
svn cp makefile.config.template                     $EXPORTDIR/
svn cp makefile.config.default                      $EXPORTDIR/
svn cp makefile.project                             $EXPORTDIR/
svn cp makefile.selection.opensource                $EXPORTDIR/makefile.selection.default
svn cp CMakeLists.project.txt                       $EXPORTDIR/
svn cp CMakeLists.txt                               $EXPORTDIR/
svn cp cmake.selection.opensource                   $EXPORTDIR/cmake.selection.default
svn cp cmake.includeTest.cpp.in                     $EXPORTDIR/
svn cp cmake.mgw.bat.in                             $EXPORTDIR/

# copy license file
svn cp util/cddl.txt                                $EXPORTDIR/

cp util/addShortLicenseToFiles.py $EXPORTDIR/util/

cd $EXPORTDIR

# remove unused modules from documentation
sed -i 's/<li>\\ref surfdoc "Using module surf" <br>Needs modules %aol, %multigrid, %quoc<\/li>//' modules/aol/mainpage.dox
sed -i 's/<li>\\ref bemdoc "Using module bem" <br>Needs modules %aol, %quoc<\/li>//' modules/aol/mainpage.dox
sed -i 's/<li>\\ref DTGriddoc "Using the Dynamic Tubular Grid (dtgrid)" <br>Needs modules %aol, %multigrid, %quoc, %tpcfe<\/li>//' modules/aol/mainpage.dox

# modify documentation
svn rm doc/DocOOProgramming.tex
svn rm doc/DocPrinciplesSVN.tex
svn rm doc/DocSelectedProblems.tex
svn rm doc/DocBoundaryElements.tex
svn rm doc/DocNewtonMethod.tex
svn rm doc/DocTimestepping.tex
svn rm doc/DocGradientflow.tex

sed -i 's/\\input{DocOOProgramming.tex}//' doc/manual.tex
sed -i 's/\\input{DocPrinciplesSVN.tex}//' doc/manual.tex
sed -i 's/\\input{DocSelectedProblems.tex}//' doc/manual.tex
sed -i 's/\\input{DocBoundaryElements.tex}//' doc/manual.tex
sed -i 's/\\input{DocNewtonMethod.tex}//' doc/manual.tex
sed -i 's/\\input{DocTimestepping.tex}//' doc/manual.tex
sed -i 's/\\input{DocGradientflow.tex}//' doc/manual.tex

util/addShortLicenseToFiles.py
rm util/addShortLicenseToFiles.py

# the following contain third-party code which we may redistribute and should not receive our license
cd $SOURCEDIR
svn cp external/bz2                                 $EXPORTDIR/external
svn cp external/dirent                              $EXPORTDIR/external
svn cp external/vtkfox                              $EXPORTDIR/external
# cmake modules OK for release
svn cp cmake/FindSUITESPARSE.cmake                  $EXPORTDIR/cmake/
cd $EXPORTDIR

sed -i 's/ | perl ${QUOCPATH}\/internal\/util\/colorize.pl//' util/cmakeParseError.sh
sed -i 's/file exceptionTest.cpp, line 8/file exceptionTest.cpp, line 18/' selfTest/aol/exceptionTest.txt

# files potentially to be modified (touch reminders that should have svn status "?") -- this needs to be done manually

touch makefile.selection.default.modify
echo "modify makefile.selection.default (remove comments)!"
touch cmake.selection.default.modify
echo "modify cmake.selection.default (remove comments)!"
touch CMakeLists.txt.modify
echo "modify CMakeLists.txt (remove cuda, triangle related stuff)!"


# adapt QUOCMESH_VERSION string in CMakeLists.txt, set -DUSE_MERCURIAL=0 and -DUSE_SUBVERSION=0 s.t. it is not overwritten by cmake
# generate documentation (make doclib, make docall)
