#!/bin/bash 

# Call astyle and remove trailing whitespace on all c/cpp/h/hpp in current folder and subfolders

EXE=`dirname $0`/standardizeSourceFormat.sh

for i in `find . -name "*.c"` `find . -name "*.cpp"` `find . -name "*.h"` `find . -name "*.hpp"`; do
  ${EXE} $i;
done

# additionally, remove tab stops from python files
for i in `find . -name "*.py"`; do
    sed -i 's/\t/  /g' $i;
done
