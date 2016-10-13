#!/bin/bash

# Call astyle (must be in path) and remove trailing whitespace on individual file

EXEA=astyle
EXET=`dirname $0`/removeTrailingWhitespace.py

echo ${EXEA}

if [ ! -e $1 ]; then
  echo "$1 does not exist...";
  exit 1;
fi

${EXEA} -M80 -m0 -o -O --mode=c --pad-paren --convert-tabs --pad-oper --indent=spaces=2 --brackets=attach --indent-switches $1
${EXET} $1

