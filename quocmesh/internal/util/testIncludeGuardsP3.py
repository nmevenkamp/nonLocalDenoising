#!/usr/bin/python

import os
import sys

# replace one or more expressions
def checkIncludeGuards(arg, drctry, files):
  os.chdir(drctry)
  drctry = drctry.split("quocmesh")
  drctry = drctry[-1]
  for datei in files:
    if os.path.isfile(datei) and (datei.lower()).endswith('.h'):
      f = open ( datei, 'r' )
      test = f.readlines()
      # Unify the line endings, this way a Windows checkout should pass the test under Linux.
      test[0] = test[0].replace("\r\n","\n")
      test[1] = test[1].replace("\r\n","\n")
      firstLineShouldBe = "#ifndef __" + datei.upper() + "\n"
      firstLineShouldBe = firstLineShouldBe.replace(r".","_")
      secondLineShouldBe = "#define __" + datei.upper() + "\n"
      secondLineShouldBe = secondLineShouldBe.replace(r".","_")
      if ( firstLineShouldBe != test[0] ):
        print ( firstLineShouldBe )
        print ( test[0] )
        print ( datei + " has incorrect include guard" )
        sys.exit(1)
      if ( secondLineShouldBe != test[1] ):
        print ( datei + " has incorrect include guard" )
        sys.exit(1)


# the main program: start the walk
for dirpath, dirnames, filenames in os.walk(os.getcwd()):
  checkIncludeGuards(4, dirpath, dirnames + filenames)
