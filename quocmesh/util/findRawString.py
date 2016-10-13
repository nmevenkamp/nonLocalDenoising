#!/usr/bin/python

# findRawString.py
# script traverses the current directory tree and calls 'findRawString'
# on every .h or .cpp-file. If the expression that is specified here in
# the source code is found, the name of the file is printed.

# The expression is specified as a raw string, i.e. you can use any character in it,
# it is not interpreted as some regular expression.

from os import *
from os.path import *
from re import *


# searches for the expression in file arg in directory drctry
def findRawString(arg, drctry, files):
  chdir(drctry)
  for datei in files:
    if isfile(datei) and ( (datei.lower()).endswith('.cpp') or (datei.lower()).endswith('.h')):
      f = file(datei, 'r')        
      text = f.read()             # load file into a string
      f.close()
      # ATTENTION: Don't change and commit this file. To avoid an accidental commit
      # you can for example copy this script to your quocmesh-directory and only edit
      # this local copy.
      if text.find(r"Complicated Example String *\n \t[M]") > 0:
        print "Found expression in file " + drctry + "/" + datei + "..."
      

# the main program: start the walk
print "Searching..."
walk(getcwd(), findRawString, 4)


print "\n\nready"
