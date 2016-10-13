#!/usr/bin/python

from os import *
from os.path import *
from re import *

slic = ["Copyright (c) 2001-2014 AG Rumpf, INS, Universitaet Bonn                     ",
        "                                                                             ",
        "The contents of this file are subject to the terms of the Common Development ",
        "and Distribution License Version 1.0 (the \"License\"); you may not use      ",
        "this file except in compliance with the License. You may obtain a copy of    ",
        "the License at http://www.opensource.org/licenses/CDDL-1.0                   ",
        "                                                                             ",
        "Software distributed under the License is distributed on an \"AS IS\" basis, ",
        "WITHOUT WARRANTY OF ANY KIND, either expressed or implied.                   "];

def addLicenseBlock(arg, drctry, files):
  
  chdir(drctry)
  drctry = drctry.split("quocmesh")
  drctry = drctry[-1]
  for datei in files:
    if isfile(datei) and ( (datei.lower()).endswith('.cpp') or (datei.lower()).endswith('.h')):
      f = file(datei, 'r')
      text = f.read()               # load file into a string
      f.close()

      f = file(datei, 'w')          # save edited file

      f.write('/*\n')
      for tline in slic:
        f.write ( " * " + tline + " *\n" )
      f.write(' */\n\n')

      f.write(text)
      f.close()    

# the main program:
files = ["CMakeLists.txt",
         "CMakeLists.project.txt"]

for datei in files:
  f = file(datei, 'r')
  text = f.read()               # load file into a string
  f.close()

  f = file(datei, 'w')          # save edited file

  for tline in slic:
    f.write ( "# " + tline + "\n" )
  f.write ('\n')

  f.write(text)
  f.close()    

walk(getcwd(), addLicenseBlock, 4)
