#!/usr/bin/python

# replaceExpressions.py
# script traverses the current directory tree and calls 'replaceExpressions'
# on every .h or .cpp-file. If one of the expressions that are specified here in
# the source code is found, it will be replaced.

from os import *
from os.path import *
from re import *

# replace one or more expressions
def replaceExpressions(arg, drctry, files):
  chdir(drctry)
  drctry = drctry.split("quocmesh")
  drctry = drctry[-1]
  for datei in files:
    if isfile(datei) and ( (datei.lower()).endswith('.cpp') or (datei.lower()).endswith('.h')):
      print "Testing '" + drctry + "/" + datei + "'...                                 \r",
      f = file(datei, 'r')          
      text = f.read()               # load file into a string
      f.close()
      # ATTENTION: Don't change and commit this file. To avoid an accidental commit
      # you can for example copy this script to your quocmesh-directory and only edit
      # this local copy.
      textNeu = text.replace(r"oldExpression1",r"")
      textNeu = textNeu.replace(r"oldExpression2",r"newExpression2")
      #textNeu = textNeu.replace(r"oldExpression3",r"newExpression3")
     
      if textNeu != text:     
        print "Modified file '" + datei + "'...                                              "
        f = file(datei, 'w')          # save edited file
        f.write(textNeu)
        f.close()    
        #uncomment the following lines if you want backup-files
        #f = file(datei+".orig", 'w')          # save original file
        #f.write(text)
        #f.close()    

# the main program: start the walk
walk(getcwd(), replaceExpressions, 4)


print "                                                                          \r", 
print "\nready"
