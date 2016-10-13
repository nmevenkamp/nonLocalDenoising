#!/usr/bin/python
import sys

for fname in sys.argv[1:]:
    infile = open ( fname, mode='r' )
    text = infile.read();
    infile.close()
    
    infile = open ( fname, mode='r' )     # open file again
    lines = []
    changemade = 0
    for line in infile.readlines():
        lines.append( line.rstrip() ) # rstrip removes trailing whitespace
        if ( lines[-1]+'\n' != line ):
            changemade = 1
    infile.close()

    if ( changemade == 1 ):
        print "Modified file '" + fname + "'...                                           "
        backupfile = open ( fname+".bak", mode='w' )      # save backup
        backupfile.write ( text )
        backupfile.close()

        outfile = open ( fname, mode='w' )    # save to original filename
        outfile.seek ( 0 )
        for line in lines:
            outfile.write ( line+'\n' )
        outfile.close()
