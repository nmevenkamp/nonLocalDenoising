#!/usr/bin/python

# Skript that checks, if source directories that would be included in the opensource release are either
# listed in cmake.selection.opensource or explicitly excluded

import argparse
import re
import os

failed = False

# Set up command line parser
parser = argparse.ArgumentParser(description='Check code in open source directories against cmake open source selection and exclusion files.')
parser.add_argument('directories', metavar='directories', nargs='+', help='The directories that should be checked for code')
parser.add_argument('--selection', metavar='path', default='cmake.selection.opensource', help='The location of the open source selection file')
parser.add_argument('--exclusion', metavar='path', default='exclusion.opensource', help='The location of the open source exclusion file')

args = parser.parse_args()

# Open and read selection and exclusion lists
selectionfile = open(args.selection, 'r')
selection = selectionfile.read()

exclusionfile = open(args.exclusion, 'r')
exclusion = exclusionfile.read()

# Get the directories to check (subdirectories will be included)
dirs = args.directories

# Regex pattern that matches paths below directories in dirs
pathpattern = re.compile(r'(?:' + '|'.join(dirs) + ')' + '/\S+')
# Regex pattern that matches paths below directories in dirs, only at the start of a line
blpathpattern = re.compile(r'^(?:' + '|'.join(dirs) + ')' + '/\S+', re.MULTILINE)

# Remove comments in selection and exclusion
selection = re.sub(r'#.*', '', selection)
exclusion = re.sub(r'#.*', '', exclusion)

# Extract paths from lists
selection = pathpattern.findall(selection)
exclusion = blpathpattern.findall(exclusion)

# Pattern that matches source files
sfpattern = re.compile(r'\S+\.(c|cpp|cu|cl|h|hpp)$')

# Walk dirs
for d in dirs:
    for root, dirs, files in os.walk(d):
        # If a source file exists in a subdirectory
        if any([sfpattern.match(f) for f in files]):
            # Windows uses a different directory separator. Replace it to match the one used in the exclusion and selection files.
            root = root.replace("\\","/")
            # Check, if the directory has been included or excluded
            if not ((root in selection) or (root in exclusion)):
                # If not, print an error
                print root, 'is NOT listed in', args.selection, 'or', args.exclusion + ', but contains source files!'
                failed = True

# Exit with correct exit code
if failed:
    exit(1)

exit(0)
