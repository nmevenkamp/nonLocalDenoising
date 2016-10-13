#!/bin/sh
TABLEFILE=`grep -l INSERT_PROJECT_TABLE_HERE doc/all/*.html`
perl -0777pe "s/\s*INSERT_PROJECT_TABLE_HERE.*//s" $TABLEFILE > projectTable.html
echo >> projectTable.html
echo >> projectTable.html
echo >> projectTable.html
`dirname $0`/tableOfProjects.sh `echo $@ | sed "s/[^ ]*selfTest[^ ]* *//g"` >> projectTable.html
echo >> projectTable.html
echo >> projectTable.html
echo >> projectTable.html
`dirname $0`/detailedTableOfProjects.sh `echo $@ | sed "s/[^ ]*\(selfTest\|projects\|finishedProjects\)[^ ]* *//g"` >> projectTable.html
perl -0777pe "s/.*INSERT_PROJECT_TABLE_HERE\s*//s" $TABLEFILE >> projectTable.html
cp projectTable.html $TABLEFILE
rm -f projectTable.html 
