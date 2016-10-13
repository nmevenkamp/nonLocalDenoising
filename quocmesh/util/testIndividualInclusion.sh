#!/bin/tcsh

# call this script from the quocmesh main directory

# For each header file found in the modules directory, this script
# writes one test program including the header (and only it).
# Moreover, three cpp files (two become only object files, one an
# executable) are generated.

mkdir selfTest/include

foreach i( `find modules -name "*.h"` )
	setenv FNAME `basename $i`
	setenv OFNAME selfTest/include/include_`basename $i .h`_h.cpp
	echo '#include<'$FNAME'>' >> selfTest/include/include_all.cpp
	echo '#include<'$FNAME'>' >> selfTest/include/include_all_one.cpp
	echo '#include<'$FNAME'>' >> selfTest/include/include_all_two.cpp
	echo '#include<'$FNAME'>' >> $OFNAME
	echo 'int main ( int, char** ) {' >> $OFNAME
	echo '  return 0;'  >> $OFNAME
	echo '}'  >> $OFNAME
	end

cd modules
foreach	i( *)
	foreach j($i/*.h)
		set FNROOT=`basename $j`
		set NOTFOUND=`grep -L $FNROOT ../selfTest/$i/selfTest.cpp`
		if ( $NOTFOUND != "" ) then 
			echo $j 'was not found in corresponding selfTest'
			exit(1)
		endif
	end
end
cd ..

echo 'int main ( int, char** ) {' >> selfTest/include/include_all.cpp
echo '  return 0;'  >> selfTest/include/include_all.cpp
echo '}'  >> selfTest/include/include_all.cpp

#echo ''
#echo 'now add selfTest/include to the PROJECTS section in your makefile.selection and run make'
#echo 'afterwards, remove selfTest/include from makefile.selection and delete that directory'   

cp makefile.selection makefile.selection.script_backup
echo 'PROJECTS += selfTest/include' >> makefile.selection
make projects
mv -f makefile.selection.script_backup makefile.selection
rm -rf selfTest/include
