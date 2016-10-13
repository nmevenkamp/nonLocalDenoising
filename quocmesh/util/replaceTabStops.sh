#!/bin/tcsh
foreach i(`find -name "*.c"` `find -name "*.cpp"` `find -name "*.h"` `find -name "*.hpp"` `find -name "*.py"` `find -name "*.cu"` )
    sed -i 's/\t/  /g' $i
end
