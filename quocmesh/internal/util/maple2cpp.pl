#!/usr/bin/perl -w
# script for converting maples "C" matrix output into block format
# e.g.: cat mat | maple2cpp.pl where mat countains the pasted maple output

while ( <> ) {
  chomp;
  $line .= $_;
}

@entries = split( /;/, $line );

print "{";
$k = 0;
for ( $i = 0; $i < 8; $i++ ) {
  print " { ";
  for ( $j = 0; $j < 8; $j++ ) {
    $entries[ $k++ ] =~ /=(.+)$/;
    print $1; 
    if ( $j < 7 ) {
      print ", ";
    }
  }
  if ( $i < 7 ) {
    print " },\n ";
  } else {
    print " }};";
  }
}

