#!/usr/bin/perl

$threshold = 1.0;

if($ARGV [0] =~ /--help/ || $ARGV [0] =~ /-h/ ) {
  print"Usage: diff_l2.pl file1 file2\n";
  print"Computes l2 distance between two number sequences\n";
  print"read from two ASCII files\n";
  print"return true if dist < $threshold\n";
  die("\n");
}

$correctfile = $ARGV [0];
$testingfile = $ARGV [1];

open (CORRECT, $correctfile);
open (TESTING, $testingfile);

$difference = 0;

sub isanumber
{
    $string = shift @_;
    return $string =~ m/^[\+\-]?\d*[\.\d]\d*([eE][\+\-]?\d+)?$/;
}

while (!eof (CORRECT) || !eof (TESTING) || @corrlist || @testlist) {

    while (!@corrlist && !eof (CORRECT)) { @corrlist = split (' ', <CORRECT>); }
    while (!@testlist && !eof (TESTING)) { @testlist = split (' ', <TESTING>); }

    $corr = shift @corrlist; $test = shift @testlist;

    if (isanumber ($corr) && isanumber ($test)) {
	$delta = $corr - $test;
	$difference += $delta * $delta;
    }
    else { if ($corr ne $test) { die "Format mismatch.\n"; } }
}

$difference = sqrt ($difference);
print "Difference: $difference  ...  ";
if ($difference < $threshold) { print "Ok.\n"; }
else { print "Failed.\n"; die "\n"; }

    
