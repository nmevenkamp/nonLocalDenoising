#!/usr/bin/perl -w

if($ARGV[0] =~ /--help/ || $ARGV[0] =~ /-h/ ) {
  print"Usage: within ERRPARSE pipe\n";
  print"Suppresses error messages or warnings listed in suppress.conf\n";
  die("\n");
}


$0 =~ s/\.pl$/.conf/;

open PREFS, $0 or die "Could not open configuration file $0, perhaps call me with absolute path?";

$expr = <PREFS>;
chomp $expr;
$expr = "($expr)";

while ($line = <PREFS>) {
    chomp $line;
    $expr .= "|($line)";
}

while (<>) {
    if (!m/$expr/) {
	print;
    }
}
