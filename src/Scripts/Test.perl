## Removing leading spaces

 $String="   This is a test";
 print $String."\n";
 $String =~ s/^\s+//;
 print $String."\n";
