#!/usr/bin/perl
#
##
##  Move NRg data files from ./run$ii to ./Data 
##
require('GetParamsFromNRGin.perl');

sub Print_Help(){
  print "Usage: $0 -i dir0 -f dirN --datadir=\"Dir in Data\" \n";
  exit(0);
}

use Getopt::Long;

if (@ARGV > 0){
  GetOptions('h|help' => \$help,
             'i=i' => \$dir0,
             'f=i' => \$dirN,
             'datadir=s' => \$DataDir);

}else{Print_Help();}
if( (!defined($dir0))||
    (!defined($dirN))||
    (!defined($DataDir)) ){Print_Help();}

my $DirNo=30;
for ($DirNo=$dir0; $DirNo<=$dirN; $DirNo++){ 
  my $extension=&GetParamsFromNRGin($DirNo);
  $extension=~ s/_$//;

  print "$extension \n";
  if ($DirNo==$dir0){ ## Prompt the first time
    system("./MoveNRGresults.sh -d $DirNo -N $extension -D $DataDir");
  }
  else{
    system("./MoveNRGresults.sh -d $DirNo -N $extension -D $DataDir -p");
  }
}

