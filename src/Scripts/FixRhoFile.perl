#!/usr/bin/perl
##
## Fix Rho Files and save it in the standard format
## Energy rho rho rho rho int
##
###############################
sub Print_Help() {
  print " Usage $0 --rhofile=RhoFilename \n"; 
#  print " Options:  \n";
#  print "    --chainfile=ChainFilename (SuscepChain.dat by default) \n";
 exit;

}
###############################

########################
##   Main program     ##
########################

use Getopt::Long;


if ( @ARGV > 0 ){
GetOptions('f|rhofile=s'=>\$RhoFilename,
          'h|help'=>\$help);
}
##else{Print_Help();}
if (!defined($RhoFilename)){Print_Help();}
if (defined($help)){Print_Help();}

if (-e $RhoFilename) {print "Found $RhoFilename . Keep going... \n";}
else{ print " $Rhofilename not found! Exiting... \n"; exit;}

$OutFilename=$RhoFilename;
##
##  Add an OLD 
##

$SansDatName=$RhoFilename;
$SansDatName =~ s/.dat//;
$OldFilename=$SansDatName."_OLD.dat";

print "Saving to $OldFilename \n ..";
$CpCommand="cp -p ".$RhoFilename." ".$OldFilename;
system($CpCommand);

################################################
##   Read Old File and Save to $RhoFilename   ##
################################################
my @LineCols=();
my $NewLine="";

 open(OUTFILE,"> ".$RhoFilename);
 open(INFILE,$OldFilename);
   while (<INFILE>){
     $TheLine = $_;
     chomp($TheLine);
     @LineCols=split(/ +/,$TheLine);
     $NewLine=
     $TheLine.=" 0 \n";
     printf(OUTFILE "%20.20e %20.20e %20.20e %20.20e %20.20e %i\n",$LineCols[0],$LineCols[1],$LineCols[1],$LineCols[1],$LineCols[1],0);
   }
 close(INFILE);
 close(OUTFILE);


print "...Fix Rho Done!\n";
#######
# END #
#######
