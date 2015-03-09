#!/usr/bin/perl
##
## Simple loop in z
##
###############################
sub Print_Help() {
  print " Usage $0 --z0=z0 --zF=zF --zstep=ztep\n"; 
  print " Options:  \n";
  print "    --band=BandType (\"SideDot\" is default, \"Const\" is an option \n";
#  print "    --dirF=N1: If N1 is not specified, N0:=N1  \n";
  exit;

}
###############################

use Getopt::Long;

my $z0=1.0;
my $zF=1.0;
my $zstep=0.0;
my $BandType="SideDot";

if ( @ARGV > 2 ){
GetOptions('z0=f'=>\$z0,
           'zF=f'=>\$zF,
           'zstep=f'=>\$zstep,
           'band=s'=>\$BandType,
           'h|help'=>\$help);
}
else{Print_Help();}
if (defined($help)){Print_Help();}

print "Z0 = $z0, ZF = $zF , Zstep = $zstep \n";


for ($zz=$z0;$zz<=$zF;$zz+=$zstep){
#  $CommandStr="nice ./NRG_main -b SideDot ";
  $CommandStr="nice ./NRG_main -b ".$BandType." ";
#  $OutputFile="output_SideDot_Anderson_zEQ".$zz.".txt";
  $OutputFile="output_".$BandType."_Anderson_zEQ".$zz.".txt";
  $CommandStr.="-z ".$zz." > ".$OutputFile." \& \n";
  print "Executing : \n".$CommandStr;
  system($CommandStr);
}
