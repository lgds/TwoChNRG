#!/usr/bin/perl
##
## Read HybFunc.dat in dir and extract 
## Gamma0 and gamma from
##  Gamma(w)=Gamma_0|w|+gamma
##
###############################

##
sub Print_Help() {
  print " Usage $0 -i dir0 -f dirN\n"; 
  print " Fits HybFuc.dat around w=0 with A|w|+B \n";
  print " Output cols: idir  A  B  A+  A-\n";

#  print "    --shiftchain=nshift : Shift chain values by nshift(>0). Does not read chain file.  \n";
  exit;

}
###############################



########################
##   Main program     ##
########################

use Getopt::Long;


if ( @ARGV > 0 ){
GetOptions('i=i'=>\$dir0,
           'f=i'=>\$dirF,
           'h|help'=>\$help);
}
##else{Print_Help();}
if (defined($help)){Print_Help();}
if ( (!defined($dir0))||(!defined($dirF)) ){Print_Help();}


##
## Loop in dirs
##

my $idir=$dir0;
my $HybFuncFileName="HybFunc.dat";
my $OutFileName="Gamma0_gamma.dat";

my @Gamma0Array=();
my @gammaArray=();


open(OUTFILE,"> ".$OutFileName);

while ($idir<=$dirF){

  my $nodedirname='./run'.$idir;
  chdir($nodedirname) or die "Could not change dir to $nodedirname \n";
  my $auxEn=-10000.0;
  my $auxGamma=0.0;

  my $Gamma0=0.0;
  my $Gamma0pos=0.0;
  my $Gamma0neg=0.0;
  my $smallgamma=0.0;
  my $FoundZero=0;

  open(INFILE,$HybFuncFileName);
  while (<INFILE>){
    my $oldEn=$auxEn;
    my $oldGamma=$auxGamma;
     $TheLine = $_;
    chomp($TheLine);
    ($auxEn,$auxGamma)=split(/ +/,$TheLine);
    if ($FoundZero==1){ ## Right after zero
      $Gamma0pos=($auxGamma-$smallgamma)/abs($auxEn);
      $Gamma0=0.5*($Gamma0pos+$Gamma0neg);
#      print "Gamma0pos = $Gamma0pos, Gamma0neg = $Gamma0neg \n";
      $FoundZero=0;
    }
    if (abs($auxEn)<1e-3){
      $FoundZero=1;
      $smallgamma=$auxGamma;
      $Gamma0neg=($oldGamma-$smallgamma)/abs($oldEn);
    }
    ## if found zero

  }
  close(INFILE);
  
  print "dir: $idir - Gamma0 = $Gamma0, gamma = $smallgamma \n";

  print OUTFILE $idir."   ".$Gamma0." ".$smallgamma."  ".$Gamma0neg."  ".$Gamma0pos."\n";


#  push(@gammaArray,$smallgamma);
#  push(@Gamma0Array,$Gamma0);

  chdir('../');
  $idir++;

}
## end loop in dirs

close(OUTFILE);
