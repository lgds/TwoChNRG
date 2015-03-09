#!/usr/bin/perl
##
##  Greps data from output*.txt in several directories and 
##  saves in 
##    OpAvg_T_param.dat in dir
## 
##    OpAvg_vsParam.dat in ./Runs
##
##
sub Print_Help(){
  print "Usage: \n"; 
  print "$0 -p choiceparam --grep=pattern -i nodei -f nodeF \n";
  print " (in development) \n";
  print "  Options: \n";
  print "   --chiloc  : calculates local suscep and saves in file SuscepImp_25_727.dat \n";
  print "   --indir   : run from within dir ./run(node) \n";
  print "   --Tzero   : also saves the T=0 value in file (pattern)_x_(dir).dat \n";
  print "  Gets params from FileParam(=input_nrg.dat by default):  \n";
  print "  use choiceparam as: \n";
  print "   General: \n";
  print "       0 - Nsites , 1- Nstates ; 2 - U ; 3 - Gamma; 4 - ed  \n";
  print "       5 - Lambda , 6- Dband  ; 7 - () ; 8 - () ; 9 - calcdens  \n";
  print "   CMphonon: \n";
  print "       11- w0 ; 12 - lambda ; 13 - alpha ; 14 - Nph \n";
  print "   SMM (and DQD up to #13): \n";
  print "       11 - U2;      12 - Gamma2;    13 - ed2 ; 14 - J12 \n";       
  print "       15 - BmagPar; 16 - BmagPerp;  17 - Dz ;  18 - B2  \n";       
#  print " Other Options: \n";
#  print "  --file_pattern FilePattern \n";
  exit(0);
}

use Getopt::Long;

##
##  Get Parameters
##
$InputFileName="input_nrg.dat";
if ( @ARGV > 0 ) {
GetOptions('file_pattern=s'=>\$FilePattern,
           'grep=s' => \$GrepPattern,
           'i|node0=i' => \$nodenum0,
           'f|nodeF=i' => \$nodenumF,
           'p|choiceparam=i'=>\$choiceparam,
           'chiloc'=>\$ChiLoc,
           'Tzero'=>\$Tzero,
           'indir'=>\$InDir,
           'h|help' => \$help);
}
else {Print_Help();}

if (defined($ChiLoc)){$GrepPattern="Szdot";}
if (!defined($GrepPattern)){Print_Help();}
if (defined($InDir)){
  if (!defined($nodenum0)){$nodenum0=1;}
  $nodenumF=$nodenum0;
}

$nodeStep=1;
if ($help){Print_Help();}
$param=-1;
$GrepCommand="cat output*.txt | grep $GrepPattern "; 

## Setting OpName: 
$OpName=$GrepPattern;
$OpName =~ s/=//;
##print " Command: $GrepCommand \n";


for ($nodenum=$nodenum0; $nodenum<= $nodenumF; $nodenum+=$nodeStep){

  if (!defined($InDir)){
    $nodedirname='./run'.$nodenum;
    chdir($nodedirname) or die "Could not change dir \n";
  }
  open(INFILE, $InputFileName) or die "Could not open the file $InputFileName \
n";
     @AllFile=<INFILE>; 
  close(INFILE);
  $param=@AllFile[$choiceparam];
  my $hfield=@AllFile[11];
## Needed for chiloc
  chop($param);
  chop($hfield);
##  print "Dir is $nodedirname ; param is $param ; field is $hfield \n";
  if (defined($ChiLoc)){
    $OutFileName="SuscepImp_25_727.dat";
  }
  else{
    $OutFileName=$OpName."_Temp_Param_".$param.".dat";
  }
##  print "File Name is $OutFileName \n";

  open (GREPDATA,"$GrepCommand |");
  open (OUTFILE,"> $OutFileName");
   while (<GREPDATA>){
     $TheLine=$_;
     chomp($TheLine);
#     if (defined($RemoveString)){$TheLine =~ s/$RemoveString//;}
     @OutCols=split(/\s+/,$TheLine);
##     print $TheLine."\n";
##     print OUTFILE $OutCols[$colx-1]." ".$OutCols[$coly-1]."\n";
##     print $OutCols[2]."   ".$OutCols[4]."\n";
     if (defined($ChiLoc)){
       my $TChiLoc=($OutCols[2])*(0.5*$OutCols[4])/$hfield;
       print OUTFILE $OutCols[2]."  ".$TChiLoc." ".$TChiLoc." ".$TChiLoc."\n";
     }else{
       print OUTFILE $OutCols[2]."   ".$OutCols[4]."\n";
     }
     ## end if ChiLoc
   }
## end while
  close (OUTFILE);
  close(GREPDATA);
  if (defined($Tzero)){
    $OutFileName=$OpName."_x_".$nodenum.".dat";
    open (OUTFILE,"> $OutFileName");
    print OUTFILE "$param   $OutCols[4] \n";
    close (OUTFILE);
  } else{print "$param   $OutCols[4] \n";}

  if (!defined($InDir)){chdir('../');}

}
