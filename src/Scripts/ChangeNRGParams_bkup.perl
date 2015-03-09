#!/usr/bin/perl

#
# Changes parameters in NRG input files in a loop 
# renames output file. 
#
# First lines of a typical lanc.in file:
#
# whichbandtype (0)
# e2            (1) 
# Delta2        (2)
# Delta1        (3)
# lambda        (4)
# MagFlux       (5)
#
########
sub Print_Help(){
  print "Usage: \n";
  print "$0 -i,--node0=ni -f,--nodeF=nf (--nodeStep=nstep) --choiceparam=(0,1,2,...) --file=(Filename)\n";
  print "Additional options: \n";
  print "  --file=(Filename) \n";
  print "  Filename=input_nrg.dat (default):  use choiceparam as: \n";
  print "   General: \n";
  print "       0 - Nsites , 1- Nstates ; 2 - U ; 3 - Gamma; 4 - ed  \n";
  print "       5 - Lambda , 6- Dband  ; 7 - () ; 8 - () ; 9 - calcdens  \n";
  print "   CMphonon: \n";
  print "       11- w0 ; 12 - lambda ; 13 - alpha ; 14 - Nph \n";
  print "   SMM (and DQD up to #13): \n";
  print "       11 - U2;      12 - Gamma2;    13 - ed2 ; 14 - J12 \n";       
  print "       15 - BmagPar; 16 - BmagPerp;  17 - Dz ;  18 - B2  \n";       
  print "  Filename=Input_Phonon.dat (older TwoChQS code). In this case, choiceparam is: \n";
  print "       0 - Nph , 1- w0 ; 2 - lambda ; 3 - tp; 4 - alpha \n";
  print "  --p0=p0 \n";
  print "  --pF=pF \n";
  print "  --pstep=pstep \n";
  print "  --askparam  : Ask for parameter value each time \n"; 
  print "  --justcheck : Just show the param and not change it \n"; 
#  print "  --runonly : Keep STMparams.in file (no change) \n"; 
#  print "  --rename  : Rename STM.dat file \n"; 
#  print "  --renamerho : All it does is to rename Rho_wNeven.dat and RhowBulla.dat files \n"; 

  exit(0);
}

##
## Main script
##

use Getopt::Long;

$nodeStep=1;
#$InputFileName="Input_Phonon.dat";
$InputFileName="input_nrg.dat";
if (@ARGV > 0){
  GetOptions('h|help' => \$help,
           'i|node0=i' => \$nodenum0,
           'f|nodeF=i' => \$nodenumF,
           'nodeStep=i' => \$nodeStep,
           'choiceparam=i' => \$choiceparam,
           'file=s' => \$InputFileName,
           'p0=f'=>\$param0,
           'pF=f'=>\$paramF,
           'pstep=f' => \$paramstep,
           'askparam' => \$AskParam,
           'justcheck' => \$JustCheck);
  if (($help)||(!$nodenum0)||(!$nodenumF)){Print_Help();}
}
else{Print_Help();}

 while ((!defined($choiceparam))||($choiceparam<0)||($choiceparam>20) ){
  print "Which parameter do you want to changei in $InputFileName? (0-20)? \n";
## Change later
  for ($InputFileName){
   if (/Input_phonon.dat/){
    print "  0 - Nph , 1- w0 ; 2 - lambda ; 3 - tp; 4 - alpha \n";
   } 
   elsif (/input_nrg.dat/){
    print "   General: \n";
    print "       0 - Nsites , 1- Nstates ; 2 - U ; 3 - Gamma; 4 - ed  \n";
    print "       5 - Lambda , 6- Dband  ; 7 - () ; 8 - () ; 9 - calcdens  \n";
    print "   CMphonon: \n";
    print "       11- w0 ; 12 - lambda ; 13 - alpha ; 14 - Nph \n";
    print "   SMM (and DQD up to #13): \n";
    print "       11 - U2;      12 - Gamma2;    13 - ed2 ; 14 - J12 \n";       
    print "       15 - BmagPar; 16 - BmagPerp;  17 - Dz ;  18 - B2  \n";
   }
  }       
  chomp($choiceparam = <STDIN>);
 }
 if (!defined($AskParam)){
   if (!defined($param0)){
    print "Enter the initial parameter value (ex: 0.01)? \n";
    chomp($param0 = <STDIN>);
   }
   if (!defined($paramstep)){
    print "Enter the step in parameter ? \n";
    chomp($paramstep = <STDIN>);
   }
 }
# end if askparam not defined
#
#
#
##
##  LOOP in the "node" directories
##

$param1=$param0;
for ($nodenum=$nodenum0; $nodenum<= $nodenumF; $nodenum+=$nodeStep){

  
  $nodedirname='./run'.$nodenum;
  chdir($nodedirname) or die "Could not change dir \n";
##
## Change Filename 
##
   open(INFILE, $InputFileName) or die "Could not open the file $InputFileName \n";
     @AllFile=<INFILE>; 
   close(INFILE);

   if (defined($AskParam)){
     print "Current value of parameter in node $nodenum is @AllFile[$choiceparam] ";
     print "Enter param value : \n";
     chomp($param1 = <STDIN>);
   }
   print "Dir $nodenum : Changing line $choiceparam in file $InputFileName from : @AllFile[$choiceparam] ";
   if(!defined($JustCheck)){@AllFile[$choiceparam]="$param1\n";}
   print "to : @AllFile[$choiceparam] \n";
   if(defined($JustCheck)){print "(would be to $param1)\n";}
##   $TestName="test.out";
##   open(OUTFILE, ">$TestName");
   open(OUTFILE, ">$InputFileName");
   print OUTFILE @AllFile;
   close(OUTFILE);


   chdir('../');
   $param1+=$paramstep;
   if (abs($param1)<1e-10){ $param1=0;}
## Round off at 10 places if not integer
   if ($param1!=int($param1)){$param1 = sprintf("%.10g",$param1);}   


 }
###END OF FOR LOOP
