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
  print "$0 -i,--node0=ni -f,--nodeF=nf (--nodeStep=nstep) --choiceparam=(0,1,2,...) --p0=p0 --pF=pF   --pstep=pstep \n";
  print "Additional options: \n";
  print "  --file=(Filename) \n";
  print "  --askparam  : Ask for parameter value each time \n"; 
  print "  -j (--justcheck) : Just show the param and not change it \n";
  print "  --usefile=Filename (--propconst=A) : gets parameters from external file (ex DFT) and multiplies by A\n";
  print "  --indir : Runs from within directory \n"; 
  print "  --options : List all choiceparam options for all possible input files \n";
#  print "  --rename  : Rename STM.dat file \n"; 
#  print "  --renamerho : All it does is to rename Rho_wNeven.dat and RhowBulla.dat files \n"; 

  exit(0);
}


sub ListOptions(){
  print " All choiceparam Options \n"; 
  print "  Filename=input_nrg.dat (default):  use choiceparam as: \n";
  print "   General: \n";
  print "       0 - Nsites , 1- Nstates ; 2 - U ; 3 - Gamma; 4 - ed  \n";
  print "       5 - Lambda , 6- ChemPot ; 7 - () ; 8 - UseCFS ; 9 - calcdens  \n";
  print "   CMphonon: \n";
  print "       11- w0 ; 12 - lambda ; 13 - alpha ; 14 - Nph \n";
  print "   SMM (and DQD up to #13): \n";
  print "       11 - U2;      12 - Gamma2;    13 - ed2 ; 14 - J12 \n";       
  print "       15 - BmagPar; 16 - BmagPerp;  17 - Dz ;  18 - B2  \n";       
  print "   Majorana : \n";
  print "       11 - Bmag ; 12 - t1; 13 - t2 ; 14 - phi_mag; 15 - em  \n";     
  print "  Filename=Input_Phonon.dat (older TwoChQS code). In this case, choiceparam is: \n";
  print "       0 - Nph , 1- w0 ; 2 - lambda ; 3 - tp; 4 - alpha \n";
  print "  Filename=lanc.in :  use choiceparam as: \n";
  print "       0 - whichbandtype , 1 - e2, 2 - Gamma2, 3 - Gamma1, 4 - lambda, 5 - MagFlux \n";
  print "       (Pseudogap) G(w)=G0*|w-w0|^r + gamma : 1 - r , 2 - G0, 3 - gamma , 4 - w0 \n";
  print "       whichbandtype: 1 -> Pseudogap 4 -> SideDot, 5 -> Parallel Dot (see lanc.in for a description) \n";

  exit(0);
}

##
## Main script
##

use Getopt::Long;
require('GetYgivenX.perl');

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
           'j|justcheck' => \$JustCheck,
           'usefile=s'=>\$UseFilename,
           'indir'=>\$InDir,
           'propconst=f'=>\$PropConst,
           'options'=>\$ListOptions);
  if ($ListOptions){ListOptions();}
  if (($help)||(!$nodenum0)||(!$nodenumF)){Print_Help();}
}
else{Print_Help();}

 while ((!defined($choiceparam))||($choiceparam<0)||($choiceparam>20) ){
  print "Which parameter do you want to change in $InputFileName? (0-20)? \n";
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
 if ( (defined($UseFilename))&&(!defined($PropConst)) ){
   $PropConst=1;
   print "Parameters from $UseFilename: \n";
   print "NRGparam = A * F(param) (F(X) defined in file)\n"; 
   print " Enter A: "; 
   chomp($PropConst = <STDIN>);
 }
# end if askparam not defined
#
#
#
##
##  LOOP in the "node" directories
##

if(defined($InDir)){$nodenumF=$nodenum0;}
$param1=$param0;
for ($nodenum=$nodenum0; $nodenum<= $nodenumF; $nodenum+=$nodeStep){

  
  $nodedirname='./run'.$nodenum;
  if(!defined($InDir)){chdir($nodedirname) or die "Could not change dir \n";}
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
   my $WouldBeParam=$param1;
    
   my $rhoJ=0;
   if(defined($UseFilename)){
    my $nDFT=GetYgivenX($UseFilename,$param1);
#### Linear
##     $rhoJ=$PropConst*$nDFT;
#### Sqrt
##     $rhoJ=$PropConst*sqrt($nDFT);
#### Log: -1/(log()+log(B/D))
#    $rhoJ=-1.0/(log($nDFT)+log($PropConst));
#### rhoJ0+B*nDFT
    my $rhoJ0=0.15;
    $rhoJ=$rhoJ0+$PropConst*$nDFT;
#    print "nDFT = $nDFT, B=$PropConst, rhoJ = $rhoJ , J = ".2.0*$rhoJ." log(0.15)=".log(0.15)." \n"; 
    if(!defined($JustCheck)){@AllFile[$choiceparam]=2.0*$rhoJ."\n";}
    else{$WouldBeParam=2.0*$rhoJ;}
   }
   print "to : @AllFile[$choiceparam] \n";
   if(defined($JustCheck)){print "(would be to $WouldBeParam)\n";}
##   $TestName="test.out";
##   open(OUTFILE, ">$TestName");
   open(OUTFILE, ">$InputFileName");
   print OUTFILE @AllFile;
   close(OUTFILE);


   if(!defined($InDir)){chdir('../');}
   $param1+=$paramstep;
   if (abs($param1)<1e-10){ $param1=0;}
## Round off at 10 places if not integer
   if ($param1!=int($param1)){$param1 = sprintf("%.10g",$param1);}   


 }
###END OF FOR LOOP
