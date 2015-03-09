#!/usr/bin/perl
##
## Script to start the NRG code in the UT/Winnbago clusters 
##
##using getopt
##
##use Getopt::Std;
##%options=();
##getopts("",\%options);

sub Print_Help(){
  print "Usage $0 -i=d1(--dir0=d1)  -f=d2(--dirF=d2)  (-C,--command=\"Command -m Option\") \n";
  print "Options : \n";
  print "\t--help\t\t: prints this message \n";
  print "\t--addtoname=string \t: adds \"string\" to the output file names extension\n";
  print "\t--queue=queuename \t: sends to queue \"queuename\"\n";
  exit;

}
use Getopt::Long;

GetOptions("queue=s"=>\$QueueName,
           "i|dir0=i" =>\$Dir0,
           "f|dirF=i" =>\$DirF,
           "C|command=s" => \$Command,
           "addtoname=s" => \$AddToName,
           "h|help" => \$HelpOpt);

Print_Help() if defined $HelpOpt;
Print_Help() if ( (!defined($Dir0))||(!defined($DirF)) );

##
## Check hostname
##
                                                                              
   chomp($HostName=`/bin/hostname --long`);
   for ($HostName)
    {
      if (/winnebago/){print "Winnebago cluster \n"; }
      elsif (/correlated/){print "Correlated cluster. \n";}
      elsif (/wilson/){print "Wilson cluster. \n";}
      elsif (/glyph/){print "Glyph cluster. \n";}
      else{print "Cant recognize cluster. Exiting... \n"; exit(0);}
    }
                                                                              

$test1="Running NRG ";
$test1.="in cluster $HostName ";
$test1.="in queue $QueueName " if defined $QueueName;
$test1.="in dirs $Dir0 through $DirF ";
##$test1.=" and avoiding $AvoidNode" if defined $AvoidNode;
$test1.=" NOW! \n";
print $test1;

##
##  Choose code:
##

my $UseAllCodes=0;
if (!defined($Command)){
  $ChoiceCode=3;
  print "Choose code: \n";
  print "      1     - OneChQS \n";
  print "      2     - TwoChQS \n";
  print "      3     - NRG_main (default) \n";
  print "        3.1 - NRG_main with ./DM_NRG and ./Conductance \n";
  print "      4     - DM_NRG \n";
  print "        4.1 - DM_NRG with ./Conductance \n";
  print "      5     - Conductance \n \n";
  print "  choice : ";
##  chomp($ChoiceCode = <STDIN>);
  chomp($test1 = <STDIN>);
  if ( ($test1 ne '')&&($test1>0)&&($test1<=5) ){
    $ChoiceCode=$test1;
  }else{print "Invalid choice. Using default (=$ChoiceCode) \n";}

  for ($ChoiceCode){
   if (/^1/){
    $CodeName="OneChQS";
   }
   elsif (/^2/){
    $CodeName="TwoChQS";
   }
   elsif (/^3/){
    $CodeName="NRG_main";
    if ($ChoiceCode=~ m/^3.1/){
        $ChoiceCode=3;
	$UseAllCodes=1;
	print "Using ALL codes...\n";
    }
   }
   elsif (/^4/){
    $CodeName="DM_NRG";
    if ($ChoiceCode=~ m/^4.1/){
        $ChoiceCode=4;
	$UseAllCodes=1;
	print "Adding Conductance code...\n";
    }
   }
   elsif (/^5/){
    $CodeName="Conductance";
   }
   else{print "Code not valid. Exiting... \n";exit; }
  }
  print " --- Code = $CodeName \n";
##
##  Choose model:
##
  if (($ChoiceCode==4)||($ChoiceCode==5)){
    $ChoiceModel=7;
  }
  else{
    $ChoiceModel=0;
    print "Choose model : \n";
    print "      0     - Anderson QS basis (default) \n";
    print "        0.1 - Anderson Q basis (NRG_main only) \n";
    print "        0.2 - Anderson (Q Sz) basis (NRG_main only) \n";
    print "      1     - Kondo \n";
    print "      2     - Phonon (1ch) \n";
    print "      3     - Chain \n";
    print "      4     - CM phonon (2ch) \n";
    print "        4.1 - CM phonon w/ parity (QSP basis, NRG_main only) \n";
    print "      5     - SMM \n";
    print "      6     - Double Quantum Dot (DQD) 1chQS \n";
    print "        6.1 - DQD with Zeeman in both dots (1chQSz) (NRG_main only) \n";
    print "      7     - DM_NRG/Conductance calculation (no model) \n";
    print "  choice : ";
#    chomp($ChoiceModel = <STDIN>);
    chomp($test1 = <STDIN>);
    if ($test1 ne ''){
      $ChoiceModel=$test1;
    }else{print "Using default (=$ChoiceModel) \n";}
  }
  for ($ChoiceModel){
   if (/^0/){
    $ModelName="Anderson";
    if ($ChoiceModel =~ m/^0.1/){$ModelName="1chQ_Anderson";}
    if ($ChoiceModel =~ m/^0.2/){$ModelName="1chQSz_Anderson";}
   }
   elsif (/^1/){
    $ModelName="Kondo";
   }
   elsif (/^2/){
    $ModelName="Phonon";
   }
   elsif (/^3/){
    $ModelName="Chain";
   }
   elsif (/^4/){
    $ModelName="CMphonon";
    if ($ChoiceModel =~ m/^4.1/){$ModelName="2chQSP_CMphonon";}
   }
   elsif (/^5/){
    $ModelName="SMM";
   }
   elsif (/^6/){
    $ModelName="1chQS_DQD";
    if ($ChoiceModel =~ m/^6.1/){$ModelName="1chQSz_DQD";}
   }
   elsif (/^7/){
    $ModelName="";
   }
     else{print "Model not valid. Exiting... \n";exit; }
  }
 
  print " --- Model: $ModelName \n";
##
##  Choose band :
##
  if ($ChoiceCode==4){ 
    $ChoiceBand=0;
    $Mtemp=1000;
    $betabar=0.727;
    $FiniteT=0;
    $NewBBar=0;
    $LoopMtemp=0;
    print "Mtemp = 1000 (T=0) and betabar=0.727. Change? (y/n) ";
    chomp($ChoiceYN = <STDIN>);
    if ($ChoiceYN eq "y"){
      print "Mtemp = "; chomp($Input = <STDIN> );
      if ( ($Input ne '') && ($Input>0) && ($Input<1000) ){
        $Mtemp=$Input;
        $FiniteT=1;
        print "Mtemp = $Mtemp \n";
        print "Loop Mtemp in each dir? (Mtemp0=$Mtemp ) (y/n) ";
        chomp($ChoiceLoopMtemp = <STDIN> );
        if ($ChoiceLoopMtemp eq "y"){
          $LoopMtemp=1;
##          print "MtempFinal = "; chomp($MtempFinal = <STDIN> );
          print "MtempStep = "; chomp($MtempStep = <STDIN> );
        } ## end if LoopMtemp
      } else {print "Not valid. Keeping Mtemp = $Mtemp;"}
      print "(default value: 0.727) betabar = "; chomp($Input = <STDIN> );
      if ( ($Input ne '') && ($Input>0) && ($Input<1.0) ){
        $betabar=$Input;
        $NewBBar=1;
        print "betabar = $betabar \n";
      } else {print "Not valid. Keeping betabar = $betabar\n";}
    }## end if choice=y
  } ## end if ChoiceCode=DM_NRG
  elsif ($ChoiceCode==5){
    $ChoiceBand=0;
  }## ChoiceCode=Conductance 
  else{
    $ChoiceBand=0;
    print "Choose host band type : \n";
    print "      0     - Constant (square) DoS [rho(e)=rho_0 -D < e < D] (default) \n";
    print "      1     - Side Dot (needs lanc.in) \n";
    print "      2     - Const (to use z-trick, for instance) \n";
    print "      3     - PowerLaw \n";
    print "      4     - FromFile (not implemented yet) \n";
    print "  choice : ";
#    chomp($ChoiceBand = <STDIN>);
    chomp($test1 = <STDIN>);
    if ($test1 ne ''){
      $ChoiceBand=$test1;
    }else{print "Using default (=$ChoiceBand) \n";}
  }
  for ($ChoiceBand){
   if (/^0/){
    $BandName="SquareWilson";
   }
   elsif (/^1/){
    $BandName="SideDot";
   }
   elsif (/^2/){
    $BandName="Const";
   }
   elsif (/^3/){
    $BandName="PowerLaw";
   }
   else{print "Band not implemented. Using Default \n";
     $ChoiceBand=0;
     $BandName="SquareWilson";
   }
  }
# end for ChoiceBand
##
##  Ztrick loop? 
##
  print "Use Ztrick (z=0.75,1,1.25,1.5) ? (y/n) :";
  chomp($UseZtrickYN = <STDIN>);
  if ($UseZtrickYN eq "y"){$UseZtrick=1;}else{$UseZtrick=0;}
}
# end if !defined($Command)
else{
## Strip "./" from Command
 $Command =~ s/\.//;
 $Command =~ s/\///;
 print "Code/Option is $Command \n";
}

##
## Loop on dirs 
##
##
## Executing
##

for ($Dir=$Dir0;
    $Dir<=$DirF;
    $Dir+=1){

   $TwoChDir=$ENV{'TWOCHDIR'};
   $LocalDir=$TwoChDir."/Runs/run$Dir";
   print "Root Dir: $TwoChDir\n"; 
##   system("cd $TwoChDir/Runs/run$Dir/;ls");
##   chomp($LocalDir=`/bin/pwd`);
   print "Running from : $LocalDir \n";
   $User=$ENV{'USER'};
##   $Extension="NRG_Code_".$CodeName."_Model_".$ModelName."_Dir".$Dir;
   $Extension="NRG_Code_".$CodeName;
   if ($ModelName ne ""){$Extension.="_Model_".$ModelName;}
   if ($ChoiceBand != 0){$Extension.="_Band_".$BandName;}
   $Extension.="_Dir".$Dir;
   if(defined($AddToName)){$Extension.=$AddToName;}
   if(defined($Command)){$ExecName=$Command;}
   else{
     $ExecName=$CodeName;
     if ($ChoiceCode==4){
       if ($FiniteT==1){
         $ExecName.=" -M ".$Mtemp;
         if ($LoopMtemp==1){print "Running Mtemp = $Mtemp \n";$Mtemp+=$MtempStep;}
       }
       if ($NewBBar==1){$ExecName.=" -b ".$betabar;}
     } elsif ($ChoiceCode!=5) {
##     } else {
       if ($ModelName ne ""){$ExecName.=" -m ".$ModelName;}
       if ($ChoiceBand != 0){$ExecName.=" -b ".$BandName;}
     } # end if ChoiceCode=4
   }
###
### Preparing submission
###

# Build Script file
#
   $ScriptFileName="sub_script".$Extension;
   $String="\#!/bin/sh \n";

   if ($HostName =~ m/winnebago/)
     {
      $ScriptFileName="SLURM_".$ScriptFileName;
      $String.="\#SLURM -s\n\#\n\#SLURM -J job".$Extension."\n\#\n";
      $String.="\#SLURM -e $LocalDir/job".$Extension.".error\n\#\n";
      $String.="\#SLURM -o $LocalDir/job".$Extension.".output\n\#\n";
##      $String.="\#SLURM -D ".$LocalDir."\n\#\n";
##      $SubCommand="srun -b ";
      $SubCommand="sbatch ";
      $SubCommand.="-p $QueueName " if defined $QueueName;
      $SubCommand.=" $LocalDir/$ScriptFileName \n";
     }
   elsif ($HostName =~ m/correlated/)
     {
      $ScriptFileName="LAVA_".$ScriptFileName;
      $String.="\#BSUB -J job".$Extension."\n\#\n";
      $String.="\#BSUB -e job".$Extension.".error\n\#\n";
      $String.="\#BSUB -o job".$Extension.".output\n\#\n";
##   $String.="\#BSUB -D ".$LocalDir."\n\#\n";
      $SubCommand="bsub ";
      $SubCommand.="-q $QueueName " if defined $QueueName;
      $SubCommand.=" < $LocalDir/$ScriptFileName \n";
     }
   elsif ( ($HostName =~ m/wilson/)||($HostName =~ m/glyph/) ) 
     {
      $ScriptFileName="TORQUE_".$ScriptFileName;
      $String.="\#PBS -N job".$Extension."\n\#\n";
      $String.="\#PBS -e job".$Extension.".error\n\#\n";
      $String.="\#PBS -o job".$Extension.".output\n\#\n";
##   $String.="\#PBS -d ".$LocalDir."\n\#\n";
      $SubCommand="qsub ";
      $SubCommand.="-q $QueueName " if defined $QueueName;
      $SubCommand.=" < $LocalDir/$ScriptFileName \n";
     }
   else {print "Which cluster is this? \n"; exit(0);}
   $String.="cd ".$LocalDir."\n";
## Adding z-trick
   if ($UseZtrick==1){
     $z0=0.75;
     $zF=1.51;
     $zstep=0.25;
     for ($zz=$z0;$zz<=$zF;$zz+=$zstep){
       $ExecNameZ=$ExecName." -z ".$zz;
       $ExtensionZ=$Extension."_zEQ".$zz;
       $String.="./$ExecNameZ > output_$ExtensionZ.txt \n";
     }
     # end loop in Z
   }
   else{
    $String.="./$ExecName > output_$Extension.txt \n";
   }

   print "UseAllCodes = $UseAllCodes \n" ;
   if ($UseAllCodes==1){
      print "Adding DM_NRG and/or Conductance \n";
      if($ChoiceCode != 4){$String.="./DM_NRG > output_DMNRG_Dir$Dir.txt \n";}
      if($ChoiceCode != 5){$String.="./Conductance > output_Conductance_Dir$Dir.txt \n";}
   }
   
   $String.="echo Finished > $LocalDir/Finish_$Extension \n";
   print "Building $ScriptFileName ... ";

  open (OUT,"> $LocalDir/".$ScriptFileName);
    print OUT $String;
  close(OUT);

  print "... Done.\n";
#
# Run PB/SLURM...
#
  print "Executing $SubCommand ...\n";
  system($SubCommand);
  print "... Done!\n";

}
## End Loop in Dirs
