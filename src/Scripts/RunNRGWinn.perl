#!/usr/bin/perl
##
## Script to start the NRG code in the Winnebago cluster 
##
##using getopt
##
##use Getopt::Std;
##%options=();
##getopts("",\%options);

sub Print_Help(){
  print "Usage $0 --dir0=d1 --dirF=d2 ";
  print "Options : \n";
  print "\t--help\t\t: prints this message \n";
  print "\t--queue=queuename \t: sends to queue \"queuename\"\n";
  exit;

}
use Getopt::Long;

GetOptions("queue=s"=>\$QueueName,
           "dir0=i" =>\$Dir0,
           "dirF=i" =>\$DirF,
           "help" => \$HelpOpt);

Print_Help() if defined $HelpOpt;
Print_Help() if ( (!defined($Dir0))||(!defined($DirF)) );

$test1="Running NRG ";
$test1.="in queue $QueueName " if defined $QueueName;
$test1.=" NOW! \n";
print $test1;

##
##  Choose code:
##
$ChoiceCode=1;
print "Choose code: \n";
print "      1     - OneChQS \n";
print "      2     - TwoChQS \n \n";
print "  choice : ";
chomp($ChoiceCode = <STDIN>);

for ($ChoiceCode){
 if (1){
  $CodeName="OneChQS";
 }
 elsif (2){
  $CodeName="TwoChQS";
 }
 else{print "Code not valid. Exiting... \n";exit; }
}

##
##  Choose model:
##

$ChoiceModel=0;
print "Choose model : \n";
print "      0     - Anderson \n";
print "      1     - Kondo \n";
print "      2     - Phonon (1ch) \n";
print "      3     - Chain \n";
print "      4     - CM phonon (2ch) \n";
print "  choice : ";
chomp($ChoiceModel = <STDIN>);

for ($ChoiceModel){
 if (0){
  $ModelName="Anderson";
 }
 elsif (1){
  $ModelName="Kondo";
 }
 elsif (2){
  $ModelName="Phonon";
 }
 elsif (3){
  $ModelName="Chain";
 }
 elsif (4){
  $ModelName="CMPhonon";
 }
 else{print "Model not valid. Exiting... \n";exit; }
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
   $Extension="NRG_Code_".$CodeName."_Model_".$ModelName."_Dir".$Dir;
   $ExecName=$CodeName." -m ".$ModelName;
###
### Preparing submission
###


# Build Script file
#

   $ScriptFileName="srun_script".$Extension;
   print "Building $ScriptFileName ... ";

   $String="\#!/bin/sh \n";
   $String.="\#SLURM -s\n\#\n\#SLURM -J job".$Extension."\n\#\n";
   $String.="\#SLURM -e job".$Extension.".error\n\#\n";
   $String.="\#SLURM -o job".$Extension.".output\n\#\n";
   $String.="\#SLURM -D ".$LocalDir."\n\#\n";
   $String.="./$ExecName > output_$Extension \n";
   $String.="echo Finished > $LocalDir/srun_finish_$Extension \n";

  open (OUT,"> $LocalDir/".$ScriptFileName);
    print OUT $String;
  close(OUT);

  print "... Done.\n";
#
# Run PBS
#

##  $String="srun -b $ScriptFileName \n";
  $String="srun -b ";
  $String.="-p $QueueName " if defined $QueueName;
  $String.="$ScriptFileName \n";
#
#
  print "Executing $String ...";

#  system($String);

}
## End Loop in Dirs
