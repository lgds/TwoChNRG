#!/usr/bin/perl
##
##  Calculates conductance from NRG phase shifts
##
#Requires GrepParamFromFile.perl
$HomeDir=`echo \$HOME`;
chomp($HomeDir);
$ScriptsDir=$HomeDir."/DMRGtest/Runs/";
push(@INC,$ScriptsDir);
require('GrepParamFromFile.perl');

sub Print_Help(){
  print "Usage : $0 --file=[file pattern] --param=[param] -N [Nsite] \n";
  print "  [file pattern] are either :";
  print "      phaseshift_*.dat files (default) ";
  print "   OR output_* files (use --calcphaseshift then)";
  print " Optional: \n";
  print "  --save    : Saves to G_vs_[param].dat \n";
  print "  --nosave  : Does not save (default) \n";
  print "  --outname=[name] : Adds _[name] for output file name. \n";
  print "  --calcphaseshift : [file pattern] are output_* files and the code calculates everything from scratch.\n";
  exit(0);
}

use Getopt::Long;

if (@ARGV>0){
  GetOptions('h|help' => \$help,
             'f|file=s' => \$FilePattern,
             'param=s' => \$ParamName,
             'save!' => \$SaveFile,
             'outname=s' => \$AddName,
             'calciphaseshift' => \$CalcPhaseShift,
             'N=i' => \$Nsite);
}else{Print_Help();}

if (defined($help)){Print_Help();}
if ( (!defined($Nsite))||(!defined($FilePattern))||(!defined($ParamName)) ){
  print "Please define all options below: \n";
  Print_Help();
}

@MatchingFiles=glob($FilePattern);
if (defined($AddName)){$AddName="_".$AddName;}else{$AddName="";}
$OutputFilename="G_vs_".$ParamName.$AddName.".dat";
##if ($SaveFile){open (OUTFILE,"> G_vs_".$ParamName.$AddName.".dat");}
my $TheLine="";
if ($SaveFile){open (OUTFILE,"> ".$OutputFilename);}
foreach $FileName (@MatchingFiles){
  my $x_param=GetParamFromFile($FileName,$ParamName);
  if (defined($CalcPhaseShift)){
   $TheLine=`CalcPhaseShiftLevels -f $FileName --N0=$Nsite --Nf=$Nsite --E1chain=0.800048 --nosave | tail -n 1`;
  }else{
   $TheLine=`cat $FileName | grep \"^$Nsite \"`;
  }
  my @SplitLine=split(/ +/,$TheLine);
  my $ind=0;
  my $Cond=$SplitLine[4];
  if (defined($CalcPhaseShift)){
    foreach $Token (@SplitLine){
      if ($Token =~m/G:/){$Cond=$SplitLine[$ind+1];}
      else{$ind++;}
    }
  } #end if 
  print $x_param."  ".$Cond."\n";
  if ($SaveFile){ print OUTFILE $x_param."  ".$Cond."\n"; }
}
if ($SaveFile){close(OUTFILE);}

if ($SaveFile){
  system("sort -g $OutputFilename > G_test.dat; mv G_test.dat $OutputFilename"); 
}
# Done!
