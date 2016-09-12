#!/usr/bin/perl
##
## Read Thermo files in same dir (or different dirs) 
## and compute an average.
##
##   Interpolation: use the zEQ1 file as the "starting" point
##
##

###############################

sub ReadSusFile{

my $TempArray=shift(@_);
my $SusArray=shift(@_);
my $SusChainArray=shift(@_);
my $SusImpArray=shift(@_);
my $FileName=shift(@_);
my $nt=0;
my $LowTemp=10;

@$TempArray=();
@$SusArray=();
@$SusChainArray=();
@$SusImpArray=();

 open(INFILE,$FileName);
   while (<INFILE>){
     $TheLine = $_;
     chomp($TheLine);
     ($auxtemp,$auxsuschainpimp,$auxsuschain,$auxsusimp)=split(/ +/,$TheLine);
     if (abs($auxtemp)<$LowTemp){$LowTemp=abs($auxtemp);}
     push(@$TempArray,$auxtemp);
     push(@$SusArray,$auxsuschainpimp);
     push(@$SusChainArray,$auxsuschain);
     push(@$SusImpArray,$auxsusimp);
     $nt++;
   }
 close(INFILE);

 return(($nt,$LowTemp));
}

###############################

sub SaveSusFile{

my $NoArgs=@_;
my $OddEven="n";

my $TempArray=shift(@_);
my $SusArray=shift(@_);
my $SusChainArray=shift(@_);
my $SusImpArray=shift(@_);
my $FileName=shift(@_);
if ($NoArgs==6){$OddEven=shift(@_);}
my $nt=@$TempArray;

my $i0=0;
my $istep=1;

for ($OddEven){
 if (/o/){$i0=1;$istep=2;}
 elsif(/e/){$i0=0;$istep=2;}
}

 open(OUTFILE,"> ".$FileName);
#  for ($ii=0;$ii<$nt;$ii++){
  for ($ii=$i0;$ii<$nt;$ii+=$istep){
    printf(OUTFILE "%20.20e %20.20e %20.20e %20.20e\n",@$TempArray[$ii],@$SusArray[$ii],@$SusChainArray[$ii],@$SusImpArray[$ii]);
  }
 close(OUTFILE);

 return;
}
# End Save SusFile

##########
##########

sub Interpol{

## Input: X_array, Y_array, Xvalue
## Output: Yvalue
##
## Remember: Assumes that X_array is in DECREASING order
##
## OPS! Can't assume anything.

  my $XArray=shift(@_);
  my $YArray=shift(@_);
  my $xval=shift(@_);

  my $A1=0.0;
  my $B1=0.0;
  my $yval=0.0;

## How to get this one??

  my $ArraySize=@$XArray;

  my $DecreasingXarray=1;
  my $MaxX=@$XArray[0];
  my $YMaxX=@$YArray[0];
  my $MinX=@$XArray[$ArraySize-1];
  my $YMinX=@$YArray[$ArraySize-1];
  ## assumes X0>X_1>...>XN

## Check if X0<X_1<...<XN instead
  if (@$XArray[0]<@$XArray[1]){
   ## Detected increasing X values (spectral functions...)
   $DecreasingXarray=0;
   $MaxX=@$XArray[$ArraySize-1];
   $YMaxX=@$YArray[$ArraySize-1];
   $MinX=@$XArray[0];
   $YMinX=@$YArray[0];
  }
  ##

## Local $xval in $XArray



## Debugging
  my $PrintNow=0;
  ##if ((abs($xval+1.66e-2)<0.0001)||(abs($xval-1.66e-2)<0.001)){$PrintNow=1;}

  if ($xval>=$MaxX){return($YMaxX);}
  if ($xval<=$MinX){return($YMinX);}

  my $iix=0;
## Looks for iix such that x_iix<= xval < x_iix+1
## Watch out for the EQUAL sign!
  if ($DecreasingXarray){
    if ($PrintNow){print "Decreasing X array \n";}
    while ( ($xval<@$XArray[$iix])&&((abs((@$XArray[$iix]-$xval)/@$XArray[$iix]))>1e-8)&&($iix<$ArraySize) ){$iix++;}
  }else{
    if ($PrintNow){print "Increasing X array \n";}
    while ( ($xval>@$XArray[$iix])&&((abs((@$XArray[$iix]-$xval)/@$XArray[$iix]))>1e-8)&&($iix<$ArraySize) ){$iix++;}
    $iix--; ## Needed in this case
  }
## end if Decreasing
  if ((abs(@$XArray[$iix]-@$XArray[$iix+1])/abs(@$XArray[$iix]))<1e-8){$iix++;}

  if ($PrintNow){
    print "iix = $iix , |x0-x1|/x0=".(abs((@$XArray[$iix]-@$XArray[$iix+1])/@$XArray[$iix]))." \n";
    print "x(".$iix.")= ".@$XArray[$iix]." , x(".($iix+1).")= ".@$XArray[$iix+1]." \n";
  }
  ## end PrintNow

  $A1=((@$YArray[$iix]-@$YArray[$iix+1]))/((@$XArray[$iix]-@$XArray[$iix+1]));

  $B1=@$YArray[$iix]-$A1*@$XArray[$iix];

  $yval=$A1*$xval+$B1;

  if ($PrintNow){
    print "y(".$iix.")= ".@$YArray[$iix]." , y(".($iix+1).")= ".@$YArray[$iix+1]." \n";
    print "A1 = ".$A1." B1 = ".$B1."\n";
    print "xval = $xval, yval = $yval \n";
  }

  return($yval);
}
## End Interpol

###############################
sub Print_Help() {
  print " Usage $0 (--type=String)\n";
  print " Options:  \n";
  print "    --type=String (--type=SuscepImp is default): Files begining with \"String\" and containing \"zEQ\" will be used for z-averaging \n";
  print "    --verbose/-v : Prints the values on the screen \n";
  print "    --noask  : Does not confirm the files to use (use with care!) \n";
  print "    --getext  : Gets extension from output*zEQ1*.txt  \n";
  exit;

}
###############################
########################
##   Main program     ##
########################

use Getopt::Long;

$PrefixFile="SuscepImp";

if ( @ARGV > 0 ){
GetOptions('dir0=i'=>\$dir0,
           'dirF=i'=>\$dirF,
	   'type=s'=>\$PrefixFile,
           'v|verbose'=>\$PrintStuff,
           'noask'=>\$NoAsk,
           'getext'=>\$GetExt,
           'h|help'=>\$help);
}
##else{Print_Help();}
if (defined($help)){Print_Help();}

##
## List of files
##
# strip dIdV
##$FilePattern =~ s/dIdV//;
$FilePattern="zEQ";

$OutputExt="";
if (defined($GetExt)){
  my @OutputFiles=();
  @OutputFiles=glob("output*zEQ1.txt");
  $OutputExt=$OutputFiles[0];
  $OutputExt =~ s/output_//;
  $OutputExt =~ s/zEQ1.txt//;
  print " Extension: $OutputExt \n"; 
}

print "Looking for files ".$PrefixFile."\*".$FilePattern."\*"."dat\n";
my $CurrentDir= $ENV{'PWD'};
print "in current dir: ".$CurrentDir."\n";
##chdir $CurrentDir;
my @MatchingFiles=();
##@MatchingFiles=glob("$PrefixFile*$FilePattern*.dat");
@MatchingFiles=glob "${PrefixFile}*${FilePattern}*dat";

##foreach $FileName (@MatchingFiles){print $FileName."\n";}

##exit(0);

$noFiles=0;
$izeq1=0;
foreach $FileName (@MatchingFiles){
 print " File $noFiles is: $FileName \n";
 if ($FileName =~ m/(zEQ1|zEQ1\.00)\.dat/){
   $izeq1=$noFiles;
   $FileNameZeq1=$FileName;
   print "   --- File $izeq1 is the z=1 file \n";
 }
 $noFiles++;
}

## Checks if OutFileName is here...
$OutFileName=$FileNameZeq1;
##$OutFileName =~ s/zEQ1/zEQAVG/;
$OutFileName =~ s/zEQ1/zEQAVG/;
if (defined($GetExt)){
  $OutputExt.="zEQAVG.dat";
  $OutFileName =~ s/zEQAVG.dat/$OutputExt/;
}
##if (-e $OutFileName){
if ((-e $OutFileName)&&(!defined($NoAsk))){
 print "$OutFileName exists in this dir. Remove (y/n)? -> ";
 chomp($contYN = <STDIN>);
 if ($contYN eq "y"){print "Removing $OutFileName \n";
   system("rm $OutFileName");
   @MatchingFiles=glob("$PrefixFile*$FilePattern*.dat");
   $noFiles=0;
   $izeq1=0;
   foreach $FileName (@MatchingFiles){
     print " File $noFiles is: $FileName \n";
     if ($FileName =~ m/(zEQ1|zEQ1\.00)\.dat/){
      $izeq1=$noFiles;
      $FileNameZeq1=$FileName;
      print "   --- File $izeq1 is the z=1 file \n";
     }
     $noFiles++;
   }
## end second loop in Matching Files

 } 
## end if Remove file
} 
## end check in Output file exists
if (!defined($NoAsk)){
  print "\n $noFiles files will be processed. Output file: $OutFileName  \n Continue (y/n)? -> ";
  chomp($contYN = <STDIN>);
  if ($contYN ne "y"){print "Exiting... \n";exit(0);}
}
print "Ok, lets do it. \n";



## 4 columns: Temp SusChainPImp SusChain SusImp for EACH file
##
@Temps=();
@SusChainPImp=();
@SusChain=();
@SusImp=();
@TempArray=();
@NtempsFile=();
my $ntemps=0;
my $ntempsChain=0;

##########################
##   Read Files         ##
##########################

$ifile=0;
my $LowTempCutOff=1e-100;
foreach $FileName (@MatchingFiles){
## Add an empty Array to each
  @tmp=();

  push(@Temps,[ @tmp ] );
  push(@SusChainPImp,[ @tmp ] );
  push(@SusChain,[ @tmp ] );
  push(@SusImp,[ @tmp ] );

  ## When using array of arrays, parameters to function call are like THIS:
  ($ntemps,$LowTemp)=ReadSusFile($Temps[$ifile],$SusChainPImp[$ifile],$SusChain[$ifile],$SusImp[$ifile],$FileName);
  if ($LowTemp>$LowTempCutOff){$LowTempCutOff=$LowTemp;}
  ## LowTempCutOff will be the lowest energy present in ALL files.
  push(@NtempsFile, $ntemps);
  $ifile++;

}

my $NumFiles=$ifile;
## NOT assuming ntemps is the same for all data files...

#print " Num files = $NumFiles = ".($#Temps+1)."; ntemps_Zeq1 = ".$NtempsFile[$izeq1]." \n";

print "Low Energy Cut off = $LowTempCutOff \n";

my @BaseXtemps=();
for ($ii=0;$ii<=$NtempsFile[$izeq1];$ii++){
  push(@BaseXtemps,$Temps[$izeq1][$ii]);
}


#####################
##   Z-Averaging   ##
#####################

my @YAvg=();

for ($ii=0;$ii<=$#BaseXtemps;$ii++){
## Tests
##for ($ii=70;$ii<=70;$ii++){
  my $y1=0.0;
  my $x1test=$BaseXtemps[$ii];
  my $yaux=0.0;
  if (defined($PrintStuff)){print $x1test." -- ";}
   for ($ifile=0;$ifile<$NumFiles;$ifile++){
##  for ($ifile=0;$ifile<1;$ifile++){
   $y1=Interpol($Temps[$ifile],$SusImp[$ifile],$x1test);
   if (defined($PrintStuff)){print $y1."  ";}
   $yaux+=$y1;
  }
  $yaux/=$NumFiles;
  if (defined($PrintStuff)){print " : ".$yaux."\n";}
  push(@YAvg,$yaux);
}


############################
##   Save in a New File   ##
############################

print "Saving in $OutFileName \n";

open(OUTFILE,"> ".$OutFileName);
##  for ($ii=0;$ii<=$#BaseXtemps-20;$ii++){
  for ($ii=0;$ii<=$#BaseXtemps-1;$ii++){
    if (abs($BaseXtemps[$ii])>=$LowTempCutOff){
      if ($PrefixFile =~ m/rho/){
       ## 5 cols and integer in rho
    printf(OUTFILE "%20.20e %20.20e %20.20e %20.20e %20.20e %i\n",$BaseXtemps[$ii],$YAvg[$ii],$YAvg[$ii],$YAvg[$ii],$YAvg[$ii],0);
       } else {
       ## 4 cols in Suscep
    printf(OUTFILE "%20.20e %20.20e %20.20e %20.20e\n",$BaseXtemps[$ii],$YAvg[$ii],$YAvg[$ii],$YAvg[$ii]); 
      } ## end if PrefixFile is rho (density)
    }
    ## only save energies above than cut-off
  }
close(OUTFILE);

print "...Done!\n";
#######
# END #
#######
