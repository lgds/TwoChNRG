#!/usr/bin/perl
##
## Make a histogram out of input data
##
###############################

sub ReadRawFile{

## $Nys=ReadRawFile(\@XArray,\@YArray,$RawDataFilename,$ycol);


my $XArray=shift(@_);
my $YArray=shift(@_);
my $FileName=shift(@_);
my $ycol=shift(@_);

my $Nys=0;

@$XArray=();
@$YArray=();


my $FirstLine=0;

my @ThisLine=();
my $auxY=0;
my $auxX=0;

 open(INFILE,$FileName);
   while (<INFILE>){
     $TheLine = $_;
     chomp($TheLine);
     @ThisLine=split(/ +/,$TheLine);
##     print " No entries in line : ".($#ThisLine+1)."\n";
     if ($#ThisLine>=0){ ## no empty lines
       $auxX=$ThisLine[0];
       $auxY=$ThisLine[$ycol];
       push(@$XArray,$auxX);
       push(@$YArray,$auxY);
##       print " Y = $auxY \n";
       $Nys+=1;
     }
   }
 close(INFILE);

 return($Nys);
}

###############################


###############################

sub CalcHist{

##CalcHist(\@YArray,\@HistArray,\@xHistArray, $Nys,$Nbins,$Y0,$YF,$Ystep,$UseLog,$Lambda,$AddExt);


use List::Util qw( min max );


my $NoArgs=@_;

my $YArray=shift(@_);
my $HistArray=shift(@_);
my $xHistArray=shift(@_);
my $Nys=shift(@_);
my $Nbins=shift(@_);
my $YMin=shift(@_);
my $YMax=shift(@_);
my $Ystep=shift(@_);
my $UseLog=shift(@_);
my $Lambda=shift(@_);
my $Extension=shift(@_);

my $icount=0;

my $Nhist=0;
my $Area=0;
##my $HistFileName="HistTsq_n".$Ndisc.".dat";

##my $YMax=max @$YArray;
##my $YMin=min @$YArray;
my $ThisY=0.0;

print "YMin = $YMin, YMax = $YMax, Nbins= $Nbins, UseLog = $UseLog Lambda = $Lambda \n";

@$HistArray=();

my $Ylow=1.0;
my $Yhigh=1.0;

## Need to work on this/ YMin and YMax are FIXED.
## Nbins will then be fixed too.
if ($UseLog==1){
##  $Ylow=$YMax*$Lambda**(-$Nbins);
  $Nbins=int(log($YMax/$YMin)/log($Lambda))+1;
  print " New Nbin = $Nbins \n";
  $Ylow=$YMin;
  $Yhigh=$Ylow*$Lambda;
##  if ($YMin>$Ylow){
##    print " Ops, Ylow = $Ylow < YMin \n";}
}
else{
  $Ylow=0.9999*$YMin;
  $Yhigh=$YMin+$Ystep;
}


##open(OUTFILE,"> histogram".$Extension.".dat");

my $OldFrq=0;
my $Frq=0;
my $OldY=0;


for ($icount=0;$icount<=$Nbins;$icount+=1){

  $OldFrq=$Frq;
  $Frq=0;

  foreach $ThisY (@$YArray){

##    print "$Ylow <  $ThisY <= $Yhigh ?\n";

    if ( ($ThisY>$Ylow)&&($ThisY<=$Yhigh) ){$Frq++;}


  }
  #end foreach
  $Nhist+=$Frq;
  $Area+=0.5*($Frq+$OldFrq)*($Ylow-$OldY);

##  my $xHist=0.5*($Yhigh+$Ylow);
##  if ($UseLog==1){$xHist=$Ylow*sqrt($Lambda);}
  my $xHist=$Ylow; ## P(y): no of entries between y and y+Dy OR
                   ## P(y): no of entries between y and y*Lambda


  push(@$HistArray,$Frq);
  push(@$xHistArray,$xHist);


##  print " $xHist   $Frq \n";
##  print OUTFILE " $xHist   $Frq\n";

  $OldY=$Ylow;
  $Ylow=$Yhigh;

  if ($UseLog==1){$Yhigh*=$Lambda;}
  else{$Yhigh+=$Ystep;}

}
#end loop over bins

##close(OUTFILE);

return($Nhist,$Area);
}
# End CalcHist
#
#

sub SaveHist{

##SaveHist(\@HistArray,\@xHistArray,$Nhist,$AddExt,$Area);

use List::Util qw( min max );

my $NoArgs=@_;

my $HistArray=shift(@_);
my $xHistArray=shift(@_);
my $Nhist=shift(@_);
my $Extension=shift(@_);
my $Area=shift(@_);

my $ArraySize=$#$HistArray+1;

if ($Nhist==$ArraySize){
   print "SaveHIst: ok\n";
}else{
   print "SaveHist: Not ok\n";
   print "ArraySize=$ArraySize, Nhist = $Nhist\n";
}


open(OUTFILE,"> histogram".$Extension.".dat");
my $icount=0;
my $Freq=0;

my $OldFrq=0;
my $xHist=@$xHistArray[0];
my $OldX=0;
my $ThisArea=0;
my $TotArea=0;

foreach $Freq (@$HistArray){
  $OldX=$xHist;
  $xHist=@$xHistArray[$icount];

  $ThisArea=0.5*($Freq/$Area+$OldFrq/$Area)*($xHist-$OldX);
  $TotArea+=$ThisArea;

  print "$xHist  $Freq  ".($Freq/$Nhist)." ".($Freq/$Area)."\n";
  print OUTFILE " $xHist   $Freq  ".($Freq/$Nhist)." ".($Freq/$Area)."\n";
  $OldFrq=$Freq;
  $icount++;
}


close(OUTFILE);
my $NormFrq=$Frq/$Nhist;

print "Area = $Area , NormArea = $TotArea \n";

}
# End SaveHist
#
#



###
##
sub Print_Help() {
  print " Usage $0 --file=dataxy.dat --nb=nbins \n"; 
  print " Options:  \n";
  print "    --y0=y0          : default: min value in file \n";
  print "    --yF=yF          : default: max value in file \n";
  print "    --ystep=ystep    : instead of nbins \n";
  print "    --col=(0,1,2...) : 1 by default \n";
  print "    --uselog=Lambda  : (needs nbins)  yF*Lambda^{-n+1}<y<yF*Lambda^{-n} down to yF*Lambda^{-nbins}\n";
  print "    --addext=text    : Adds \"text\" to outputfile name \n";
  print "    --special        : Does something with the data (edit it in code) \n";
  exit;

}
###############################




########################
##   Main program     ##
########################

use Getopt::Long;

use List::Util qw( min max );

my $ycol=1;
my $UseLog=0;
my $Normalize=0;


if ( @ARGV > 0 ){
GetOptions('file=s'=>\$RawDataFilename,
	   'col=i' => \$ycol,
	   'uselog=f' => \$Lambda,
	   'nb=i' => \$Nbins,
	   'y0=f' => \$Y0,
	   'yF=f' => \$YF,
	   'ystep=f' => \$Ystep,
	   'addext=s' => \$AddExt,
	   'norm' => \$Norm,
	   'special' => \$Special,
           'h|help'=>\$help);
}
##else{Print_Help();}
if (defined($help)){Print_Help();}
if ( !defined($RawDataFilename) ){Print_Help();}

if (defined($Lambda)){$UseLog=1;}
if (defined($Norm)){$Normalize=1;}

if ( (!defined($Nbins))&&(!defined($Ystep)) ){Print_Help();}

if (defined($AddExt)){$AddExt="_".$AddExt;}

my @XArray=();
my @YArray=();
my @HistArray=();
my @xHistArray=();
my $nrealiz=1;
my $Nys=-1;
my $Nhist=0;
my $Area=0;

###########################
##   Read RawDataFile    ##
###########################

$Nys=ReadRawFile(\@XArray,\@YArray,$RawDataFilename,$ycol);

print " No entries = $Nys \n";

##if ( (!defined($YF))||($YF<(max @YArray)) ){$YF= max @YArray;}
##if ( (!defined($Y0))||($Y0>(min @YArray)) ){$Y0= min @YArray;}

if ( (!defined($YF)) ){$YF= max @YArray;}
if ( (!defined($Y0)) ){$Y0= min @YArray;}

if (defined($Special)){
## Nov 2013: take the log
  foreach my $ThisY (@YArray){
    $ThisY=log($ThisY);
  }
  $Y0=log($Y0);
  $YF=log($YF);
}
#end special

if ($YF<(max @YArray)){print "Warning: YF<max data \n";}
if ($Y0>(min @YArray)){print "Warning: Y0>min data: Y0 = $Y0 min data =".(min @YArray)." \n";}

###########################
##  Set Nbins or $Ystep  ##
###########################

if (!defined($Nbins)){$Nbins= int(($YF-$Y0)/$Ystep);}

if (!defined($Ystep)){$Ystep=($YF-$Y0)/$Nbins;}


###########################
##   Make Histogram      ##
###########################


($Nhist,$Area)=CalcHist(\@YArray,\@HistArray,\@xHistArray, $Nys,$Nbins,$Y0,$YF,$Ystep,$UseLog,$Lambda,$AddExt);


print "Nhist = $Nhist Area =$Area \n";

SaveHist(\@HistArray,\@xHistArray,$Nhist,$AddExt,$Area);

print "...Done!\n";
#######
# END #
#######


##############
# TRASH BIN  #
##############
