#!/usr/bin/perl
##
## Split disorder realizations into HybDeltas file; update input_nrg.
##
###############################

sub ReadRawFile{

my $EnArray=shift(@_);
my $tnSqArray=shift(@_);
my $E0Array=shift(@_);
my $FileName=shift(@_);
my $nrealiz=1;
my $Nens=-1;

@$EnArray=();
@$tnSqArray=();
@$E0Array=();

my $FirstLine=0;

my $auxEn=0;
my $auxtn=0;

 open(INFILE,$FileName);
   while (<INFILE>){
     $TheLine = $_;
     chomp($TheLine);
## Remove leading spaces
     $TheLine =~ s/^\s+//;
     ($auxEn,$auxtn)=split(/ +/,$TheLine);
     if ($FirstLine==0){
       $Nens=$auxtn;
##       print "En =  $auxEn \n";
       push(@$E0Array,$auxEn);
       $FirstLine=1;
     }else{
       if ($auxtn==$Nens){
##	 print "En =  $auxEn \n";
	 push(@$E0Array,$auxEn);
	 $nrealiz++;
       }else{
	 push(@$EnArray,$auxEn);
	 push(@$tnSqArray,$auxtn);
       }
       ## end if Separator

     }
## end if FirstLine
   }
 close(INFILE);

 return($nrealiz,$Nens);
}

###############################

sub RescaleEns{

my $EnArray=shift(@_);
my $tnSqArray=shift(@_);
my $E0Array=shift(@_);
my $nrealiz=shift(@_);
my $Nens=shift(@_);
my $Dband=shift(@_);
my $mu=shift(@_);
my $vzero=shift(@_);

my $icount=0;
my $EnMax=-10;
my $auxEnPlus=0;
my $auxEnMinus=0;



for (my $ir=0;$ir<$nrealiz;$ir++){
  $auxEnPlus=@$EnArray[$icount];
#  print " Real: $ir ; EnPlus= $auxEnPlus \n";
  if (abs($auxEnPlus)>$EnMax){$EnMax=abs($auxEnPlus);}
  $icount+=$Nens-1;
  $auxEnMinus=@$EnArray[$icount];
#  print " Real: $ir ; EnMinus= $auxEnMinus \n";
  if (abs($auxEnMinus)>$EnMax){$EnMax=abs($auxEnMinus);}
  $icount++;
}
#loop in realiz

print " EnMax = $EnMax \n";
##$Dband=$EnMax;

if ($EnMax<$Dband){print "Problems: Enmax < Dband \n"; exit(0);}

$icount=0;
for (my $ir=0;$ir<$nrealiz;$ir++){
##  print " Old: @$E0Array[$ir] ";
  @$E0Array[$ir]/=$Dband;
##
##  Changing the coupling Udis->sqrt(V_0)*Udis
##
  @$E0Array[$ir]*=sqrt($vzero); ## SHOULD BE sqrt(v0)!!!
  @$E0Array[$ir]-=$mu;
##  print " New: @$E0Array[$ir] \n";
  for (my $ien=0;$ien<$Nens;$ien++){
    @$EnArray[$icount]/=$Dband;
    @$EnArray[$icount]-=$mu;
##
##  NOTE: EnArray will also change with v0: the off-diagonal terms in H_PP will increase!! 
##
##
## We need to SUBTRACT mu beause we are changing the energy axis: 
##   e->e1=e-mu (-1-mu < e1 < 1-mu)
## BUT 
## if we only have the expression for f(e) THEN
##  f(e)=f(e1+mu) 
## so, for a given e1, we need to calculate f(e1+mu)
##  If we have (x_n, f_n), then we are ok: just do x_n-=mu, keeping f_n.

    @$tnSqArray[$icount]/=($Dband*$Dband); ## Is this correct?
##
## Yes, for the following reason: say all t's are equal. Then
## Gamma(e)=pi.t^2.rho(e) -> int_{e} Gamma(e) de = pi.t^2.Nst
## Transf: e'=e/D
##
## int_{e'} Gamma(e') de' = (1/D) int_{e} Gamma(e) de
##
## then
##
## int_{e'} Gamma(e') de' = pi.(t^2/D).Nst 
## or
## Gamma(e')=pi.(t^2/D).rho(e')
##
##  BUT, I want to have Gamma in units of D as well. So
## 
##  Gamma(e')/D=pi.(t^2/D^2).rho(e')
## so I need t^2/D^2 up there
##
##  Changing the coupling Udis->sqrt(V_0)*Udis
##
    @$tnSqArray[$icount]*=$vzero; 

    $icount++;
  }
## loop in ens
  $auxEnPlus=@$EnArray[$ir*$Nens];
##  print " Real: $ir ; EnPlus= $auxEnPlus \n";
  $auxEnMinus=@$EnArray[($ir+1)*$Nens-1];
##  print " Real: $ir ; EnMinus= $auxEnMinus \n";

}
#loop in realiz

return;
}
###############################

sub SaveBandFiles{

##   SaveBandFiles(\@EnArray,\@tnSqArray,\@BinPositive,\@BinNegative,\@E0Array,$nrealiz,$Nens,$dir0,$dirF,$ir0,$DeltaE,$SavePlot,$LambdaNRG,$ztwist,$mu);
##

my $NoArgs=@_;

my $EnArray=shift(@_);
my $tnSqArray=shift(@_);
my $BinPositive=shift(@_);
my $BinNegative=shift(@_);
my $E0Array=shift(@_);

my $nrealiz=shift(@_);
my $Nens=shift(@_);
my $dir0=shift(@_);
my $dirF=shift(@_);
my $ir0=shift(@_);
my $DeltaE=shift(@_);
my $SavePlot=shift(@_);
my $LambdaNRG=shift(@_);
my $ztwist=shift(@_);
my $mu=shift(@_);


my $EnMax=-10;
my $auxEnPlus=0;
my $auxEnMinus=0;

my $ir=$ir0;
my $idir=$dir0;

my $OutFileName="HybDeltas.dat";

my $Npos=@$BinPositive;
my $Nneg=@$BinNegative;



#for (my $ii=0; $ii<4; $ii++){

#  my $irand=int(rand($Npos));
#  print "irand = $irand Positive = ".@$BinPositive[$irand]."\n";
#  my $irand=int(rand($Nneg));
#  print "irand = $irand Negative = ".@$BinNegative[$irand]."\n";

#}


my $elim0=$LambdaNRG**(-$Ndisc-$ztwist);
my $elim1=$LambdaNRG**(-$Ndisc-$ztwist+1);


##my $icount=0;
my $icount=$ir0*$Nens;

while ( ($ir<$nrealiz)&&($idir<=$dirF) ){

##
## changes E0 to E0[ir]
##

  my $Ed=@$E0Array[$ir];
  my $ChgCommand="./ChangeNRGParams.perl -i $idir  -f $idir --choiceparam=4 --p0=$Ed --pF=$Ed   --pstep=0 \n";

  ##print "Executing $ChgCommand \n";
  system($ChgCommand);
  $ChgCommand="./ChangeNRGParams.perl -i $idir  -f $idir --choiceparam=6 --p0=$mu --pF=$mu   --pstep=0.0 \n";
  system($ChgCommand);



##
## cd into dir
##

  my $nodedirname='./run'.$idir;
  chdir($nodedirname) or die "Could not change dir to $nodedirname \n";

  my $PlotFileName="HybPlot_".$idir.".dat";
  my $PlotFileNameNewEns="HybPlotAddEns_".$idir.".dat";


##  my $DeltaE=0.1; ## input???


##
## Saves HybFunc.dat
##


  CalcGamma(\@$EnArray,\@$tnSqArray,$nrealiz,$Nens,$DeltaE,$ir,1,$mu);



  open(OUTFILE,"> ".$OutFileName);
  if ($SavePlot){open(PLOTFILE,"> ".$PlotFileName);}
  if ($SavePlot){open(PLOTFILE2,"> ".$PlotFileNameNewEns);}
  ## end if SavePlot
  for (my $ii=0;$ii<$Nens;$ii++){
    my $lastN=int(-log(abs(@$EnArray[$icount]))/log($LambdaNRG))+2-$ztwist;


    printf(OUTFILE "%20.20e  %20.20e\n",@$EnArray[$icount],@$tnSqArray[$icount]);

    ##
    ## Check for change of sign (moved to trash down below)
    ##

    if ($SavePlot){
      printf(PLOTFILE "%20.20e  %20.20e\n",@$EnArray[$icount],0.0);
      printf(PLOTFILE "%20.20e  %20.20e\n",@$EnArray[$icount],@$tnSqArray[$icount]);
      printf(PLOTFILE "%20.20e  %20.20e\n",@$EnArray[$icount],0.0);
    }
    ## end if SavePlot
    $icount++;
  }
  if ($SavePlot){close(PLOTFILE);close(PLOTFILE2);}
  close(OUTFILE);

  chdir('../');

  $ir++;
  $idir++;
}

print "SaveBandFiles: ir0=$ir0 Last ir=".($ir-1)." \n";



return;

}
##
## End Save Band files
##

###############################

sub CalcPn{

## CalcPn(\@EnArray,\@tnSqArray,$nrealiz,$Nens,$Ndisc,$LambdaNRG)

use List::Util qw( min max );


my $NoArgs=@_;

my $EnArray=shift(@_);
my $tnSqArray=shift(@_);
##my $E0Array=shift(@_);
my $nrealiz=shift(@_);
my $Nens=shift(@_);
my $Ndisc=shift(@_);
my $LambdaNRG=shift(@_);

#my $elim0=$LambdaNRG**(-$Ndisc-1);
#my $elim1=$LambdaNRG**(-$Ndisc);

my $ztwist=1;
my $elim0=$LambdaNRG**(-$Ndisc-$ztwist);
my $elim1=$LambdaNRG**(-$Ndisc-$ztwist+1);


my $ir=0;
my $icount=0;

my $Nhist=0;
my @HistTsqArray=();
my $HistFileName="HistTsq_n".$Ndisc.".dat";


print "Calculating Pn: n= $Ndisc, e0 = $elim0, e1 = $elim1, Nr = $nrealiz Nens = $Nens\n";

while ($ir<$nrealiz){
  for (my $ii=0;$ii<$Nens;$ii++){
    if ( (@$EnArray[$icount]>=$elim0)&&(@$EnArray[$icount]<=$elim1) ){
      push(@HistTsqArray,@$tnSqArray[$icount]);
      $Nhist++;
      ##print "Found: e=".@$EnArray[$icount]." Nhist = ".$Nhist."\n";

    }
    ## end if elim0<e<elim1
    $icount++;
  }
  $ir++;
}


my $MaxTsq=max @HistTsqArray;


##$DeltaHist=($MaxTsq/100);

$DeltaHist=1e-06;

print "Ndisc = $Ndisc Max t2 = $MaxTsq, Nhist = $Nhist, Delta = $DeltaHist \n";

if ($Nhist>0){
  SaveHistogram(\@HistTsqArray,$Nhist,$DeltaHist,$MaxTsq,$HistFileName);
} else {print "Ops: Nhist=$Nhist \n";}


return($Nhist);
}
# End CalcPn
#
#

sub SaveHistogram{

my $HistArray=shift(@_);
my $Nhist=shift(@_);
my $Deltax=shift(@_);
my $MaxTsq=shift(@_);
my $HistFileName=shift(@_);

my $ThisTsq=0;
my $xhist=0;


open(OUTFILE,"> ".$HistFileName);


for ($xhist=0;$xhist<=$MaxTsq;$xhist+=$Deltax){

  my $Frq=0;

  foreach $ThisTsq (@$HistArray){
    if ( ($ThisTsq>=$xhist)&&($ThisTsq<=$xhist+$Deltax) ){$Frq++;}
  }
  #end foreach

  printf(OUTFILE "%20.20e  %20.20e %5d \n",0.5*($xhist+$Deltax),$Frq/$Nhist,$Frq);
  ##printf("%20.20e  %20.20e %5d \n",0.5*($xhist+$Deltax),$Frq/$Nhist,$Frq);


}
##end hist


close(OUTFILE);

return;
}
# End SaveHistogram
#
#

###############################

sub GetBin{
##
## Adds missing energies
##
## GetBin(\@EnArray,\@tnSqArray,\@BinPositive,\@BinNegative,$nrealiz,$Nens,$Ndisc,$LambdaNRG)

use List::Util qw( min max );

my $NoArgs=@_;

my $EnArray=shift(@_);
my $tnSqArray=shift(@_);
my $BinPositive=shift(@_);
my $BinNegative=shift(@_);
my $nrealiz=shift(@_);
my $Nens=shift(@_);
my $Ndisc=shift(@_);
my $LambdaNRG=shift(@_);



#my $elim0=$LambdaNRG**(-$Ndisc-1);
#my $elim1=$LambdaNRG**(-$Ndisc);

my $ztwist=1;
my $elim0=$LambdaNRG**(-$Ndisc-$ztwist);
my $elim1=$LambdaNRG**(-$Ndisc-$ztwist+1);


@$BinPositive=();
@$BinNegative=();


my $ir=0;
my $count=0;
my $Npos=0;
my $Nneg=0;


while ($ir<$nrealiz){
  for (my $ii=0;$ii<$Nens;$ii++){

    # Positive energies
    if ( (@$EnArray[$icount]>=$elim0)&&(@$EnArray[$icount]<=$elim1) ){
      push(@$BinPositive,@$tnSqArray[$icount]);
      ##print "Found: e=".@$EnArray[$icount]." Nhist = ".$Nhist."\n";
      $Npos++;
    }
    ## end if elim0<e<elim1

    # Negative energies
    if( (@$EnArray[$icount]>=(-1)*$elim1)&&(@$EnArray[$icount]<=(-1)*$elim0) ){
      push(@$BinNegative,@$tnSqArray[$icount]);
      ##print "Found: e=".@$EnArray[$icount]." Nhist = ".$Nhist."\n";
      $Nneg++;

    }
    ## end if -elim1<e<-elim0

    $icount++;

  }
  $ir++;
}
## end loop

print "Npos = $Npos ; Nneg = $Nneg \n";



return;
}
# End GetBin
#
#

sub CalcGamma{

##CalcGamma(\@EnArray,\@tnSqArray,$nrealiz,$Nens,$DeltaE,$ir,$InDir);

use Math::Trig;

my $NoArgs=@_;

my $EnArray=shift(@_);
my $tnSqArray=shift(@_);
my $nrealiz=shift(@_);
my $Nens=shift(@_);
my $DeltaE=shift(@_);
my $ir=shift(@_);
my $InDir=shift(@_);
## number from 0 to nrealiz-1
my $mu=shift(@_);


my $countE=0;
##my $En=-1.0;
my $En=-(1.0+$mu);

if ($ir<0){

  my $CountFileName="N_E_Delta".$DeltaE.".dat";

  my @EnArraySorted=sort {$a <=> $b} @$EnArray; ## ascending


  my $NensTot=@EnArraySorted;
  my $ThisEbeta=0;


  print "Nens = $Nens  NensTot= $NensTot =(?) ".($Nens*$nrealiz)."\n";
  open(COUNTFILE,">".$CountFileName);

  foreach $ThisEbeta (@EnArraySorted) {
    #for (my $ii=0;$ii<11;$ii++){
    #  $ThisEbeta=@EnArraySorted[$ii];
    #  print "Ebeta = $ThisEbeta En=$En En+Delta=".($En+$DeltaE)."\n";

    if ( ($ThisEbeta>=$En)&&($ThisEbeta<=$En+$DeltaE) ) {
      $countE++;
    } else {
      ##    print "E= ".($En+$DeltaE/2)." N_E= ".$countE." \n";
      print COUNTFILE ($En+$DeltaE/2)." ".($countE/$DeltaE)." ".($countE/($DeltaE*$NensTot))."\n";
      $countE=0;
      $En+=$DeltaE;
    }
    ## end if

  }
  ## end foreach
## Need the last one too.
  print COUNTFILE ($En+$DeltaE/2)." ".($countE/$DeltaE)." ".($countE/($DeltaE*$NensTot))."\n";
  close(COUNTFILE);

} else {

  if ($InDir==1){
    $CountFileName="HybFunc.dat";}
  else{
    $CountFileName="N_E_Delta".$DeltaE."_realiz".$ir.".dat";}

  open(COUNTFILE,">".$CountFileName);

##
## Need to reverse...
## 

##  my $icount=$ir*$Nens;
  my $icount=($ir+1)*$Nens-1;
  $countE=0;
##  $En=-1.0-$DeltaE/2;
  $En=-(1.0+$mu)-$DeltaE/2;

  my $tsq_avg=0;

  for (my $ii=0;$ii<$Nens;$ii++) {

    ##print "Ei = ".@$EnArray[$icount]." Ei-1 = ".@$EnArray[$icount-1]." E= ".$En." E+Delta = ".($En+$DeltaE)." \n";

    if ( (@$EnArray[$icount]>=$En)&&(@$EnArray[$icount]<=$En+$DeltaE) ) {
      $countE++;
      ## also calculate |t|^2(E)
      $tsq_avg+=@$tnSqArray[$icount];
    } elsif (@$EnArray[$icount]>$En+$DeltaE){
      ##print "E= ".($En+$DeltaE/2)." N_E= ".$countE." \n";
      if ($countE!=0) {
	$tsq_avg/=$countE;
      }
      if ($InDir==1){
	print COUNTFILE ($En+$DeltaE/2)." ".(pi*$tsq_avg*$countE/$DeltaE)."\n";}
      else{
	print COUNTFILE ($En+$DeltaE/2)." ".($countE/$DeltaE)." ".$countE/($DeltaE*$Nens)." ".$tsq_avg."  ".(pi*$tsq_avg*$countE/$DeltaE)."\n";}
      $countE=0;
      $tsq_avg=0;
      $En+=$DeltaE;
    }
    ## end if
    ##    $icount++;
    $icount--;
  }
  ## end for

  ## Need the last one too.
  if ($countE!=0) {
    $tsq_avg/=$countE;
  }
  if ($InDir==1){
    print COUNTFILE ($En+$DeltaE/2)." ".(pi*$tsq_avg*$countE/$DeltaE)."\n";}
  else{
    print COUNTFILE ($En+$DeltaE/2)." ".($countE/$DeltaE)." ".$countE/($DeltaE*$Nens)." ".$tsq_avg."  ".(pi*$tsq_avg*$countE/$DeltaE)."\n";}

  close(COUNTFILE);

}
## end if ir<0

return;
}
# End CalcGamma 
#
#


###
##
sub Print_Help() {
  print " Usage $0 --file=params_kondo_file -i dir0 -f dirN\n"; 
  print " Options:  \n";
  print "    -D , --Dband=Dband (in units of t. D=3.1 by default) \n";
  print "    --saveplot    : Save \"grass\" data in plottable format \n";
  print "    --hist        : Calculate histograms of |t|^2 \n";
  print "    --dos=deltaE  : Calulate Hyb with spacing deltaE(default=0.1) \n";
  print "    --ir=ir0      : Start from realization ir0 (default =0) \n";
  print "    --nosave      : Does not save in dirs (Histograms/DoS only).  \n";
  print "    --mu=mu       : displace energies E->E-mu   \n";
  print "    --vzero=v0    : t^2->t^2*v0 (default v0=1)   \n";
#  print "    --shiftchain=nshift : Shift chain values by nshift(>0). Does not read chain file.  \n";
  exit;

}
###############################




########################
##   Main program     ##
########################

use Getopt::Long;

my $Dband=3;
my $ir0=0;
my $DeltaE=0.1;
my $mu=0.0;
my $vzero=1.0;

if ( @ARGV > 0 ){
GetOptions('file=s'=>\$RawDataFilename,
           'D|Dband=f'=>\$Dband,
           'i=i'=>\$dir0,
           'f=i'=>\$dirF,
           'saveplot'=>\$SavePlot,
           'hist'=>\$CalcHists,
           'dos=f'=>\$DeltaE,
           'mu=f'=>\$mu,
           'vzero=f'=>\$vzero,
           'ir=i'=>\$ir0,
#           'overwrite'=>\$OvrtFile,
           'nosave'=>\$NoSave,
           'h|help'=>\$help);
}
##else{Print_Help();}
if (defined($help)){Print_Help();}
if ( (!defined($NoSave))&&( (!defined($dir0))||(!defined($dirF)) )  ){Print_Help();}


my @EnArray=();
my @tnSqArray=();
my @E0Array=();
my $nrealiz=1;
my $Nens=-1;



###########################
##   Read RawDataFile    ##
###########################

($nrealiz,$Nens)=ReadRawFile(\@EnArray,\@tnSqArray,\@E0Array,$RawDataFilename);

print " No Realiz =  $nrealiz; No Ens = $Nens \n";

#################################
##   Rescale energies          ##
#################################

RescaleEns(\@EnArray,\@tnSqArray,\@E0Array,$nrealiz,$Nens,$Dband,$mu,$vzero);

#################################
##   Get representative bin    ##
#################################

my $LambdaNRG=2.5;
my $ztwist=1;
my $Ndisc=0;


GetBin(\@EnArray,\@tnSqArray,\@BinPositive,\@BinNegative,$nrealiz,$Nens,$Ndisc,$LambdaNRG);



#################################
##   Histograms (NEW)          ##
#################################


my $Nhist=100;

if (defined($CalcHists)){

 while ($Nhist>0){
   $Nhist=CalcPn(\@EnArray,\@tnSqArray,$nrealiz,$Nens,$Ndisc,$LambdaNRG);
##   print "Ndisc = $Ndisc, Nhist = $Nhist \n";
   $Ndisc++;
 }
### end while
}
## end if CalcHists

#################################
##   Density of states         ##
#################################

if (defined($DeltaE)){

 CalcGamma(\@EnArray,\@tnSqArray,$nrealiz,$Nens,$DeltaE,$ir0,0,$mu);
}

#################################
##   Save in different files   ##
#################################

if (!defined($NoSave)){
  SaveBandFiles(\@EnArray,\@tnSqArray,\@BinPositive,\@BinNegative,\@E0Array,$nrealiz,$Nens,$dir0,$dirF,$ir0,$DeltaE,$SavePlot,$LambdaNRG,$ztwist,$mu);
}

print "...Done!\n";
#######
# END #
#######


##############
# TRASH BIN  #
##############


    ##
    ## Check for change of sign
    ##
#    if ( (@$EnArray[$icount]*@$EnArray[$icount+1]<0)&&(@$EnArray[$icount]>0) ){
#      print "changed sign: E = ".@$EnArray[$icount],"\n";
#      print " E = ".@$EnArray[$icount]." is on bin (".($lastN-1)." - ".$lastN.")\n";
##
##    Picking new |t|^2: positive energies
##
#      for (my $iN=$lastN+1;$iN<=50;$iN++){
#        for (my $ii=1; $ii<=1; $ii++){
#	  my $irand=int(rand($Npos));
#	  my $En=$LambdaNRG**(-($iN+$ztwist-1+$ii*0.5));
#	  ##print "iN = $iN, ii=$ii New E = $En t2 = ".@$BinPositive[$irand]."\n";
###	  printf(OUTFILE "%20.20e  %20.20e\n",$En,@$BinPositive[$irand]);
#	  if ($SavePlot){
#	    printf(PLOTFILE2 "%20.20e  %20.20e\n",$En,0.0);
#	    printf(PLOTFILE2 "%20.20e  %20.20e\n",$En,@$BinPositive[$irand]);
#	    printf(PLOTFILE2 "%20.20e  %20.20e\n",$En,0.0);
#	  }
#	  ## end if SavePlot

#	}
#      }
#      ## end loop on remaining bins

#      my $lastNneg=int(-log(abs(@$EnArray[$icount+1]))/log($LambdaNRG))+2-$ztwist;

#      print " E = ".@$EnArray[$icount+1]." is on bin (".($lastNneg-1)." - ".$lastNneg.")\n";
###
###    Picking new |t|^2: negative energies
###
#      for (my $iN=50;$iN>$lastNneg;$iN--){
#        for (my $ii=1; $ii>=1; $ii--){
#	  my $irand=int(rand($Nneg));
#	  my $En=-$LambdaNRG**(-($iN+$ztwist-1+$ii*0.5));
###	  print "iN = $iN, ii=$ii New Eneg = $En t2 = ".@$BinNegative[$irand]."\n";
####	  printf(OUTFILE "%20.20e  %20.20e\n",$En,@$BinNegative[$irand]);
#	  if ($SavePlot){
#	    printf(PLOTFILE2 "%20.20e  %20.20e\n",$En,0.0);
#	    printf(PLOTFILE2 "%20.20e  %20.20e\n",$En,@$BinNegative[$irand]);
#	    printf(PLOTFILE2 "%20.20e  %20.20e\n",$En,0.0);
#	  }
#	  ## end if SavePlot


#	}
#      }
#      ## end loop on remaining bins


#    }
### end if changed sign
