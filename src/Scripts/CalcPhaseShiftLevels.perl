#!/usr/bin/perl
##
## Plots energy levels using GetEnLevels.sh
##

sub print_help(){

  print "Usage: ";
  print "$0 -f file --N0=N0 --Nf=Nf --E1chain=f (--plot --Nlevels=Nlevels) \n";
  print " Options: \n";
  print "  -v, --verbose : display messages \n";
  print "  --nosave      : does not save file \n";
  print "  --outname     : saves in \"phase_shift_[outname].dat\" \n";
  exit(0);

}
#######################

sub dEqual{
## Checks if A=B (can be real numbers)

my $NoArgs=@_;
my $NumA=shift(@_);
my $NumB=shift(@_);
my $Tolerance=1.0E-15;
if ($NoArgs>2){$Tolerance=shift(@_);}

##print "NoArgs = $NoArgs, tol = $Tolerance \n";
if ($NumA==$NumB){return(1);}
else{
  if ( abs($NumA-$NumB)<$Tolerance){return(1);}
  else{return(0);}
}

}


#######################

sub GetSinglePLev{

  my $HoleEnLevels=shift(@_);
  my $SinglePLev=shift(@_);

## No of levels: 
  my $Nlevels=$#HoleEnLevels+1;
  my $PrevLevel=-100;
  my $CurrentLevel=0;
  my $idegen=1;

###
###  1 - Filter Degeneracies
###

  for (my $ii=0; $ii<50;$ii++){
   $CurrentLevel=@$HoleEnLevels[$ii];
##   if ((abs($CurrentLevel-$PrevLevel)<0.01)||
##       (abs($CurrentLevel)<0.0001)){
   if (abs($CurrentLevel-$PrevLevel)<0.01){
     $idegen++;
##     print "Found GS/deg level: deg = $idegen \n" ;
   }
   else{$idegen=1;push(@$SinglePLev,$CurrentLevel);}
   $PrevLevel=$CurrentLevel;
  }

## Watch the syntax!! Makes you crazy...
  my $NSPlevels=$#$SinglePLev+1;

##  foreach $en (@$SinglePLev){
##    print "En = $en \n";
##  }

###
###  2 - Filter multiples of same level
###

  my $BaseLevel=-100;
  for (my $ii=0; $ii<$NSPlevels;$ii++){
    $BaseLevel=@$SinglePLev[$ii];

    for (my $jj=$ii+1; $jj<$NSPlevels;$jj++){

      $CurrentLevel=@$SinglePLev[$jj];

      my $Ratio=$CurrentLevel/$BaseLevel;
      if ( (($Ratio-int($Ratio))<0.001)||
           (($Ratio-int($Ratio))>0.999) )
	{
## Remove/Insert elements from list: 
##  splice(@Nums, pos, no_el_to_replace, New list)
	splice(@$SinglePLev,$jj,1,() );
	$NSPlevels--;
	}
##      else {print "No! \n";}
    }
### end loop in jj

  }
## end loop in ii


###
###  3 - Filter more complicated combinations of different levels
###  (still to come)
###

  return(0);

}

#######################
##
## main
##

use Getopt::Long;
use constant PI    => 4*atan2(1, 1);

if ( @ARGV > 0 ){

GetOptions('f|file=s'=>\$FileName,
	   'N0=i' =>\$N0,
	   'Nf=i' =>\$Nf,
           'plot' =>\$plot,
           'Nlevels=i'=>\$NlevelsMax,
	   'E1chain=f'=>\$E1chain,
           'v|verbose'=>\$verbose,
	   'nosave' => \$nosave,
	   'outname=s' =>\$outname,
	   'h|help' =>\$help);

if ((!$N0)||(!$Nf)){print_help();}
if (!$FileName){$FileName="output2ch.txt";}
if (!$NlevelsMax){$NlevelsMax=60;}

if (!defined($E1chain)){$E1chain=0.5;}

}
else{print_help();}

if ($help){print_help();}

print "File = $FileName \n";
print "N0 = $N0 \n";
print "Nf = $Nf \n";

my @En_N=(); #List of lists
my $NumNs=0;
##my $NlevelsMax=50;

##
##  Get Extension
##
if (defined($outname)){$Ext=$outname;}
else{
 $Ext=$FileName;
 $Ext =~ s/^output_//;
 $Ext =~ s/.txt$//;
}
if (!defined($nosave)){open(SAVEPHASE,"> phase_shift_$Ext.dat");}

for ($Ns=$N0;$Ns<=$Nf;$Ns+=2){
  my @HoleEnLevels=();
  my @ElecEnLevels=();
  my @HolePm1EnLevels=();
  my @ElecPm1EnLevels=();
  my @SPholeLevs=();
  my @SPelecLevs=();
  my @SPholePm1Levs=();
  my @SPelecPm1Levs=();
  my $nlevels=0;
  my @LineData=();
  my $Nqns=0;
  my $Qgs=0;
  my $Sgs=0;
  my $Pgs=0;

  my $NumGS=0;
  my $TotNlev=0;
  my @QGS=();
  my @SGS=();
  my @PGS=();


  print "Obtaining levels: N = $Ns ; ";


##  push(@HoleEnLevels,$Ns); ## set Ns
##  $GetEnCommand="GetEnLevels -f $FileName -N $Ns | grep \"| -1  0\" | awk \'{print \$1}\'";
  $GetEnCommand="GetEnLevels -f $FileName -N $Ns ";
##  print "Command = $GetEnCommand \n";
  open (GREPDATA,"$GetEnCommand |");
  while (<GREPDATA>){
    $TheLine=$_;
    chomp($TheLine);
    @LineData=split(/ +/,$TheLine);
###   Get No QNs:
    $Nqns=@LineData-4;
### Set "parity=1" if onle Q,S
    if ($Nqns==2){$LineData[4]=1;} 

##
##  Find Ground state(s) for different parities
##
    if (dEqual($LineData[0],0.0,0.0001)){
      push (@QGS,$LineData[2]);
      if ($Nqns>1){push(@SGS,$LineData[3]);}
      if ($Nqns>2){push(@PGS,$LineData[4]);}
      $NumGS++;
      if (defined($verbose)){
	print "GS found: \n";
	print " QNs : ";
	for (my $ii=0;$ii<$Nqns;$ii++){
	  print $LineData[2+$ii]."  ";
	}
	print "\n";
	print "Num Gs = $NumGS \n";
	print "Qgs : ".$QGS[$NumGS-1]." \n";
	print "Sgs : ".$SGS[$NumGS-1]." \n";
	print "Pgs : ".$PGS[$NumGS-1]." \n";
      } ## verbose
    }
    ## end if E_i =0
##    if (dEqual($LineData[0],0.0)){
##       $Qgs=$LineData[2];
##       if ($Nqns>1){$Sgs=$LineData[3];}
##       if ($Nqns>2){$Pgs=$LineData[4];}
##       if (defined($verbose)){
##         print "GS found: Q = $Qgs \n";
##         print " QNs : ";
##         for (my $ii=0;$ii<$Nqns;$ii++){
##           print $LineData[2+$ii]."  ";
##         }
##         print "\n";
##       } ## verbose
##     }
    $nlevels++;
  }
  close(GREPDATA);
  $TotNlev=$nlevels;
  print "Nlevels = ".$TotNlev." \n";
  print "N GS = ".$NumGS." \n";

##
##  Once the GS have been identified, get the excitation spectrum
##

  $nlevels=0;
  open (GREPDATA,"$GetEnCommand |");
  while (<GREPDATA>){
    $TheLine=$_;
    chomp($TheLine);
    @LineData=split(/ +/,$TheLine);

    ## Loop in GS
    for (my $igs=0;$igs<$NumGS;$igs++){
##
##  Hole excitations (P=1)
##

      if ( (!dEqual($LineData[0],0.0,0.0001))&&
	   (dEqual($LineData[2],$QGS[$igs]-1.0))&&
	   ( (dEqual($LineData[3],abs($SGS[$igs]-0.5)))||
	     (dEqual($LineData[3],abs($SGS[$igs]+0.5))) )&&
           (dEqual($LineData[4],1)) ){
        if ($nlevels<=$NlevelsMax){push(@HoleEnLevels,$LineData[0]);}
      }

##
##  Electron excitations (P=1)
##
      if ( (!dEqual($LineData[0],0.0,0.0001))&&
	   (dEqual($LineData[2],$QGS[$igs]+1.0))&&
	   ( (dEqual($LineData[3],abs($SGS[$igs]-0.5)))||
	     (dEqual($LineData[3],abs($SGS[$igs]+0.5))) )&&
           (dEqual($LineData[4],1)) ){
        if ($nlevels<=$NlevelsMax){push(@ElecEnLevels,$LineData[0]);}
      }


      if ($Nqns==3){
##
##  Hole excitations (P=-1)
##
        if ( (!dEqual($LineData[0],0.0,0.0001))&&
	     (dEqual($LineData[2],$QGS[$igs]-1.0))&&
             ( (dEqual($LineData[3],abs($SGS[$igs]-0.5)))||
	       (dEqual($LineData[3],abs($SGS[$igs]+0.5))) )&&
             (dEqual($LineData[4],-1)) ){
          if ($nlevels<=$NlevelsMax){push(@HolePm1EnLevels,$LineData[0]);}
        }

##
##  Electron excitations (P=-1)
##
        if ( (!dEqual($LineData[0],0.0,0.0001))&&
	     (dEqual($LineData[2],$QGS[$igs]+1.0))&&
             ( (dEqual($LineData[3],abs($SGS[$igs]-0.5)))||
	       (dEqual($LineData[3],abs($SGS[$igs]+0.5))) )&&
             (dEqual($LineData[4],-1)) ){
          if ($nlevels<=$NlevelsMax){push(@ElecPm1EnLevels,$LineData[0]);}
        }
      }
## end if nqns=3

    }
##end loop in GS
    $nlevels++;
  }
## end while grep data
  close(GREPDATA);


  GetSinglePLev(\@HoleEnLevels,\@SPholeLevs);
  GetSinglePLev(\@ElecEnLevels,\@SPelecLevs);
  if ($Nqns==3){
   GetSinglePLev(\@HolePm1EnLevels,\@SPholePm1Levs);
   GetSinglePLev(\@ElecPm1EnLevels,\@SPelecPm1Levs);
  }
  if (defined($verbose)){
     print "Hole single-particle P=1 levels \n";
   
     for (my $ii=0;$ii<6;$ii++){
       $en=$SPholeLevs[$ii];
       print "Eh_SP = $en , delta/pi = ".$en/(2.0*$E1chain)." \n";
     }
     print "Electron single-particle P=1 levels \n";
   
     for (my $ii=0;$ii<6;$ii++){
       $en=$SPelecLevs[$ii];
       print "Ee_SP = $en \n";
     }
   
     print "Hole single-particle P=-1 levels \n";
   
     for (my $ii=0;$ii<6;$ii++){
       $en=$SPholePm1Levs[$ii];
       print "EhP=-1_SP = $en , delta/pi = ".$en/(2.0*$E1chain)." \n";
     }
     print "Electron single-particle P=-1 levels \n";
   
     for (my $ii=0;$ii<6;$ii++){
       $en=$SPelecPm1Levs[$ii];
       print "EeP=-1_SP = $en \n";
     }
   
  } ##verbose
##
## Calc phase shift
##

  my $PhaseShift=0.0;
  my $Conductance=0.0;
  my $ps_even=0.0;
  my $ps_odd=0.0;
  my $UseHoleLevs=0;  
  if (($SPelecLevs[0]-$SPholeLevs[0])>0.001){
   $UseHoleLevs=1;
   if (defined($verbose)){print "Using hole phase shift. \n"}
  } 
 
  $ps_even=$SPelecLevs[0]/(2.0*$E1chain);
  if ($UseHoleLevs){$ps_even=$SPholeLevs[0]/(2.0*$E1chain);}
  if ($Nqns==3){
    $ps_odd=$SPelecPm1Levs[0]/(2.0*$E1chain);
    if ($UseHoleLevs){$ps_odd=$SPholePm1Levs[0]/(2.0*$E1chain);}
##    $PhaseShift=abs(($SPelecLevs[0]-$SPelecPm1Levs[0])/(2.0*$E1chain));
    $PhaseShift=abs($ps_even-$ps_odd);
    $Conductance=(sin(PI*$PhaseShift))**2;
  }


  print "ps: $Ns even: $ps_even  odd : $ps_odd  |even-odd|: $PhaseShift \n";
  printf ("N : %i e: %7.5f o: %7.5f |o-e| : %7.5f G: %7.5f \n", $Ns, $ps_even,  $ps_odd , $PhaseShift, $Conductance);
if (!defined($nosave)){  print SAVEPHASE "$Ns  $ps_even $ps_odd  $PhaseShift $Conductance \n";}

##  push(@En_N,\@HoleEnLevels);

  $NumNs++;
}
##end Loop in Ns
if (!defined($nosave)){close(SAVEPHASE);}
