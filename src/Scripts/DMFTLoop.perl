#!/usr/bin/perl
##
## DMFT Loop: start from HybFunc.dat, calculate rho_0_0_OmegaRhow.dat
##  and then do:
##
##  Gamma(w)=phi.V^2.rho(w)
##
## and run it again.
##
###############################

sub CreateHybfuncFile {

  use constant PI    => 4*atan2(1, 1);

  my $V0=shift(@_);
  my $SpecFuncFileName=shift(@_);
  my $HybFuncFileName=shift(@_);
  my $facPrevious=shift(@_);
  my $Diff=0.0;

  my @LineData=();
  my @EnArray=();
  my @RhoArray=();
  my @HybEnArray=();
  my @PrevHybArray=();


##
## Read Spectral function
##
  open(INFILE,$SpecFuncFileName);
  while (<INFILE>){
     $TheLine = $_;
    chomp($TheLine);
    @LineData=split(/ +/,$TheLine);
    my $En=$LineData[0];
    my $rho=$LineData[1];
    push(@EnArray,$En);
    push(@RhoArray,$rho);
   # end if -1<=en<=1
  }
  close(INFILE);
##
## Read Previous Hybfunction
##
  open(INFILE,$HybFuncFileName);
  while (<INFILE>){
     $TheLine = $_;
    chomp($TheLine);
    @LineData=split(/ +/,$TheLine);
    my $En=$LineData[0];
    my $PrevGammaw=$LineData[1];
    push(@HybEnArray,$En);
    push(@PrevHybArray,$PrevGammaw);
  }
  close(INFILE);

  $Nens=@EnArray-1;
  $NensHyb=@HybEnArray-1;

  print "NenSpec = $Nens , NenHyb= $NensHyb , Factor= $facPrevious\n";

  if ($Nens!=$NensHyb){
    print " Problems: energies in the files do not match! NenSpec neq NenHyb \n";
    exit;
  }

  print "Opening $HybFuncFileName \n";
  open(OUTFILE,"> ".$HybFuncFileName);
  for (my $ii=0;$ii<=$Nens;$ii++){
   my $En=@EnArray[$ii];
   my $rho=@RhoArray[$ii];
   my $PrevHyb=@PrevHybArray[$ii];
   my $NewHyb=PI*($V0**2)*$rho;
   $Diff+=abs($NewHyb-$PrevHyb);
   ## Mix old and new
   $NewHyb=(1.0-$facPrevious)*PI*($V0**2)*$rho + $facPrevious*$PrevHyb;
   print OUTFILE $En."  ".$NewHyb."\n";
##   print $En."  ".$NewHyb."\n";
  }
  close(OUTFILE);

 return($Diff);

}

###############################
sub Print_Help() {
  print " Usage $0 (-f factor) (-M Mtemp) > output.txt \n";
  print "  DMFT loop: Reads V0 from input_nrg.dat constructs Delta(w)=V0^2.Gamma(w) \n";  
  print " Options:  \n";
  print "  -f factor : (default factor=0) uses (1-factor)NewHyb+factor.OldHyb as input to the next iteration. \n";
  print "  -M Mtemp  : Finite-temperature calculation. Mtemp is the NRG site corresponding to the scale T_M.\n";
  print "  -n NwEachBin : (default Nw=1) Calculates Nw omega values for each NRG bin (note: this might produce artifacts in the band). \n";
  exit;
}
###############################



########################
##   Main program     ##
########################

use Getopt::Long;

my $Factor=0.0;
my $Nw=1;
if ( @ARGV > 0 ){
GetOptions('f=f'=>\$Factor,
           'M=i'=>\$Mtemp,
           'n=i'=>\$Nw,
           'h|help'=>\$help);
}
##else{Print_Help();}
if (defined($help)){Print_Help();}
## Reading V0 from input_nrg.dat
##if ( (!defined($V0)) ){Print_Help();}

use constant PI    => 4*atan2(1, 1);

##
## DMFT Loop
##
if (-e "input_nrg.dat" ){
  print "Reading V0 from input_nrg.dat. \n";
}else{
  print "Ops!! input_nrg.dat not found!! exiting...\n";
  exit;
} 
my $V0=`cat input_nrg.dat | head -n 4 | tail -n 1`;
print "V0=".$V0."\n";

my $HybFuncFileName="HybFunc.dat";
my $HybFuncNm1FileName="HybFunc_Nm1.dat";
my $HybFuncNm2FileName="HybFunc_Nm2.dat";

my $SpecFuncFileName="rho_0_0_OmegaRhow.dat";

## Clean bin files 
system("rm -f *.bin \n");

##
##  Initialize: first HybFunc.dat gess and first NRG run
##

### Next steps: read these from file
##my $Lambda=2.5;
my $Lambda=`cat input_nrg.dat | head -n 6 | tail -n 1`;
chomp($Lambda);
my $Nsites=`cat input_nrg.dat | head -n 1`;
chomp($Nsites);
my $factorWN=$Lambda**(1.25);
## Finite temp
if (defined($Mtemp)){$Nsites=$Mtemp+1};

print "Reading NRG paramaters: Lambda= $Lambda and Nsites = $Nsites \n";

if (-e $HybFuncFileName ){
  print "File HybFunc.dat exits. Using it as a first guess. \n";
}else{
##  system("cp -p HybFunc0_TEST.dat HybFunc.dat");
## NEED A BETTER INITIAL FILE, WITH THE CORRECT ENERGIES
##
  print "Generating Hybfunc.dat \n";
  open(OUTFILE,"> ".$HybFuncFileName);
  for (my $Nsh=0; $Nsh<$Nsites; $Nsh++){
##     my $en=-$factorWN*0.5*(1.0 + (1.0/$Lambda))*$Lambda**(-($Nsh - 1)/2.0);
     my $DN=0.5*(1.0 + (1.0/$Lambda))*$Lambda**(-($Nsh - 1)/2.0);
     ## Fine structure in each bin (negative omega):
     ## DN(Nsh+1)=DN(Nsh)*Lambda^-1/2
     ## Thus, we want omega=-DN(Nsh)*Lambda^-(iw/2Nw) iw=0,1,Nw-1
     for (my $iw=0;$iw<$Nw;$iw++){
       my $NwFactor=$Lambda**(-$iw/(2*$Nw));
       my $en=-$factorWN*$NwFactor*$DN;
       my $Gammaw=PI*($V0**2)*(1.0/2.0);
       print "en= ".$en." Gammaw= ".$Gammaw."\n";
       print OUTFILE $en." ".$Gammaw."\n";
     }
     # end for in each bin
  }
  # end for
  for (my $Nsh=$Nsites-1; $Nsh>=0; $Nsh--){
##     my $en=$factorWN*0.5*(1.0 + (1.0/$Lambda))*$Lambda**(-($Nsh - 1)/2.0);
     my $DN=0.5*(1.0 + (1.0/$Lambda))*$Lambda**(-($Nsh - 1)/2.0);
     ## Fine structure in each bin (positive omega):
     ## DN(Nsh-1)*Lambda^-1/2=DN(Nsh)
     ## Thus, we want omega=+DN(Nsh)*Lambda^(-iw/2Nw) iw=Nw-1,Nw-2,...,0
     for (my $iw=$Nw-1;$iw>=0;$iw--){
       my $NwFactor=$Lambda**(-$iw/(2*$Nw));
       my $en=$factorWN*$NwFactor*$DN;
       my $Gammaw=PI*($V0**2)*(1.0/2.0);
       print "en= ".$en." Gammaw= ".$Gammaw."\n";
       print OUTFILE $en." ".$Gammaw."\n";
     }
     # end for in each bin
  }
  # end for
  close(OUTFILE);
############
}
# end if

## Set HybFuncNm1
system("cp -f ".$HybFuncFileName."  ".$HybFuncNm1FileName."\n");

##
##  NRG calculation commands (feel free to edit).
##

$NRGCommand="nice ./NRG_main -m Anderson -b FromFile > out_NRG_Anderson.txt \n";

## Setting DM_NRG command
$DMNRGCommand="nice ./DM_NRG";

#if (defined($Mtemp)){$DMNRGCommand="nice ./DM_NRG -M ".$Mtemp." > out_DMNRG.dat \n";}
#else{$DMNRGCommand="nice ./DM_NRG > out_DMNRG.dat \n";}

## Non-zero temperature
if (defined($Mtemp)){$DMNRGCommand.=" -M ".$Mtemp;}
## -n Nw option
if ($Nw>1){$DMNRGCommand.=" -n ".$Nw;}
$DMNRGCommand.=" > out_DMNRG.dat \n";

##
##
## DMFT Loop
##
##

my $Nloop=0;

my $ThisCommand="";
my $Tolerance=0.0001;
my $Diff=100;

my $Nmax=20;  # max number of iterations

while ( ($Nloop<$Nmax)&&($Diff>$Tolerance) ){

##
##
##  Remove .bin files, move SpecFuncFileName and run NRG
##
  $ThisCommand="rm -f *.bin \n";
  $ThisCommand.="mv -f ".$SpecFuncFileName." rhoDMFT_Nm1.dat \n";
##

  print "Loop: $Nloop : Running NRG... \n";

  system($ThisCommand);
  system($NRGCommand);

  print "... done.  Now running DM-NRG... \n";

  system($DMNRGCommand);

  if (!(-e $SpecFuncFileName )){
    print "Ops, something went wrong: File $SpecFuncFileName not found. \n";
    exit;
  }

  $ThisCommand="cp -f ".$HybFuncNm1FileName."  ".$HybFuncNm2FileName."\n";
  $ThisCommand.="cp -f ".$HybFuncFileName."  ".$HybFuncNm1FileName."\n";
##  print $ThisCommand."\n";
  system($ThisCommand);
##  $Diff=CreateHybfuncFile($V0,$SpecFuncFileName,$HybFuncFileName,0.0);
  $Diff=CreateHybfuncFile($V0,$SpecFuncFileName,$HybFuncFileName,$Factor);
  print "Diff = $Diff \n";

  $Nloop++;
}
#end loop
