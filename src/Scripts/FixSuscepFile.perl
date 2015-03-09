#!/usr/bin/perl
##
## Fix Suscep, Entropy Files
##
###############################
sub Print_Help() {
  print " Usage $0 --susfile=SuscepFilename \n"; 
  print " Options:  \n";
  print "    --chainfile=ChainFilename (SuscepChain.dat by default) \n";
  print "    --entropy          : Entropy files used \n";
  print "    --addconst=Const   : Add constant to Imp values (ex 0.25) \n";
  print "    --split            : Split file into odd and even (assumes first N=0) \n";
  print "    --overwrite        : Overwrite original file, saving a copy marked \"OLD\" \n";
  print "    --nochain          : Does not read chain file.  \n";
  print "    --shiftchain=nshift : Shift chain values by nshift(>0). Does not read chain file.  \n";
  exit;

}
###############################

sub ReadSusFile{

my $TempArray=shift(@_);
my $SusArray=shift(@_);
my $SusChainArray=shift(@_);
my $SusImpArray=shift(@_);
my $FileName=shift(@_);
my $nt=0;

@$TempArray=();
@$SusArray=();
@$SusChainArray=();
@$SusImpArray=();

 open(INFILE,$FileName);
   while (<INFILE>){
     $TheLine = $_;
     chomp($TheLine);
###     ($auxtemp,$auxsuschainpimp,$auxsuschain,$auxsusimp)=split(/ /,$TheLine);
     ($auxtemp,$auxsuschainpimp,$auxsuschain,$auxsusimp)=split(/ +/,$TheLine);
     push(@$TempArray,$auxtemp);
     push(@$SusArray,$auxsuschainpimp);
     push(@$SusChainArray,$auxsuschain);
     push(@$SusImpArray,$auxsusimp);
     $nt++;
   }
 close(INFILE);

 return($nt);
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
########################
##   Main program     ##
########################

use Getopt::Long;

$PrefixFile="Suscep";
$Nshift=0;

if ( @ARGV > 0 ){
GetOptions('f|susfile=s'=>\$SuscepFilename,
           'chainfile=s'=>\$ChainFilename,
           'entropy'=>\$EntropyFile,
           'addconst=f'=>\$AddConst,
           'split'=>\$Split,
           'overwrite'=>\$OvrtFile,
           'nochain'=>\$NoChain,
           'shiftchain=i'=>\$Nshift,
           'h|help'=>\$help);
}
##else{Print_Help();}
if (defined($help)){Print_Help();}
if (defined($EntropyFile)){$PrefixFile="Entropy";}


if (defined($OvrtFile)){$OvrtFile=1;}else{$OvrtFile=0;}
if (!defined($SuscepFilename)){$SuscepFilename=$PrefixFile."Imp2Ch_25_726.dat";}
if (!defined($ChainFilename)){$ChainFilename=$PrefixFile."Chain2Ch.dat";}
$ReadChain=1;
if (defined($AddConst)||defined($NoChain)||defined($Split)){$ReadChain=0;}
if ($Nshift>0){$ReadChain=0;}


## 4 columns: Temp SusChainPImp SusChain SusImp
##
@Temps=();
@SusChainPImp=();
@SusChain=();
@SusImp=();
@TempArray=();
my $ntemps=0;
my $ntempsChain=0;

##########################
##   Read SuscepFile    ##
##########################

$ntemps=ReadSusFile(\@Temps,\@SusChainPImp,\@SusChain,\@SusImp,$SuscepFilename);

#New
for ($ishift=0;$ishift<$Nshift;$ishift++){
   shift(@SusChain);
   pop(@Temps);
   pop(@SusChainPImp);
   pop(@SusImp);
   $ntemps--;
}

############################
##   Read SuscepChain     ##
############################

if ($ReadChain){
 print " Reading chain file $ChainFilename \n";
 $ntempsChain=ReadSusFile(\@AuxArray,\@SusChain,\@AuxArray,\@AuxArray,$ChainFilename);


 if ($ntemps!=$ntempsChain){
   print "Mismatch in Chain and Imp files: $ntemps neq $ntempsChain\n";
   if ($ntempsChain>$ntemps){print "Using first $ntemps values in Chain file \n";}
   else{ print "Missing chain sites. Exiting...\n"; exit;}
 }
}else{$ntempsChain=$ntemps;}
#End if ReadChain
############################
##   Combine the two      ##
############################
print " Ntemps = $ntemps = ".($#Temps+1)." , NChain = $ntempsChain \n";
for ($ii=0;$ii<=$#Temps;$ii++){
## print "ii = $ii Temp: $Temps[$ii] Sus : $SusChainPImp[$ii] SusChain $SusChain[$ii] \n";
## printf("ii = %d Temp: %5.3e  Sus : %5.5e SusChain: %5.5e \n",$ii,$Temps[$ii],$SusChainPImp[$ii],$SusChain[$ii]);
## printf("ii = %d Temp: %5.3e  Susi-SusChain= %5.5e SusImp: %5.5e \n",$ii,$Temps[$ii],$SusChainPImp[$ii]-$SusChain[$ii],$SusImp[$ii]);
 if (defined($AddConst)){$SusChainPImp[$ii]+=$AddConst;} 
 $SusImp[$ii]=$SusChainPImp[$ii]-$SusChain[$ii];
}

############################
##   Save in a New File   ##
############################

##if (defined($OvrtFile)){
if ($OvrtFile==1){
  $SansDatName=$SuscepFilename;
  $SansDatName =~ s/.dat//;
  $OldFilename=$SansDatName."_OLD.dat";
  if (defined($Split)){
   print "Saving to ".$SansDatName."_Even.dat ...\n";
   SaveSusFile(\@Temps,\@SusChainPImp,\@SusChain,\@SusImp,$SansDatName."_Even.dat","e");
   print "Saving to ".$SansDatName."_Odd.dat ...\n";
   SaveSusFile(\@Temps,\@SusChainPImp,\@SusChain,\@SusImp,$SansDatName."_Odd.dat","o");
  }
  else{
   print "Copying $SuscepFilename to  $OldFilename  and rewriting ...\n";
   system("cp $SuscepFilename $OldFilename");
   SaveSusFile(\@Temps,\@SusChainPImp,\@SusChain,\@SusImp,$SuscepFilename);
  }
}
else
{
  if (defined($Split)){
    print "Saving to ".$PrefixFile."Out_Even.dat ...\n";
    SaveSusFile(\@Temps,\@SusChainPImp,\@SusChain,\@SusImp,$PrefixFile."Out_Even.dat","e");
    print "Saving to ".$PrefixFile."Out_Odd.dat ...\n";
    SaveSusFile(\@Temps,\@SusChainPImp,\@SusChain,\@SusImp,$PrefixFile."Out_Odd.dat","o");
  }
  else{
    print "Saving to ".$PrefixFile."Out.dat ...\n";
    SaveSusFile(\@Temps,\@SusChainPImp,\@SusChain,\@SusImp,$PrefixFile."Out.dat");
  }
}

print "...Done!\n";
#######
# END #
#######
