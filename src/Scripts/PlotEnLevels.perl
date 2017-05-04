#!/usr/bin/perl
##
## Plots energy levels using GetEnLevels.sh
##

sub print_help(){

  print "Usage: ";
  print "$0 -f file --N0=N0 --Nf=Nf (--plot --Nlevels=Nlevels --spec) \n";
  print "  --spec : Generates an additional file with the output of GetEnLevels.sh for each N \n";
  exit(0);

}

##
## main
##

use Getopt::Long;

if ( @ARGV > 0 ){

GetOptions('f|file=s'=>\$FileName,
	   'N0=i' =>\$N0,
	   'Nf=i' =>\$Nf,
           'plot' =>\$plot,
           'Nlevels=i'=>\$NlevelsMax,
           'spec'=>\$GenSpec,
	   'h|help' =>\$help);

if (!defined($N0)||!defined($Nf)){print_help();}
if (!$FileName){$FileName="output.txt";}
if (!$NlevelsMax){$NlevelsMax=50;}

}
else{print_help();}

if ($help){print_help();}

print "File = $FileName \n";
print "N0 = $N0 \n"; 
print "Nf = $Nf \n"; 

##
## Get energy levels from GetEnLevels.sh
##

my @En_N=(); #List of lists
my $NumNs=0;
##my $NlevelsMax=50;
for ($Ns=$N0;$Ns<=$Nf;$Ns+=2){
  my @EnLevels=();
  my $nlevels=0;
  push(@EnLevels,$Ns); ## set Ns
##  $GetEnCommand="GetEnLevels -f $FileName -N $Ns | awk \'{print \$1}\'";
  $GetEnCommand="GetEnLevels.sh -f $FileName -N $Ns | awk \'{print \$1}\'";
##  print "Command = $GetEnCommand \n";
  open (GREPDATA,"$GetEnCommand |");
  while (<GREPDATA>){
    $TheLine=$_;
    chomp($TheLine);
    if ($nlevels<=$NlevelsMax){push(@EnLevels,$TheLine);}
    $nlevels++;
  }
  close(GREPDATA);
  push(@En_N,\@EnLevels);
  print "Obtaining levels: N = $Ns ; ";
  print "Nlevels = ".$nlevels." \n";


  $NumNs++;
}
##end Loop in Ns

##for ($ii=0;$ii<$NumNs;$ii++){
##  print "N = $En_N[$ii][0] En_N = $En_N[$ii][1]   $En_N[$ii][2]  $En_N[$ii][50]\n";
##}

##
##  Constructing the data file
##

open (OUTPUTFILE,"> En_vs_N.dat");

for ($ilev=1;$ilev<=$NlevelsMax;$ilev++){
  for ($ii=0;$ii<$NumNs;$ii++){
   print OUTPUTFILE "$En_N[$ii][0] $En_N[$ii][$ilev]\n";
  }
  for ($ii=$NumNs-2;$ii>=0;$ii--){
   print OUTPUTFILE "$En_N[$ii][0] $En_N[$ii][$ilev]\n";
  }
}
## Loop in levels
close(OUTPUTFILE);

##
## Generate gnuplot
##
if ($plot){
print "Plotting... \n";
open (GPLOT,"| gnuplot -persist") or die "No Gnuplot here \n";
$GnuPlotCommand.="set ylabel \"Levels\" \n set xlabel \"N\" \n";
$GnuPlotCommand.="set xrange[$N0:*] \n";
$GnuPlotCommand.="plot \"En_vs_N.dat\" with lp \n";
print GPLOT $GnuPlotCommand;
close(GPLOT);
print "... done.\n";
} # end if plot

if (defined($GenSpec)){
  print "Creating EnLevels.dat ... \n"; 
  open (OUTPUTFILE,"> EnLevels.dat");
  for ($Ns=$N0;$Ns<=$Nf;$Ns+=2){
    $RawGetEnCommand="GetEnLevels.sh -f $FileName -N $Ns ";
    open (GREPDATA,"$RawGetEnCommand |");
     while (<GREPDATA>){
     $TheLine=$_;
     chomp($TheLine);
     print OUTPUTFILE "$TheLine \n";
     }
    close(GREPDATA);
  } # end for
  close(OUTPUTFILE);
  print "... done creating EnLevels.dat  \n"; 
} # end if spec


