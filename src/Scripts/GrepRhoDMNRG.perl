#!/usr/bin/perl
##
##  Greps DM-NRG spectral function data from output*DMNRG*.txt in several directories and 
##  Plots it 
##
sub Print_Help(){
  print "Usage: \n"; 
  print "$0 -c choiceRho -i nodei -f nodeF (--6col)\n";
  print " (in development) \n";
  print "  Greps RhoEven(Odd,Interpol) from output*_DM_NRG_* in directories ./run[i]. \n";
  print "  use choiceRho as: \n";
  print "    interp  -- greps RhoInterpol (default) \n";
  print "    even    -- greps RhoEven \n";
  print "    odd     -- greps RhoOdd \n";
  print "  --6col : saves in \%f \%f \%f \%f \%f \%d format (compatible with ./ConducT code) \n";
  exit(0);
}

use Getopt::Long;

##
##  Get Parameters
##
$InputFileName="input_nrg.dat";
$choiceRho="interp";
if ( @ARGV > 0 ) {
GetOptions('c=s'=>\$choiceRho,
           'i|node0=i' => \$nodenum0,
           'f|nodeF=i' => \$nodenumF,
           '6col' => \$sixcol,
           'h|help' => \$help);
}
else {Print_Help();}

for ($choiceRho){
 if (/interp/){
  $GrepPattern="RhoInterpol";
 }
 elsif (/even/) {
  $GrepPattern="RhoEven";
 }
 elsif (/odd/) {
  $GrepPattern="RhoOdd";
 }
 else{"-c $choiceRho: Invalid option. \n";Print_Help();}
}
$nodeStep=1;
if ($help){Print_Help();}



## Setting OpName: 
$OpName=$GrepPattern;
$OpName =~ s/=//;

for ($nodenum=$nodenum0; $nodenum<= $nodenumF; $nodenum+=$nodeStep){

  $nodedirname='./run'.$nodenum;
  chdir($nodedirname) or die "Could not change dir \n";

  my @MatchingFiles=glob("out*DM*NRG*");
  my $ifilechoice=0;
  my $NoFiles=$#MatchingFiles+1;
  print "There are $NoFiles out*DM*NRG* files in dir $nodedirname: \n";
  my $icount=0;
  foreach $FileName (@MatchingFiles) {
    print "$icount  -- $MatchingFiles[$icount] \n";
   
    $GrepCommand="cat $MatchingFiles[$icount] | grep $GrepPattern";
    print " Command: $GrepCommand \n";
    $OutFileName=$MatchingFiles[$icount];
    $OutFileName =~ s/output/$GrepPattern/;
    $OutFileName =~ s/NRG_Code_//;
    $OutFileName =~ s/\.txt/\.dat/;
    print " OutFileName : $OutFileName \n";

    open (GREPDATA,"$GrepCommand |");
    open (OUTFILE,"> $OutFileName");
    while (<GREPDATA>){
       $TheLine=$_;
       chomp($TheLine);
       @OutCols=split(/\s+/,$TheLine);
#     print $OutCols[1]."   ".$OutCols[3]."\n";
       if (defined($sixcol)){
         print OUTFILE $OutCols[1]."   ".$OutCols[3]."   ".$OutCols[3]."   ".$OutCols[3]."   ".$OutCols[3]."  0 \n";
       } else {
         print OUTFILE $OutCols[1]."   ".$OutCols[3]."\n";
       }
    }
    close (OUTFILE);
    close(GREPDATA);
#    system("GPlotFiles -f $OutFileName"); 

    $icount++;
  }
  # loop in matching files
  chdir('../');
}
# loop in dirs
