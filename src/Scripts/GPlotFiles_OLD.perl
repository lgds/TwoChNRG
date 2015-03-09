#!/usr/bin/perl
##
##  Subroutine to plot files in Gnuplot
##
sub Print_Help(){

  print "Usage: \n";
  print " $0 -f|--file Filepattern (--xcol=xc --ycol==yc) \n";
  print "   Options: \n";
  print "   -f|--file      : Matching files (mandatory). Ex: Sus\\*.dat \n";
  print "   --xcol=xc, --ycol=yc : gnuplots using specified cols ( u xc:yc). xc=1,yc=2 by default \n";
  print "   --setlog=(x,y)       : sets log (x,y) \n";
  print "   --every=(2,3,...)    : sets every (1 by default) \n";
  print "   --xmax=xmax          : set xrange[*:xmax] \n";
  print "   --xmin=x0,--xmax=x1  : set xrange[x0:x1] \n";
  print "   --ymin=y0,--ymax=y1  : set yrange[y0:y1] \n";
  print "   --save               : saves as GPlot.plt\n";
  print "   --plotps             : saves as GPlot.ps\n";
  exit(0);
}

##           ##
## Main Code ##
##           ##

use Getopt::Long;

$xc=1;
$yc=2;
$every=1;
$xmin="*";
$xmax="*";
$ymin="*";
$ymax="*";
if ( @ARGV>0 ){
   GetOptions('f|file=s'=>\$FilePattern,
              'xcol=s'=>\$xc,
              'ycol=s'=>\$yc,
              'setlog=s' => \$setlog,
##              'every=i' => \$every,
              'every=s' => \$every,
              'save' => \$save,
              'xmin=s' => \$xmin,
              'xmax=s' => \$xmax,
              'ymin=s' => \$ymin,
              'ymax=s' => \$ymax,
              'plotps' => \$plotps,
              'h|help' =>\$help);
}
else{Print_Help();}

if ($help){Print_Help();}
if (!defined($FilePattern)){Print_Help();}

##
##  Get Files: curly brackets work!!
##
##my @MatchingFiles=glob("$FilePattern"); #Incomplete
my @MatchingFiles=glob($FilePattern);
#my @MatchingFiles=<${FilePattern}>;
## Test
##my $NoFiles=@MatchingFiles;
##for (my $if=0;$if<$NoFiles;$if++){print $MatchingFiles[$if]."\n";}
##exit(0);
#
open (GP, "|gnuplot -persist") or die "no Gnuplot";
##
##  Gnuplots 
##
if (defined($plotps)){
  $GnuplotCommand.="set terminal postscript color enhanced \n";
  $GnuplotCommand.="set output \"GPlot.ps\"\n";
}
if (defined($setlog)){$GnuplotCommand.="set log $setlog \n";}
$GnuplotCommand.="set xrange[$xmin:$xmax] \n";
$GnuplotCommand.="set yrange[$ymin:$ymax] \n";
$GnuplotCommand.="plot ";
  foreach $Filename (@MatchingFiles)
  {
    $GnuplotCommand.="\"$Filename\" every $every u $xc:$yc with lp,"
  }
  $GnuplotCommand =~ s/,$//;
  $GnuplotCommand.="\n";
  if (defined($plotps)){
    $GnuplotCommand.="set terminal x11\n";
    $GnuplotCommand.="set output\n";
  }
  if (defined($save)){$GnuplotCommand.="save \"GPlot.plt\" \n";}
  print "$GnuplotCommand \n";
  print GP "$GnuplotCommand \n";

close(GP);
