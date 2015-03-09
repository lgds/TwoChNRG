#!/usr/bin/perl
##
## Subroutine GetYgivenX.perl
## Parameters
##   Filename: data file with X Y columns
##   X       : X value 
##
sub GetYgivenX{

  my $Filename=shift(@_);
  my $Xvalue=shift(@_);

  my $ThisX=0; my $PrevX=10;
  my $ThisY=0; my $PrevY=10;
  my $Yfound=0.0; 
  my $Alin=0; my $Blin=0;
 
  open(DATAFILE,"cat $Filename | sort -g |") or die "Cant open $Filename \n";
  while (<DATAFILE>) {
     my $TheLine = $_;
     chomp($TheLine);
     ($ThisX, $ThisY)=split(/ +/, $TheLine);
     if ( ($Xvalue>=$PrevX)&&($Xvalue<=$ThisX) )
      {
#        print "Found : $PrevX < $Xvalue < $ThisX \n";
        $Alin=($ThisY-$PrevY)/($ThisX-$PrevX);
        $Blin=$ThisY-$Alin*$ThisX;
        $Yfound=$Alin*$Xvalue+$Blin;
      }
     $PrevX=$ThisX;
     $PrevY=$ThisY;
  }
  close(DATAFILE);

  return($Yfound);
}
# Returns a true value
1;

#########
# Test  #
#########
#my $Xval=3.0;
#my $Filename="./Data/STM/Curve1.dat";
#my $Yvalue=GetYgivenX($Filename,$Xval);
#print "X = $Xval, Y = $Yvalue \n";

