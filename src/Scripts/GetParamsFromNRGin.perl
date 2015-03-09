#!/usr/bin/perl
##
## Returns a string with all parameters listed in NRG_in.in in dir ./run$DirNo
##
sub GetParamsFromNRGin(){

## Input dir number
my $DirNo=shift(@_);

my $ParamString="";

##
## Open NRG_in
##

  $DirName='./run'.$DirNo;
  chdir($DirName) or die "Could not change dir \n";
##
## Change Filename
##
   open(INFILE, "NRG_in.txt") or die "Could not open the file NRG_in.txt \n";
     while(<INFILE>){
      $TheLine=$_;
      chomp($TheLine);
      if ($TheLine =~ m/ = /){
# Split line
        @ParamVal=split(/\s+=\s+/,$TheLine);
# Get rid of ALL spaces
        foreach $String (@ParamVal){
#trailing only 
#         $String =~ s/^\s+//;
#         $String =~ s/\s+$//;
         $String =~ s/\s//g;
        }
##      end foreach
#      print $ParamVal[0]." = ".$ParamVal[1]."\n";       
      $ParamString.=$ParamVal[0]."EQ".$ParamVal[1]."_";       
      }
 ## end if
     }
#     @AllFile=<INFILE>;
   close(INFILE);

  chdir('../') or die "Could not chdir back \n";

return($ParamString);

}
## Needs to return a true value!
1;

