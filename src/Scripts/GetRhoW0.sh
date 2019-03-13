#!/bin/bash
##
##  Script to extract the rho(w=0) values from the pre-saved 
##  spectral function data in a directory 
## 
## 
##
##
## ./Conductance needs the following files to run:
##  - ThisCodePars.dat
##  - rho_0_0_OmegaRhow.dat 
##  -   OR rho_0_0_OmegaRhow_zEQAVG.00.dat (z-trick used) 
##  - input_nrg.dat
##  - lanc.in (if NRG uses it)
##
##  To be used with ./ChangeNRGparams.perl

print_help()
{
  echo "Usage $0 -d (dir) [ -M Mtemp] "
  echo " -d dir : local dir where rho_0_0 are. Example: -d ./SpecDensFiles)"
##  echo " -f : first time: get correct links etc. "
##  echo " -r : Use FixRhoFile.perl if rho_0_0_OmegaRhow[].dat is in the old format (no integer at the end of the line) "
  echo " -M Mtemp : Finite temperature calculations. Will look for rhoMtemp[Mtemp]_0_0_OmegaRhow* files."
  echo "   Right now, using e2 in lanc.in as the running parameter "
  echo "   Please edit script to parse the correct parameter "
##  echo " -c : does NOT use DM-NRG (runs ConducT only)"
  exit 1
}


conductonly=0
Mtemp=1000
FixRhoFile=0
while getopts "d:p: M: h" option; do
  case $option in
    d) DirName=$OPTARG ;;
    p) ChoiceParam=$OPTARG ;;
    M) Mtemp=$OPTARG ;;
    h) print_help ;;
  esac
done
##echo "There were $# arguments"
if [ $# -lt 2 ]; then
  print_help
fi
if [ ! -d $DirName ]; then
  echo "Directory $DirName not found "
  print_help
fi



CurrentDir=`pwd`;
##cd $DirName;


if [ $Mtemp -lt 1000 ]; then
  FilesInDir=`ls ${DirName}/rhoMtemp${Mtemp}_0_0_OmegaRhow*`;
  betabar=0.727; 
  Lambda=2.5;
  DN=`perl -e "print 0.5*(1+1/$Lambda)*($Lambda**(-($Mtemp-1)/2))"`; 
  Temp=`perl -e "print $DN/$betabar"`; 
##  echo $Temp;
##  echo $Temp"  "$Mtemp > Temp_Mtemp.dat;
else
  FilesInDir=`ls ${DirName}/rho_0_0_OmegaRhow*`;
fi

##FilesInDir=`ls ${DirName}/rho_0_0_OmegaRhow_zAVG_e2-0.0738977.dat`;
#echo "Rho files in directory $DirName:" 
#echo $FilesInDir

ChoiceParam=1

for ThisFile in $FilesInDir; do

  DirNameNoDotsSlashes=`echo $DirName | sed "s/\.//" | sed "s/\///"`

  ThisFile=`echo $ThisFile | sed "s/\.\/$DirNameNoDotsSlashes\/\(.*\)/\1/"`;

#  echo "File: $ThisFile"
#  echo "File with Dir: $DirName/$ThisFile"


  if [ $Mtemp -lt 1000 ]; then
    ThisParam=`echo $ThisFile | sed "s/rhoMtemp${Mtemp}_0_0_OmegaRhow\(.*\)e2\(.*\).dat/\2/"`
  else
   ThisParam=`echo $ThisFile | sed "s/rho_0_0_OmegaRhow\(.*\)e2\(.*\).dat/\2/"`
  fi


  ThisRho0=`grep 4\.1076435 $DirName/$ThisFile | tail -n 1 | awk '{print \$2}'`

  echo "$ThisParam $ThisRho0" 


done

#  if [ $Mtemp -lt 1000 ]; then
#    cat output_Conductance_e2* | grep GoverG0 | awk '{print $3 "  "  $15}' | sort -g > Conductance_vse2_Mtemp${Mtemp}.dat
#  else
#    cat output_Conductance_e2* | grep GoverG0 | awk '{print $3 "  "  $15}' | sort -g > Conductance_vse2.dat
#  fi
##cd $CurrentDir;
