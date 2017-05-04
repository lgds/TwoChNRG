#!/bin/bash
##
##  Script to calculate the T=0 conductance with pre-saved 
##  spectral function data in a directory 
## 
##  Needs a symbolic link to the Conductance code.
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
  echo "Usage $0 -d (dir where rho_0_0 are) -p (choiceparam)"
##  echo " -c : does NOT use DM-NRG (runs ConducT only)"
  exit 1
}


conductonly=0
while getopts "d:p: h" option; do
  case $option in
    d) DirName=$OPTARG ;;
    p) ChoiceParam=$OPTARG ;;
    h) print_help ;;
  esac
done
echo "There were $# arguments"
if [ $# -lt 2 ]; then
  print_help
fi
if [ ! -f ChangeNRGparams.perl ]; then
  echo " Needs ./ChangeNRGparams.perl in the same dir !" 
  print_help
fi
if [ ! -d $DirName ]; then
  echo "Directory $DirName not found "
  print_help
fi

CurrentDir=`pwd`;
##cd $DirName;
FilesInDir=`ls ${DirName}/rho_0_0_OmegaRhow*`;
##FilesInDir=`ls ${DirName}/rho_0_0_OmegaRhow_zAVG_e2-0.0738977.dat`;
##echo "Rho files in directory $DirName:" 
##echo $FilesInDir

rm -f output_Conductance_e2*.txt

for ThisFile in $FilesInDir; do

ThisFile=`echo $ThisFile | sed "s/$DirName\/\(.*\)/\1/"`;

echo "File: $ThisFile"

ThisParam=`echo $ThisFile | sed "s/rho_0_0_OmegaRhow\(.*\)e2\(.*\).dat/\2/"`

echo "Param: $ThisParam"

perl ChangeNRGParams.perl -i 99 -f 99 --file="lanc.in" --choiceparam=$ChoiceParam --p0=$ThisParam --pstep=0 --indir

cp -f ${DirName}/$ThisFile ./rho_0_0_OmegaRhow_zEQAVG.00.dat

echo "Calculating Conductance..."
./Conductance > output_Conductance_e2${ThisParam}.txt
echo "... done!"


done

##cd $CurrentDir;