#!/bin/bash
##
##  Script to calculate finte-T conductances using
##  ConducT.c with DM-NRG data!
##
print_help()
{
  echo "Usage $0 -i dir 0 -f dirF"
  echo " -c : does NOT use DM-NRG (runs ConducT only)"
  exit 1
}

conductonly=0
while getopts "i:f: c h" option; do
  case $option in
    i) nodir0=$OPTARG ;;
    f) nodirF=$OPTARG ;;
    c) conductonly=1;;
    h) print_help ;;
  esac
done
echo "There were $# arguments"
if [ $# -lt 4 ]; then
  print_help
fi

if [ $conductonly -eq 0 ]; then
##  ./GrepRhoDMNRG.perl -c even -i $nodir0 -f $nodirF --6col
## Use Interpolated data instead of Neven data as input
  ./GrepRhoDMNRG.perl -c interpol -i $nodir0 -f $nodirF --6col
  ./GrepRhoDMNRG.perl -c odd -i $nodir0 -f $nodirF --6col
fi

DirNumber=1
param=3 ## e2
DataCol=3

for ii in `seq $nodir0 $nodirF`; do 
  echo "entering ./run$ii"
  cd ./run$ii
  if [ $conductonly -eq 0 ]; then
##    mv RhoEven*.dat RhoT_wNeven.dat
    mv RhoInterpol*.dat RhoT_wNeven.dat
    mv RhoOdd*.dat RhoT_wNodd.dat
  fi

  if test $DirNumber -lt 10;
  then
    DirNumberCond="0$DirNumber"
##    echo "number = " $DirNumberCond
  else
    DirNumberCond="$DirNumber"
##    echo "number = " $DirNumberCond
  fi
  ./ConducT $param $DirNumberCond $DataCol
  cd ../
  DirNumber=$[$DirNumber+1]
done

echo "Result:"
./RunLoop.sh $nodir0 $nodirF "cat GT_x_*.dat" --silent
echo "Saving in GT.dat"
./RunLoop.sh $nodir0 $nodirF "cat GT_x_*.dat" --silent > GT.dat
