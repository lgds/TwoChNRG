#!/bin/bash
##
##  Script to calculate finte-T conductances using
##  ConducT.c with DM-NRG data!
##
print_help()
{
  echo "Usage $0 -i dir 0 -f dirF"
  exit 1
}

while getopts "i:f: h" option; do
  case $option in
    i) nodir0=$OPTARG ;;
    f) nodirF=$OPTARG ;;
    h) print_help ;;
  esac
done
echo "There were $# arguments"
if [ $# -lt 4 ]; then
  print_help
fi

#./GrepRhoDMNRG.perl -c even -i 3 -f 3
#./GrepRhoDMNRG.perl -c odd -i 3 -f 3
./GrepRhoDMNRG.perl -c even -i $nodir0 -f $nodirF --6col
./GrepRhoDMNRG.perl -c odd -i $nodir0 -f $nodirF --6col

DirNumber=1
param=3 ## e2
DataCol=3

for ii in `seq $nodir0 $nodirF`; do 
  echo "entering ./run$ii"
  cd ./run$ii
  mv RhoEven*.dat RhoT_wNeven.dat
  mv RhoOdd*.dat RhoT_wNodd.dat

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

./RunLoop.sh $nodir0 $nodirF "cat GT_x_*.dat" --silent
