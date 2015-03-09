#!/bin/bash
##
##  Run DM-NRG and conductance for several temps in a directory 
##
print_help()
{
  echo "Usage $0 -i M0 -f MF -d dir (-R) (-b broad)" 
  echo " Runs DM-NRG and Conductance (-R) in ./run(dir) for temps M0 to MF "
  echo "  -b broad : runs ./DM_NRG -b broad"
  echo "  -B broadtemp : runs ./DM_NRG -B broadtemp"
  echo "  -w twindow : runs ./DM_NRG -w twindow "
  echo "  -C : runs ./DM_NRG -C (CFS method) "
  exit 1 
}

resistivity=0
RunCFS=0
broad=""
broadtemp=""
twindow=""
while getopts "i:f:d:b:B:w: R C h" option; do
  case $option in
    i) M0=$OPTARG ;;
    f) MF=$OPTARG ;;
    d) dir=$OPTARG;;
    b) broad=$OPTARG;;
    B) broadtemp=$OPTARG;;
    w) twindow=$OPTARG;;
    R) resistivity=1;;
    C) RunCFS=1;;
    h) print_help ;;
  esac
done

if ( ( [ -z $M0 ] )||( [ -z $MF ] )||( [ -z $dir ] ) ); then
  print_help
fi
if ( [ $M0 -gt $MF ] ); then
  print_help
fi

for Mtemp in `seq $M0 2 $MF`; do 
  echo " Calculating M=$Mtemp ..."
  ./CleanDirs.sh -i $dir -f $dir -d -n
  cd ./run$dir
  Command="nice ./DM_NRG -M $Mtemp"
  if !( [ -z $broad ] ); then
    Command=${Command}" -b "$broad
  fi
  if !( [ -z $broadtemp ] ); then
    Command=${Command}" -B "$broadtemp
  fi
  if !( [ -z $twindow ] ); then
    Command=${Command}" -w "$twindow
  fi
  if [ $RunCFS -eq 1 ]; then
    Command=${Command}" -C "
  fi
  echo "Command is $Command"
  $Command > output_DMNRG.txt
#  if ( [ -z $broad ] ); then
#    nice ./DM_NRG -M $Mtemp > output_DMNRG.txt
#  else
#    nice ./DM_NRG -M $Mtemp -b $broad > output_DMNRG.txt
#  fi
  if [ $resistivity -eq 1 ]; then
    nice ./Conductance -R > resist_M$Mtemp.txt
  else
    nice ./Conductance > cond_M$Mtemp.txt
  fi
  cd ../
  echo "... done"
done

