#!/bin/bash
##
##
##
print_help()
{
  echo "Usage $0 -i dir0 -f dirF -m M0 -M MF (-R) (-b broad)"
  echo " Runs ./RunGvsT.sh with the above parameters in dirs ./run(dir)"
  echo "  using SLURM " 
  echo " ./RunGvsT.sh: Runs DM-NRG and Conductance (-R) in ./run(dir) for temps M0 to MF "
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
while getopts "i:f:m:M:b:B:w: R C h" option; do
  case $option in
    i) dir0=$OPTARG ;;
    f) dirF=$OPTARG ;;
    m) M0=$OPTARG ;;
    M) MF=$OPTARG ;;
    b) broad=$OPTARG;;
    B) broadtemp=$OPTARG;;
    w) twindow=$OPTARG;;
    R) resistivity=1;;
    C) RunCFS=1;;
    h) print_help ;;
  esac
done

if ( ( [ -z $M0 ] )||( [ -z $MF ] )||( [ -z $dir0 ] )||( [ -z $dirF ] ) ); then
  print_help
fi
if ( [ $M0 -gt $MF ] ); then
  print_help
fi


for dir in `seq $dir0 $dirF`; do
  Command="./RunGvsT.sh -i $M0 -f $MF -d $dir"
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
  if [ $resistivity -eq 1 ]; then
    Command=${Command}" -R "
  fi
  echo "Command is $Command"
  SLURMrun.sh "$Command" $TWOCHDIR/Runs/ $dir
done
