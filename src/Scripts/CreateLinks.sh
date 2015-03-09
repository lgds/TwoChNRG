#!/bin/sh
##
## Create links i dirs
##
print_help()
{
  echo "Usage: $0 -i ni -f nf"
  echo " Creates links in directories run(ni) through run(nf)"
  exit 1
}

while getopts "i:f: h" option; do
  case $option in
     h) print_help ;;
     i) ni=$OPTARG;; 
     f) nf=$OPTARG;;
  esac
done
if ( ( [ -z $ni ] )||( [ -z $nf ] ) ); then
  print_help
fi

for ii in `seq $ni $nf`; do
  if [ -d run$ii ]; then
    echo "Creating links to executable in $TWOCHDIR/Runsi/run$ii"
    ln -s $TWOCHDIR/src/Main/NRG_main ./run$ii/NRG_main 
    ln -s $TWOCHDIR/src/DM_NRG/DM_NRG ./run$ii/DM_NRG 
    ln -s $TWOCHDIR/src/DM_NRG/Conductance ./run$ii/Conductance 
##    ln -sf $TWOCHDIR/src/GetEnLevels.sh  ./run$ii/GetEnLevels
##    ln -sf $NRGANDDIR/Dats/CalcTK/CalcTK ./run$ii/CalcTK
  else
    echo "run$ii does not exist. Create first"
  fi
done
