#!/bin/sh
##
## Create dir run$ii
##
print_help()
{
  echo "Usage: $0 -i ni -f nf"
  echo " Creates directories run(ni) through run(nf)"
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
    echo "run$ii exists"
  else
    echo "Creating run$ii in $TWOCHDIR/Runs"
    mkdir run$ii
##    ln -sf $TWOCHDIR/src/TwoChQS/TwoChQS ./run$ii/TwoChQS
##    ln -sf $TWOCHDIR/src/OneChQS/OneChQS ./run$ii/OneChQS
    ln -sf $TWOCHDIR/src/Main/NRG_main ./run$ii/NRG_main 
    ln -sf $TWOCHDIR/src/DM_NRG/DM_NRG ./run$ii/DM_NRG 
    ln -sf $TWOCHDIR/src/DM_NRG/Conductance ./run$ii/Conductance 
##    ln -sf $TWOCHDIR/Runs/run2/SuscepChain2Ch.dat ./run$ii/SuscepChain2Ch.dat 
##    ln -sf $TWOCHDIR/Runs/run2/EntropyChain2Ch.dat ./run$ii/EntropyChain2Ch.dat 
##    ln -sf $TWOCHDIR/Runs/run2/SuscepChain1Ch.dat ./run$ii/SuscepChain1Ch.dat 
##    cp ./run3/nrg_input_TwoCh.dat ./run$ii
##    cp ./run3/Input_Phonon.dat  ./run$ii
##    cp ./run1/nrg_input_OneChQS.dat ./run$ii
    cp ./run3/input_nrg.dat ./run$ii
    cp ./run3/lanc.in ./run$ii
#    ln -sf $TWOCHDIR/src/GetEnLevels.sh  ./run$ii/GetEnLevels
    ln -sf $NRGANDDIR/Dats/CalcTK/CalcTK ./run$ii/CalcTK
  fi
done
