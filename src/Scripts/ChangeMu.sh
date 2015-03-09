#!/bin/bash
print_help()
{
  echo "Usage $0 -i dir0 -f dirF -m mu0 -s mustep (-v V0)" 
  exit 1 
}

v0="1.0"
while getopts "i:f:m:s:v:  h" option; do
  case $option in
    i) dir0=$OPTARG ;;
    f) dirF=$OPTARG ;;
    m) mu0=$OPTARG ;;
    v) v0=$OPTARG ;;
    s) mustep=$OPTARG ;;
    h) print_help ;;
  esac
done
if ( ( [ -z $dir0 ] )||( [ -z $dirF ] ) ); then
  print_help
fi
if ( ( [ -z $mu0 ] )||( [ -z $mustep ] ) ); then
  print_help
fi
if ( [ $dir0 -gt $dirF ] ); then
  print_help
fi

mu=$mu0 
for idir in `seq $dir0 $dirF`; do 
echo "mu = " $mu
./SplitDisorderRealiz.perl --file=param_kondo_novo -i $idir -f $idir --mu=$mu --vzero=$v0
##mu=$(echo "$mu + 0.01" | bc -l)
mu=$(echo "$mu + $mustep " | bc -l)
mu=$(echo $mu | awk '{printf "%08f", $1}')
done

./ChangeNRGParams.perl --choiceparam=9 --p0=3 --pstep=0 -i $dir0  -f $dirF 
