#!/bin/bash
##
## $1 - dirI
## $2 - dirF
## $3 - which temp
print_help()
{
  echo "Usage $0 -i dir0 -f dirF -n nTemp" 
  exit 1 
}
while getopts "i:f:n:  h" option; do
  case $option in
    i) dir0=$OPTARG ;;
    f) dirF=$OPTARG ;;
    n) nT=$OPTARG ;;
    h) print_help ;;
  esac
done
if ( ( [ -z $dir0 ] )||( [ -z $dirF ] )||( [ -z $nT ] ) ); then
  print_help
fi



for ii in `seq $dir0 $dirF` 
do 
 cd ./run${ii}
 U=`head -n 3 input_nrg.dat | tail -n 1`
 Gamma=`head -n 4 input_nrg.dat | tail -n 1`
 ed=`head -n 5 input_nrg.dat | tail -n 1`
## cp -p GT_OmegaRhow.dat ../Data/ImpSolvers/GvsT/GvsT_U${U}_Gamma_${Gamma}_ed${ed}.dat 
 Temp=`head -n $nT GT_OmegaRhow.dat | tail -n 1 | awk '{print $1}'`
 GT=`head -n $nT  GT_OmegaRhow.dat | tail -n 1 | awk '{print $2}'` 
 echo $ed $GT $Temp $U $Gamma
 cd ../;
done

## Use together with
## for ii in `seq 1 51`; do ./ProcGTFiles.sh -i 130 -f 210 -n $ii > ./Data/ImpSolvers/GTvsEd/GTvsed_U0.25_Gamma0.025_nTemp${ii}.dat; done
##
