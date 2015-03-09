#!/bin/bash
##
## Simple script to Move data files
##
print_help()
{
  echo "Usage: $0 -i dir0 -f dirf -a alpha_beg -s alpha_step[ed_step] -N Nph -e ed_beg (-l lambda )"
  echo "  Runs \"MoveNRGresults\" in dir0 through dirf "
  echo "  -E : loops over ed"
  exit 1
}

Nph=9
ed=-0.25
alpha=0.0
alphastep=0.0
lambda=0.0
edcalc=0
precision=3
while getopts "i:f:a:s:N:e:l: E h" option; do
  case $option in
    i) nodir0=$OPTARG ;;
    f) nodirF=$OPTARG ;;
    a) alpha=$OPTARG ;;
    s) alphastep=$OPTARG ;;
    N) Nph=$OPTARG ;;
    e) ed=$OPTARG ;;
    E) edcalc=1 ;;
    l) lambda=$OPTARG ;;
    p) precision=$OPTARG ;;
    h) print_help ;;
  esac
done

if [ $# -lt 4 ]; then
  print_help
fi

for ii in `seq $nodir0 $nodirF`; do
# alpha=`printf %5.3f $alpha`
 alpha=`printf %5.4f $alpha`
# alpha=`printf %6.5f $alpha`
 ed=`printf %4.2f $ed`
 FullName="lambda"$lambda"_alpha"$alpha"_ed"$ed"_Nph"$Nph"_2000"

 echo $ii
 echo $FullName

 ./MoveNRGresults.sh -d $ii -N $FullName

 if [ $edcalc -eq 1 ]; then
   ed=`echo $ed+$alphastep | bc -l`
 else
   alpha=`echo $alpha+$alphastep | bc -l`
 fi 

done 
