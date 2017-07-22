#!/bin/bash
##
##  Clean dir
##
print_help()
{
  echo "Usage $0 -i dir0 -f dir1" 
  echo "Removes the following files: "
  echo " job* SLURM* Finish_* NRG_end.txt"
  echo " Data*bin Data*dat  " 
  echo " -o : removes output_*.txt *Temp*.dat SuscepImp*.dat EntropyImp*.dat Rho* as well"
  echo " -d : removes ONLY rho*.bin, output*DM*NRG*.txt output*Conductance*.txt Rho*.dat GT_x_*.dat rho*OmegaRhow*dat files"
  echo " -n : does not prompt before cleaning"
  exit 1 
}

rmoutput=0
dmnrgonly=0
noprompt=0
while getopts "i:f: d o n h" option; do
  case $option in
    i) dir0=$OPTARG ;;
    f) dirF=$OPTARG ;;
    d) dmnrgonly=1;;
    o) rmoutput=1 ;;
    n) noprompt=1 ;;
    h) print_help ;;
  esac
done

if ( ( [ -z $dir0 ] )||( [ -z $dirF ] ) ); then
  print_help
fi
if ( [ $dir0 -gt $dirF ] ); then
  print_help
fi

echo "This will remove the following files: "
if [ $dmnrgonly -eq 1 ]; then

  echo "rho*.bin output*DM*NRG*.txt output*Conductance*.txt Rho*.dat GT_x_*.dat rho_*OmegaRhow*dat"

else

  echo " Data*bin Data*dat Abasis*bin Acut*bin Mat*bin rhoDM*bin" 
  echo " job* SLURM* Finish_* NRG_end.txt "
  if [ $rmoutput -eq 1 ]; then
    echo " Also, it will remove output*.txt *Temp*.dat SuscepImp[Chain]*.dat EntropyImp[Chain]*.dat Rho*.dat rho_*OmegaRhow*dat"
  fi

fi
echo "Clean dirs from $dir0 through $dirF ?"
if [ $noprompt -eq 0 ]; then
 read optionYN
else
 optionYN="y"
fi

if [ "$optionYN" == "y" ]; then

for dir in `seq $dir0 $dirF`; do 
  echo "Cleaning ./run$dir/..."
  if [ $dmnrgonly -eq 1 ]; then
    rm  ./run$dir/rhoDM*bin
    rm ./run$dir/output*DM*NRG*.txt
    rm ./run$dir/output*Conductance*.txt  
    rm ./run$dir/output*Conductance*.dat  
    rm ./run$dir/Rho*.dat 
    rm ./run$dir/GT_x*.dat 
    rm ./run$dir/rho*OmegaRhow*dat
  else
    rm  ./run$dir/Data*bin
    rm  ./run$dir/Data*dat
    rm  ./run$dir/Acut*bin
    rm  ./run$dir/Abasis*bin
    rm  ./run$dir/rhoDM*bin
    rm  ./run$dir/rho*OmegaRhow*dat
    rm  ./run$dir/rho_*SubGap*dat
    rm  ./run$dir/Mat*bin
    rm  ./run$dir/job*
    rm  ./run$dir/SLURM*
    rm  ./run$dir/PBS*
    rm  ./run$dir/Finish_*
    rm  ./run$dir/NRG_end.txt
    rm  ./run$dir/TK_*.dat 
    rm  ./run$dir/Ndot_*.dat 

    if [ $rmoutput -eq 1 ]; then
       echo "Removing output_*.txt *Temp*.dat SuscepImp[Chain]*.dat EntropyImp[Chain]*.dat files"
       rm ./run$dir/*Temp*.dat
       rm ./run$dir/output*.txt
       rm ./run$dir/SuscepImp*.dat 
       rm ./run$dir/EntropyImp*.dat 
       rm ./run$dir/SuscepChain*.dat 
       rm ./run$dir/EntropyChain*.dat 
       rm ./run$dir/Rho*.dat 
       rm ./run$dir/Rho*.txt 
       rm ./run$dir/HybFunc.dat ./run$dir/HybDeltas.dat
       rm ./run$dir/En_vs_N.dat
       rm ./run$dir/EnLevels.dat
    fi
  fi
  echo "... done"
done

fi
