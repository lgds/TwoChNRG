#!/bin/bash
##
## Bash script to plot energy levels
##
print_help()
{
  echo "Usage $0 (-f filename) -N Nsites"
  echo " Options"
  echo "  -f filename (output.txt by default)"
  echo "  -N Nsites"
  echo "  -u : unsorted output"
  echo "  -h : prints this message"
}
#################

File="output.txt"
while getopts ":f:N: u h" option; do
  case $option in
    f) File=$OPTARG;;
    N) Nsites=$OPTARG;;
    u) NoSort=1;;
    h) print_help
       exit 1 ;;
   esac 
done

if [ $# -eq 0 ]; then
 print_help
 exit 1
fi
if [ -z "$NoSort" ]; then
  cat $File | grep "En" | grep "_N= $Nsites$" | sed 's/En = //' | sort -g 
else
  cat $File | grep "En" | grep "_N= $Nsites$" | sed 's/En = //'  
fi
