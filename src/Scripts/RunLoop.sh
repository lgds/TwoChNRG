#!/bin/bash
## Usage: (example) 
## %./RunLoop 2 12 "ls -l NRG*.txt"
##  Runs command ls -l NRG*.txt in directories "nodes2" through "node12"  
##  Note: it will store "2" in $1 and "12" in $2 and "ls -l NRG*.txt" in $3
##
if [ $# -lt 3 ]; then
  echo "Usage $0 node_i node_f \"Command\" "
  exit
fi

for np in `seq $1 $2`;
do
#
#
 if [ "$4" != "--silent" ]; then
   echo Node $np : executing $3
 fi 
 cd $TWOCHDIR/Runs/run$np 
 $3
done

