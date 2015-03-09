#!/bin/bash
##
## Move SuscepImp EntropyImp somewhere
##
print_help()
{
  echo "Usage $0 -d dir_no -N name (-D dest_dir)"
  echo "  Copy files from ./node(dir_no) [ or ./run(dir_no)] to respective subdirs in ./Data/dest_dir"
  echo "  dest_dir=Phonon by default"  
  echo "  Options: "
  echo "   -p : no prompt"
  echo "   -F : saves files with full name + extra name "
  exit 1
}

destdir="Phonon"
noprompt=0;
fullname=0;
while getopts "d:N:D: p F h" option; do
  case $option in
    d) dir1=./run$OPTARG
       nodir=$OPTARG ;;
    N) extraname=$OPTARG ;;
    D) destdir=$OPTARG ;;
    p) noprompt=1;;
    F) fullname=1;;
    h) print_help ;;
  esac
done

if ( ( [ -z $dir1 ] )|| ( [ -z $extraname ] ) ); then
  print_help
fi

echo "Dir : $dir1"
echo "Name : $extraname "
echo " Dest Dir : $destdir " 
##
##  Check if destination dir exists otherwise create it
##
if ( !( [ -d ./Data/$destdir/ ] ) ); then
  echo "Creating ./Data/$destdir ..."
  mkdir  ./Data/$destdir/
  mkdir ./Data/$destdir/CtrolFiles/
  mkdir ./Data/$destdir/Outputs/
  mkdir ./Data/$destdir/Thermo/
  mkdir ./Data/$destdir/SpecFuncs/
  echo "...Done" 
fi


echo "Dir$nodir " 
extension=` ls $dir1/output*Dir$nodir.txt | sed "s/\(.*\)output\(.*\)Dir$nodir.txt/\2/"`
if ( [ $fullname -ne  1 ] ); then
 extension="_"
fi
echo "extension: " $extension
echo
  echo "Executing the following commands:"
  ls $dir1/output*Dir$nodir.txt | sed "s/\(.*\)output\(.*\)Dir$nodir.txt/cp -i \1output\2Dir$nodir.txt .\/Data\/$destdir\/Outputs\/output$extension$extraname.txt/"
  echo
  CommandMoveNRGin="cp -p -i $dir1/NRG_in.txt $TWOCHDIR/Runs/Data/$destdir/CtrolFiles/NRG_in$extension$extraname.txt"
  echo $CommandMoveNRGin
  echo
  ls $dir1/[S,E]*Imp*.dat | sed "s/\(.*\)\/\(.*\)Imp\(.*\).dat/cp -i & .\/Data\/$destdir\/Thermo\/\2Imp$extension$extraname.dat/"
  ls $dir1/rho_*.dat | sed "s/\(.*\)\/rho_\(.*\).dat/cp -i & .\/Data\/$destdir\/SpecFuncs\/rho_\2$extension$extraname.dat/"
if ( [ $noprompt -eq 0 ] ); then 
  echo "Are you sure?"
  read ExecYN
else
  ExecYN="y"
fi
if ( [ "$ExecYN" == "y" ] ); then
  ls $dir1/output*Dir$nodir.txt | sed "s/\(.*\)output\(.*\)Dir$nodir.txt/cp -p  -i \1output\2Dir$nodir.txt .\/Data\/$destdir\/Outputs\/output$extension$extraname.txt/" | sh
$CommandMoveNRGin  
  ls $dir1/[S,E]*Imp*.dat | sed "s/\(.*\)\/\(.*\)Imp\(.*\).dat/cp -p -i & .\/Data\/$destdir\/Thermo\/\2Imp$extension$extraname.dat/" | sh
  ls $dir1/rho_*.dat | sed "s/\(.*\)\/rho_\(.*\).dat/cp -p -i & .\/Data\/$destdir\/SpecFuncs\/rho_\2$extension$extraname.dat/" | sh
fi

#ls $dir1/SuscepImp*
#ls $dir1/EntropyImp*
#ls $dir1/output*
