#!/bin/bash
##
## Move SuscepImp EntropyImp somewhere
##
print_help()
{
  echo "Usage $0 -d dir -N name"
  exit 1
}

while getopts "d:N:f: h" option; do
  case $option in
    d) dir1=./run$OPTARG
       nodir=$OPTARG ;;
    N) extraname=$OPTARG ;;
    h) print_help ;;
  esac
done

if ( ( [ -z $dir1 ] )|| ( [ -z $extraname ] ) ); then
  print_help
fi

echo "Dir : $dir1"
echo "Name : $extraname "

echo "Dir$nodir " 
extension=` ls $dir1/output*Dir$nodir.txt | sed "s/\(.*\)output\(.*\)Dir$nodir.txt/\2/"`
echo "extension: " $extension
echo
echo "Execute the Following commands?"
 ls $dir1/output*Dir$nodir.txt | sed "s/\(.*\)output\(.*\)Dir$nodir.txt/mv -i \1output\2Dir$nodir.txt .\/Data\/Phonon\/Outputs\/output\2$extraname.txt/"
echo
CommandMoveNRGin="mv -i $dir1/NRG_in.txt $TWOCHDIR/Runs/Data/Phonon/CtrolFiles/NRG_in$extension$extraname.txt"
echo $CommandMoveNRGin
echo
ls $dir1/[S,E]*Imp*.dat | sed "s/\(.*\)\/\(.*\)Imp\(.*\).dat/mv -i & .\/Data\/Phonon\/Thermo\/\2Imp\3$extension$extraname.dat/"
echo "?"
read ExecYN
if ( [ "$ExecYN" == "y" ] ); then
 ls $dir1/output*Dir$nodir.txt | sed "s/\(.*\)output\(.*\)Dir$nodir.txt/mv -i \1output\2Dir$nodir.txt .\/Data\/Phonon\/Outputs\/output\2$extraname.txt/" | sh
$CommandMoveNRGin  
ls $dir1/[S,E]*Imp*.dat | sed "s/\(.*\)\/\(.*\)Imp\(.*\).dat/mv -i & .\/Data\/Phonon\/Thermo\/\2Imp\3$extension$extraname.dat/" | sh
fi

#ls $dir1/SuscepImp*
#ls $dir1/EntropyImp*
#ls $dir1/output*
