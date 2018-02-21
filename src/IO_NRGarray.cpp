#include <iostream> 
#include <fstream>

#include "NRGclasses.hpp"
 
using namespace std;


////////////////////////
void SaveNRGarrayBin(CNRGarray *pArray,char arqname[]){

  ofstream FileOut(arqname, ios::out | ios::binary);

  if (!FileOut){cout << "SaveNRGArrayBin: Cannot save data in " << arqname << endl; return;}


  FileOut.write((const char *)pArray, sizeof(CNRGarray));
  FileOut.close();

  if (!FileOut.good()){cout << "SaveNRGArrayBin: error saving" << endl;}

}
// end SaveNRGArrayBin

////////////////////////


// DOES NOT WORK!!!
void ReadNRGarrayBin(CNRGarray *pArray,char arqname[]){

  ifstream FileIn(arqname, ios::in | ios::binary);

  if (!FileIn){cout << "InNRGArrayBin: Cannot open " << arqname << endl; return;}


  FileIn.read((char *)pArray, sizeof(CNRGarray));
  FileIn.close();

  if (!FileIn.good()){cout << "ReadNRGArrayBin: error reading" << endl;}

}
// end ReadNRGArrayBin

