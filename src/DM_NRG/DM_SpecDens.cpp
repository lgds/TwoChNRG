
#include <iostream>
//#include <iomanip>
//#include <fstream>

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include <vector>
#include <cmath>
#include <cstring>

//#include <unistd.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
//#include "NRG_main.hpp"
//#include "TwoChQS.hpp"


#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif


int main (int argc, char* argv[]){


  CNRGCodeHandler ThisCode;

  CNRGbasisarray Acut;
  CNRGbasisarray Abasis;

  CNRGmatrix* MatArray;

  ThisCode.NumNRGmats=0;
  ThisCode.pAcut=&Acut;
  ThisCode.pAbasis=&Abasis;
  //ThisCode.MatArray=MatArray;


  ThisCode.ReadGenPars("");

  cout << "Lambda = " << ThisCode.Lambda << endl;
  cout << "Nsitesmax = " << ThisCode.Nsitesmax << endl;

  Acut.Nshell=0;
  Abasis.Nshell=0;
  ThisCode.ReadArrays("test");

  Acut.PrintAll();


}
// end MAIN
