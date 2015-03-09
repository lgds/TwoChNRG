
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
#include "DM_NRG.hpp"
#include "NRGOpMatRules.hpp"


#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif


int main (int argc, char* argv[]){


  CNRGCodeHandler ThisCode;

  //CNRGbasisarray AcutN;
  CNRGbasisarray* AcutN;
  CNRGbasisarray AbasisN;

  //CNRGbasisarray AcutNp1;
  CNRGbasisarray AbasisNp1;

  CNRGmatrix** OpArrayN;

  CNRGmatrix* RhoN;
  vector<double> ParamsTemp;
  double betabar, DM,TM;
  // 0 - betabar; 1 - TM; 2 - DM ?

  char CNsites[8];
  char CNmat[8];
  char arqname[32];
  char ext[32];

  int NshellMax=3;
  int NFermiOps=2;

  ThisCode.NumNRGmats=0;
  ThisCode.pAcut=&AcutN[NshellMax];
  ThisCode.pAbasis=&AbasisN;
  //ThisCode.MatArray=MatArray;

  // Set betabar
  betabar=0.727;
  ParamsTemp.push_back(betabar);

  // Tricky

  OpArrayN=new CNRGmatrix* [NFermiOps];
  for (int iop=0;iop<NFermiOps;iop++){
    OpArrayN[iop]=new CNRGmatrix [NshellMax+1];
  }

  AcutN = new CNRGbasisarray [NshellMax+1];
  RhoN = new CNRGmatrix [NshellMax+1];

  ThisCode.ReadGenPars("");

  // This for now
  ThisCode.Nsitesmax=NshellMax;

  cout << "Lambda = " << ThisCode.Lambda << endl;
  cout << "Nsitesmax = " << ThisCode.Nsitesmax << endl;

  // Read last Nshell

  AcutN[NshellMax].Nshell=NshellMax;
  AbasisNp1.Nshell=NshellMax;
  //ThisCode.ReadArrays("test"); // ReadsInto pAcutN
  // Not too handy...

  sprintf(CNsites,"%d",NshellMax);
  strcpy(ext,"test_N");
  strcat(ext,CNsites);
  strcat(ext,".bin");

  // Read Abasis
  strcpy(arqname,"Abasis_");
  strcat(arqname,ext);
  AbasisNp1.ReadBin(arqname);

  // Read Acut
  strcpy(arqname,"Acut_");
  strcat(arqname,ext);
  AcutN[NshellMax].ReadBin(arqname);

  // Read Operators 
  for (int iop=0;iop<NFermiOps;iop++){
    sprintf(CNmat,"%d",iop);
    strcpy(arqname,"Mat");
    strcat(arqname,CNmat);
    strcat(arqname,"_");
    strcat(arqname,ext);
    OpArrayN[iop][NshellMax].ReadBin(arqname);
    // Need a better way to do this but for now it will do:
    // Actually, from Mat block it works!
    //OpArrayN[iop][NshellMax].CheckForMatEl=OneChQSz_cd_check;

  }
  // end read operators

  //AcutN[NshellMax].PrintAll();

  // Set density matrix at the LAST NRG iteration

  DM_NRG_SetRhoNmax(ParamsTemp,&AcutN[NshellMax],&RhoN[NshellMax]);

  //RhoN[NshellMax].PrintAllBlocks();

  // Calculate reduced density matrices

  for (int Nshell=NshellMax-1;Nshell>=0;Nshell--){

    AcutN[Nshell].ClearAll();
    AbasisN.ClearAll();
    AcutN[Nshell].Nshell=Nshell;
    AbasisN.Nshell=Nshell;
    //ThisCode.ReadArrays("test"); // ReadsInto pAcutN
    // Not too handy...

    sprintf(CNsites,"%d",Nshell);
    strcpy(ext,"test_N");
    strcat(ext,CNsites);
    strcat(ext,".bin");

    // Read Abasis
    strcpy(arqname,"Abasis_");
    strcat(arqname,ext);
    AbasisN.ReadBin(arqname);

    // Read Acut
    strcpy(arqname,"Acut_");
    strcat(arqname,ext);
    AcutN[Nshell].ReadBin(arqname);


    // Set ChildStates in AcutN
    DM_NRG_SetChildSt(&AcutN[Nshell],&AbasisNp1);

    // Set RhoN from RhoNp1

    DM_NRG_CalcRhoN(ParamsTemp,
		    &AcutN[Nshell],&AcutN[Nshell+1],&AbasisNp1,
		    &RhoN[Nshell],
		    &RhoN[Nshell+1]);

    //RhoN[Nshell].PrintAllBlocks();

    // Update AbasisNp1
    AbasisNp1.ClearAll();
    AbasisNp1=AbasisN;

    // Read Operators 
    for (int iop=0;iop<NFermiOps;iop++){
      sprintf(CNmat,"%d",iop); // Get all SavedMatrices
      strcpy(arqname,"Mat");
      strcat(arqname,CNmat);
      strcat(arqname,"_");
      strcat(arqname,ext);
      OpArrayN[iop][Nshell].ReadBin(arqname);
      // Need a better way to do this but for now it will do:
      // Not needed!
      //OpArrayN[iop][Nshell].CheckForMatEl=OneChQSz_cd_check;

    }
    // end read operators



  }
  // end loop in Nshell



  // Given AcutN, rhoN and the Operators, calculate the spectral density


  DM_NRG_CalcSpecFuncs(&ThisCode,AcutN,RhoN,OpArrayN);



  delete[] RhoN;
  delete[] AcutN;
  for (int iop=NFermiOps-1;iop>=0;iop--){
    delete[] OpArrayN[iop];
  }
  delete[] OpArrayN;

}
// end MAIN
