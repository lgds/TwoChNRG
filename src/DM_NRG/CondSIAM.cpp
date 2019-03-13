
#include <iostream>


#include <boost/timer.hpp>

#include <vector>
#include <cmath>
#include <cstring>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "SpecFuncClass.hpp"
//#include "CondSIAM.hpp"
#include "NRGOpMatRules.hpp"


#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif


int main (int argc, char* argv[]){

  CNRGCodeHandler ThisCode;

  CNRGbasisarray* AcutN;
  CNRGbasisarray AbasisN;

  CNRGbasisarray SingleSite;

  CNRGmatrix** OpArrayN;

//   vector<double> ParamsTemp;
  double betabar=0.727;
  double DM,TM;
  // 0 - betabar; 1 - TM; 2 - DM ?
  double Temp; // In units of the bandwidth (fixed!)
               // Temp/DN = TempBar = 1/betabar 

  char CNsites[8];
  char CNmat[8];
  char arqname[32];
  char ext[32];

  int NshellMax=3;
  int NFermiOps=2;

  // Check time
  boost::timer MyTime;
  double time_elapsed;


  int UseCFS=0;
  int Nw=1; // No of omegas in each shell

  int c;
  opterr = 0;

  char cmNtemp[8],cbetabar[8];

  while ((c = getopt (argc, argv, "h b: n:")) != -1)
    switch (c){
    case 'b':
      strcpy(cbetabar,optarg);
      betabar=atof(cbetabar);
      break;
    case 'n':
      strcpy(cmNtemp,optarg);
      Nw = atoi(cmNtemp);
      break;
    case 'h':
      cout << "Usage: " << argv[0] << " (-b betabar -n NtempsEachN)" << endl;
      cout << " Default values: betabar=0.727 ; NtempsEachN=1 " << endl;
      // 	cout << " Additional options: " << endl;
      // 	cout << " -C : uses CFS method " << endl;
      exit(0);
      break;
    case '?':
      cout << "Usage: " << argv[0] << " -b betabar -n NtempsEachN" << endl;
      cout << " Default values: betabar=0.727 ; NtempsEachN=1 " << endl;
      exit(0);
      break;
    default:
      break;
    }
  // end CommandLine



  // Read code parameters, including z_twist
  ThisCode.ReadGenPars(true); 



  // Set Temp
  cout << " betabar = " << betabar 
       << " NomegasEachShell = " << Nw << " (if =1 then use standard interpolation) "
       << endl; 


  cout << "Lambda    = " << ThisCode.Lambda << endl;
  cout << "Nsites0   = " << ThisCode.Nsites0 << endl;
  cout << "Nsitesmax = " << ThisCode.Nsitesmax << endl;
  cout << "NFermiOps = " << ThisCode.NopsSaved << endl;
  cout << "Ext arq   = " << ThisCode.SaveArraysFileName << endl;
  cout << "Symmetry  = " << ThisCode.SymNo << endl;
  cout << "Model     = " << ThisCode.ModelNo << endl;
  if ( (ThisCode.ModelNo != 0)  ){
    cout << " CalcSIAM: Only conductance in the Anderson model supported at this point." << endl;
    cout << " CalcSIAM: Proceed at your own risk ... " << endl;
  }
  if (ThisCode.totalS){cout << " SU(2) symmetry detected. " << endl;}
  cout << "Oliveira z = " << ThisCode.chain.z_twist << endl;

  ThisCode.ReadParams((char *)"input_nrg.dat",1); // Needs Gamma

  NshellMax=ThisCode.Nsitesmax-1;
  NFermiOps= ThisCode.NopsSaved;



  // TotS QNs in ThisCode
  ThisCode.SetTotS();
  // Setting SingleSite (depends only on SymNo)
  ThisCode.SetSingleSite(&SingleSite);

  // Allocate matrices: Tricky

  OpArrayN=new CNRGmatrix* [NFermiOps];
  for (int iop=0;iop<NFermiOps;iop++){
    OpArrayN[iop]=new CNRGmatrix [NshellMax+1];
  }

  AcutN = new CNRGbasisarray [NshellMax+1];


  // SIAM Conductance calculation

  // Read Acut files and calculates the conductance

  //for (int Nshell=NshellMax;Nshell>=ThisCode.Nsites0;Nshell--){
  for (int Nshell=ThisCode.Nsites0;Nshell<=NshellMax;Nshell++){

    cout << "CondSIAM: Reading files for Nshell = " << Nshell << endl;

    AcutN[Nshell].ClearAll();
    AcutN[Nshell].Nshell=Nshell;

    sprintf(CNsites,"%d",Nshell);
    //strcpy(ext,"test_N");
    strcpy(ext,ThisCode.SaveArraysFileName);
    strcat(ext,"_N");
    strcat(ext,CNsites);
    strcat(ext,".bin");

    // Read Acut
    strcpy(arqname,"Acut_");
    strcat(arqname,ext);
    AcutN[Nshell].ReadBin(arqname);

    // Read Operators 
    for (int iop=0;iop<NFermiOps;iop++){
      sprintf(CNmat,"%d",iop); // Get all SavedMatrices
      strcpy(arqname,"Mat");
      strcat(arqname,CNmat);
      strcat(arqname,"_");
      strcat(arqname,ext);
      OpArrayN[iop][Nshell].ReadBin(arqname);
    }
    // end read operators
    
  }
  // end loop in Nshell


  // Define spectral function

  CSpecFunction spec1;

  spec1.Lambda=ThisCode.Lambda;
  spec1.z_twist=ThisCode.chain.z_twist;
  spec1.NshellMax=ThisCode.Nsitesmax;
  spec1.NshellMin=ThisCode.Nsites0;
  spec1.AcutN=AcutN;

  //spec1.Temp=Temp;
  //spec1.Mtemp=Mtemp;
  spec1.Betabar=betabar;

  //spec1.TwindowFac=twindow;
  //spec1.dBroadTemp=broadtemp;

  spec1.NonDiagGF=false;
  spec1.Op1N=OpArrayN[0];

  spec1.CalcCondSIAM_ManyTs(Nw);

  // Need to multiply it by  Gamma

  // File Name
  strcpy(spec1.Name,"GT");

  // Will add OmegaRhow, z-value and .dat
  spec1.SaveOmegaRhow();

// Need to multiply it by Gamma...


// de-allocate matrices. Watch out.
  delete[] AcutN;
  for (int iop=NFermiOps-1;iop>=0;iop--){
    delete[] OpArrayN[iop];
  }
  delete[] OpArrayN;

}
// end MAIN
