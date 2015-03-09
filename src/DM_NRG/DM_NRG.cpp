
#include <iostream>
//#include <iomanip>
//#include <fstream>

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include <boost/timer.hpp>

#include <vector>
#include <cmath>
#include <cstring>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "SpecFuncClass.hpp"
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

  CNRGbasisarray SingleSite;

  CNRGmatrix** OpArrayN;

  CNRGmatrix* RhoN;
  vector<double> ParamsTemp;
  double betabar=0.727;
    double DM,TM;
  // 0 - betabar; 1 - TM; 2 - DM ?
  double Temp; // In units of the bandwidth (fixed!)
               // Temp/DN = TempBar = 1/betabar 
  int Mtemp=1001;

  double twindow=2.0;

  double broadtemp=0.8; // Good for CFS

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
  int Nw=0; // No of omegas in each shell

  // Read code parameters, including z_twist
  ThisCode.ReadGenPars(true); 

  double bbroad=0.5*log(ThisCode.Lambda);

  // Command-line: Set Temp
  DM_NRG_CommandLineRead(argc,argv,Mtemp,betabar,twindow,broadtemp,bbroad,UseCFS,Nw);

  // Set Temp
  cout << " Mtemp = " << Mtemp 
       << " betabar = " << betabar 
       << " twindow = " << twindow 
       << " dBroad = " << bbroad 
       << " dBroadTemp = " << broadtemp << endl
       << " NomegasEachShell = " << Nw << " (if =0 then use standard interpolation) "
       << endl; 
  //betabar=0.727;
  Temp=0.0;  // Real Temperature
  //ParamsTemp.push_back(betabar); // Need to fix this today...
  if (Mtemp<ThisCode.Nsitesmax-1) {
    ThisCode.Nsitesmax=Mtemp+1;
    DM=CalcDN(ThisCode.Lambda,Mtemp);
    Temp=DM/betabar;
  }
  else{ // OK, Let me try this. Instead of 0.0, "Temp" will be DM(Nsites+200)/betabar.
    //Temp=0.0;
    Temp=CalcDN(ThisCode.Lambda,ThisCode.Nsitesmax+199)/betabar;
    DM=CalcDN(ThisCode.Lambda,ThisCode.Nsitesmax-1);
  }
  cout << " Temp = "<<Temp
       << " DMtemp = " << DM
       << endl;
  //DM=CalcDN(ThisCode.Lambda,NshellMax);
  ParamsTemp.push_back(Temp/DM);
  // Save in Temp_Mtemp.dat
  ofstream OutFile;
  strcpy(arqname,"Temp_Mtemp.dat");
  OutFile.open(arqname);
  OutFile << Temp << " " << Mtemp << endl;
  OutFile.close();


  cout << "Lambda    = " << ThisCode.Lambda << endl;
  cout << "Nsites0   = " << ThisCode.Nsites0 << endl;
  cout << "Nsitesmax = " << ThisCode.Nsitesmax << endl;
  cout << "NFermiOps = " << ThisCode.NopsSaved << endl;
  cout << "Ext arq   = " << ThisCode.SaveArraysFileName << endl;
  cout << "Symmetry  = " << ThisCode.SymNo << endl;
  if ( (ThisCode.SymNo != 2)&&(ThisCode.SymNo != 3) ){
    cout << " DM-NRG: Only 1chQSz and 1chQ symmetries supported at this point. Testing QS symmetry now..." << endl;
  }
  if (ThisCode.totalS){cout << " SU(2) symmetry detected. " << endl;}
  cout << "Oliveira z = " << ThisCode.chain.z_twist << endl;


  NshellMax=ThisCode.Nsitesmax-1;
  NFermiOps= ThisCode.NopsSaved;

  // April 2010
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
  RhoN = new CNRGmatrix [NshellMax+1];


  // Read last Nshell

  cout << " Reading bin Files for NshellMax = " << NshellMax << endl;
  AcutN[NshellMax].Nshell=NshellMax;
  AbasisNp1.Nshell=NshellMax;
  //ThisCode.ReadArrays("test"); // ReadsInto pAcutN
  // Not too handy...

  sprintf(CNsites,"%d",NshellMax);
  //strcpy(ext,"test_N");
  strcpy(ext,ThisCode.SaveArraysFileName);
  strcat(ext,"_N");
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
  cout << " NshellMax Files read "<< endl;

  // Set density matrix at the LAST NRG iteration

  // Either CalcRhoN AND save it OR read it from file
  // Read RhoNmax
//    sprintf(CNmat,"%d",NshellMax);
//    strcpy(arqname,"rhoDM_Nmax");
//    strcat(arqname,CNmat);
//    strcat(arqname,"_");
  strcpy(arqname,"rhoDM_");
  strcat(arqname,ext);

  if (ThisCode.CheckFileExists(arqname)){
    cout << " DM_NRG: Found file " << arqname << endl;
    RhoN[NshellMax].ReadBin(arqname);
  }else{
    DM_NRG_SetRhoNmax(ParamsTemp,&AcutN[NshellMax],&RhoN[NshellMax]);
    RhoN[NshellMax].SaveBin(arqname);
  }


  cout << " Rho Nshell = " << NshellMax << endl;
  //if (NshellMax==48) 
  //RhoN[NshellMax].PrintAllBlocks();
  // Debugging
  double qnums[2];
  int iBl=0;
  
  // Calculate reduced density matrices

  //for (int Nshell=NshellMax-1;Nshell>=0;Nshell--){
  for (int Nshell=NshellMax-1;Nshell>=ThisCode.Nsites0;Nshell--){

    cout << "DM-NRG: working on Nshell = " << Nshell << endl;

    AcutN[Nshell].ClearAll();
    AbasisN.ClearAll();
    AcutN[Nshell].Nshell=Nshell;
    AbasisN.Nshell=Nshell;
    //ThisCode.ReadArrays("test"); // ReadsInto pAcutN
    // Not too handy...

    sprintf(CNsites,"%d",Nshell);
    //strcpy(ext,"test_N");
    strcpy(ext,ThisCode.SaveArraysFileName);
    strcat(ext,"_N");
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

    MyTime.restart();

    // Either CalcRhoN AND save it OR read it from file

    // Read RhoN
    strcpy(arqname,"rhoDM_");
    strcat(arqname,ext);

    if (ThisCode.CheckFileExists(arqname)){
      cout << " DM_NRG: Found file " << arqname << endl;
      RhoN[Nshell].ReadBin(arqname);
    }else{
      if (ThisCode.totalS){
	DM_NRG_CalcRhoN_withSU2(&AcutN[Nshell],&AcutN[Nshell+1],
			  &AbasisNp1,&SingleSite,
			  &RhoN[Nshell],&RhoN[Nshell+1]);

      }
      else{
	DM_NRG_CalcRhoN(&AcutN[Nshell],&AcutN[Nshell+1],&AbasisNp1,
			&RhoN[Nshell],
			&RhoN[Nshell+1]);
      }
      // if SU(2) symmetry
      RhoN[Nshell].SaveBin(arqname);
    }

    time_elapsed=MyTime.elapsed();
    cout << " Rho Nshell = " << Nshell 
	 << " completed in " << time_elapsed << " secs " << endl;

    // Debugging
     if (Nshell==48){
       //RhoN[Nshell].PrintAllBlocks();
//     qnums[0]=-2.000;
//     qnums[1]=1.000;
//     iBl=RhoN[Nshell].GetBlockFromQNumbers(qnums);
//     RhoN[Nshell].PrintMatBlock(iBl,iBl);
//     qnums[0]=2.000;
//     qnums[1]=1.000;
//     iBl=RhoN[Nshell].GetBlockFromQNumbers(qnums);
//     RhoN[Nshell].PrintMatBlock(iBl,iBl);
     }

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
    
//     cout << " Op1 : " << endl;
//     OpArrayN[0][Nshell].PrintAllBlocks();
//     cout << " Op2 : " << endl;
//     OpArrayN[1][Nshell].PrintAllBlocks();

    // Debug
//     if ((Nshell==0)||(Nshell==1)){
//       cout << " DM-NRG: Abasis(N=0): " << endl;
//       AbasisN.PrintBasisAll();
//       cout << " Op1 (Nshell="<< Nshell<<"): " << endl;
//       OpArrayN[0][Nshell].PrintAllBlocks();
//     }
    // end debug


  }
  // end loop in Nshell


  // Save RhoN matrices

  // Read RhoN matrices

  // Given AcutN, rhoN and the Operators, calculate the spectral density

  // Calculate rho_0_0 and rho_Costi. Good for debugging
//   DM_NRG_CalcSpecFuncs(&ThisCode,AcutN,RhoN,OpArrayN,0,0);
//   DM_NRG_CalcSpecFuncs(&ThisCode,AcutN,RhoN,OpArrayN,1,1);
//   DM_NRG_CalcSpecFuncs(&ThisCode,AcutN,RhoN,OpArrayN,2,2);


   CSpecFunction spec1;

   spec1.Lambda=ThisCode.Lambda;
   spec1.z_twist=ThisCode.chain.z_twist;
   spec1.NshellMax=ThisCode.Nsitesmax;
   spec1.NshellMin=ThisCode.Nsites0;
   spec1.AcutN=AcutN;
   spec1.RhoN=RhoN;

   // Spec Dens
   spec1.BDelta=BroadDelta;
   //spec1.dBroad=0.5*log(ThisCode.Lambda);
   spec1.dBroad=bbroad;

   // Finite Temp stuff (need to test this!)
   if (UseCFS==1) spec1.BDeltaTemp=LorentzDeltaAnders;
   else spec1.BDeltaTemp=LorentzDelta;

   //spec1.Temp=0.0;
   //spec1.Mtemp=1000;
   spec1.Temp=Temp;
   spec1.Mtemp=Mtemp;
   spec1.Betabar=betabar;

   spec1.TwindowFac=twindow;
   spec1.dBroadTemp=broadtemp;


   // Calculate ALL spectral functions!!
   for (int iop=0; iop<NFermiOps; iop++){
     if (UseCFS==1){
       DM_NRG_CalcSpecFunc_ij(&spec1,OpArrayN,iop,iop,UseCFS,Nw);
     } else {
       for (int jop=0; jop<NFermiOps; jop++){
	 DM_NRG_CalcSpecFunc_ij(&spec1,OpArrayN,iop,jop,UseCFS,Nw);
       }
     }
     // end if CFS
   }
   // end loop in Fermi Ops

   //int iop=0;
   //int jop=0;
   //CFS_CalcSpecFunc_ij(&spec1,OpArrayN,iop,jop);
   


  delete[] RhoN;
  delete[] AcutN;
  for (int iop=NFermiOps-1;iop>=0;iop--){
    delete[] OpArrayN[iop];
  }
  delete[] OpArrayN;

}
// end MAIN
