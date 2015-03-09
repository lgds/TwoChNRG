#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>
#include <cmath>
#include <cstring>

#include <unistd.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "OneChQS.hpp"

#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif

#ifndef _CHIN_
#define _CHIN_


double chiN(int Nsites, double Lambda)
{

  double daux[3];
  daux[0]=1.0-pow( Lambda,-((double)(Nsites)+1.0) );
  daux[1]=sqrt( 1.0-pow(Lambda,-(2.0*(double)(Nsites)+1.0)) );
  daux[2]=sqrt( 1.0-pow(Lambda,-(2.0*(double)(Nsites)+3.0)) );  

  return(daux[0]/(daux[1]*daux[2]));

}

#endif


int main (int argc, char* argv[]){


  // Parameters for command-line passing (GetOpt)

  //char ModelOption[]="Anderson";
  //int ModelNo=0;

#include"ModelOpt.cpp"


  // NRG objects

  CNRGarray Aeig(2);

  CNRGbasisarray AeigCut(2);

  CNRGbasisarray Abasis(2);

  CNRGbasisarray SingleSite(3);

  // STL vector

  CNRGmatrix HN;
  CNRGmatrix Qm1fNQ;

  //CNRGmatrix MQQp1;

  CNRGmatrix* MatArray;
  int NumNRGarrays=3;
  // MatArray 0 is nd
  // MatArray 1 is cd
  // MatArray 2 is a test

  int KeepSz=0;

  double U,ed,Gamma;
  double Lambda;
  double HalfLambdaFactor;
  double Dband;
  int auxIn;
  int calcdens=0;
  vector <double> Params;
  vector <int> Indexes;

  double chi_m1,chi_N,daux[4];
  int Nsites;
  int Nsitesmax=5;
  int Nsites0=0;
  double DN=0.0;


  // Thermodynamics

  CNRGthermo Suscep;
  CNRGthermo Entropy;


  double TM=0.0;
  double betabar=0.727;
  double Temp=0.0;
  double Sus=0.0;
  double SusDot=1.0/8.0;
  vector<double> SuscepChain;

  double ndot=0.0;
  double Sz2=0.0;

  int Ncutoff=700;
  int UpdateBefCut=1;

  char arqSus[32],arqname[32];

  int ii,jj,i1,i2;

  // STL iterator:

  vector<double>::iterator diter;

  // outstream
  ofstream OutFile;
  // instream
  ifstream InFile;


  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  ////                    READ PARAMETERS                      ////
  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////


  U=0.5;
  ed=-0.5*U;
  Gamma=0.0282691;
  Lambda=2.5;


  InFile.open("nrg_input_OneChQS.dat");
  if (InFile.is_open())
    {
      InFile >> Nsitesmax;
      InFile >> Ncutoff;
      InFile >> U;
      InFile >> Gamma;
      InFile >> ed;
      InFile >> Lambda;
      InFile >> Dband;
      InFile >> auxIn;
      InFile >> UpdateBefCut;
      InFile >> calcdens;
    }
  else
    {
      cout << "can't open nrg_input_OneChQS.dat" << endl;
      exit(0);
   }

  InFile.close();

  // NRG chain
  CNRGchain chain1(Lambda,Nsitesmax);


  // New stuff
  strcpy(Suscep.ArqName,"SuscepImp1Ch_25_726.dat");
  strcpy(Suscep.ChainArqName,"SuscepChain1Ch.dat");
  Suscep.Calc=CalcSuscep;
  
  strcpy(Entropy.ArqName,"EntropyImp1Ch_25_726.dat");
  strcpy(Entropy.ChainArqName,"EntropyChain1Ch.dat");
  Entropy.Calc=CalcEntropy;

  // Check for SuscepChain.dat file. If it is there, read Temp, Sus
  // If not, change to calcdens=3, U=0, ed=0 Gamma=1
  if ( (calcdens==1)&&(ModelNo==1) )
    {
      cout << "Kondo model: Can't calculate spectral function" << endl;
      exit(0);
    }
  if (ModelNo==3) calcdens=3; // Chain calculation
  if (calcdens==2)
    {
      Suscep.CalcChain=false;
      Entropy.CalcChain=false;
      Suscep.Nsite0Chain=0; // Anderson chains
      Entropy.Nsite0Chain=0; // Anderson chains

      if ( (!FileExists(Suscep.ChainArqName))||
	   (!FileExists(Entropy.ChainArqName)) )
	{
	  cout << " Can't find chain files: " << endl 
	  << "   " << Suscep.ChainArqName << endl
	  << "   " << Entropy.ChainArqName << endl;
	  cout << " Exiting... " << endl;
	  exit(0);
	}
      else
	  cout << " Found files " 
	       << Suscep.ChainArqName << ", " 
	       << Entropy.ChainArqName 
	       << endl;



      //int ReadChainStatus=Suscep.ReadChain();
      //int ReadChainStatus2=Entropy.ReadChain();
      //if (ReadChainStatus==-1)calcdens=3;

    }
  if (calcdens==3)
    {
      Suscep.CalcChain=true;
      Entropy.CalcChain=true;
      U=0.0;ed=0.0;Gamma=0.0;
    }


  HalfLambdaFactor=0.5*(1.0+(1.0/Lambda));
  chi_m1=sqrt(2.0*Gamma/pi)/(sqrt(Lambda)*HalfLambdaFactor);

  double U_tilde=ed/(Lambda*HalfLambdaFactor);
  cout << "e1 = " << 0.5*U/(Lambda*HalfLambdaFactor) << endl;
  cout << "e2 = " << (ed+0.5*U)/(Lambda*HalfLambdaFactor) << endl;
  cout << "e3 = " << (2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor) << endl;
  cout << "chi_m1  = " << chi_m1 << endl;


  cout << " Nsitesmax    = " << Nsitesmax-1 << endl;
  cout << " Ncutoff      = " << Ncutoff << endl;
  if (strcmp(ModelOption,"Anderson")==0){
    cout << " U            = " << U << endl;
    cout << " Gamma        = " << Gamma << endl;
    cout << " ed           = " << ed << endl;
  }
  if (strcmp(ModelOption,"Kondo")==0){
    cout << " JK            = " << Gamma << endl;
  }
  cout << " Lambda       = " << Lambda << endl;
  cout << " Dband        = " << Dband << endl;
  cout << " UpdateBefCut = " << UpdateBefCut << endl;
  cout << " calcdens     = " << calcdens << endl;

  /////////////////////////////////////////////////
  OutFile.open("NRG_in.txt");
  OutFile << "Begin NRG calculation. Model :" << ModelOption << endl;
  OutFile << " Nsitesmax    = " << Nsitesmax-1 << endl;
  OutFile << " Ncutoff      = " << Ncutoff << endl;
  if (strcmp(ModelOption,"Anderson")==0){
    OutFile << " U            = " << U << endl;
    OutFile << " Gamma        = " << Gamma << endl;
    OutFile << " ed           = " << ed << endl;
  }
  if (strcmp(ModelOption,"Kondo")==0){
    OutFile << " JK            = " << Gamma << endl;
  }
  OutFile << " Lambda       = " << Lambda << endl;
  OutFile << " Dband        = " << Dband << endl;
  OutFile << " UpdateBefCut = " << UpdateBefCut << endl;
  OutFile << " calcdens     = " << calcdens << endl;
  OutFile.close();
  
  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  ////                   END READ PARAMETERS                   ////
  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////

  // Set SingleSite

  OneChQS_SetSingleSite(SingleSite);


  // Allocate MatArray (hope this works!)
  MatArray = new CNRGmatrix [NumNRGarrays];


  // Test Kondo model
  
  Params.clear();


  switch (ModelNo)
    {
    case 0 :
      // Set initial CNRG array (N=-1)
      Params.push_back(U);
      Params.push_back(ed);
      Params.push_back(Lambda);
      Params.push_back(HalfLambdaFactor);
      //OneChQS_SetAndersonHm1(Params,&Aeig, &Qm1fNQ, &MQQp1);
      OneChQS_SetAndersonHm1(Params,&Aeig, &Qm1fNQ, MatArray);

      Nsites0=0; // Adding N=0 site
      chi_N=chi_m1; // tunneling to site N=0 (hi_N connects sites N and N+1)
      SusDot=1.0/8.0;
      // Local Susceptibility,Entropy
      Suscep.dImpValue=SusDot;
      Entropy.dImpValue=2.0*log(2.0);
      break;
    case 1 :
      Params.push_back(Gamma); // JK
      OneChQS_SetKondoH0(Params,&Aeig,&Qm1fNQ,&SingleSite );
      Nsites0=1; // Adding N=1 site
      chi_N=chiN(0,Lambda);
      SusDot=1.0/4.0; // check this difference
      // Local Susceptibility,Entropy
      Suscep.dImpValue=SusDot;
      Entropy.dImpValue=log(2.0);
      break;
    case 2 :
      cout << " Not implemented yet..." << endl;
      exit(0);
      break;
    case 3 :
      OneChQS_SetH0Chain(&Aeig, &Qm1fNQ);
      Nsites0=1; // Adding N=1 site
      chi_N=chiN(0,Lambda); // tunneling to site N=1 
                            //(chi_N connects sites N and N+1)
      // Local Susceptibility,Entropy
      Suscep.CalcChain=true;
      Entropy.CalcChain=true;
      Suscep.dImpValue=0.0;
      Entropy.dImpValue=0.0;
      break;
    case 5 :
      OneChQS_SetSMM_Hm1(Params,&Aeig, &Qm1fNQ,MatArray);
      Nsites0=0; // N=0 site only
      chi_N=chiN(0,Lambda); // tunneling to site N=1 
                            //(chi_N connects sites N and N+1)
      KeepSz=1;
      calcdens=0;
      break;
    default :
      cout << " Model not implemented. Exiting... " << endl;
      exit(0);
      break;
    }
  //end switch models


  // Set AeigCut
  AeigCut.ClearAll();
  AeigCut=CutStates(&Aeig, Ncutoff);
  
 // Set operators to update


  // Will this work??? Yes!!

//   MatArray[0].NeedOld=false;
//   MatArray[0].CheckForMatEl=OneChQS_cd_check;
//   MatArray[0].CalcMatEl=OneChQS_fN_MatEl;



  //////////////
  Nsites=Nsites0-1;
  DN=HalfLambdaFactor*pow(Lambda,(-(Nsites-1)/2.0) );
  TM=DN/betabar;
  cout << "DN = " << DN << "TM = " << TM << endl;

  if (calcdens==1)
    {
      MatArray[0].NeedOld=true;
      MatArray[0].CheckForMatEl=OneChQS_nd_check;
      MatArray[0].CalcMatEl=OneChQS_nd_MatEl;

      MatArray[1].NeedOld=true;
      MatArray[1].CheckForMatEl=OneChQS_cd_check;
      MatArray[1].CalcMatEl=OneChQS_cd_MatEl;

      // Sep 08: calculating <Sz2>
      MatArray[2].NeedOld=true;
      MatArray[2].CheckForMatEl=OneChQS_nd_check;
      MatArray[2].CalcMatEl=OneChQS_nd_MatEl;



      Params.clear();
      Params.push_back(betabar);  
      ndot=CalcOpAvg(Params,&AeigCut,&MatArray[0],true,1);
      cout << " ed = " << ed << " Initial ndot = " << ndot << endl;
      Sz2=CalcOpAvg(Params,&AeigCut,&MatArray[2],true,1);
      cout << " Sz2 = " << Sz2 << " Initial Sz2 = " << Sz2 << endl;


      for (int ibl=0;ibl<MatArray[1].NumMatBlocks();ibl++)
	MatArray[1].PrintMatBlock(ibl);

      // Interface with old functions.
//       AeigCut.SaveQSParameters();
//       MatArray[1].SaveInOldFormat();

    }
  if ( (ModelNo==3)&&(calcdens==3) )
    {
      TM=DN/betabar;
      Params.clear();
      Params.push_back(betabar);
      // New Stuff
      Suscep.ReadNChainValue(Nsites,Nsites0); // Chain starts at N=0
      Suscep.AddValue(Params,&Aeig,1,true,TM);
      Suscep.SaveNValue(Nsites,Nsites0);
      Entropy.ReadNChainValue(Nsites,Nsites0);
      Entropy.AddValue(Params,&Aeig,1,true,TM);
      Entropy.SaveNValue(Nsites,Nsites0);
    }

  // Jul 09: testing CNRGchain

  //chain1.SetChainWilson(Nsitesmax);

  // Loop on Nsites

  // Entering calculation of H_Nsites0
  Nsites=Nsites0;
  while (Nsites<Nsitesmax)
    {
      cout << "Nsites = " << Nsites << endl;
      cout << "Old Nshell = " << Aeig.Nshell << endl;
      cout << "BEG Nshell = " << Aeig.Nshell+1 << endl;
      cout << "chi_(N-1) = " << chi_N << endl;
      DN=HalfLambdaFactor*pow(Lambda,(-(Nsites-1)/2.0) );
      TM=DN/betabar;
      cout << "DN = " << DN << "TM = " << TM << endl;

      // Build Abasis

      /////QS_BuildBasis(&AeigCut,&Abasis,&SingleSite,UpdateBefCut);
      QS_BuildBasis(&AeigCut,&Abasis,&SingleSite,UpdateBefCut,KeepSz);

      //Abasis.PrintAll();
      
      cout << "Basis Nstates = " << Abasis.Nstates() << endl;
      
      // Build and diagonalize H_N+1 (general): Build Aeigv
      // 1 - Get old matrix elements
      // 2 - Build and diagonalize H_N+1

      cout << "Diagonalizing HN... " << endl;

      OneChQS_DiagHN(Qm1fNQ,Abasis,SingleSite,Aeig,0.0,chi_N,Lambda);

      cout << "... done diagonalizing HN. " << endl;

      Aeig.PrintEn();
//       for (int ibl=0;ibl<Aeig.NumBlocks();ibl++)
// 	Aeig.PrintBlock(ibl);

      if (UpdateBefCut==1)
	{
	  cout << "Updating matrices before cutting... " << endl;
	  OneChQS_UpdateQm1fQ(&Qm1fNQ,&Aeig,&Abasis,&SingleSite);
	  cout << "... done updating matrices. " << endl;
	}


      // Calculate Susceptibility

      if ( (calcdens==2)||(calcdens==3) )
	{
	  TM=DN/betabar;
	  Params.clear();
	  Params.push_back(betabar);
	  // New Stuff
	  Suscep.ReadNChainValue(Nsites,Nsites0); // Chain starts at N=0
	  Suscep.AddValue(Params,&Aeig,1,true,TM);
	  Suscep.SaveNValue(Nsites,Nsites0);
	  Entropy.ReadNChainValue(Nsites,Nsites0);
	  Entropy.AddValue(Params,&Aeig,1,true,TM);
	  Entropy.SaveNValue(Nsites,Nsites0);
	}

      // Eliminate states, update matrices, calculate stuff 

      // Try this here!
      cout << "Cutting states..." << endl;
      
      AeigCut.ClearAll();
      AeigCut=CutStates(&Aeig, Ncutoff);

      cout << "... done cutting states." << endl;

      // Calculate new matrix elements using the CUT basis:
      // update Qm1fNQ, Qm1cdQ, etc.

      if (UpdateBefCut==0)
	{
	  cout << "Updating matrices after cutting... " << endl;
	  //OneChQS_UpdateMatrixAfterCutting(&Qm1fNQ,&MQQp1,
	  //          &AeigCut,&Abasis,&SingleSite);
	  OneChQS_UpdateMatrixAfterCutting(&Qm1fNQ,&MatArray[1],
					   &AeigCut,&Abasis,&SingleSite);
	  cout << "... done updating matrices. " << endl;

	  // New Update

	  // Calculate Spectral density, <ndot>
	  if (calcdens==1)
	    {	      
// 	      UpdateMatrices(&SingleSite,&AeigCut, 
// 			     &Abasis,MatArray, 2);
// Calc Sz2 too now: update 3 mats
	      UpdateMatrices(&SingleSite,&AeigCut, 
			     &Abasis,MatArray, 3);

// 	      for (int ibl=0;ibl<MatArray[1].NumMatBlocks();ibl++)
// 		MatArray[1].PrintMatBlock(ibl);

	      if (Nsites==0)
		for (int ist=0;ist<MatArray[1].Nstates();ist++)
		  {
		    for (int jst=0;jst<MatArray[1].Nstates();jst++) 
		      cout << MatArray[1].GetMatEl(ist,jst) << "  ";
		    cout << endl;
		  }

	      // Changed to NumMats to 2
	      TM=DN/betabar;
	      Params.clear();
	      Params.push_back(betabar);
	      ndot=CalcOpAvg(Params,&AeigCut,&MatArray[0],true,1);
	      	      
	      
	      cout << " ed = " << ed; 
	      cout.precision(10);
	      cout << scientific << " Temp = " << TM;
	      cout << resetiosflags (ios_base::floatfield);
	      cout << " ndot = " << ndot << endl;

	      Sz2=CalcOpAvg(Params,&AeigCut,&MatArray[2],true,1);
	      cout << " ed = " << ed; 
	      cout.precision(10);
	      cout << scientific << " Temp = " << TM;
	      cout << resetiosflags (ios_base::floatfield);
	      cout << " Sz2 = " << Sz2 << endl;


	      // Interface with old functions.
	      AeigCut.SaveQSParameters();
	      MatArray[1].SaveInOldFormat();
	    }
	  // end if calcdens==1
	  	  
	}
      // end if UpdateCefCut==0

      cout << "END Nshell = " << Aeig.Nshell << endl;

     // Update chi_N

//       daux[0]=(double)( 1.0-pow(Lambda,(-(Nsites+1))) );
//       daux[1]=(double)sqrt( 1.0-pow(Lambda,-(2*(Nsites+1)-1)) );
//       daux[2]=(double)sqrt( 1.0-pow(Lambda,-(2*(Nsites+1)+1)) );  
//       daux[3]=0.5*(1.0+(1.0/Lambda))*(double)sqrt(Lambda);

//       chi_N=daux[0]/(daux[1]*daux[2]);
// new (Jul 09)
      chi_N=chain1.GetChin(Nsites);

      cout << "chi_N = " << chi_N << endl;
//       cout << "Chi(N="<< Nsites <<") = " << chiN(Nsites,Lambda) << endl; 

      // Update sites
      Nsites++;
    }
  // end NRG loop

  // De-allocate MatArray 
  delete[] MatArray;


  cout << "=== Calculation Finished! ==== "<< endl;
  OutFile.open("NRG_end.txt");
  OutFile << "END NRG calculation" << endl;
  OutFile.close();



}
// END code



////////////////////////////////
///                          ///
///       Trash can          ///
///                          ///
////////////////////////////////


