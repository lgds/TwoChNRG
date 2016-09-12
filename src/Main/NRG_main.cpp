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
#include "NRGOpMatRules.hpp"
#include "NRG_main.hpp"
#include "TwoChQS.hpp"


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
  
  CNRGCodeHandler ThisCode;

#include"ModelOptMain.cpp"


  ThisCode.SaveData=false;

  // NRG objects
  
  CNRGarray Aeig;
  CNRGbasisarray AeigCut;
  CNRGbasisarray Abasis;
  CNRGbasisarray SingleSite;

  // STL vector

  CNRGmatrix HN;
  CNRGmatrix Qm1fNQ;

  //CNRGmatrix MQQp1;

  CNRGmatrix* MatArray;
  int NumNRGarrays=4;
  // Jul 09: Will this work???
  vector<CNRGmatrix> STLMatArray;
  CNRGmatrix auxNRGMat;
  // MatArray 0 is f_ch1
  // MatArray 1 is f_ch2
  // MatArray 2 is Sz
  // MatArray 3 is Sz2



  // STL vectors
  vector <double> Params;
  vector <double> ParamsHN;
  vector <double> ParamsBetabar;
  vector<int> CommonQNs; 
  vector<int> totSpos;

  // Thermodynamics

  CNRGthermo Suscep;
  CNRGthermo Entropy;

  CNRGthermo *ThermoArray;
  int NumThermoArrays=2;

  double TM=0.0;
  //double DN=0.0;
  //double betabar=0.727;
  ThisCode.Nsites=0;
  ThisCode.betabar=0.727;
  // Add more than one betabar in the code... Done!
  //ThisCode.betabar=0.6; 

  // instream
  ifstream InFile;
  // outstream
  ofstream OutFile;


  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  ///                                              ///
  ///                Main code                     ///
  ///                                              ///
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////


  // Read stuff (separate routine)

// Allocate MatArray (hope this works!)
// Yes, but, for some reason, needs to be BEFORE everything!


   MatArray = new CNRGmatrix [NumNRGarrays];
   ThermoArray = new CNRGthermo [NumThermoArrays];

   for (int imat=0;imat<NumNRGarrays;imat++){STLMatArray.push_back(auxNRGMat);}

   ThisCode.Nsites0=0;
   ThisCode.NumChannels=1;
   ThisCode.NumNRGmats=NumNRGarrays;
   ThisCode.NumThermoMats=NumThermoArrays;
   ThisCode.InitialSetUp();
   //InFile.open("input_nrg.dat"); // Needs this here 
                                   // new MatArray after!! Why????

   ThisCode.SetSingleSite(&SingleSite);

   // Copy pointers for saving/reading.
   ThisCode.pAbasis=&Abasis;
   //ThisCode.MatArray=MatArray;
   ThisCode.MatArray=&STLMatArray[0];
   ThisCode.pAcut=&AeigCut;

  // Set H0
  // Set initial Aeig and matrices (such as Qm1fNQ[])/Operators
  //  - Model dependent functions (hardest part)
  //  - Set quantum numbers, etc,etc,
  //  - HN. will 


  // Thermobasics

   strcpy(ThermoArray[0].ArqName,"SuscepImp_Main.dat");
   strcpy(ThermoArray[0].ChainArqName,"SuscepChain2Ch.dat");
   ThermoArray[0].Calc=CalcSuscep;
  
   strcpy(ThermoArray[1].ArqName,"EntropyImp_Main.dat");
   strcpy(ThermoArray[1].ChainArqName,"EntropyChain2Ch.dat");
   ThermoArray[1].Calc=CalcEntropy;

   for (int ithermo=0;ithermo<NumThermoArrays;ithermo++){
     if (ThisCode.calcdens==3)
       ThermoArray[ithermo].CalcChain=true;
     else
       ThermoArray[ithermo].CalcChain=false;
   }
   // HN basics
   HN.NeedOld=false;
   HN.UpperTriangular=true;
   HN.CheckForMatEl=Diag_check;


   // To do: send ALL THIS to codeHandler!!
   /// For all codes
   double Gamma=0.0282691;
   double Lambda=2.5;
   if (ThisCode.dInitParams.size()>1)
     Gamma=ThisCode.dInitParams[1];
   if (ThisCode.Lambda>0.0)
     Lambda=ThisCode.Lambda;
   //double HalfLambdaFactor=0.5*(1.0+(1.0/Lambda)); // also defined in CNRGcodehandler
   double HalfLambdaFactor=ThisCode.HalfLambdaFactor; // fix this double counting later

   double chi_m1=sqrt(2.0*Gamma/pi)/(sqrt(Lambda)*HalfLambdaFactor);

   // added: eN for chains away from phs (getting en for N=0)
   double eN=ThisCode.chain.GetEn(0);

   cout << " Fsq = " << ThisCode.chain.Fsq 
	<< " chi_m1 = " << chi_m1 
	<< " eN(0) = " << eN << endl;

   // Adding A_Lambda factor (except if using Campo-Oliveira)
   if (ThisCode.chain.DiscScheme!=1){
     chi_m1*=sqrt(0.5*log(Lambda)*(Lambda+1)/(Lambda-1));
   }

   cout  << " sqrt(ALambda)*chi_m1 = " << chi_m1  << endl;


   // Sep 08: 1ChQ chain calculations diverted to Anderson (ModelNo=0)
   // Aug 09: Why?? Is this correct?

   if ( (ThisCode.ModelNo==3)&&(ThisCode.SymNo==3) ){
     ThisCode.ZeroParams();
     ThisCode.ModelNo=0;
   }

   double DeltaSC=0.0;
   // Set 2015: added SC leads
   if ( (ThisCode.SymNo==7)||(ThisCode.SymNo==8) ){
     if (ThisCode.ModelNo==0){
       DeltaSC=ThisCode.dInitParams[3];
       cout << " DeltaSC = " << DeltaSC << endl;
       }
   }
   // end add DeltaSC

   /***************/
   // Model specific hamiltonians //
   /***************/	

   // to do: models for which Nsites0=1 need eN as input here

   ThisCode.ModelSwitch(CommonQNs,totSpos,&Aeig,&SingleSite,
			&HN,STLMatArray,ThermoArray,chi_m1);

   CutStates(Aeig, AeigCut, ThisCode.Ncutoff);

 
   // watch out for the "Nsites0-2"
   ThisCode.DN=HalfLambdaFactor*pow(Lambda,(-(ThisCode.Nsites0-2)/2.0) );
   // Aug 2010: Z-trick factor added into HalfLambdaFactor
   TM=ThisCode.DN/ThisCode.betabar;
   cout << "DN = " << ThisCode.DN << "  TM = " << TM << endl;

   // 
   for (int ithermo=0;ithermo<NumThermoArrays;ithermo++){
     ThisCode.ThermoSTLArray.push_back(ThermoArray[ithermo]);
   }

   if (ThisCode.Nsites0==1) // Calculates Thermo from N=0 on (not N=-1)
     ThisCode.CalcThermo(ThermoArray,&Aeig);

   ThisCode.CalcStuff(&Aeig,&AeigCut);


   // Save Params?
   if (ThisCode.SaveData){ThisCode.SaveGenPars();}

   // Ok, with that we can go into the main loop

   // Parameters for HN (always the same?)
   ParamsHN.push_back(Lambda);
   for (int ich=0;ich<ThisCode.NumChannels;ich++)
     ParamsHN.push_back(chi_m1);
   // This is not good. Depends on setting chi_m1 at ModelSwitch!!
   // Also: if Nsites0 is NOT -1, we need something like this
   //        ParamsHN[ich+1]=ThisCode.chain.GetChin(ThisCode.Nsites);
   // What if the chains are not equivalent?

   // added: chains away from phs
   ParamsHN.push_back(eN);
   // added: superconducting leads
   ParamsHN.push_back(DeltaSC);


   // Entering calculation of H_Nsites0
   ThisCode.Nsites=ThisCode.Nsites0;

   while (ThisCode.Nsites<ThisCode.Nsitesmax){
     ThisCode.DN=HalfLambdaFactor*pow(Lambda,(-(ThisCode.Nsites-1)/2.0) );
     // Aug 2010: Z-trick factor added into HalfLambdaFactor
     TM=ThisCode.DN/ThisCode.betabar;
     cout << "DN = " << ThisCode.DN << " TM = " << TM << endl;
     // Debugging...
     //bool disp=(ThisCode.Nsites==3?true:false);
     //bool disp=true;
     bool disp=false;

     // Build Basis

     //SingleSite.PrintQNumbers();
     //Aeig.PrintQNumbers();


     // If UpdateBefCut=1, StCameFrom is the "old" (Before Cutting) one! 
     // Matrices will be in sync with Aeig, not AeigCut 

     BuildBasis(CommonQNs, totSpos,&AeigCut,&Abasis, 
		&SingleSite,ThisCode.UpdateBefCut);


     // Diagonalize HN

     cout << "Diagonalizing H(N="<<ThisCode.Nsites<<")... " << endl;
     cout << " chi_N = " << ParamsHN[1] << endl;

     // Need to adapt ALL models to use STLMatArray... 
     // in the meantime, this will work.
     if ( (ThisCode.ModelNo==4) ){
       HN.DiagHN(ParamsHN,&Abasis,&SingleSite,MatArray,&Aeig);
     }else{
       HN.DiagHN(ParamsHN,&Abasis,&SingleSite,&STLMatArray[0],&Aeig,disp);
     }

     cout << "... done diagonalizing HN. " << endl;
    

     Aeig.SetKept(ThisCode.Ncutoff);
     Aeig.PrintEn();

     //       for (int ibl=0;ibl<Aeig.NumBlocks();ibl++)
     // 	Aeig.PrintBlock(ibl);

     // Calculate Susceptibility/Entropy (does not need update)

     ThisCode.CalcThermo(ThermoArray,&Aeig);


     // NEW (updates after cutting)

     if (ThisCode.UpdateBefCut==0){
       // Update AFTER CUTTING (faster)
       CutStates(Aeig, AeigCut, ThisCode.Ncutoff);
     }else{
       // Update BEFORE CUTTING
       //Aeig.SetKept(ThisCode.Ncutoff);
       CutStates(Aeig, AeigCut, 2*Aeig.Nstates());
       // AeigCut is essentially a copy of Aeig
     }
     // end if UpdateBefCut


     // Update matrices (in either case)

     if (ThisCode.Nsites<ThisCode.Nsitesmax){
       // Need to adapt ALL models to use STLMatArray... 
       // in the meantime, this will work.
       if ( (ThisCode.ModelNo==4) ){
	 UpdateMatrices(&SingleSite,&AeigCut, 
			&Abasis,MatArray, ThisCode.NumNRGmats);

       }else{
	 UpdateMatrices(&SingleSite,&AeigCut, 
			&Abasis,&STLMatArray[0],STLMatArray.size(),disp);
       }
     }
     // end update matrices

     // Save stuff to files


     if (ThisCode.SaveData){ThisCode.SaveArrays();}

 
     // Calculate Other things that need updated matrices
     //
     // The idea is to eliminate this and put into CalcStuff
     // (need to test with OneChQ!!
     if (ThisCode.SymNo==3){
       ParamsBetabar.clear();
       ParamsBetabar.push_back(ThisCode.betabar);
       double dSz=CalcOpAvg(ParamsBetabar,&AeigCut,&STLMatArray[2],false,0);
       double dSz2=CalcOpAvg(ParamsBetabar,&AeigCut,&STLMatArray[3],false,0);
       //double dSz=CalcOpAvg(ParamsBetabar,&AeigCut,MatArray,false,0);
       //double dSz2=CalcOpAvg(ParamsBetabar,&AeigCut,MatArray,false,0);
       
       cout << " T = " << TM << " Sz = " << dSz << " Sz2 = " << dSz2 
	    << "  T_M chi = " << dSz2-dSz*dSz << endl;
     }

     ThisCode.CalcStuff(&Aeig,&AeigCut);


     if (ThisCode.UpdateBefCut==1){
       // Now we cut for real
       CutStates(Aeig, AeigCut, ThisCode.Ncutoff);
     }

     cout << " N = " << ThisCode.Nsites << endl;

     // Update ParamsHN
     //ParamsHN[1]=chiN(ThisCode.Nsites,Lambda); // Old
     for (int ich=0;ich<ThisCode.NumChannels;ich++)
       //ParamsHN[ich+1]=chiN(ThisCode.Nsites,Lambda);
       // Try this...
       ParamsHN[ich+1]=ThisCode.chain.GetChin(ThisCode.Nsites);

     
     // Update sites
     ThisCode.Nsites++;

     // Adding eN: is this before or after Nsites is updated? After!
     if (ParamsHN.size()>ThisCode.NumChannels+2){
       eN=ThisCode.chain.GetEn(ThisCode.Nsites);
       ParamsHN[ThisCode.NumChannels+1]=eN;
       cout << " eN_" << ThisCode.Nsites << " = " << eN << endl;
       // Adding DeltaSC RENORMALIZED
       double ScaleFactorNm1=ThisCode.chain.ScaleFactor(ThisCode.Nsites-1);
       ParamsHN[ThisCode.NumChannels+2]=ScaleFactorNm1*DeltaSC;
       cout << " ScaleFactorN_" << ThisCode.Nsites-1 << " = " << ScaleFactorNm1
	    << " DeltaSC_scaled = " << " DeltaSC_scaled = " << ScaleFactorNm1*DeltaSC
	    << endl;

     }
     // apparently this is ALWAYS true!

     
   }
   // end Nsites loop

   // De-allocate MatArray 
   delete[] ThermoArray;
   delete[] MatArray;

   ThisCode.WrapUp();

}
// end main


//
// Trash
//
//////////////////////////


