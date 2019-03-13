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

//   char ModelOption[]="Anderson";
//   int ModelNo=5;

  // Read Model, etc. (not now)

#include"ModelOptMain.cpp"


  strcpy(ThisCode.Symmetry,ModelSymmetry);
  strcpy(ThisCode.ModelOption,ModelOption);
  ThisCode.ModelNo=ModelNo;
  ThisCode.SymNo=SymNo;
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

  // Integers
  //int Nsites;

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
   //InFile.open("input_nrg.dat"); // Needs this here is I 
                                      // new MatArray after!! Why????

   ThisCode.SetSingleSite(&SingleSite);

   // Copy pointers for saving/reading.
   ThisCode.pAbasis=&Abasis;
   //ThisCode.MatArray=MatArray;
   ThisCode.MatArray=&STLMatArray[0];
   ThisCode.pAcut=&AeigCut;



  // Set H0
  // HN.SetH0(Params,&SingleSite,&Aeig,&Abasis,MatArray);
  //(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0,CNRGmatrix* NRGMats)

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


   /// For all codes
   double Gamma=0.0282691;
   double Lambda=2.5;
   if (ThisCode.dInitParams.size()>1)
     Gamma=ThisCode.dInitParams[1];
   if (ThisCode.Lambda>0.0)
     Lambda=ThisCode.Lambda;
   double HalfLambdaFactor=0.5*(1.0+(1.0/Lambda));
   double chi_m1=sqrt(2.0*Gamma/pi)/(sqrt(Lambda)*HalfLambdaFactor);

   // Jul 09: NRG chain
   CNRGchain chain1(ThisCode.Lambda,ThisCode.Nsitesmax);


   // Sep 08: 1ChQ chain calculations diverted to Anderson (ModelNo=0)
   // Aug 09: Why?? Is this correct?

   if ( (ThisCode.ModelNo==3)&&(ThisCode.SymNo==3) ){
     ThisCode.ZeroParams();
     ThisCode.ModelNo=0;
   }

 
 
///////////////////////////////////
////   Model specific hamiltonians ///
///////////////////////////////////

   switch(ThisCode.ModelNo){
   case 0: // Anderson
     switch(ThisCode.SymNo){
     case 2: // Anderson model: OneChQSz
       
       // Initialize matrices (OneChQ routines work here!)
       // Jul 09: using STMatArray

       STLMatArray[0].NeedOld=false;	  
       STLMatArray[0].CheckForMatEl=OneChQ_cd_check;
       STLMatArray[0].CalcMatEl=OneChQ_fNup_MatEl;
       STLMatArray[0].SaveMatYN=false;

       STLMatArray[1].NeedOld=false;
       STLMatArray[1].CheckForMatEl=OneChQ_cd_check;
       STLMatArray[1].CalcMatEl=OneChQ_fNdn_MatEl;
       STLMatArray[1].SaveMatYN=false;


       // Param for H0
       Params.push_back(ThisCode.dInitParams[0]); // U
       Params.push_back(ThisCode.dInitParams[2]); // ed
       Params.push_back(Lambda);
       Params.push_back(HalfLambdaFactor);
       Params.push_back(ThisCode.dInitParams[3]);        // Mag Field 

       //OneChQSz_SetAnderson_Hm1(Params,&Aeig,MatArray);
       OneChQSz_SetAnderson_Hm1(Params,&Aeig,STLMatArray);

       // cd_up and cd_dn (OneChQSz_SetAnderson equals 
       // MatArray[2,3] to MatArray[0,1]
       STLMatArray[2].NeedOld=true;
       STLMatArray[2].CheckForMatEl=OneChQSz_cdup_check;
       STLMatArray[2].CalcMatEl=OneChQSz_cdup_MatEl;
       STLMatArray[2].SaveMatYN=true;

       STLMatArray[3].NeedOld=true;
       STLMatArray[3].CheckForMatEl=OneChQSz_cddn_check;
       STLMatArray[3].CalcMatEl=OneChQSz_cddn_MatEl;
       STLMatArray[3].SaveMatYN=true;

       cout << " N = -1" << endl;
       // 	cout << " fN_up : " << endl;
       // 	MatArray[0].PrintAllBlocks();
       // 	cout << " fN_dn : " << endl;
       // 	MatArray[1].PrintAllBlocks();
       cout << " cdN_up : " << endl;
       STLMatArray[2].PrintAllBlocks();
       cout << " cdN_dn : " << endl;
       STLMatArray[3].PrintAllBlocks();


       ////////////////////////

       // BuildBasis params
       CommonQNs.push_back(2); // No of common QNs
       CommonQNs.push_back(0); // pos of QN 1 in old
       CommonQNs.push_back(1); // pos of QN 2 in old

       CommonQNs.push_back(0); // pos of QN 1 in SingleSite
       CommonQNs.push_back(1); // pos of QN 1 in SingleSite

       // No total S variables. Leave totSpos empty
       HN.CalcHNMatEl=OneChQSz_HN_MatEl;

       ThisCode.NumChannels=1;
       ThisCode.Nsites0=0;
       //ThisCode.NumNRGmats=4;//2
       ThisCode.NumNRGmats=STLMatArray.size();
       ThisCode.MatArray=&STLMatArray[0]; // Need to do this after
                                          // changes in STLMatArray!

       ThisCode.SaveData=true;

       break;
	
     case 3: // Anderson model: OneChQ

       // OLD!!

//        MatArray[0].NeedOld=false;	  
//        MatArray[0].CheckForMatEl=OneChQ_cd_check;
//        MatArray[0].CalcMatEl=OneChQ_fNup_MatEl;
	  
//        MatArray[1].NeedOld=false;
//        MatArray[1].CheckForMatEl=OneChQ_cd_check;
//        MatArray[1].CalcMatEl=OneChQ_fNdn_MatEl;

//        // Sz and Sz2
	  
//        MatArray[2].NeedOld=true;
//        MatArray[2].UpperTriangular=false;
//        MatArray[2].CheckForMatEl=Diag_check;
//        MatArray[2].CalcMatEl=ImpOnly_MatEl;

//        MatArray[3].NeedOld=true;
//        MatArray[3].UpperTriangular=false;
//        MatArray[3].CheckForMatEl=Diag_check;
//        MatArray[3].CalcMatEl=ImpOnly_MatEl;


       // Initialize matrices (uncomment later)

       STLMatArray[0].NeedOld=false;	  
       STLMatArray[0].CheckForMatEl=OneChQ_cd_check;
       STLMatArray[0].CalcMatEl=OneChQ_fNup_MatEl;
	  
       STLMatArray[1].NeedOld=false;
       STLMatArray[1].CheckForMatEl=OneChQ_cd_check;
       STLMatArray[1].CalcMatEl=OneChQ_fNdn_MatEl;

       // Sz and Sz2
       STLMatArray[2].NeedOld=true;
       STLMatArray[2].UpperTriangular=false;
       STLMatArray[2].CalcAvg=true;  // CalcAvg
       STLMatArray[2].CheckForMatEl=Diag_check;
       STLMatArray[2].CalcMatEl=ImpOnly_MatEl;
       strcpy(STLMatArray[2].MatName,"Sz");

       STLMatArray[3].NeedOld=true;
       STLMatArray[3].UpperTriangular=false;
       STLMatArray[3].CalcAvg=true;  // CalcAvg
       STLMatArray[3].CheckForMatEl=Diag_check;
       STLMatArray[3].CalcMatEl=ImpOnly_MatEl;
       strcpy(STLMatArray[3].MatName,"Sz2");

       if (ThisCode.calcdens==1){
	 // All operators have the same Checks/MatEls in |Q> basis
	 // what changes is the initial set-up
	 auxNRGMat.NeedOld=true;
	 auxNRGMat.UpperTriangular=false;
	 auxNRGMat.CheckForMatEl=OneChQ_cd_check;
	 auxNRGMat.CalcMatEl=OneChQ_cd_MatEl;
	 auxNRGMat.SaveMatYN=true;

	 STLMatArray.push_back(auxNRGMat); // 4 - cd1_up
	 STLMatArray.push_back(auxNRGMat); // 5 - cd1_dn
	 NumNRGarrays=STLMatArray.size();
	 ThisCode.NumNRGmats=NumNRGarrays;
	 ThisCode.MatArray=&STLMatArray[0]; //
	 // changes in STLMatArray!
	 ThisCode.SaveData=true;
	 strcpy(ThisCode.SaveArraysFileName,"1chQAnd");
       }
       else{
	 // Set Thermodynamics
	 // Entropy calculation only
	 strcpy(ThermoArray[0].ArqName,"EntropyImp1Ch_Anderson_Q.dat");
	 strcpy(ThermoArray[0].ChainArqName,"EntropyChain1Ch_Q.dat");
	 ThermoArray[0].Calc=CalcEntropy;
	 ThermoArray[0].dImpValue=2.0*log(2.0);
	 //ThermoArray[0].CalcChain=false;
	 ThisCode.NumThermoMats=1;
       }
       // calculate spec function or thermo?


       // Param for H0
       Params.push_back(ThisCode.dInitParams[0]); // U
       Params.push_back(ThisCode.dInitParams[2]); // ed
       Params.push_back(Lambda);
       Params.push_back(HalfLambdaFactor);
       OneChQ_SetAnderson_Hm1(Params,&Aeig,STLMatArray);
       //OneChQ_SetAnderson_Hm1_old(Params,&Aeig,MatArray);


       ////////////////////////

       // BuildBasis params
       CommonQNs.push_back(1); // No of common QNs
       CommonQNs.push_back(0); // pos of QN 1 in old
       CommonQNs.push_back(0); // pos of QN 1 in SingleSite
       // No total S variables. Leave totSpos empty
       HN.CalcHNMatEl=OneChQ_HN_MatEl;

       ThisCode.NumChannels=1;
       ThisCode.Nsites0=0;


       break;
     default:
       cout << " Symmetry " << ThisCode.Symmetry 
	    << " not implemented for " << ThisCode.ModelOption << endl;
       exit(0);
     }
     break;
   case 3: // Chain only
     ThisCode.ZeroParams();
     ThisCode.Nsites0=1;
     ThisCode.NumThermoMats=2;
     switch(ThisCode.SymNo){
     case 0: // Chain model: OneChQS
       // fN (reduced)
       STLMatArray[0].NeedOld=false;
       STLMatArray[0].UpperTriangular=false;
       STLMatArray[0].CheckForMatEl=OneChQS_cd_check;
       STLMatArray[0].CalcMatEl=OneChQS_fN_MatEl;
       STLMatArray.pop_back(); // Only ONE array
       STLMatArray.pop_back(); // Only ONE array
       STLMatArray.pop_back(); // Only ONE array
       // Set ChainH0: SingleSite is already set, just get Aeig and STLMatArray.
       OneChQS_SetChainH0(&Aeig,STLMatArray);
       ThisCode.NumNRGmats=STLMatArray.size();
       ThisCode.MatArray=&STLMatArray[0]; // Need to do this after
       ThisCode.SaveData=false;
       // BuildBasis params
       CommonQNs.push_back(2); // Q and S and commont QNs
       CommonQNs.push_back(0); // position of Q in old basis
       CommonQNs.push_back(1); // position of S in old basis
       CommonQNs.push_back(0); // position of Q in SingleSite
       CommonQNs.push_back(1); // position of S in SingleSite
       totSpos.push_back(1);   // SU(2) symmetry in position 1
       // HN params
       HN.CalcHNMatEl=OneChQS_HN_MatEl;
       ThisCode.NumChannels=1;
       chi_m1=chiN(0,Lambda); // is actually chi0
       strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QS.dat");
       strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QS.dat");
       break;
//      case 1: // Chain model: TwoChQS
//        TwoChQS_SetChainH0(&Aeig,STLMatArray);
//        break;
//      case 2: // Chain model: OneChQSz
//        OneChQSz_SetChainH0(&Aeig,STLMatArray);
//        break;
//      case 3: // Chain model: OneChQ
//        OneChQ_SetChainH0(&Aeig,STLMatArray);
//        break;
     default:
       cout << " Symmetry " << ThisCode.Symmetry 
	    << " not implemented for " << ThisCode.ModelOption << endl;
       exit(0);
     }
     //cout << " Chain calculation still under implementation" << endl;
     //exit(0);
     break;
   case 5: // SMM model, symmetry is OneChQ
     // Params for H0
     Params.push_back(Lambda);  // Lambda
     Params.push_back(HalfLambdaFactor); // Lambda factor


     if (ThisCode.dInitParams.size()>10){
       Params.push_back(ThisCode.dInitParams[0]); // U1
       Params.push_back(ThisCode.dInitParams[2]); // ed1

       Params.push_back(ThisCode.dInitParams[3]); // U2
       Params.push_back(ThisCode.dInitParams[5]); // ed2
 
       Params.push_back(ThisCode.dInitParams[6]); // J12

       Params.push_back(ThisCode.dInitParams[7]); // BmagPar      
       Params.push_back(ThisCode.dInitParams[8]); // BmagPerp

       Params.push_back(ThisCode.dInitParams[9]); // Dz      
       Params.push_back(ThisCode.dInitParams[10]); // B2 anisot

       double chi2_m1=sqrt(2.0*ThisCode.dInitParams[4]/pi)/(sqrt(Lambda)*HalfLambdaFactor);
       Params.push_back(chi_m1); // gamma1
       Params.push_back(chi2_m1); // gamma2
     }
     else{
       double U=0.5;
       double ed=-0.5*U;
       Params.push_back(U); // U1
       Params.push_back(ed); // ed1

       Params.push_back(0.0); // U2
       Params.push_back(0.0); // ed2
 
       Params.push_back(0.0); // J12

       Params.push_back(0.0); // BmagPar      
       Params.push_back(0.0); // BmagPerp

       Params.push_back(0.0); // Dz      
       Params.push_back(0.0); // B2 anisot

       Params.push_back(chi_m1); // gamma1
       Params.push_back(0.0); // gamma2
     }


     // Set rules for MatArray
     // Jul 09: using STMatArray

     // fN_up
     STLMatArray[0].NeedOld=false;
     STLMatArray[0].CheckForMatEl=OneChQ_cd_check;
     STLMatArray[0].CalcMatEl=OneChQ_fNup_MatEl;

     // fN_dn
     STLMatArray[1].NeedOld=false;
     STLMatArray[1].CheckForMatEl=OneChQ_cd_check;
     STLMatArray[1].CalcMatEl=OneChQ_fNdn_MatEl;


     // Sz and Sz2

     STLMatArray[2].NeedOld=true;
     STLMatArray[2].UpperTriangular=false;
     STLMatArray[2].CalcAvg=true;  // CalcAvg
     STLMatArray[2].CheckForMatEl=Diag_check;
     STLMatArray[2].CalcMatEl=ImpOnly_MatEl;
     strcpy(STLMatArray[2].MatName,"Sz");


     STLMatArray[3].NeedOld=true;
     STLMatArray[3].UpperTriangular=false;
     STLMatArray[2].CalcAvg=true;
     STLMatArray[3].CheckForMatEl=Diag_check;
     STLMatArray[3].CalcMatEl=ImpOnly_MatEl;
     strcpy(STLMatArray[3].MatName,"Sz2"); 

     // cd1_up, cd1_dn, cd2_up, cd2_dn  
     // MatArray is an STL vector!!!
     if (ThisCode.calcdens==1){
       // All operators have the same Checks/MatEls in |Q> basis
       // what changes is the initial set-up
       auxNRGMat.NeedOld=true;
       auxNRGMat.CheckForMatEl=OneChQ_cd_check;
       auxNRGMat.CalcMatEl=OneChQ_cd_MatEl;
       auxNRGMat.SaveMatYN=true;

       STLMatArray.push_back(auxNRGMat); // 4 - cd1_up
       STLMatArray.push_back(auxNRGMat); // 5 - cd1_dn
       STLMatArray.push_back(auxNRGMat); // 6 - cd2_up
       STLMatArray.push_back(auxNRGMat); // 7 - cd2_dn
       NumNRGarrays=STLMatArray.size();
       ThisCode.NumNRGmats=NumNRGarrays;
       ThisCode.MatArray=&STLMatArray[0]; // Need to do this after
                                          // changes in STLMatArray!
       ThisCode.SaveData=true;
       strcpy(ThisCode.SaveArraysFileName,"SMM");
      
     }
     else{
       // Set Thermodynamics
	
       // Entropy calculation only
      
       strcpy(ThermoArray[0].ArqName,"EntropyImp1Ch_Anderson_Q.dat");
       strcpy(ThermoArray[0].ChainArqName,"EntropyChain1Ch_Q.dat");
       ThermoArray[0].Calc=CalcEntropy;
       ThermoArray[0].dImpValue=4.0*log(2.0);
       ThermoArray[0].CalcChain=false;
       ThisCode.NumThermoMats=1;
     }
     // calculate spec function or thermo?


     // Set H0: output pAeig, MatArray



     OneChQ_SetSMM_H0(Params,&Aeig,&SingleSite,STLMatArray);


     // BuildBasis params
     CommonQNs.push_back(1); // No of common QNs
     CommonQNs.push_back(0); // pos of QN 1 in old
     CommonQNs.push_back(0); // pos of QN 1 in SingleSite
     // No total S variables. Leave totSpos empty
     HN.CalcHNMatEl=OneChQ_HN_MatEl;
     ThisCode.Nsites0=1;
     ThisCode.NumChannels=1;

     chi_m1=chiN(0,Lambda);

     break;

   case 4: // CM Phonons - 2channel
     switch(ThisCode.SymNo){
     case 4:{ // TwoChQSP
       Params.push_back(Lambda);  // Lambda
       Params.push_back(HalfLambdaFactor); // Lambda factor

       Params.push_back(ThisCode.dInitParams[0]); // U1
       Params.push_back(ThisCode.dInitParams[2]); // ed1

       Params.push_back(chi_m1); // sqrt(2gamma1/Pi)/sqrt(L)*HalfLambdaFactor

       Params.push_back(ThisCode.dInitParams[3]); // w0
       Params.push_back(ThisCode.dInitParams[4]); // lambda_ph
       Params.push_back(ThisCode.dInitParams[5]); // alpha
       Params.push_back(ThisCode.dInitParams[6]); // Nph      


       // Update rules for f_ch1 and fch2

       // Set rules for MatArray

       // fN_ch1 (reduced)
       MatArray[0].NeedOld=false;
       MatArray[0].CheckForMatEl=TwoChQS_cd_check;
       MatArray[0].CalcMatEl=TwoChQS_cd_ich1_Phonon_MatEl;

       // fN_ch2 (reduced)
       MatArray[1].NeedOld=false;
       MatArray[1].CheckForMatEl=TwoChQSP_fA_check; // mixes parities
       MatArray[1].CalcMatEl=TwoChQS_cd_ich2_Phonon_MatEl;

       //CNRGarray MyTest1;
       

       // Set H0, get pAeig and Update MatArray
       TwoChQSP_SetH0CMphonon(Params, &Aeig,&SingleSite,MatArray);


       // No phonons needed from now on
       MatArray[0].CalcMatEl=TwoChQS_fNch1_MatEl;
       MatArray[1].CalcMatEl=TwoChQS_fNch2_MatEl;


       // BuildBasis params
       CommonQNs.push_back(3); // No of common QNs
       CommonQNs.push_back(0); // pos of QN 1 in old
       CommonQNs.push_back(1); // pos of QN 2 in old
       CommonQNs.push_back(2); // pos of QN 3 in old

       CommonQNs.push_back(0); // pos of QN 1 in SingleSite
       CommonQNs.push_back(1); // pos of QN 2 in SingleSite
       CommonQNs.push_back(3); // pos of QN 3 in SingleSite

       CommonQNs.push_back(2); // Pos of parity in old
       totSpos.push_back(1); // Pos of S in old

       // Hamiltonian/Code params
       HN.CalcHNMatEl=TwoChQS_HN_MatEl;
       ThisCode.Nsites0=1;
       chi_m1=chiN(0,Lambda);
       ThisCode.NumNRGmats=2;
       ThisCode.NumChannels=2;
       ThisCode.NumThermoMats=2;

       ThermoArray[0].dImpValue=0.0;
       ThermoArray[1].dImpValue=0.0; // Check this!!

       // Thermobasics

       strcpy(ThermoArray[0].ArqName,"SuscepImp2Ch_CMphonon_QSP.dat");
       strcpy(ThermoArray[0].ChainArqName,"SuscepChain2Ch.dat");
       ThermoArray[0].Calc=CalcSuscep;
  
       strcpy(ThermoArray[1].ArqName,"EntropyImp2Ch_CMphonon_QSP.dat");
       strcpy(ThermoArray[1].ChainArqName,"EntropyChain2Ch.dat");
       ThermoArray[1].Calc=CalcEntropy;

       cout << " Implementing this model NOW " << endl;

       cout << " fN_ch1 : " << endl;
       MatArray[0].PrintAllBlocks();
       cout << " fN_ch2 : " << endl;
       MatArray[1].PrintAllBlocks();
       break;
     }// Need curly braces in "case" is declaring anything in it
     default:
       cout << " Symmetry " << ThisCode.Symmetry 
	    << " not implemented for " << ThisCode.ModelOption << endl;
       exit(0);
     }
     break;
   case 6: // Double Quantum dot 
     switch(ThisCode.SymNo){
     case 0:{ // OneChQS
       cout << " Implementing DQD case... " << endl;

       double chi2_m1=sqrt(2.0*ThisCode.dInitParams[4]/pi)/(sqrt(Lambda)*HalfLambdaFactor);

       Params.push_back(Lambda);
       Params.push_back(HalfLambdaFactor);
       Params.push_back(ThisCode.dInitParams[0]); // U1
       Params.push_back(ThisCode.dInitParams[2]); // ed1
       Params.push_back(chi_m1); // gamma1
       Params.push_back(ThisCode.dInitParams[3]); // U2
       Params.push_back(ThisCode.dInitParams[5]); // ed2
       Params.push_back(chi2_m1); // gamma2

       // Question: Can these things be implemented within the subroutine??

       // fN (reduced)
       STLMatArray[0].NeedOld=false;
       STLMatArray[0].UpperTriangular=false;
       STLMatArray[0].CheckForMatEl=OneChQS_cd_check;
       STLMatArray[0].CalcMatEl=OneChQS_fN_MatEl;


       switch (ThisCode.calcdens){
       case 0: // Calculates Op averages only
	 auxNRGMat.NeedOld=true;
	 auxNRGMat.UpperTriangular=false; // update procedure only works for 'false'
	 auxNRGMat.CalcAvg=true;  // CalcAvg
	 auxNRGMat.CheckForMatEl=Diag_check;
	 auxNRGMat.CalcMatEl=ImpOnly_MatEl;
	 STLMatArray[1]=auxNRGMat; // 1 - ndot (there already)
	 STLMatArray[2]=auxNRGMat; // 2 - ndot^2 (there already)
	 STLMatArray[3]=auxNRGMat; // 3 - Sdot^2 (there already)
	 STLMatArray.push_back(auxNRGMat); // 4 - S1 dot S2
	 strcpy(STLMatArray[1].MatName,"Ndqd");
	 strcpy(STLMatArray[2].MatName,"NdqdSq");
	 strcpy(STLMatArray[3].MatName,"SdqdSq");
	 strcpy(STLMatArray[4].MatName,"S1dotS2");
	 NumNRGarrays=STLMatArray.size();
	 ThisCode.NumNRGmats=NumNRGarrays;
	 ThisCode.MatArray=&STLMatArray[0]; // Need to do this after
	 // changes in STLMatArray!
	 ThisCode.SaveData=false;
	 break;
       case 1: // Calc Spectral functions
	 // All operators have the same Checks/MatEls in |Q> basis
	 // what changes is the initial set-up
	 auxNRGMat.NeedOld=true;
	 auxNRGMat.CheckForMatEl=OneChQS_cd_check;
	 auxNRGMat.CalcMatEl=OneChQS_cd_MatEl;
	 auxNRGMat.SaveMatYN=true;
	 STLMatArray[1]=auxNRGMat; // 1 - cd1
	 STLMatArray[2]=auxNRGMat; // 2 - cd2
	 strcpy(STLMatArray[1].MatName,"cdot1");
	 strcpy(STLMatArray[2].MatName,"cdot2");
	 STLMatArray.pop_back(); // 3 elements only
	 NumNRGarrays=STLMatArray.size();
	 ThisCode.NumNRGmats=NumNRGarrays;
	 ThisCode.MatArray=&STLMatArray[0]; // Need to do this after
	 // changes in STLMatArray!
	 ThisCode.SaveData=true;
	 strcpy(ThisCode.SaveArraysFileName,"DQD");
	 break;
       case 2: // Thermodynamics
	 STLMatArray.pop_back(); // 1 Matrix only
 	 STLMatArray.pop_back(); // 1 Matrix only
 	 STLMatArray.pop_back(); // 1 Matrix only
      	 strcpy(ThermoArray[0].ArqName,"SuscepImp1ChQS_DQD.dat");
	 strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QS.dat");
      	 strcpy(ThermoArray[1].ArqName,"EntropyImp1ChQS_DQD.dat");
	 strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QS.dat");
	 ThermoArray[0].Calc=CalcSuscep;
	 ThermoArray[0].dImpValue=1.0/4.0; // check
	 ThermoArray[1].Calc=CalcEntropy;
	 ThermoArray[1].dImpValue=4.0*log(2.0);
	 ThisCode.NumThermoMats=2;
	 break;
       default:
	 cout << " Calculating levels only " << endl;

       }
       // end switch calcdens


       OneChQS_SetHm1DQD(Params,&Aeig,&SingleSite,STLMatArray);

       // BuildBasis params
       CommonQNs.push_back(2); // Q and S and commont QNs
       CommonQNs.push_back(0); // position of Q in old basis
       CommonQNs.push_back(1); // position of S in old basis
       CommonQNs.push_back(0); // position of Q in SingleSite
       CommonQNs.push_back(1); // position of S in SingleSite
       totSpos.push_back(1);   // SU(2) symmetry in position 1
       // HN params
       HN.CalcHNMatEl=OneChQS_HN_MatEl;
       ThisCode.Nsites0=1;
       ThisCode.NumChannels=1;
       chi_m1=chiN(0,Lambda);
       break;
     }
     default:
       cout << " Symmetry " << ThisCode.Symmetry 
	    << " not implemented for " << ThisCode.ModelOption << endl;
       exit(0);
     }
     break;
   default:
     cout << " Model not implemented. Exiting... " << endl;
     exit(0);
   }
   // end switch ModelNo 


   // Set Initial AeigCut
   //CNRGbasisarray AeigCut=CutStates(&Aeig, ThisCode.Ncutoff);
   // ThisCode.pAcut=&AeigCut;

   CutStates(Aeig, AeigCut, ThisCode.Ncutoff);
 
  //AeigCut=CutStates(&Aeig, ThisCode.Ncutoff);
   cout << "Got here... I" << endl;

   // watch out for the "Nsites0-2"
   ThisCode.DN=HalfLambdaFactor*pow(Lambda,(-(ThisCode.Nsites0-2)/2.0) );
   TM=ThisCode.DN/ThisCode.betabar;
   cout << "DN = " << ThisCode.DN << "  TM = " << TM << endl;

   //
   // The idea is to eliminate this and put into CalcStuff
   //
//    if (ThisCode.SymNo==3){
//      ParamsBetabar.clear();
//      ParamsBetabar.push_back(ThisCode.betabar);
//      //double dSz=CalcOpAvg(ParamsBetabar,&AeigCut,&MatArray[2],false,0);
//      //double dSz2=CalcOpAvg(ParamsBetabar,&AeigCut,&MatArray[3],false,0);
//      double dSz=CalcOpAvg(ParamsBetabar,&AeigCut,&STLMatArray[2],false,0);
//      double dSz2=CalcOpAvg(ParamsBetabar,&AeigCut,&STLMatArray[3],false,0);

//      cout << " T = " << TM << " Sz = " << dSz << " Sz2 = " << dSz2 
// 	  << "  T_M chi = " << dSz2-dSz*dSz << endl;
//    }

   // 
   for (int ithermo=0;ithermo<NumThermoArrays;ithermo++){
     ThisCode.ThermoSTLArray.push_back(ThermoArray[ithermo]);
   }

   if (ThisCode.Nsites0=1) // Calculates Thermo from N=0 on (not N=-1)
     ThisCode.CalcThermo(ThermoArray,&Aeig);

   ThisCode.CalcStuff(&Aeig,&AeigCut);


   // Save Params?
   if (ThisCode.SaveData){ThisCode.SaveGenPars();}

   // Ok, with that we can go into the main loop

   // Parameters for HN (always the same?)
   ParamsHN.push_back(Lambda);
   for (int ich=0;ich<ThisCode.NumChannels;ich++)
     ParamsHN.push_back(chi_m1);

   // Entering calculation of H_Nsites0
   ThisCode.Nsites=ThisCode.Nsites0;

   while (ThisCode.Nsites<ThisCode.Nsitesmax){
     ThisCode.DN=HalfLambdaFactor*pow(Lambda,(-(ThisCode.Nsites-1)/2.0) );
     TM=ThisCode.DN/ThisCode.betabar;
     cout << "DN = " << ThisCode.DN << "TM = " << TM << endl;
     // Debugging...
     //bool disp=(ThisCode.Nsites==3?true:false);
     bool disp=false;

     // Build Basis

     //SingleSite.PrintQNumbers();
     //Aeig.PrintQNumbers();

     BuildBasis(CommonQNs, totSpos,&AeigCut,&Abasis, 
		&SingleSite,ThisCode.UpdateBefCut);

     //Abasis.PrintAll();

     // Diagonalize HN

     cout << "Diagonalizing HN... " << endl;
     cout << " chi_N = " << ParamsHN[1] << endl;

     // Need to adapt ALL models to use STLMatArray... 
     // in the meantime, this will work.
//      if ( (ThisCode.ModelNo==5)||(ThisCode.ModelNo==0)||(ThisCode.ModelNo==6) ){
//        HN.DiagHN(ParamsHN,&Abasis,&SingleSite,&STLMatArray[0],&Aeig,disp);
//      }else{
//        HN.DiagHN(ParamsHN,&Abasis,&SingleSite,MatArray,&Aeig);
//      }
     if ( (ThisCode.ModelNo==4) ){
       HN.DiagHN(ParamsHN,&Abasis,&SingleSite,MatArray,&Aeig);
     }else{
       HN.DiagHN(ParamsHN,&Abasis,&SingleSite,&STLMatArray[0],&Aeig,disp);
     }

     cout << "... done diagonalizing HN. " << endl;
    
     Aeig.PrintEn();

  

     //       for (int ibl=0;ibl<Aeig.NumBlocks();ibl++)
     // 	Aeig.PrintBlock(ibl);

     // Calculate Susceptibility/Entropy (does not need update)

     ThisCode.CalcThermo(ThermoArray,&Aeig);

     // Update matrices (if UpdateBefCut)
     // Cut and Update AC (If it does not need updated matrices, this should be cut). 
     // OLD
     //AeigCut.ClearAll();
     //AeigCut=CutStates(&Aeig, ThisCode.Ncutoff);

     // NEW (updates after cutting)
     CutStates(Aeig, AeigCut, ThisCode.Ncutoff);

      
     if (ThisCode.Nsites<ThisCode.Nsitesmax){
       // Need to adapt ALL models to use STLMatArray... 
       // in the meantime, this will work.
       if ( (ThisCode.ModelNo==4) ){
	 UpdateMatrices(&SingleSite,&AeigCut, 
			&Abasis,MatArray, ThisCode.NumNRGmats);

       }else{
//  	 UpdateMatrices(&SingleSite,&AeigCut, 
//  			&Abasis,&STLMatArray[0],STLMatArray.size());
         // Debugging 
	 UpdateMatrices(&SingleSite,&AeigCut, 
			&Abasis,&STLMatArray[0],STLMatArray.size(),disp);
       }
     }
     // end update matrices
     // Q Sz tests (remove later)
     cout << " N = " << ThisCode.Nsites << endl;

     // Save stuff to files

     //// Read/Write test (remove later) ///////////

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


     // Update ParamsHN
     //ParamsHN[1]=chiN(ThisCode.Nsites,Lambda); // Old
     for (int ich=0;ich<ThisCode.NumChannels;ich++)
       //ParamsHN[ich+1]=chiN(ThisCode.Nsites,Lambda);
       // Try this...
       ParamsHN[ich+1]=chain1.GetChin(ThisCode.Nsites);
     
     // Update sites
     ThisCode.Nsites++;
     
   }
   // end Nsites loop


   // De-allocate MatArray 
   delete[] ThermoArray;
   delete[] MatArray;

   ThisCode.WrapUp();

}
// end main
