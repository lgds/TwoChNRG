
#include <iostream>
#include <vector>
#include <math.h>


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "NRGOpMatRules.hpp"
#include "TwoChQS.hpp"
#include "NRG_main.hpp"


#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif

///////////////////
// Oct 2013.
// These routines are being called by CodeHandler_ModelSwitch
// In principle, they should be moved out of "Main" and into the
// src dir. Of course, this is a major move so I will need a couple
// of days to fix this.
//
///////////////////

///////////////////////////
///                     ///
///  1ch Anderson Hm1   ///
///                     ///
///  | Q > basis        ///
///                     ///
///////////////////////////



void OneChQ_SetAnderson_Hm1_old(vector<double> Params,
				CNRGarray* pAeig, 
				CNRGmatrix* NRGMats){
 // Set initial CNRG array (N=-1)

  double U=Params[0];
  double ed=Params[1];
  double Lambda=Params[2];
  double HalfLambdaFactor=Params[3];

  
  // First Test: the usual Anderson model

  pAeig->NQNumbers=1;
  pAeig->Nshell=-1;

  // new thing
  pAeig->totalS=false;
  

  pAeig->QNumbers.push_back(-1.0);
  pAeig->iDegen.push_back(0); // Sz=0  


  // Two states in this one
  pAeig->QNumbers.push_back(0.0);
  pAeig->iDegen.push_back(5); // Sz=1/2
  pAeig->iDegen.push_back(-5); // Sz=-1/2


  pAeig->QNumbers.push_back(1.0);
  pAeig->iDegen.push_back(0); // Sz=0

  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor));


  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);

  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(1.0);


  pAeig->dEigVec.push_back(1.0);

  pAeig->PrintEn();

  // Set Matrix elements

  // fN_up

  NRGMats[0].SyncNRGarray(*pAeig);
  NRGMats[0].UpperTriangular=false;

 
  NRGMats[0].MatEl.push_back(1.0);
  NRGMats[0].MatEl.push_back(0.0);
  NRGMats[0].MatBlockMap.push_back(0);
  NRGMats[0].MatBlockMap.push_back(1);
  NRGMats[0].MatBlockBegEnd.push_back(0);
  NRGMats[0].MatBlockBegEnd.push_back(1);

  NRGMats[0].MatEl.push_back(0.0);
  NRGMats[0].MatEl.push_back(1.0);
  NRGMats[0].MatBlockMap.push_back(1);
  NRGMats[0].MatBlockMap.push_back(2);
  NRGMats[0].MatBlockBegEnd.push_back(2);
  NRGMats[0].MatBlockBegEnd.push_back(3);


  // fN_dn
  NRGMats[1].SyncNRGarray(*pAeig);

  NRGMats[1].UpperTriangular=false;
 
  NRGMats[1].MatEl.push_back(0.0);
  NRGMats[1].MatEl.push_back(1.0);
  NRGMats[1].MatBlockMap.push_back(0);
  NRGMats[1].MatBlockMap.push_back(1);
  NRGMats[1].MatBlockBegEnd.push_back(0);
  NRGMats[1].MatBlockBegEnd.push_back(1);

  NRGMats[1].MatEl.push_back(-1.0);
  NRGMats[1].MatEl.push_back(0.0);
  NRGMats[1].MatBlockMap.push_back(1);
  NRGMats[1].MatBlockMap.push_back(2);
  NRGMats[1].MatBlockBegEnd.push_back(2);
  NRGMats[1].MatBlockBegEnd.push_back(3);


  cout << " OLD: Setting Sz... " << endl;

  NRGMats[2].SyncNRGarray(*pAeig);
  NRGMats[3].SyncNRGarray(*pAeig);

  int i1=0;
  for (int ibl=0; ibl<pAeig->NumBlocks(); ibl++)
    {
      int ist0=pAeig->GetBlockLimit(ibl,0);
      int ist1=pAeig->GetBlockLimit(ibl,1);

      double Qi=pAeig->GetQNumber(ibl,0);

      cout << " Q = " << Qi
	   << " ist0 = " << ist0
	   << " ist1 = " << ist1 << endl;

      // Block diagonal
      NRGMats[2].MatBlockMap.push_back(ibl);
      NRGMats[2].MatBlockMap.push_back(ibl);
      NRGMats[3].MatBlockMap.push_back(ibl);
      NRGMats[3].MatBlockMap.push_back(ibl);
      // Position in MatBlock follows ist
      NRGMats[2].MatBlockBegEnd.push_back(i1);
      NRGMats[3].MatBlockBegEnd.push_back(i1);
      for (int ist=ist0;ist<=ist1;ist++)
	{
	  double Szi=(double)(pAeig->iDegen[ist])/10.0;
	  for (int jst=ist0;jst<=ist1;jst++)
	    {
	      // Sz and Sz2 diagonal only (waste of space but...)
	      if (ist==jst)
		{
		  NRGMats[2].MatEl.push_back(Szi);
		  NRGMats[3].MatEl.push_back(Szi*Szi);
		  // Test in Sep 08: calc ndot instead WORKS!
// 		  NRGMats[2].MatEl.push_back(Qi+1.0);
// 		  NRGMats[3].MatEl.push_back(Qi+1.0);
		}
	      else
		{
		  NRGMats[2].MatEl.push_back(0.0);
		  NRGMats[3].MatEl.push_back(0.0);
		}
	      i1++;
	    } // loop in jst
	} //ist
      NRGMats[2].MatBlockBegEnd.push_back(i1-1);
      NRGMats[3].MatBlockBegEnd.push_back(i1-1);

    }
  // end Set-up of Sz,Sz2


  //NRGMats[2].PrintAllBlocks();
  //NRGMats[3].PrintAllBlocks();

  // Suscep test
    
  CNRGbasisarray AeigCut=CutStates(pAeig, 1000);

  vector<double> ParamsBetabar;
  double betabar=0.727;
  ParamsBetabar.clear();
  ParamsBetabar.push_back(betabar);
  
  double dSz0=CalcOpAvg(ParamsBetabar,&AeigCut,&NRGMats[2],false,0);
  double dSz20=CalcOpAvg(ParamsBetabar,&AeigCut,&NRGMats[3],false,0);
  
  cout << " N=-1: Sz = " << dSz0 << " Sz2 = " << dSz20 
       << "    T_M chi = " << dSz20-dSz0*dSz0 << endl;


}
/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////


void OneChQ_SetAnderson_Hm1(vector<double> Params,
			    CNRGarray* pAeig, 
			    vector<CNRGmatrix> &STLNRGMats){

  // Set initial CNRG array (N=-1)

  double U=Params[0];
  double ed=Params[1];
  double Lambda=Params[2];
  double HalfLambdaFactor=Params[3];

  
  // First Test: the usual Anderson model

  pAeig->NQNumbers=1;
  pAeig->Nshell=-1;

  // new thing
  pAeig->totalS=false;
  

  pAeig->QNumbers.push_back(-1.0);
  pAeig->iDegen.push_back(0); // Sz=0  


  // Two states in this one
  pAeig->QNumbers.push_back(0.0);
  pAeig->iDegen.push_back(5); // Sz=1/2
  pAeig->iDegen.push_back(-5); // Sz=-1/2


  pAeig->QNumbers.push_back(1.0);
  pAeig->iDegen.push_back(0); // Sz=0

  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor));


  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);

  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(1.0);


  pAeig->dEigVec.push_back(1.0);

  pAeig->PrintEn();

  // Set Matrix elements

  // fN_up

  STLNRGMats[0].SyncNRGarray(*pAeig);
  STLNRGMats[0].UpperTriangular=false;

 
  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatEl.push_back(0.0);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);

  STLNRGMats[0].MatEl.push_back(0.0);
  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockMap.push_back(2);
  STLNRGMats[0].MatBlockBegEnd.push_back(2);
  STLNRGMats[0].MatBlockBegEnd.push_back(3);


  // fN_dn
  STLNRGMats[1].SyncNRGarray(*pAeig);

  STLNRGMats[1].UpperTriangular=false;
 
  STLNRGMats[1].MatEl.push_back(0.0);
  STLNRGMats[1].MatEl.push_back(1.0);
  STLNRGMats[1].MatBlockMap.push_back(0);
  STLNRGMats[1].MatBlockMap.push_back(1);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);

  STLNRGMats[1].MatEl.push_back(-1.0);
  STLNRGMats[1].MatEl.push_back(0.0);
  STLNRGMats[1].MatBlockMap.push_back(1);
  STLNRGMats[1].MatBlockMap.push_back(2);
  STLNRGMats[1].MatBlockBegEnd.push_back(2);
  STLNRGMats[1].MatBlockBegEnd.push_back(3);

  if (STLNRGMats.size()>5){
    cout << " Setting up cd_up and cd_dn... ";
    STLNRGMats[4].CopyData(&STLNRGMats[0]); // cd_up
    STLNRGMats[5].CopyData(&STLNRGMats[1]); // cd_dn
    cout << " ... done." << endl;

  }


  cout << " Setting Sz... " << endl;

  STLNRGMats[2].SyncNRGarray(*pAeig);
  STLNRGMats[3].SyncNRGarray(*pAeig);

  int i1=0;
  for (int ibl=0; ibl<pAeig->NumBlocks(); ibl++)
    {
      int ist0=pAeig->GetBlockLimit(ibl,0);
      int ist1=pAeig->GetBlockLimit(ibl,1);

      double Qi=pAeig->GetQNumber(ibl,0);

      cout << " Q = " << Qi
	   << " ist0 = " << ist0
	   << " ist1 = " << ist1 << endl;

      // Block diagonal
      STLNRGMats[2].MatBlockMap.push_back(ibl);
      STLNRGMats[2].MatBlockMap.push_back(ibl);
      STLNRGMats[3].MatBlockMap.push_back(ibl);
      STLNRGMats[3].MatBlockMap.push_back(ibl);
      // Position in MatBlock follows ist
      STLNRGMats[2].MatBlockBegEnd.push_back(i1);
      STLNRGMats[3].MatBlockBegEnd.push_back(i1);
      for (int ist=ist0;ist<=ist1;ist++)
	{
	  double Szi=(double)(pAeig->iDegen[ist])/10.0;
	  for (int jst=ist0;jst<=ist1;jst++)
	    {
	      // Sz and Sz2 diagonal only (waste of space but...)
	      if (ist==jst)
		{
		  STLNRGMats[2].MatEl.push_back(Szi);
		  STLNRGMats[3].MatEl.push_back(Szi*Szi);
		  // Test in Sep 08: calc ndot instead WORKS!
// 		  STLNRGMats[2].MatEl.push_back(Qi+1.0);
// 		  STLNRGMats[3].MatEl.push_back(Qi+1.0);
		}
	      else
		{
		  STLNRGMats[2].MatEl.push_back(0.0);
		  STLNRGMats[3].MatEl.push_back(0.0);
		}
	      i1++;
	    } // loop in jst
	} //ist
      STLNRGMats[2].MatBlockBegEnd.push_back(i1-1);
      STLNRGMats[3].MatBlockBegEnd.push_back(i1-1);

    }
  // end Set-up of Sz,Sz2


  //STLNRGMats[2].PrintAllBlocks();
  //STLNRGMats[3].PrintAllBlocks();

  // Suscep test
    
  CNRGbasisarray AeigCut=CutStates(pAeig, 1000);

  vector<double> ParamsBetabar;
  double betabar=0.727;
  ParamsBetabar.clear();
  ParamsBetabar.push_back(betabar);
  
  double dSz0=CalcOpAvg(ParamsBetabar,&AeigCut,&STLNRGMats[2],false,0);
  double dSz20=CalcOpAvg(ParamsBetabar,&AeigCut,&STLNRGMats[3],false,0);
  
  cout << " N=-1: Sz = " << dSz0 << " Sz2 = " << dSz20 
       << "    T_M chi = " << dSz20-dSz0*dSz0 << endl;



}
///////////////////
///////////////////
///////////////////
///////////////////

///////////////////////////
///                     ///
///  1ch Anderson Hm1   ///
///                     ///
///  | Q Sz > basis     ///
///                     ///
///////////////////////////


////
//// | Q Sz > basis
////


void OneChQSz_SetAnderson_Hm1(vector<double> Params,
			    CNRGarray* pAeig, 
			    vector<CNRGmatrix> &STLNRGMats){


  // Set initial CNRG array (N=-1)

  double U=Params[0];
  double ed=Params[1];
  double Lambda=Params[2];
  double HalfLambdaFactor=Params[3];

  double hz=0.5*Params[4]; // hz=g mu_B Bz: Zeeman term is:  - hz S_z 
                           // 0.5 added in Feb 2014
  // First Test: the usual Anderson model

  pAeig->NQNumbers=2;
  pAeig->Nshell=-1;
  // new thing
  pAeig->totalS=false;
  

  // |0> = |Q=-1, Sz=0>
  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(0.0);
  pAeig->iDegen.push_back(0); // Sz=0  

  // |up> = |Q=0, Sz=1/2>
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.5);
  pAeig->iDegen.push_back(0); // Sz=1/2

  // |dn> = |Q=0, Sz=-1/2>
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(-0.5);
  pAeig->iDegen.push_back(0); // Sz=1/2

  // |up dn> = |Q=1, Sz=0>
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(0.0);
  pAeig->iDegen.push_back(0); // Sz=0

  // One state per block
  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U-hz)/(Lambda*HalfLambdaFactor)); // up
  pAeig->dEn.push_back((ed+0.5*U+hz)/(Lambda*HalfLambdaFactor)); // dn
  pAeig->dEn.push_back((2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor));


  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);

  pAeig->PrintEn();

  // Set Matrix elements

  // fN_up

  STLNRGMats[0].SyncNRGarray(*pAeig);
  STLNRGMats[0].UpperTriangular=false;

  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);

  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(2);
  STLNRGMats[0].MatBlockMap.push_back(3);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);

  // fN_dn
  STLNRGMats[1].SyncNRGarray(*pAeig);
  STLNRGMats[1].UpperTriangular=false;
 
  STLNRGMats[1].MatEl.push_back(1.0);
  STLNRGMats[1].MatBlockMap.push_back(0);
  STLNRGMats[1].MatBlockMap.push_back(2);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);

  STLNRGMats[1].MatEl.push_back(-1.0);
  STLNRGMats[1].MatBlockMap.push_back(1);
  STLNRGMats[1].MatBlockMap.push_back(3);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);

  // Add stuff here to calculate Sz and nd
  // Needs to differentiate calcdens (HOW??)
  // By matrix names

  if (	(strcmp(STLNRGMats[2].MatName,"Ndot")==0)&&
	(strcmp(STLNRGMats[3].MatName,"Szdot")==0) ){

    //////////////////

    cout << " Setting Nd and Sz... " << endl;

    STLNRGMats[2].SyncNRGarray(*pAeig);
    STLNRGMats[3].SyncNRGarray(*pAeig);

    int i1=0;
    for (int ibl=0; ibl<pAeig->NumBlocks(); ibl++){
      double Qi=pAeig->GetQNumber(ibl,0);
      double Szi=pAeig->GetQNumber(ibl,1);

      cout << " Q = " << Qi
	   << " Sz = " << Szi 
	   << endl;

      // Block diagonal
      STLNRGMats[2].MatBlockMap.push_back(ibl);
      STLNRGMats[2].MatBlockMap.push_back(ibl);
      STLNRGMats[3].MatBlockMap.push_back(ibl);
      STLNRGMats[3].MatBlockMap.push_back(ibl);
      // Position in MatBlock follows ist
      STLNRGMats[2].MatBlockBegEnd.push_back(i1);
      STLNRGMats[3].MatBlockBegEnd.push_back(i1);
      

      // Only one state per block !
      STLNRGMats[2].MatEl.push_back(Qi+1.0); // <Nd>
      STLNRGMats[3].MatEl.push_back(Szi);    // <Sz>
      i1++;

      STLNRGMats[2].MatBlockBegEnd.push_back(i1-1);
      STLNRGMats[3].MatBlockBegEnd.push_back(i1-1);

    }
    // end loop in blocks (Set-up of Nd,Sz)

    STLNRGMats[2].PrintAllBlocks();
    STLNRGMats[3].PrintAllBlocks();

  }else if (strcmp(STLNRGMats[2].MatName,"Szomega")==0){
      cout << " Dynamical spin susceptibiliy: setting up <Sz>_N=-1  ";
      // Needs work here

      STLNRGMats[2].SyncNRGarray(*pAeig);
      STLNRGMats[3].SyncNRGarray(*pAeig);
      STLNRGMats[4].SyncNRGarray(*pAeig);
      STLNRGMats[5].SyncNRGarray(*pAeig);

      int i1=0;
      for (int ibl=0; ibl<pAeig->NumBlocks(); ibl++){
	double Qi=pAeig->GetQNumber(ibl,0);
	double Szi=pAeig->GetQNumber(ibl,1);

	cout << " Q = " << Qi
	     << " Sz = " << Szi 
	     << endl;

	// Sz dynamical and Nd dynamical
	for (int imat=2;imat<=3; imat++){
	// Block diagonal
	STLNRGMats[imat].MatBlockMap.push_back(ibl);
	STLNRGMats[imat].MatBlockMap.push_back(ibl);
	// Position in MatBlock follows ist
	// Only one state per block !
	STLNRGMats[imat].MatBlockBegEnd.push_back(i1);
      	STLNRGMats[imat].MatBlockBegEnd.push_back(i1);
	}
	STLNRGMats[2].MatEl.push_back(Szi); // <Sz>
	STLNRGMats[3].MatEl.push_back(Qi+1.0); // <Nd>
	i1++;
      }
      // end loop in blocks (Set-up of Nd,Sz,Sz2)

      //Setting Sz2 and Nd2
      for (int imat=4; imat<=5; imat++){
      STLNRGMats[imat].CopyData(&STLNRGMats[imat-2]);
      for (int ii=0;ii<STLNRGMats[imat].MatEl.size();ii++)
	STLNRGMats[imat].MatEl[ii]*=STLNRGMats[imat].MatEl[ii];
      }
      // end setting Sz2 and Nd2

      cout << " ... done." << endl;


//       STLNRGMats[2].PrintAllBlocks();
//       STLNRGMats[3].PrintAllBlocks();

  }else if ((STLNRGMats.size()==6)&&(strcmp(STLNRGMats[4].MatName,"cdupNdn"))==0){
      cout << " Spectral functions: Self-energy trick  ";
      // Needs work here
      STLNRGMats[2].CopyData(&STLNRGMats[0]);  //c1_up
      STLNRGMats[3].CopyData(&STLNRGMats[1]);  //c1_dn
      //cd_up Ndn
      STLNRGMats[4].CopyData(&STLNRGMats[0]); //c1_up Ndn
      STLNRGMats[4].MatEl[0]=0.0;  // <0|cd_up Ndn|up>=0
      //cd_dn
      STLNRGMats[5].CopyData(&STLNRGMats[1]);  //c1_dn Nup
      STLNRGMats[5].MatEl[0]=0.0;  // <0|cd_dn Nup|dn>=0

  }else {
  // end if Mat[1]=Szomega
    switch(STLNRGMats.size()){
    case 2: //fd_up fd_dn -> calcdens=2
      cout << " Thermodynamics only." << endl;
      break;
    case 4: // calcdens=1
      // This should work
      // cd_up
      STLNRGMats[2].CopyData(&STLNRGMats[0]);  //c1_up
      //cd_dn
      STLNRGMats[3].CopyData(&STLNRGMats[1]);  //c1_dn
      break;
    default:
      cout << " Can't figure out what to do with the STLNRGMats... Exiting." << endl;
      exit(0);
    }
    // end switch
  }
  // if calcdens
}
// end set H0 for Anderson with Q Sz symmetry

///////////////////


///////////////////////////
///                     ///
///  1chQSz DQD H0      ///
///                     ///
///  | Q Sz > basis      ///
///                     ///
///////////////////////////

void OneChQSz_SetH0_DQD(vector<double> Params, 
		       CNRGarray* pAeig,
			CNRGbasisarray* pSingleSite, // Do we need this???
		       vector<CNRGmatrix> &STLNRGMats){

  // Set initial CNRG array (N=0)
  // Actually, it will set H0 at the end of the day

  // Includes lambda AND Zeeman in BOTH dots.

  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];


  double U1tilde=Params[2]/(sqrt(Lambda)*HalfLambdaFactor);
  double ed1tilde=Params[3]/(sqrt(Lambda)*HalfLambdaFactor);
  double gamma1tilde=Params[4];
//  double hz1=Params[5]; // hz=g mu_B Bz: Zeeman term is:  - hz S_z 
  double hz1=0.5*Params[5]/(sqrt(Lambda)*HalfLambdaFactor);; // hz=g mu_B Bz: Zeeman term is:  - hz S_z 

  // Remember:
  // gamma~=(2*Gamma/pi)^1/2/(sqrt(Lambda)*HalfLambdaFactor)
  // Why? see Wilson's paper (E. 2.18 with N=-1)
  // Also: Multiplied U1~ and ed~ by sqrt(Lambda) as compared to previous codes
  // (i.e. divide by sqrt(Lambda) instead of Lambda)
  // Thus, no lambda^{1/2} factor in setting up H0.

  double U2tilde=Params[6]/(sqrt(Lambda)*HalfLambdaFactor);
  double ed2tilde=Params[7]/(sqrt(Lambda)*HalfLambdaFactor);
  double gamma2tilde=Params[8];
//  double hz2=Params[9]; // hz=g mu_B Bz: Zeeman term is:  - hz S_z 
  double hz2=0.5*Params[9]/(sqrt(Lambda)*HalfLambdaFactor); // hz=g mu_B Bz: Zeeman term is:  - hz S_z 

  // Is the factor right???
  double small_lambda=Params[10]/(sqrt(Lambda)*HalfLambdaFactor);

  // The Basis: a two Channel-like basis in QSz basis 
  // with Sz entering iDegen (Sz=10*iDegen)


  cout << " U1~= " << U1tilde
       << " ed1~= " << ed1tilde
       << " t1~= " << gamma1tilde
       << " hz1~= " << hz1 << endl
       << " U2~= " << U2tilde
       << " ed2~= " << ed2tilde
       << " t2~= " << gamma2tilde
       << " hz2~= " << hz2
       << " lambda= " << small_lambda
       << endl;

    
  CNRGbasisarray AbasisHm1; // will be AbasisHm1


  // Set AbasisHm1 as a 2-ch QSz site
  TwoChQSz_SetSingleSite(&AbasisHm1);

  AbasisHm1.Nshell=-1;
  // Set energies (diagonal)

  //|0>i = 0.5*Ui
  //|up>i = ei+0.5*Ui-hzi*Sz
  //|dn>i = ei+0.5*Ui+hzi*Sz
  //|up dn>i = 2ei+1.5*Ui

  AbasisHm1.dEn.push_back(0.5*U1tilde+0.5*U2tilde); // | 0 0 >

  AbasisHm1.dEn.push_back(0.5*U1tilde+ed2tilde+0.5*U2tilde-hz2); // | 0 up>
  AbasisHm1.dEn.push_back(0.5*U2tilde+ed1tilde+0.5*U1tilde-hz1); // | up 0>

  AbasisHm1.dEn.push_back(0.5*U1tilde+ed2tilde+0.5*U2tilde+hz2); // | 0 dn>
  AbasisHm1.dEn.push_back(0.5*U2tilde+ed1tilde+0.5*U1tilde+hz1); // | dn 0>

  //|up dn> 
  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde-hz1+ed2tilde+0.5*U2tilde+hz2); 
  AbasisHm1.dEn.push_back(0.5*U1tilde+2.0*ed2tilde+1.5*U2tilde); //|0  up dn>  
  AbasisHm1.dEn.push_back(0.5*U2tilde+2.0*ed1tilde+1.5*U1tilde); //|up dn  0>  
  //|dn up>
  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+hz1+ed2tilde+0.5*U2tilde-hz2); 

  //|up up>
  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde-hz1+ed2tilde+0.5*U2tilde-hz2);  

  //|dn dn>
  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+hz1+ed2tilde+0.5*U2tilde+hz2);  
  
  //|up  updn>
  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde-hz1+2.0*ed2tilde+1.5*U2tilde);
  //|updn  up>
  AbasisHm1.dEn.push_back(ed2tilde+0.5*U2tilde-hz2+2.0*ed1tilde+1.5*U1tilde);

  //|dn  updn>
  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+hz1+2.0*ed2tilde+1.5*U2tilde);
  //|updn  dn>
  AbasisHm1.dEn.push_back(ed2tilde+0.5*U2tilde+hz2+2.0*ed1tilde+1.5*U1tilde);

  //|updn  updn >  
  AbasisHm1.dEn.push_back(2.0*ed1tilde+1.5*U1tilde+2.0*ed2tilde+1.5*U2tilde); 

  AbasisHm1.SetEigVecToOne();

  // Set up cd1, cd2

  for (int imat=0;imat<STLNRGMats.size();imat++)
    STLNRGMats[imat].SyncNRGarray(AbasisHm1);


  vector<CNRGmatrix> AuxMatArray; // Arrays needed to set up H0 (c1, c2)
  CNRGmatrix auxMat;

  auxMat=STLNRGMats[0]; // Same rules as fN_up
  auxMat.SetZeroMatrix(); // Build block structure
  AuxMatArray.push_back(auxMat); // cd1_up
  auxMat=STLNRGMats[1]; // Same rules as fN_dn
  auxMat.SetZeroMatrix(); // Build block structure
  AuxMatArray.push_back(auxMat); // cd1_dn
  auxMat=STLNRGMats[0]; // Same rules as fN_up
  auxMat.SetZeroMatrix(); // Build block structure
  AuxMatArray.push_back(auxMat); // cd2_up
  auxMat=STLNRGMats[1]; // Same rules as fN_dn
  auxMat.SetZeroMatrix(); // Build block structure
  AuxMatArray.push_back(auxMat); // cd2_dn


  double auxEl[4]={0.0};

  // Need to work on this...

  for (int ist=0;ist<16;ist++){
    for (int jst=0;jst<16;jst++){
      auxEl[0]=TwoCh_fd_table(1, 1, jst, ist); // c1_up (not dagger)
      auxEl[1]=TwoCh_fd_table(1, -1, jst, ist); // c1_dn (not dagger)
      auxEl[2]=TwoCh_fd_table(2, 1, jst, ist); // c2_up (not dagger)
      auxEl[3]=TwoCh_fd_table(2, -1, jst, ist); // c2_dn (not dagger)

      AuxMatArray[0].PushMatEl(auxEl[0],ist,jst); 
      AuxMatArray[1].PushMatEl(auxEl[1],ist,jst); 
      AuxMatArray[2].PushMatEl(auxEl[2],ist,jst); 
      AuxMatArray[3].PushMatEl(auxEl[3],ist,jst); 
    }
  }
  // end push Mat else

  //c1_up
//    AuxMatArray[0].PrintAllBlocks();
//   AuxMatArray[0].PrintMatBlock(0,1);
  //c1_dn
//    AuxMatArray[0].PrintAllBlocks();
//   AuxMatArray[0].PrintMatBlock(0,1);


  // Add lambda term in Hm1, diagonalize it and set AeigHm1 

  // This should be ok.

  CNRGarray AeigHm1; // will be AeigHm1

  auxMat.ClearAll(); //  Hm1 (re-using auxMat)
  auxMat.NeedOld=false;
  auxMat.UpperTriangular=true;
  auxMat.CheckForMatEl=Diag_check;
  auxMat.CalcHNMatEl=TwoDotQSz_Hm1_MatEl;

  vector<double> ParamsHm1;
  ParamsHm1.push_back(small_lambda);

  // and diagonalize the Hm1 Hamiltonian

  auxMat.DiagHN(ParamsHm1,&AbasisHm1,pSingleSite,&AuxMatArray[0],&AeigHm1);
  
  AeigHm1.PrintEn();


  //c1_up
  //AuxMatArray[0].PrintAllBlocks();
  // AuxMatArray[0].PrintMatBlock(3,6);

  // c1_dn
  //AuxMatArray[1].PrintAllBlocks();
  //AuxMatArray[1].PrintMatBlock(1,3);

  // Rotate c1,c2 (in the UnCut basis!)

  for (int imat=0;imat<AuxMatArray.size();imat++){
    //RotateMatrix((&AuxMatArray[imat]),(&AeigHm1),&auxMat);
    RotateMatrix_NoCut((&AuxMatArray[imat]),(&AeigHm1),&auxMat,1);
    AuxMatArray[imat].CopyData(&auxMat);
  }
  // end loop in matrices

  //c1_up
  //AuxMatArray[0].PrintAllBlocks();
//   AuxMatArray[0].PrintMatBlock(3,6);

  // c1_dn
  //AuxMatArray[1].PrintAllBlocks();
//   AuxMatArray[1].PrintMatBlock(3,7);


  //c2_up
  //AuxMatArray[2].PrintAllBlocks();
  //AuxMatArray[2].PrintMatBlock(3,6);

  // c2_dn
  //AuxMatArray[3].PrintAllBlocks();
  //AuxMatArray[3].PrintMatBlock(3,7);



//   for (int imat=0;imat<AuxMatArray.size();imat++){
//      cout << " ---- cd" << imat+1 << "  ---- " << endl; 
//      AuxMatArray[imat].PrintAllBlocks();
//   }

  
  // Set up operators in the Hm1 basis
  // ndot, Sz, S1 dot S2 if calcdens==0

  switch(STLNRGMats.size()){
  case 2: //fd_up fd_dn
    cout << " Thermodynamics only." << endl;
    break;
  case 6: //fd_up fd_dn c1_up c1_dn c2_up c2_dn
    // cd1, cd2 if calcdens==1 (just use the ops defined above)
    STLNRGMats[2].CopyData(&AuxMatArray[0]);  //c1_up
    STLNRGMats[3].CopyData(&AuxMatArray[1]);  //c1_dn
    STLNRGMats[4].CopyData(&AuxMatArray[2]);  //c2_up
    STLNRGMats[5].CopyData(&AuxMatArray[3]);  //c2_dn

    cout << " DQD: spectral functions still at implementation stage..." << endl;
    break;
  default:
    cout << " Can't figure out what to do with the STLNRGMats... Exiting." << endl;
    exit(0);
  }
  // end switch STLNRGMats.size()

  
//   for (int imat=0;imat<STLNRGMats.size();imat++){
//      cout << " ---- Matrix " << imat << "  ---- " << endl; 
//      STLNRGMats[imat].PrintAllBlocks();
//   }



  // Connect to first site of the chain: Create new basis and so on. 
  // Build basis for DQD + 1 site

  vector<int> CommonQNs;
  vector<int> totSpos;
  CNRGbasisarray AbasisH0;
  CommonQNs.push_back(2); // Q and Sz are common QNs
  CommonQNs.push_back(0); // position of Q in old basis
  CommonQNs.push_back(1); // position of Sz in old basis
  CommonQNs.push_back(0); // position of Q in SingleSite
  CommonQNs.push_back(1); // position of Sz in SingleSite

  // AbasisHm1 becomes "AcutHm1"
  AbasisHm1.FalseCut(&AeigHm1);

  // Use "AcutHm1" to build the H0 basis
  BuildBasis(CommonQNs,totSpos,&AbasisHm1,&AbasisH0,pSingleSite,0);


  // Diagonalize H0: Aeig will be the new vector

  auxMat.ClearAll(); //  H0 (re-using auxMat)
  auxMat.NeedOld=false;
  auxMat.UpperTriangular=true;
  auxMat.CheckForMatEl=Diag_check;
  auxMat.CalcHNMatEl=OneChQSz_H0DQD_MatEl;

  vector<double> ParamsH0;
  ParamsH0.push_back(gamma1tilde);
  ParamsH0.push_back(gamma2tilde);

  auxMat.DiagHN(ParamsH0,&AbasisH0,pSingleSite,&AuxMatArray[0],pAeig);
  
  pAeig->PrintEn();


  // Update all operators return all.

  // AbasisHm1 becomes "Acut"
  AbasisHm1.FalseCut(pAeig);
  
  UpdateMatrices(pSingleSite,&AbasisHm1, 
		 &AbasisH0,&STLNRGMats[0],STLNRGMats.size());

//   for (int imat=0;imat<STLNRGMats.size();imat++){
//      cout << " ---- Matrix " << imat << "  ---- " << endl; 
//      STLNRGMats[imat].PrintAllBlocks();
//   }

//  for (int imatbl=0;imatbl<5;imatbl++){
//   for (int imat=0;imat<2;imat++){
//      cout << " ---- Matrix " << imat << "  ---- " << endl;
//      STLNRGMats[imat].PrintMatBlock(imatbl);
//   }
//   //end loop in matrices
//  }
//  //end loop in blocks

  cout << " ...  1chQSz_DQD H0 done. " << endl;


// Oct 11 - Stopped here


}

/////////////////






///////////////////
///////////////////
///////////////////




///////////////////////////
///                     ///
///  1ch Q SMM H0       ///
///                     ///
///                     ///
///////////////////////////




///////////////////
// void OneChQ_SetSMM_H0(vector<double> Params,
// 		       CNRGarray* pAeig,
// 		       CNRGbasisarray* pSingleSite,
// 		       CNRGmatrix* NRGMats){

void OneChQ_SetSMM_H0(vector<double> Params,
		       CNRGarray* pAeig,
		       CNRGbasisarray* pSingleSite,
		       vector<CNRGmatrix> &STLNRGMats){


  // Set initial CNRG array (N=-1)
  

  // The Basis: a two Channel-like basis in QS basis 
  // with Sz entering iDegen (Sz=10*iDegen)

  CNRGbasisarray AbasisH0; // will be BOTH AbasisHm1 and AbasisH0

  // Set AbasisHm1
  TwoChQSNoSz_SetSingleSite(&AbasisH0);
  AbasisH0.Nshell=-1;

  // Set Sz, Sz2 (Sz is given in iDegen in AbasisH0)

  cout << " Setting Sz and Sz2... " << endl;

  STLNRGMats[2].SyncNRGarray(AbasisH0);
  STLNRGMats[3].SyncNRGarray(AbasisH0);

  int i1=0;
  for (int ibl=0; ibl<AbasisH0.NumBlocks(); ibl++){
    int ist0=AbasisH0.GetBlockLimit(ibl,0);
    int ist1=AbasisH0.GetBlockLimit(ibl,1);

    // Block diagonal
    STLNRGMats[2].MatBlockMap.push_back(ibl);
    STLNRGMats[2].MatBlockMap.push_back(ibl);
    STLNRGMats[3].MatBlockMap.push_back(ibl);
    STLNRGMats[3].MatBlockMap.push_back(ibl);
    // Position in MatBlock follows ist
    STLNRGMats[2].MatBlockBegEnd.push_back(i1);
    STLNRGMats[3].MatBlockBegEnd.push_back(i1);
    for (int ist=ist0;ist<=ist1;ist++){
      double Szi=(double)(AbasisH0.iDegen[ist])/10.0;
      for (int jst=ist0;jst<=ist1;jst++){
	// Sz and Sz2 diagonal only (waste of space but...)
	if (ist==jst){
	  STLNRGMats[2].MatEl.push_back(Szi);
	  STLNRGMats[3].MatEl.push_back(Szi*Szi);
	}
	else{
	  STLNRGMats[2].MatEl.push_back(0.0);
	  STLNRGMats[3].MatEl.push_back(0.0);
	}
	i1++;
      } // loop in jst
    } //ist
    STLNRGMats[2].MatBlockBegEnd.push_back(i1-1);
    STLNRGMats[3].MatBlockBegEnd.push_back(i1-1);
  }
  // end Set-up of Sz,Sz2

  // Set cd1_up, cd1_dn, cd2_up, cd2_dn

  if (STLNRGMats.size()>4){
    cout << " Setting up cd1 and cd2... " << endl;
    
    // 16 state QS basis (Sz is set in iDegen)
    for (int imat=4;imat<STLNRGMats.size();imat++){
      STLNRGMats[imat].SyncNRGarray(AbasisH0);
      STLNRGMats[imat].NeedOld=true; // just to make sure
      STLNRGMats[imat].CalcMatEl=AlwaysOne_MatEl; // for now
    }

    int icount[4]={0};
    for (int ibl=0; ibl<AbasisH0.NumBlocks(); ibl++){
      int ist0=AbasisH0.GetBlockLimit(ibl,0);
      int ist1=AbasisH0.GetBlockLimit(ibl,1);
      
      for (int iblp=0; iblp<AbasisH0.NumBlocks(); iblp++){
	int istp0=AbasisH0.GetBlockLimit(iblp,0);
	int istp1=AbasisH0.GetBlockLimit(iblp,1);
	

       	double auxMatEl=0.0;
	if ( TwoChQS_cd_check(&AbasisH0,ibl,iblp) ){ // same for all ops.
	  for (int iop=4;iop<STLNRGMats.size();iop++){
	    STLNRGMats[iop].MatBlockMap.push_back(ibl);
	    STLNRGMats[iop].MatBlockMap.push_back(iblp);
	    STLNRGMats[iop].MatBlockBegEnd.push_back(icount[iop-4]);
	    int ich=1;
	    int isigma=1;
	    switch(iop){
	    case 4:
	      ich=1;isigma=1;
	      break;
	    case 5:
	      ich=1;isigma=-1;
	      break;
	    case 6:
	      ich=2;isigma=1;
	      break;
	    case 7:
	      ich=2;isigma=-1;
	      break;
	    }
	    // end switch
	    // <ist| c_{idot sigma} | istp>=<istp| c^{dagger}_{idot sigma} | ist>
	    for (int ist=ist0;ist<=ist1;ist++){
	    for (int istp=istp0;istp<=istp1;istp++){
	      auxMatEl=TwoChQS_fd_table(ich,isigma, istp, ist); // c_{ich isigma}
	      STLNRGMats[iop].MatEl.push_back(auxMatEl);
	      icount[iop-4]++;
// 	      if (iop==7){
// 		cout << " icount = " << icount[iop-4]
// 		     << " ist = " << ist
// 		     << " istp = " << istp
// 		     << " matel = " << auxMatEl
// 		     << endl;
// 	      }
	    }	    
	    }
	    // end loop in block states
	    STLNRGMats[iop].MatBlockBegEnd.push_back(icount[iop-4]-1);
	  }
	  // end loop in operators 		
	}
	// end if there is a matrix element there
      }
      // end loop in iblp (AbasisH0 blocks)
    }
    // end loop in ibl (AbasisH0 blocks)

    cout << " ... done " << endl;
    //STLNRGMats[7].PrintAllBlocks();

  }
  // end if calculating dynamics, setting cd1_up, cd1_dn, cd2_up, cd2_dn

//   int ii=11;
//   int jj=11;
//   cout << " <"<<ii<<"|H0|"<<jj<<"> = " 
//        << H0.CalcHNMatEl(Params,&AbasisH0,pEmptySite,pEmptyArray,ii,jj)
//        << endl;


  // Set Up SMM Hamiltonian

  CNRGmatrix H0(AbasisH0); // will be BOTH Hm1 and H0

  CNRGmatrix* pEmptyArray;
  CNRGbasisarray* pEmptySite;

  // Set up and diag Hm1
  H0.NeedOld=false;
  H0.UpperTriangular=true;
  H0.CheckForMatEl=Diag_check;
  H0.CalcHNMatEl=TwoChQSNoSz_Hm1SMM_MatEl;


  H0.DiagHN(Params,&AbasisH0,pEmptySite,pEmptyArray,pAeig);

  pAeig->PrintEn();

  // Next step:
  // Build the Hamiltonian for N=0. 
  // There will be FOUR operators, spin up and down for each site.

  CNRGmatrix* fHm1Array; 

  // Allocate Matrices 
  fHm1Array = new CNRGmatrix [4];

  //fm1_dot1up
  fHm1Array[0].CheckForMatEl=TwoChQS_cd_check;
  fHm1Array[0].CalcMatEl=TwoChQSNoSz_fm1_dot1up_MatEl;
  fHm1Array[0].NeedOld=false;
  //fm1_dot1dn
  fHm1Array[1].CheckForMatEl=TwoChQS_cd_check;
  fHm1Array[1].CalcMatEl=TwoChQSNoSz_fm1_dot1dn_MatEl;
  fHm1Array[1].NeedOld=false;
  //fm1_dot2up
  fHm1Array[2].CheckForMatEl=TwoChQS_cd_check;
  fHm1Array[2].CalcMatEl=TwoChQSNoSz_fm1_dot2up_MatEl;
  fHm1Array[2].NeedOld=false;
  //fm1_dot2dn
  fHm1Array[3].CheckForMatEl=TwoChQS_cd_check;
  fHm1Array[3].CalcMatEl=TwoChQSNoSz_fm1_dot2dn_MatEl;
  fHm1Array[3].NeedOld=false;


  // Update matrices to close the Hm1 case
    
  CNRGbasisarray AeigCut=CutStates(pAeig, 1000);
  
  UpdateMatrices(pEmptySite,&AeigCut, 
		 &AbasisH0,fHm1Array, 4);

  // Update Sz, Sz2 in Hm1

  cout << endl;
  cout << " Updating NRGMats for the 1st time... " << endl;
  cout << endl;

//   UpdateMatrices(pEmptySite,&AeigCut, 
// 		 &AbasisH0,NRGMats, 4);

  UpdateMatrices(pEmptySite,&AeigCut, 
		 &AbasisH0,&STLNRGMats[0],STLNRGMats.size());
  cout << " ... done " << endl;


  //STLNRGMats[4].PrintAllBlocks();

  if (STLNRGMats.size()>4){
    STLNRGMats[7].PrintAllBlocks();
    for (int imat=4;imat<STLNRGMats.size();imat++){
      STLNRGMats[imat].CalcMatEl=OneChQ_cd_MatEl; // from now on
      STLNRGMats[imat].SaveMatYN=true;
    }
  }


  // Suscep test
  vector<double> ParamsBetabar;
  double betabar=0.727;
  ParamsBetabar.clear();
  ParamsBetabar.push_back(betabar);
  
  double dSz0=CalcOpAvg(ParamsBetabar,&AeigCut,&STLNRGMats[2],false,0);
  double dSz20=CalcOpAvg(ParamsBetabar,&AeigCut,&STLNRGMats[3],false,0);
  
  cout << " N=0: Sz = " << dSz0 << " Sz2 = " << dSz20 
       << "    T_M chi = " << dSz20-dSz0*dSz0 << endl;



  // New basis N=0
  // BuildBasis params
  vector <int> CommonQNs;
  vector<int> totSpos;

  CommonQNs.push_back(1); // No of common QNs
  CommonQNs.push_back(0); // pos of QN 1 in old
  CommonQNs.push_back(0); // pos of QN 1 in SingleSite
  // No total S variables. Leave totSpos empty
 
  BuildBasis(CommonQNs, totSpos,&AeigCut,&AbasisH0, 
	     pSingleSite,1);

  // Set up H_(N=0)

  H0.CalcHNMatEl=OneChQ_H0SMM_MatEl;
  H0.DiagHN(Params,&AbasisH0,pSingleSite,fHm1Array,pAeig);


  // De-allocate MatArray 
  delete[] fHm1Array;

  pAeig->PrintEn();


  // Set-up Matrix elements fN0_up, fN0_dn

  AeigCut.ClearAll();
  AeigCut=CutStates(pAeig,1000);

//   UpdateMatrices(pSingleSite,&AeigCut, 
// 		 &AbasisH0,NRGMats, 4);
  UpdateMatrices(pSingleSite,&AeigCut, 
		 &AbasisH0,&STLNRGMats[0],STLNRGMats.size());

//   cout << " Got here !!" << endl;
//   STLNRGMats[2].PrintAllBlocks();

}
///////////////////////////
///////////////////////////
///////////////////////////
///////////////////////////
///////////////////////////


///////////////////////////
///                     ///
///  2ch QSP CMphonons  ///
///                     ///
///                     ///
///////////////////////////


void TwoChQSP_SetH0CMphonon(vector<double> Params, 
			    CNRGarray* pAeig,
			    CNRGbasisarray* pSingleSite,
			    CNRGmatrix* NRGMats){

// void TwoChQSP_SetH0CMphonon(vector<double> Params, 
// 			    CNRGarray* pAeig,
// 			    CNRGbasisarray* pSingleSite,
// 			    CNRGmatrix* NRGMats){



  // Parameters:

  // Get Params
  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];

  // Watch out for normalization factors

  double U_norm=0.5*Params[2]/(sqrt(Lambda)*HalfLambdaFactor); 
  double ed_norm=Params[3]/(sqrt(Lambda)*HalfLambdaFactor);

  double chi_N[2]={Params[4],Params[7]*Params[4]};

  double w0=Params[5];
  double lambda=Params[6];
  double alpha=Params[7];

  int Nph=(int)Params[8];

  // Gamma_tilde = (2 Gamma/pi)/(1/2*(1+Lambda^-1))^2
  // Params[2] = sqrt(Gamma_tilde/Lambda) = chi_S
  // Change here: chi_S = Params[2], chi_A = alpha*Params[2]
  //double chi_N[2]={Params[2],Params[3]};

  double auxEl=0.0;

  cout << "e1 = " << U_norm << endl;
  cout << "e2 = " << ed_norm+U_norm << endl;
  cout << "e3 = " << 2.0*ed_norm+3.0*U_norm << endl;


  cout << "chi_N(S,A) = " << chi_N[0] << "  " << chi_N[1] << endl;

  cout << " w0 = " << w0 << endl;
  cout << " lambda = " << lambda << endl;
  cout << " alpha = " << alpha << endl;
  cout << " Nph = " << Nph << endl;


  // Set H0

  pAeig->ClearAll();

  CNRGarray AeigHimp(3);

  // Set basis for Hm1
  AeigHimp.Nshell=-1;
  AeigHimp.NQNumbers=3; //(Q,S,P)

  // new thing
  AeigHimp.totalS=true;
  AeigHimp.Sqnumbers.push_back(1);


  // |0>  = |-1 0 P> 
  AeigHimp.QNumbers.push_back(-1.0);
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(1.0);
  AeigHimp.dEn.push_back(U_norm);
  AeigHimp.BlockBegEnd.push_back(0);
  AeigHimp.BlockBegEnd.push_back(0);
  // |up> = |0 0.5 P> Simple as that (or |dn>)
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(0.5);
  AeigHimp.QNumbers.push_back(1.0);
  AeigHimp.BlockBegEnd.push_back(1);
  AeigHimp.BlockBegEnd.push_back(1);
  AeigHimp.dEn.push_back(ed_norm+U_norm);


  // |up dn> = |1 0 P> 
  AeigHimp.QNumbers.push_back(1.0);
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(1.0);
  AeigHimp.BlockBegEnd.push_back(2);
  AeigHimp.BlockBegEnd.push_back(2);
  AeigHimp.dEn.push_back(2.0*ed_norm+3.0*U_norm);
  //

  
  ////
  ////  Set the matrices Hm1 and <||d||>
  ////

  CNRGmatrix MatArrayH0[2];

  for (int imat=0; imat<2; imat++)
    MatArrayH0[imat].SyncNRGarray(AeigHimp);

  //
  // MatArrayH0[0]: Hm1 in the OLD basis
  //

  MatArrayH0[0].MatEl.push_back(U_norm);
  MatArrayH0[0].MatEl.push_back(ed_norm+U_norm);
  MatArrayH0[0].MatEl.push_back(2.0*ed_norm+3.0*U_norm);

  // one state per block, diagonal
  for (int ii=0;ii<=2;ii++)
    {
      MatArrayH0[0].MatBlockMap.push_back(ii);
      MatArrayH0[0].MatBlockMap.push_back(ii);
      MatArrayH0[0].MatBlockBegEnd.push_back(ii);
      MatArrayH0[0].MatBlockBegEnd.push_back(ii);
    }



  // MatArrayH0[1]: reduced <||d||> elements in the OLD basis

  MatArrayH0[1].MatEl.push_back(1.0);
  MatArrayH0[1].MatBlockMap.push_back(0);
  MatArrayH0[1].MatBlockMap.push_back(1);
  MatArrayH0[1].MatBlockBegEnd.push_back(0);
  MatArrayH0[1].MatBlockBegEnd.push_back(0);

  MatArrayH0[1].MatEl.push_back(-sqrt(2.0));
  MatArrayH0[1].MatBlockMap.push_back(1);
  MatArrayH0[1].MatBlockMap.push_back(2);
  MatArrayH0[1].MatBlockBegEnd.push_back(1);
  MatArrayH0[1].MatBlockBegEnd.push_back(1);

  //
  // Add a site: N=0 basis 
  //



  //CNRGbasisarray ACut=CutStates(&AeigHimp, 100);
  CNRGbasisarray ACut;
  CutStates(AeigHimp, ACut, 100);
  ACut.PrintAll();

  // build basis

  // Common QNs

  vector<int> CommonQNs;
  vector<int> totSpos;

  CommonQNs.push_back(3); // NoCommon Q,S and P

  CommonQNs.push_back(0);
  CommonQNs.push_back(1); // Pos in old
  CommonQNs.push_back(2);

  CommonQNs.push_back(0);
  CommonQNs.push_back(1); // Pos in site
  CommonQNs.push_back(3);

  CommonQNs.push_back(2); // Pos of parity in old

  totSpos.push_back(1); // Pos of S in old



  CNRGbasisarray AbasisH0;


  pSingleSite->PrintAll();



  BuildBasis(CommonQNs,totSpos,&ACut,&AbasisH0,pSingleSite,1);


  //AbasisH0.PrintQNumbers();


  // Add Phonons: simply copy blocks

  for (int ibl=0;ibl<AbasisH0.NumBlocks();ibl++)
    AbasisH0.CopyBlock(ibl,Nph);



  // Set parity

  vector<double> dAux;

  for (int ist=0;ist<AbasisH0.Nstates();ist++)
//      ACut.iDegen[ist]=(int)pow(-1.0,ACut.iDegen[ist]);
    dAux.push_back(pow(-1.0,AbasisH0.iDegen[ist])*AbasisH0.GetQNumberFromSt(ist,2));
 
  AbasisH0.SetLastQNumber(dAux);

  //AbasisH0.PrintAll();

  // Add vectors to block basis
  // 3xNph states



  CNRGmatrix H0(AbasisH0);

  H0.NeedOld=false;
  H0.UpperTriangular=true;
  H0.CheckForMatEl=Diag_check;
  H0.CalcHNMatEl=TwoChQSP_H0ph_MatEl;


  vector<double> ParamsH0;
  ParamsH0.push_back(chi_N[0]);
  ParamsH0.push_back(chi_N[1]);
  ParamsH0.push_back(w0);
  ParamsH0.push_back(lambda);
  ParamsH0.push_back(alpha);

  int ist=2;
  int jst=4;
  double testEl=H0.CalcHNMatEl(ParamsH0,&AbasisH0,pSingleSite,MatArrayH0,ist,jst);


  cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 

//    for (ist=8;ist<=19;ist++)
//      {
//        for (jst=8;jst<=19;jst++)
// 	 {
//     testEl=H0.CalcHNMatEl(ParamsH0,pAbasisH0,pSingleSite,MatArrayH0,ist,jst);
//     cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 
// 	 }
//      }


  H0.DiagHN(ParamsH0,&AbasisH0,pSingleSite,MatArrayH0,pAeig);
  
  pAeig->PrintEn();


  // Set-up Matrix elements fN_ch1, fN_ch2

  ACut.ClearAll();
  ACut=CutStates(pAeig,2000);

  UpdateMatrices(pSingleSite,&ACut, 
		 &AbasisH0,NRGMats, 2);



}

////////////////////////////////////////////

///////////////////////////
///                     ///
///  1chQS DQD H0       ///
///                     ///
///  | Q S > basis      ///
///                     ///
///////////////////////////

void OneChQS_SetH0_DQD(vector<double> Params, 
		       CNRGarray* pAeig,
		       CNRGbasisarray* pSingleSite,
		       vector<CNRGmatrix> &STLNRGMats){

  // Set initial CNRG array (N=0)
  // Actually, it will set H0 at the end of the day

  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];
  // Testing:
  //Lambda=1.0;
  //HalfLambdaFactor=1.0;


  double U1tilde=Params[2]/(sqrt(Lambda)*HalfLambdaFactor);
  double ed1tilde=Params[3]/(sqrt(Lambda)*HalfLambdaFactor);
  double gamma1tilde=Params[4];

  // Remember:
  // gamma~=(2*Gamma/pi)^1/2/(sqrt(Lambda)*HalfLambdaFactor)
  // Why? see Wilson's paper (E. 2.18 with N=-1)
  // Also: Multiplied U1~ and ed~ by sqrt(Lambda) as compared to previous codes
  // (i.e. divide by sqrt(Lambda) instead of Lambda)
  // Thus, no lmabda^{1/2} factor in setting up H0.

  double U2tilde=Params[5]/(sqrt(Lambda)*HalfLambdaFactor);
  double ed2tilde=Params[6]/(sqrt(Lambda)*HalfLambdaFactor);
  double gamma2tilde=Params[7];

  // The Basis: a two Channel-like basis in QS basis 
  // with Sz entering iDegen (Sz=10*iDegen)

  cout << " U1~= " << U1tilde
       << " ed1~= " << ed1tilde
       << " t1~= " << gamma1tilde
       << " U2~= " << U2tilde
       << " ed2~= " << ed2tilde
       << " t2~= " << gamma2tilde
       << endl;
    
  CNRGbasisarray AbasisHm1; // will be AbasisHm1

  // Set AbasisHm1
  //TwoChQS_SetSingleSite(&AbasisHm1);
  TwoDotQS_SetInitialSite(&AbasisHm1);

  AbasisHm1.Nshell=-1;
  // Set energies (diagonal)

  AbasisHm1.dEn.push_back(0.5*U1tilde+0.5*U2tilde); // | 0 0 >

  AbasisHm1.dEn.push_back(0.5*U1tilde+ed2tilde+0.5*U2tilde); // | 0 up>
  AbasisHm1.dEn.push_back(0.5*U2tilde+ed1tilde+0.5*U1tilde); // | up 0>
  //AbasisHm1.dEn.push_back(0.5*U1tilde+ed2tilde+0.5*U2tilde); // | 0 dn>
  //AbasisHm1.dEn.push_back(0.5*U2tilde+ed1tilde+0.5*U1tilde); // | dn 0>

  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+ed2tilde+0.5*U2tilde); //|up dn> - |dn up>/sqrt2 
  AbasisHm1.dEn.push_back(0.5*U1tilde+2.0*ed2tilde+1.5*U2tilde); //|0  up dn>  
  AbasisHm1.dEn.push_back(0.5*U2tilde+2.0*ed1tilde+1.5*U1tilde); //|up dn  0>  

  //AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+ed2tilde+0.5*U2tilde); //|dn dn> 
  //AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+ed2tilde+0.5*U2tilde); //|up dn> + |dn up>/sqrt2 
  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+ed2tilde+0.5*U2tilde); //|up up> 

  AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+2.0*ed2tilde+1.5*U2tilde);//|up  updn>
  AbasisHm1.dEn.push_back(ed2tilde+0.5*U2tilde+2.0*ed1tilde+1.5*U1tilde);//|updn  up>
  //AbasisHm1.dEn.push_back(ed1tilde+0.5*U1tilde+2.0*ed2tilde+1.5*U2tilde);//|dn  updn>
  //AbasisHm1.dEn.push_back(ed2tilde+0.5*U2tilde+2.0*ed1tilde+1.5*U1tilde);//|updn  dn>

  AbasisHm1.dEn.push_back(2.0*ed1tilde+1.5*U1tilde+2.0*ed2tilde+1.5*U2tilde); 
  //|updn  updn >  

  AbasisHm1.SetEigVecToOne();

  for (int imat=0;imat<STLNRGMats.size();imat++)
    STLNRGMats[imat].SyncNRGarray(AbasisHm1);

  // Set up cd1, cd2
  vector<CNRGmatrix> AuxMatArray; // Arrays needed to set up H0
  CNRGmatrix auxMat;

  auxMat=STLNRGMats[0]; // Same rules as fN
  auxMat.SetZeroMatrix();
  AuxMatArray.push_back(auxMat); // cd1
  AuxMatArray.push_back(auxMat); // cd2

  double auxEl[2]={0.0};


  for (int ist=0;ist<10;ist++){
    for (int jst=0;jst<10;jst++){
      auxEl[0]=TwoDotQS_fd_reduced_table(1,jst,ist); // cd1 (not dagger)
      auxEl[1]=TwoDotQS_fd_reduced_table(2,jst,ist); // cd2 (not dagger)
      //cout << " <" << ist << "|cd1|"<< jst <<"> = " << auxEl[0];
      //cout << "   <" << ist << "|cd2|"<< jst <<"> = " << auxEl[1] << endl;
      //cout << " Putting into cd1: " << endl;
      AuxMatArray[0].PushMatEl(auxEl[0],ist,jst); 
      //cout << " Putting into cd2: " << endl;
      AuxMatArray[1].PushMatEl(auxEl[1],ist,jst); 
    }
  }
  // end push Mat else


//   for (int imat=0;imat<AuxMatArray.size();imat++){
//      cout << " ---- cd" << imat+1 << "  ---- " << endl; 
//      AuxMatArray[imat].PrintAllBlocks();

//   }


  // Set up operators in the Hm1 basis
  // ndot, Sz, S1 dot S2 if calcdens==0

  switch(STLNRGMats.size()){
  case 1:
    cout << " Thermodynamics only." << endl;
    break;
  case 3:
    // cd1, cd2 if calcdens==1 (just use the ops defined above)
//     STLNRGMats[1]=AuxMatArray[0];
//     STLNRGMats[2]=AuxMatArray[1];
    STLNRGMats[1].CopyData(&AuxMatArray[0]);
    STLNRGMats[2].CopyData(&AuxMatArray[1]);
    cout << " DQD: spectral functions still at implementation stage..." << endl;
    //exit(0);
    break;
  case 5:
    // ndot    - imat=1 
    // ndot^2 - imat=2
    // Sdot^2   - imat=3
    // S1 dot S2 - imat=4
    // Setting all to Identity (all diagonal)
    for (int imat=1;imat<STLNRGMats.size();imat++){
      STLNRGMats[imat].SetDiagonalMatrix(0.0);
    }
    // Setting ndot, Sz, Sz2
    // Loop over diagonal
    for (int ibl=0;ibl<AbasisHm1.NumBlocks();ibl++){
      double ndot=AbasisHm1.GetQNumber(ibl,0)+2.0; // nel=Q+2
      double S=AbasisHm1.GetQNumber(ibl,1); // S is set on each state.
      for (int istbl=0;istbl<AbasisHm1.GetBlockSize(ibl);istbl++){
	STLNRGMats[1].PushBlockMatEl(ndot,ibl,ibl,istbl,istbl); // n
	STLNRGMats[2].PushBlockMatEl(ndot*ndot,ibl,ibl,istbl,istbl); //n^2
	STLNRGMats[3].PushBlockMatEl(S*(S+1.0),ibl,ibl,istbl,istbl); // S^2
      }
      // end loop over diag states
    }
    // end loop over blocks

    // S1 dot S2 - only 2 nonzero mat els
    STLNRGMats[4].PushMatEl(-0.75,3,3);
    STLNRGMats[4].PushMatEl(0.25,6,6);
    break;
  default:
    cout << " Can't figure out what to do with the STLNRGMats... Exiting." << endl;
    exit(0);
  }
  // end switch STLNRGMats.size()

  
//   for (int imat=0;imat<STLNRGMats.size();imat++){
//      cout << " ---- Matrix " << imat << "  ---- " << endl; 
//      STLNRGMats[imat].PrintAllBlocks();
//   }



  // Connect to first site of the chain: Create new basis and so on. 
  // Build basis for DQD + 1 site

  vector<int> CommonQNs;
  vector<int> totSpos;
  CNRGbasisarray AbasisH0;
  CommonQNs.push_back(2); // Q and S and commont QNs
  CommonQNs.push_back(0); // position of Q in old basis
  CommonQNs.push_back(1); // position of S in old basis
  CommonQNs.push_back(0); // position of Q in SingleSite
  CommonQNs.push_back(1); // position of S in SingleSite

  totSpos.push_back(1);   // SU(2) symmetry in position 1

  BuildBasis(CommonQNs,totSpos,&AbasisHm1,&AbasisH0,pSingleSite,0);


  // Diagonalize H0: Aeig will be the new vector

  auxMat.ClearAll(); //  H0 (re-using auxMat)
  auxMat.NeedOld=false;
  auxMat.UpperTriangular=true;
  auxMat.CheckForMatEl=Diag_check;
  auxMat.CalcHNMatEl=OneChQS_H0DQD_MatEl;

  vector<double> ParamsH0;
  ParamsH0.push_back(gamma1tilde);
  ParamsH0.push_back(gamma2tilde);

//   int ist=5;
//   int jst=9;
//   double testEl=auxMat.CalcHNMatEl(ParamsH0,&AbasisH0,pSingleSite,&AuxMatArray[0],ist,jst);

//   cout << "<" << ist <<"|H0|"<<jst<<">= " << testEl 
//        << ";  Gamma1 = " << gamma1tilde 
//        << ";  Gamma2 = " << gamma2tilde 
//        << endl;

  auxMat.DiagHN(ParamsH0,&AbasisH0,pSingleSite,&AuxMatArray[0],pAeig);
  
  pAeig->PrintEn();


  // Update all operators return all.

  AbasisHm1.ClearAll();
  AbasisHm1=CutStates(pAeig,2000);
  
  UpdateMatrices(pSingleSite,&AbasisHm1, 
		 &AbasisH0,&STLNRGMats[0],STLNRGMats.size());
//   for (int imat=0;imat<STLNRGMats.size();imat++){
//      cout << " ---- Matrix " << imat << "  ---- " << endl; 
//      STLNRGMats[imat].PrintAllBlocks();
//   }

  cout << " ... H0 done. " << endl;


}

/////////////////


//////////////////////////////////
///                            ///
///  1chQS Anderson model Hm1  ///
///                            ///
///  | Q S > basis             ///
///                            ///
//////////////////////////////////

void OneChQS_SetAnderson_Hm1(vector<double> Params, 
			     CNRGarray* pAeig,
			     vector<CNRGmatrix> &STLNRGMats){

  //			  CNRGbasisarray* pSingleSite,


  // Set initial CNRG array (N=-1)

  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];
  // Testing:
  //Lambda=1.0;
  //HalfLambdaFactor=1.0;

  double U=Params[2];
  double ed=Params[3];


//   double U1tilde=Params[2]/(sqrt(Lambda)*HalfLambdaFactor);
//   double ed1tilde=Params[3]/(sqrt(Lambda)*HalfLambdaFactor);
//   double gamma1tilde=Params[4];

  pAeig->NQNumbers=2;
  pAeig->Nshell=-1;
  // new
  pAeig->totalS=true;
  pAeig->Sqnumbers.push_back(1);
  

  // |0> = |Q=-1, S=0>
  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(0.0);
  //pAeig->iDegen.push_back(0); // S=0  // Do I need iDegen???

  // |up> = |Q=0, S=1/2>
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.5);
  //pAeig->iDegen.push_back(0); // S=1/2


  // |up dn> = |Q=1, Sz=0>
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(0.0);
  //pAeig->iDegen.push_back(0); // S=0

  // One state per block
  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor)); // up
  pAeig->dEn.push_back((2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor));


  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);

  pAeig->PrintEn();

  // Set Matrix elements (reduced!)

  // fN (check -sqrt(2) and whatnot...

  //STLNRGMats[0].SyncNRGarray(*pAeig);
  // Sync all
  for (int imat=0;imat<STLNRGMats.size();imat++)
    STLNRGMats[imat].SyncNRGarray(*pAeig);

  STLNRGMats[0].UpperTriangular=false;

  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);

  STLNRGMats[0].MatEl.push_back(-sqrt(2.0));
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockMap.push_back(2);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);

  // This should work
  // cd_red
//   if (STLNRGMats.size()>1){
//       cout << " Setting up cd_red and cd_dn... ";
//       //STLNRGMats[1]=STLNRGMats[0];
//       STLNRGMats[1].CopyData(&STLNRGMats[0]);
//       cout << " ... done." << endl;
//   }
  // end if STLNRGMats.size()>1

  // Changing the way things are being done (June 2013)
  if (strcmp(STLNRGMats[1].MatName,"Ndot")==0){

    cout << " Setting up n_d, nd^2 and S^2... " << endl;
    // ndot    - imat=1 
    // ndot^2 - imat=2
    // Sdot^2   - imat=3
    // Setting all to Identity (all diagonal)
    for (int imat=1;imat<STLNRGMats.size();imat++){
      STLNRGMats[imat].SetDiagonalMatrix(0.0);
    }
    // Setting ndot, Sz, Sz2
    // Loop over diagonal
    for (int ibl=0;ibl<pAeig->NumBlocks();ibl++){
      double ndot=pAeig->GetQNumber(ibl,0)+1.0; // nel=Q+1
      double S=pAeig->GetQNumber(ibl,1); // S is set on each state.
      for (int istbl=0;istbl<pAeig->GetBlockSize(ibl);istbl++){
	STLNRGMats[1].PushBlockMatEl(ndot,ibl,ibl,istbl,istbl); // n
	STLNRGMats[2].PushBlockMatEl(ndot*ndot,ibl,ibl,istbl,istbl); //n^2
	STLNRGMats[3].PushBlockMatEl(S*(S+1.0),ibl,ibl,istbl,istbl); // S^2
      }
      // end loop over diag states
    }
    // end loop over blocks
    // Test
    //cout << " <0|nd|0> = " << STLNRGMats[1].GetMatEl(0,0) << endl;
  }else if (strcmp(STLNRGMats[1].MatName,"Szomega")==0){
    cout << " Dynamical spin susceptibiliy: setting up <Sz>_N=-1  ";
    // Already Sync'ed above... 
    int i1=0;
    for (int ibl=0; ibl<pAeig->NumBlocks(); ibl++){
      double Qi=pAeig->GetQNumber(ibl,0);
      double Si=pAeig->GetQNumber(ibl,1);

      cout << " Q = " << Qi
	   << " S = " << Si 
	   << endl;

      // Sz dynamical and Nd dynamical
      for (int imat=1;imat<=2; imat++){
	// Block diagonal
	STLNRGMats[imat].MatBlockMap.push_back(ibl);
	STLNRGMats[imat].MatBlockMap.push_back(ibl);
	// Position in MatBlock follows ist
	// Only one state per block !
	STLNRGMats[imat].MatBlockBegEnd.push_back(i1);
      	STLNRGMats[imat].MatBlockBegEnd.push_back(i1);
      }
      STLNRGMats[1].MatEl.push_back(Si*sqrt(3.0)); // <Sz> reduced!!!
      STLNRGMats[2].MatEl.push_back(Qi+1.0); // <Nd> 
      i1++;
    }
    // end loop in blocks (Set-up of Nd,Sz,Sz2)

    //Setting Sz2 and Nd2
    for (int imat=3; imat<=4; imat++){
      STLNRGMats[imat].CopyData(&STLNRGMats[imat-2]);
      for (int ii=0;ii<STLNRGMats[imat].MatEl.size();ii++)
	STLNRGMats[imat].MatEl[ii]*=STLNRGMats[imat].MatEl[ii];
    }
    // end setting Sz2 and Nd2

    STLNRGMats[1].PrintAllBlocks();
    STLNRGMats[2].PrintAllBlocks();

    cout << " ... done." << endl;
  } else if (strcmp(STLNRGMats[1].MatName,"cdot1")==0){
    cout << " Spectral function: setting up cd_red and cd_dn... ";
    STLNRGMats[1].CopyData(&STLNRGMats[0]);
    cout << " ... done." << endl;
  }
  // end if Mat[1]=Ndot


  cout << " ... H_{-1} done. " << endl;


}
// end set H_-1 for Anderson with Q S symmetry

/////////////////

//////////////////////////////////
///                            ///
///  1chQS Kond model H0       ///
///                            ///
///  | Q S > basis             ///
///                            ///
//////////////////////////////////

void OneChQS_SetKondoH0(vector<double> Params, 
			CNRGarray* pAeig,
			CNRGbasisarray* pSingleSite,
			vector<CNRGmatrix> &STLNRGMats){

  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];
  // Testing:
  //Lambda=1.0;
  //HalfLambdaFactor=1.0;

  // is Given by chi_m1! Should get the pseudogap case right.
  // Calculate Jtilde??? I guess I should...
  // Ok, here's the deal:
  // - rho_0 J_0 enters as "Gamma" in the code (or is calculated as 0.5*Fsq*(rho_0 J_0)
  //
  // - Params[2](="chi_m1") is sqrt(Fsq*rho_0*J_0/pi)/(sqrt(Lambda)*HalfLambdaFactor); 
  // 
  // - Jtilde is Fsq*rho_0*J_0/(sqrt(Lambda)*HalfLambdaFactor)
  //
  // Therefore:

  double Jtilde=Params[2]*Params[2]*pi*(sqrt(Lambda)*HalfLambdaFactor);

  cout << "Jtilde  = " << Jtilde << endl;

  double auxEn=0.0;

  pAeig->ClearAll();
  STLNRGMats[0].ClearAll();

  pAeig->NQNumbers=2;
  pAeig->Nshell=0;
  // new
  pAeig->totalS=true;
  pAeig->Sqnumbers.push_back(1);


  // QNumbers
  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(0.5);

  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.0);

  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(1.0);

  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(0.5);


  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  auxEn=Jtilde*Sdots_totalS(0.5, 0.5, 0.5, 0.0);
  pAeig->dEn.push_back(auxEn);
  auxEn=Jtilde*Sdots_totalS(0.0, 0.0, 0.5, 0.5);
  pAeig->dEn.push_back(auxEn);
  auxEn=Jtilde*Sdots_totalS(1.0, 1.0, 0.5, 0.5);
  pAeig->dEn.push_back(auxEn);
  auxEn=Jtilde*Sdots_totalS(0.5, 0.5, 0.5, 0.0);
  pAeig->dEn.push_back(auxEn);


  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);


  // Set REDUCED Matrix elements <Q S||f||Q' S'> = <Q' S'||f+||Q S> 


  // Test Kondo model
  vector <int> Indexes;
  double MatEl=0.0;

  STLNRGMats[0].SyncNRGarray(*pAeig);
 
  STLNRGMats[0].MatBlockMap.push_back(0); // <-1 0.5||f0||0 0>=<0 0||f+0||-1 0.5>
  STLNRGMats[0].MatBlockMap.push_back(1); // |0 0>=|0 0.5>x|0 0.5 -0.5> (stsite=2)
                                     // |-1 0.5>=|0 0.5>x|-1 0 0> (stsite=0)

  Params.clear();
  Indexes.clear(); 
  Params.push_back(0.0);   // STp
  Params.push_back(0.5); // Soldp
  Params.push_back(0.5);    // ST
  Params.push_back(0.5);  // Sold

  Indexes.push_back(2); // sitest'
  Indexes.push_back(0); // sitest
  Indexes.push_back(1); // no channels

  MatEl=QfdQm1_totS(Params,Indexes, pSingleSite);
  STLNRGMats[0].MatEl.push_back(MatEl);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);

  STLNRGMats[0].MatBlockMap.push_back(0); // <-1 0.5||f0||0 1>=<0 1||f+0||-1 0.5>
  STLNRGMats[0].MatBlockMap.push_back(2); // |0 1>=|0 0.5>x|0 0.5 0.5> (stsite=1)
                                     // |-1 0.5>=|0 0.5>x|-1 0> (stsite=0)
  Params.clear();
  Indexes.clear(); 
  Params.push_back(1.0);   // STp
  Params.push_back(0.5); // Soldp
  Params.push_back(0.5);    // ST
  Params.push_back(0.5);  // Sold

  Indexes.push_back(1); // sitest'
  Indexes.push_back(0); // sitest
  Indexes.push_back(1); // no channels

  MatEl=QfdQm1_totS(Params,Indexes, pSingleSite);
  STLNRGMats[0].MatEl.push_back(MatEl);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);

  STLNRGMats[0].MatBlockMap.push_back(1); // <0 0||f0||1 0.5>=<1 0.5||f+0||0 0>
  STLNRGMats[0].MatBlockMap.push_back(3); // |1 0.5>=|0 0.5>x|1 0 0> (stsite=3)
                                     // |0 0>=|0 0.5>x|0 0.5 -0.5> (stsite=2)
  
  Params.clear();
  Indexes.clear(); 
  Params.push_back(0.5);   // STp
  Params.push_back(0.5); // Soldp
  Params.push_back(0.0);    // ST
  Params.push_back(0.5);  // Sold

  Indexes.push_back(3); // sitest'
  Indexes.push_back(2); // sitest
  Indexes.push_back(1); // no channels

  MatEl=QfdQm1_totS(Params,Indexes, pSingleSite);
  STLNRGMats[0].MatEl.push_back(MatEl);
  STLNRGMats[0].MatBlockBegEnd.push_back(2);
  STLNRGMats[0].MatBlockBegEnd.push_back(2);


  STLNRGMats[0].MatBlockMap.push_back(2); // <0 1||f0||1 0.5>=<1 0.5||f+0||0 1>
  STLNRGMats[0].MatBlockMap.push_back(3); // |1 0.5>=|0 0.5>x|1 0 0> (stsite=3)
                                     // |0 1>=|0 0.5>x|0 0.5 0.5> (stsite=1)
                                     
  Params.clear();
  Indexes.clear(); 
  Params.push_back(0.5);   // STp
  Params.push_back(0.5); // Soldp
  Params.push_back(1.0);    // ST
  Params.push_back(0.5);  // Sold

  Indexes.push_back(3); // sitest'
  Indexes.push_back(1); // sitest
  Indexes.push_back(1); // no channels

  MatEl=QfdQm1_totS(Params,Indexes, pSingleSite);
  STLNRGMats[0].MatEl.push_back(MatEl);
  STLNRGMats[0].MatBlockBegEnd.push_back(3);
  STLNRGMats[0].MatBlockBegEnd.push_back(3);



//   for (int imbl=0;imbl<STLNRGMats[0].NumMatBlocks();imbl++)
//     {
//       STLNRGMats[0].PrintMatBlock(imbl);
//     }




}
// end set H_0 for Kondo with Q S symmetry



///////////////////////////
///                     ///
///    2ch Kondo        ///
///                     ///
///  | Q S > basis      ///
///                     ///
///////////////////////////


///////////////////////////////////////
////  Two-channel Kondo model   ///////
///////////////////////////////////////

void TwoChQS_SetKondoH0(vector<double> Params,
			CNRGbasisarray* pSingleSite,
			CNRGarray* pAeig,
			vector<CNRGmatrix> &STLNRGMats){
  //
  // Set the 2-channel Kondo hamiltonian in the Q,S basis.
  // 

  // Do I need to return AbasisH0 ?? I don't think so...
  //			CNRGbasisarray* pAbasisH0,


  CNRGbasisarray AbasisH0;

  double J1=Params[0];
  double J2=Params[1];

  double auxEl=0.0;


  cout << "J1 = " << J1 << endl;
  cout << "J2 = " << J2 << endl;

  ////////////////////////
  // 1 - Diagonalize impurity Hamiltonian: H_(-1)
  ////////////////////////

  CNRGarray AeigHimp(2);

  AeigHimp.Nshell=-1;
  AeigHimp.NQNumbers=2;

  // Only 1 state

  // |up> = |0 0.5> Simple as that (or |dn>)
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(0.5);

  AeigHimp.dEn.push_back(0.0);

  AeigHimp.BlockBegEnd.push_back(0);AeigHimp.BlockBegEnd.push_back(0);


  ////////////////////////
  // 2 - Add one site: Build basis
  ////////////////////////


  //CNRGbasisarray AbasisH0(3);

  AbasisH0.NQNumbers=2;

  CNRGbasisarray ACut=CutStates(&AeigHimp, 100);
  QS_BuildBasis(&ACut,&AbasisH0,pSingleSite,1);

  cout << "No blocks = " << AbasisH0.NumBlocks() << endl;
  cout << "No states = " << AbasisH0.Nstates() << endl;

  AbasisH0.PrintAll();


  ////////////////////////
  // 3 - Diagonalize Himp + Hcoupling.
  ////////////////////////


  //CNRGmatrix H0(*pAbasisH0);
  CNRGmatrix H0(AbasisH0);

  H0.NeedOld=false;
  H0.UpperTriangular=false;
  H0.CheckForMatEl=Diag_check;
  H0.CalcHNMatEl=TwoChQS_H0Kondo_MatEl; // Does not actually use STLNRGMats


  vector<double> ParamsH0;
  ParamsH0.push_back(J1);
  ParamsH0.push_back(J2);



   int ist=8;
   int jst=8;
   double testEl=H0.CalcHNMatEl(ParamsH0,&AbasisH0,pSingleSite,&STLNRGMats[0],ist,jst);

   cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 

   // Actually, H0 is diagonal in this basis... what the hell, let's do it.

   H0.DiagHN(ParamsH0,&AbasisH0,pSingleSite,&STLNRGMats[0],pAeig); // Does not actually use STLNRGMats



 // Set REDUCED Matrix elements <Q S||f1||Q' S'> = <Q' S'||f1+||Q S> 
//                              <Q S||f2||Q' S'> = <Q' S'||f2+||Q S> 



  // 
  // Ok Need to write the following according to the TwoChannel case.
  // Better to write a loop.
  // 
  // There are 2 x 16 = 32 states in AbasisH0, pAeig.

  STLNRGMats[0].SyncNRGarray(AbasisH0);
  STLNRGMats[1].SyncNRGarray(AbasisH0);
 
  // The following is STILL IN THE ONE-CHANNEL CASE!!
  // Need a loop... will do tomorrow.
  // 
  // for imat=0,1
  // Loop: states in pAbasisH0

  // check if nonzero with STLMatArray[0].CheckForMatEl (=TwoChQS_cd_check)
  // If YES, calculate with STLMatArray[0].CalcMatEl (=TwoChQS_fNch1(2)_MatEl;)
  // Actually, I made this much simpler now.

  for (int imat=0; imat<=1;imat++)
    STLNRGMats[imat].SetMatrix(&AbasisH0,pSingleSite);
  // end loop in imat

  // Ok, the eigenvectors of Aeig are NOT in a diagonal form.
  // Let's rotate the STLmatrices in the Aeig basis.

  // Using AbasisH0 for that...

  cout << " f0_ch1 (in BASIS! Need to \"rotate\") : " << endl;
  //STLNRGMats[0].PrintAllBlocks();
  STLNRGMats[0].PrintMatBlock(1,3);

  CNRGbasisarray AeigCut=CutStates(pAeig,2000);
  
  UpdateMatrices(pSingleSite,&AeigCut, 
		 &AbasisH0,&STLNRGMats[0],STLNRGMats.size());

  cout << " f0_ch1 : " << endl;
  //STLNRGMats[0].PrintAllBlocks();
  STLNRGMats[0].PrintMatBlock(1,3);


}
// end setH0Kondo


///////////////////////////////
///                         ///
///  1ch | Nup Pdn > basis  ///
///    (Oct 2013)           ///
///                         ///
///////////////////////////////



void OneChNupPdn_SetH0_AndersonMajorana(vector<double> Params,
					CNRGbasisarray* pSingleSite,
					CNRGarray* pAeig, 
					vector<CNRGmatrix> &STLNRGMats){

  // Set initial CNRG array (N=-1)


  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];

  double Utilde=Params[2]/(Lambda*HalfLambdaFactor);
  double edtilde=Params[3]/(Lambda*HalfLambdaFactor);
  double gammatilde=Params[4];
  double hz=0.5*Params[5]/(Lambda*HalfLambdaFactor);; // hz=g mu_B Bz: Zeeman term is:  - hz S_z 

  // Remember:
  // gamma~=(2*Gamma/pi)^1/2/(sqrt(Lambda)*HalfLambdaFactor)
  // Why? see Wilson's paper (E. 2.18 with N=-1)

  // Note that "Lambda" (and not sqrt(Lambda)) appears in the
  // denominator of U~ and ed~.
  // This is different than the DQD code above since we are
  // generating only Hm1 here: H0 will be set up later by the NRG engine
  // and there will be a \sqrt(Lambda) multiplicating factor later on. 

  // Add other parameters here! t1, t2, phi, etc.
  
  double t1=Params[6]/(Lambda*HalfLambdaFactor);
  double t2=Params[7]/(Lambda*HalfLambdaFactor);
  double phi_mag=Params[8];
  double em=Params[9]/(Lambda*HalfLambdaFactor);

  
  cout  << " Lambda= " << Lambda
	<< " HalfLambdaFactor= " << HalfLambdaFactor
	<< endl
	<< " U~= " << Utilde
	<< " ed~= " << edtilde
	<< " t~= " << gammatilde
	<< " hz~= " << hz << endl
	<< " t1~= " << t1
	<< " t2~= " << t2
	<< " phimag= " << phi_mag
	<< " em~= " << em
	<< endl;


  // Build H_m1 (N=-1) 

  CNRGbasisarray AbasisHm1(2);

  // Set basis for Hm1
  AbasisHm1.Nshell=-1;
  AbasisHm1.NQNumbers=2; //(Nup,Pdn=(-1)^Ndn)
  AbasisHm1.totalS=false;

  // |0>d |dn>f , |dn>d |0>f : |Nup=0 Pdn=-1>
  AbasisHm1.QNumbers.push_back(0.0);
  AbasisHm1.QNumbers.push_back(-1.0);
  AbasisHm1.dEn.push_back(0.5*Utilde+em);
  AbasisHm1.dEn.push_back(edtilde+0.5*Utilde-hz-em);
  AbasisHm1.BlockBegEnd.push_back(0);AbasisHm1.BlockBegEnd.push_back(1);

  // |0>d |0>f , |dn>d |dn>f : |Nup=0 Pdn=+1>
  AbasisHm1.QNumbers.push_back(0.0);
  AbasisHm1.QNumbers.push_back(1.0);
  AbasisHm1.dEn.push_back(0.5*Utilde-em);
  AbasisHm1.dEn.push_back(edtilde+0.5*Utilde-hz+em);
  AbasisHm1.BlockBegEnd.push_back(2);AbasisHm1.BlockBegEnd.push_back(3);

  // |up>d |dn>f , |up dn>d |0>f : |Nup=1 Pdn=-1>
  AbasisHm1.QNumbers.push_back(1.0);
  AbasisHm1.QNumbers.push_back(-1.0);
  AbasisHm1.dEn.push_back(edtilde+0.5*Utilde+hz+em);
  AbasisHm1.dEn.push_back(2.0*edtilde+1.5*Utilde-em);
  AbasisHm1.BlockBegEnd.push_back(4);AbasisHm1.BlockBegEnd.push_back(5);

  // |up>d |0>f , |up dn>d |dn>f : |Nup=1 Pdn=+1>
  AbasisHm1.QNumbers.push_back(1.0);
  AbasisHm1.QNumbers.push_back(1.0);
  AbasisHm1.dEn.push_back(edtilde+0.5*Utilde+hz-em);
  AbasisHm1.dEn.push_back(2.0*edtilde+1.5*Utilde+em);
  AbasisHm1.BlockBegEnd.push_back(6);AbasisHm1.BlockBegEnd.push_back(7);

  AbasisHm1.SetEigVecToOne();

  // Set up cd and fd: easier in this basis. 
  // Will Rotate later.
  for (int imat=0;imat<STLNRGMats.size();imat++)
    STLNRGMats[imat].SyncNRGarray(AbasisHm1);


  //AbasisHm1.PrintAll();


  vector<CNRGmatrix> AuxMatArray; // Arrays needed to set up H0 
                       //(cd_up, cd_dn and f_Maj_dn)
  CNRGmatrix auxMat;

  auxMat=STLNRGMats[0]; // Same rules as fN_up
  auxMat.SetZeroMatrix(); // Build block structure
  AuxMatArray.push_back(auxMat); // cd_up
  auxMat=STLNRGMats[1]; // Same rules as fN_dn
  auxMat.SetZeroMatrix(); // Build block structure
  AuxMatArray.push_back(auxMat); // cd_dn

  // c_up|up>|0(dn)>=+|0>|0(dn)>  
  // c_up|up dn>|0(dn)>=+|dn>|0(dn)>  
  
  AuxMatArray[0].PushMatEl(1.0,2,6); 
  AuxMatArray[0].PushMatEl(1.0,0,4);
  AuxMatArray[0].PushMatEl(1.0,1,5); 
  AuxMatArray[0].PushMatEl(1.0,3,7); 
   
  // c_dn|dn>|0(dn)>=+|0>|0(dn)>  
  // c_dn|up dn>|0(dn)>=-|up>|0(dn)>  
  AuxMatArray[1].PushMatEl(1.0,2,1); 
  AuxMatArray[1].PushMatEl(1.0,0,3); 
  AuxMatArray[1].PushMatEl(-1.0,6,5); 
  AuxMatArray[1].PushMatEl(-1.0,4,7); 

 
  if (	(strcmp(STLNRGMats[2].MatName,"cd_up")==0)&&
	(strcmp(STLNRGMats[3].MatName,"cd_dn")==0) ){

    AuxMatArray.push_back(auxMat); // [2] f_Maj_dn

    // cd1, cd2 if calcdens==1 (just use the ops defined above)
    // Define c_up and c_dn and f_Maj_dn matrices in this basis!
    
    cout << " Majorana: spectral functions still at implementation stage..." << endl;

    // f_dn|up (dn)>|dn>=-|up (dn)>|0>  
    // f_dn|up dn (0) >|dn>=+|up dn (0)>|0>  
    AuxMatArray[2].PushMatEl(1.0,2,0); 
    AuxMatArray[2].PushMatEl(-1.0,1,3); 
    AuxMatArray[2].PushMatEl(1.0,5,7); 
    AuxMatArray[2].PushMatEl(-1.0,6,4); 
    
    //   cout << " f_dn: " << endl;
    //   AuxMatArray[2].PrintAllBlocks();

  }else if ( (strcmp(STLNRGMats[2].MatName,"Ndot")==0)&&
	     (strcmp(STLNRGMats[3].MatName,"Szdot")==0) ){
    //////////////////
    
    cout << " Setting up 1-body observables..." << endl;

    auxMat=STLNRGMats[2]; // Same rules as Ndot
    auxMat.SetZeroMatrix(); // Build block structure
  
    AuxMatArray.push_back(auxMat); // [2] Ndot
    AuxMatArray.push_back(auxMat); // [3] Szdot
    AuxMatArray.push_back(auxMat); // [4] NMaj
  
//     for (int imat=2; imat<=4; imat++){
//       AuxMatArray[imat]=STLNRGMats[imat]; // Same rules as ndot, szdot, nmaj
//       AuxMatArray[imat].SetZeroMatrix(); // Build block structure
//     }
 
    // ndot=nup+ndown
    AuxMatArray[2].PushMatEl(0.0,0,0); 
    AuxMatArray[2].PushMatEl(1.0,1,1);
    AuxMatArray[2].PushMatEl(0.0,2,2); 
    AuxMatArray[2].PushMatEl(1.0,3,3); 
    AuxMatArray[2].PushMatEl(1.0,4,4); 
    AuxMatArray[2].PushMatEl(2.0,5,5);
    AuxMatArray[2].PushMatEl(1.0,6,6); 
    AuxMatArray[2].PushMatEl(2.0,7,7); 
   
    // szdot
    AuxMatArray[3].PushMatEl(0.0,0,0); 
    AuxMatArray[3].PushMatEl(-0.5,1,1);
    AuxMatArray[3].PushMatEl(0.0,2,2); 
    AuxMatArray[3].PushMatEl(-0.5,3,3); 
    AuxMatArray[3].PushMatEl(0.5,4,4); 
    AuxMatArray[3].PushMatEl(0.0,5,5);
    AuxMatArray[3].PushMatEl(0.5,6,6); 
    AuxMatArray[3].PushMatEl(0.0,7,7); 

    // nMaj
    AuxMatArray[4].PushMatEl(1.0,0,0); 
    AuxMatArray[4].PushMatEl(0.0,1,1);
    AuxMatArray[4].PushMatEl(0.0,2,2); 
    AuxMatArray[4].PushMatEl(1.0,3,3); 
    AuxMatArray[4].PushMatEl(1.0,4,4); 
    AuxMatArray[4].PushMatEl(0.0,5,5);
    AuxMatArray[4].PushMatEl(0.0,6,6); 
    AuxMatArray[4].PushMatEl(1.0,7,7); 
  
  } else { 
    switch(STLNRGMats.size()){
    case 2: //fd_up fd_dn
      cout << " Thermodynamics only." << endl;
      break;

    default:
      cout << " Can't figure out what to do with the STLNRGMats... Exiting." << endl;
      exit(0);
    }
    // end switch STLNRGMats.size()
  }
  // end if MatName==Ndot, Sz

//   cout << " Ndot: " << endl;
//   AuxMatArray[2].PrintAllBlocks();
//   cout << " NMaj: " << endl;
//   AuxMatArray[4].PrintAllBlocks();

  /////
  // set up and diagonalize Hm1 8x8 matrix
  /////
  auxMat.CheckForMatEl=Diag_check;
  auxMat.CalcHNMatElCplx=OneChNupPdn_Hm1_Majorana_MatEl;
  auxMat.IsComplex=true;
  //  Block Structure
  auxMat.SyncNRGarray(AbasisHm1);
  auxMat.NeedOld=false;
  auxMat.UpperTriangular=true;

  //CNRGarray AeigHm1; // will be AeigHm1
  vector<double> ParamsHm1;
  ParamsHm1.push_back(t1);
  ParamsHm1.push_back(t2);
  ParamsHm1.push_back(phi_mag);
  ParamsHm1.push_back(em);

  //auxMat.DiagHN(ParamsHm1,&AbasisHm1,pSingleSite,&AuxMatArray[0],&AeigHm1);
  //AeigHm1.PrintEn();

  auxMat.DiagHN(ParamsHm1,&AbasisHm1,pSingleSite,&AuxMatArray[0],pAeig,true);
  
  pAeig->PrintEn();


  // Rotate c1,c2 (in the UnCut basis!)
  // Matrices are complex!

  // Check the loop below. Getting there!
  for (int imat=0;imat<AuxMatArray.size();imat++){
    //RotateMatrix((&AuxMatArray[imat]),(&AeigHm1),&auxMat);
    //RotateMatrix_NoCut((&AuxMatArray[imat]),(&AeigHm1),&auxMat,1);
    RotateMatrix_NoCut((&AuxMatArray[imat]),pAeig,&auxMat,1);
    AuxMatArray[imat].CopyData(&auxMat);
  }
  // end loop in matrices

//    cout << " f_{-1 up}, c_up: " << endl;
//   cout << " Ndot: " << endl;
//   AuxMatArray[2].PrintAllBlocks();

  // Set up operators in the Hm1 basis
  // ndot, Sz, S1 dot S2 if calcdens==0

  STLNRGMats[0].CopyData(&AuxMatArray[0]);  //f_{-1 up}
  STLNRGMats[1].CopyData(&AuxMatArray[1]);  //f_{-1 dn}

  if (	(strcmp(STLNRGMats[2].MatName,"cd_up")==0)&&
	(strcmp(STLNRGMats[3].MatName,"cd_dn")==0) ){
    // spectral functions
    STLNRGMats[2].CopyData(&AuxMatArray[0]); //c_up
    STLNRGMats[3].CopyData(&AuxMatArray[1]); //c_dn
    STLNRGMats[4].CopyData(&AuxMatArray[2]);  //f_dn
  } else if ( (strcmp(STLNRGMats[2].MatName,"Ndot")==0)&&
	     (strcmp(STLNRGMats[3].MatName,"Szdot")==0) ){
    STLNRGMats[2].CopyData(&AuxMatArray[2]); //Ndot  (rotated)
    STLNRGMats[3].CopyData(&AuxMatArray[3]); //Szdot (rotated)
    STLNRGMats[4].CopyData(&AuxMatArray[4]); //NMaj (rotated)
 
  }

  // Ok up to here.
 
  // Finally, build and diagonalize H_0 (N=0) and take it from there. 
  // Perhaps this can be done in the code itself.
  

}
// end OneChNupPdn_SetAndersonMajorana_H0
///////////////////
///////////////////
///////////////////
///////////////////

//////////////////////////////////
///                            ///
///  1chS Anderson model Hm1   ///
///                            ///
///  | S > basis               ///
///                            ///
//////////////////////////////////

void OneChS_SetAnderson_Hm1(vector<double> Params, 
			     CNRGarray* pAeig,
			     vector<CNRGmatrix> &STLNRGMats){


  // Set initial CNRG array (N=-1)

  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];
  // Testing:
  //Lambda=1.0;
  //HalfLambdaFactor=1.0;

  double U=Params[2];
  double ed=Params[3];

  pAeig->NQNumbers=1;
  pAeig->Nshell=-1;
  // new
  pAeig->totalS=true;
  pAeig->Sqnumbers.push_back(0);
  

  // |0>, |up dn>: |S=0>
  pAeig->QNumbers.push_back(0.0);

  // |up>, |dn>: |S=1/2>
  pAeig->QNumbers.push_back(0.5);
  //pAeig->iDegen.push_back(0); // S=1/2

  // Two states in block 1; one state in block 2
  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor)); // up


  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(1.0);

  pAeig->dEigVec.push_back(1.0);

  pAeig->PrintEn();

  // Set Matrix elements (reduced!)

  // fN (check -sqrt(2) and whatnot...

  //STLNRGMats[0].SyncNRGarray(*pAeig);
  // Sync all
  for (int imat=0;imat<STLNRGMats.size();imat++)
    STLNRGMats[imat].SyncNRGarray(*pAeig);

  STLNRGMats[0].UpperTriangular=false;

  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatEl.push_back(0.0);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);

  STLNRGMats[0].MatEl.push_back(0.0);
  STLNRGMats[0].MatEl.push_back(-sqrt(2.0));
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(2);
  STLNRGMats[0].MatBlockBegEnd.push_back(3);


  // Changing the way things are being done (June 2013)
  if (strcmp(STLNRGMats[1].MatName,"Ndot")==0){

    cout << " Setting up n_d, nd^2 and S^2... " << endl;
    // ndot    - imat=1 
    // ndot^2  - imat=2
    // Sdot^2  - imat=3
    // Setting all to Identity (all diagonal)
    for (int imat=1;imat<STLNRGMats.size();imat++){
      STLNRGMats[imat].SetDiagonalMatrix(0.0);
    }
    // Setting ndot, Sz, Sz2
    // Loop over diagonal
    for (int ibl=0;ibl<pAeig->NumBlocks();ibl++){
      double S=pAeig->GetQNumber(ibl,0); // S is set on each state.
      double ndot=0.0;
      if (ibl==1) ndot=1.0; // nel=1 for S=1/2 block
      for (int istbl=0;istbl<pAeig->GetBlockSize(ibl);istbl++){
	if (ibl==0) ndot=2*istbl; // nel=0, 2 in S=0 block
	STLNRGMats[1].PushBlockMatEl(ndot,ibl,ibl,istbl,istbl); // n
	STLNRGMats[2].PushBlockMatEl(ndot*ndot,ibl,ibl,istbl,istbl); //n^2
	STLNRGMats[3].PushBlockMatEl(S*(S+1.0),ibl,ibl,istbl,istbl); // S^2
      }
      // end loop over diag states
    }
    // end loop over blocks
    // Test
    //cout << " <0|nd|0> = " << STLNRGMats[1].GetMatEl(0,0) << endl;
  } else {
    if (strcmp(STLNRGMats[1].MatName,"cdot1")==0){
      cout << " Spectral function: setting up cd_red  ";
      STLNRGMats[1].CopyData(&STLNRGMats[0]);
      cout << " ... done." << endl;
    }
    // end if Mat[1]=cdot1
  }
  // end if Mat[1]=Ndot


  STLNRGMats[1].PrintAllBlocks();

  cout << " ... H_{-1} done. " << endl;

}
// end set H_-1 for Anderson with S symmetry

/////////////////

//////////////////////////////////
///                            ///
///  1chSz Anderson model Hm1  ///
///                            ///
///  | Sz > basis              ///
///                            ///
//////////////////////////////////

void OneChSz_SetAnderson_Hm1(vector<double> Params, 
			     CNRGarray* pAeig,
			     vector<CNRGmatrix> &STLNRGMats){


  // Set initial CNRG array (N=-1)

  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];
  // Testing:
  //Lambda=1.0;
  //HalfLambdaFactor=1.0;

  double U=Params[2];
  double ed=Params[3];

  double hz=0.5*Params[4]; // hz=g mu_B Bz: Zeeman term is:  - hz S_z 
                           // 0.5 added in Feb 2014

  pAeig->NQNumbers=1;
  pAeig->Nshell=-1;
  pAeig->totalS=false;
  
  // |0>, |up dn>: |Sz=0>
  pAeig->QNumbers.push_back(0.0);
  // |up>: |Sz=1/2>
  pAeig->QNumbers.push_back(0.5);
  // |dn>: |Sz=-1/2>
  pAeig->QNumbers.push_back(-0.5);

  // Two states in block 1; one state in blocks 2 and 3
  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U-hz)/(Lambda*HalfLambdaFactor)); // up
  pAeig->dEn.push_back((ed+0.5*U+hz)/(Lambda*HalfLambdaFactor)); // dn


  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(1.0);

  pAeig->dEigVec.push_back(1.0);

  pAeig->dEigVec.push_back(1.0);

  pAeig->PrintEn();

  // Set Matrix elements (reduced!)

  // fN (check -sqrt(2) and whatnot...

  //STLNRGMats[0].SyncNRGarray(*pAeig);
  // Sync all
  for (int imat=0;imat<STLNRGMats.size();imat++)
    STLNRGMats[imat].SyncNRGarray(*pAeig);


  // fN_up
  STLNRGMats[0].UpperTriangular=false;

  STLNRGMats[0].MatEl.push_back(1.0); //<0|fup|up>=1
  STLNRGMats[0].MatEl.push_back(0.0); //<up dn|fup|up>=0
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);

  STLNRGMats[0].MatEl.push_back(0.0); //<dn|fup|0>=0
  STLNRGMats[0].MatEl.push_back(1.0); //<dn|fup|up dn>=+1
  STLNRGMats[0].MatBlockMap.push_back(2); 
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(2);
  STLNRGMats[0].MatBlockBegEnd.push_back(3);


  // fN_dn
  STLNRGMats[1].UpperTriangular=false;

  STLNRGMats[1].MatEl.push_back(1.0); //<0|fdn|dn>=1
  STLNRGMats[1].MatEl.push_back(0.0); //<up dn|fdn|dn>=0
  STLNRGMats[1].MatBlockMap.push_back(0);
  STLNRGMats[1].MatBlockMap.push_back(2);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);

  STLNRGMats[1].MatEl.push_back(0.0); //<up|fdn|0>=0
  STLNRGMats[1].MatEl.push_back(-1.0); //<up|fdn| up dn>=-1
  STLNRGMats[1].MatBlockMap.push_back(1);
  STLNRGMats[1].MatBlockMap.push_back(0);
  STLNRGMats[1].MatBlockBegEnd.push_back(2);
  STLNRGMats[1].MatBlockBegEnd.push_back(3);


  // Changing the way things are being done (June 2013)
  if (strcmp(STLNRGMats[2].MatName,"Ndot")==0){
    cout << " Setting up n_d, Sz and Sz^2... " << endl;
    // ndot    - imat=2
    // ndot^2  - imat=3
    // Sdot^2  - imat=4
    // Setting all to Identity (all diagonal)
    for (int imat=2;imat<STLNRGMats.size();imat++){
      STLNRGMats[imat].SetDiagonalMatrix(0.0);
    }
    // Setting ndot, Sz, Sz2
    // Loop over diagonal
    for (int ibl=0;ibl<pAeig->NumBlocks();ibl++){
      double Sz=pAeig->GetQNumber(ibl,0); // Sz is set on each state.
      double ndot=0.0;
      if ( (ibl==1)||(ibl==2) ) ndot=1.0; // nel=1 for Sz=+-1/2 blocks
      for (int istbl=0;istbl<pAeig->GetBlockSize(ibl);istbl++){
	if (ibl==0) ndot=2*istbl; // nel=0, 2 in S=0 block
	STLNRGMats[2].PushBlockMatEl(ndot,ibl,ibl,istbl,istbl); // n
	STLNRGMats[3].PushBlockMatEl(Sz,ibl,ibl,istbl,istbl); // Sz
	STLNRGMats[4].PushBlockMatEl(Sz*Sz,ibl,ibl,istbl,istbl); // Sz^2
      }
      // end loop over diag states
    }
    // end loop over blocks
    // Test
    //cout << " <0|nd|0> = " << STLNRGMats[1].GetMatEl(0,0) << endl;
  } else if (strcmp(STLNRGMats[2].MatName,"cdotup")==0){
      cout << " Spectral function: setting up cd_up and cd_dn...  ";
      STLNRGMats[2].CopyData(&STLNRGMats[0]);
      STLNRGMats[3].CopyData(&STLNRGMats[1]);
      cout << " ... done." << endl;
    // end if Mat[1]=cdot1

  }else if (strcmp(STLNRGMats[2].MatName,"Szomega")==0){
      cout << " Dynamical spin susceptibiliy: setting up <Sz>_N=-1  ";
      // Needs work here


      for (int imat=2;imat<=5; imat++){
	STLNRGMats[imat].SyncNRGarray(*pAeig);
	// Block diagonal
	STLNRGMats[imat].SetDiagonalMatrix(0.0);
      }
      // end loop in matrices

      for (int ibl=0; ibl<pAeig->NumBlocks(); ibl++){
	double Szi=pAeig->GetQNumber(ibl,0);
	double ndot=0.0;

	if ( (ibl==1)||(ibl==2) ) ndot=1.0; // nel=1 for Sz=+-1/2 blocks

	cout << " Sz = " << Szi << endl;
	// Sz dynamical and Nd dynamical
	for (int istbl=0;istbl<pAeig->GetBlockSize(ibl);istbl++){
	  if (ibl==0) ndot=2*istbl; // nel=0, 2 in S=0 block
	  STLNRGMats[2].PushBlockMatEl(Szi,ibl,ibl,istbl,istbl);// <Sz>
	  STLNRGMats[3].PushBlockMatEl(ndot,ibl,ibl,istbl,istbl); // <Nd>
	}
      }
      // end loop in blocks (Set-up of Nd,Sz,Sz2)

      //Setting Sz2 and Nd2
      for (int imat=4; imat<=5; imat++){
      STLNRGMats[imat].CopyData(&STLNRGMats[imat-2]);
      for (int ii=0;ii<STLNRGMats[imat].MatEl.size();ii++)
	STLNRGMats[imat].MatEl[ii]*=STLNRGMats[imat].MatEl[ii];
      }
      // end setting Sz2 and Nd2

      cout << " ... done." << endl;

//       STLNRGMats[2].PrintAllBlocks();
//       STLNRGMats[3].PrintAllBlocks();


  }
  // end if Mat[1]=Ndot

  cout << " ... H_{-1} done. " << endl;

}
// end set H_-1 for Anderson with Sz symmetry

/////////////////


///////////////////////////////
///                         ///
///  1ch | Pup Pdn > basis  ///
///    (Jul 2017)           ///
///                         ///
///////////////////////////////



void OneChPupPdn_SetHm1_AndersonMajorana(vector<double> Params,
					CNRGbasisarray* pSingleSite,
					CNRGarray* pAeig, 
					vector<CNRGmatrix> &STLNRGMats){
  // In construction (05-Jul-2017)

  // Set initial CNRG array (N=-1)


  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];

  double Utilde=Params[2]/(Lambda*HalfLambdaFactor);
  double edtilde=Params[3]/(Lambda*HalfLambdaFactor);
  double gammatilde=Params[4];
  double hz=0.5*Params[5]/(Lambda*HalfLambdaFactor);; // hz=g mu_B Bz: Zeeman term is:  - hz S_z 

  // Remember:
  // gamma~=(2*Gamma/pi)^1/2/(sqrt(Lambda)*HalfLambdaFactor)
  // Why? see Wilson's paper (E. 2.18 with N=-1)

  // Note that "Lambda" (and not sqrt(Lambda)) appears in the
  // denominator of U~ and ed~.
  // This is different than the DQD code above since we are
  // generating only Hm1 here: H0 will be set up later by the NRG engine
  // and there will be a \sqrt(Lambda) multiplicating factor later on. 

  // Add other parameters here! t1, t2, phi, etc.
  
  double t1=Params[6]/(Lambda*HalfLambdaFactor);
  double t2=Params[7]/(Lambda*HalfLambdaFactor);
  double phi_mag=Params[8];
//   double em=Params[9]/(Lambda*HalfLambdaFactor);

  
  cout  << " Lambda= " << Lambda
	<< " HalfLambdaFactor= " << HalfLambdaFactor
	<< endl
	<< " U~= " << Utilde
	<< " ed~= " << edtilde
	<< " t~= " << gammatilde
	<< " hz~= " << hz << endl
	<< " t1~= " << t1
	<< " t2~= " << t2
	<< " phimag/Pi= " << phi_mag
	<< endl;


  // Build H_m1 (N=-1) 

  CNRGbasisarray AbasisHm1(2);

  // Diagonal energies
  double ediag[4];
  ediag[0]=0.5*Utilde;
  ediag[1]=edtilde+0.5*Utilde-hz;
  ediag[2]=edtilde+0.5*Utilde+hz;
  ediag[3]=2.0*edtilde+1.5*Utilde;

  // Set basis for Hm1
  AbasisHm1.Nshell=-1;
  AbasisHm1.NQNumbers=2; //(Pup,Pdn=(-1)^Ndn)
  AbasisHm1.totalS=false;

  // |0>d |0>f , |up>d |up>f, |dn>d |dn>f, |up dn>d |up dn>f : |Pup=+1 Pdn=+1>
  AbasisHm1.QNumbers.push_back(1.0);
  AbasisHm1.QNumbers.push_back(1.0);
  for (int ii=0;ii<4;ii++){
    AbasisHm1.dEn.push_back(ediag[ii]);
  } // end set diagonal energies
  AbasisHm1.BlockBegEnd.push_back(0);AbasisHm1.BlockBegEnd.push_back(3);

  // |0>d |up dn>f , |up>d |dn>f, |dn>d |up>f, |up dn>d |0>f : |Pup=-1 Pdn=-1>
  AbasisHm1.QNumbers.push_back(-1.0);
  AbasisHm1.QNumbers.push_back(-1.0);
  for (int ii=0;ii<4;ii++){
    AbasisHm1.dEn.push_back(ediag[ii]);
  } // end set diagonal energies
  AbasisHm1.BlockBegEnd.push_back(4);AbasisHm1.BlockBegEnd.push_back(7);

  // |0>d |up>f , |up>d |0>f, |dn>d |up dn>f, |up dn>d |dn>f : |Pup=-1 Pdn=+1>
  AbasisHm1.QNumbers.push_back(-1.0);
  AbasisHm1.QNumbers.push_back(1.0);
  for (int ii=0;ii<4;ii++){
    AbasisHm1.dEn.push_back(ediag[ii]);
  } // end set diagonal energies
  AbasisHm1.BlockBegEnd.push_back(8);AbasisHm1.BlockBegEnd.push_back(11);

  // |0>d |dn>f , |up>d |0>f, |dn>d |0>f, |up dn>d |up>f : |Pup=+1 Pdn=-1>
  AbasisHm1.QNumbers.push_back(1.0);
  AbasisHm1.QNumbers.push_back(-1.0);
  for (int ii=0;ii<4;ii++){
    AbasisHm1.dEn.push_back(ediag[ii]);
  } // end set diagonal energies
  AbasisHm1.BlockBegEnd.push_back(12);AbasisHm1.BlockBegEnd.push_back(15);

  AbasisHm1.SetEigVecToOne();

  // Set up cd and fd: easier in this basis. 
  // Will Rotate later.
  for (int imat=0;imat<STLNRGMats.size();imat++)
    STLNRGMats[imat].SyncNRGarray(AbasisHm1);


  AbasisHm1.PrintAll();


  vector<CNRGmatrix> AuxMatArray; // Arrays needed to set up H0 
                       //(cd_up, cd_dn, f_Maj_up and f_Maj_dn)
  CNRGmatrix auxMat;

  auxMat=STLNRGMats[0]; // Same rules as fN_up
  auxMat.SetZeroMatrix(); // Build block structure
  AuxMatArray.push_back(auxMat); // cd_up
  auxMat=STLNRGMats[1]; // Same rules as fN_dn
  auxMat.SetZeroMatrix(); // Build block structure
  AuxMatArray.push_back(auxMat); // cd_dn

  // c_up|up>|k>=+|0>|k>  
  // c_up|up dn>|k>=+|dn>|k>  
  
  AuxMatArray[0].PushMatEl(1.0,8,1); 
  AuxMatArray[0].PushMatEl(1.0,12,5);
  AuxMatArray[0].PushMatEl(1.0,0,9); 
  AuxMatArray[0].PushMatEl(1.0,4,13); 
  AuxMatArray[0].PushMatEl(1.0,10,3); 
  AuxMatArray[0].PushMatEl(1.0,14,7);
  AuxMatArray[0].PushMatEl(1.0,2,11); 
  AuxMatArray[0].PushMatEl(1.0,6,15); 
   
  // c_dn|dn>|k>=+|0>|k>  
  // c_dn|up dn>|k>=-|up>|k>  
  AuxMatArray[1].PushMatEl(1.0,12,2); 
  AuxMatArray[1].PushMatEl(1.0,8,6); 
  AuxMatArray[1].PushMatEl(1.0,4,10); 
  AuxMatArray[1].PushMatEl(1.0,0,14); 
  AuxMatArray[1].PushMatEl(-1.0,13,3); 
  AuxMatArray[1].PushMatEl(-1.0,9,7); 
  AuxMatArray[1].PushMatEl(-1.0,5,11); 
  AuxMatArray[1].PushMatEl(-1.0,1,15); 

 
  if (	(strcmp(STLNRGMats[2].MatName,"cd_up")==0)&&
	(strcmp(STLNRGMats[3].MatName,"cd_dn")==0) ){

    AuxMatArray.push_back(auxMat); // [2] f_Maj_up
    AuxMatArray.push_back(auxMat); // [3] f_Maj_dn

    // cd1, cd2 if calcdens==1 (just use the ops defined above)
    // Define c_up and c_dn and f_Maj_dn matrices in this basis!
    
    cout << " Majorana: spectral functions still at implementation stage..." << endl;

    // f_up|up dn (0) >|up>=+|up dn (0)>|0>  
    // f_up|up dn (0) >|up dn>=+|up dn (0)>|dn>  
    // f_up|up (dn)>|up>=-|up (dn)>|0>  
    // f_up|up (dn)>|up dn>=-|up (dn)>|dn>  
    AuxMatArray[2].PushMatEl(1.0,0,8); 
    AuxMatArray[2].PushMatEl(1.0,7,15); 
    AuxMatArray[2].PushMatEl(1.0,11,3); 
    AuxMatArray[2].PushMatEl(1.0,12,4); 
    AuxMatArray[2].PushMatEl(-1.0,9,1); 
    AuxMatArray[2].PushMatEl(-1.0,2,10); 
    AuxMatArray[2].PushMatEl(-1.0,14,6); 
    AuxMatArray[2].PushMatEl(-1.0,5,13); 


    // f_dn|up dn (0) >|dn>=+|up dn (0)>|0>  
    // f_dn|up (dn)>|up dn>=+|up (dn)>|up>  
    // f_dn|up (dn)>|dn>=-|up (dn)>|0>  
    // f_dn|up dn (0) >|up dn>=-|up dn (0)>|up>  
    AuxMatArray[3].PushMatEl(1.0,7,11); 
    AuxMatArray[3].PushMatEl(1.0,0,12); 
    AuxMatArray[3].PushMatEl(1.0,6,10); 
    AuxMatArray[3].PushMatEl(1.0,1,13); 
    AuxMatArray[3].PushMatEl(-1.0,15,3); 
    AuxMatArray[3].PushMatEl(-1.0,8,4); 
    AuxMatArray[3].PushMatEl(-1.0,14,2); 
    AuxMatArray[3].PushMatEl(-1.0,9,5); 
    
//     cout << " f_dn: " << endl;
//     AuxMatArray[3].PrintAllBlocks();

    // Stopped here (6/7/17)

  }else if ( (strcmp(STLNRGMats[2].MatName,"Ndot")==0)&&
	     (strcmp(STLNRGMats[3].MatName,"Szdot")==0) ){
    //////////////////
    
    cout << " Setting up 1-body observables..." << endl;

    auxMat=STLNRGMats[2]; // Same rules as Ndot
    auxMat.SetZeroMatrix(); // Build block structure
  
    AuxMatArray.push_back(auxMat); // [2] Ndot
    AuxMatArray.push_back(auxMat); // [3] Szdot
    AuxMatArray.push_back(auxMat); // [4] NMaj
  
//    for (int imat=2; imat<=4; imat++){
//      AuxMatArray[imat]=STLNRGMats[imat]; // Same rules as ndot, szdot, nmaj
//      AuxMatArray[imat].SetZeroMatrix(); // Build block structure
//     }
 
    // ndot=nup+ndown [2]
    // ndot=0 for ist=0,4,8,12 : ist mod 4 = 0
    // ndot=1 for ist=1,5,9,13 : ist mod 4 = 1
    // ndot=1 for ist=2,6,10,14: ist mod 4 = 2
    // ndot=2 for ist=3,7,11,15: ist mod 4 = 3
    // szdot          [3]
    // szdot=0   for ist=0,4,8,12 : ist mod 4 = 0
    // szdot=0.5 for ist=1,5,9,13 : ist mod 4 = 1
    // szdot=-0.5 for ist=2,6,10,14: ist mod 4 = 2
    // szdot=0 for ist=3,7,11,15: ist mod 4 = 3

    int ist=0;
    for (int ist=0;ist<16;ist++){
      int istmod4=ist % 4;
      switch (istmod4){
      case 0:
	AuxMatArray[2].PushMatEl(0.0,ist,ist); 
	AuxMatArray[3].PushMatEl(0.0,ist,ist); 
	break;
      case 1:
	AuxMatArray[2].PushMatEl(1.0,ist,ist); 
	AuxMatArray[3].PushMatEl(0.5,ist,ist); 
	break;
      case 2:
	AuxMatArray[2].PushMatEl(1.0,ist,ist); 
	AuxMatArray[3].PushMatEl(-0.5,ist,ist); 
	break;
      case 3:
	AuxMatArray[2].PushMatEl(2.0,ist,ist); 
	AuxMatArray[3].PushMatEl(0.0,ist,ist); 
	break;
      default:
	AuxMatArray[2].PushMatEl(0.0,ist,ist); 
	AuxMatArray[3].PushMatEl(0.0,ist,ist); 
      } // end switch
    } //end loop in ist
   
    // nMaj
    AuxMatArray[4].PushMatEl(0.0,0,0); 
    AuxMatArray[4].PushMatEl(1.0,1,1);
    AuxMatArray[4].PushMatEl(1.0,2,2); 
    AuxMatArray[4].PushMatEl(2.0,3,3); 
    AuxMatArray[4].PushMatEl(2.0,4,4); 
    AuxMatArray[4].PushMatEl(1.0,5,5);
    AuxMatArray[4].PushMatEl(1.0,6,6); 
    AuxMatArray[4].PushMatEl(0.0,7,7); 
    AuxMatArray[4].PushMatEl(1.0,8,8); 
    AuxMatArray[4].PushMatEl(0.0,9,9); 
    AuxMatArray[4].PushMatEl(2.0,10,10); 
    AuxMatArray[4].PushMatEl(1.0,11,11); 
    AuxMatArray[4].PushMatEl(1.0,12,12); 
    AuxMatArray[4].PushMatEl(2.0,13,13); 
    AuxMatArray[4].PushMatEl(0.0,14,14); 
    AuxMatArray[4].PushMatEl(1.0,15,15); 
  
  } else { 
    switch(STLNRGMats.size()){
    case 2: //fd_up fd_dn
      cout << " Thermodynamics only." << endl;
      break;

    default:
      cout << " Can't figure out what to do with the STLNRGMats... Exiting." << endl;
      exit(0);
    }
    // end switch STLNRGMats.size()
  }
  // end if MatName==Ndot, Sz

//    cout << " Ndot: " << endl;
//    AuxMatArray[2].PrintAllBlocks();
//    cout << " NMaj: " << endl;
//    AuxMatArray[4].PrintAllBlocks();

  /////
  // set up and diagonalize Hm1 16x16 matrix
  /////
  auxMat.CheckForMatEl=Diag_check;
  auxMat.CalcHNMatElCplx=OneChPupPdn_Hm1_Majorana_MatEl;
  auxMat.IsComplex=true;
  //  Block Structure
  auxMat.SyncNRGarray(AbasisHm1);
  auxMat.NeedOld=false;
  auxMat.UpperTriangular=true;

  //CNRGarray AeigHm1; // will be AeigHm1
  vector<double> ParamsHm1;
  ParamsHm1.push_back(t1);
  ParamsHm1.push_back(t2);
  //ParamsHm1.push_back(phi_mag*pi); //actual phi 
  ParamsHm1.push_back(phi_mag); //phi/Pi 
  //ParamsHm1.push_back(em);

  //auxMat.DiagHN(ParamsHm1,&AbasisHm1,pSingleSite,&AuxMatArray[0],&AeigHm1);
  //AeigHm1.PrintEn();

  auxMat.DiagHN(ParamsHm1,&AbasisHm1,pSingleSite,&AuxMatArray[0],pAeig,true);
  
  pAeig->PrintEn();

//   cout << " Main_SetHm1 c_up BEFORE rotation: " << endl;
//   AuxMatArray[0].PrintAllBlocks();
//   cout << " Main_SetHm1 c_dn BEFORE rotation: " << endl;
//   AuxMatArray[1].PrintAllBlocks();

  // Rotate c1,c2 (in the UnCut basis!)
  // Matrices are complex!

  // Check the loop below. Getting there!
  for (int imat=0;imat<AuxMatArray.size();imat++){
    //RotateMatrix((&AuxMatArray[imat]),(&AeigHm1),&auxMat);
    //RotateMatrix_NoCut((&AuxMatArray[imat]),(&AeigHm1),&auxMat,1);
    RotateMatrix_NoCut((&AuxMatArray[imat]),pAeig,&auxMat,1);
    AuxMatArray[imat].CopyData(&auxMat);
  }
  // end loop in matrices

//    cout << " Main_SetHm1 c_up rotated: " << endl;
//    AuxMatArray[0].PrintAllBlocks();
//    cout << " Main_SetHm1 c_dn rotated: " << endl;
//    AuxMatArray[1].PrintAllBlocks();


  // Set up operators in the Hm1 basis
  // ndot, Sz, S1 dot S2 if calcdens==0

  STLNRGMats[0].CopyData(&AuxMatArray[0]);  //f_{-1 up}
  STLNRGMats[1].CopyData(&AuxMatArray[1]);  //f_{-1 dn}

  if (	(strcmp(STLNRGMats[2].MatName,"cd_up")==0)&&
	(strcmp(STLNRGMats[3].MatName,"cd_dn")==0) ){
    // spectral functions
    STLNRGMats[2].CopyData(&AuxMatArray[0]); //c_up
    STLNRGMats[3].CopyData(&AuxMatArray[1]); //c_dn
    STLNRGMats[4].CopyData(&AuxMatArray[2]);  //f_up
    STLNRGMats[5].CopyData(&AuxMatArray[3]);  //f_dn
  } else if ( (strcmp(STLNRGMats[2].MatName,"Ndot")==0)&&
	     (strcmp(STLNRGMats[3].MatName,"Szdot")==0) ){
    STLNRGMats[2].CopyData(&AuxMatArray[2]); //Ndot  (rotated)
    STLNRGMats[3].CopyData(&AuxMatArray[3]); //Szdot (rotated)
    STLNRGMats[4].CopyData(&AuxMatArray[4]); //NMaj (rotated)
 
  }

  // Ok up to here.
 
  // Finally, build and diagonalize H_0 (N=0) and take it from there. 
  // Perhaps this can be done in the code itself.
  

}
// end OneChPupPdn_SetAndersonMajorana_Hm1
///////////////////
///////////////////
///////////////////
///////////////////







/////////////////////////
//  Chain hamiltonians //
/////////////////////////

////
// OneChQS
////

void OneChQS_SetChainH0(CNRGarray* pAeig, 
			vector<CNRGmatrix> &STLNRGMats){

  // First Test: the usual Anderson model

  pAeig->NQNumbers=2;
  pAeig->Nshell=0;
  // new
  pAeig->totalS=true;
  pAeig->Sqnumbers.push_back(1);
  
  pAeig->QNumbers.push_back(-1.0); // |0>
  pAeig->QNumbers.push_back(0.0);


  pAeig->QNumbers.push_back(0.0);  // | up >
  pAeig->QNumbers.push_back(0.5);

  pAeig->QNumbers.push_back(1.0); // |0>
  pAeig->QNumbers.push_back(0.0);


  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);

  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);

  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);

  // Set Matrix elements

  // fN reduced

  STLNRGMats[0].SyncNRGarray(*pAeig);
  STLNRGMats[0].UpperTriangular=false;
 
  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);

  STLNRGMats[0].MatEl.push_back(-sqrt(2.0));
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockMap.push_back(2);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);

 

}


/////////////////



////
// OneChQSz
////

void OneChQSz_SetChainH0(CNRGarray* pAeig, 
			vector<CNRGmatrix> &STLNRGMats){

  pAeig->NQNumbers=2;
  pAeig->Nshell=0;
  pAeig->totalS=false;
  
  // |0> = |Q=-1, Sz=0>
  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(0.0);
  pAeig->iDegen.push_back(0); // Sz=0  

  // |up> = |Q=0, Sz=1/2>
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.5);
  pAeig->iDegen.push_back(0); // Sz=1/2

  // |dn> = |Q=0, Sz=-1/2>
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(-0.5);
  pAeig->iDegen.push_back(0); // Sz=1/2

  // |up dn> = |Q=1, Sz=0>
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(0.0);
  pAeig->iDegen.push_back(0); // Sz=0

  // One state per block
  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);

  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);

  // Set Matrix elements

  // fN_up

  STLNRGMats[0].SyncNRGarray(*pAeig);
  STLNRGMats[0].UpperTriangular=false;
  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);

  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(2);
  STLNRGMats[0].MatBlockMap.push_back(3);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);


  // fN_dn
  STLNRGMats[1].SyncNRGarray(*pAeig);
  STLNRGMats[1].UpperTriangular=false;
  STLNRGMats[1].MatEl.push_back(1.0);
  STLNRGMats[1].MatBlockMap.push_back(0);
  STLNRGMats[1].MatBlockMap.push_back(2);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);

  STLNRGMats[1].MatEl.push_back(-1.0);
  STLNRGMats[1].MatBlockMap.push_back(1);
  STLNRGMats[1].MatBlockMap.push_back(3);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);

}
// end OneChQSz_SetChain

/////////////////


///////////////////
// TwoChQS chain //
///////////////////

void TwoChQS_SetH0Chain(CNRGarray* pAeig, 
			CNRGbasisarray* pSingleSite, 
			vector<CNRGmatrix> &STLNRGMats){


 
  // Aeig will need
  // - Energies
  // - Eigenvectors


  CNRGbasisarray AbasisH0;
  AbasisH0.ClearAll();
  pAeig->ClearAll();


//   // Set AbasisH0
//   TwoChQSNoSz_SetSingleSite(&AbasisH0);

//   // Need to do this.
//   for (int ii=0;ii<AbasisH0.StCameFrom.size();ii++)
//     AbasisH0.StCameFrom[ii]=0;

  AbasisH0.Nshell=0;
  AbasisH0.NQNumbers=2;

  // |-2 0>  = |0 0> 
  AbasisH0.QNumbers.push_back(-2.0);
  AbasisH0.QNumbers.push_back(0.0);
  AbasisH0.BlockBegEnd.push_back(0);AbasisH0.BlockBegEnd.push_back(0);
  AbasisH0.iType.push_back(0);


  // |-1 1/2>  = |0 up>, |up 0> 
  AbasisH0.QNumbers.push_back(-1.0);
  AbasisH0.QNumbers.push_back(0.5);
  AbasisH0.BlockBegEnd.push_back(1);AbasisH0.BlockBegEnd.push_back(2);
  AbasisH0.iType.push_back(1);
  AbasisH0.iType.push_back(2);


  // |0 0>  = |singlet>, |0  up dn>. |up dn  0> 
  AbasisH0.QNumbers.push_back(0.0);
  AbasisH0.QNumbers.push_back(0.0);
  AbasisH0.BlockBegEnd.push_back(3);AbasisH0.BlockBegEnd.push_back(5);
  AbasisH0.iType.push_back(5);
  AbasisH0.iType.push_back(6);
  AbasisH0.iType.push_back(7);


  // |0 1>  = |triplet > = |up up>
  AbasisH0.QNumbers.push_back(0.0);
  AbasisH0.QNumbers.push_back(1.0);
  AbasisH0.BlockBegEnd.push_back(6);AbasisH0.BlockBegEnd.push_back(6);
  AbasisH0.iType.push_back(10);


  // |1 1/2>  = |up  up dn>. |up dn up> 
  AbasisH0.QNumbers.push_back(1.0);
  AbasisH0.QNumbers.push_back(0.5);
  AbasisH0.BlockBegEnd.push_back(7);AbasisH0.BlockBegEnd.push_back(8);
  AbasisH0.iType.push_back(11);
  AbasisH0.iType.push_back(12);


  // |2 0>  = |up dn  up dn>
  AbasisH0.QNumbers.push_back(2.0);
  AbasisH0.QNumbers.push_back(0.0);
  AbasisH0.BlockBegEnd.push_back(9);AbasisH0.BlockBegEnd.push_back(9);
  AbasisH0.iType.push_back(15);



  // Aeig is the same

  *pAeig=AbasisH0;
  
  for (int ibl=0;ibl<pAeig->NumBlocks();ibl++){
    for (int ii=pAeig->GetBlockLimit(ibl,0);
	 ii<=pAeig->GetBlockLimit(ibl,1);ii++){
      pAeig->dEn.push_back(0); // all the same
      AbasisH0.StCameFrom.push_back(0);
      for (int jj=pAeig->GetBlockLimit(ibl,0);
	   jj<=pAeig->GetBlockLimit(ibl,1);jj++){
	if (ii==jj)
	  pAeig->dEigVec.push_back(1.0);
	else
	  pAeig->dEigVec.push_back(0.0);
      }
    }
    //loop in ii
  }
  // Loop in blocks
  // Type labels the state

  STLNRGMats[0].SyncNRGarray(AbasisH0);
  STLNRGMats[1].SyncNRGarray(AbasisH0);

  for (int imat=0; imat<=1;imat++)
    STLNRGMats[imat].SetMatrix(&AbasisH0,pSingleSite);

  // No need to Rotate (Eigenvectors form a unitary matrix)
  cout << " f0_ch1 : " << endl;
  //STLNRGMats[0].PrintAllBlocks();
  STLNRGMats[0].PrintMatBlock(1,3);

}

/////////////////


////
// OneChNupPdn
////

void OneChNupPdn_SetChainH0(CNRGarray* pAeig, 
			vector<CNRGmatrix> &STLNRGMats){

  pAeig->NQNumbers=2;
  pAeig->Nshell=0;
  pAeig->totalS=false;
  
  // |0> = |Nup=0, Pdn=1>
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(1.0);
  pAeig->iDegen.push_back(0); // Sz=0  

  // |up> = |Nup=1, Pdn=1>
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(1.0);
  pAeig->iDegen.push_back(0); // Sz=1/2

  // |dn> = |Nup=0, Pdn=-1>
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(-1.0);
  pAeig->iDegen.push_back(0); // Sz=1/2

  // |up dn> = |Nup=1, Pdn=-1>
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(-1.0);
  pAeig->iDegen.push_back(0); // Sz=0

  // One state per block
  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);

  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);

  // Set Matrix elements

  // fN_up

  STLNRGMats[0].SyncNRGarray(*pAeig);
  STLNRGMats[0].UpperTriangular=false;
  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);

  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(2);
  STLNRGMats[0].MatBlockMap.push_back(3);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);


  // fN_dn
  STLNRGMats[1].SyncNRGarray(*pAeig);
  STLNRGMats[1].UpperTriangular=false;
  STLNRGMats[1].MatEl.push_back(1.0);
  STLNRGMats[1].MatBlockMap.push_back(0);
  STLNRGMats[1].MatBlockMap.push_back(2);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);

  STLNRGMats[1].MatEl.push_back(-1.0);
  STLNRGMats[1].MatBlockMap.push_back(1);
  STLNRGMats[1].MatBlockMap.push_back(3);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);

}
// end OneChNupPdn_SetChain


////
// OneChPupPdn
////

void OneChPupPdn_SetChainH0(CNRGarray* pAeig, 
			vector<CNRGmatrix> &STLNRGMats){

  pAeig->NQNumbers=2;
  pAeig->Nshell=0;
  pAeig->totalS=false;
  
  // |0> = |Pup=1, Pdn=1>
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(1.0);
  pAeig->iDegen.push_back(0);   

  // |up> = |Pup=-1, Pdn=1>
  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(1.0);
  pAeig->iDegen.push_back(0);

  // |dn> = |Pup=1, Pdn=-1>
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(-1.0);
  pAeig->iDegen.push_back(0);

  // |up dn> = |Pup=-1, Pdn=-1>
  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(-1.0);
  pAeig->iDegen.push_back(0); 

  // One state per block
  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);

  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);

  // Set Matrix elements

  // fN_up

  STLNRGMats[0].SyncNRGarray(*pAeig);
  STLNRGMats[0].UpperTriangular=false;
  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(0);
  STLNRGMats[0].MatBlockMap.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);
  STLNRGMats[0].MatBlockBegEnd.push_back(0);

  STLNRGMats[0].MatEl.push_back(1.0);
  STLNRGMats[0].MatBlockMap.push_back(2);
  STLNRGMats[0].MatBlockMap.push_back(3);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);
  STLNRGMats[0].MatBlockBegEnd.push_back(1);


  // fN_dn
  STLNRGMats[1].SyncNRGarray(*pAeig);
  STLNRGMats[1].UpperTriangular=false;
  STLNRGMats[1].MatEl.push_back(1.0);
  STLNRGMats[1].MatBlockMap.push_back(0);
  STLNRGMats[1].MatBlockMap.push_back(2);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);
  STLNRGMats[1].MatBlockBegEnd.push_back(0);

  STLNRGMats[1].MatEl.push_back(-1.0);
  STLNRGMats[1].MatBlockMap.push_back(1);
  STLNRGMats[1].MatBlockMap.push_back(3);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);
  STLNRGMats[1].MatBlockBegEnd.push_back(1);

}
// end OneChPupPdn_SetChain




/////////////////
// Trash

//   // Complex number test
//   cout << "Complex Number test" << endl;
//   complex<double> A1 (1.0,0.0);
//   cout << "A1 = " << A1 << endl;



//   auxMat.SetZeroMatrix(); // Build block structure
//   cout << "Got here II " << endl;

//   cout << "Mat(0) = " << auxMat.MatElCplx[0] << endl;
//   cout << "Mat(1) = " << auxMat.MatElCplx[1] << endl;
//   cout << "Mat(2) = " << auxMat.MatElCplx[2] << endl;

 
//   A1.real()=1.1234;
//   cout << "A1 = " << A1 << endl;
//   double A2=1.234;
//   auxMat.TemplateTest(A2);
//   auxMat.TemplateTest(A1);
//   A1.real()=1.0;
//   auxMat.PushMatEl(A1,0,0);
//   A1=OneImC;
//   auxMat.PushMatEl(A1,0,1);
//   A1=ZeroC; A1.real()=2.0;
//   auxMat.PushMatEl(A1,1,1);

//   auxMat.PrintAllBlocks();

//   // Use zheev instead!

//   vector<double> denergies;
//   vector<complex<double> > dzvec;

//   auxMat.DiagBlock(0,denergies,dzvec);


//   cout << " e1 = " << denergies[0] << endl
//        << " e2 = " << denergies[1] << endl;

//   cout << " z1a = " << dzvec[0]
//        << "  z1b = " << dzvec[1] << endl
//        << " z2a = " << dzvec[2] 
//        << "  z2b = " << dzvec[3] << endl;


//   cout << "Got here III " << endl;

//   auxMat.DiagBlock(0,pAeig->dEn,pAeig->cEigVec);

//   pAeig->PrintAll();
