#include <iostream>
//#include <iomanip>
//#include <fstream>

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include <vector>
//#include <cmath>
//#include <cstring>

//#include <unistd.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "NRGOpMatRules.hpp"
#include "NRG_main.hpp"
#include "TwoChQS.hpp"


#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif


void CNRGCodeHandler::ModelSwitch(  vector<int> &CommonQNs,
				    vector<int> &totSpos,   
				    CNRGarray* pAeig,
				    CNRGbasisarray* pSingleSite, 
				    CNRGmatrix* pHN, 
				    vector<CNRGmatrix> &STLMatArray,
				    CNRGthermo* ThermoArray,
				    double &chi_m1){

  // Ok, this can be changed later to reduce the number of parameters 
  // Leave as it is for now.

  // List of Non CNRGcodehandler objects that needs to be set/defined

  // ALL of these have to be DEFINED as fields of CNRGCodeHandler
  // STLMatArray -> MatArray (but routines need an STL vector... 
  //                          SHOULD DEFINE MatArray as an STL vector in 
  //                          CNRGcodehandler 
  //                          and adapt in the other routines. 
  //                          But leave it for now)
  //                MatArray enters with 4 elements. Add if you need more.
  // CommonQNs   -> (input)
  // totSpos     -> (input)
  // HN          -> (input)
  // NumNRGarrays (not really, just set it as NumNRGmats)
  // ThermoArray -> ThermoSTLArray (ok)

  // Defined here:
  
  CNRGmatrix auxNRGMat;
  vector <double> Params;

  double chi2_m1=0.0;

  // To do: getrid of ThermoArray and STLMatArray: 
  // define them in codehandler class.
  // get rid of all parameters?



  ///////////////////////////////////
  ////   Model specific hamiltonians ///
  ///////////////////////////////////

// SymNo:
//  0 - OneChQS
//  1 - TwoChQS
//  2 - OneChQSz
//  3 - OneChQ
//  4 - TwoChQSP
//  5 - TwoChQSz
//  6 - OneChNupPdn
//  7 - OneChS
//  8 - OneChSz
//  9 - OneChPupPdn

// ModelNo:
//  0 - Anderson
//  1 - Kondo
//  2 - Anderson w Holstein phonons
//  3 - Chain only
//  4 - CM phonons (2ch only)
//  5 - SMM (1chQ only)
//  6 - DQD 
//  7 - Anderson+Local Majorana (1ch NupPdn only)
//

// NEW: Add zEQ(z) to SuscepImp names.
//

    char zstring[30];
    char zvalue[8];
    strcpy(zstring,"_zEQ");
    if (dEqual(chain.z_twist,1.0))
      strcpy(zvalue,"1");
    else
      sprintf(zvalue,"%4.2f",chain.z_twist);
    strcat(zstring,zvalue);
    cout << "Zstring = " << zstring << endl;

    switch(SymNo){
    case 0: // OneChQS
      ////////////////////////////
      /////                   ////
      ///// Symmetry: OneChQS ////
      /////                   ////
      ////////////////////////////
      /// Models so far:       ///
      ///  0 - Anderson        ///
      ///  1 - Kondo           ///
      ///  3 - Chain           ///
      ///  6 - DQD             ///
      ////////////////////////////
      switch(ModelNo){
      case 0: // Single-impurity Anderson model Q,S symmetry

	// Initialize matrices (OneChQ routines work here!)
	// Jul 09: using STMatArray

	//NumThermoMats=2;
      
	// fN (reduced)
	STLMatArray[0].NeedOld=false;
	STLMatArray[0].UpperTriangular=false;
	STLMatArray[0].CheckForMatEl=OneChQS_cd_check;
	STLMatArray[0].CalcMatEl=OneChQS_fN_MatEl;
	STLMatArray[0].SaveMatYN=false;

	switch (calcdens){
	case 0: // OpAvg: n_d(T), n^2_d(T), S^2(T) - NEW (2012)
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false; // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;
	  STLMatArray[1]=auxNRGMat; // 1 - ndot (there already)
	  STLMatArray[2]=auxNRGMat; // 2 - ndot^2 (there already)
	  STLMatArray[3]=auxNRGMat; // 3 - Sdot^2 (there already)
	  strcpy(STLMatArray[1].MatName,"Ndot");
	  strcpy(STLMatArray[2].MatName,"NdSq");
	  strcpy(STLMatArray[3].MatName,"SdSq");
	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; // Need to do this after
	  // changes in STLMatArray!
	  SaveData=false;
// 	  STLMatArray.pop_back(); // 1 Matrices only
// 	  STLMatArray.pop_back(); // 1 Matrices only
// 	  STLMatArray.pop_back(); // 1 Matrices only
	  // No Thermo
	  NumThermoMats=0;
	  break;
	case 1: //Spectral density
	  STLMatArray.pop_back(); // 2 Matrices
	  STLMatArray.pop_back(); // 2 Matrices
	  NumThermoMats=0;
	  // cd operator 
	  STLMatArray[1].NeedOld=true;
	  STLMatArray[1].UpperTriangular=false;
	  STLMatArray[1].CheckForMatEl=OneChQS_cd_check;
	  STLMatArray[1].CalcMatEl=OneChQS_cd_MatEl;
	  STLMatArray[1].SaveMatYN=true;
	  strcpy(STLMatArray[1].MatName,"cdot1");
	  NumNRGmats=STLMatArray.size();
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"1chQSAnderson");
	  break;
	case 2: // Thermodynamics ONLY
//	  STLMatArray.pop_back(); // 1 Matrix only
//	  STLMatArray.pop_back(); // 1 Matrix only
//	  STLMatArray.pop_back(); // 1 Matrix only
// If you want <n_d>(T), use calcdens=0
// June 2013: I'm changing this. Want to calculate ndot toghether with chi_imp
//

 	  auxNRGMat.NeedOld=true;
 	  auxNRGMat.UpperTriangular=false;
	  // update procedure only works for 'false'
 	  auxNRGMat.CalcAvg=true;  // CalcAvg
 	  auxNRGMat.CheckForMatEl=Diag_check;
 	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;
 	  STLMatArray[1]=auxNRGMat; // 1 - ndot (there already)
 	  STLMatArray[2]=auxNRGMat; // 2 - ndot^2 (not right now)
 	  STLMatArray[3]=auxNRGMat; // 3 - Sdot^2 (not right now)
 	  strcpy(STLMatArray[1].MatName,"Ndot");
 	  strcpy(STLMatArray[2].MatName,"NdSq");
 	  strcpy(STLMatArray[3].MatName,"SdSq"); 
////////////
	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; // Need to do this after
	  // changes in STLMatArray!
	  SaveData=false;
	  // Thermo
	  // Z trick in file name.
	  // Changinf chain arq names for 1chQS too (Sep 2010)
	  strcpy(ThermoArray[0].ArqName,"SuscepImp1ChQS_Anderson");
	  strcat(ThermoArray[0].ArqName,zstring);
	  strcat(ThermoArray[0].ArqName,".dat");
	  //
	  strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QS");
	  strcat(ThermoArray[0].ChainArqName,zstring);
	  strcat(ThermoArray[0].ChainArqName,".dat");
	  //
	  strcpy(ThermoArray[1].ArqName,"EntropyImp1ChQS_Anderson");
	  strcat(ThermoArray[1].ArqName,zstring);
	  strcat(ThermoArray[1].ArqName,".dat");
	  //
	  strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QS");
	  strcat(ThermoArray[1].ChainArqName,zstring);
	  strcat(ThermoArray[1].ChainArqName,".dat");
	  ThermoArray[0].Calc=CalcSuscep;
	  ThermoArray[0].dImpValue=1.0/8.0; 
	  ThermoArray[1].Calc=CalcEntropy;
	  ThermoArray[1].dImpValue=2.0*log(2.0);
	  NumThermoMats=2;
	  break;
	case 4: // Dynamical spin susceptibility and <Sz> and <Sz^2> (2015)
	  // 5 Matrices
	  STLMatArray.push_back(auxNRGMat);
	  NumThermoMats=0;
	  // <Sz> and <Sz2> operator 
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false;
	  // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;

	  STLMatArray[1]=auxNRGMat; // 1 - Sz dynamical
	  STLMatArray[1].SaveMatYN=true;
	  STLMatArray[1].CalcAvg=false;
	  STLMatArray[1].CalcMatEl=OneChQS_Sz_MatEl;
	  strcpy(STLMatArray[1].MatName,"Szomega");
	  STLMatArray[1].WignerEckartL=1.0; // Needed by DM_NRG later


	  STLMatArray[2]=auxNRGMat; // 2 - Nd dynamical
	  STLMatArray[2].SaveMatYN=true;
	  STLMatArray[2].CalcAvg=false;
	  strcpy(STLMatArray[2].MatName,"Ndomega");

	  STLMatArray[3]=auxNRGMat; // 3 - <Sz2>(T) static
	  STLMatArray[3].SaveMatYN=false;
	  strcpy(STLMatArray[3].MatName,"Sz2dot");

	  STLMatArray[4]=auxNRGMat; // 4 - <Nd2>(T) static
	  STLMatArray[4].SaveMatYN=false;
	  strcpy(STLMatArray[4].MatName,"Nd2dot");

	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; // Need to do this after
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"1chQSAnderson");
	  break;
	default:
	  cout << " calcdens = " << calcdens << " not implemented. Exiting. " << endl;
	  exit(0);

	}
	// end switch calcdens

	// Param for H0
	Params.push_back(Lambda);
	Params.push_back(HalfLambdaFactor);
	Params.push_back(dInitParams[0]); // U
	Params.push_back(dInitParams[2]); // ed
	Params.push_back(chi_m1); // gamma

	OneChQS_SetAnderson_Hm1(Params,pAeig,STLMatArray);


	// BuildBasis params
	CommonQNs.push_back(2); // Q and S and commont QNs
	CommonQNs.push_back(0); // position of Q in old basis
	CommonQNs.push_back(1); // position of S in old basis
	CommonQNs.push_back(0); // position of Q in SingleSite
	CommonQNs.push_back(1); // position of S in SingleSite
	totSpos.push_back(1);   // SU(2) symmetry in position 1
	// HN params
	pHN->CalcHNMatEl=OneChQS_HN_MatEl;
	NumChannels=1;




	break;
	// end OneChQS Anderson model set up

      case 1: // Single-impurity Kondo model Q,S symmetry


	Nsites0=1;

	// fN (reduced)
	STLMatArray[0].NeedOld=false;
	STLMatArray[0].UpperTriangular=false;
	STLMatArray[0].CheckForMatEl=OneChQS_cd_check;
	STLMatArray[0].CalcMatEl=OneChQS_fN_MatEl;
	STLMatArray[0].SaveMatYN=false;

	STLMatArray.pop_back(); // 1 Matrix only
	STLMatArray.pop_back(); // 1 Matrix only
	STLMatArray.pop_back(); // 1 Matrix only

	switch (calcdens){
	case 0: // Levels only
	  NumThermoMats=0;
	  break;
	case 2: // Thermodynamics
	  // Z trick in file name.
	  strcpy(ThermoArray[0].ArqName,"SuscepImp1ChQS_Kondo");
	  strcat(ThermoArray[0].ArqName,zstring);
	  strcat(ThermoArray[0].ArqName,".dat");
	  strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QS");
	  strcat(ThermoArray[0].ChainArqName,zstring);
	  strcat(ThermoArray[0].ChainArqName,".dat");
	  strcpy(ThermoArray[1].ArqName,"EntropyImp1ChQS_Kondo");
	  strcat(ThermoArray[1].ArqName,zstring);
	  strcat(ThermoArray[1].ArqName,".dat");
	  strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QS");
	  strcat(ThermoArray[1].ChainArqName,zstring);
	  strcat(ThermoArray[1].ChainArqName,".dat");
	  // Values for the single-impurity Kondo model
	  ThermoArray[0].Calc=CalcSuscep;
	  ThermoArray[0].dImpValue=1.0/4.0; 
	  ThermoArray[1].Calc=CalcEntropy;
	  ThermoArray[1].dImpValue=log(2.0);
	  NumThermoMats=2;
	  break;
      
	default:
	  cout << " Calculating levels only " << endl;

	}
	// end switch calcdens

	// Param for H0
	Params.push_back(Lambda);
	Params.push_back(HalfLambdaFactor);
	Params.push_back(chi_m1); // will become Jtilde

	//"chi_m1"= sqrt(Fsq*rho_0*J_0/pi)/(sqrt(Lambda)*HalfLambdaFactor); 

	cout << " Fsq = " << chain.Fsq 
	     << " rho0 J0  = " << Lambda*HalfLambdaFactor*HalfLambdaFactor*chi_m1*chi_m1*pi/chain.Fsq << endl;


	OneChQS_SetKondoH0(Params,pAeig,pSingleSite,STLMatArray);

	// BuildBasis params
	CommonQNs.push_back(2); // Q and S and commont QNs
	CommonQNs.push_back(0); // position of Q in old basis
	CommonQNs.push_back(1); // position of S in old basis
	CommonQNs.push_back(0); // position of Q in SingleSite
	CommonQNs.push_back(1); // position of S in SingleSite
	totSpos.push_back(1);   // SU(2) symmetry in position 1
	// HN params
	pHN->CalcHNMatEl=OneChQS_HN_MatEl;
	NumChannels=1;
	// Set chi_0 to enter HN!
	chi_m1=chain.GetChin(0); // Should be EQUAL to chiN(0,Lambda) if Square
	cout << " chi_0= " << chi_m1 
	     << " Square band chi_0 : " << chiN(0,Lambda)  
	     << endl;

	break;
	// end OneChQS Kondo model set up


      case 3: // Chain only

	ZeroParams();
	Nsites0=1;
	NumThermoMats=2;

	// fN (reduced)
	STLMatArray[0].NeedOld=false;
	STLMatArray[0].UpperTriangular=false;
	STLMatArray[0].CheckForMatEl=OneChQS_cd_check;
	STLMatArray[0].CalcMatEl=OneChQS_fN_MatEl;
	STLMatArray.pop_back(); // Only ONE array
	STLMatArray.pop_back(); // Only ONE array
	STLMatArray.pop_back(); // Only ONE array
	// Set ChainH0: SingleSite is already set, just get Aeig and STLMatArray.
	OneChQS_SetChainH0(pAeig,STLMatArray);
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	SaveData=false;
	// BuildBasis params
	CommonQNs.push_back(2); // Q and S and commont QNs
	CommonQNs.push_back(0); // position of Q in old basis
	CommonQNs.push_back(1); // position of S in old basis
	CommonQNs.push_back(0); // position of Q in SingleSite
	CommonQNs.push_back(1); // position of S in SingleSite
	totSpos.push_back(1);   // SU(2) symmetry in position 1
	// HN params
	pHN->CalcHNMatEl=OneChQS_HN_MatEl;
	NumChannels=1;
	if (BandNo==0){
	  chi_m1=chiN(0,Lambda); // is actually chi0 // chiN is a codehandler func
	// Valid ONLY for SquareBand, z_twist=1! Check side dot calcs...
	}
	else{chi_m1=chain.GetChin(0);}
	// Need to change chain files for ALL Models!!
	//strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QS.dat");
	strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QS");
	strcat(ThermoArray[0].ChainArqName,zstring);
	strcat(ThermoArray[0].ChainArqName,".dat");
	//strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QS.dat");
	strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QS");
	strcat(ThermoArray[1].ChainArqName,zstring);
	strcat(ThermoArray[1].ChainArqName,".dat");
	break;

      case 6: // 1chQS Double Quantum dot

	cout << " Implementing DQD case... " << endl;

	chi2_m1=sqrt(2.0*dInitParams[4]/pi)/(sqrt(Lambda)*HalfLambdaFactor);

	Params.push_back(Lambda);
	Params.push_back(HalfLambdaFactor);
	Params.push_back(dInitParams[0]); // U1
	Params.push_back(dInitParams[2]); // ed1
	Params.push_back(chi_m1); // gamma1
	Params.push_back(dInitParams[3]); // U2
	Params.push_back(dInitParams[5]); // ed2
	Params.push_back(chi2_m1); // gamma2

	// Question: Can these things be implemented within the subroutine??

	// fN (reduced)
	STLMatArray[0].NeedOld=false;
	STLMatArray[0].UpperTriangular=false;
	STLMatArray[0].CheckForMatEl=OneChQS_cd_check;
	STLMatArray[0].CalcMatEl=OneChQS_fN_MatEl;


	switch (calcdens){
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
	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; // Need to do this after
	  // changes in STLMatArray!
	  SaveData=false;
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
	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; // Need to do this after
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"DQD");
	  break;
	case 2: // Thermodynamics
	  STLMatArray.pop_back(); // 1 Matrix only
	  STLMatArray.pop_back(); // 1 Matrix only
	  STLMatArray.pop_back(); // 1 Matrix only
	  //strcpy(ThermoArray[0].ArqName,"SuscepImp1ChQS_DQD.dat");
	  //strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QS.dat");
	  strcpy(ThermoArray[0].ArqName,"SuscepImp1ChQS_DQD");
	  strcat(ThermoArray[0].ArqName,zstring);
	  strcat(ThermoArray[0].ArqName,".dat");
	  strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QS");
	  strcat(ThermoArray[0].ChainArqName,zstring);
	  strcat(ThermoArray[0].ChainArqName,".dat");
	  //strcpy(ThermoArray[1].ArqName,"EntropyImp1ChQS_DQD.dat");
	  //strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QS.dat");
	  strcpy(ThermoArray[1].ArqName,"EntropyImp1ChQS_DQD");
	  strcat(ThermoArray[1].ArqName,zstring);
	  strcat(ThermoArray[1].ArqName,".dat");
	  strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QS");
	  strcat(ThermoArray[1].ChainArqName,zstring);
	  strcat(ThermoArray[1].ChainArqName,".dat");

	  ThermoArray[0].Calc=CalcSuscep;
	  ThermoArray[0].dImpValue=1.0/4.0; // check. Irrelevant, actually.
	  ThermoArray[1].Calc=CalcEntropy;
	  ThermoArray[1].dImpValue=4.0*log(2.0);
	  NumThermoMats=2;
	  break;
      
	default:
	  cout << " Calculating levels only " << endl;

	}
	// end switch calcdens


	OneChQS_SetH0_DQD(Params,pAeig,pSingleSite,STLMatArray);

	// BuildBasis params
	CommonQNs.push_back(2); // Q and S and commont QNs
	CommonQNs.push_back(0); // position of Q in old basis
	CommonQNs.push_back(1); // position of S in old basis
	CommonQNs.push_back(0); // position of Q in SingleSite
	CommonQNs.push_back(1); // position of S in SingleSite
	totSpos.push_back(1);   // SU(2) symmetry in position 1
	// HN params
	pHN->CalcHNMatEl=OneChQS_HN_MatEl;
	Nsites0=1;
	NumChannels=1;
	chi_m1=chiN(0,Lambda);
      
	break;
   
      default:
	cout << " Model not implemented for this symmetry. Exiting... " << endl;
	exit(0);
      }
      // end switch ModelNo
      break;

    case 1: // TwoChQS
      ////////////////////////////
      /////                   ////
      ///// Symmetry: TwoChQS ////
      /////                   ////
      ////////////////////////////
      /// Models so far:       ///
      ///  1 - Kondo           ///
      ///  3 - Chain           ///
      ////////////////////////////
      switch(ModelNo){
      case 1: 
	// Single-impurity Kondo model Q,S symmetry
	cout << " Implementing the two-Channel Kondo model " << endl;

	Nsites0=1;
	NumThermoMats=2;

	Params.clear();

	//Params.push_back(Lambda);
	//Params.push_back(HalfLambdaFactor); Do we need this??

	chi2_m1=sqrt(2.0*dInitParams[4]/pi)/(sqrt(Lambda)*HalfLambdaFactor);

	cout << " FOR NOW, not including the Lambda factor in J1,J2 " << endl;

// 	Params.push_back(chi_m1); // J1
// 	Params.push_back(chi2_m1); // J2

	Params.push_back(dInitParams[1]); // J1
	Params.push_back(dInitParams[4]); // J2


	switch (calcdens){
	case 0: // Levels only
	  NumThermoMats=0;
	  break;
	case 2:
	  // Z trick in file name.
	  strcpy(ThermoArray[0].ArqName,"SuscepImp2ChQS_Kondo");
	  strcat(ThermoArray[0].ArqName,zstring);
	  strcat(ThermoArray[0].ArqName,".dat");
	  strcpy(ThermoArray[0].ChainArqName,"SuscepChain2Ch_QS");
	  strcat(ThermoArray[0].ChainArqName,zstring);
	  strcat(ThermoArray[0].ChainArqName,".dat");
	  strcpy(ThermoArray[1].ArqName,"EntropyImp2ChQS_Kondo");
	  strcat(ThermoArray[1].ArqName,zstring);
	  strcat(ThermoArray[1].ArqName,".dat");
	  strcpy(ThermoArray[1].ChainArqName,"EntropyChain2Ch_QS");
	  strcat(ThermoArray[1].ChainArqName,zstring);
	  strcat(ThermoArray[1].ChainArqName,".dat");
	  // Values for the single-impurity Kondo model
	  ThermoArray[0].Calc=CalcSuscep;
	  ThermoArray[0].dImpValue=1.0/4.0; 
	  ThermoArray[1].Calc=CalcEntropy;
	  ThermoArray[1].dImpValue=log(2.0);
	  NumThermoMats=2;
	  break;
// 	case 3:
// 	  break;
	default:
	  cout << " Please select calcdens=0 or 2 for TwoChQS_Kondo " << endl;
	  exit(0);
	}
	// end switch calcdens

	STLMatArray.pop_back(); // 2 Matrices
	STLMatArray.pop_back(); // 2 Matrices

	// fN_ch1 (reduced)
	STLMatArray[0].NeedOld=false;
	STLMatArray[0].UpperTriangular=false; // WHY? Not symmetric
	STLMatArray[0].SaveMatYN=false;
	STLMatArray[0].CheckForMatEl=TwoChQS_cd_check;
	STLMatArray[0].CalcMatEl=TwoChQS_fNch1_MatEl;


	// fN_ch2 (reduced)
	STLMatArray[1].NeedOld=false;
	STLMatArray[1].UpperTriangular=false;  // WHY? Not symmetric
	STLMatArray[1].SaveMatYN=false;
	STLMatArray[1].CheckForMatEl=TwoChQS_cd_check;
	STLMatArray[1].CalcMatEl=TwoChQS_fNch2_MatEl;


	TwoChQS_SetKondoH0(Params,pSingleSite,pAeig,STLMatArray);
	// Do I need pAbasis?? Nope...

	// BuildBasis params
	CommonQNs.push_back(2); // Q and S and commont QNs
	CommonQNs.push_back(0); // position of Q in old basis
	CommonQNs.push_back(1); // position of S in old basis
	CommonQNs.push_back(0); // position of Q in SingleSite
	CommonQNs.push_back(1); // position of S in SingleSite
	totSpos.push_back(1);   // SU(2) symmetry in position 1
	// Hamiltonian/Code params
	// Need to check these things. (May 8/2011)
	pHN->CalcHNMatEl=TwoChQS_HN_MatEl;
	Nsites0=1;
	NumNRGmats=2;
	NumChannels=2;

	// Is this ok? What if there is a pseudogap??
	chi_m1=chiN(0,Lambda);
	chi2_m1=chiN(0,Lambda);
	if (BandNo==0){
	  chi_m1=chiN(0,Lambda); // is actually chi0 // chiN is a codehandler func
	  chi2_m1=chiN(0,Lambda);
	}
	// Valid ONLY for SquareBand, z_twist=1! Check side dot calcs...
	else{
	  chi_m1=chain.GetChin(0);
	  chi2_m1=chain.GetChin(0);
	}
	///////////


	// Debugging: Print and Exit
	pAeig->PrintEn();
	//exit(0);
	break;
      case 3: // TwoChQS Chain only

	ZeroParams();
	Nsites0=1;
	NumNRGmats=2;
	NumChannels=2;
	NumThermoMats=2;

	strcpy(ThermoArray[0].ChainArqName,"SuscepChain2Ch_QS");
	strcat(ThermoArray[0].ChainArqName,zstring);
	strcat(ThermoArray[0].ChainArqName,".dat");
	strcpy(ThermoArray[1].ChainArqName,"EntropyChain2Ch_QS");
	strcat(ThermoArray[1].ChainArqName,zstring);
	strcat(ThermoArray[1].ChainArqName,".dat");


	STLMatArray.pop_back(); // 2 Matrices
	STLMatArray.pop_back(); // 2 Matrices
	// fN_ch1 (reduced)
	STLMatArray[0].NeedOld=false;
	STLMatArray[0].UpperTriangular=false; // WHY? Not symmetric
	STLMatArray[0].SaveMatYN=false;
	STLMatArray[0].CheckForMatEl=TwoChQS_cd_check;
	STLMatArray[0].CalcMatEl=TwoChQS_fNch1_MatEl;
	// fN_ch2 (reduced)
	STLMatArray[1].NeedOld=false;
	STLMatArray[1].UpperTriangular=false;  // WHY? Not symmetric
	STLMatArray[1].SaveMatYN=false;
	STLMatArray[1].CheckForMatEl=TwoChQS_cd_check;
	STLMatArray[1].CalcMatEl=TwoChQS_fNch2_MatEl;

	TwoChQS_SetH0Chain(pAeig,pSingleSite,STLMatArray);
	// Debugging: Print and Exit
	pAeig->PrintEn();

	// BuildBasis params
	CommonQNs.push_back(2); // Q and S and commont QNs
	CommonQNs.push_back(0); // position of Q in old basis
	CommonQNs.push_back(1); // position of S in old basis
	CommonQNs.push_back(0); // position of Q in SingleSite
	CommonQNs.push_back(1); // position of S in SingleSite
	totSpos.push_back(1);   // SU(2) symmetry in position 1
	// Hamiltonian/Code params
	pHN->CalcHNMatEl=TwoChQS_HN_MatEl;
	if (BandNo==0){
	  chi_m1=chiN(0,Lambda); // is actually chi0 // chiN is a codehandler func
	  chi2_m1=chiN(0,Lambda);
	}
	// Valid ONLY for SquareBand, z_twist=1! Check side dot calcs...
	else{
	  chi_m1=chain.GetChin(0);
	  chi2_m1=chain.GetChin(0);
	}
	///////////
	break;
      default:
	cout << " Model not implemented for this symmetry. Exiting... " << endl;
	exit(0);
      }
      // end switch ModelNo
      break;


    case 2: // OneChQSz
      ////////////////////////////
      /////                    ////
      ///// Symmetry: OneChQSz ////
      /////                    ////
      //////////////////////////////////////
      /// Models so far:                 ///
      ///  0 - Anderson (Bfield)         ///
      ///  3 - Chain                     ///
      ///  6 - DQD (Zeeman in both dots) ///
      //////////////////////////////////////

      switch(ModelNo){
      case 0: // Anderson

	// Initialize matrices (OneChQ routines work here!)
	// Jul 09: using STMatArray

	STLMatArray[0].NeedOld=false;	  
	STLMatArray[0].CheckForMatEl=OneChQ_cd_check; //WHY???
	STLMatArray[0].CalcMatEl=OneChQ_fNup_MatEl;
	STLMatArray[0].SaveMatYN=false;

	STLMatArray[1].NeedOld=false;
	STLMatArray[1].CheckForMatEl=OneChQ_cd_check; //WHY???
	STLMatArray[1].CalcMatEl=OneChQ_fNdn_MatEl;
	STLMatArray[1].SaveMatYN=false;

	// Param for H0
	Params.push_back(dInitParams[0]); // U
	Params.push_back(dInitParams[2]); // ed
	Params.push_back(Lambda);
	Params.push_back(HalfLambdaFactor);
	Params.push_back(dInitParams[3]);        // Mag Field 

	switch (calcdens){
	case 0: // OpAvg: n_d(T), S_z(T) -(2013)
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false; // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;
	  STLMatArray[2]=auxNRGMat; // 2 - ndot (there already)
	  STLMatArray[3]=auxNRGMat; // 3 - Sz (there already)
	  strcpy(STLMatArray[2].MatName,"Ndot");
	  strcpy(STLMatArray[3].MatName,"Szdot");
	  // changes in STLMatArray!
	  SaveData=false;
	  // No Thermo
	  NumThermoMats=0;
	  break;
 	case 1: //Spectral density
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

	  strcpy(STLMatArray[2].MatName,"cdup");
	  strcpy(STLMatArray[3].MatName,"cddn");

	  cout << " N = -1" << endl;
	  cout << " cdN_up : " << endl;
	  STLMatArray[2].PrintAllBlocks();
	  cout << " cdN_dn : " << endl;
	  STLMatArray[3].PrintAllBlocks();

	  SaveData=true;

 	  break;
 	case 2: // Thermodynamics ONLY

	  STLMatArray.pop_back(); // 2 Matrices only
	  STLMatArray.pop_back(); // 2 Matrices only
	  strcpy(ThermoArray[0].ArqName,"SuscepImp1ChQSz");
	  strcat(ThermoArray[0].ArqName,zstring);
	  strcat(ThermoArray[0].ArqName,".dat");
	  strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QSz");
	  strcat(ThermoArray[0].ChainArqName,zstring);
	  strcat(ThermoArray[0].ChainArqName,".dat");

	  strcpy(ThermoArray[1].ArqName,"EntropyImp1ChQSz");
	  strcat(ThermoArray[1].ArqName,zstring);
	  strcat(ThermoArray[1].ArqName,".dat");
	  strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QSz");
	  strcat(ThermoArray[1].ChainArqName,zstring);
	  strcat(ThermoArray[1].ChainArqName,".dat");

	  ThermoArray[0].Calc=CalcSuscep;
	  ThermoArray[0].dImpValue=1.0/8.0; // check
	  ThermoArray[1].Calc=CalcEntropy;
	  ThermoArray[1].dImpValue=2.0*log(2.0);
	  NumThermoMats=2;
	  SaveData=false;

 	  break;

	case 4: // Dynamical spin susceptibility and <Sz> and <Sz^2> (2013)
	  // 6 Matrices
	  STLMatArray.push_back(auxNRGMat);
	  STLMatArray.push_back(auxNRGMat);
	  NumThermoMats=0;
	  // <Sz> and <Sz2> operator 
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false;
	  // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;

	  STLMatArray[2]=auxNRGMat; // 2 - Sz dynamical
	  STLMatArray[2].SaveMatYN=true;
	  STLMatArray[2].CalcAvg=false;


	  STLMatArray[3]=auxNRGMat; // 3 - Nd dynamical
	  STLMatArray[3].SaveMatYN=true;
	  STLMatArray[3].CalcAvg=false;

	  STLMatArray[4]=auxNRGMat; // 4 - <Sz2>(T) static
	  STLMatArray[4].SaveMatYN=false;

	  STLMatArray[5]=auxNRGMat; // 5 - <Nd2>(T) static
	  STLMatArray[5].SaveMatYN=false;

	  strcpy(STLMatArray[2].MatName,"Szomega");
	  strcpy(STLMatArray[3].MatName,"Ndomega");
	  strcpy(STLMatArray[4].MatName,"Sz2dot");
	  strcpy(STLMatArray[5].MatName,"Nd2dot");

	  NumNRGmats=STLMatArray.size();
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"1chQSzAnderson");
	  break;

 	case 5: //Spectral density with Self-energy trick
	  // 6 Matrices
	  STLMatArray.push_back(auxNRGMat);
	  STLMatArray.push_back(auxNRGMat);
	  // cd_up and cd_dn
	  STLMatArray[2].NeedOld=true;
	  STLMatArray[2].CheckForMatEl=OneChQSz_cdup_check;
	  STLMatArray[2].CalcMatEl=OneChQSz_cdup_MatEl;
	  STLMatArray[2].SaveMatYN=true;

	  STLMatArray[3].NeedOld=true;
	  STLMatArray[3].CheckForMatEl=OneChQSz_cddn_check;
	  STLMatArray[3].CalcMatEl=OneChQSz_cddn_MatEl;
	  STLMatArray[3].SaveMatYN=true;

	  // <cd_up N_dn> and <cd_dn N_up>
	  STLMatArray[4].NeedOld=true;
	  STLMatArray[4].CheckForMatEl=OneChQSz_cdup_check;
	  STLMatArray[4].CalcMatEl=OneChQSz_cdup_MatEl; // Same Fermi Factor
	  STLMatArray[4].SaveMatYN=true;

	  STLMatArray[5].NeedOld=true;
	  STLMatArray[5].CheckForMatEl=OneChQSz_cddn_check;
	  STLMatArray[5].CalcMatEl=OneChQSz_cddn_MatEl; // Same Fermi Factor
	  STLMatArray[5].SaveMatYN=true;

	  strcpy(STLMatArray[2].MatName,"cdup");
	  strcpy(STLMatArray[3].MatName,"cddn");
	  strcpy(STLMatArray[4].MatName,"cdupNdn");
	  strcpy(STLMatArray[5].MatName,"cddnNup");

	  SaveData=true;

 	  break;

	default:
	  cout << " calcdens = " << calcdens << " not implemented. Exiting. " << endl;
	  exit(0);

	}
	// end switch calcdens


	OneChQSz_SetAnderson_Hm1(Params,pAeig,STLMatArray);

	////////////////////////

	// BuildBasis params
	CommonQNs.push_back(2); // No of common QNs
	CommonQNs.push_back(0); // pos of QN 1 in old
	CommonQNs.push_back(1); // pos of QN 2 in old

	CommonQNs.push_back(0); // pos of QN 1 in SingleSite
	CommonQNs.push_back(1); // pos of QN 1 in SingleSite

	// No total S variables. Leave totSpos empty
	pHN->CalcHNMatEl=OneChQSz_HN_MatEl;

	// Wrap up
	NumChannels=1;
	Nsites0=0;
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	// changes in STLMatArray!


	break;

      case 3: // OneChQSz Chain only

	ZeroParams();
	Nsites0=1;
	NumThermoMats=2;
	NumChannels=1;

	// fN_up and fN_dn
	// fd_up
	STLMatArray[0].NeedOld=false;	  
	STLMatArray[0].UpperTriangular=false;
	STLMatArray[0].CheckForMatEl=OneChQSz_cdup_check; //This should work
	STLMatArray[0].CalcMatEl=OneChQSz_fNup_MatEl;
	STLMatArray[0].SaveMatYN=false;

	// fd_dn
	STLMatArray[1].NeedOld=false;
	STLMatArray[1].UpperTriangular=false;
	STLMatArray[1].CheckForMatEl=OneChQSz_cddn_check; //This should work
	STLMatArray[1].CalcMatEl=OneChQSz_fNdn_MatEl;
	STLMatArray[1].SaveMatYN=false;

	STLMatArray.pop_back(); // Only TWO mats
	STLMatArray.pop_back(); // Only TWO mats
	// Set ChainH0: SingleSite is already set, just get Aeig and STLMatArray.
	OneChQSz_SetChainH0(pAeig,STLMatArray);
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	SaveData=false;
	// BuildBasis params
	CommonQNs.push_back(2); // Q and Sz and commont QNs
	CommonQNs.push_back(0); // position of Q in old basis
	CommonQNs.push_back(1); // position of Sz in old basis
	CommonQNs.push_back(0); // position of Q in SingleSite
	CommonQNs.push_back(1); // position of Sz in SingleSite

	// No total S variables. Leave totSpos empty
	pHN->CalcHNMatEl=OneChQSz_HN_MatEl;

	if (BandNo==0){
	  chi_m1=chiN(0,Lambda); // is actually chi0 // chiN is a codehandler func
	// Valid ONLY for SquareBand, z_twist=1! Check side dot calcs...
	}
	else{chi_m1=chain.GetChin(0);}

	NumNRGmats=STLMatArray.size(); // DO I STILL NEED THIS???
	MatArray=&STLMatArray[0]; // Need to do this after
	// changes in STLMatArray! 

	// Need to change chain files for ALL Models!!
	strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QSz");
	strcat(ThermoArray[0].ChainArqName,zstring);
	strcat(ThermoArray[0].ChainArqName,".dat");
	strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QSz");
	strcat(ThermoArray[1].ChainArqName,zstring);
	strcat(ThermoArray[1].ChainArqName,".dat");
	break;
	// end OneChQSz chain only

      case 6: // OneChQSz DQD (Zeeman in both dots)

	// Initialize matrices (OneChQ routines work here!)
	// Jul 09: using STMatArray
	// fd_up
	STLMatArray[0].NeedOld=false;	  
	STLMatArray[0].CheckForMatEl=OneChQSz_cdup_check; //This should work
	STLMatArray[0].CalcMatEl=OneChQSz_fNup_MatEl;
	STLMatArray[0].SaveMatYN=false;

	// fd_dn
	STLMatArray[1].NeedOld=false;
	STLMatArray[1].CheckForMatEl=OneChQSz_cddn_check; //This should work
	STLMatArray[1].CalcMatEl=OneChQSz_fNdn_MatEl;
	STLMatArray[1].SaveMatYN=false;


	// Param for H0 (11 parameters)
	Params.push_back(Lambda);
	Params.push_back(HalfLambdaFactor);

	Params.push_back(dInitParams[0]); // U1
	Params.push_back(dInitParams[2]); // ed1
	Params.push_back(chi_m1); // gamma1
	Params.push_back(dInitParams[7]); // Bmag1      

	Params.push_back(dInitParams[3]); // U2
	Params.push_back(dInitParams[5]); // ed2
	chi2_m1=sqrt(2.0*dInitParams[4]/pi)/(sqrt(Lambda)*HalfLambdaFactor);
	Params.push_back(chi2_m1); // gamma2
	Params.push_back(dInitParams[8]); // Bmag2
	  
	Params.push_back(dInitParams[6]); // lambda


	switch (calcdens){
	case 1: // Calculates spectral densities
	  // cd_up and cd_dn (OneChQSz_SetAnderson equals 
	  // c1_up
	  STLMatArray[2].NeedOld=true;
	  STLMatArray[2].CheckForMatEl=OneChQSz_cdup_check;
	  STLMatArray[2].CalcMatEl=OneChQSz_cdup_MatEl;
	  STLMatArray[2].SaveMatYN=true;
	  strcpy(STLMatArray[2].MatName,"cdot1_up");
	  // c1_dn
	  STLMatArray[3].NeedOld=true;
	  STLMatArray[3].CheckForMatEl=OneChQSz_cddn_check;
	  STLMatArray[3].CalcMatEl=OneChQSz_cddn_MatEl;
	  STLMatArray[3].SaveMatYN=true;
	  strcpy(STLMatArray[3].MatName,"cdot1_dn");
	  // c2_up
	  auxNRGMat=STLMatArray[2];
	  STLMatArray.push_back(auxNRGMat);
	  strcpy(STLMatArray[4].MatName,"cdot2_up");
	  // c2_dn
	  auxNRGMat=STLMatArray[3];
	  STLMatArray.push_back(auxNRGMat);
	  strcpy(STLMatArray[5].MatName,"cdot2_dn");
	  // changes in STLMatArray
	  SaveData=true;
	  strcpy(SaveArraysFileName,"DQD_QSz");

	  break;
	case 2: // Thermodynamics for DQD+Zeeman
	  STLMatArray.pop_back(); // 2 Matrices only
	  STLMatArray.pop_back(); // 2 Matrices only
	  //strcpy(ThermoArray[0].ArqName,"SuscepImp1ChQSz_DQD.dat");
	  strcpy(ThermoArray[0].ArqName,"SuscepImp1ChQSz_DQD");
	  strcat(ThermoArray[0].ArqName,zstring);
	  strcat(ThermoArray[0].ArqName,".dat");
	  strcpy(ThermoArray[0].ChainArqName,"SuscepChain1Ch_QSz");
	  strcat(ThermoArray[0].ChainArqName,zstring);
	  strcat(ThermoArray[0].ChainArqName,".dat");
	  //strcpy(ThermoArray[1].ArqName,"EntropyImp1ChQSz_DQD.dat");
	  strcpy(ThermoArray[1].ArqName,"EntropyImp1ChQSz_DQD");
	  strcat(ThermoArray[1].ArqName,zstring);
	  strcat(ThermoArray[1].ArqName,".dat");
	  strcpy(ThermoArray[1].ChainArqName,"EntropyChain1Ch_QSz");
	  strcat(ThermoArray[1].ChainArqName,zstring);
	  strcat(ThermoArray[1].ChainArqName,".dat");

	  ThermoArray[0].Calc=CalcSuscep;
	  ThermoArray[0].dImpValue=1.0/4.0; // check
	  ThermoArray[1].Calc=CalcEntropy;
	  ThermoArray[1].dImpValue=4.0*log(2.0);
	  NumThermoMats=2;
	  SaveData=false;

	  break;

	default:
	  cout << " Please select calcdens=1,2 or 3 for 1chQSz_DQD " << endl;
	  exit(0);
	}
	// end switch calcdens

	OneChQSz_SetH0_DQD(Params,pAeig,pSingleSite,STLMatArray);

	// Wrap up
	NumChannels=1;
	Nsites0=1;
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	if (BandNo==0)
	  chi_m1=chiN(0,Lambda); // is actually chi0
	// Valid ONLY for SquareBand, z_twist=1! Check side dot calcs...
	else
	  chi_m1=chain.GetChin(0);
	///////////
	// BuildBasis params
	CommonQNs.push_back(2); // No of common QNs
	CommonQNs.push_back(0); // pos of QN 1 in old
	CommonQNs.push_back(1); // pos of QN 2 in old
	CommonQNs.push_back(0); // pos of QN 1 in SingleSite
	CommonQNs.push_back(1); // pos of QN 1 in SingleSite
	// No total S variables. Leave totSpos empty
	pHN->CalcHNMatEl=OneChQSz_HN_MatEl;

	break;
	// end DQD 1chQSz

      default:
	cout << " Model not implemented for this symmetry. Exiting... " << endl;
	exit(0);
      }
      // end switch ModelNo
      break;
    case 3: // OneChQ
      ///// Symmetry: OneChQ ////
      ////////////////////////////
      /////                  ////
      ///// Symmetry: OneChQ ////
      /////                  ////
      ////////////////////////////
      /// Models so far:       ///
      ///  0 - Anderson        ///
      ///  5 - SMM             ///
      ////////////////////////////

      switch(ModelNo){
      case 0: // Anderson

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

	if (calcdens==1){
	  // All operators have the same Checks/MatEls in |Q> basis
	  // what changes is the initial set-up
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false;
	  auxNRGMat.CheckForMatEl=OneChQ_cd_check;
	  auxNRGMat.CalcMatEl=OneChQ_cd_MatEl;
	  auxNRGMat.SaveMatYN=true;

	  STLMatArray.push_back(auxNRGMat); // 4 - cd1_up
	  STLMatArray.push_back(auxNRGMat); // 5 - cd1_dn
	  //NumNRGarrays=STLMatArray.size(); // do we need NumNRGarrays later?
	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; //
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"1chQAnd");
	}
	else{
	  // Set Thermodynamics
	  // Entropy calculation only
	  strcpy(ThermoArray[0].ArqName,"EntropyImp1Ch_Anderson_Q.dat");
	  strcpy(ThermoArray[0].ChainArqName,"EntropyChain1Ch_Q.dat");
	  ThermoArray[0].Calc=CalcEntropy;
	  ThermoArray[0].dImpValue=2.0*log(2.0);
	  //ThermoArray[0].CalcChain=false;
	  NumThermoMats=1;
	}
	// calculate spec function or thermo?


	// Param for H0
	Params.push_back(dInitParams[0]); // U
	Params.push_back(dInitParams[2]); // ed
	Params.push_back(Lambda);
	Params.push_back(HalfLambdaFactor);
	OneChQ_SetAnderson_Hm1(Params,pAeig,STLMatArray);
	//OneChQ_SetAnderson_Hm1_old(Params,&Aeig,MatArray);


	////////////////////////

	// BuildBasis params
	CommonQNs.push_back(1); // No of common QNs
	CommonQNs.push_back(0); // pos of QN 1 in old
	CommonQNs.push_back(0); // pos of QN 1 in SingleSite
	// No total S variables. Leave totSpos empty
	pHN->CalcHNMatEl=OneChQ_HN_MatEl;

	NumChannels=1;
	Nsites0=0;

	break;
      case 5: // SMM model

	// Params for H0
	Params.push_back(Lambda);  // Lambda
	Params.push_back(HalfLambdaFactor); // Lambda factor


	if (dInitParams.size()>10){
	  Params.push_back(dInitParams[0]); // U1
	  Params.push_back(dInitParams[2]); // ed1

	  Params.push_back(dInitParams[3]); // U2
	  Params.push_back(dInitParams[5]); // ed2
 
	  Params.push_back(dInitParams[6]); // J12

	  Params.push_back(dInitParams[7]); // BmagPar      
	  Params.push_back(dInitParams[8]); // BmagPerp

	  Params.push_back(dInitParams[9]); // Dz      
	  Params.push_back(dInitParams[10]); // B2 anisot

	  double chi2_m1=sqrt(2.0*dInitParams[4]/pi)/(sqrt(Lambda)*HalfLambdaFactor);
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
	if (calcdens==1){
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
	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; // Need to do this after
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"SMM");
      
	}
	else{
	  // Set Thermodynamics
	
	  // Entropy calculation only
      
	  strcpy(ThermoArray[0].ArqName,"EntropyImp1Ch_Anderson_Q.dat");
	  strcpy(ThermoArray[0].ChainArqName,"EntropyChain1Ch_Q.dat");
	  ThermoArray[0].Calc=CalcEntropy;
	  ThermoArray[0].dImpValue=4.0*log(2.0); // Seems to be irrelevant...
	  ThermoArray[0].CalcChain=false;
	  NumThermoMats=1;
	}
	// calculate spec function or thermo?


	// Set H0: output pAeig, MatArray



	OneChQ_SetSMM_H0(Params,pAeig,pSingleSite,STLMatArray);


	// BuildBasis params
	CommonQNs.push_back(1); // No of common QNs
	CommonQNs.push_back(0); // pos of QN 1 in old
	CommonQNs.push_back(0); // pos of QN 1 in SingleSite
	// No total S variables. Leave totSpos empty
	pHN->CalcHNMatEl=OneChQ_HN_MatEl;
	Nsites0=1;
	NumChannels=1;

	chi_m1=chiN(0,Lambda);

	break;
      default:
	cout << " Model not implemented for this symmetry. Exiting... " << endl;
	exit(0);
      }
      // end switch ModelNo
      break;

    case 4: // TwoChQSP
      switch(ModelNo){
      ////////////////////////////
      /////                     ///
      ///// Symmetry: TwoChQSP  ///
      /////                     ///
      ////////////////////////////
      /// Models so far:       ///
      ///  4 - CM Phonons (2ch)///
      ////////////////////////////
      case 4: // CM Phonons - 2channel
	Params.push_back(Lambda);  // Lambda
	Params.push_back(HalfLambdaFactor); // Lambda factor

	Params.push_back(dInitParams[0]); // U1
	Params.push_back(dInitParams[2]); // ed1

	Params.push_back(chi_m1); // sqrt(2gamma1/Pi)/sqrt(L)*HalfLambdaFactor

	Params.push_back(dInitParams[3]); // w0
	Params.push_back(dInitParams[4]); // lambda_ph
	Params.push_back(dInitParams[5]); // alpha
	Params.push_back(dInitParams[6]); // Nph      


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
	TwoChQSP_SetH0CMphonon(Params, pAeig,pSingleSite,MatArray);


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
	pHN->CalcHNMatEl=TwoChQS_HN_MatEl;
	Nsites0=1;
	chi_m1=chiN(0,Lambda);
	NumNRGmats=2;
	NumChannels=2;
	NumThermoMats=2;

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
      default:
	cout << " Model not implemented for this symmetry. Exiting... " << endl;
	exit(0);
      }  
      // end switch ModelNo
      break;

    case 6: // OneChNupPdn
      ////////////////////////////////
      /////                        ///
      ///// Symmetry: OneChNupPdn  ///
      /////                        ///
      ////////////////////////////////

      // BuildBasis params
      CommonQNs.push_back(2); // No of common QNs
      CommonQNs.push_back(0); // pos of QN 1 in old
      CommonQNs.push_back(1); // pos of QN 2 in old

      CommonQNs.push_back(0); // pos of QN 1 in SingleSite
      CommonQNs.push_back(1); // pos of QN 2 in SingleSite

      CommonQNs.push_back(1); // pos of Parity QN
      // No total S variables. Leave totSpos empty

      switch(ModelNo){
      ///////////////////////////////////
      /// Models so far:              ///
      ///  3 - Chain only             ///
      ///  7 - Anderson+Local Majorana///
      ///////////////////////////////////
      case 3: // Chain only

	ZeroParams();
	Nsites0=1; // IS THIS CORRECT? Yes. H_0 is set so we start NRG from H_1
	NumThermoMats=1; // Entropy only

	// Initialize matrices (real matrices)
	// fd_up
	STLMatArray[0].NeedOld=false;	  
 	STLMatArray[0].CheckForMatEl=OneChNupPdn_cdup_check;
  	STLMatArray[0].CalcMatEl=OneChNupPdn_fNup_MatEl;
 	//STLMatArray[0].CalcMatElCplx=OneChNupPdn_fNup_MatElCplx;
	STLMatArray[0].SaveMatYN=false;
	STLMatArray[0].IsComplex=false;

	// fd_dn
	STLMatArray[1].NeedOld=false;
 	STLMatArray[1].CheckForMatEl=OneChNupPdn_cddn_check;
  	STLMatArray[1].CalcMatEl=OneChNupPdn_fNdn_MatEl;
	STLMatArray[1].SaveMatYN=false;
	STLMatArray[1].IsComplex=false;

	STLMatArray.pop_back(); // Only TWO arrays
	STLMatArray.pop_back(); // Only TWO arrays

	OneChNupPdn_SetChainH0(pAeig,STLMatArray);

	// Set pHN
	pHN->CalcHNMatEl=OneChNupPdn_HN_MatEl; // REAL
	// Wrap up
	NumChannels=1;
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	// changes in STLMatArray!

	if (BandNo==0){
	  chi_m1=chiN(0,Lambda); // is actually chi0 // chiN is a codehandler func
	// Valid ONLY for SquareBand, z_twist=1! Check side dot calcs...
	}
	else{chi_m1=chain.GetChin(0);}
	// Need to change chain files for ALL Models!!
	strcpy(ThermoArray[0].ChainArqName,"EntropyChain1Ch_NupPdn");
	strcat(ThermoArray[0].ChainArqName,zstring);
	strcat(ThermoArray[0].ChainArqName,".dat");
	ThermoArray[0].Calc=CalcEntropy;
	ThermoArray[0].CalcChain=true; // Just to make sure

	break;

      case 7: 
	cout << " Majorana + QD effective model " << endl;

	Params.push_back(Lambda);  // Lambda
	Params.push_back(HalfLambdaFactor); // Lambda factor

	Params.push_back(dInitParams[0]); // U1
	Params.push_back(dInitParams[2]); // ed1
	Params.push_back(chi_m1); // sqrt(2gamma1/Pi)/sqrt(L)*HalfLambdaFactor

	Params.push_back(dInitParams[3]); // Mag Field 

	Params.push_back(dInitParams[4]); // t1
	Params.push_back(dInitParams[5]); // t2
	Params.push_back(dInitParams[6]); // phi_mag
	Params.push_back(dInitParams[7]); // em      

	// Initialize matrices 
	// fd_up
	STLMatArray[0].NeedOld=false;	  
 	STLMatArray[0].CheckForMatEl=OneChNupPdn_cdup_check;
//  	STLMatArray[0].CalcMatEl=OneChNupPdn_fNup_MatEl;
 	STLMatArray[0].CalcMatElCplx=OneChNupPdn_fNup_MatElCplx;
	STLMatArray[0].SaveMatYN=false;
	STLMatArray[0].IsComplex=true;

	// fd_dn
	STLMatArray[1].NeedOld=false;
 	STLMatArray[1].CheckForMatEl=OneChNupPdn_cddn_check;
 	STLMatArray[1].CalcMatElCplx=OneChNupPdn_fNdn_MatElCplx;
	STLMatArray[1].SaveMatYN=false;
	STLMatArray[1].IsComplex=true;


	switch (calcdens){
	case 0: // OpAvg: n_d(T), S_z(T) -(2014)
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false; // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatElCplx=ImpOnly_MatElCplx;
	  auxNRGMat.IsComplex=true;

	  STLMatArray[2]=auxNRGMat; // 2 - ndot (there already)
	  STLMatArray[3]=auxNRGMat; // 3 - Sz (there already)
	  strcpy(STLMatArray[2].MatName,"Ndot");
	  strcpy(STLMatArray[3].MatName,"Szdot");
	  STLMatArray.push_back(auxNRGMat);
	  strcpy(STLMatArray[4].MatName,"Nf"); // 4 -Nf
	  // changes in STLMatArray!
	  SaveData=false;
	  // No Thermo
	  NumThermoMats=0;
	  break;
 	case 1: //Spectral density
	  // cd_up and cd_dn (OneChQSz_SetAnderson equals 
	  // MatArray[2,3] to MatArray[0,1]
	  STLMatArray[2].NeedOld=true;
	  STLMatArray[2].CheckForMatEl=OneChNupPdn_cdup_check;
	  STLMatArray[2].CalcMatElCplx=OneChNupPdn_cdup_MatElCplx;
	  STLMatArray[2].SaveMatYN=true;
	  STLMatArray[2].IsComplex=true;
	  strcpy(STLMatArray[2].MatName,"cd_up");

	  STLMatArray[3].NeedOld=true;
	  STLMatArray[3].CheckForMatEl=OneChNupPdn_cddn_check;
	  STLMatArray[3].CalcMatElCplx=OneChNupPdn_cddn_MatElCplx;
	  STLMatArray[3].SaveMatYN=true;
	  STLMatArray[3].IsComplex=true;
	  strcpy(STLMatArray[3].MatName,"cd_dn");

	  // f_Maj (this looks ok)
	  auxNRGMat=STLMatArray[3];
	  STLMatArray.push_back(auxNRGMat);
	  strcpy(STLMatArray[4].MatName,"f_Maj");

	  // TODO: Add other spectral functions??

	  cout << " N = -1: Nothing here" << endl;
// 	  cout << " fdN_up : " << endl;
// 	  STLMatArray[0].PrintAllBlocks();
// 	  cout << " fdN_dn : " << endl;
// 	  STLMatArray[1].PrintAllBlocks();
// 	  cout << " cdN_up : " << endl;
// 	  STLMatArray[2].PrintAllBlocks();
// 	  cout << " cdN_dn : " << endl;
// 	  STLMatArray[3].PrintAllBlocks();
// 	  cout << " f_Maj : " << endl;
// 	  STLMatArray[4].PrintAllBlocks();

	  SaveData=true;

 	  break;
	case 2:
	  cout << "Entropy calculation (implementing it...)" << endl;

	  STLMatArray.pop_back(); // Only TWO arrays
	  STLMatArray.pop_back(); // Only TWO arrays

	  // Need to change chain files for ALL Models!!
	  strcpy(ThermoArray[0].ArqName,"EntropyImp1ChNupPdn_Majorana");
	  strcat(ThermoArray[0].ArqName,zstring);
	  strcat(ThermoArray[0].ArqName,".dat");
	  strcpy(ThermoArray[0].ChainArqName,"EntropyChain1Ch_NupPdn");
	  strcat(ThermoArray[0].ChainArqName,zstring);
	  strcat(ThermoArray[0].ChainArqName,".dat");
	  ThermoArray[0].Calc=CalcEntropy;
	  ThermoArray[0].dImpValue=3.0*log(2.0); // log(8)
	  ThermoArray[0].CalcChain=false;
	  NumThermoMats=1;
	  SaveData=false;

	  break;
	default:
	  cout << " calcdens = " << calcdens << " not implemented. Exiting. " << endl;
	  exit(0);
	}
	// end switch calcdens

	// actually SetHm1 !
	OneChNupPdn_SetH0_AndersonMajorana(Params,pSingleSite,pAeig,STLMatArray);

// 	// BuildBasis params (already out of here)
	// Set pHN
	pHN->IsComplex=true;
	pHN->CalcHNMatElCplx=OneChNupPdn_HN_MatElCplx;

	// Wrap up
	NumChannels=1;
	Nsites0=0;
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	// changes in STLMatArray!

	cout << " N=-1 ok. Going to N=0... " << endl;
	
	cout << " pAeig is complex ? " << pAeig->CheckComplex() << endl;


	//exit(0);
	break;
      default:
	cout << " Model not implemented for this symmetry. Exiting... " << endl;
	exit(0);
      }  
      // end switch ModelNo for OneChNupPdn
      break;

    case 7: // OneChS
      /////////////////////////////
      /////                    ////
      ///// Symmetry: OneChS   ////
      /////                    ////
      /////////////////////////////
      // BuildBasis params
      CommonQNs.push_back(1); // No of common QNs
      CommonQNs.push_back(0); // pos of QN 1 in old
      CommonQNs.push_back(0); // pos of QN 1 in SingleSite
      
      totSpos.push_back(0);   // SU(2) symmetry in position 1

      switch(ModelNo){
      //////////////////////////////////////
      /// Models so far:                 ///
      ///  0 - Anderson (SC leads)       ///
      //////////////////////////////////////
      case 0: // Single-impurity Anderson model S symmetry with SC leads

	// Initialize matrices 
	// Jul 09: using STMatArray

	// fN (reduced)
	STLMatArray[0].NeedOld=false;
	STLMatArray[0].UpperTriangular=false;
	STLMatArray[0].CheckForMatEl=OneChS_cd_check;
	STLMatArray[0].CalcMatEl=OneChS_fN_MatEl;
	STLMatArray[0].SaveMatYN=false;

	switch (calcdens){
	case 0: // OpAvg: n_d(T), n^2_d(T), S^2(T) - NEW (2012)
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false; // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;
	  STLMatArray[1]=auxNRGMat; // 1 - ndot (there already)
	  STLMatArray[2]=auxNRGMat; // 2 - ndot^2 (there already)
	  STLMatArray[3]=auxNRGMat; // 3 - Sdot^2 (there already)
	  strcpy(STLMatArray[1].MatName,"Ndot");
	  strcpy(STLMatArray[2].MatName,"NdSq");
	  strcpy(STLMatArray[3].MatName,"SdSq");
	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; // Need to do this after
	  // changes in STLMatArray!
	  SaveData=false;
// 	  STLMatArray.pop_back(); // 1 Matrices only
// 	  STLMatArray.pop_back(); // 1 Matrices only
// 	  STLMatArray.pop_back(); // 1 Matrices only
	  // No Thermo
	  NumThermoMats=0;
	  break;
	case 1: //Spectral density
	  STLMatArray.pop_back(); // 2 Matrices
	  STLMatArray.pop_back(); // 2 Matrices
	  NumThermoMats=0;
	  // cd operator 
	  STLMatArray[1].NeedOld=true;
	  STLMatArray[1].UpperTriangular=false;
	  STLMatArray[1].CheckForMatEl=OneChS_cd_check;
	  STLMatArray[1].CalcMatEl=OneChS_cd_MatEl;
	  STLMatArray[1].SaveMatYN=true;
	  strcpy(STLMatArray[1].MatName,"cdot1");
	  NumNRGmats=STLMatArray.size();
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"1chSAndersonSC");
	  break;
	case 4: // Dynamical spin susceptibility and <Sz> and <Sz^2> (2015)
	  // 6 Matrices
	  STLMatArray.push_back(auxNRGMat);
	  STLMatArray.push_back(auxNRGMat);
	  NumThermoMats=0;
	  // <Sz> and <Sz2> operator 
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false;
	  // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;

	  STLMatArray[2]=auxNRGMat; // 2 - Sz dynamical
	  STLMatArray[2].SaveMatYN=true;
	  STLMatArray[2].CalcAvg=false;

	  STLMatArray[3]=auxNRGMat; // 3 - Nd dynamical
	  STLMatArray[3].SaveMatYN=true;
	  STLMatArray[3].CalcAvg=false;

	  STLMatArray[4]=auxNRGMat; // 4 - <Sz2>(T) static
	  STLMatArray[4].SaveMatYN=false;

	  STLMatArray[5]=auxNRGMat; // 5 - <Nd2>(T) static
	  STLMatArray[5].SaveMatYN=false;

	  strcpy(STLMatArray[2].MatName,"Szomega");
	  strcpy(STLMatArray[3].MatName,"Ndomega");
	  strcpy(STLMatArray[4].MatName,"Sz2dot");
	  strcpy(STLMatArray[5].MatName,"Nd2dot");

	  NumNRGmats=STLMatArray.size();
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"1chSAndersonSC");
	  break;

	default:
	  cout << " calcdens = " << calcdens << " not implemented. Exiting. " << endl;
	  exit(0);

	}
	// end switch calcdens

	// Param for H0
	Params.push_back(Lambda);
	Params.push_back(HalfLambdaFactor);
	Params.push_back(dInitParams[0]); // U
	Params.push_back(dInitParams[2]); // ed

	OneChS_SetAnderson_Hm1(Params,pAeig,STLMatArray);

	// HN params
	pHN->CalcHNMatEl=OneChS_HNsc_MatEl;

	// Wrap up
	NumChannels=1;
	Nsites0=0;
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	// changes in STLMatArray!

	break;
	// end OneChS Anderson model set up

      default:
	cout << " Model not implemented for OneChS symmetry. Exiting... " << endl;
	exit(0);
      }  
      // end switch ModelNo for OneChS
      break;

      ///////////////////////

    case 8: // OneChSz
      /////////////////////////////
      /////                    ////
      ///// Symmetry: OneChSz   ////
      /////                    ////
      /////////////////////////////
      // BuildBasis params
      CommonQNs.push_back(1); // No of common QNs
      CommonQNs.push_back(0); // pos of QN 1 in old
      CommonQNs.push_back(0); // pos of QN 1 in SingleSite
      
      // No total S variables. Leave totSpos empty

      switch(ModelNo){
      //////////////////////////////////////
      /// Models so far:                 ///
      ///  0 - Anderson (SC leads)       ///
      //////////////////////////////////////
      case 0: // Single-impurity Anderson model S symmetry with SC leads

	// Initialize matrices 
	// Jul 09: using STMatArray

	// fNup and fNdn
	STLMatArray[0].NeedOld=false;
	STLMatArray[0].UpperTriangular=false;
// 	STLMatArray[0].CheckForMatEl=OneChS_cd_check; // Should work
	STLMatArray[0].CheckForMatEl=OneChSz_cdup_check; // More efficient
	STLMatArray[0].CalcMatEl=OneChSz_fNup_MatEl;
	STLMatArray[0].SaveMatYN=false;

	STLMatArray[1].NeedOld=false;
	STLMatArray[1].UpperTriangular=false;
// 	STLMatArray[1].CheckForMatEl=OneChS_cd_check; // Should work
	STLMatArray[1].CheckForMatEl=OneChSz_cddn_check; // Should work
	STLMatArray[1].CalcMatEl=OneChSz_fNdn_MatEl;
	STLMatArray[1].SaveMatYN=false;



	switch (calcdens){
	case 0: // OpAvg: n_d(T), Sz(T), Sz(T) - NEW (2012)
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false; // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;
	  STLMatArray.push_back(auxNRGMat); // 5 matrices
	  STLMatArray[2]=auxNRGMat; // 1 - ndot (there already)
	  STLMatArray[3]=auxNRGMat; // 2 - ndot^2 (there already)
	  STLMatArray[4]=auxNRGMat; // 3 - Sdot^2 (there already)
	  strcpy(STLMatArray[2].MatName,"Ndot");
	  strcpy(STLMatArray[3].MatName,"NdSq");
	  strcpy(STLMatArray[4].MatName,"SdSq");
	  NumNRGmats=STLMatArray.size();
	  MatArray=&STLMatArray[0]; // Need to do this after
	  // changes in STLMatArray!
	  SaveData=false;
	  // No Thermo
	  NumThermoMats=0;
	  break;
	case 1: //Spectral density - 4 matrices (no pop-up)
	  NumThermoMats=0;
	  // cdup operator 
	  STLMatArray[2].NeedOld=true;
	  STLMatArray[2].UpperTriangular=false;
// 	  STLMatArray[2].CheckForMatEl=OneChS_cd_check;
	  STLMatArray[2].CheckForMatEl=OneChSz_cdup_check;
	  STLMatArray[2].CalcMatEl=OneChSz_cdup_MatEl;
	  STLMatArray[2].SaveMatYN=true;
	  strcpy(STLMatArray[2].MatName,"cdotup");
	  // cddn operator 
	  STLMatArray[3].NeedOld=true;
	  STLMatArray[3].UpperTriangular=false;
// 	  STLMatArray[3].CheckForMatEl=OneChS_cd_check;
	  STLMatArray[3].CheckForMatEl=OneChSz_cddn_check;
	  STLMatArray[3].CalcMatEl=OneChSz_cddn_MatEl;
	  STLMatArray[3].SaveMatYN=true;
	  strcpy(STLMatArray[3].MatName,"cdotdn");


	  NumNRGmats=STLMatArray.size();
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"1chSzAndersonSC");
	  break;
	case 4: // Dynamical spin susceptibility and <Sz> and <Sz^2> (2015)
	  /// Stopped here
	  // 6 Matrices
	  STLMatArray.push_back(auxNRGMat);
	  STLMatArray.push_back(auxNRGMat);
	  NumThermoMats=0;
	  // <Sz> and <Sz2> operator 
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false;
	  // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatEl=ImpOnly_MatEl;

	  STLMatArray[2]=auxNRGMat; // 2 - Sz dynamical
	  STLMatArray[2].SaveMatYN=true;
	  STLMatArray[2].CalcAvg=false;

	  STLMatArray[3]=auxNRGMat; // 3 - Nd dynamical
	  STLMatArray[3].SaveMatYN=true;
	  STLMatArray[3].CalcAvg=false;

	  STLMatArray[4]=auxNRGMat; // 4 - <Sz2>(T) static
	  STLMatArray[4].SaveMatYN=false;

	  STLMatArray[5]=auxNRGMat; // 5 - <Nd2>(T) static
	  STLMatArray[5].SaveMatYN=false;

	  strcpy(STLMatArray[2].MatName,"Szomega");
	  strcpy(STLMatArray[3].MatName,"Ndomega");
	  strcpy(STLMatArray[4].MatName,"Sz2dot");
	  strcpy(STLMatArray[5].MatName,"Nd2dot");

	  NumNRGmats=STLMatArray.size();
	  // changes in STLMatArray!
	  SaveData=true;
	  strcpy(SaveArraysFileName,"1chSzAndersonSC");
	  break;

	default:
	  cout << " calcdens = " << calcdens << " not implemented. Exiting. " << endl;
	  exit(0);

	}
	// end switch calcdens

	// Param for H0
	Params.push_back(Lambda);
	Params.push_back(HalfLambdaFactor);
	Params.push_back(dInitParams[0]); // U
	Params.push_back(dInitParams[2]); // ed
	Params.push_back(dInitParams[4]); // Mag Field.
       //Note that dInitParams[3] is DeltaSC 

	OneChSz_SetAnderson_Hm1(Params,pAeig,STLMatArray);

	// HN params
	pHN->CalcHNMatEl=OneChSz_HNsc_MatEl;

	NumChannels=1;

	// Wrap up
	NumChannels=1;
	Nsites0=0;
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	// changes in STLMatArray!

	//exit(0);
	break;
	// end OneChSz Anderson model set up


      default:
	cout << " Model not implemented for OneChSz symmetry. Exiting... " << endl;
	exit(0);
      }  
      // end switch ModelNo for OneChSz
      break;
    case 9: // OneChPuPdn
      //////////////////////////////////
      /////                         ////
      ///// Symmetry: OneChPupPdn   ////
      /////                         ////
      //////////////////////////////////
      // BuildBasis params
      CommonQNs.push_back(2); // No of common QNs
      CommonQNs.push_back(0); // pos of QN 1 in old
      CommonQNs.push_back(1); // pos of QN 2 in old

      CommonQNs.push_back(0); // pos of QN 1 in SingleSite
      CommonQNs.push_back(1); // pos of QN 2 in SingleSite

      CommonQNs.push_back(0); // pos of Parity QN
      CommonQNs.push_back(1); // pos of Parity QN

      // No total S variables. Leave totSpos empty

      switch(ModelNo){
      ///////////////////////////////////////
      /// Models so far:                  ///
      ///  3 - Chain only                 ///
      ///  7 - Anderson+ 2 Local Majoranas///
      ///////////////////////////////////////
      case 3: // Chain only

	ZeroParams();
	Nsites0=1; // IS THIS CORRECT? Yes. H_0 is set so we start NRG from H_1
	NumThermoMats=1; // Entropy only

	// Initialize matrices (real matrices)
	// fd_up
	STLMatArray[0].NeedOld=false;	  
 	STLMatArray[0].CheckForMatEl=OneChPupPdn_cdup_check;
  	STLMatArray[0].CalcMatEl=OneChPupPdn_fNup_MatEl;
 	//STLMatArray[0].CalcMatElCplx=OneChNupPdn_fNup_MatElCplx;
	STLMatArray[0].SaveMatYN=false;
	STLMatArray[0].IsComplex=false;

	// fd_dn
	STLMatArray[1].NeedOld=false;
 	STLMatArray[1].CheckForMatEl=OneChPupPdn_cddn_check;
  	STLMatArray[1].CalcMatEl=OneChPupPdn_fNdn_MatEl;
	STLMatArray[1].SaveMatYN=false;
	STLMatArray[1].IsComplex=false;

	STLMatArray.pop_back(); // Only TWO arrays
	STLMatArray.pop_back(); // Only TWO arrays

	OneChPupPdn_SetChainH0(pAeig,STLMatArray);

	// Set pHN
	pHN->CalcHNMatEl=OneChPupPdn_HN_MatEl; // REAL
	// Wrap up
	NumChannels=1;
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	// changes in STLMatArray!

	if (BandNo==0){
	  chi_m1=chiN(0,Lambda); // is actually chi0 // chiN is a codehandler func
	// Valid ONLY for SquareBand, z_twist=1! Check side dot calcs...
	}
	else{chi_m1=chain.GetChin(0);}
	// Need to change chain files for ALL Models!!
	strcpy(ThermoArray[0].ChainArqName,"EntropyChain1Ch_PupPdn");
	strcat(ThermoArray[0].ChainArqName,zstring);
	strcat(ThermoArray[0].ChainArqName,".dat");
	ThermoArray[0].Calc=CalcEntropy;
	ThermoArray[0].CalcChain=true; // Just to make sure

	break;

// Working on this.
      case 7: 
	cout << " 2 Majoranas + QD effective model " << endl;

	Params.push_back(Lambda);  // Lambda
	Params.push_back(HalfLambdaFactor); // Lambda factor

	Params.push_back(dInitParams[0]); // U1
	Params.push_back(dInitParams[2]); // ed1
	Params.push_back(chi_m1); // sqrt(2gamma1/Pi)/sqrt(L)*HalfLambdaFactor

	Params.push_back(dInitParams[3]); // Mag Field 

	Params.push_back(dInitParams[4]); // lambdaL (or tL)
	Params.push_back(dInitParams[5]); // lambdaR
	Params.push_back(dInitParams[6]); // phi/pi

	// Initialize matrices 
	// fd_up
	STLMatArray[0].NeedOld=false;	  
 	STLMatArray[0].CheckForMatEl=OneChPupPdn_cdup_check;
 	STLMatArray[0].CalcMatElCplx=OneChPupPdn_fNup_MatElCplx;
	STLMatArray[0].SaveMatYN=false;
	STLMatArray[0].IsComplex=true;

	// fd_dn
	STLMatArray[1].NeedOld=false;
 	STLMatArray[1].CheckForMatEl=OneChPupPdn_cddn_check;
 	STLMatArray[1].CalcMatElCplx=OneChPupPdn_fNdn_MatElCplx;
	STLMatArray[1].SaveMatYN=false;
	STLMatArray[1].IsComplex=true;


	switch (calcdens){
	case 0: // OpAvg: n_d(T), S_z(T) -(2014)
	  auxNRGMat.NeedOld=true;
	  auxNRGMat.UpperTriangular=false; // update procedure only works for 'false'
	  auxNRGMat.CalcAvg=true;  // CalcAvg
	  auxNRGMat.CheckForMatEl=Diag_check;
	  auxNRGMat.CalcMatElCplx=ImpOnly_MatElCplx;
	  auxNRGMat.IsComplex=true;

	  STLMatArray[2]=auxNRGMat; // 2 - ndot (there already)
	  STLMatArray[3]=auxNRGMat; // 3 - Sz (there already)
	  strcpy(STLMatArray[2].MatName,"Ndot");
	  strcpy(STLMatArray[3].MatName,"Szdot");
	  STLMatArray.push_back(auxNRGMat);
	  strcpy(STLMatArray[4].MatName,"Nf"); // 4 -Nf
	  // changes in STLMatArray!
	  SaveData=false;
	  // No Thermo
	  NumThermoMats=0;
	  break;
 	case 1: //Spectral density
	  // cd_up and cd_dn (OneChQSz_SetAnderson equals 
	  // MatArray[2,3] to MatArray[0,1]
	  STLMatArray[2].NeedOld=true;
	  STLMatArray[2].CheckForMatEl=OneChPupPdn_cdup_check;
	  STLMatArray[2].CalcMatElCplx=OneChPupPdn_cdup_MatElCplx;
	  STLMatArray[2].SaveMatYN=true;
	  STLMatArray[2].IsComplex=true;
	  strcpy(STLMatArray[2].MatName,"cd_up");


	  STLMatArray[3].NeedOld=true;
	  STLMatArray[3].CheckForMatEl=OneChPupPdn_cddn_check;
	  STLMatArray[3].CalcMatElCplx=OneChPupPdn_cddn_MatElCplx;
	  STLMatArray[3].SaveMatYN=true;
	  STLMatArray[3].IsComplex=true;
	  strcpy(STLMatArray[3].MatName,"cd_dn");

	  // f_Maj_up (this looks ok)
	  auxNRGMat=STLMatArray[2];
	  STLMatArray.push_back(auxNRGMat);
	  strcpy(STLMatArray[4].MatName,"f_MajUp");

	  // f_Maj_dn (this looks ok)
	  auxNRGMat=STLMatArray[3];
	  STLMatArray.push_back(auxNRGMat);
	  strcpy(STLMatArray[5].MatName,"f_MajDn");

	  // TODO: Add other spectral functions??

	  cout << " N = -1: Nothing here" << endl;
// 	  cout << " fdN_up : " << endl;
// 	  STLMatArray[0].PrintAllBlocks();
// 	  cout << " fdN_dn : " << endl;
// 	  STLMatArray[1].PrintAllBlocks();
// 	  cout << " cdN_up : " << endl;
// 	  STLMatArray[2].PrintAllBlocks();
// 	  cout << " cdN_dn : " << endl;
// 	  STLMatArray[3].PrintAllBlocks();
// 	  cout << " f_MajUp : " << endl;
// 	  STLMatArray[4].PrintAllBlocks();
// 	  cout << " f_MajDn : " << endl;
// 	  STLMatArray[5].PrintAllBlocks();

	  SaveData=true;

 	  break;
	case 2:
	  cout << "Entropy calculation (implementing it...)" << endl;

	  STLMatArray.pop_back(); // Only TWO arrays
	  STLMatArray.pop_back(); // Only TWO arrays

	  // Need to change chain files for ALL Models!!
	  strcpy(ThermoArray[0].ArqName,"EntropyImp1ChPupPdn_Majorana");
	  strcat(ThermoArray[0].ArqName,zstring);
	  strcat(ThermoArray[0].ArqName,".dat");
	  strcpy(ThermoArray[0].ChainArqName,"EntropyChain1Ch_PupPdn");
	  strcat(ThermoArray[0].ChainArqName,zstring);
	  strcat(ThermoArray[0].ChainArqName,".dat");
	  ThermoArray[0].Calc=CalcEntropy;
	  ThermoArray[0].dImpValue=4.0*log(2.0); // log(16) (check)
	  ThermoArray[0].CalcChain=false;
	  NumThermoMats=1;
	  SaveData=false;

	  break;
	default:
	  cout << " calcdens = " << calcdens << " not implemented. Exiting. " << endl;
	  exit(0);
	}
	// end switch calcdens

	// actually SetHm1 !
	OneChPupPdn_SetHm1_AndersonMajorana(Params,pSingleSite,pAeig,STLMatArray);

// 	// BuildBasis params (already out of here)
	// Set pHN
	pHN->IsComplex=true;
	pHN->CalcHNMatElCplx=OneChPupPdn_HN_MatElCplx;

	// Wrap up
	NumChannels=1;
	Nsites0=0;
	NumNRGmats=STLMatArray.size();
	MatArray=&STLMatArray[0]; // Need to do this after
	                          // changes in STLMatArray!

	cout << " N=-1 ok. Going to N=0... " << endl;
	
	cout << " pAeig is complex ? " << pAeig->CheckComplex() << endl;


	//exit(0);
	break;
      default:
	cout << " Model not implemented for this symmetry. Exiting... " << endl;
	exit(0);
      }  
      // end switch ModelNo for OneChPupPdn
      break;
/////////////////////////////    
    default:
      cout << " Symmetry not implemented. Exiting... " << endl;
      exit(0);    
    }
    // end switch SymNo

}
// end ModelSwitch
