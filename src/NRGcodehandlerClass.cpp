
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "NRGbandtypes.hpp"
#include "OneChQS.hpp"
#include "TwoChQS.hpp"
#include "NRG_main.hpp"
using namespace std;



//////////////////////////

void CNRGCodeHandler::DoesNothing(){

  // instream
  ifstream InFile;
  
  //InFile.open("input_nrg.dat");

   CodeInFile.open("input_nrg.dat");
   CodeInFile.close();

}


//////////////////////////

//void CNRGCodeHandler::InitialSetUp(){
// Trying this
void CNRGCodeHandler::InitialSetUp(bool ReadParamsOnly){

  strcpy(ParamFileName,"input_nrg.dat");
  strcpy(BeginFileName,"NRG_in.txt");
  strcpy(EndFileName,"NRG_end.txt");
  strcpy(SaveArraysFileName,"test");
  strcpy(LancInFileName,"lanc.in");

  NoInputParamsDouble=3;
 
  switch(ModelNo){
  case 0: // 0 - Anderson, 1 - Kondo , 2 - Phonon
    if (SymNo==2){NoInputParamsDouble=4;} // Add Bmag in Q Sz
    if (SymNo==7){NoInputParamsDouble=4;} // Add DeltaSC in S
    if (SymNo==8){NoInputParamsDouble=5;} // Add DeltaSC and B in Sz
    break;
  case 1:
    if (SymNo==1){NoInputParamsDouble=5;} // Add J1 in TwoCh_Kondo
    break;
  case 3: // Chain
    NoInputParamsDouble=0;
    break;
  case 4: // 2ch Phonon
    NoInputParamsDouble=9;
    break;
  case 5: // SMM
    NoInputParamsDouble=11;
    break;
  case 6: // DQD
    NoInputParamsDouble=6; // U1 ed1 gamma1 U2 ed2 gamma2
    if (SymNo==2){NoInputParamsDouble=9;} // Add lambda, B1, B2
    break;
  case 7: // Majorana
    NoInputParamsDouble=10;
    break;

  }

  ReadParams();

  if (!ReadParamsOnly){
    if ( (calcdens==3)&&(ModelNo!=3) ){
      cout << " Chain calculation detected: Changing ModelNo from " 
	   << ModelNo <<" to 3 (chain)" << endl;
      ModelNo=3;
      ZeroParams();
    }
    
    SetChain();
    
    SaveBegFile();
    SetTotS();
  }
  // if not ReadParamsOnly

}

//////////////////////////

void CNRGCodeHandler::ReadParams(){

  double auxD;
  int auxI=0;

  char Nothing[32];
  int iNothing;

  // instream
  ifstream InFile;

  CodeInFile.open(ParamFileName,ifstream::in);
  if (CodeInFile.is_open())
    {
      CodeInFile >> Nsitesmax;
      CodeInFile >> Ncutoff;
      CodeInFile >> auxD;
      dInitParams.push_back(auxD); // U
      CodeInFile >> auxD;
      dInitParams.push_back(auxD); // Gamma
      CodeInFile >> auxD;
      dInitParams.push_back(auxD); // ed
      CodeInFile >> Lambda;
//       CodeInFile >> Dband;
      CodeInFile >> ChemPot;
      CodeInFile >> auxI;
      CodeInFile >> UpdateBefCut;
      CodeInFile >> calcdens;
      CodeInFile >> iNothing; // separator
      for (int ii=3;ii<NoInputParamsDouble;ii++)
	 {CodeInFile >> auxD;dInitParams.push_back(auxD);}
    }
  else
    {
      cout << "can't open " << ParamFileName << endl;
      exit(0);
    }
  CodeInFile.close();

  // added May 2010
  if (Lambda>0.0)
    //HalfLambdaFactor=0.5*(1.0+(1.0/Lambda));
    HalfLambdaFactor=0.5*(1.0+(1.0/Lambda))*pow(Lambda,-(chain.z_twist)+1.0);
  else
    HalfLambdaFactor=0.0;
}


//////////////////////////

void CNRGCodeHandler::ReadParams(char FileName[], int NoAdditionalParams){

  double auxD;
  int auxI=0;

  char Nothing[32];
  int iNothing;

  // instream
  ifstream InFile;

  InFile.open(FileName,ifstream::in);
  if (InFile.is_open()){
    InFile >> Nsitesmax;
    InFile >> Ncutoff;
    InFile >> auxD;
    dInitParams.push_back(auxD); // U
    InFile >> auxD;
    dInitParams.push_back(auxD); // Gamma
    InFile >> auxD;
    dInitParams.push_back(auxD); // ed
    InFile >> Lambda;
//     InFile >> Dband;
    InFile >> ChemPot;
    InFile >> auxI;
    InFile >> UpdateBefCut;
    InFile >> calcdens;
    InFile >> iNothing; // separator
    for (int ii=0;ii<NoAdditionalParams;ii++)
      {InFile >> auxD;dInitParams.push_back(auxD);}
  }
  else{
    cout << "can't open " << FileName << endl;
    exit(0);
  }
  InFile.close();

  // added May 2010
//   if (Lambda>0.0)
//     //HalfLambdaFactor=0.5*(1.0+(1.0/Lambda));
//     HalfLambdaFactor=0.5*(1.0+(1.0/Lambda))*pow(Lambda,-(chain.z_twist)+1.0);
//   else
//     HalfLambdaFactor=0.0;
}


//////////////////////////

void CNRGCodeHandler::SetChain(){

  vector<double> dParams;
  chain.Lambda=Lambda;
  chain.HybFuncParams=&dParams;
  // new (Jun 2013)
  chain.BandOffset=ChemPot;

  double daux=0.0;

  int Nsiteschain=Nsitesmax>100?Nsitesmax:100;

  switch(BandNo){
  case 0: // Square but no
    chain.SetChainWilson(Nsiteschain);
    break;
  case 4: // Lorentzian
    // Read Params from lanc.in
    chain.HybFunction=HybFunc_Lorentzian;
    chain.HybFuncWithEn=HybFunc_Lorentzian_timesEn;
    chain.HybFuncParams=&dParams;
    chain.ReadParams(LancInFileName,6);
    chain.SetChainLanczos(Nsiteschain);
    break;
  case 41: // Cavity: multi_peak LorenztianLorentzian
    // Read Params from lanc.in
    chain.HybFunction=HybFunc_Cavity;
    chain.HybFuncWithEn=HybFunc_Cavity_timesEn;
    chain.HybFuncParams=&dParams;
    chain.ReadParams(LancInFileName,8);
    cout << " Testing ..." << endl;
    cout << " Delta(-0.3)= " << chain.GetHyb_w(-0.3) << endl;
    cout << " Delta(-0.2)= " << chain.GetHyb_w(-0.2) << endl;
    cout << " Delta(-0.1)= " << chain.GetHyb_w(-0.1) << endl;
    chain.SetChainLanczos(Nsiteschain);
    break;
  case 11: // Const DoS/ConsHyb
    chain.HybFunction=HybFunc_Const;
    switch (chain.DiscScheme){
    case 0:
      chain.HybFuncWithEn=HybFunc_Const_timesEn;
      break;
    case 1: // Campo-Oliveira
      chain.HybFuncWithEn=HybFunc_Const_divEn;
      break;
    default:
      chain.HybFuncWithEn=HybFunc_Const_timesEn;
    }
    // end switch disc scheme
    chain.HybFuncParams=&dParams; // HybFuncParams points to dParams.
    if ( (dInitParams.size()>1)&&(dNEqual(dInitParams[1],0.0)) )
      dParams.push_back(dInitParams[1]); // Parameter is Gamma
    else
      dParams.push_back(1.0);
    // z_twist comes from command line
    cout << " z_twist = " << chain.z_twist << endl;
    chain.SetChainLanczos(Nsiteschain); // Should use z_twist...
    break;
  case 12: // PowerLaw
    // Read Params from lanc.in
    chain.HybFunction=HybFunc_PowerLaw;
    chain.HybFuncWithEn=HybFunc_PowerLaw_timesEn;
    chain.HybFuncParams=&dParams; // HybFuncParams points to dParams.
    chain.ReadParams(LancInFileName,5);
    cout << " Power-law band: Gamma(w)=Gamma0*|w-w0|^r + small_gamma " << endl;
    cout << " whichbandtype = " << dParams[0]
	 << " r = " << dParams[1]
	 << " Gamma0 = " << dParams[2]
	 << " small_gamma = " << dParams[3]
	 << " w0 = " << dParams[4]
	 << endl;
    // Setting HybRef (will always be Gamma0). Needs to change!!
//     chain.HybRef=dParams[2];
//     chain.HybRefIsSet=true;
    daux=chain.GetHyb_w(0.0);
    if (dEqual(daux,0.0)){
      cout << "  SetChain: Hyb(0)=0. Setting HybRef to " << dParams[2] << endl;
      chain.HybRef=dParams[2];
    }else{chain.HybRef=daux;}
    chain.HybRefIsSet=true;
    // Set Nsiteschain according with Delta!!
    chain.SetChainLanczos(Nsiteschain);
    break;
  case 13: // From File
    // Read Params from lanc.in
    chain.HybFunction=HybFunc_FromFile;
//     chain.HybFuncWithEn=HybFunc_FromFile_timesEn;
    switch (chain.DiscScheme){
    case 0:
      chain.HybFuncWithEn=HybFunc_FromFile_timesEn;
      break;
    case 1: // Campo-Oliveira
      chain.HybFuncWithEn=HybFunc_FromFile_divEn;
      break;
    default:
      chain.HybFuncWithEn=HybFunc_FromFile_timesEn;
    }
    // end switch disc scheme
    chain.HybFuncParams=&dParams; 
    // We need HybFuncParams to be a single STL 
    // vector otherwise the GSL function won't work... 
    char HybFuncDatFileName[64];
    strcpy(HybFuncDatFileName,"HybFunc.dat");
    chain.ReadParams(HybFuncDatFileName); // Reads the whole thing.
    daux=chain.GetHyb_w(0.0);
    cout << " Hyb(0.0) = " << daux << endl;
    // NOTE: this file is the FULL hybridization:
    // Gamma(e)=pi*t(e)^2*rho(e)
    // If Gamma0=0.0, set HybRef.
    if ( dEqual(daux,0.0) ){
      cout << "  SetChain: Hyb(0)=0. Setting HybRef to 1.0 " << endl;
      chain.HybRef=1.0;
      chain.HybRefIsSet=true;
    }
    // checking z-trick
    cout << " z_twist = " << chain.z_twist << endl;
    chain.SetChainLanczos(Nsiteschain);
    //chain.PrintChain(20);
    //exit(0);
    break;
  case 14: // Sum of Deltas 
    // Read Params from lanc.in
    chain.HybFunction=HybFunc_FromFile;
    chain.HybFuncWithEn=HybFunc_FromFile_timesEn;
    chain.HybFuncParams=&dParams; 
    // We need HybFuncParams to be a single STL 
    // vector otherwise the GSL function won't work... 
    //char HybFuncDatFileName[64];
    strcpy(HybFuncDatFileName,"HybDeltas.dat");
    chain.ReadParams(HybFuncDatFileName); // Reads the whole thing.
    daux=0.0;
    //daux=chain.GetHyb_w(0.0);
    // If Gamma0=0.0, set HybRef.
    if ( dEqual(daux,0.0) ){
      cout << "  SetChain: Hyb(0)=0. Setting HybRef to 1.0 " << endl;
      chain.HybRef=1.0;
      chain.HybRefIsSet=true;
    }
    // Testing graphene
    chain.HybFuncIsSumDeltas=true;
    chain.SetChainLanczos(Nsiteschain);
    chain.PrintChain(20);
    exit(0);
    break;
  default:
    cout << "BandNo = " << BandNo << "not implemented. Exiting..." << endl;
    exit(0);
    //chain.SetChainWilson(Nsiteschain);
  }
  // end switch band type

  if (BandNo!=0){
    cout << " Special band used: " << BandType << endl;
    if (dInitParams.size()>1){
      //cout << " Changed Gamma to 1/2*Fsq*Gamma(0) = " << chain.HybRef << endl;
      dInitParams[1]=0.5*chain.Fsq*chain.HybRef;
      cout << " Changed Gamma to 1/2*Fsq*Gamma(0) = " << dInitParams[1] << endl;

    }
    chain.PrintChain(20);
    // Debugging. Remove later.
    //if (BandNo==12) exit(0);
  }
  // end if BandNo!=0



}


//////////////////////////
void CNRGCodeHandler::SaveBegFile()
{

  // outstream
  ofstream OutFile;


  OutFile.open(BeginFileName);
  OutFile << "Begin NRG calculation. Model :" << ModelOption 
  	  << " - Symmetry :" << Symmetry 
	  << " - Band Type: " << BandType  << endl;
  OutFile << " Nsitesmax    = " << Nsitesmax-1 << endl;
  OutFile << " Ncutoff      = " << Ncutoff << endl;
  OutFile << " U            = " << dInitParams[0] << endl;
  OutFile << " Gamma        = " << dInitParams[1] << endl;
  OutFile << " ed           = " << dInitParams[2] << endl;
  OutFile << " Lambda       = " << Lambda << endl;
//   OutFile << " Dband        = " << Dband << endl;
  OutFile << " ChemPot      = " << ChemPot << endl;
  OutFile << " UpdateBefCut = " << UpdateBefCut << endl;
  OutFile << " calcdens     = " << calcdens << endl;
  OutFile << " Oliveira z   = " << chain.z_twist << endl;
  OutFile << " Additional parameters : " << endl;
  for (int ii=3;ii<NoInputParamsDouble;ii++)
    OutFile << " Param " << ii << "  = " << dInitParams[ii] << endl;
  OutFile.close();


  // couts as well

  cout << " NRG: Model : " << ModelOption 
       << " - Symmetry : " << Symmetry 
       << " - Band Type: " << BandType  << endl;
  cout << " Nsitesmax    = " << Nsitesmax-1 << endl;
  cout << " Ncutoff      = " << Ncutoff << endl;
  cout << " U            = " << dInitParams[0] << endl;
  cout << " Gamma        = " << dInitParams[1] << endl;
  cout << " ed           = " << dInitParams[2] << endl;
  if (strcmp(ModelOption,"Kondo")==0){
    cout << " JK            = " <<  dInitParams[1] << endl;
  }
  cout << " Lambda       = " << Lambda << endl;
//   cout << " Dband        = " << Dband << endl;
  cout << " ChemPot      = " << ChemPot << endl;
  cout << " UpdateBefCut = " << UpdateBefCut << endl;
  cout << " calcdens     = " << calcdens << endl;
  cout << " Oliveira z   = " << chain.z_twist << endl;
  cout << " Additional parameters : " << endl;
  for (int ii=3;ii<NoInputParamsDouble;ii++)
    cout << " Param " << ii << "  = " << dInitParams[ii] << endl;


}



//////////////////////////

void CNRGCodeHandler::SetSingleSite(CNRGbasisarray* pSingleSite){

  cout << " SetSingleSite begin... " << endl; 

  switch (SymNo){
  case 0:
    cout << " Got OneChQS " << endl;
    OneChQS_SetSingleSite(*pSingleSite);
    break;
  case 1:
    cout << " Got TwoChQS " << endl;
    TwoChQS_SetSingleSite(pSingleSite);
    break;
  case 2:
    OneChQSz_SetSingleSite(*pSingleSite);
    break;
  case 3:
    cout << " Got OneChQ " << endl;
    pSingleSite->NQNumbers=1;
    pSingleSite->Nshell=0;
    // Q=-1 (one state, S=0)
    pSingleSite->QNumbers.push_back(-1.0);
    pSingleSite->BlockBegEnd.push_back(0);
    pSingleSite->BlockBegEnd.push_back(0);
    // Q=0 (two states, S=1/2)
    pSingleSite->QNumbers.push_back(0.0);
    pSingleSite->BlockBegEnd.push_back(1);
    pSingleSite->BlockBegEnd.push_back(2);
    // Q=1 (One state, S=0)
    pSingleSite->QNumbers.push_back(1.0);
    pSingleSite->BlockBegEnd.push_back(3);
    pSingleSite->BlockBegEnd.push_back(3);
    for (int ii=0;ii<4;ii++){
      // Type labels the state
      pSingleSite->iType.push_back(ii);
    }
    break;
  case 4:
    TwoChQSP_SetSingleSite(pSingleSite); // old
    ////TwoChQSP_SetSingleSite(*pSingleSite); // new
    break;
  case 6:
    cout << " Got OneChNupPdn " << endl;
    pSingleSite->NQNumbers=2;
    pSingleSite->Nshell=0;
    // |0>: Nup=0 Pdn=+1 
    pSingleSite->QNumbers.push_back(0.0);
    pSingleSite->QNumbers.push_back(1.0);
    pSingleSite->BlockBegEnd.push_back(0);
    pSingleSite->BlockBegEnd.push_back(0);
    // |up>: Nup=1 Pdn=+1 
    pSingleSite->QNumbers.push_back(1.0);
    pSingleSite->QNumbers.push_back(1.0);
    pSingleSite->BlockBegEnd.push_back(1);
    pSingleSite->BlockBegEnd.push_back(1);
    // |dn>: Nup=0 Pdn=-1 
    pSingleSite->QNumbers.push_back(0.0);
    pSingleSite->QNumbers.push_back(-1.0);
    pSingleSite->BlockBegEnd.push_back(2);
    pSingleSite->BlockBegEnd.push_back(2);
    // |up dn>: Nup=1 Pdn=-1 
    pSingleSite->QNumbers.push_back(1.0);
    pSingleSite->QNumbers.push_back(-1.0);
    pSingleSite->BlockBegEnd.push_back(3);
    pSingleSite->BlockBegEnd.push_back(3);
    for (int ii=0;ii<4;ii++){
      // Type labels the state
      pSingleSite->iType.push_back(ii);
    }
    break;
  case 7:
    cout << " Got OneChS " << endl;
    pSingleSite->NQNumbers=2;
    pSingleSite->Nshell=0;
    pSingleSite->totalS=true;
    pSingleSite->Sqnumbers.push_back(0);    
    // |0>, |up dn>: S=0 Sz=0 
    pSingleSite->QNumbers.push_back(0.0);
    pSingleSite->QNumbers.push_back(0.0);
    pSingleSite->BlockBegEnd.push_back(0);
    pSingleSite->BlockBegEnd.push_back(1);
    // |up>: S=0.5 Sz=0.5 
    pSingleSite->QNumbers.push_back(0.5);
    pSingleSite->QNumbers.push_back(0.5);
    pSingleSite->BlockBegEnd.push_back(2);
    pSingleSite->BlockBegEnd.push_back(2);
    // |dn>: S=0.5 Sz=-0.5 
    pSingleSite->QNumbers.push_back(0.5);
    pSingleSite->QNumbers.push_back(-0.5);
    pSingleSite->BlockBegEnd.push_back(3);
    pSingleSite->BlockBegEnd.push_back(3);
    for (int ii=0;ii<4;ii++){
      // Type labels the state
      pSingleSite->iType.push_back(ii);
    }
    break;
  case 8:
    cout << " Got OneChSz " << endl;
    pSingleSite->NQNumbers=1;
    pSingleSite->Nshell=0;
    pSingleSite->totalS=false;
    // |0>, |up dn>: Sz=0 
    pSingleSite->QNumbers.push_back(0.0);
    pSingleSite->BlockBegEnd.push_back(0);
    pSingleSite->BlockBegEnd.push_back(1);
    // |up>: S=0.5 Sz=0.5 
    pSingleSite->QNumbers.push_back(0.5);
    pSingleSite->BlockBegEnd.push_back(2);
    pSingleSite->BlockBegEnd.push_back(2);
    // |dn>: S=0.5 Sz=-0.5 
    pSingleSite->QNumbers.push_back(-0.5);
    pSingleSite->BlockBegEnd.push_back(3);
    pSingleSite->BlockBegEnd.push_back(3);
    for (int ii=0;ii<4;ii++){
      // Type labels the state
      pSingleSite->iType.push_back(ii);
    }
    break;
  default:
    cout << Symmetry << " symmetry not implemented. Exiting... " <<endl;
    exit(0);
  }


  cout << " .. SetSingleSite done!" << endl; 
}

//////////////////////////

void CNRGCodeHandler::CalcThermo(CNRGthermo* ThermoArray, CNRGarray* pAeig){

  // Needs Nsites, Nsites0(!), DN and betabar set 

  vector<double> Params;

  // Added a loop over betabar values:
  // (For now, let's do four betabar values per iteration)
  // betabar0:=ThisCode.betabar;
  // betabar=betabar0/Lambda^{n/8})

  int Nmax=4; //=8;

// Need nc=0 to compare same temps and diff Nsites
// BUT let's try with nc!=Nmax for now
  //for (int nc=0; nc<=Nmax; nc++){ // Compare with Kevin's 
  for (int nc=1; nc<=Nmax; nc++){ // No repeated TM values 
    //for (int nc=0; nc<Nmax; nc++){ // Old stuff
    //double betabarN=betabar*pow(Lambda,nc/(2.0*Nmax)); // Old stuff
    // THIS one works!!
    double betabarN=betabar*pow(Lambda,(nc-Nmax)/(2.0*Nmax));

    if ( (calcdens==2)||(calcdens==3) ){
      double TM=DN/betabarN;
      Params.clear();
      Params.push_back(betabarN);
      cout  << " CalcThermo: Nsites = " << Nsites
	    << " Nsites0 = " << Nsites0 
	    << " nc = " << nc << " of " << Nmax
	    << " TM = " << TM
	//	    << " Nvalue = " << Nsites*(Nmax+1)+nc
	    << " Nvalue = " << Nsites*Nmax+nc-1
	    << endl;

      for (int ithermo=0;ithermo<NumThermoMats;ithermo++){
	//ThermoArray[ithermo].ReadNChainValue(Nsites,0);
	//  OLD Now there are Nmax+1 temps/values for each Nsites!
	//ThermoArray[ithermo].ReadNChainValue(Nsites*(Nmax+1)+nc,0);
	//  NEW Now there are Nmax temps/values for each Nsites!
	ThermoArray[ithermo].ReadNChainValue(Nsites*Nmax+nc-1,0);
	ThermoArray[ithermo].AddValue(Params,pAeig,Sqnumber,totalS,TM);
	//ThermoArray[ithermo].SaveNValue(Nsites,0);
	//ThermoArray[ithermo].SaveNValue(Nsites*(Nmax+1)+nc,0);
	ThermoArray[ithermo].SaveNValue(Nsites*Nmax+nc-1,0);

      }
    }
    else
      cout << "CalcThermo: Not calculating thermodynamics " << endl;
    //end if calcdens=2 or 3

  }
  // end loop in betabar values


}
//////////////////////////

// Still working on this
void CNRGCodeHandler::CalcStuff(CNRGarray* pAeig, CNRGbasisarray* pAeigCut){

  vector <double> ParamsBetabar;

  ParamsBetabar.clear();
  ParamsBetabar.push_back(betabar);

  // HOW ABOUT THIS INSTEAD??
  // June 2013: taking the calcdens constraint.
  // Now we can "CalcStuff" for any calcdens
  //if (calcdens==0){
  for (int imat=0;imat<NumNRGmats;imat++){
    //       cout << "imat = " << imat
    // 	   << " CalcAvg = " << MatArray[imat].CalcAvg
    // 	   << endl;
    if (MatArray[imat].CalcAvg){
      double MatAvg=CalcOpAvg(ParamsBetabar,pAeigCut,&MatArray[imat],totalS,Sqnumber);
      cout << " T= " << DN/betabar << " " << MatArray[imat].MatName << "= " << MatAvg  << endl;
    }else{
      cout << " CalcStuff: Not Calculating stuff for Mat " << imat << endl;
    }
  }
  // end loop in mats

}
// end CalcStuff

//////////////////////////

void CNRGCodeHandler::WrapUp(){


  // outstream
  ofstream OutFile;

  cout << "=== Calculation Finished! ==== "<< endl;
  OutFile.open(EndFileName);
  OutFile << "END NRG calculation" << endl;
  OutFile.close();


}


///////////////////////////

void CNRGCodeHandler::ZeroParams(){

  vector<double>::iterator it;

  for (it=dInitParams.begin();it<dInitParams.end();it++)
    *it=0.0;
  

}


///////////////////////////


void CNRGCodeHandler::SetTotS(){


  switch (SymNo){
  case 0: // OneChQS
    totalS=true;
    Sqnumber=1;
    break;
  case 1: // TwoChQS
    totalS=true;
    Sqnumber=1;
    break;
  case 2: // OneChQSz
    totalS=false;
    Sqnumber=1;
    break;
  case 3: // OneChQ
    totalS=false;
    Sqnumber=-1; // no spin whatsoever
    break;
  case 4: // TwoChQSP
    totalS=true;
    Sqnumber=1;
    break;
  case 5: // TwoChQSz
    totalS=false;
    Sqnumber=1;
    break;
  case 6: // OneChNupPdn
    totalS=false;
    Sqnumber=-1; // no spin whatsoever
    break;
  case 7: // OneChS
    totalS=true;
    Sqnumber=0;
    break;
  case 8: // OneChSz
    totalS=false;
    Sqnumber=0;
    break;
  default:
    totalS=true;
    Sqnumber=1;
  }
  // end switch


}


///////////////////////////


/// Mar. 09

//void CNRGCodeHandler::SaveGenPars(char idname[]){
void CNRGCodeHandler::SaveGenPars(){

  ofstream OutFile;
  char arqname[32];

  // Save general parameters
  
  strcpy(arqname,"ThisCodePars.dat");
  //strcat(arqname,ext);

  SetNopsSaved();

  OutFile.open(arqname, ios::out);
  if (!OutFile){cout << "SaveGenPars: Cannot save data in " << arqname << endl; return;}

  OutFile << Lambda << endl;
  OutFile << Nsites0 << endl;
  OutFile << Nsitesmax << endl;


  OutFile << NopsSaved << endl;
  OutFile << SaveArraysFileName << endl;
  OutFile << SymNo << endl;
  // NEW LINE
  OutFile << ModelNo << endl;
  OutFile << BandNo << endl;
  OutFile << chain.z_twist << endl;

  
  OutFile << endl;
  OutFile << "# Lambda used" << endl;
  OutFile << "# Nsites0 used " << endl;
  OutFile << "# Nsitesmax used " << endl;
  OutFile << "# Num Ops saved " << endl;
  OutFile << "# File extension used " << endl;
  OutFile << "# Symmetry No " << endl;
  OutFile << "# Model No  " << endl;
  OutFile << "# Band No  " << endl;
  OutFile << "# Oliveira z  " << endl;
  OutFile << "###################" << endl;
  OutFile << "# Symmetry Numbers: " << endl;
  OutFile << "#  0 - OneChQS " << endl;
  OutFile << "#  1 - TwoChQS " << endl;
  OutFile << "#  2 - OneChQSz " << endl;
  OutFile << "#  3 - OneChQ " << endl;
  OutFile << "#  4 - TwoChQSP " << endl;
  OutFile << "#  5 - TwoChQSz " << endl;
  OutFile << "#  6 - OneChNupPdn " << endl;
  OutFile << "#  7 - OneChS " << endl;
  OutFile << "#  8 - OneChSz " << endl;


//  
//  
// 
//  
//  


  OutFile.close();



}


///////////////////////////

void CNRGCodeHandler::SaveArrays(){

  //char idname[]="";
  //SaveArrays(idname);

  SaveArrays(SaveArraysFileName);

}


///////////////////////////

void CNRGCodeHandler::SaveArrays(char idname[]){

      char CNsites[8];
      char CNmat[8];
      char arqname[32];
      char ext[32];

      sprintf(CNsites,"%d",pAcut->Nshell);

      strcpy(ext,idname);
      strcat(ext,"_N");
      strcat(ext,CNsites);
      strcat(ext,".bin");

      // Save Abasis
      strcpy(arqname,"Abasis_");
      strcat(arqname,ext);
      pAbasis->SaveBin(arqname);

      // Save Acut
      strcpy(arqname,"Acut_");
      strcat(arqname,ext);
      pAcut->SaveBin(arqname);


      cout << "CodeHandler: Saving matrices " << endl;
      // Save Amatrices (not all of them)
      int isave=0;
      for (int imat=0;imat<NumNRGmats;imat++){
	//sprintf(CNmat,"%d",imat);
	cout << "imat = " << imat
	     << " SaveMatYN = " << MatArray[imat].SaveMatYN
	     << endl;
	if (MatArray[imat].SaveMatYN){
	  sprintf(CNmat,"%d",isave);
	  strcpy(arqname,"Mat");
	  strcat(arqname,CNmat);
	  strcat(arqname,"_");
	  strcat(arqname,ext);
	  MatArray[imat].SaveBin(arqname);
	  isave++;
	}
	// end loop in isave
      }
      // end loop in imat
}


///////////////////////////

// Keep this until all ThisCodePars.dat have been updated :)
//void CNRGCodeHandler::ReadGenPars(){
void CNRGCodeHandler::ReadGenPars(bool ReadModelBand){

  ifstream ReadFile;
  char arqname[32];

  // Read general parameters
  
  strcpy(arqname,"ThisCodePars.dat");
  //strcat(arqname,ext);

  // Just Lambda, Nsitesmax, Ncutoff for now...

  ReadFile.open(arqname, ios::in);
  if (!ReadFile){cout << "ReadGenPars: Cannot read data from " << arqname << endl; return;}

//   ReadFile.read((char*)&Lambda, sizeof(double));
//   ReadFile.read((char*)&Nsitesmax, sizeof(int));
  ReadFile >> Lambda;
  ReadFile >> Nsites0;
  ReadFile >> Nsitesmax;

  ReadFile >> NopsSaved;
  ReadFile >> SaveArraysFileName;
  ReadFile >> SymNo;
  // NEW LINE
  if (ReadModelBand){
    ReadFile >> ModelNo;
    ReadFile >> BandNo;
    ReadFile >> chain.z_twist;
  }
  
  ReadFile.close();


  // new
  SetTotS();

}


///////////////////////////

void CNRGCodeHandler::ReadArrays(){

  //char idname[]="";
  //ReadArrays(idname);

  ReadArrays(SaveArraysFileName);

}

///////////////////////////


void CNRGCodeHandler::ReadArrays(char idname[]){


      char CNsites[8];
      char CNmat[8];
      char arqname[32];
      char ext[32];

      sprintf(CNsites,"%d",pAcut->Nshell);

      strcpy(ext,idname);
      strcat(ext,"_N");
      strcat(ext,CNsites);
      strcat(ext,".bin");

      // Read Abasis
      strcpy(arqname,"Abasis_");
      strcat(arqname,ext);
      pAbasis->ReadBin(arqname);

      // Read Acut
      strcpy(arqname,"Acut_");
      strcat(arqname,ext);
      pAcut->ReadBin(arqname);

      // Read Amatrices
      for (int imat=0;imat<NumNRGmats;imat++){
	sprintf(CNmat,"%d",imat);
	strcpy(arqname,"Mat");
	strcat(arqname,CNmat);
	strcat(arqname,"_");
	strcat(arqname,ext);
	MatArray[imat].ReadBin(arqname);
      }


}

///////////////////////////

///////////////////////////

// double CNRGCodeHandler::CalcDN(int iNshell){

//   double HalfLambdaFactor=1.0;

//   if (Lambda>0.0)
//     HalfLambdaFactor=0.5*(1.0+(1.0/Lambda));

//   return(HalfLambdaFactor*pow(Lambda,(-(iNshell-1)/2.0) ));


// }

///////////////////////////

void CNRGCodeHandler::SetCurrentDN(){

//   DN=CalcDN(Lambda,Nsites);
  DN=CalcDN(Lambda,Nsites,chain.z_twist);

}


///////////////////////////

void CNRGCodeHandler::SetNopsSaved(){

  NopsSaved=0;
  for (int imat=0;imat<NumNRGmats;imat++){
    if (MatArray[imat].SaveMatYN){NopsSaved++;}
  }
    
}


///////////////////////////


bool CNRGCodeHandler::CheckFileExists(char arqname[]){

  // Checks whether file arqname exists or not

  // ifstream
  ifstream InFile(arqname);

  return(InFile); // casts InFile into bool

}

///////////////////////////

double CNRGCodeHandler::chiN(int Nsites, double Lambda)
{

  // Valid ONLY for z_twist=1!!
  double daux[3];
  daux[0]=1.0-pow( Lambda,-((double)(Nsites)+1.0) );
  daux[1]=sqrt( 1.0-pow(Lambda,-(2.0*(double)(Nsites)+1.0)) );
  daux[2]=sqrt( 1.0-pow(Lambda,-(2.0*(double)(Nsites)+3.0)) );  

  return(daux[0]/(daux[1]*daux[2]));

}

///////////////////////////


void CNRGCodeHandler::PrintSettings(){

  cout << " ====== Code Settings : ========= " << endl;

  cout << " ModelNo = " << ModelNo << endl
       << " SymNo = " << SymNo << endl
       << " BandNo = " << BandNo << endl;
  cout << " Nsitesmax    = " << Nsitesmax-1 << endl;
  cout << " Ncutoff      = " << Ncutoff << endl;
  cout << " U            = " << dInitParams[0] << endl;
  cout << " Gamma        = " << dInitParams[1] << endl;
  cout << " ed           = " << dInitParams[2] << endl;
  cout << " Lambda       = " << Lambda << endl;
  //  cout << " Dband        = " << Dband << endl;
  cout << " ChemPot      = " << ChemPot << endl;
  cout << " UpdateBefCut = " << UpdateBefCut << endl;
  cout << " calcdens     = " << calcdens << endl;
  cout << " Oliveira z   = " << chain.z_twist << endl;
  cout << " Additional parameters : " << endl;
  for (int ii=3;ii<NoInputParamsDouble;ii++)
    cout << " Param " << ii << "  = " << dInitParams[ii] << endl;

  cout << " NopsSaved    = " << NopsSaved << endl;


}


///////////////////////////
// Trash


