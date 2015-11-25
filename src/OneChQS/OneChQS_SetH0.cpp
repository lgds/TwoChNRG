
#include <iostream>
#include <vector>
#include <math.h>


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"



void OneChQS_SetAndersonHm1(vector<double> Params,
			    CNRGarray* pAeig, 
			    CNRGmatrix* pQm1fNQ,
			    CNRGmatrix* NRGMats){

  // Set initial CNRG array (N=-1)

  double U=Params[0];
  double ed=Params[1];
  double Lambda=Params[2];
  double HalfLambdaFactor=Params[3];


  pAeig->NQNumbers=2;
  pAeig->Nshell=-1;
  

  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.5);
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(0.0);

  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);

//   pAeig->dEn.push_back(0.0);
//   pAeig->dEn.push_back(ed/(Lambda*HalfLambdaFactor));
//   pAeig->dEn.push_back((2.0*ed+U)/(Lambda*HalfLambdaFactor));

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor));


  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);

  pAeig->PrintEn();

  // Set Matrix elements

  //fN

  pQm1fNQ->SyncNRGarray(*pAeig);

  pQm1fNQ->UpperTriangular=false;
 
  pQm1fNQ->MatEl.push_back(1.0);
  pQm1fNQ->MatBlockMap.push_back(0);
  pQm1fNQ->MatBlockMap.push_back(1);
  pQm1fNQ->MatBlockBegEnd.push_back(0);
  pQm1fNQ->MatBlockBegEnd.push_back(0);

  pQm1fNQ->MatEl.push_back(-sqrt(2.0));
  pQm1fNQ->MatBlockMap.push_back(1);
  pQm1fNQ->MatBlockMap.push_back(2);
  pQm1fNQ->MatBlockBegEnd.push_back(1);
  pQm1fNQ->MatBlockBegEnd.push_back(1);


  // cd
  // Will this work?  
  NRGMats[1]=*pQm1fNQ;
  NRGMats[1].UpperTriangular=false;

  // n_{d up}+ n_{d dn} (diagonal)


  NRGMats[0].SyncNRGarray(*pAeig);
  NRGMats[0].UpperTriangular=false;

 
  NRGMats[0].MatEl.push_back(0.0);
  NRGMats[0].MatBlockMap.push_back(0);
  NRGMats[0].MatBlockMap.push_back(0);
  NRGMats[0].MatBlockBegEnd.push_back(0);
  NRGMats[0].MatBlockBegEnd.push_back(0);

  NRGMats[0].MatEl.push_back(1.0);
  NRGMats[0].MatBlockMap.push_back(1);
  NRGMats[0].MatBlockMap.push_back(1);
  NRGMats[0].MatBlockBegEnd.push_back(1);
  NRGMats[0].MatBlockBegEnd.push_back(1);

  NRGMats[0].MatEl.push_back(2.0);
  NRGMats[0].MatBlockMap.push_back(2);
  NRGMats[0].MatBlockMap.push_back(2);
  NRGMats[0].MatBlockBegEnd.push_back(2);
  NRGMats[0].MatBlockBegEnd.push_back(2);

  // Sep 08:
  // Sz^2 (diagonal)


  NRGMats[2].SyncNRGarray(*pAeig);
  NRGMats[2].UpperTriangular=false;

 
  NRGMats[2].MatEl.push_back(0.0);
  NRGMats[2].MatBlockMap.push_back(0);
  NRGMats[2].MatBlockMap.push_back(0);
  NRGMats[2].MatBlockBegEnd.push_back(0);
  NRGMats[2].MatBlockBegEnd.push_back(0);

  NRGMats[2].MatEl.push_back(0.25);
  NRGMats[2].MatBlockMap.push_back(1);
  NRGMats[2].MatBlockMap.push_back(1);
  NRGMats[2].MatBlockBegEnd.push_back(1);
  NRGMats[2].MatBlockBegEnd.push_back(1);

  NRGMats[2].MatEl.push_back(0.0);
  NRGMats[2].MatBlockMap.push_back(2);
  NRGMats[2].MatBlockMap.push_back(2);
  NRGMats[2].MatBlockBegEnd.push_back(2);
  NRGMats[2].MatBlockBegEnd.push_back(2);



}

////////////
////////////
////////////

void OneChQS_SetKondoH0(vector<double> Params,
			CNRGarray* pAeig, CNRGmatrix* pQm1fNQ, 
			CNRGbasisarray* pSingleSite ){

  // Set initial CNRG array (N=0) for the Kondo model

  double JK=Params[0];
  //double Lambda=Params[1];
  //double HalfLambdaFactor=Params[2];

  double auxEn=0.0;

  pAeig->ClearAll();
  pQm1fNQ->ClearAll();

  pAeig->NQNumbers=2;
  pAeig->Nshell=0;
  

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

  auxEn=JK*Sdots_totalS(0.5, 0.5, 0.5, 0.0);
  pAeig->dEn.push_back(auxEn);
  auxEn=JK*Sdots_totalS(0.0, 0.0, 0.5, 0.5);
  pAeig->dEn.push_back(auxEn);
  auxEn=JK*Sdots_totalS(1.0, 1.0, 0.5, 0.5);
  pAeig->dEn.push_back(auxEn);
  auxEn=JK*Sdots_totalS(0.5, 0.5, 0.5, 0.0);
  pAeig->dEn.push_back(auxEn);


  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);


  // Set REDUCED Matrix elements <Q S||f||Q' S'> = <Q' S'||f+||Q S> 


  // Test Kondo model
  vector <int> Indexes;
  double MatEl=0.0;

  pQm1fNQ->SyncNRGarray(*pAeig);
 
  pQm1fNQ->MatBlockMap.push_back(0); // <-1 0.5||f0||0 0>=<0 0||f+0||-1 0.5>
  pQm1fNQ->MatBlockMap.push_back(1); // |0 0>=|0 0.5>x|0 0.5 -0.5> (stsite=2)
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
  pQm1fNQ->MatEl.push_back(MatEl);
  pQm1fNQ->MatBlockBegEnd.push_back(0);
  pQm1fNQ->MatBlockBegEnd.push_back(0);



  pQm1fNQ->MatBlockMap.push_back(0); // <-1 0.5||f0||0 1>=<0 1||f+0||-1 0.5>
  pQm1fNQ->MatBlockMap.push_back(2); // |0 1>=|0 0.5>x|0 0.5 0.5> (stsite=1)
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
  pQm1fNQ->MatEl.push_back(MatEl);
  pQm1fNQ->MatBlockBegEnd.push_back(1);
  pQm1fNQ->MatBlockBegEnd.push_back(1);




  pQm1fNQ->MatBlockMap.push_back(1); // <0 0||f0||1 0.5>=<1 0.5||f+0||0 0>
  pQm1fNQ->MatBlockMap.push_back(3); // |1 0.5>=|0 0.5>x|1 0 0> (stsite=3)
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
  pQm1fNQ->MatEl.push_back(MatEl);
  pQm1fNQ->MatBlockBegEnd.push_back(2);
  pQm1fNQ->MatBlockBegEnd.push_back(2);




  pQm1fNQ->MatBlockMap.push_back(2); // <0 1||f0||1 0.5>=<1 0.5||f+0||0 1>
  pQm1fNQ->MatBlockMap.push_back(3); // |1 0.5>=|0 0.5>x|1 0 0> (stsite=3)
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
  pQm1fNQ->MatEl.push_back(MatEl);
  pQm1fNQ->MatBlockBegEnd.push_back(3);
  pQm1fNQ->MatBlockBegEnd.push_back(3);



//   for (int imbl=0;imbl<pQm1fNQ->NumMatBlocks();imbl++)
//     {
//       pQm1fNQ->PrintMatBlock(imbl);
//     }



}


///////////////
///////////////

void OneChQS_SetH0Chain(CNRGarray* pAeig, 
			CNRGmatrix* pQm1fNQ){


  pAeig->NQNumbers=2;
  pAeig->Nshell=0;

  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.5);
  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(0.0);

  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(1);
  pAeig->BlockBegEnd.push_back(2);pAeig->BlockBegEnd.push_back(2);

  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);
  pAeig->dEn.push_back(0.0);

  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(1.0);


  //fN

  pQm1fNQ->SyncNRGarray(*pAeig);
 
  pQm1fNQ->MatEl.push_back(1.0);
  pQm1fNQ->MatBlockMap.push_back(0);
  pQm1fNQ->MatBlockMap.push_back(1);
  pQm1fNQ->MatBlockBegEnd.push_back(0);
  pQm1fNQ->MatBlockBegEnd.push_back(0);

  pQm1fNQ->MatEl.push_back(-sqrt(2.0));
  pQm1fNQ->MatBlockMap.push_back(1);
  pQm1fNQ->MatBlockMap.push_back(2);
  pQm1fNQ->MatBlockBegEnd.push_back(1);
  pQm1fNQ->MatBlockBegEnd.push_back(1);



}

///////////////////////////


void OneChQS_SetSMM_Hm1(vector<double> Params,
		       CNRGarray* pAeig, 
		       CNRGmatrix* pQm1fNQ,
		       CNRGmatrix* NRGMats){

  // Set initial CNRG array (N=-1)

  double U=Params[0];
  double ed=Params[1];
  double Lambda=Params[2];
  double HalfLambdaFactor=Params[3];

  double Dz=0.0; //Params[4];
  double B2=0.0; //Params[5];
  double BmagPar=0.0; //Params[6];
  double BmagPerp=0.0; //Params[7];
  
  double U2=0.0; //Params[8];
  double ed2=0.0; //Params[9];
  double JFM=0.0; //Params[10];
  

  // First Test: the usual Anderson model

  pAeig->NQNumbers=2;
  pAeig->Nshell=-1;
  

  pAeig->QNumbers.push_back(-1.0);
  pAeig->QNumbers.push_back(0.0);

  // Two states in this one
  pAeig->QNumbers.push_back(0.0);
  pAeig->QNumbers.push_back(0.5);

  pAeig->QNumbers.push_back(1.0);
  pAeig->QNumbers.push_back(0.0);

  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

//   pAeig->dEn.push_back(0.0);
//   pAeig->dEn.push_back(ed/(Lambda*HalfLambdaFactor));
//   pAeig->dEn.push_back((2.0*ed+U)/(Lambda*HalfLambdaFactor));

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed-BmagPerp+0.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+BmagPerp+0.5*U)/(Lambda*HalfLambdaFactor));
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

  //fN

  pQm1fNQ->SyncNRGarray(*pAeig);

  pQm1fNQ->UpperTriangular=false;
 
  pQm1fNQ->MatEl.push_back(1.0);
  pQm1fNQ->MatEl.push_back(0.0);
  pQm1fNQ->MatBlockMap.push_back(0);
  pQm1fNQ->MatBlockMap.push_back(1);
  pQm1fNQ->MatBlockBegEnd.push_back(0);
  pQm1fNQ->MatBlockBegEnd.push_back(1);

  pQm1fNQ->MatEl.push_back(-1.0);
  pQm1fNQ->MatEl.push_back(0.0);
  pQm1fNQ->MatBlockMap.push_back(1);
  pQm1fNQ->MatBlockMap.push_back(2);
  pQm1fNQ->MatBlockBegEnd.push_back(2);
  pQm1fNQ->MatBlockBegEnd.push_back(3);



  // fN_up
  // Will this work?  
  NRGMats[0]=*pQm1fNQ;
  NRGMats[0].UpperTriangular=false;

  // fN_dn
  NRGMats[1].SyncNRGarray(*pAeig);

  NRGMats[1].UpperTriangular=false;
 
  NRGMats[1].MatEl.push_back(0.0);
  NRGMats[1].MatEl.push_back(1.0);
  NRGMats[1].MatBlockMap.push_back(0);
  NRGMats[1].MatBlockMap.push_back(1);
  NRGMats[1].MatBlockBegEnd.push_back(0);
  NRGMats[1].MatBlockBegEnd.push_back(1);

  NRGMats[1].MatEl.push_back(0.0);
  NRGMats[1].MatEl.push_back(1.0);
  NRGMats[1].MatBlockMap.push_back(1);
  NRGMats[1].MatBlockMap.push_back(2);
  NRGMats[1].MatBlockBegEnd.push_back(2);
  NRGMats[1].MatBlockBegEnd.push_back(3);


}
///////////////////
