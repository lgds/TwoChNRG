#include <iostream>

#include <vector>
using namespace std;

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

void DM_NRG_SetRhoNmax(vector<double> ParamsTemp,
		       CNRGbasisarray* pAcutNp1,
		       CNRGmatrix* pRhoNmax){

  cout << " DM-NRG: Calculating Rho Nmax ..." << endl;


  pRhoNmax->ClearAll();

  pRhoNmax->SyncNRGarray(*pAcutNp1);
  pRhoNmax->UpperTriangular=true;
  // Check for complex matrices
  if (pAcutNp1->CheckComplex()){pRhoNmax->IsComplex=true;}

  double TempBar=ParamsTemp[0];
  double betabar=0.0;

  bool ZeroTemp=false;

  double ThisZN[2]={0.0,0.0};

  if (dEqual(TempBar,0.0)){
    ZeroTemp=true;
    betabar=1.0e10;
  }
  else
    betabar=1.0/TempBar;

  // Partition function only for Q Sz symmetry

  double ZN=ParamsTemp[1];


  // PartitionFuncTeq0() does not work well... Use the other one.
//   if (ZeroTemp)
//     ZN=pAcutNp1->PartitionFuncTeq0();
//   else
//     ZN=pAcutNp1->PartitionFunc(betabar);

  ThisZN[0]=pAcutNp1->PartitionFunc(betabar);


  cout << " DM_NRG_SetRhoNmax: TempBar = " << TempBar
       << " betabar = " << betabar 
       << " ZN = " << ZN << " at Nmax." 
       << " Matrix is complex? " << pRhoNmax->IsComplex 
       << " Zero Temperature? " << ZeroTemp 
       << endl;

  // Set Matrix elements
  
  int i1=0;
  for (int ibl=0;ibl<pAcutNp1->NumBlocks();ibl++){
    // Rho is block diagonal
    pRhoNmax->MatBlockMap.push_back(ibl);
    pRhoNmax->MatBlockMap.push_back(ibl);

    pRhoNmax->MatBlockBegEnd.push_back(i1);
    
    // NO NEED! WRONG!
//     double Si=0.0;
//     if (pAcutNp1->totalS){
//       Si=pAcutNp1->GetQNumber(ibl,pAcutNp1->Sqnumbers[0]); 
//       // only a single SU(2) for now
//     }
//     //Debug!!
//     Si=0.0; //Should be wrong... But, apparently, it's not!!!

    int ist0=pAcutNp1->GetBlockLimit(ibl,0);
    int ist1=pAcutNp1->GetBlockLimit(ibl,1);

    for (int ist=ist0;ist<=ist1;ist++){
      double auxExp=0.0;
      auxExp=exp(-betabar*pAcutNp1->dEn[ist]);
      ThisZN[1]+=auxExp;

      // THIS IS WRONG!!
      //auxExp=(2.0*Si+1.0)*exp(-betabar*pAcutNp1->dEn[ist]);

      // debug
      // if ( (dNEqual(auxExp,0.0)) ){
      // 	cout << " DMNRG_SetRhoNmax: non-zero element "
      // 	     << " ibl= " << ibl
      // 	     << " ist= " << ist
      // 	     << " betabar= " << betabar
      // 	     << " Ei= " << pAcutNp1->dEn[ist]
      // 	     << " auxExp= " << auxExp
      // 	     << endl;
      // }
      // end debug

      if (pRhoNmax->IsComplex)
	pRhoNmax->MatElCplx.push_back(auxExp/ZN);
      else
	pRhoNmax->MatEl.push_back(auxExp/ZN);

      i1++;
      for (int jst=ist+1;jst<=ist1;jst++){ 
	//pRhoNmax->MatEl.push_back(0.0);
	if (pRhoNmax->IsComplex)
	  pRhoNmax->MatElCplx.push_back(ZeroC);
	else
	  pRhoNmax->MatEl.push_back(0.0);

	i1++;
      }
      // loop over jst
    }
    // loop over ist
    pRhoNmax->MatBlockBegEnd.push_back(i1-1);
  }
  // end loop over blocks


  cout << " ... Rho Nmax done. ThisZN[1] = " << ThisZN[1] 
       << " should be " << ThisZN[0] << endl;

}
// end subroutine

