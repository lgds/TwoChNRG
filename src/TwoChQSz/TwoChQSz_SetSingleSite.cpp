
#include "NRGclasses.hpp"


int TwoChQSz_SetSingleSite(CNRGbasisarray* pSingleSite){


  // two channels: Q,Sz basis:


  pSingleSite->NQNumbers=2;
  pSingleSite->Nshell=0;
  
  // 16 states:
  // |0, 0>, |0, up>, |0, dn>, |0, up dn>

  pSingleSite->QNumbers.push_back(-2.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(-0.5);

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);

  // |up, 0>, |up, up>, |up, dn>, |up, up dn>
  
  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);

  // |dn, 0>, |dn, up>, |dn, dn>, |dn, up dn>
  
  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(-0.5);

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(-0.5);


  // |up dn, 0>, |up dn, up>, |up dn, dn>, |up dn, up dn>
  
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(-0.5);

  pSingleSite->QNumbers.push_back(2.0);
  pSingleSite->QNumbers.push_back(0.0);



  // One state per block
  for (int ii=0;ii<16;ii++){
    pSingleSite->BlockBegEnd.push_back(ii);
    pSingleSite->BlockBegEnd.push_back(ii);
    // Type labels the state
    pSingleSite->iType.push_back(ii);
  }


  return(0);

}
