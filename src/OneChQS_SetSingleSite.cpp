#include "NRGclasses.hpp"
#include "OneChQS.hpp"


int OneChQS_SetSingleSite(CNRGbasisarray &SingleSite){

  // one channel, Q,S,Sz basis:


  SingleSite.NQNumbers=3;
  SingleSite.Nshell=0;
  // new (April 2010)
  SingleSite.totalS=true;
  SingleSite.Sqnumbers.push_back(1);

  
  // QS blocks
  SingleSite.QNumbers.push_back(-1.0);
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.0);

  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(0.5);

  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(-0.5);

  SingleSite.QNumbers.push_back(1.0);
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.0);

  // One state per block
  for (int ii=0;ii<4;ii++){
    SingleSite.BlockBegEnd.push_back(ii);
    SingleSite.BlockBegEnd.push_back(ii);
    // Type labels the state
    SingleSite.iType.push_back(ii);
  }

  // Important: SingleSite has Q,S,Sz block structure, one state per block

}
