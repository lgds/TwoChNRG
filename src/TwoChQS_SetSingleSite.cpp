
#include "NRGclasses.hpp"


int TwoChQS_SetSingleSite(CNRGbasisarray* pSingleSite){


  // two channels: Q,S,Sz basis:


  pSingleSite->NQNumbers=3;
  pSingleSite->Nshell=0;
  // new (April 2010)
  pSingleSite->totalS=true;
  pSingleSite->Sqnumbers.push_back(1);

  
  // 16 states divided in 10 blocks

  // -2 0 0 : | 0 0>
  pSingleSite->QNumbers.push_back(-2.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);


  pSingleSite->BlockBegEnd.push_back(0);
  pSingleSite->BlockBegEnd.push_back(0);

  // -1 1/2 1/2 : | 0 up> and | up 0> 

  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->BlockBegEnd.push_back(1);
  pSingleSite->BlockBegEnd.push_back(2);

  // -1 1/2 -1/2 : | 0 dn> and | dn 0> 


  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);

  pSingleSite->BlockBegEnd.push_back(3);
  pSingleSite->BlockBegEnd.push_back(4);

  // 0 0 0 : 1/sqrt(2)(| up dn>-| dn up>); | 0 up dn>; | up dn 0>; 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0); // includes Singlet
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->BlockBegEnd.push_back(5);
  pSingleSite->BlockBegEnd.push_back(7);


  // 0 1 -1 : | dn dn> 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0); // Triplet
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(8);
  pSingleSite->BlockBegEnd.push_back(8);


  // 0 1 0 : 1/sqrt(2)(| up dn>+| dn up>) 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0); // Triplet
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->BlockBegEnd.push_back(9);
  pSingleSite->BlockBegEnd.push_back(9);


  // 0 1 1 : | up up> 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0); // Triplet
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(10);
  pSingleSite->BlockBegEnd.push_back(10);

  // 1 1/2 1/2 : | up  up dn>; | up dn   up> 

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->BlockBegEnd.push_back(11);
  pSingleSite->BlockBegEnd.push_back(12);

  // 1 1/2 -1/2 : | dn  up dn>; | up dn   dn> 

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);

  pSingleSite->BlockBegEnd.push_back(13);
  pSingleSite->BlockBegEnd.push_back(14);


  // 2 0 0 : | up dn  up dn>

  pSingleSite->QNumbers.push_back(2.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->BlockBegEnd.push_back(15);
  pSingleSite->BlockBegEnd.push_back(15);


  // Type labels the state
  for (int ii=0;ii<16;ii++){
    pSingleSite->iType.push_back(ii);
  }


  return(0);

}


////////////////////////////////////////////////////

int TwoChQSP_SetSingleSite(CNRGbasisarray* pSingleSite){


  // two channels: Q,S,Sz and P (parity) basis:


  pSingleSite->NQNumbers=4;
  pSingleSite->Nshell=0;
  
  // 16 states divided in 10 blocks

  // -2 0 0 P=+1 : | 0 0>
  pSingleSite->QNumbers.push_back(-2.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);


  pSingleSite->BlockBegEnd.push_back(0);
  pSingleSite->BlockBegEnd.push_back(0);

  // -1 1/2 1/2 P=-1 : | 0 up> 

  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(1);
  pSingleSite->BlockBegEnd.push_back(1);


  // -1 1/2 1/2 P=+1 : | up 0> 

  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(2);
  pSingleSite->BlockBegEnd.push_back(2);


  // -1 1/2 -1/2 P=-1 : | 0 dn> 


  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(3);
  pSingleSite->BlockBegEnd.push_back(3);

  // -1 1/2 -1/2 P=+1 : | dn 0> 


  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(4);
  pSingleSite->BlockBegEnd.push_back(4);


  // 0 0 0 P=-1 : 1/sqrt(2)(| up dn>-| dn up>); 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0); // Singlet
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(5);
  pSingleSite->BlockBegEnd.push_back(5);


  // 0 0 0 P=1 : | 0 up dn>; | up dn 0>;


  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(6);
  pSingleSite->BlockBegEnd.push_back(7);
 


  // 0 1 -1 P=-1: | dn dn> 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0); // Triplet
  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(8);
  pSingleSite->BlockBegEnd.push_back(8);


  // 0 1 0 P=-1: 1/sqrt(2)(| up dn>+| dn up>) 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0); // Triplet
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(9);
  pSingleSite->BlockBegEnd.push_back(9);


  // 0 1 1 P=-1: | up up> 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0); // Triplet
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(10);
  pSingleSite->BlockBegEnd.push_back(10);

  // 1 1/2 1/2 P=1: | up  up dn>;

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(11);
  pSingleSite->BlockBegEnd.push_back(11);

  // 1 1/2 1/2 P=-1: | up dn   up> 

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(12);
  pSingleSite->BlockBegEnd.push_back(12);


  // 1 1/2 -1/2 P=1: | dn  up dn>

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(13);
  pSingleSite->BlockBegEnd.push_back(13);

  // 1 1/2 -1/2 P=-1: | up dn   dn> 

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(14);
  pSingleSite->BlockBegEnd.push_back(14);

  // 2 0 0 P=+1: | up dn  up dn>

  pSingleSite->QNumbers.push_back(2.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(15);
  pSingleSite->BlockBegEnd.push_back(15);


  // Type labels the state
  for (int ii=0;ii<16;ii++){
    pSingleSite->iType.push_back(ii);
  }


  return(0);

}
