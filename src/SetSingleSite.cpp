#include "NRGclasses.hpp"

//
// Note: as of Sept. 09, this file is yet not being called by
// the main code (NRG_main.cpp) and does not have a header file either.
//
// Currently, SetSingleSite routines are in Main_SetSingleFile.cpp
//
// The idea is to eventually concentrate all SetSingleSite routines here.
//
// Ops, DM-NRG uses this!
//

int TwoChQSNoSz_SetSingleSite(CNRGbasisarray* pSingleSite){

  pSingleSite->NQNumbers=2;
  pSingleSite->Nshell=0;
  // new
  pSingleSite->totalS=true;
  pSingleSite->Sqnumbers.push_back(1);


  // 16 states divided in 6 Q,S blocks
  // with Sz entering in iDegen (Sz=10*iDegen)


  // -2 0 : | 0 0>
  pSingleSite->QNumbers.push_back(-2.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->BlockBegEnd.push_back(0);
  pSingleSite->BlockBegEnd.push_back(0);

  pSingleSite->iDegen.push_back(0); // Sz=0

  // -1 1/2 : | 0 up> and | up 0> ; | 0 dn> and | dn 0> 

  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->BlockBegEnd.push_back(1);
  pSingleSite->BlockBegEnd.push_back(4);

  pSingleSite->iDegen.push_back(5); // Sz=0.5  
  pSingleSite->iDegen.push_back(5); // Sz=0.5
  pSingleSite->iDegen.push_back(-5); // Sz=-0.5  
  pSingleSite->iDegen.push_back(-5); // Sz=-0.5


  // 0 0 : 1/sqrt(2)(| up dn>-| dn up>); | 0 up dn>; | up dn 0>; 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0); // includes Singlet

  pSingleSite->BlockBegEnd.push_back(5);
  pSingleSite->BlockBegEnd.push_back(7);

  pSingleSite->iDegen.push_back(0); // Sz=0
  pSingleSite->iDegen.push_back(0); // Sz=0
  pSingleSite->iDegen.push_back(0); // Sz=0


  // 0 1 : | dn dn> ;  1/sqrt(2)(| up dn>+| dn up>) ;  | up up> 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0); // Triplet

  pSingleSite->BlockBegEnd.push_back(8);
  pSingleSite->BlockBegEnd.push_back(10);

  pSingleSite->iDegen.push_back(-10); // Sz=0
  pSingleSite->iDegen.push_back(0); // Sz=0
  pSingleSite->iDegen.push_back(10); // Sz=0


  // 1 1/2 : | up  up dn>; | up dn   up> ; | dn  up dn>; | up dn   dn> 

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->BlockBegEnd.push_back(11);
  pSingleSite->BlockBegEnd.push_back(14);

  pSingleSite->iDegen.push_back(5); // Sz=0.5  
  pSingleSite->iDegen.push_back(5); // Sz=0.5
  pSingleSite->iDegen.push_back(-5); // Sz=-0.5  
  pSingleSite->iDegen.push_back(-5); // Sz=-0.5


  // 2 0 : | up dn  up dn>

  pSingleSite->QNumbers.push_back(2.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->BlockBegEnd.push_back(15);
  pSingleSite->BlockBegEnd.push_back(15);

  pSingleSite->iDegen.push_back(0); // Sz=0


  // Type labels the state
  for (int ii=0;ii<16;ii++){
    pSingleSite->iType.push_back(ii);
    // New and useful...
    pSingleSite->StCameFrom.push_back(ii);
  }


  return(0);
}
////

/////////////////////////////////////
//                                 //
// Two channel:                    //
// SU(2)iso + SU(2)spin + Parity   //
//                                 //
//   I,Iz,S,Sz,P                   //
//                                 //
// (The whole enchilada)           //
//                                 //
/////////////////////////////////////

int TwoChISP_SetSingleSite(CNRGbasisarray* pSingleSite){

  pSingleSite->NQNumbers=5;
  pSingleSite->Nshell=0;

  // new. Needs check.
  pSingleSite->totalS=true;
  pSingleSite->Sqnumbers.push_back(0);
  pSingleSite->Sqnumbers.push_back(2);


  //
  // 16 states divided in 6 I,S,P blocks
  //
  // (Iz=Q/2)



  //
  // I=1 S=0 P=+1 (3 states)
  //

  // | 0 0> (I=1 Iz=-1 S=0 Sz=0 P=+1)
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(0);
  pSingleSite->BlockBegEnd.push_back(0);


  // 1/sqrt2(| up dn  0> + |0  up dn> (I=1 Iz=0 S=0 Sz=0 P=+1)
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(1);
  pSingleSite->BlockBegEnd.push_back(1);

  
  // | up dn   up dn >  (I=1 Iz=1 S=0 Sz=0 P=+1)
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(2);
  pSingleSite->BlockBegEnd.push_back(2);


  //
  // I=1/2 S=1/2 P=-1 (4 states)
  //

  // |up dn  up >  (I=1/2 Iz=+1/2 S=1/2 Sz=+1/2 P=-1) 

  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(3);
  pSingleSite->BlockBegEnd.push_back(3);

  // | 0  up >  (I=1/2 Iz=-1/2 S=1/2 Sz=+1/2 P=-1) 

  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(4);
  pSingleSite->BlockBegEnd.push_back(4);

  // |up dn  dn >  (I=1/2 Iz=+1/2 S=1/2 Sz=-1/2 P=-1) 

  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(5);
  pSingleSite->BlockBegEnd.push_back(5);

  // | 0  dn >  (I=1/2 Iz=-1/2 S=1/2 Sz=-1/2 P=-1) 

  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(6);
  pSingleSite->BlockBegEnd.push_back(6);


  //
  // I=1/2 S=1/2 P=+1 (4 states)
  //

  // |up  up dn >  (I=1/2 Iz=+1/2 S=1/2 Sz=+1/2 P=+1) 

  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(7);
  pSingleSite->BlockBegEnd.push_back(7);

  // | up  0 >  (I=1/2 Iz=-1/2 S=1/2 Sz=+1/2 P=+1) 

  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(8);
  pSingleSite->BlockBegEnd.push_back(8);

  // |dn  up dn >  (I=1/2 Iz=+1/2 S=1/2 Sz=-1/2 P=1) 

  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(9);
  pSingleSite->BlockBegEnd.push_back(9);

  // | dn  0 >  (I=1/2 Iz=-1/2 S=1/2 Sz=-1/2 P=1) 

  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(0.5);
  pSingleSite->QNumbers.push_back(-0.5);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(10);
  pSingleSite->BlockBegEnd.push_back(10);


  //
  // I=0 S=1 P=-1 (3 states)
  //

  // | dn dn > (I=0 Iz=0 S=1 Sz=-1 P=-1)

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(11);
  pSingleSite->BlockBegEnd.push_back(11);



  // 1/sqrt2(|up dn> + |dn up>) (I=0 Iz=0 S=1 Sz=0 P=-1)

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(12);
  pSingleSite->BlockBegEnd.push_back(12);


  // | up up > (I=0 Iz=0 S=1 Sz=1 P=-1)

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(13);
  pSingleSite->BlockBegEnd.push_back(13);


  //
  // I=0 S=0 P=-1 (1 state)
  //
  // 1/sqrt2(|up dn> - |dn up>) (I=0 Iz=0 S=0 Sz=0 P=-1)

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(-1.0);

  pSingleSite->BlockBegEnd.push_back(14);
  pSingleSite->BlockBegEnd.push_back(14);

  //
  // I=0 S=0 P=+1 (1 state)
  //
  // 1/sqrt2(| up dn  0> - |0  up dn> (I=0 Iz=0 S=0 Sz=0 P=+1)
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0);

  pSingleSite->BlockBegEnd.push_back(15);
  pSingleSite->BlockBegEnd.push_back(15);


  // Type labels the state
  for (int ii=0;ii<16;ii++){
    pSingleSite->iType.push_back(ii);
    // New and useful...
    pSingleSite->StCameFrom.push_back(ii);
  }


  return(0);
}
///////

int OneChQSz_SetSingleSite(CNRGbasisarray &SingleSite){

  // one channel, Q,S,Sz basis:


  SingleSite.NQNumbers=2;
  SingleSite.Nshell=0;
  // new
  SingleSite.totalS=false;

  
  // QSz blocks
  SingleSite.QNumbers.push_back(-1.0);
  SingleSite.QNumbers.push_back(0.0);

  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.5);

  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(-0.5);

  SingleSite.QNumbers.push_back(1.0);
  SingleSite.QNumbers.push_back(0.0);

  // One state per block
  for (int ii=0;ii<4;ii++){
    SingleSite.BlockBegEnd.push_back(ii);
    SingleSite.BlockBegEnd.push_back(ii);
    // Type labels the state
    SingleSite.iType.push_back(ii);
  }

  // Important: SingleSite has Q,Sz block structure, one state per block


  return(0);
}
/////////////

// NEW (similar to the one in TwoChQS dir, but better)

void TwoChQSP_SetSingleSite(CNRGbasisarray &SingleSite){


  // two channels: Q,S,Sz and P (parity) basis:


  SingleSite.NQNumbers=4;
  SingleSite.Nshell=0;
  // new
  SingleSite.totalS=true;
  SingleSite.Sqnumbers.push_back(1);

  
  // 16 states divided in 10 blocks

  // -2 0 0 P=+1 : | 0 0>
  SingleSite.QNumbers.push_back(-2.0);
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(1.0);


  SingleSite.BlockBegEnd.push_back(0);
  SingleSite.BlockBegEnd.push_back(0);

  // -1 1/2 1/2 P=-1 : | 0 up> 

  SingleSite.QNumbers.push_back(-1.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(-1.0);

  SingleSite.BlockBegEnd.push_back(1);
  SingleSite.BlockBegEnd.push_back(1);


  // -1 1/2 1/2 P=+1 : | up 0> 

  SingleSite.QNumbers.push_back(-1.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(1.0);

  SingleSite.BlockBegEnd.push_back(2);
  SingleSite.BlockBegEnd.push_back(2);


  // -1 1/2 -1/2 P=-1 : | 0 dn> 


  SingleSite.QNumbers.push_back(-1.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(-0.5);
  SingleSite.QNumbers.push_back(-1.0);

  SingleSite.BlockBegEnd.push_back(3);
  SingleSite.BlockBegEnd.push_back(3);

  // -1 1/2 -1/2 P=+1 : | dn 0> 


  SingleSite.QNumbers.push_back(-1.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(-0.5);
  SingleSite.QNumbers.push_back(1.0);

  SingleSite.BlockBegEnd.push_back(4);
  SingleSite.BlockBegEnd.push_back(4);


  // 0 0 0 P=-1 : 1/sqrt(2)(| up dn>-| dn up>); 

  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.0); // Singlet
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(-1.0);

  SingleSite.BlockBegEnd.push_back(5);
  SingleSite.BlockBegEnd.push_back(5);


  // 0 0 0 P=1 : | 0 up dn>; | up dn 0>;


  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(1.0);

  SingleSite.BlockBegEnd.push_back(6);
  SingleSite.BlockBegEnd.push_back(7);
 


  // 0 1 -1 P=-1: | dn dn> 

  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(1.0); // Triplet
  SingleSite.QNumbers.push_back(-1.0);
  SingleSite.QNumbers.push_back(-1.0);

  SingleSite.BlockBegEnd.push_back(8);
  SingleSite.BlockBegEnd.push_back(8);


  // 0 1 0 P=-1: 1/sqrt(2)(| up dn>+| dn up>) 

  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(1.0); // Triplet
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(-1.0);

  SingleSite.BlockBegEnd.push_back(9);
  SingleSite.BlockBegEnd.push_back(9);


  // 0 1 1 P=-1: | up up> 

  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(1.0); // Triplet
  SingleSite.QNumbers.push_back(1.0);
  SingleSite.QNumbers.push_back(-1.0);

  SingleSite.BlockBegEnd.push_back(10);
  SingleSite.BlockBegEnd.push_back(10);

  // 1 1/2 1/2 P=1: | up  up dn>;

  SingleSite.QNumbers.push_back(1.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(1.0);

  SingleSite.BlockBegEnd.push_back(11);
  SingleSite.BlockBegEnd.push_back(11);

  // 1 1/2 1/2 P=-1: | up dn   up> 

  SingleSite.QNumbers.push_back(1.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(-1.0);

  SingleSite.BlockBegEnd.push_back(12);
  SingleSite.BlockBegEnd.push_back(12);


  // 1 1/2 -1/2 P=1: | dn  up dn>

  SingleSite.QNumbers.push_back(1.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(-0.5);
  SingleSite.QNumbers.push_back(1.0);

  SingleSite.BlockBegEnd.push_back(13);
  SingleSite.BlockBegEnd.push_back(13);

  // 1 1/2 -1/2 P=-1: | up dn   dn> 

  SingleSite.QNumbers.push_back(1.0);
  SingleSite.QNumbers.push_back(0.5);
  SingleSite.QNumbers.push_back(-0.5);
  SingleSite.QNumbers.push_back(-1.0);

  SingleSite.BlockBegEnd.push_back(14);
  SingleSite.BlockBegEnd.push_back(14);

  // 2 0 0 P=+1: | up dn  up dn>

  SingleSite.QNumbers.push_back(2.0);
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(0.0);
  SingleSite.QNumbers.push_back(1.0);

  SingleSite.BlockBegEnd.push_back(15);
  SingleSite.BlockBegEnd.push_back(15);


  // Type labels the state
  for (int ii=0;ii<16;ii++){
    SingleSite.iType.push_back(ii);
  }

  SingleSite.ChildStates.push_back( vector<int>());
  SingleSite.ChildStates.clear();  

  return;

}

/////////////

void TwoDotQS_SetInitialSite(CNRGbasisarray* pSingleSite){


  // two dots: Q,S basis:
  // (only states with Sz=S are included!)

  pSingleSite->NQNumbers=2;
  pSingleSite->Nshell=0;
  // new
  pSingleSite->totalS=true;
  pSingleSite->Sqnumbers.push_back(1);

  
  // 10 states divided in 6 blocks

  // -2 0 0 : | 0 0>
  pSingleSite->QNumbers.push_back(-2.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->BlockBegEnd.push_back(0);
  pSingleSite->BlockBegEnd.push_back(0);

  // -1 1/2 1/2 : | 0 up> and | up 0> 
  // -1 1/2 -1/2 : | 0 dn> and | dn 0> not included

  pSingleSite->QNumbers.push_back(-1.0);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->BlockBegEnd.push_back(1);
  pSingleSite->BlockBegEnd.push_back(2);


  // 0 0 0 : 1/sqrt(2)(| up dn>-| dn up>); | 0 up dn>; | up dn 0>; 

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(0.0); // includes Singlet
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->BlockBegEnd.push_back(3);
  pSingleSite->BlockBegEnd.push_back(5);


  // 0 1 1 : | up up> 

  // 0 1 -1 : | dn dn> not included
  // 0 1 0 : 1/sqrt(2)(| up dn>+| dn up>) not included

  pSingleSite->QNumbers.push_back(0.0);
  pSingleSite->QNumbers.push_back(1.0); // Triplet

  pSingleSite->BlockBegEnd.push_back(6);
  pSingleSite->BlockBegEnd.push_back(6);

  // 1 1/2 1/2 : | up  up dn> and | up dn   up> 
  // 1 1/2 -1/2 : | dn  up dn>; | up dn   dn> not included

  pSingleSite->QNumbers.push_back(1.0);
  pSingleSite->QNumbers.push_back(0.5);

  pSingleSite->BlockBegEnd.push_back(7);
  pSingleSite->BlockBegEnd.push_back(8);

  // 2 0 0 : | up dn  up dn>

  pSingleSite->QNumbers.push_back(2.0);
  pSingleSite->QNumbers.push_back(0.0);

  pSingleSite->BlockBegEnd.push_back(9);
  pSingleSite->BlockBegEnd.push_back(9);


  // Type labels the state
  for (int ii=0;ii<10;ii++){
    pSingleSite->iType.push_back(ii);
  }


}

///////////////////////////
