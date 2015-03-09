#include <iostream>

#include <vector>
using namespace std;

#include <boost/timer.hpp>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "DM_NRG.hpp"


void DM_NRG_CalcRhoN_old(CNRGbasisarray* pAcutN,
		     CNRGbasisarray* pAcutNp1,
		     CNRGbasisarray* pAbasisNp1,
		     CNRGmatrix* pRhoN,
		     CNRGmatrix* pRhoNp1){

  // Poor old routine, which does not use uBLAS... (slow old man)

  vector< vector<int> > ChildSt_SameType;
  vector< vector<int> >::iterator arrayit;
  vector<int>::iterator iit;

  // Check time
  boost::timer MyTime;
  double time_elapsed;



  cout << " DM-NRG: Calculating Rho N ..." << endl;

  // Initial set up
  pRhoN->ClearAll();
  pRhoN->SyncNRGarray(*pAcutN);
  pRhoN->UpperTriangular=true;

  // Calculate matrix elements
  int i1=0;
  for (int ibl=0;ibl<pAcutN->NumBlocks();ibl++){
  // Testing in a few blocks
  //for (int ibl=0;ibl<=8;ibl+=8){
    // Rho is block diagonal
    pRhoN->MatBlockMap.push_back(ibl);
    pRhoN->MatBlockMap.push_back(ibl);
    pRhoN->MatBlockBegEnd.push_back(i1);
    
    int ist0=pAcutN->GetBlockLimit(ibl,0);
    int ist1=pAcutN->GetBlockLimit(ibl,1);

    // testing...


    int BlSize=pAcutN->GetBlockSize(ibl);
    cout << " Setting Block " << ibl << " of " << pAcutN->NumBlocks()-1
	 << "  Size : " << BlSize << " ... "; 

    MyTime.restart();
    // Loop in block states (and calculate rhoN(ist,jst))
    for (int ist=ist0;ist<=ist1;ist++){
      for (int jst=ist;jst<=ist1;jst++){

	double auxMatEl=0.0;
	// Debugging
// 	if ( (ist==61)&&(jst==61) ){
// 	  cout << " ist = " << ist << " jst = " << jst << endl;
// 	}

	DM_NRG_SetChildSameType(ist,jst,pAcutN,pAbasisNp1,ChildSt_SameType);

// 	if ( (ist==61)&&(jst==61) ){
// 	  cout << "bl  ist_Np1  jst_Np1 " << endl;
// 	  for (arrayit=ChildSt_SameType.begin();
// 	       arrayit<ChildSt_SameType.end();arrayit++){
// 	    for(iit=(*arrayit).begin();iit<(*arrayit).end();iit++){
// 	      cout << (*iit) << " ";
// 	    }
// 	    cout << endl;
// 	  }
// 	}
	// debug

	// Loop over blocks in ChildSt_SameType (N+1)

	for (arrayit=ChildSt_SameType.begin();
	     arrayit<ChildSt_SameType.end();arrayit++){
	  iit=(*arrayit).begin();
	  int iblbasis_Np1=(*iit);
	  int istbasis_Np1=(*(iit+1));
	  int jstbasis_Np1=(*(iit+2));

// 	  if ( (ist==61)&&(jst==61) ){
// 	    cout << " iblbasis_Np1 = " << iblbasis_Np1
// 		 << " istbasis_Np1 = "<< istbasis_Np1 
// 		 << " jstbasis_Np1 = "<< jstbasis_Np1 
// 		 << endl;
// 	  }
	  //
	  // Find corresponding block in AcutNp1 (compare QNumbers)
	  //
	  double* qnums=new double [pAbasisNp1->NQNumbers];
	  for (int iqn=0;iqn<pAbasisNp1->NQNumbers;iqn++){
	    qnums[iqn]=pAbasisNp1->GetQNumber(iblbasis_Np1,iqn);
	  }
	  int ibl_Np1=pAcutNp1->GetBlockFromQNumbers(qnums);
// 	  if (ibl_Np1==-1){
// 	    cout << " Block  not found in Acut N+1: ";
// 	    cout << " QNs = " ;
// 	    for (int iqn=0;iqn<pAbasisNp1->NQNumbers;iqn++){
// 	      cout << qnums[iqn] << " ";
// 	    } 
// 	    cout << endl;
// 	  }
	  // if block not found
	  delete [] qnums;

// 	  if ( (ist==61)&&(jst==61) ){
// 	    cout << " ibl_Np1 = " << ibl_Np1 << endl;
// 	    if (ibl_Np1>0)
// 	      cout << " states: from " << pAcutNp1->GetBlockLimit(ibl_Np1,0)
// 		   << " to " << pAcutNp1->GetBlockLimit(ibl_Np1,1) << endl;
// 	  }
	  // debug
	  // Loop over states in the block
	  if (ibl_Np1>0){ // Block exists in cut basis (kept)
	    for(int ist_Np1=pAcutNp1->GetBlockLimit(ibl_Np1,0);
		ist_Np1<=pAcutNp1->GetBlockLimit(ibl_Np1,1);
		ist_Np1++){
	      for(int jst_Np1=pAcutNp1->GetBlockLimit(ibl_Np1,0);
		  jst_Np1<=pAcutNp1->GetBlockLimit(ibl_Np1,1);
		  jst_Np1++){
	       double U_wi=pAcutNp1->GetEigVecComponent(ist_Np1,istbasis_Np1);
	       double U_wpj=pAcutNp1->GetEigVecComponent(jst_Np1,jstbasis_Np1);
	       double rho_wwp=pRhoNp1->GetMatEl(ist_Np1,jst_Np1);

	       auxMatEl+=rho_wwp*U_wi*U_wpj;

// 	      if ( (ist==61)&&(jst==61) ){
// 		cout << " ist_Np1 = " << ist_Np1 
// 		     << " jst_Np1 = " << jst_Np1 
// 		     << endl; 
// 		cout << " U_wi= " << U_wi
// 		     << " U_wpj= " << U_wpj 
// 		     << " rho_wwp= " << rho_wwp
// 		     << endl; 
// 		cout << " auxMatEl = " << auxMatEl << endl;
// 	      }
	      // debug
	      }
	      // end loop in j_states in ibl_Np1
	    }
	    // end loop in i_states in ibl_Np1
	  } 
	  // end if block exists in Acut

// 	  if ( (ist==61)&&(jst==61) ){
// 	    cout << " auxMatEl = " << auxMatEl << endl;
// 	  }

	}
	// end loop in ChildSt_SameType

	pRhoN->MatEl.push_back(auxMatEl);
	i1++;
      }
      // loop over jst    
    }
    // loop over ist
    pRhoN->MatBlockBegEnd.push_back(i1-1);

    time_elapsed=MyTime.elapsed();
    cout << " ... done in " << time_elapsed << " secs." << endl;

  }
  // end loop in blocks

  cout << " ... Rho N done." << endl;


}
// end subroutine

/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////


void DM_NRG_CalcRhoN(CNRGbasisarray* pAcutN,
		     CNRGbasisarray* pAcutNp1,
		     CNRGbasisarray* pAbasisNp1,
		     CNRGmatrix* pRhoN,
		     CNRGmatrix* pRhoNp1){

  // New routine, which rotates the matrix first with uBLAS (!)
  // FASTER... like a Ferrari compared to an old Volkswagen beetle.

  vector< vector<int> > ChildSt_SameType;
  vector< vector<int> >::iterator arrayit;
  vector<int>::iterator iit;

  // Check time
  boost::timer MyTime;
  double time_elapsed;

  bool debug=false;

  CNRGmatrix RhoNp1basis;

  cout << " DM-NRG: Calculating Rho N ..." << endl;

  cout << " Rotating Rho N... " << endl;

  MyTime.restart();
  RotateMatrix(pRhoNp1,pAcutNp1,(&RhoNp1basis));
  time_elapsed=MyTime.elapsed();

  cout << " ... done in " << time_elapsed << " secs. " << endl;

  int ThisBl=0;
//   if (pAcutN->Nshell==4){
//     cout << " Rho_N=5: " << endl;
//     pRhoNp1->PrintAllBlocks();
// //     pRhoNp1->PrintMatBlock(11,11);
// //     RhoNp1basis.PrintMatBlock(11,11);
//     ThisBl=11;
//     cout << " ibl = " << ThisBl << endl;
//     cout << " pAcutNp1 block structure: " << endl;
//     pAcutNp1->PrintBlockQNumbers(ThisBl);
//     cout << " Block " << ThisBl << " goes from ist = " 
// 	 << pAcutNp1->GetBlockLimit(ThisBl,0) 
// 	 << " to "
// 	 << pAcutNp1->GetBlockLimit(ThisBl,1)
// 	 << endl;
//     cout << " RhoNp1basis block structure: " << endl;
//     RhoNp1basis.PrintBlockQNumbers(ThisBl);
//     cout << " Block " << ThisBl << " goes from ist = " 
// 	 << RhoNp1basis.GetBlockLimit(ThisBl,0) 
// 	 << " to "
// 	 << RhoNp1basis.GetBlockLimit(ThisBl,1)
// 	 << endl;
//     ThisBl=18;
//     cout << " pAbasisNp1 block structure: " << endl;
// //     pAbasisNp1->PrintBlockQNumbers(ThisBl);
// //     cout << " Block " << ThisBl << " goes from ist = " 
// // 	 << pAbasisNp1->GetBlockLimit(ThisBl,0) 
// // 	 << " to "
// // 	 << pAbasisNp1->GetBlockLimit(ThisBl,1)
// // 	 << endl;
// //     pAbasisNp1->PrintBlockBasis(ThisBl);
//   }
  
  // Initial set up
  pRhoN->ClearAll();
  pRhoN->SyncNRGarray(*pAcutN);
  pRhoN->UpperTriangular=true;
  pRhoN->IsComplex=pRhoNp1->IsComplex;

  // Calculate matrix elements
  int i1=0;
  for (int ibl=0;ibl<pAcutN->NumBlocks();ibl++){


    // Rho is block diagonal
    pRhoN->MatBlockMap.push_back(ibl);
    pRhoN->MatBlockMap.push_back(ibl);
    pRhoN->MatBlockBegEnd.push_back(i1);
    
    int ist0=pAcutN->GetBlockLimit(ibl,0);
    int ist1=pAcutN->GetBlockLimit(ibl,1);

    // testing...

    // Get Sold from ibl



    int BlSize=pAcutN->GetBlockSize(ibl);
    if (debug)
      cout << " Setting Block " << ibl << " of " << pAcutN->NumBlocks()-1
	   << "  Size : " << BlSize << " ... "; 

    MyTime.restart();
    // Loop in block states (and calculate rhoN(ist,jst))
    for (int ist=ist0;ist<=ist1;ist++){
      for (int jst=ist;jst<=ist1;jst++){


	// Ok, watch out for the " type" thing for SU(2) symmetry

	DM_NRG_SetChildSameType(ist,jst,pAcutN,pAbasisNp1,ChildSt_SameType);

	double auxMatEl=0.0;
	complex<double> cMatEl=ZeroC;
	// Loop over blocks in ChildSt_SameType (N+1)
	for (arrayit=ChildSt_SameType.begin();
	     arrayit<ChildSt_SameType.end();arrayit++){
	  iit=(*arrayit).begin();
	  int iblbasis_Np1=(*iit);
	  int istbasis_Np1=(*(iit+1));
	  int jstbasis_Np1=(*(iit+2));

	  // Get S_Np1 from iblbasis_Np1
	  // Get S_type= from istbasis_Np1
	  // Calculate Clebsh-Gordan coefs.


	  //
	  // No need to find corresponding block in AcutNp1:
	  //
	  // if block is not there, then 
	  // RhoNp1basis.GetMatEl(istbasis_Np1,jstbasis_Np1)
	  // will return 0.0 (in principle!)
	  // 

	  if (!(RhoNp1basis.IsComplex))
	    auxMatEl+=RhoNp1basis.GetMatEl(istbasis_Np1,jstbasis_Np1);
	  else
	    cMatEl+=RhoNp1basis.cGetMatEl(istbasis_Np1,jstbasis_Np1);
	    

//  	  if ( (debug)&&(iblbasis_Np1==18)&&(ist==jst) ){
// 	    cout << " ibl = " << ibl  
// 		 << " ist = " << ist
// 		 << " jst = " << jst
// 		 << " iblbasis_Np1 = " << iblbasis_Np1
// 		 << " istbasis_Np1 = " << istbasis_Np1
// 		 << " jstbasis_Np1 = " << jstbasis_Np1
// 		 << " mat el = " << RhoNp1basis.GetMatEl(istbasis_Np1,jstbasis_Np1)
// 		 << " auxMatEl = " << auxMatEl
// 		 << endl; 
// 	  }


	}
	// end loop in Child states
	if (!(pRhoN->IsComplex))	  
	  pRhoN->MatEl.push_back(auxMatEl);
	else
	  pRhoN->MatElCplx.push_back(cMatEl);
	  
	i1++;
      }
      // end loop in jst
    }
    // end loop in ist
    pRhoN->MatBlockBegEnd.push_back(i1-1);


    time_elapsed=MyTime.elapsed();
    if (debug)
     cout << " ... done in " << time_elapsed << " secs." << endl;

  }
  // end loop in pAcutN blocks

  cout << " ... Rho N done." << endl;

//   if (pAcutN->Nshell==4){
//     cout << " Rho_N=4: " << endl;
//     //pRhoN->PrintMatBlock(10,10);
//     ThisBl=10;
//     cout << " ibl = " << ThisBl << endl;
//     cout << " pRhoN block structure: " << endl;
//     pRhoN->PrintBlockQNumbers(ThisBl);
//     cout << " Block " << ThisBl << " goes from ist = " 
// 	 << pRhoN->GetBlockLimit(ThisBl,0) 
// 	 << " to "
// 	 << pRhoN->GetBlockLimit(ThisBl,1) 
// 	 << endl;
//     //pRhoN->PrintMatBlock(13,13);
//     ThisBl=13;
//     cout << " ibl = " << ThisBl << endl;
//     cout << " pRhoN block structure: " << endl;
//     pRhoN->PrintBlockQNumbers(ThisBl);
//     cout << " Block " << ThisBl << " goes from ist = " 
// 	 << pRhoN->GetBlockLimit(ThisBl,0) 
// 	 << " to "
// 	 << pRhoN->GetBlockLimit(ThisBl,1) 
// 	 << endl;
//   }

}
// end set RhoN new


/////////////////////////////////////////
/////   Version with SU(2) symmetry /////
/////////////////////////////////////////


void DM_NRG_CalcRhoN_withSU2(CNRGbasisarray* pAcutN,
		     CNRGbasisarray* pAcutNp1,
		     CNRGbasisarray* pAbasisNp1,
		     CNRGbasisarray* pSingleSite,
		     CNRGmatrix* pRhoN,
		     CNRGmatrix* pRhoNp1){

  // New routine, which rotates the matrix first with uBLAS (!)
  // FASTER... like a Ferrari compared to an old Volkswagen beetle.

  vector< vector<int> > ChildSt_SameType;
  vector< vector<int> >::iterator arrayit;
  vector<int>::iterator iit;

  // Check time
  boost::timer MyTime;
  double time_elapsed;

  bool debug=false;

  CNRGmatrix RhoNp1basis;

  cout << " DM-NRG: Calculating Rho N (w SU(2))..." << endl;

  cout << " Rotating Rho N... " << endl;

  MyTime.restart();
  RotateMatrix(pRhoNp1,pAcutNp1,(&RhoNp1basis));
  time_elapsed=MyTime.elapsed();

  cout << " ... done in " << time_elapsed << " secs. " << endl;

  int ThisBl=0;
  
  // Initial set up
  pRhoN->ClearAll();
  pRhoN->SyncNRGarray(*pAcutN);
  pRhoN->UpperTriangular=true;

  // Calculate matrix elements
  int i1=0;
  for (int ibl=0;ibl<pAcutN->NumBlocks();ibl++){


    // Rho is block diagonal
    pRhoN->MatBlockMap.push_back(ibl);
    pRhoN->MatBlockMap.push_back(ibl);
    pRhoN->MatBlockBegEnd.push_back(i1);
    
    int ist0=pAcutN->GetBlockLimit(ibl,0);
    int ist1=pAcutN->GetBlockLimit(ibl,1);

    // testing...

    // Get Sold from ibl
    double Sold=0.0;
    if (pAcutN->totalS){
      Sold=pAcutN->GetQNumber(ibl,pAcutN->Sqnumbers[0]); 
      // only a single SU(2) for now
    }


    int BlSize=pAcutN->GetBlockSize(ibl);
    if (debug)
      cout << " Setting Block " << ibl << " of " << pAcutN->NumBlocks()-1
	   << "  Size : " << BlSize << " ... "; 

    MyTime.restart();
    // Loop in block states (and calculate rhoN(ist,jst))
    for (int ist=ist0;ist<=ist1;ist++){
      for (int jst=ist;jst<=ist1;jst++){

	// Ok, watch out for the " type" thing for SU(2) symmetry

	DM_NRG_SetChildSameType(ist,jst,pAcutN,pAbasisNp1,ChildSt_SameType);

	double auxMatEl=0.0;
	complex<double> cMatEl=ZeroC;
	// Loop over blocks in ChildSt_SameType (N+1)
	for (arrayit=ChildSt_SameType.begin();
	     arrayit<ChildSt_SameType.end();arrayit++){
	  iit=(*arrayit).begin();
	  int iblbasis_Np1=(*iit);
	  int istbasis_Np1=(*(iit+1));
	  int jstbasis_Np1=(*(iit+2));

	  // Get S_Np1 from iblbasis_Np1
	  double Si=0.0;
	  if (pAbasisNp1->totalS){
	    Si=pAbasisNp1->GetQNumber(iblbasis_Np1,pAbasisNp1->Sqnumbers[0]); 
	    // only a single SU(2) for now
	  }
	  // Get S_type= from istbasis_Np1
	  int itype_Np1=pAbasisNp1->iType[istbasis_Np1];
	  double Stilde=0.0;
	  if (pSingleSite->totalS){
	    //Stilde=pSingleSite->GetQNumber(ibls,pSingleSite->Sqnumbers[0]);
	    Stilde=pSingleSite->GetQNumberFromSt(itype_Np1,pSingleSite->Sqnumbers[0]);
	    // only a single SU(2) for now
	  }

	  // Calculate Clebsh-Gordan coefs.

	  double CGfactor=0.0;
	  for (double Sz=-Si;Sz<=Si;Sz+=1.0){
	    double auxCG=CGordan(Sold,Sold,Stilde,Sz-Sold,Si,Sz);
	    CGfactor+=auxCG*auxCG; // Square it!! (oct 2010)
	  }
	  // end calcCG

	  //
	  // No need to find corresponding block in AcutNp1:
	  //
	  // if block is not there, then 
	  // RhoNp1basis.GetMatEl(istbasis_Np1,jstbasis_Np1)
	  // will return 0.0 (in principle!)
	  // 
	  //auxMatEl+=CGfactor*RhoNp1basis.GetMatEl(istbasis_Np1,jstbasis_Np1);
	  if (!(RhoNp1basis.IsComplex))
	    auxMatEl+=CGfactor*RhoNp1basis.GetMatEl(istbasis_Np1,jstbasis_Np1);
	  else
	    cMatEl+=CGfactor*RhoNp1basis.cGetMatEl(istbasis_Np1,jstbasis_Np1);

	}
	// end loop in Child states
	if (pRhoN->IsComplex)
	  pRhoN->MatElCplx.push_back(cMatEl);
	else
	  pRhoN->MatEl.push_back(auxMatEl);

	i1++;
      }
      // end loop in jst
    }
    // end loop in ist
    pRhoN->MatBlockBegEnd.push_back(i1-1);


    time_elapsed=MyTime.elapsed();
    if (debug)
     cout << " ... done in " << time_elapsed << " secs." << endl;

  }
  // end loop in pAcutN blocks

  cout << " ... Rho N done." << endl;

//   if (pAcutN->Nshell==4){
//     cout << " Rho_N=4: " << endl;
//     //pRhoN->PrintMatBlock(10,10);
//     ThisBl=10;
//     cout << " ibl = " << ThisBl << endl;
//     cout << " pRhoN block structure: " << endl;
//     pRhoN->PrintBlockQNumbers(ThisBl);
//     cout << " Block " << ThisBl << " goes from ist = " 
// 	 << pRhoN->GetBlockLimit(ThisBl,0) 
// 	 << " to "
// 	 << pRhoN->GetBlockLimit(ThisBl,1) 
// 	 << endl;
//     //pRhoN->PrintMatBlock(13,13);
//     ThisBl=13;
//     cout << " ibl = " << ThisBl << endl;
//     cout << " pRhoN block structure: " << endl;
//     pRhoN->PrintBlockQNumbers(ThisBl);
//     cout << " Block " << ThisBl << " goes from ist = " 
// 	 << pRhoN->GetBlockLimit(ThisBl,0) 
// 	 << " to "
// 	 << pRhoN->GetBlockLimit(ThisBl,1) 
// 	 << endl;
//   }

}
// end set RhoN new
