#include <iostream>

#include <vector>
using namespace std;

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"



void DM_NRG_SetChildSt(CNRGbasisarray* pAcutN,
		   CNRGbasisarray* pAbasisNp1){

  vector<int>::iterator iit;

  cout << " DM-NRG: Setting child states... " << endl;

  pAcutN->ChildStates.clear();

  // Add Nstates empty vectors to AcutN.ChildStates (lines)

  for (int ist=0;ist<pAcutN->Nstates();ist++){
    pAcutN->ChildStates.push_back( vector<int>() );
  }

  // Add Columns to AcutN.ChildStates (lines)
  
  int icount=0;
  for (iit=pAbasisNp1->StCameFrom.begin();
       iit<pAbasisNp1->StCameFrom.end();iit++){
    int stcf=*iit;
    pAcutN->ChildStates[stcf].push_back(icount);
    icount++;
  }

  cout << " ...Set child states done." << endl;

}
// end routine

///
/// Next: Given ist and jst in N and the N+1 
///       child vectors, set up a 
///       STL array of the for ibl_child ist_child jst_child
///       for pairt (ist_child, jst_child) with the SAME type

void DM_NRG_SetChildSameType(int ist_N, int jst_N,
			     CNRGbasisarray* pAcutN, 
			     CNRGbasisarray* pAbasisNp1,
			     vector< vector<int> > &ChildSt_SameType){

  vector< vector<int> >::iterator arrayit;
  vector<int>::iterator iit,iit2;

  int ist_Np1,jst_Np1;
  int itype_Np1,jtype_Np1;
  int ibl_Np1,jbl_Np1;
  int pos_bl;
  

  ChildSt_SameType.clear();

  if( pAcutN->ChildStates.size()<(ist_N>jst_N?ist_N:jst_N) ){
    cout << "ChildStates too short" << endl;
    return;
  }

  int ifound=0;
  for(iit=pAcutN->ChildStates[ist_N].begin();
      iit<pAcutN->ChildStates[ist_N].end();iit++){
    ist_Np1=*iit;
    itype_Np1=pAbasisNp1->iType[ist_Np1];
    ibl_Np1=pAbasisNp1->GetBlockFromSt(ist_Np1,pos_bl);

//     cout << " ibl= " << ibl_Np1
//       	 << " itype_Np1= " << itype_Np1
// 	 << " ist_Np1= " << ist_Np1
// 	 << endl;


    for(iit2=pAcutN->ChildStates[jst_N].begin();
	iit2<pAcutN->ChildStates[jst_N].end();iit2++){
      jst_Np1=*iit2;
      jtype_Np1=pAbasisNp1->iType[jst_Np1];
      jbl_Np1=pAbasisNp1->GetBlockFromSt(jst_Np1,pos_bl);
      
//       cout << " jbl= " << jbl_Np1
// 	   << " jtype_Np1= " << jtype_Np1
// 	   << " jst_Np1= " << jst_Np1
// 	   << endl;

      
      if ( (itype_Np1==jtype_Np1)&&(ibl_Np1==jbl_Np1) ){
// 	cout << " FOUND! ifound = " << ifound << endl;
// 	cout << " ibl= " << ibl_Np1
// 	      << " jbl= " << jbl_Np1
// 	      << " ist_Np1= " << ist_Np1
// 	      << " jst_Np1= " << jst_Np1
// 	      << endl;
	ChildSt_SameType.push_back( vector<int>() );
	ChildSt_SameType[ifound].push_back(ibl_Np1);
	ChildSt_SameType[ifound].push_back(ist_Np1);
	ChildSt_SameType[ifound].push_back(jst_Np1);
	ifound++;
      }
      // end if found
    }
    // end loop in children from jst
  }
  // end loop in children from ist


}
// end subroutine
////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
