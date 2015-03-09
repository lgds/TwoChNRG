
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include <iostream>
#include <vector>

// Oct 07 : In construction

int QS_BuildBasis(CNRGbasisarray* pAeigCut,  CNRGbasisarray* pAbasis, 
                  CNRGbasisarray* pSingleSite, int BefCut, int KeepSz){

  //int QS_BuildBasis(CNRGbasisarray AeigCut,  CNRGbasisarray &Abasis, 
  //               CNRGbasisarray SingleSite, int Ncutoff, int BefCut){

  // Removing Ncutoff from argument list (March 2008)
  // Switch args to pointers
  int Nqns=pAeigCut->NQNumbers;
  int iq,is;
  int iqold,isold;
  double Qold,Sold,Qnew,Snew;
  double EnOld; 
  int istnew;


  //CNRGbasisarray AeigCut=CutStates(&Aeig, Ncutoff);


  pAbasis->Nshell=pAeigCut->Nshell+1;
  pAbasis->NQNumbers=pAeigCut->NQNumbers;

  pAbasis->ClearAll();

  cout << "Build Basis Nshell = " << pAbasis->Nshell 
       << " NQNums = " << pAbasis->NQNumbers << endl;


  
  // Build Block structure

  //  for (int ibl=0;ibl<Aeig.NumBlocks();ibl++)
  for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++)
    {

//       Qold=Aeig.GetQNumber(ibl,0);
//       Sold=Aeig.GetQNumber(ibl,1);
      Qold=pAeigCut->GetQNumber(ibl,0);
      Sold=pAeigCut->GetQNumber(ibl,1);
      //cout << "Qold Sold  = " << Qold << " " << Sold << endl; 

      for (int ibls=0; ibls<pSingleSite->NumBlocks(); ibls++)
	{
	  double Qsite=pSingleSite->GetQNumber(ibls,0);
	  double Ssite=pSingleSite->GetQNumber(ibls,1);
	  Snew=Sold+Ssite;
	  //cout << "Block site = " << ibls << endl; 
	  //cout << "Q,S site = " << Qsite << " " << Ssite << endl; 
	  while (Snew>=fabs(Sold-Ssite))
	    {
	      //cout << "Qnew = " << Qold+Qsite << " Snew = " << Snew << endl;
	      pAbasis->QNumbers.push_back(Qold+Qsite);
	      pAbasis->QNumbers.push_back(Snew);
	      Snew-=1.0;
	    }
	}
      
    }

  //pAbasis->PrintQNumbers();

  pAbasis->FilterQNumbers();

  //pAbasis->PrintQNumbers();
 
  // Set Type, StCameFrom, dEn, BlockBegEnd
  // Watch out!! Not necessarily one state per block in SingleSite!!

  istnew=0;
  for (int ibl=0;ibl<pAbasis->NumBlocks();ibl++)
    {
      // Block begin
      pAbasis->BlockBegEnd.push_back(istnew);

      Qnew=pAbasis->GetQNumber(ibl,0);
      Snew=pAbasis->GetQNumber(ibl,1);


//       for (int ibl_old=0;ibl_old<Aeig.NumBlocks();ibl_old++)
      for (int ibl_old=0;ibl_old<pAeigCut->NumBlocks();ibl_old++)
	{

// 	  Qold=Aeig.GetQNumber(ibl_old,0);
// 	  Sold=Aeig.GetQNumber(ibl_old,1);
	  Qold=pAeigCut->GetQNumber(ibl_old,0);
	  Sold=pAeigCut->GetQNumber(ibl_old,1);

	  for (int ibls=0; ibls<pSingleSite->NumBlocks(); ibls++)
	    {
	      double Qsite=pSingleSite->GetQNumber(ibls,0);
	      double Ssite=pSingleSite->GetQNumber(ibls,1);
	      double Szsite=pSingleSite->GetQNumber(ibls,2);
	      // end loop in SingleSite block states
	      // Extra loop (March 2008)
// 	      cout << " Qnew = " << Qnew 
// 		   << " Snew = " << Snew
// 		   << "   Qold = " << Qold 
// 		   << " Sold = " << Sold
// 		   << endl
// 		   << " Qsite = " << Qsite 
// 		   << " Ssite = " << Ssite
// 		   << " Szsite = " << Szsite
// 		   << endl;

// 	      cout << " Triangle Ineq : " << TriangIneq(Sold,Ssite,Snew) << endl;

	      for (int ists=pSingleSite->GetBlockLimit(ibls,0);
		   ists<=pSingleSite->GetBlockLimit(ibls,1);ists++)
		{

		  bool CheckQNs=(dEqual(Qnew,Qold+Qsite))&&
		       (dEqual(Snew,Sold+Szsite))&&
		       TriangIneq(Sold,Ssite,Snew);

		  if (KeepSz==1)
		    CheckQNs=(dEqual(Qnew,Qold+Qsite))&&
		      TriangIneq(Sold,Ssite,Snew);

// 		  if ( (dEqual(Qnew,Qold+Qsite))&&
// 		       (dEqual(Snew,Sold+Szsite))&&
// 		       TriangIneq(Sold,Ssite,Snew) )
		  if ( CheckQNs )
		{
		  for (int istbl=pAeigCut->BlockBegEnd[2*ibl_old];
		       istbl<=pAeigCut->BlockBegEnd[2*ibl_old+1];istbl++)
		    {
		  ////pAbasis->iType.push_back(pSingleSite->iType[ibls]);
			  // New (March 2008)
		      pAbasis->iType.push_back(pSingleSite->iType[ists]);
		      EnOld=pAeigCut->dEn[istbl];
		      pAbasis->dEn.push_back(EnOld);
// Has to be in the uncut basis if BefCut=1!
		      if (BefCut==1)
		      pAbasis->StCameFrom.push_back(pAeigCut->StCameFrom[istbl]);
		      else
		      pAbasis->StCameFrom.push_back(istbl);
		      istnew++;
		    }
		  // end loop in EigCut block states

		}
	      // end IF match

		}
	      // end loop in SingleSite block states
	    }
	  //END loop in SingleSite blocks

	}
      // END for in old structure
      // Block end
      pAbasis->BlockBegEnd.push_back(istnew-1);

    }
  //Loop in Abasis blocks
  
//   for (int ibl=0;ibl<pAbasis->NumBlocks();ibl++)
//    pAbasis->PrintBlockBasis(ibl);

  return(0);
}
