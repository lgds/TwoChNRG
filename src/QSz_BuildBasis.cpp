
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include <iostream>
#include <vector>

// Dec 07 : In construction

int QSz_BuildBasis(CNRGbasisarray* pAeigCut,  CNRGbasisarray* pAbasis, 
                  CNRGbasisarray* pSingleSite, int BefCut){



  int Nqns=pAeigCut->NQNumbers;
  int iq,is;
  int iqold,isold;
  double Qold,Sold,Qnew,Snew;
  double EnOld; 
  int istnew;



  pAbasis->Nshell=pAeigCut->Nshell+1;
  pAbasis->NQNumbers=pAeigCut->NQNumbers;

  pAbasis->ClearAll();

  cout << "Build Basis Nshell = " << pAbasis->Nshell 
       << " NQNums = " << pAbasis->NQNumbers << endl;


  
  // Build Block structure

  for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++)
    {

      Qold=pAeigCut->GetQNumber(ibl,0);
      Sold=pAeigCut->GetQNumber(ibl,1);
      //      cout << "Qold Sold  = " << Qold << " " << Sold << endl; 

      for (int ibls=0; ibls<pSingleSite->NumBlocks(); ibls++)
	{
	  double Qsite=pSingleSite->GetQNumber(ibls,0);
	  double Ssite=pSingleSite->GetQNumber(ibls,1);
	  //  cout << "Block site = " << ibls << endl; 
	  //  cout << "Q,S site = " << Qsite << " " << Ssite << endl; 
	  pAbasis->QNumbers.push_back(Qold+Qsite);
	  pAbasis->QNumbers.push_back(Sold+Ssite);
	}
      
    }

  //pAbasis->PrintQNumbers();

  pAbasis->FilterQNumbers();

  //pAbasis->PrintQNumbers();
 
  // Set Type, StCameFrom, dEn, BlockBegEnd

  istnew=0;
  for (int ibl=0;ibl<pAbasis->NumBlocks();ibl++)
    {
      // Block begin
      pAbasis->BlockBegEnd.push_back(istnew);

      Qnew=pAbasis->GetQNumber(ibl,0);
      Snew=pAbasis->GetQNumber(ibl,1);

      for (int ibl_old=0;ibl_old<pAeigCut->NumBlocks();ibl_old++)
	{

	  Qold=pAeigCut->GetQNumber(ibl_old,0);
	  Sold=pAeigCut->GetQNumber(ibl_old,1);

	  for (int ibls=0; ibls<pSingleSite->NumBlocks(); ibls++)
	    {
	      double Qsite=pSingleSite->GetQNumber(ibls,0);
	      double Szsite=pSingleSite->GetQNumber(ibls,1);
	      
	      if ( (dEqual(Qnew,Qold+Qsite))&&(dEqual(Snew,Sold+Szsite)) )
		{
		  for (int istbl=pAeigCut->BlockBegEnd[2*ibl_old];
		       istbl<=pAeigCut->BlockBegEnd[2*ibl_old+1];istbl++)
		    {
		      pAbasis->iType.push_back(pSingleSite->iType[ibls]);
		      EnOld=pAeigCut->dEn[istbl];
		      pAbasis->dEn.push_back(EnOld);
// Has to be in the uncut basis if BefCut=1!
		      if (BefCut==1)
		      pAbasis->StCameFrom.push_back(pAeigCut->StCameFrom[istbl]);
		      else
		      pAbasis->StCameFrom.push_back(istbl);
		      istnew++;
		      //pAbasis->StCameFrom.push_back(pAeigCut->StCameFrom[istbl]);
		      //istnew++;
		    }

		}
	    }
	  //END loop in SingleSite blocks

	}
      // END for in old structure
      // Block end
      pAbasis->BlockBegEnd.push_back(istnew-1);

    }
  //Loop in Abasis blocks
  
  //for (int ibl=0;ibl<pAbasis->NumBlocks();ibl++)
  // pAbasis->PrintBlockBasis(ibl);

  return(0);
}

