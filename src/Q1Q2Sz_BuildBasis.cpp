

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include <iostream>
#include <vector>


int Q1Q2Sz_BuildBasis(CNRGbasisarray* pAeigCut, CNRGbasisarray* pAbasis, 
                  CNRGbasisarray* pSingleSite, int BefCut){


  //CNRGbasisarray AeigCut=CutStates(pAeig, Ncutoff); // This should work as is...


  pAbasis->Nshell=pAeigCut->Nshell+1;
  pAbasis->NQNumbers=pAeigCut->NQNumbers;

  pAbasis->ClearAll();


  // Build Block structure

  for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++)
    {
      
      double Q1old=pAeigCut->GetQNumber(ibl,0);
      double Q2old=pAeigCut->GetQNumber(ibl,1);
      double Szold=pAeigCut->GetQNumber(ibl,2);



      for (int ibls=0; ibls<pSingleSite->NumBlocks(); ibls++)
	{
	  double Q1site=pSingleSite->GetQNumber(ibls,0);
	  double Q2site=pSingleSite->GetQNumber(ibls,1);
	  double Szsite=pSingleSite->GetQNumber(ibls,2);

	  pAbasis->QNumbers.push_back(Q1old+Q1site);
	  pAbasis->QNumbers.push_back(Q2old+Q2site);
	  pAbasis->QNumbers.push_back(Szold+Szsite);

	}
      // END loop in SingleSite states

    }
  // END loop in Acut blocks

  // Filter blocks
  pAbasis->FilterQNumbers();

  // Build new block structure
  int istnew=0;
  for (int ibl_new=0;ibl_new<pAbasis->NumBlocks();ibl_new++)
    {

      double Q1new=pAbasis->GetQNumber(ibl_new,0);
      double Q2new=pAbasis->GetQNumber(ibl_new,1);
      double Sznew=pAbasis->GetQNumber(ibl_new,2);

      // Mark Block begin
      pAbasis->BlockBegEnd.push_back(istnew);
      for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++)
	{
      
	  double Q1old=pAeigCut->GetQNumber(ibl,0);
	  double Q2old=pAeigCut->GetQNumber(ibl,1);
	  double Szold=pAeigCut->GetQNumber(ibl,2);

	  for (int ibls=0; ibls<pSingleSite->NumBlocks(); ibls++)
	    {
	      double Q1site=pSingleSite->GetQNumber(ibls,0);
	      double Q2site=pSingleSite->GetQNumber(ibls,1);
	      double Szsite=pSingleSite->GetQNumber(ibls,2);

	      if (dEqual(Q1new,Q1old+Q1site)&&
		  dEqual(Q2new,Q2old+Q2site)&&
		  dEqual(Sznew,Szold+Szsite)  )
		{
		  // Loop in the Acut states in that block
		  for (int istbl=pAeigCut->BlockBegEnd[2*ibl];
		       istbl<=pAeigCut->BlockBegEnd[2*ibl+1];istbl++)
		    {
		      pAbasis->iType.push_back(pSingleSite->iType[ibls]);
		      pAbasis->dEn.push_back(pAeigCut->dEn[istbl]);
		      // StCameFrom in the UN-CUT basis!!
		      if (BefCut==1)
		      pAbasis->StCameFrom.push_back(pAeigCut->StCameFrom[istbl]);
		      else
		      pAbasis->StCameFrom.push_back(istbl);
		      istnew++;
		    }
		}
	      // END check if blocks match

	    }
	  // END loop in SingleSite states

	}
      // END loop in Acut blocks


      // Mark Block end
      pAbasis->BlockBegEnd.push_back(istnew-1);

    }
  // END loop in Abasis blocks



  return(0);
}
