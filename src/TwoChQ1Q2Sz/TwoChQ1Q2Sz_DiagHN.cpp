
#include <iostream>
#include <vector>
#include <math.h>


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

#include "TwoChQSz.hpp"



void TwoChQSz_DiagHN(vector<double> Params,
		     CNRGbasisarray* pAbasis,CNRGbasisarray* pSingleSite,
		     CNRGmatrix* Qm1fNQ, CNRGarray* pAeig){


  CNRGmatrix HN(*pAbasis);
  double OldEl[2];

  double chi_N[2]={Params[0],Params[1]};
  double Lambda=Params[2];
  double auxEl;

  int icount=0;
  for (int ibl=0;ibl<pAbasis->NumBlocks();ibl++)
    {
      // HN is block diagonal (ibl,ibl) so we are setting the blocks
      HN.MatBlockMap.push_back(ibl);
      HN.MatBlockMap.push_back(ibl);

      // Each MatBlock is Nbl x Nbl in lenght
      HN.MatBlockBegEnd.push_back(icount);
      icount+=pAbasis->GetBlockSize(ibl)*pAbasis->GetBlockSize(ibl);
      HN.MatBlockBegEnd.push_back(icount-1);

      cout << "  Setting up H_(N = " << pAbasis->Nshell << ") Block : " << ibl ;
      cout << " of " << pAbasis->NumBlocks();
      cout << "  size: " << pAbasis->GetBlockSize(ibl) << endl;

      double Q1i=pAbasis->GetQNumber(ibl,0);
      double Q2i=pAbasis->GetQNumber(ibl,1);
      double Szi=pAbasis->GetQNumber(ibl,2);


      // Calculate matrix elements
      for (int ist=pAbasis->GetBlockLimit(ibl,0);
	   ist<=pAbasis->GetBlockLimit(ibl,1);ist++)
	{
	  int type_i=pAbasis->iType[ist];

	  for (int jst=pAbasis->GetBlockLimit(ibl,0);
	       jst<=pAbasis->GetBlockLimit(ibl,1);jst++)
	    {
	      int type_j=pAbasis->iType[jst];

	      if (ist==jst)   // Diagonal terms
		{
		  auxEl=sqrt(Lambda)*(pAbasis->dEn[ist]);
		  HN.MatEl.push_back(auxEl);
		}
	      else 		  //Off-diagonal terms
		{
		  // Loop in channels
		  auxEl=0.0;
		  for (int ich=1;ich<=2;ich++)
		    {
		      OldEl[ich-1]=Qm1fNQ[ich-1].GetMatEl(pAbasis->StCameFrom[ist],pAbasis->StCameFrom[jst]);
		      int typep=type_i;
		      int type=type_j;

		      // if zero, try the h.c term
		      if (dEqual(fabs(OldEl[ich-1]),0.0))
			{
			  OldEl[ich-1]=Qm1fNQ[ich-1].GetMatEl(pAbasis->StCameFrom[jst],
						       pAbasis->StCameFrom[ist]);
			  typep=type_j;
			  type=type_i;
			}

		      // Check Fermi Sign
		      double Q1tilde=pSingleSite->GetQNumber(type,0);
		      double Q2tilde=pSingleSite->GetQNumber(type,1);
		      double FermiSign=1.0;
		      // will be -1 if only one of them is zero (A XOR B) 
		      if (dEqual(Q1tilde,0.0)!=dEqual(Q2tilde,0.0)) FermiSign=-1.0;
		      //Loop in spins
		      for (int sigma=-1;sigma<=1;sigma+=2)
			auxEl+=chi_N[ich-1]*OldEl[ich-1]*fd_table(ich,sigma,typep,type)*FermiSign;
		    }
		  //end loop in channels

		  HN.MatEl.push_back(auxEl);
		}
	      // END if ist=jst
	    }
	  //END loop in jst
	}
      // END loop in ist
    }
  // END Loop in blocks (ibl)


  cout << "Updating Aeig " << endl;
  pAeig->ClearAll();

  // Syncronize with HN
  *pAeig=HN;

  // Set dEn,dEigVec
  pAeig->dEn.clear();
  pAeig->dEigVec.clear();

  for (int ibl=0;ibl<HN.NumBlocks();ibl++)
    {
      //pAbasis->PrintBlockBasis(ibl);
      //HN.PrintMatBlock(ibl,ibl);
      HN.DiagBlock(ibl,pAeig->dEn,pAeig->dEigVec);
    }

  pAeig->SetE0zero();



}
