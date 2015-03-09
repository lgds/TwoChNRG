

#include "NRGclasses.hpp"


void CNRGmatrix::DiagHN(vector<double> ParamsHN,
			CNRGbasisarray* pAbasis, 
			CNRGbasisarray* pSingleSite,
			CNRGmatrix* MatArray,
			CNRGarray* pAeig,bool display)

{

  // New (July 08) // display is dummy, can be omitted
  ClearAll();

  SyncNRGarray(*pAbasis);

  int icount=0;
  for (int ibl=0;ibl<pAbasis->NumBlocks();ibl++)
    //for (int ibl=0;ibl<=2;ibl++)
    {
      cout << " Setting up H_(N =" << pAbasis->Nshell << " : ibl = " << ibl
	   << " size = " << pAbasis->GetBlockSize(ibl)
	   << " (Nbls = " << pAbasis->NumBlocks() << ")" << endl;
      cout << " Limits : " <<  pAbasis->GetBlockLimit(ibl,0) << "  " 
	   << pAbasis->GetBlockLimit(ibl,1) << endl;
      // HN is block diagonal (ibl,ibl) so we are setting the blocks
      MatBlockMap.push_back(ibl);
      MatBlockMap.push_back(ibl);

      // Each MatBlock is Nbl x Nbl in lenght
      MatBlockBegEnd.push_back(icount);
      //icount+=pAbasis->GetBlockSize(ibl)*pAbasis->GetBlockSize(ibl);
      // New
      int sizebl=pAbasis->GetBlockSize(ibl);
      if (UpperTriangular){
 	icount+=(sizebl*(sizebl+1))/2;
	//cout << "Got here " << endl;
      }
      else
	icount+=sizebl*sizebl;
      MatBlockBegEnd.push_back(icount-1);

      for (int ist=pAbasis->GetBlockLimit(ibl,0);
	   ist<=pAbasis->GetBlockLimit(ibl,1);ist++)
	{
// 	  for (int jst=pAbasis->GetBlockLimit(ibl,0);
// 	       jst<=pAbasis->GetBlockLimit(ibl,1);jst++)
	  // NEW THING! only calculates half of matrix els!
	  int j0;
	  if (UpperTriangular) j0=ist; 
	  else j0=pAbasis->GetBlockLimit(ibl,0);
	  for (int jst=j0;
	           jst<=pAbasis->GetBlockLimit(ibl,1);jst++){
	    if (!IsComplex){
	      double auxEl=CalcHNMatEl(ParamsHN,
				     pAbasis, 
				     pSingleSite,
				     MatArray,
				       ist,jst);
	      MatEl.push_back(auxEl);
	    }else{
	      complex<double> auxEl=CalcHNMatElCplx(ParamsHN,
				    pAbasis, 
				    pSingleSite,
				    MatArray,
				    ist,jst);
	      MatElCplx.push_back(auxEl);

	    }
	    // end if is complex
	    }
	  // end loop in ist

	}
      // end loop in ist

    }
  // end loop in blocks

  cout << " Done setting up HN. " << endl;
  cout << "Updating Aeig... " << endl;
  pAeig->ClearAll();

  // Syncronize with HN
  pAeig->Nshell=CNRGarray::Nshell;
  pAeig->NQNumbers=CNRGarray::NQNumbers;
  pAeig->QNumbers=CNRGarray::QNumbers;
  pAeig->BlockBegEnd=CNRGarray::BlockBegEnd;
  pAeig->totalS=CNRGarray::totalS;
  pAeig->Sqnumbers=CNRGarray::Sqnumbers;


  // Set dEn,dEigVec
  pAeig->dEn.clear();
  pAeig->dEigVec.clear();
  pAeig->cEigVec.clear();

  for (int ibl=0;ibl<NumBlocks();ibl++)
    {
      if (display){
	cout << " ******************* " << endl;
	pAbasis->PrintBlockBasis(ibl);
	PrintMatBlock(ibl,ibl);
      }
      cout << " Diagonalizing block " << ibl << " of " << NumBlocks()-1 << endl;
      if(!IsComplex) DiagBlock(ibl,pAeig->dEn,pAeig->dEigVec);
      else DiagBlock(ibl,pAeig->dEn,pAeig->cEigVec);
//       if (display){pAeig->PrintBlockEn(ibl);}
      if (display){pAeig->PrintBlock(ibl);}

    }

  //uset this later
  pAeig->SetE0zero();


}
// end subroutine
