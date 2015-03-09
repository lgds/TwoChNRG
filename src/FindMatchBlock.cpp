
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

int FindMatchBlock(CNRGarray* pAeigCut, int iblock, 
		    CNRGbasisarray* pAbasis){

  //
  // Given a block in pAeigCut, find the corresponding block in pAbasis 
  //

  int Nqns=pAeigCut->NQNumbers;
  
  double* qnums = new double [Nqns];

  for (int iqn=0;iqn<Nqns;iqn++)
    qnums[iqn]=pAeigCut->GetQNumber(iblock,iqn);

  int iblockBC=pAbasis->GetBlockFromQNumbers(qnums);

  delete[] qnums;

  return(iblockBC);
}
