
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
 
#include "TwoChQS.hpp"
 


////////////////////////////////////////////
///////////     MatElChecks    /////////////
////////////////////////////////////////////
 
 
bool TwoChQS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2)
{
 
  double Qi=pAeigCut->GetQNumber(iblock1,0);
  double Si=pAeigCut->GetQNumber(iblock1,1);
 
  double Qj=pAeigCut->GetQNumber(iblock2,0);
  double Sj=pAeigCut->GetQNumber(iblock2,1);
 
   
  if (  (dEqual(Qj,(Qi+1.0)))&&
        ( (dEqual(Sj,(Si+0.5)))||(dEqual(Sj,(Si-0.5))) )
        )
    return(true);
  else
    return(false);
  
}

////////////////////////
////////////////////////
////////////////////////
                                                                                                      

