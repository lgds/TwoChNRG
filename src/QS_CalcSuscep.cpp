#include <iostream>
#include <vector>
#include <math.h>


#include "NRGclasses.hpp"



double QS_CalcSuscep(vector<double> Params, CNRGarray* pAeig){

  double betabar=Params[0];

  double sum=0.0;
  double sum2=0.0;

  for (int ibl=0;ibl < pAeig->NumBlocks();ibl++)
    {
      double Si=pAeig->GetQNumber(ibl,1);

      for (int ist=pAeig->GetBlockLimit(ibl,0);
	   ist<=pAeig->GetBlockLimit(ibl,1);ist++)
	{
         double auxS=(2.0*Si+1.0);
         double auxExp=exp(-betabar*pAeig->dEn[ist]);
         sum+=auxS*auxExp;
         sum2+=(1.0/12.0)*auxS*(auxS*auxS-1.0)*auxExp;
	}
    }

  return(sum2/sum);

}
