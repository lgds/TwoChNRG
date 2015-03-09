
#include <iostream>
#include <vector>
#include <math.h>
//#include <cstring>


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"



double CalcOpAvg(vector<double> Params, 
		 CNRGbasisarray* pAeig, CNRGmatrix* pAop, 
		 bool totalS,int sqnumber){

  // Given a local diagonal operator A (values in pAop)
  // calculates
  // <A>(T)=\sum_i exp(-betabar*Ei) <i|A|i>

  //  Params[0]: betabar
  //  totalS   : true if full SU(2) symmetry used
  //  sqnumber : position of S (if totalS) 
  //             in the QNumber vector. For Q,S or Q,Sz is 1.
  //             (needed for degeneracy factor)

  double betabar=Params[0];

  double sum=0.0;
  double ZN=0.0;
  double DegFactor=1.0;

  for (int ibl=0;ibl < pAeig->NumBlocks();ibl++){
    if (totalS)
      DegFactor=2.0*(pAeig->GetQNumber(ibl,sqnumber))+1.0;
    else
      DegFactor=1.0;

//     cout << " OpAvg: bl = " << ibl << " of " << pAeig->NumBlocks()
// 	   << "  Si = " << (pAeig->GetQNumber(ibl,sqnumber))
// 	   << "  2S+1 = " << DegFactor
//  	   << endl;

    for (int ist=pAeig->GetBlockLimit(ibl,0);
	 ist<=pAeig->GetBlockLimit(ibl,1);ist++){
      double auxExp=exp(-betabar*(pAeig->dEn[ist]));
// 	double MatEl=pAop->GetMatEl(ist,ist);
      double MatEl=0.0;
      complex<double> MatElCplx=ZeroC;
      if (pAop->IsComplex){
	MatElCplx=pAop->cGetMatEl(ist,ist);
	MatEl=MatElCplx.real(); // diagonal element
      }
      else
	MatEl=pAop->GetMatEl(ist,ist);

	  

      // 	  cout << " ist = " << ist  
      // 	       << " En = " << pAeig->dEn[ist] 
      // 	       << " e^(-betabar*Ei) = " << auxExp 
      // 	       << " MatEl = "  << pAop->GetMatEl(ist,ist) 
      // 	       << endl;

      sum+=DegFactor*MatEl*auxExp;
      ZN+=DegFactor*auxExp;
    }
    // end loop in block states
  }
  // end loop in blocks

  //   cout << " OpAvg: Sum = " << sum
  //        << " ZN = " << ZN
  //        << " Sum/ZN = " << sum/ZN
  //        << endl;

  if (dNEqual(ZN,0.0)) sum=sum/ZN;


  return(sum);



}
