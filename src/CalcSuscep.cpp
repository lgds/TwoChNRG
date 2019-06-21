#include <iostream>
#include <vector>
#include <math.h>
#include <cstring>


#include "NRGclasses.hpp"


////
//// Susceptibility
////



double CalcSuscep(vector<double> Params, CNRGarray* pAeig, 
		  int sqnumber,bool totalS){

  //  Params[0]: betabar
  //  sqnumber : position of S (if totalS) or Sz (if not totalS)  
  //             in the QNumber vector. For Q,S or Q,Sz is 1.
  //  totalS   : true if full SU(2) symmetry used

  double betabar=Params[0];

  double sum[3]={0.0,0.0,0.0};
  double SzSq=0.0;
  double Sz=0.0;
  double auxS=0.0;

  for (int ibl=0;ibl < pAeig->NumBlocks();ibl++)
    {
      double Si=pAeig->GetQNumber(ibl,sqnumber);

      for (int ist=pAeig->GetBlockLimit(ibl,0);
	   ist<=pAeig->GetBlockLimit(ibl,1);ist++)
	{

	  //// 
	  //// Get <Sz^2>, <Sz>
	  ////

	  if (totalS)
	    {
	      auxS=(2.0*Si+1.0);
	      SzSq=(1.0/12.0)*auxS*(auxS*auxS-1.0);
	      Sz=0.0;
	    }
	  else
	    {
	      auxS=1.0;
	      SzSq=Si*Si;
	      Sz=Si;
	    }
	  double auxExp=exp(-betabar*pAeig->dEn[ist]);
         sum[0]+=auxS*auxExp;
         sum[1]+=SzSq*auxExp;
	 sum[2]+=Sz*auxExp;
	}
    }

  cout << "Tr Sz^2 = " << sum[1] << endl;
  cout << "Tr Sz = " << sum[2] << endl;
  cout << "ZN = " << sum[0] << endl;
  cout << "[(Tr Sz^2)-(Tr Sz)^2]/ZN = " << (sum[1]-sum[2]*sum[2])/sum[0] << endl;

  return((sum[1]-0.0*sum[1]*sum[2])/sum[0]);

}

////
//// Entropy
////

double CalcEntropy(vector<double> Params, CNRGarray* pAeig, 
		  int sqnumber,bool totalS){

  //  Params[0]: betabar
  //  sqnumber : position of S (if totalS) or Sz (if not totalS)  
  //             in the QNumber vector. For Q,S or Q,Sz is 1.
  //  totalS   : true if full SU(2) symmetry used

  double betabar=Params[0];

  double ZN=0.0;
  double DegFactor=1.0;
  double sum=0.0;


  for (int ibl=0;ibl < pAeig->NumBlocks();ibl++)
    {
      double Si=pAeig->GetQNumber(ibl,sqnumber);

      if (totalS)
        DegFactor=2.0*Si+1.0;
      else
	DegFactor=1.0;


      for (int ist=pAeig->GetBlockLimit(ibl,0);
	   ist<=pAeig->GetBlockLimit(ibl,1);ist++)
	{

	  double En=pAeig->dEn[ist];
	  double auxExp=exp(-betabar*En);
	  sum+=En*auxExp;
	  ZN+=auxExp;
	}
    }
  
  double SN=betabar*(sum/ZN)+log(ZN);

  return(SN);

}
