
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include <iostream>
#include <vector>

// Oct 07 : In construction

int QS_BuildBasis(CNRGarray Aeig,  CNRGbasisarray &Abasis, 
                  CNRGbasisarray SingleSite, int Ncutoff){


  int Nqns=Aeig.NQNumbers;
  int iq,is;
  int iqold,isold;
  double Qold,Sold,Qnew,Snew;
  double EnOld; 
  int istnew;


  CNRGbasisarray AeigCut=CutStates(&Aeig, Ncutoff);


  //    if (Aeig.Nshell==1)
  //    {
  //     cout << "Original : " << endl;
  //      Aeig.PrintAll();
  //      CNRGbasisarray AeigCut2=CutStates(&Aeig, 3);
  //      cout << "Removed States below Ecutoff 2: " << endl;
  //      AeigCut2.CNRGbasisarray::PrintAll();
  //    }
//   Abasis.Nshell=Aeig.Nshell+1;
//   Abasis.NQNumbers=Aeig.NQNumbers;
  Abasis.Nshell=AeigCut.Nshell+1;
  Abasis.NQNumbers=AeigCut.NQNumbers;

  Abasis.ClearAll();

  cout << "Build Basis Nshell = " << Abasis.Nshell 
       << " NQNums = " << Abasis.NQNumbers << endl;


  
  // Build Block structure

  //  for (int ibl=0;ibl<Aeig.NumBlocks();ibl++)
  for (int ibl=0;ibl<AeigCut.NumBlocks();ibl++)
    {

//       Qold=Aeig.GetQNumber(ibl,0);
//       Sold=Aeig.GetQNumber(ibl,1);
      Qold=AeigCut.GetQNumber(ibl,0);
      Sold=AeigCut.GetQNumber(ibl,1);
      //      cout << "Qold Sold  = " << Qold << " " << Sold << endl; 

      for (int ibls=0; ibls<SingleSite.NumBlocks(); ibls++)
	{
	  double Qsite=SingleSite.GetQNumber(ibls,0);
	  double Ssite=SingleSite.GetQNumber(ibls,1);
	  Snew=Sold+Ssite;
	  //  cout << "Block site = " << ibls << endl; 
	  //  cout << "Q,S site = " << Qsite << " " << Ssite << endl; 
	  while (Snew>=fabs(Sold-Ssite))
	    {
	      // cout << "Qnew = " << Qold+Qsite << " Snew = " << Snew << endl;
	      Abasis.QNumbers.push_back(Qold+Qsite);
	      Abasis.QNumbers.push_back(Snew);
	      Snew-=1.0;
	    }
	}
      
    }

  //Abasis.PrintQNumbers();

  Abasis.FilterQNumbers();

  //Abasis.PrintQNumbers();
 
  // Set Type, StCameFrom, dEn, BlockBegEnd

  istnew=0;
  for (int ibl=0;ibl<Abasis.NumBlocks();ibl++)
    {
      // Block begin
      Abasis.BlockBegEnd.push_back(istnew);

      Qnew=Abasis.GetQNumber(ibl,0);
      Snew=Abasis.GetQNumber(ibl,1);


//       for (int ibl_old=0;ibl_old<Aeig.NumBlocks();ibl_old++)
      for (int ibl_old=0;ibl_old<AeigCut.NumBlocks();ibl_old++)
	{

// 	  Qold=Aeig.GetQNumber(ibl_old,0);
// 	  Sold=Aeig.GetQNumber(ibl_old,1);
	  Qold=AeigCut.GetQNumber(ibl_old,0);
	  Sold=AeigCut.GetQNumber(ibl_old,1);

	  for (int ibls=0; ibls<SingleSite.NumBlocks(); ibls++)
	    {
	      double Qsite=SingleSite.GetQNumber(ibls,0);
	      double Ssite=SingleSite.GetQNumber(ibls,1);
	      double Szsite=SingleSite.GetQNumber(ibls,2);
	      
	      if ( (Qnew==Qold+Qsite)&&(Snew==Sold+Szsite) )
		{
		  for (int istbl=AeigCut.BlockBegEnd[2*ibl_old];
		       istbl<=AeigCut.BlockBegEnd[2*ibl_old+1];istbl++)
		    {
		      Abasis.iType.push_back(SingleSite.iType[ibls]);
		      EnOld=AeigCut.dEn[istbl];
		      Abasis.dEn.push_back(EnOld);
// Has to be in the uncut basis!
		      Abasis.StCameFrom.push_back(AeigCut.StCameFrom[istbl]);
		      istnew++;
		    }

		}
	    }
	  //END loop in SingleSite blocks

	}
      // END for in old structure
      // Block end
      Abasis.BlockBegEnd.push_back(istnew-1);

    }
  //Loop in Abasis blocks
  
  //for (int ibl=0;ibl<Abasis.NumBlocks();ibl++)
  // Abasis.PrintBlockBasis(ibl);

  return(0);
}

