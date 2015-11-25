
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"


#include <iostream>
#include <vector>


int OneChQS_BuildBasis(CNRGarray Aeig,  CNRGbasisarray &Abasis, int Ncutoff){


  int Nqns=Aeig.NQNumbers;
  int iq,is;
  int iqold,isold;
  double Qold,Sold,Qnew,Snew;
  double EnOld; 
  int istnew;



  //   CNRGbasisarray AeigCut(Aeig);
  // AeigCut.SetSCFfromdEn();

  CNRGbasisarray AeigCut=CutStates(&Aeig, Ncutoff);


//    if (Aeig.Nshell==1)
//      {
//        cout << "Original : " << endl;
//        Aeig.PrintAll();
//        CNRGbasisarray AeigCut2=CutStates(&Aeig, 3);
//        cout << "Removed States below Ecutoff 2: " << endl;
//        AeigCut2.CNRGbasisarray::PrintAll();
//      }
//   Abasis.Nshell=Aeig.Nshell+1;
//   Abasis.NQNumbers=Aeig.NQNumbers;
  Abasis.Nshell=AeigCut.Nshell+1;
  Abasis.NQNumbers=AeigCut.NQNumbers;

  Abasis.ClearAll();

  cout << "Build Basis Nshell = " << Abasis.Nshell << endl;


  
  // Build Block structure

  //  for (int ibl=0;ibl<Aeig.NumBlocks();ibl++)
  for (int ibl=0;ibl<AeigCut.NumBlocks();ibl++)
    {

//       Qold=Aeig.GetQNumber(ibl,0);
//       Sold=Aeig.GetQNumber(ibl,1);
      Qold=AeigCut.GetQNumber(ibl,0);
      Sold=AeigCut.GetQNumber(ibl,1);
      
//       cout << "Block : " << ibl << ": ";
//       cout << "Qold = " << Qold  << "  "
// 	   << "Sold = " << Sold  << endl;

      Abasis.QNumbers.push_back(Qold-1.0);
      Abasis.QNumbers.push_back(Sold);
      if (Sold>=0.5) 
	{
	  Abasis.QNumbers.push_back(Qold);
	  Abasis.QNumbers.push_back(Sold-0.5);
	}
      Abasis.QNumbers.push_back(Qold);
      Abasis.QNumbers.push_back(Sold+0.5);

      Abasis.QNumbers.push_back(Qold+1.0);
      Abasis.QNumbers.push_back(Sold);
    }


  Abasis.FilterQNumbers();

 
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


	  if ( (Qold==Qnew+1.0)&&(Sold==Snew) )
	    {
// 	      for (int istbl=Aeig.BlockBegEnd[2*ibl_old];
// 		   istbl<=Aeig.BlockBegEnd[2*ibl_old+1];istbl++)
	      for (int istbl=AeigCut.BlockBegEnd[2*ibl_old];
		   istbl<=AeigCut.BlockBegEnd[2*ibl_old+1];istbl++)
		{
		  Abasis.iType.push_back(1);
// 		  EnOld=Aeig.dEn[istbl];
		  EnOld=AeigCut.dEn[istbl];
		  Abasis.dEn.push_back(EnOld);
// 		  Abasis.StCameFrom.push_back(istbl); 
// Has to be in the uncut basis!
		  Abasis.StCameFrom.push_back(AeigCut.StCameFrom[istbl]);
		  istnew++;
		}
	    }
	  if ( (Qold==Qnew)&&(Sold==Snew-0.5) )
	    {
// 	      for (int istbl=Aeig.BlockBegEnd[2*ibl_old];
// 		   istbl<=Aeig.BlockBegEnd[2*ibl_old+1];istbl++)
	      for (int istbl=AeigCut.BlockBegEnd[2*ibl_old];
		   istbl<=AeigCut.BlockBegEnd[2*ibl_old+1];istbl++)
		{
		  Abasis.iType.push_back(2);
// 		  EnOld=Aeig.dEn[istbl];
		  EnOld=AeigCut.dEn[istbl];
		  Abasis.dEn.push_back(EnOld);
// 		  Abasis.StCameFrom.push_back(istbl); 
// Has to be in the uncut basis!
		  Abasis.StCameFrom.push_back(AeigCut.StCameFrom[istbl]);
		  istnew++;
		}
	    }
	  if ( (Qold==Qnew)&&(Sold==Snew+0.5) )
	    {
// 	      for (int istbl=Aeig.BlockBegEnd[2*ibl_old];
// 		   istbl<=Aeig.BlockBegEnd[2*ibl_old+1];istbl++)
	      for (int istbl=AeigCut.BlockBegEnd[2*ibl_old];
		   istbl<=AeigCut.BlockBegEnd[2*ibl_old+1];istbl++)
		{
		  Abasis.iType.push_back(3);
// 		  EnOld=Aeig.dEn[istbl];
		  EnOld=AeigCut.dEn[istbl];
		  Abasis.dEn.push_back(EnOld);
// 		  Abasis.StCameFrom.push_back(istbl); 
// Has to be in the uncut basis!
		  Abasis.StCameFrom.push_back(AeigCut.StCameFrom[istbl]);
		  istnew++;
		}
	    }
	  if ( (Qold==Qnew-1.0)&&(Sold==Snew) )
	    {
// 	      for (int istbl=Aeig.BlockBegEnd[2*ibl_old];
// 		   istbl<=Aeig.BlockBegEnd[2*ibl_old+1];istbl++)
	      for (int istbl=AeigCut.BlockBegEnd[2*ibl_old];
		   istbl<=AeigCut.BlockBegEnd[2*ibl_old+1];istbl++)
		{
		  Abasis.iType.push_back(4);
// 		  EnOld=Aeig.dEn[istbl];
		  EnOld=AeigCut.dEn[istbl];
		  Abasis.dEn.push_back(EnOld);
// 		  Abasis.StCameFrom.push_back(istbl); 
// Has to be in the uncut basis!
		  Abasis.StCameFrom.push_back(AeigCut.StCameFrom[istbl]);
		  istnew++;
		}
	    }
	}
      // END for in old structure
      // Block end
      Abasis.BlockBegEnd.push_back(istnew-1);

    }
  //Loop in Abasis blocks
  
//   for (int ibl=0;ibl<Abasis.NumBlocks();ibl++)
//     Abasis.PrintBlockBasis(ibl);



  return(0);
}

