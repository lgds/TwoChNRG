
#include <iostream>
#include <vector>
#include <math.h>

#include <boost/timer.hpp>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

#include "TwoChQS.hpp"



void TwoChQS_DiagHN(vector<double> Params,
		     CNRGbasisarray* pAbasis,CNRGbasisarray* pSingleSite,
		     CNRGmatrix* Qm1fNQ, CNRGarray* pAeig){


  CNRGmatrix HN(*pAbasis);
  double OldEl[2];

  double chi_N[2]={Params[0],Params[1]};
  double Lambda=Params[2];
  double auxEl;

  // Check time
  boost::timer t;

  // New approach
  HN.UpperTriangular=true;

  int icount=0;
  for (int ibl=0;ibl<pAbasis->NumBlocks();ibl++)
    {
      // HN is block diagonal (ibl,ibl) so we are setting the blocks
      HN.MatBlockMap.push_back(ibl);
      HN.MatBlockMap.push_back(ibl);

      // Each MatBlock is Nbl x Nbl in lenght
      HN.MatBlockBegEnd.push_back(icount);
      //icount+=pAbasis->GetBlockSize(ibl)*pAbasis->GetBlockSize(ibl);
      // New
      int sizebl=pAbasis->GetBlockSize(ibl);
      if (HN.UpperTriangular)
 	icount+=(sizebl*(sizebl+1))/2;
      else
	icount+=sizebl*sizebl;

      HN.MatBlockBegEnd.push_back(icount-1);

      cout << "  Setting up H_(N = " << pAbasis->Nshell << ") Block : " << ibl ;
      cout << " of " << pAbasis->NumBlocks();
      cout << "  size: " << pAbasis->GetBlockSize(ibl) << endl;
      t.restart();

      double Q1i=pAbasis->GetQNumber(ibl,0);
      double Si=pAbasis->GetQNumber(ibl,1);


      // Calculate matrix elements
      for (int ist=pAbasis->GetBlockLimit(ibl,0);
	   ist<=pAbasis->GetBlockLimit(ibl,1);ist++)
	{
	  int type_i=pAbasis->iType[ist];


// 	  for (int jst=pAbasis->GetBlockLimit(ibl,0);
//  	       jst<=pAbasis->GetBlockLimit(ibl,1);jst++)
//       NEW: calculate only half of the matrix elements
//
	  int j0;
	  // NEW THING! only calculates half of matrix els!
	  if (HN.UpperTriangular) j0=ist; 
	  else j0=pAbasis->GetBlockLimit(ibl,0);
	  for (int jst=j0;
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

		      // Will this work?
		      // if still zero, next.
		      //if (dEqual(fabs(OldEl[ich-1]),0.0)) break;
		      

		      // Get SingleSite QNumbers
		      int pos=0;
		      int iblssp=pSingleSite->GetBlockFromSt(typep,pos);
		      int iblss=pSingleSite->GetBlockFromSt(type,pos);

		      double Qtildep=pSingleSite->GetQNumber(iblssp,0);
		      double Stildep=pSingleSite->GetQNumber(iblssp,1);
		      double Sztildep=pSingleSite->GetQNumber(iblssp,2);


		      double Qtilde=pSingleSite->GetQNumber(iblss,0);
		      double Stilde=pSingleSite->GetQNumber(iblss,1);
		      double Sztilde=pSingleSite->GetQNumber(iblss,2);

		      double Soldp=Si-Sztildep;
		      double Sold=Si-Sztilde;


		      // Check Fermi Sign
		      double FermiSign=1.0;
		      // will be -1 if only one of them is zero (A XOR B)
		      // Check this one. 
		      // Fermi sign now depends on which type the thing is!!!
                      // In other words, it will be -1 if the state has only one 
		      // electron in only one of the channels
		      // Single site states +-1,+-1/2!!
                      // 
		      if ( (dEqual(Sztilde,0.5)||dEqual(Sztilde,-0.5)) ) FermiSign=-1.0;

		      double siteqnumsp[]={Qtildep,Stildep,Sztildep};
		      double siteqnums[]={Qtilde,Stilde,Sztilde};

		      //Loop in spins
		      for (int sigma=-1;sigma<=1;sigma+=2)
			{
			  double dSigma=0.5*(double)sigma;
			  double FullMatEl=0.0;
			  double auxCG[]={0.0,0.0,0.0,0.0};

			  // Loop in Szold
			  for (double Szold=Sold;Szold>=-Sold;Szold-=1.0)
			    {
			      double Szoldp=Szold-dSigma;
			      Sztilde=Si-Szold;
			      Sztildep=Sztilde+dSigma;
			      // Changing dSigma changes Sztildep
			      // Need an extra loop in the "p" block

			      // Site matrix element: finds the block
			      siteqnums[2]=Sztilde;
			      siteqnumsp[2]=Sztildep;
			      int siteblock=pSingleSite->GetBlockFromQNumbers(siteqnums);
			      int siteblockp=pSingleSite->GetBlockFromQNumbers(siteqnumsp);
			      
// 			      if ((ist==302)&&( (jst==304)||(jst==305) ))
// 				cout << " ist = " << ist
// 				     << " jst = " << jst
// 				     << " ich = " << ich
// 				     << " sigma = " << sigma
// 				     << " Szold = " << Szold
// 				     << " Sztildep = " << Sztildep
// 				     << " Sztilde = " << Sztilde
// 				     << endl;



			      if ( (dLEqual(fabs(Sztildep),Stildep))&&
				   (dLEqual(fabs(Szold),Sold))&&
				   (dLEqual(fabs(Szoldp),Soldp)) )
				{
				  // CG coefs
				  auxCG[0]=CGordan(Sold,Szold,
						   Stilde,Sztilde,
						   Si,Si);				 

				  auxCG[1]=CGordan(Soldp,Szoldp,
						   Stildep,Sztildep,
						   Si,Si);

				  auxCG[2]=CGordan(Soldp,Szoldp,
						   0.5,dSigma,
						   Sold,Szold);
				  
				  //Loop in site blocks!! Necessary!
				  for (int sitestatep=pSingleSite->GetBlockLimit(siteblockp,0);sitestatep<=pSingleSite->GetBlockLimit(siteblockp,1);sitestatep++)
				    {
				      for (int sitestate=pSingleSite->GetBlockLimit(siteblock,0);sitestate<=pSingleSite->GetBlockLimit(siteblock,1);sitestate++)
					{
					  double SpSm[2]={0.0,0.0};
					  SpSm[0]=TwoChQS_SpSm_table(sitestate,type);
					  SpSm[1]=TwoChQS_SpSm_table(sitestatep,typep);
					  auxCG[3]=TwoChQS_fd_table(ich,sigma,
								    sitestatep,sitestate);
//  					  if ((ist==302)&&( (jst==304)||(jst==305) ))
// 					    cout << " ist = " << ist
// 						 << " jst = " << jst
// 						 << " ich = " << ich
// 						 << " sigma = " << sigma
// 						 << " fd_table(" 
// 						 << sitestatep << "," 
// 						 << sitestate << ") = " 
// 						 << auxCG[3] << endl
// 						 << " SpSm1(" 
// 						 << sitestate << "," << type 
// 						 << ") = " << SpSm[0] 
// 						 << endl
// 						 << " SpSm2(" 
// 						 << sitestatep << "," << typep 
// 						 << ") = " << SpSm[1] 
// 						 << endl;
				  
					  FullMatEl+=auxCG[0]*auxCG[1]*auxCG[2]*auxCG[3]*FermiSign*SpSm[0]*SpSm[1];

					}
				      // end loop in site block
				    }
				  // end loop in site blockp
				}
			      // END Calc coefs safeguard
			    }
			  // End loop in Szold

			auxEl+=chi_N[ich-1]*OldEl[ich-1]*FullMatEl;

// 			if ((ist==302)&&( (jst==304)||(jst==305) ))
// 			  cout << " ich = " << ich
// 			       << " OldEl = " << OldEl[ich-1]
// 			       << " FullMatEl = " << FullMatEl
// 			       << " auxEl = " << auxEl
// 			       << endl;

			}
		      // end loop in sigma
		    }
		  //end loop in channels

		  HN.MatEl.push_back(auxEl);
		  // NEW: add jst,ist as well
		  //HN.MatEl.push_back(auxEl);

		}
	      // END if ist=jst
	    }
	  //END loop in jst
	}
      // END loop in ist

      cout << " ...Done. Set-up Time: " << t.elapsed() << endl;
    }
  // END Loop in blocks (ibl)


  cout << "Updating Aeig " << endl;
  pAeig->ClearAll();

  // Syncronize with HN
  *pAeig=HN;

  // Set dEn,dEigVec
  pAeig->dEn.clear();
  pAeig->dEigVec.clear();

  for (int ibl=0;ibl<HN.NumBlocks();ibl++){
      //// NEW
    //       cout << " Regularizing block " << ibl << " of " << HN.NumBlocks()-1 << endl;
    //       HN.PutInRegularForm(ibl);
    //pAbasis->PrintBlockBasis(ibl);
    //HN.PrintMatBlock(ibl,ibl);
    cout << " Diagonalizing block " << ibl << " of " << HN.NumBlocks()-1 << endl;
    HN.DiagBlock(ibl,pAeig->dEn,pAeig->dEigVec);
    //pAeig->PrintBlockEn(ibl);
  }

 
  //uset this later
  pAeig->SetE0zero();



}
// end subroutine
////////////////////////
////////////////////////
////////////////////////

// Scrap code
			  // Sum in Sztilde
// 			  for (Sztilde=-Stilde;Sztilde<=Stilde;Sztilde+=1.0)
// 			    {
// 			      double Sztildep=Sztilde+dSigma;
// 			      double Szold=Si-Sztilde;
// 			      double Szoldp=Si-Sztildep;

// 			      // CG coefs
// 			      auxCG[0]=CGordan(Sold,Szold,
// 					       Stilde,Sztilde,Si,Si);
// 			      auxCG[1]=CGordan(Soldp,Szoldp,
// 					       Stildep,Sztildep,Si,Si);
// 			      auxCG[2]=CGordan(Soldp,Szoldp,
// 					       0.5,dSigma,Sold,Szold);
// 			      // Site matrix element
// 			      siteqnums[2]=Sztilde;
// 			      siteqnumsp[2]=Sztildep;
// 			      int sitestate=pSingleSite->GetBlockFromQNumbers(siteqnums);
// 			      int sitestatep=pSingleSite->GetBlockFromQNumbers(siteqnumsp);
// 			      auxCG[3]=TwoChQS_fd_table(ich,sigma,
// 							sitestatep,sitestate);

// 			      FullMatEl+=auxCG[0]*auxCG[1]*auxCG[2]*OldEl[ich-1]*auxCG[3]*FermiSign;
			      
// 			    }
// 			  // END Loop in Sztilde
