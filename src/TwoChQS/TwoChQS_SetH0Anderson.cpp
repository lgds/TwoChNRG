
#include <iostream>
#include <vector>
#include <math.h>


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

#include "TwoChQS.hpp"




int TwoChQS_SetH0Anderson(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0){

  double U_norm=Params[0];
  double ed_norm=Params[1];
  double chi_N[2]={Params[2],Params[3]};

  double auxEl=0.0;


  cout << "e1 = " << U_norm << endl;
  cout << "e2 = " << ed_norm+U_norm << endl;
  cout << "e3 = " << 2.0*ed_norm+3.0*U_norm << endl;

  cout << "chi_N = " << chi_N[0] << "  " << chi_N[1] << endl;

  // 1 - Diagonalize impurity Hamiltonian: H_(-1)

  CNRGarray AeigHimp(2);

  AeigHimp.Nshell=-1;
  AeigHimp.NQNumbers=2;

  // Only 3 states

  // |0>  = |-1 0> 
  AeigHimp.QNumbers.push_back(-1.0);
  AeigHimp.QNumbers.push_back(0.0);

  AeigHimp.dEn.push_back(U_norm);
  AeigHimp.BlockBegEnd.push_back(0);AeigHimp.BlockBegEnd.push_back(0);

  // |up> = |0 0.5> Simple as that (or |dn>)
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(0.5);

  AeigHimp.dEn.push_back(ed_norm+U_norm);

  AeigHimp.BlockBegEnd.push_back(1);AeigHimp.BlockBegEnd.push_back(1);


  // |up dn> = |1 0.0> 
  AeigHimp.QNumbers.push_back(1.0);
  AeigHimp.QNumbers.push_back(0.0);

  AeigHimp.dEn.push_back(2.0*ed_norm+3.0*U_norm);

  AeigHimp.BlockBegEnd.push_back(2);AeigHimp.BlockBegEnd.push_back(2);


  // Matrix Qm1fNQ[1] and Qm1fNQ[2]
  // Setting < ||f_{-1}|| > matrix elements.

  // These will be the REDUCED matrix elements


  CNRGmatrix Qm1fNQ; // Temporary

  Qm1fNQ.SyncNRGarray(AeigHimp);

  Qm1fNQ.MatEl.push_back(1.0);
  Qm1fNQ.MatBlockMap.push_back(0);
  Qm1fNQ.MatBlockMap.push_back(1);
  Qm1fNQ.MatBlockBegEnd.push_back(0);
  Qm1fNQ.MatBlockBegEnd.push_back(0);

  Qm1fNQ.MatEl.push_back(-sqrt(2.0));
  Qm1fNQ.MatBlockMap.push_back(1);
  Qm1fNQ.MatBlockMap.push_back(2);
  Qm1fNQ.MatBlockBegEnd.push_back(1);
  Qm1fNQ.MatBlockBegEnd.push_back(1);

  ////////////////////////
  // 2 - Add one site: Build basis
  ////////////////////////


  //CNRGbasisarray AbasisH0(3);

  pAbasisH0->NQNumbers=2;

  CNRGbasisarray ACut=CutStates(&AeigHimp, 100);
  QS_BuildBasis(&ACut,pAbasisH0,pSingleSite,1);



  cout << "No blocks = " << pAbasisH0->NumBlocks() << endl;

  cout << "No states = " << pAbasisH0->Nstates() << endl;


  //pAbasisH0->PrintAll();

  // 3 - Diagonalize Himp + Hcoupling.

  // Set Himp+Hc

  CNRGmatrix H0(*pAbasisH0);
  H0.UpperTriangular=false;

  int icount=0;
  for (int ibl=0;ibl<pAbasisH0->NumBlocks();ibl++)
    {
      double Qi=pAbasisH0->GetQNumber(ibl,0);
      double Si=pAbasisH0->GetQNumber(ibl,1);

       cout << "  Setting up HN Block : " << ibl << endl;
       cout << " Q = " << Qi << " S = " << Si;
       cout << "  size: " << pAbasisH0->GetBlockSize(ibl) << endl;


      // H0 is block diagonal (ibl,ibl) so we are setting the blocks
      H0.MatBlockMap.push_back(ibl);
      H0.MatBlockMap.push_back(ibl);

      // Each MatBlock is Nbl x Nbl in lenght
      H0.MatBlockBegEnd.push_back(icount);
      icount+=pAbasisH0->GetBlockSize(ibl)*pAbasisH0->GetBlockSize(ibl);
      H0.MatBlockBegEnd.push_back(icount-1);

      // Calculate matrix elements

      for (int ist=pAbasisH0->GetBlockLimit(ibl,0);
	   ist<=pAbasisH0->GetBlockLimit(ibl,1);ist++)
	{
	  int type_i=pAbasisH0->iType[ist];

	  for (int jst=pAbasisH0->GetBlockLimit(ibl,0);
	       jst<=pAbasisH0->GetBlockLimit(ibl,1);jst++)
	    {
	      int type_j=pAbasisH0->iType[jst];

	      if (ist==jst)   // Diagonal terms (no Lambda^1/2 factor here)
		{
		  auxEl=pAbasisH0->dEn[ist];
		  H0.MatEl.push_back(auxEl);
		}
	      else 		  //Off-diagonal terms
		{
		  // Reduced matrix element
		  double OldEl=Qm1fNQ.GetMatEl(pAbasisH0->StCameFrom[ist],
					   pAbasisH0->StCameFrom[jst]);
		  int typep=type_i;
		  int type=type_j;

		  // if zero, try the h.c term
		  if (dEqual(OldEl,0.0))
		    {
		      OldEl=Qm1fNQ.GetMatEl(pAbasisH0->StCameFrom[jst],
					       pAbasisH0->StCameFrom[ist]);
		      typep=type_j;
		      type=type_i;
		    }

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
		  // 

		  if ( (dEqual(Qtilde,1.0)||dEqual(Qtilde,-1.0))&&
		       (dEqual(Sztilde,0.5)||dEqual(Sztilde,-0.5)) ) FermiSign=-1.0;
		  // Fermi sign appears IF there is an 
		  // odd number of fermionic ops in site state 
		  // Only one electron "alone" in one channel, 
		  // e.g. |0 up> or |dn 0> (not |up dn>

		  // now get
		  // auxEl=sum_ch <fd_{ch 0 sigma} f_{-1 sigma}>
		  // REMEMBER:
		  // |Q S>_N+1 = sum <So Szo Stil Sztil|SS>|QoSo>|Qtil Stil>

// 		  cout << " ist = " << ist
// 		       << " jst = " << jst 
// 		       << " type_i = " << type_i 
// 		       << " type_j = " << type_j  
// 		       << " typep = " << typep  
// 		       << " type = " << type
// 		       << " OldEl = " << OldEl
// 		       << " FermiSign = " << FermiSign
// 		       << endl;


		  // Initialize
		  auxEl=0.0;

		  // Loop in channels
		  for (int ich=1;ich<=2;ich++)
		    {
		      double FullMatEl=0.0;
		      double auxCG[]={0.0,0.0,0.0,0.0};
		      double siteqnumsp[]={Qtildep,Stildep,Sztildep};
		      double siteqnums[]={Qtilde,Stilde,Sztilde};
		      // Loop in sigma
		      for (int sigma=-1;sigma<=1;sigma+=2)
			{
			  double dSigma=0.5*(double)sigma;

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
					  double SpSm[2];
					  SpSm[0]=TwoChQS_SpSm_table(sitestate,type);
					  SpSm[1]=TwoChQS_SpSm_table(sitestatep,typep);
					  auxCG[3]=TwoChQS_fd_table(ich,sigma,
								    sitestatep,sitestate);

					  // cout << " fd_table(" 
					  //   << sitestatep << "," << sitestate 
					  //   << ") = " << auxCG[3] << endl;
					  //  cout << " SpSm1(" 
					  //  << sitestate << "," << type 
					  //  << ") = " << SpSm[0] << endl;
					  //  cout << " SpSm2(" 
					  //    << sitestatep << "," << typep 
					  //   << ") = " << SpSm[1] << endl;
				  
					  FullMatEl+=auxCG[0]*auxCG[1]*auxCG[2]*OldEl*auxCG[3]*FermiSign*SpSm[0]*SpSm[1];

					}
				      // end loop in site block
				    }
				  // end loop in site blockp
				}
			      // END Calc coefs safeguard
			    }
			  // End loop in Szold

// 			  cout << "Partial for sigma = " 
// 			       << sigma 
// 			       << " : FullMat = " 
// 			       << FullMatEl << endl;

			}
		      //END Sum in sigma


		      auxEl+=chi_N[ich-1]*FullMatEl;

// 		      cout << " chi_ich = " << chi_N[ich-1] 
// 			   << " FullMatEl = " << FullMatEl << endl;

		    }
		  // END Loop in channels


		  H0.MatEl.push_back(auxEl);
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

  // Syncronize with H0
  *pAeig=H0;

  // Set dEn,dEigVec
  pAeig->dEn.clear();
  pAeig->dEigVec.clear();

  for (int ibl=0;ibl<H0.NumBlocks();ibl++)
    {
      pAbasisH0->PrintBlockBasis(ibl);
      H0.PrintMatBlock(ibl,ibl);
      H0.DiagBlock(ibl,pAeig->dEn,pAeig->dEigVec);
    }

  pAeig->SetE0zero();

  return(0);
}


////////////////////
////////////////////

void TwoChQS_SetH0CMphonon(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0,CNRGmatrix* NRGMats){


  // Parameters:

  double U_norm=Params[0];
  double ed_norm=Params[1];
  double chi_N[2]={Params[2],Params[3]};

  double w0=Params[4];
  double lambda=Params[5];
  double alpha=Params[6];

  int Nph=(int)Params[7];

  double auxEl=0.0;


  cout << "e1 = " << U_norm << endl;
  cout << "e2 = " << ed_norm+U_norm << endl;
  cout << "e3 = " << 2.0*ed_norm+3.0*U_norm << endl;

  cout << "chi_N = " << chi_N[0] << "  " << chi_N[1] << endl;


  // Set Hm1

  pAeig->ClearAll();

  CNRGbasisarray AeigHimp(2);

  AeigHimp.Nshell=-1;
  AeigHimp.NQNumbers=3; //(Q,S,mph)

  // 3xNph states

  int iblst=0;
  for (int mph=0;mph<=Nph;mph++)
    {
      // |0 m>  = |-1 0 m> 
      AeigHimp.QNumbers.push_back(-1.0);
      AeigHimp.QNumbers.push_back(0.0);
      AeigHimp.QNumbers.push_back(mph);
      
      AeigHimp.dEn.push_back(U_norm+w0*(double)mph);

      AeigHimp.BlockBegEnd.push_back(iblst);
      AeigHimp.BlockBegEnd.push_back(iblst);
      
      // |up> = |0 0.5 m> Simple as that (or |dn>)
      AeigHimp.QNumbers.push_back(0.0);
      AeigHimp.QNumbers.push_back(0.5);
      AeigHimp.QNumbers.push_back(mph);

      AeigHimp.dEn.push_back(ed_norm+U_norm+w0*(double)mph);

      AeigHimp.BlockBegEnd.push_back(iblst+1);
      AeigHimp.BlockBegEnd.push_back(iblst+1);


      // |up dn> = |1 0.0 m> 
      AeigHimp.QNumbers.push_back(1.0);
      AeigHimp.QNumbers.push_back(0.0);
      AeigHimp.QNumbers.push_back(mph);

      AeigHimp.dEn.push_back(2.0*ed_norm+3.0*U_norm+w0*(double)mph);

      AeigHimp.BlockBegEnd.push_back(iblst+2);
      AeigHimp.BlockBegEnd.push_back(iblst+2);

      iblst+=3;
    }
  // end loop in phonons

  // Matrix Qm1fNQ[1] and Qm1fNQ[2]
  // Setting < ||f_{-1}|| > matrix elements.

  // These will be the REDUCED matrix elements


  CNRGmatrix MatArray[2]; // Temporary
  // MatArray[0]: < |cd| >
  // MatArray[1]: < |(nd-1)cd| >

  MatArray[0].SyncNRGarray(AeigHimp);
  MatArray[1].SyncNRGarray(AeigHimp);

  int ibl=0;
  int iblpos=0;
  for (int mph=0;mph<=Nph;mph++)
    {
      // <QS|| cd ||Q'S'> 
      MatArray[0].MatEl.push_back(1.0);
      MatArray[0].MatBlockMap.push_back(ibl);
      MatArray[0].MatBlockMap.push_back(ibl+1);
      MatArray[0].MatBlockBegEnd.push_back(iblpos);
      MatArray[0].MatBlockBegEnd.push_back(iblpos);

      MatArray[0].MatEl.push_back(-sqrt(2.0));
      MatArray[0].MatBlockMap.push_back(ibl+1);
      MatArray[0].MatBlockMap.push_back(ibl+2);
      MatArray[0].MatBlockBegEnd.push_back(iblpos+1);
      MatArray[0].MatBlockBegEnd.push_back(iblpos+1);

      // <QS|| (nd-1)cd ||Q'S'> 

      MatArray[1].MatEl.push_back(-1.0);
      MatArray[1].MatBlockMap.push_back(ibl);
      MatArray[1].MatBlockMap.push_back(ibl+1);
      MatArray[1].MatBlockBegEnd.push_back(iblpos);
      MatArray[1].MatBlockBegEnd.push_back(iblpos);

      MatArray[1].MatEl.push_back(0.0);
      MatArray[1].MatBlockMap.push_back(ibl+1);
      MatArray[1].MatBlockMap.push_back(ibl+2);
      MatArray[1].MatBlockBegEnd.push_back(iblpos+1);
      MatArray[1].MatBlockBegEnd.push_back(iblpos+1);
      
      ibl+=3;
      iblpos+=2;

    }
  // end loop in phonons


  ////////////////////////
  // 2 - Add one site: Build basis
  ////////////////////////

  AeigHimp.PrintQNumbers();

  vector<int> CommonQNs;
  vector<int> totSpos;
  
  CommonQNs.push_back(2); // no of common qns

  CommonQNs.push_back(0); // pos Qold
  CommonQNs.push_back(1); // pos Sold
  //CommonQNs.push_back(2); // pos mph_old

  CommonQNs.push_back(0); // pos Qtilde
  CommonQNs.push_back(1); // pos Stilde
  //CommonQNs.push_back(-1); // (no match there)


  totSpos.push_back(1); // SU(2) symmetry

  BuildBasis(CommonQNs, totSpos,&AeigHimp,pAbasisH0,pSingleSite,0);

  //pAbasisH0->PrintEn();

  pAbasisH0->PrintAll();

  // Works!! Wow...

  ////////////////////////
  // 3 - Build Hamiltonian and diagonalize it.
  ////////////////////////

  CNRGmatrix H0(*pAbasisH0);

  H0.NeedOld=false;
  H0.UpperTriangular=false;
  H0.CheckForMatEl=Diag_check;
  H0.CalcHNMatEl=TwoChQS_H0ph_MatEl;


  vector<double> ParamsH0;
  ParamsH0.push_back(chi_N[0]);
  ParamsH0.push_back(chi_N[1]);
  ParamsH0.push_back(w0);
  ParamsH0.push_back(lambda);
  ParamsH0.push_back(alpha);

  //
  // Can I put this into a separate routine??
  //
  // Input: ParamsHN, pAbasis, pSingleSite, MatArray
  // Output: pAeig

//    double ma=1.0;
//    double mb=2.0;
//    double testPhonon=Calc_phMatEl(ma,mb,lambda+0.2,w0);

//    cout << " lambda = "  << lambda+0.2
//         << " w0 = " << w0
//         << endl;
//    cout << " <" << ma << "|Exp(...)|" << mb << "> = " << testPhonon <<endl; 
//  Working!

//    int ist=1;
//    int jst=0;

//    double testEl=H0.CalcHNMatEl(ParamsH0,pAbasisH0,pSingleSite,MatArray,ist,jst);

//    cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 


  H0.DiagHN(ParamsH0,pAbasisH0,pSingleSite,MatArray,pAeig);

  // Looking good! Continue...

}
// end subroutine



//   int icount=0;
//   for (int ibl=0;ibl<pAbasisH0->NumBlocks();ibl++)
//     {
//       // HN is block diagonal (ibl,ibl) so we are setting the blocks
//       H0.MatBlockMap.push_back(ibl);
//       H0.MatBlockMap.push_back(ibl);

//       // Each MatBlock is Nbl x Nbl in lenght
//       H0.MatBlockBegEnd.push_back(icount);
//       icount+=pAbasisH0->GetBlockSize(ibl)*pAbasisH0->GetBlockSize(ibl);
//       H0.MatBlockBegEnd.push_back(icount-1);

//       for (int ist=pAbasisH0->GetBlockLimit(ibl,0);
// 	   ist<=pAbasisH0->GetBlockLimit(ibl,1);ist++)
// 	{
// 	  for (int jst=pAbasisH0->GetBlockLimit(ibl,0);
// 	       jst<=pAbasisH0->GetBlockLimit(ibl,1);jst++)
// 	    {
// 	      double MatEl=H0.CalcHNMatEl(ParamsH0,
// 					  pAbasisH0, 
// 					  pSingleSite,
// 					  MatArray,
// 					  ist,jst);
// 	      H0.MatEl.push_back(MatEl);
// 	    }
// 	  // end loop in ist

// 	}
//       // end loop in ist

//     }
//   // end loop in blocks
