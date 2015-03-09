
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

  
			      if ( (siteblock<pSingleSite->NumBlocks())&&
				   (siteblockp<pSingleSite->NumBlocks())&&
				   (dLEqual(fabs(Sztildep),Stildep))&&
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
      //pAbasisH0->PrintBlockBasis(ibl);
      //H0.PrintMatBlock(ibl,ibl);
      H0.DiagBlock(ibl,pAeig->dEn,pAeig->dEigVec);
    }

  pAeig->SetE0zero();

  return(0);
}
/////////////////////
/////////////////////

///////////////////////////////////////////////////
//////                                       //////
//////                                       //////
///// New CM phonon model (no transformation) /////
//////                                       //////
//////                                       //////
///////////////////////////////////////////////////



void TwoChQS_SetH0CMphonon(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0,CNRGmatrix* NRGMats){


  // Parameters:

  double U_norm=Params[0];
  double ed_norm=Params[1];

  // Gamma_tilde = (2 Gamma/pi)/(1/2*(1+Lambda^-1))^2
  // Params[2] = sqrt(Gamma_tilde/Lambda) = chi_S
  // Change here: chi_S = Params[2], chi_A = alpha*Params[2]
  //double chi_N[2]={Params[2],Params[3]};
  double chi_N[2]={Params[2],Params[6]*Params[2]};

  double w0=Params[4];
  double lambda=Params[5];
  double alpha=Params[6];

  int Nph=(int)Params[7];

  double auxEl=0.0;


  cout << "e1 = " << U_norm << endl;
  cout << "e2 = " << ed_norm+U_norm << endl;
  cout << "e3 = " << 2.0*ed_norm+3.0*U_norm << endl;

  cout << "chi_N(S,A) = " << chi_N[0] << "  " << chi_N[1] << endl;

  cout << " w0 = " << w0 << endl;
  cout << " lambda = " << lambda << endl;
  cout << " alpha = " << alpha << endl;
  

  // Set H0

  pAeig->ClearAll();

  CNRGarray AeigHimp(2);

  // Set basis for Hm1
  AeigHimp.Nshell=-1;
  AeigHimp.NQNumbers=2; //(Q,S)

  // |0>  = |-1 0> 
  AeigHimp.QNumbers.push_back(-1.0);
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.dEn.push_back(0.0);
  AeigHimp.BlockBegEnd.push_back(0);
  AeigHimp.BlockBegEnd.push_back(0);
  // |up> = |0 0.5> Simple as that (or |dn>)
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(0.5);
  AeigHimp.BlockBegEnd.push_back(1);
  AeigHimp.BlockBegEnd.push_back(1);
  AeigHimp.dEn.push_back(0.0);
  // |up dn> = |1 0> 
  AeigHimp.QNumbers.push_back(1.0);
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.BlockBegEnd.push_back(2);
  AeigHimp.BlockBegEnd.push_back(2);
  AeigHimp.dEn.push_back(0.0);
  //

  ////
  ////  Set the matrices
  ////

  CNRGmatrix MatArray[3];

  for (int imat=0; imat<3; imat++)
    MatArray[imat].SyncNRGarray(AeigHimp);

  //
  // MatArray[0]: Hm1 in the OLD basis
  //

  MatArray[0].MatEl.push_back(U_norm);
  MatArray[0].MatEl.push_back(ed_norm+U_norm);
  MatArray[0].MatEl.push_back(2.0*ed_norm+3.0*U_norm);

  // one state per block, diagonal
  for (int ii=0;ii<=2;ii++)
    {
      MatArray[0].MatBlockMap.push_back(ii);
      MatArray[0].MatBlockMap.push_back(ii);
      MatArray[0].MatBlockBegEnd.push_back(ii);
      MatArray[0].MatBlockBegEnd.push_back(ii);
    }

//   cout << " Hm1 : " << endl;
//   for (int ibl=0;ibl<MatArray[0].NumMatBlocks();ibl++)
//     MatArray[0].PrintMatBlock(ibl);

  //
  // MatArray[1]: (nd_sigma-1) in the OLD basis
  //

  MatArray[1].MatEl.push_back(-1.0);
  MatArray[1].MatEl.push_back(0.0);
  MatArray[1].MatEl.push_back(1.0);

  // one state per block, diagonal
  for (int ii=0;ii<=2;ii++)
    {
      MatArray[1].MatBlockMap.push_back(ii);
      MatArray[1].MatBlockMap.push_back(ii);
      MatArray[1].MatBlockBegEnd.push_back(ii);
      MatArray[1].MatBlockBegEnd.push_back(ii);
    }

  // MatArray[2]: reduced <||d||> elements in the OLD basis

  MatArray[2].MatEl.push_back(1.0);
  MatArray[2].MatBlockMap.push_back(0);
  MatArray[2].MatBlockMap.push_back(1);
  MatArray[2].MatBlockBegEnd.push_back(0);
  MatArray[2].MatBlockBegEnd.push_back(0);

  MatArray[2].MatEl.push_back(-sqrt(2.0));
  MatArray[2].MatBlockMap.push_back(1);
  MatArray[2].MatBlockMap.push_back(2);
  MatArray[2].MatBlockBegEnd.push_back(1);
  MatArray[2].MatBlockBegEnd.push_back(1);

  //
  // Add a site: N=0 basis 
  //

  pAbasisH0->NQNumbers=2;

  CNRGbasisarray ACut=CutStates(&AeigHimp, 100);
  QS_BuildBasis(&ACut,pAbasisH0,pSingleSite,1);


  // Add Phonons: simply copy blocks

  for (int ibl=0;ibl<pAbasisH0->NumBlocks();ibl++)
    pAbasisH0->CopyBlock(ibl,Nph);

//    pAbasisH0->PrintAll();

  // Add vectors to block basis
  // 3xNph states


  CNRGmatrix H0(*pAbasisH0);

  H0.NeedOld=false;
  H0.UpperTriangular=true;
  H0.CheckForMatEl=Diag_check;
  H0.CalcHNMatEl=TwoChQS_H0ph_MatEl;


  vector<double> ParamsH0;
  ParamsH0.push_back(chi_N[0]);
  ParamsH0.push_back(chi_N[1]);
  ParamsH0.push_back(w0);
  ParamsH0.push_back(lambda);
  ParamsH0.push_back(alpha);

  
  int ist=4;
   int jst=7;
   double testEl=H0.CalcHNMatEl(ParamsH0,pAbasisH0,pSingleSite,MatArray,ist,jst);

//    cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 

//    for (ist=8;ist<=19;ist++)
//      {
//        for (jst=8;jst<=19;jst++)
// 	 {
//     testEl=H0.CalcHNMatEl(ParamsH0,pAbasisH0,pSingleSite,MatArray,ist,jst);
//     cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 
// 	 }
//      }

   H0.DiagHN(ParamsH0,pAbasisH0,pSingleSite,MatArray,pAeig);

//    int ibl=1;
//    pAbasisH0->PrintBlockBasis(ibl);
//    H0.PrintMatBlock(ibl,ibl);

//    ibl=8;
//    pAbasisH0->PrintBlockBasis(ibl);
//    H0.PrintMatBlock(ibl,ibl);
   


}
////////////////////////////


////////////////////
////////////////////


///////////////////////////////////////////////////
//////                                       //////
//////                                       //////
/////   CM phonon model (with transformation) /////
//////                                       //////
//////                                       //////
///////////////////////////////////////////////////


void TwoChQS_SetH0CMphononwTransf(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0,CNRGmatrix* NRGMats){


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
      // <QS|| cd ||Q'S'> on each pair of blocks, sub-blocks

      //int iblp=ibl;
      // for (int mphp=mph;mphp<=Nph;mphp++) one element per block
      int iblp=0;// from zero
      for (int mphp=0;mphp<=Nph;mphp++)
	{

	  MatArray[0].MatEl.push_back(1.0);
	  MatArray[0].MatBlockMap.push_back(ibl);
	  MatArray[0].MatBlockMap.push_back(iblp+1);
	  MatArray[0].MatBlockBegEnd.push_back(iblpos);
	  MatArray[0].MatBlockBegEnd.push_back(iblpos);

	  MatArray[0].MatEl.push_back(-sqrt(2.0));
	  MatArray[0].MatBlockMap.push_back(ibl+1);
	  MatArray[0].MatBlockMap.push_back(iblp+2);
	  MatArray[0].MatBlockBegEnd.push_back(iblpos+1);
	  MatArray[0].MatBlockBegEnd.push_back(iblpos+1);

	  // <QS|| (nd-1)cd ||Q'S'> 

	  MatArray[1].MatEl.push_back(-1.0);
	  MatArray[1].MatBlockMap.push_back(ibl);
	  MatArray[1].MatBlockMap.push_back(iblp+1);
	  MatArray[1].MatBlockBegEnd.push_back(iblpos);
	  MatArray[1].MatBlockBegEnd.push_back(iblpos);

	  MatArray[1].MatEl.push_back(0.0);
	  MatArray[1].MatBlockMap.push_back(ibl+1);
	  MatArray[1].MatBlockMap.push_back(iblp+2);
	  MatArray[1].MatBlockBegEnd.push_back(iblpos+1);
	  MatArray[1].MatBlockBegEnd.push_back(iblpos+1);
      
	  iblp+=3;
	  iblpos+=2;

	}
      // end 2nd loop in phonons

      ibl+=3;

    }
  // end loop in phonons

  for (int imatbl=0;imatbl<MatArray[0].NumMatBlocks();imatbl++)
    MatArray[1].PrintMatBlock(imatbl);


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

//   pAbasisH0->PrintAll();

//   pAbasisH0->PrintBlockBasis(2);

  // Works!! Wow...

  ////////////////////////
  // 3 - Build Hamiltonian and diagonalize it.
  ////////////////////////

  CNRGmatrix H0(*pAbasisH0);

  H0.NeedOld=false;
  H0.UpperTriangular=false;
  H0.CheckForMatEl=Diag_check;
  H0.CalcHNMatEl=TwoChQS_H0phwTransf_MatEl;


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

   cout << " Gamma1 = " << chi_N[0] 
	<< " Gamma2 = " << chi_N[1] 
        << " lambda = "  << lambda
        << " w0 = " << w0
        << endl;


   double ma=0.0;
   double mb=0.0;
   double testPhonon=Calc_phMatEl(ma,mb,lambda,w0);
   cout << " <" << ma << "|Exp(...)|" << mb << "> = " << testPhonon <<endl; 
   
   ma=1.0;
   mb=0.0;
   testPhonon=Calc_phMatEl(ma,mb,lambda,w0);
   cout << " <" << ma << "|Exp(...)|" << mb << "> = " << testPhonon <<endl; 


   ma=1.0;
   mb=1.0;
   testPhonon=Calc_phMatEl(ma,mb,lambda,w0);
   cout << " <" << ma << "|Exp(...)|" << mb << "> = " << testPhonon <<endl; 


//  Working!

   int ist=3;
   int jst=4;

   double testEl=H0.CalcHNMatEl(ParamsH0,pAbasisH0,pSingleSite,MatArray,ist,jst);

   cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 

   ist=2;
   jst=5;
   testEl=H0.CalcHNMatEl(ParamsH0,pAbasisH0,pSingleSite,MatArray,ist,jst);
   cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 


   H0.DiagHN(ParamsH0,pAbasisH0,pSingleSite,MatArray,pAeig);



//    ibl=1;
//    pAbasisH0->PrintBlockBasis(ibl);
//    H0.PrintMatBlock(ibl,ibl);

//    ibl=8;
//    pAbasisH0->PrintBlockBasis(ibl);
//    H0.PrintMatBlock(ibl,ibl);
   
   //pAeig->PrintEn();

   // Looking good! Continue...

}
// end subroutine

///////////////////////////////////////
////  Two-channel Kondo model   ///////
///////////////////////////////////////

void TwoChQS_SetH0Kondo(vector<double> Params,
			CNRGbasisarray* pSingleSite,
			CNRGarray* pAeig,
			CNRGbasisarray* pAbasisH0,
			CNRGmatrix* MatArray){
  //
  // Set the 2-channel Kondo hamiltonian in the Q,S basis.
  // 

  double J1=Params[0];
  double J2=Params[1];

  double auxEl=0.0;


  cout << "J1 = " << J1 << endl;
  cout << "J2 = " << J2 << endl;

  ////////////////////////
  // 1 - Diagonalize impurity Hamiltonian: H_(-1)
  ////////////////////////

  CNRGarray AeigHimp(2);

  AeigHimp.Nshell=-1;
  AeigHimp.NQNumbers=2;

  // Only 1 state

  // |up> = |0 0.5> Simple as that (or |dn>)
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(0.5);

  AeigHimp.dEn.push_back(0.0);

  AeigHimp.BlockBegEnd.push_back(0);AeigHimp.BlockBegEnd.push_back(0);


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


  ////////////////////////
  // 3 - Diagonalize Himp + Hcoupling.
  ////////////////////////


  CNRGmatrix H0(*pAbasisH0);

  H0.NeedOld=false;
  H0.UpperTriangular=false;
  H0.CheckForMatEl=Diag_check;
  H0.CalcHNMatEl=TwoChQS_H0Kondo_MatEl;


  vector<double> ParamsH0;
  ParamsH0.push_back(J1);
  ParamsH0.push_back(J2);



   int ist=8;
   int jst=8;
   double testEl=H0.CalcHNMatEl(ParamsH0,pAbasisH0,pSingleSite,MatArray,ist,jst);

   cout << " <" << ist << "|H0|" << jst << "> = " << testEl <<endl; 


   // Actually, H0 is diagonal in this basis... what the hell, let's do it.

   H0.DiagHN(ParamsH0,pAbasisH0,pSingleSite,MatArray,pAeig);


}
// end setH0Kondo

////////
//////// For Chain calculations only
////////

void TwoChQS_SetH0Chain(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0){


  // pAbasisH0 will have the (Q S Sz=S) states in pSingleSite

  // Aeig will be a copy of Abasis.

  // Aeig will need
  // - Energies
  // - Eigenvectors

  pAbasisH0->ClearAll();
  pAeig->ClearAll();

  pAbasisH0->Nshell=0;
  pAbasisH0->NQNumbers=2;

  // |-2 0>  = |0 0> 
  pAbasisH0->QNumbers.push_back(-2.0);
  pAbasisH0->QNumbers.push_back(0.0);
  pAbasisH0->BlockBegEnd.push_back(0);pAbasisH0->BlockBegEnd.push_back(0);
  pAbasisH0->iType.push_back(0);


  // |-1 1/2>  = |0 up>, |up 0> 
  pAbasisH0->QNumbers.push_back(-1.0);
  pAbasisH0->QNumbers.push_back(0.5);
  pAbasisH0->BlockBegEnd.push_back(1);pAbasisH0->BlockBegEnd.push_back(2);
  pAbasisH0->iType.push_back(1);
  pAbasisH0->iType.push_back(2);


  // |0 0>  = |singlet>, |0  up dn>. |up dn  0> 
  pAbasisH0->QNumbers.push_back(0.0);
  pAbasisH0->QNumbers.push_back(0.0);
  pAbasisH0->BlockBegEnd.push_back(3);pAbasisH0->BlockBegEnd.push_back(5);
  pAbasisH0->iType.push_back(5);
  pAbasisH0->iType.push_back(6);
  pAbasisH0->iType.push_back(7);


  // |0 1>  = |triplet > = |up up>
  pAbasisH0->QNumbers.push_back(0.0);
  pAbasisH0->QNumbers.push_back(1.0);
  pAbasisH0->BlockBegEnd.push_back(6);pAbasisH0->BlockBegEnd.push_back(6);
  pAbasisH0->iType.push_back(10);


  // |1 1/2>  = |up  up dn>. |up dn up> 
  pAbasisH0->QNumbers.push_back(1.0);
  pAbasisH0->QNumbers.push_back(0.5);
  pAbasisH0->BlockBegEnd.push_back(7);pAbasisH0->BlockBegEnd.push_back(8);
  pAbasisH0->iType.push_back(11);
  pAbasisH0->iType.push_back(12);


  // |2 0>  = |up dn  up dn>
  pAbasisH0->QNumbers.push_back(2.0);
  pAbasisH0->QNumbers.push_back(0.0);
  pAbasisH0->BlockBegEnd.push_back(9);pAbasisH0->BlockBegEnd.push_back(9);
  pAbasisH0->iType.push_back(15);


  // Aeig is the same

  *pAeig=*pAbasisH0;
  
  for (int ibl=0;ibl<pAeig->NumBlocks();ibl++)
    {
    for (int ii=pAeig->GetBlockLimit(ibl,0);
	 ii<=pAeig->GetBlockLimit(ibl,1);ii++)
      {
	pAeig->dEn.push_back(0); // all the same
      pAbasisH0->StCameFrom.push_back(0);
      for (int jj=pAeig->GetBlockLimit(ibl,0);
	   jj<=pAeig->GetBlockLimit(ibl,1);jj++)
	{
	  if (ii==jj)
	    pAeig->dEigVec.push_back(1.0);
	  else
	    pAeig->dEigVec.push_back(0.0);
	}
      }
    //loop in ii
    }
  // Loop in blocks
  // Type labels the state


}


