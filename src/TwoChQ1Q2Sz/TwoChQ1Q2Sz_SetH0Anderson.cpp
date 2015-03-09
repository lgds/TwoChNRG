
#include <iostream>
#include <vector>
#include <math.h>


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

#include "TwoChQSz.hpp"




int TwoChQSz_SetH0Anderson(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0){

  double U_norm=Params[0];
  double ed_norm=Params[1];
  double chi_N[2]={Params[2],Params[3]};

  double auxEl=0.0;


  cout << "e1 = " << U_norm << endl;
  cout << "e2 = " << ed_norm+U_norm << endl;
  cout << "e3 = " << 2.0*ed_norm+3.0*U_norm << endl;

  cout << "chi_N = " << chi_N[0] << "  " << chi_N[1] << endl;

  // 1 - Diagonalize impurity Hamiltonian: H_(-1)

  CNRGarray AeigHimp(3);

  AeigHimp.Nshell=-1;
  AeigHimp.NQNumbers=3;

  // |0> (|0,0> = |-1 -1 0>) 
  AeigHimp.QNumbers.push_back(-1.0);
  AeigHimp.QNumbers.push_back(-1.0);
  AeigHimp.QNumbers.push_back(0.0);

  AeigHimp.dEn.push_back(U_norm);
  AeigHimp.BlockBegEnd.push_back(0);AeigHimp.BlockBegEnd.push_back(0);

  // |up> (|up,0> = |0 -1 0.5>) Trying something....
  // Feb 08: perhaps the best way is to write like this:
  // |up> (|up> = |0 0 0.5>)
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(-1.0);
  //AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(0.5);

  AeigHimp.dEn.push_back(ed_norm+U_norm);

  AeigHimp.BlockBegEnd.push_back(1);AeigHimp.BlockBegEnd.push_back(1);

  // |dn> (|up,0> = |0 -1 -0.5>) 
  AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(-1.0);
  //AeigHimp.QNumbers.push_back(0.0);
  AeigHimp.QNumbers.push_back(-0.5);

  AeigHimp.dEn.push_back(ed_norm+U_norm);

  AeigHimp.BlockBegEnd.push_back(2);AeigHimp.BlockBegEnd.push_back(2);

  // |dn> (|up dn,0> = |1 -1 0.0>) 
  AeigHimp.QNumbers.push_back(1.0);
  AeigHimp.QNumbers.push_back(-1.0);
  //AeigHimp.QNumbers.push_back(1.0);
  AeigHimp.QNumbers.push_back(0.0);

  AeigHimp.dEn.push_back(2.0*ed_norm+3.0*U_norm);

  AeigHimp.BlockBegEnd.push_back(3);AeigHimp.BlockBegEnd.push_back(3);


  // Matrix Qm1fNQ[1] and Qm1fNQ[2]
  // Setting < f_{1 up} + f_{1 dn} > matrix elements: careful when retrieving
  // <Sz'| f_{1 up} |Sz> = < f_{1 up} + f_{1 dn} > if Sz = Sz'+1/2 and so on...

  CNRGmatrix Qm1fNQ; // Temporary

  Qm1fNQ.SyncNRGarray(AeigHimp);

  Qm1fNQ.MatEl.push_back(1.0);
  Qm1fNQ.MatBlockMap.push_back(0);
  Qm1fNQ.MatBlockMap.push_back(1);
  Qm1fNQ.MatBlockBegEnd.push_back(0);
  Qm1fNQ.MatBlockBegEnd.push_back(0);

  Qm1fNQ.MatEl.push_back(1.0);
  Qm1fNQ.MatBlockMap.push_back(0);
  Qm1fNQ.MatBlockMap.push_back(2);
  Qm1fNQ.MatBlockBegEnd.push_back(1);
  Qm1fNQ.MatBlockBegEnd.push_back(1);


  Qm1fNQ.MatEl.push_back(-1.0);
  Qm1fNQ.MatBlockMap.push_back(1);
  Qm1fNQ.MatBlockMap.push_back(3);
  Qm1fNQ.MatBlockBegEnd.push_back(2);
  Qm1fNQ.MatBlockBegEnd.push_back(2);


  Qm1fNQ.MatEl.push_back(1.0);
  Qm1fNQ.MatBlockMap.push_back(2);
  Qm1fNQ.MatBlockMap.push_back(3);
  Qm1fNQ.MatBlockBegEnd.push_back(3);
  Qm1fNQ.MatBlockBegEnd.push_back(3);


  ////////////////////////
  // 2 - Add one site: Build basis
  ////////////////////////


  //CNRGbasisarray AbasisH0(3);

  pAbasisH0->NQNumbers=3;

  //Q1Q2Sz_BuildBasis(&AeigHimp,&AbasisH0,pSingleSite,100);
  CNRGbasisarray ACut=CutStates(&AeigHimp, 100);
  //Q1Q2Sz_BuildBasis(&ACut,pAbasisH0,pSingleSite,1); 
// This is WRONG! Abasis needs to have 4 quantum numbers

  //Let's try this

  vector<int> CommonQNs;
  
  CommonQNs.push_back(1); // tot no of common QNs
  CommonQNs.push_back(2); // pos in old
  CommonQNs.push_back(2); // pos in pSingleSite
  

  BuildBasis_RedSymmetry(CommonQNs, &ACut,pAbasisH0, pSingleSite,1);


  //cout << "No blocks = " << pAbasisH0->NumBlocks() << endl;

  //cout << "No states = " << pAbasisH0->Nstates() << endl;


  //pAbasisH0->PrintAll();

  // 3 - Diagonalize Himp + Hcoupling.

  // Set Himp+Hc

  CNRGmatrix H0(*pAbasisH0);

  int icount=0;
  for (int ibl=0;ibl<pAbasisH0->NumBlocks();ibl++)
    {
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
	  //double Q1i=pAbasisH0->GetQNumber(ist,0);
	  //double Q2i=pAbasisH0->GetQNumber(ist,1);
	  //double Szi=pAbasisH0->GetQNumber(ist,2);
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

		  // Check Fermi Sign
		  double Q1tilde=pSingleSite->GetQNumber(type,0);
		  double Q2tilde=pSingleSite->GetQNumber(type,1);
		  double FermiSign=1.0;
		  // will be -1 if only one of them is zero (A XOR B) 
		  if (dEqual(Q1tilde,0.0)!=dEqual(Q2tilde,0.0)) FermiSign=-1.0;


		  // Loop in channels (in the real one, this comes first)

		  auxEl=0.0;
		  for (int ich=1;ich<=2;ich++)
		      for (int sigma=-1;sigma<=1;sigma+=2)
		      //Loop in spins
		  auxEl+=chi_N[ich-1]*OldEl*fd_table(ich,sigma,typep,type)*FermiSign;

		  cout << " OldEl = " << OldEl << endl;
		  cout << " Off diag = " << auxEl << endl;

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

  exit(0);

  return(0);
}

