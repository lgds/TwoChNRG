#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>


#include <vector>
#include <math.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "TwoChQSz.hpp"

void TwoChQSz_UpdateQm1fQ(CNRGbasisarray* pSingleSite,CNRGarray* pAeig, 
			  CNRGbasisarray* pAbasis,CNRGmatrix* Qm1fNQ){

  // Boost matrices
  boost::numeric::ublas::matrix<double> Zibl;
  boost::numeric::ublas::matrix<double> Zjbl;
  boost::numeric::ublas::matrix<double> fnbasis;


  boost::numeric::ublas::matrix<double> fnw;
  // Check time
  boost::timer t;
  double time_elapsed;



  //Clear all
  for (int ich=1;ich<=2;ich++)
    {
      Qm1fNQ[ich-1].ClearAll();
      Qm1fNQ[ich-1].SyncNRGarray(*pAeig);
    }

  // Note: assumes  Abasis and Aeig have 
  // EXACTLY the same block structure... should check.

  double Qi[2];
  double Qj[2];
  double Szi,Szj;

  // Loop over blocks (icount counts for each matrix)
  int icount[2]={0,0};
  for (int ibl=0;ibl<pAeig->NumBlocks();ibl++)
    {
      Qi[0]=pAeig->GetQNumber(ibl,0);
      Qi[1]=pAeig->GetQNumber(ibl,1);
      Szi=pAeig->GetQNumber(ibl,2);
      int Nstibl=pAeig->GetBlockSize(ibl);

      for (int jbl=ibl+1;jbl<pAeig->NumBlocks();jbl++)
	{
	  Qj[0]=pAeig->GetQNumber(jbl,0);
	  Qj[1]=pAeig->GetQNumber(jbl,1);
	  Szj=pAeig->GetQNumber(jbl,2);
	  int Nstjbl=pAeig->GetBlockSize(jbl);

	  //Loop over channels: check if the matrix elements are non-zero
	  for (int ich=1;ich<=2;ich++)
	    {
	      if (  dEqual(Qi[ich-1]+1.0,Qj[ich-1])&&
		   ( dEqual(Szi+0.5,Szj)||dEqual(Szi-0.5,Szj) )  )
		{
		  // Set Zi,Zj
		  cout << " BC: Nshell = " << pAeig->Nshell << endl;
		  cout << "Matrix in channel: " << ich << endl;
		  cout << "  Setting up Block i : " << ibl 
		       << " (size " <<  Nstibl 
		       << ") x Block j " << jbl
		       << " (size " << Nstjbl << ") of " << pAeig->NumBlocks() << endl;
		  Qm1fNQ[ich-1].MatBlockMap.push_back(ibl);
		  Qm1fNQ[ich-1].MatBlockMap.push_back(jbl);
		  Qm1fNQ[ich-1].MatBlockBegEnd.push_back(icount[ich-1]);
		  icount[ich-1]+=Nstibl*Nstjbl;
		  Qm1fNQ[ich-1].MatBlockBegEnd.push_back(icount[ich-1]-1);

		  cout << "    Setting up BLAS matrices..." << endl;
		  //boost::numeric::ublas::matrix<double> Zibl(Nstibl,Nstibl);
		  //boost::numeric::ublas::matrix<double> Zjbl(Nstjbl,Nstjbl);
		  Zibl.resize(Nstibl,Nstibl);
		  Zjbl.resize(Nstjbl,Nstjbl);

		  Zibl=pAeig->EigVec2BLAS(ibl);
		  Zjbl=pAeig->EigVec2BLAS(jbl);
		  cout << "    ...Zs done ..." << endl;

		  //boost::numeric::ublas::matrix<double> fnbasis (Nstibl,Nstjbl);
		  fnbasis.resize(Nstibl,Nstjbl);
		  // Set-up fn basis: Loop in the basis states
		  int istbl=0;
		  for (int ist=pAbasis->GetBlockLimit(ibl,0);
		       ist<=pAbasis->GetBlockLimit(ibl,1);ist++)
		    {
		      int typei=pAbasis->iType[ist];
		      int stcfi=pAbasis->StCameFrom[ist];

		      int jstbl=0;
		      for (int jst=pAbasis->GetBlockLimit(jbl,0);
			   jst<=pAbasis->GetBlockLimit(jbl,1);jst++)
			{
			  int typej=pAbasis->iType[jst];
			  int stcfj=pAbasis->StCameFrom[jst];
			  
			  fnbasis(istbl,jstbl)=0.0;

			  if (stcfi==stcfj)
			    {
			    fnbasis(istbl,jstbl)=fd_table(ich,-1,typej,typei)+
			      fd_table(ich,1,typej,typei);
			    }
			  jstbl++;
			}
		      // end jst loop
		      istbl++;
		    }
		  // end ist loop
		  cout << "    ...fnbasis done." << endl;

		  cout << "    Multiplying BLAS matrices... " << endl;
		  
		  fnw.resize(Nstibl,Nstjbl);
		  t.restart();
		  noalias(fnw)=prod (Zibl, 
				     boost::numeric::ublas::matrix<double>(prod(fnbasis,trans(Zjbl))) );
		  time_elapsed=t.elapsed();
		  cout << "    ...done. Elapsed time:" << time_elapsed << endl;

// 		  cout << "Z(i="<<ibl<<")  : " <<  Zibl << endl;
// 		  cout << "Z(j="<<jbl<<")  : " <<  Zjbl << endl;
// 		  cout << "fbasis   :" <<  fnbasis << endl;
// 		  cout << "Zi.fbasis.ZjT : " <<  fnw << endl;


		  // Add to Qm1fNQ[ich-1]
		  for (int ii=0;ii<fnw.size1();ii++)
		    for (int jj=0;jj<fnw.size2();jj++)
		      Qm1fNQ[ich-1].MatEl.push_back(fnw(ii,jj));
			

		}
	      //end if Q=Q'+1 etc

	    }
	  // end channel loop



	}
      // end jbl loop

    }
  // end ibl loop

}
// END subroutine

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

void TwoChQSz_UpdateMatrixAfterCutting(CNRGbasisarray* pSingleSite,
				       CNRGbasisarray* pAeigCut, 
				       CNRGbasisarray* pAbasis,
				       CNRGmatrix* Qm1fNQ, 
				       CNRGmatrix* pMQQp1){

  // Boost matrices
  boost::numeric::ublas::matrix<double> Zibl;
  boost::numeric::ublas::matrix<double> Zjbl;
  boost::numeric::ublas::matrix<double> fnbasis;


  boost::numeric::ublas::matrix<double> fnw;
  // Check time
  boost::timer t;
  double time_elapsed;

  //Clear all
  for (int ich=1;ich<=2;ich++)
    {
      Qm1fNQ[ich-1].ClearAll();
      Qm1fNQ[ich-1].SyncNRGarray(*pAeigCut);
    }

  // Note: assumes  Abasis and Aeig have 
  // EXACTLY the same block structure... should check.

  double Qi[2];
  double Qj[2];
  double Szi,Szj;

  double qnums[3];


  // Loop over blocks (icount counts for each matrix)
  int icount[2]={0,0};
  for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++)
    {
      Qi[0]=pAeigCut->GetQNumber(ibl,0);
      Qi[1]=pAeigCut->GetQNumber(ibl,1);
      Szi=pAeigCut->GetQNumber(ibl,2);
      int Nstibl=pAeigCut->GetBlockSize(ibl);

      // Get the corresponding block in pAbasis!
      qnums[0]=Qi[0];
      qnums[1]=Qi[1];
      qnums[2]=Szi;
      int ii_corr=pAbasis->GetBlockFromQNumbers(qnums);
      int NstiblBC=pAbasis->GetBlockSize(ii_corr);

      for (int jbl=ibl+1;jbl<pAeigCut->NumBlocks();jbl++)
	{
	  Qj[0]=pAeigCut->GetQNumber(jbl,0);
	  Qj[1]=pAeigCut->GetQNumber(jbl,1);
	  Szj=pAeigCut->GetQNumber(jbl,2);
	  int Nstjbl=pAeigCut->GetBlockSize(jbl);

	  // Get the corresponding block in pAbasis!
	  qnums[0]=Qj[0];
	  qnums[1]=Qj[1];
	  qnums[2]=Szj;
	  int jj_corr=pAbasis->GetBlockFromQNumbers(qnums);
	  int NstjblBC=pAbasis->GetBlockSize(jj_corr);

	  if ( (!dEqual(Qj[0],pAbasis->GetQNumber(jj_corr,0)))||
	       (!dEqual(Qj[1],pAbasis->GetQNumber(jj_corr,1)))||
	       (!dEqual(Szj,pAbasis->GetQNumber(jj_corr,2))) )
	    {
	      cout << "Ops, problems with GetBlockFromQNumbers" << endl;
	      pAbasis->PrintQNumbers();
	      cout << qnums[0] << " " << qnums[1] << " " << qnums[2] << endl;
	      cout << pAbasis->GetBlockFromQNumbers(qnums) << endl;
	    }



	  //Loop over channels: check if the matrix elements are non-zero
	  for (int ich=1;ich<=2;ich++)
	    {
	      if (  dEqual(Qi[ich-1]+1.0,Qj[ich-1])&&
		   ( dEqual(Szi+0.5,Szj)||dEqual(Szi-0.5,Szj) )  )
		{
		  // Set Zi,Zj
		  cout << " AC : Nshell = " << pAeigCut->Nshell << endl;
		  cout << "Matrix in channel: " << ich << endl;
		  cout << "  Setting up Block i : " << ibl 
		       << " (size " <<  Nstibl << " was " << NstiblBC
		       << ") x Block j " << jbl
		       << " (size " << Nstjbl << " was " << NstjblBC
		       << ") of " << pAeigCut->NumBlocks() << endl;
		  Qm1fNQ[ich-1].MatBlockMap.push_back(ibl);
		  Qm1fNQ[ich-1].MatBlockMap.push_back(jbl);
		  Qm1fNQ[ich-1].MatBlockBegEnd.push_back(icount[ich-1]);
		  icount[ich-1]+=Nstibl*Nstjbl;
		  Qm1fNQ[ich-1].MatBlockBegEnd.push_back(icount[ich-1]-1);

		  cout << "    Setting up BLAS matrices..." << endl;
		  //boost::numeric::ublas::matrix<double> Zibl(Nstibl,Nstibl);
		  //boost::numeric::ublas::matrix<double> Zjbl(Nstjbl,Nstjbl);
		  Zibl.resize(Nstibl,NstiblBC);
		  Zjbl.resize(Nstjbl,NstjblBC);

		  //Zibl=pAeigCut->EigVec2BLAS(ibl);
		  //Zjbl=pAeigCut->EigVec2BLAS(jbl);
		  Zibl=pAeigCut->EigVecCut2BLAS(ibl);
		  Zjbl=pAeigCut->EigVecCut2BLAS(jbl);
		  cout << "    ...Zs done ..." << endl;

		  fnbasis.resize(NstiblBC,NstjblBC);
		  // Set-up fn basis: Loop in the basis states
		  // Watch out here... ii_corr, not ibl
		  int istbl=0;
		  for (int ist=pAbasis->GetBlockLimit(ii_corr,0);
		       ist<=pAbasis->GetBlockLimit(ii_corr,1);ist++)
		    {
		      int typei=pAbasis->iType[ist];
		      int stcfi=pAbasis->StCameFrom[ist];

		      int jstbl=0;
		      for (int jst=pAbasis->GetBlockLimit(jj_corr,0);
			   jst<=pAbasis->GetBlockLimit(jj_corr,1);jst++)
			{
			  int typej=pAbasis->iType[jst];
			  int stcfj=pAbasis->StCameFrom[jst];
			  
			  fnbasis(istbl,jstbl)=0.0;

			  if (stcfi==stcfj)
			    {
			    fnbasis(istbl,jstbl)=fd_table(ich,-1,typej,typei)+
			      fd_table(ich,1,typej,typei);
			    }
			  jstbl++;
			}
		      // end jst loop
		      istbl++;
		    }
		  // end ist loop
		  cout << "    ...fnbasis done." << endl;

		  cout << "    Multiplying BLAS matrices... " << endl;
		  
		  fnw.resize(Nstibl,Nstjbl);
		  t.restart();
		  noalias(fnw)=prod (Zibl, 
				     boost::numeric::ublas::matrix<double>(prod(fnbasis,trans(Zjbl))) );
		  time_elapsed=t.elapsed();
		  cout << "    ...done. Elapsed time:" << time_elapsed << endl;

// 		  cout << "Z(i="<<ibl<<")  : " <<  Zibl << endl;
// 		  cout << "Z(j="<<jbl<<")  : " <<  Zjbl << endl;
// 		  cout << "fbasis   :" <<  fnbasis << endl;
// 		  cout << "Zi.fbasis.ZjT : " <<  fnw << endl;


		  // Add to Qm1fNQ[ich-1]
		  for (int ii=0;ii<fnw.size1();ii++)
		    for (int jj=0;jj<fnw.size2();jj++)
		      Qm1fNQ[ich-1].MatEl.push_back(fnw(ii,jj));
			

		}
	      //end if Q=Q'+1 etc

	    }
	  // end channel loop



	}
      // end jbl loop

    }
  // end ibl loop




}
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
