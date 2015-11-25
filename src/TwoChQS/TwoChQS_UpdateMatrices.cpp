#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>


#include <vector>
#include <math.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "TwoChQS.hpp"

void TwoChQS_UpdateQm1fQ(CNRGbasisarray* pSingleSite,CNRGarray* pAeig, 
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

  double Qi;
  double Qj;
  double Si,Sj;

  // Loop over blocks (icount counts for each matrix)
  int icount[2]={0,0};
  for (int ibl=0;ibl<pAeig->NumBlocks();ibl++)
    {
      Qi=pAeig->GetQNumber(ibl,0);
      Si=pAeig->GetQNumber(ibl,1);
      int Nstibl=pAeig->GetBlockSize(ibl);

      for (int jbl=ibl+1;jbl<pAeig->NumBlocks();jbl++)
	{
	  Qj=pAeig->GetQNumber(jbl,0);
	  Sj=pAeig->GetQNumber(jbl,1);
	  int Nstjbl=pAeig->GetBlockSize(jbl);

	  //Loop over channels: check if the matrix elements are non-zero
	  for (int ich=1;ich<=2;ich++)
	    {
	      if (  dEqual(Qi+1.0,Qj)&&
		   ( dEqual(Si+0.5,Sj)||dEqual(Si-0.5,Sj) )  )
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
			      // New approach

			      // Get SingleSite QNumbers
			      int pos=0;
			      int iblssi=pSingleSite->GetBlockFromSt(typei,pos);
			      int iblssj=pSingleSite->GetBlockFromSt(typej,pos);


			      double Qtildei=pSingleSite->GetQNumber(iblssi,0);
			      double Stildei=pSingleSite->GetQNumber(iblssi,1);
			      double Sztildei=pSingleSite->GetQNumber(iblssi,2);
			 
			      double Qtildej=pSingleSite->GetQNumber(iblssj,0);
			      double Stildej=pSingleSite->GetQNumber(iblssj,1);
			      double Sztildej=pSingleSite->GetQNumber(iblssj,2);

			      double siteqnumsi[]={Qtildei,Stildei,Sztildei};
			      double siteqnumsj[]={Qtildej,Stildej,Sztildej};

			      double Scfi=Si-Sztildei;
			      double Scfj=Sj-Sztildej; // should be the same!!

			      double auxBasis[2]={0.0,0.0};
			      double auxCG[2]={0.0,0.0};
			      double CGnorm=1.0;

			      int sigma=0;

			      // Choose which sigma
			      auxCG[0]=CGordan(Si,Si,0.5,-0.5,Sj,Sj);
			      auxCG[1]=CGordan(Si,Si,0.5,0.5,Sj,Sj);
			      if ( dNEqual(auxCG[0],0.0) )
				{sigma=-1;CGnorm=auxCG[0];}
			      else
				{
				  if ( dNEqual(auxCG[1],0.0) )
				    {sigma=1;CGnorm=auxCG[1];}
				  else fnbasis(istbl,jstbl)=0.0;
				}

// 			      cout << " UpMat: "
// 				   << " ist = " << ist 
// 				   << " jst = " << jst 
// 				   << " CG0 = " << auxCG[0]
// 				   << " CG1 = " << auxCG[1]
// 				   << " sigma = " << sigma
// 				   << endl;

			      double dSigma=0.5*sigma;
			      // Now Sj=Si+dSigma necessarily
			      // Calculate
			      if (sigma!=0)
				{
				  for (double Szcfj=Scfj;
				       Szcfj>=-Scfj;Szcfj-=1.0)
				    {
				      Sztildei=Si-Szcfj; //Szcfi=Szcfj
				      Sztildej=Sj-Szcfj; 
				      siteqnumsi[2]=Sztildei;
				      int blockstatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);
				      siteqnumsj[2]=Sztildej;
				      int blockstatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);

				      auxCG[0]=CGordan(Scfi,Szcfj,
						       Stildei,Sztildei,
						       Si,Si);
				      auxCG[1]=CGordan(Scfj,Szcfj,
						       Stildej,Sztildej,
						       Sj,Sj);

				      for (
			  int sitestatei=pSingleSite->GetBlockLimit(blockstatei,0);
			  sitestatei<=pSingleSite->GetBlockLimit(blockstatei,1);
			  sitestatei++)
					{
					  for (
			  int sitestatej=pSingleSite->GetBlockLimit(blockstatej,0);
			  sitestatej<=pSingleSite->GetBlockLimit(blockstatej,1);
			  sitestatej++)
					    {
					      auxBasis[0]+=TwoChQS_fd_table(ich,sigma,sitestatej,sitestatei)*auxCG[0]*auxCG[1]*TwoChQS_SpSm_table(sitestatei,typei)*TwoChQS_SpSm_table(sitestatej,typej);

// 					      cout << " Szcfj = " << Szcfj
// 						   << " Stildei = " << Stildei 
// 						   << " Sztildei = " << Sztildei 
// 						   << " Stildej = " << Stildej 
// 						   << " Sztildej = " << Sztildej 
// 						   << " typei = " << typei
// 						   << " sitestatei = " << sitestatei
// 						   << " typej = " << typej
// 						   << " sitestatej = " << sitestatej
// 						   << " CG0 = " << auxCG[0]
// 						   << " CG1 = " << auxCG[1]
// 						   << endl;


					    }
					  // loop in blockstatej
					}
				      // loop in blockstatei  
				    }
				  // end Loop in Szcf

				  //Calculate fbasis
				  fnbasis(istbl,jstbl)=auxBasis[0]/CGnorm;

//  				  cout << " auxBasis = " << auxBasis[0]
// 				       << " CGnorm = " << CGnorm
//  				       << endl;


				}
			      // end if (sigma!=0)

			    }
			  // end istcf=jstcf
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

//   		  cout << "Z(i="<<ibl<<")  : " <<  Zibl << endl;
//   		  cout << "Z(j="<<jbl<<")  : " <<  Zjbl << endl;
//   		  cout << "fbasis   :" <<  fnbasis << endl;
//   		  cout << "Zi.fbasis.ZjT : " <<  fnw << endl;


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

void TwoChQS_UpdateMatrixAfterCutting(CNRGbasisarray* pSingleSite,
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

  double Qi;
  double Qj;
  double Si,Sj;

  double qnums[2];


  // Loop over blocks (icount counts for each matrix)
  int icount[2]={0,0};
  for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++)
    {
      Qi=pAeigCut->GetQNumber(ibl,0);
      Si=pAeigCut->GetQNumber(ibl,1);
      int Nstibl=pAeigCut->GetBlockSize(ibl);

      // Get the corresponding block in pAbasis!
      qnums[0]=Qi;
      qnums[1]=Si;
      int ii_corr=pAbasis->GetBlockFromQNumbers(qnums);
      int NstiblBC=pAbasis->GetBlockSize(ii_corr);

      for (int jbl=ibl+1;jbl<pAeigCut->NumBlocks();jbl++)
	{
	  Qj=pAeigCut->GetQNumber(jbl,0);
	  Sj=pAeigCut->GetQNumber(jbl,1);
	  int Nstjbl=pAeigCut->GetBlockSize(jbl);

	  // Get the corresponding block in pAbasis!
	  qnums[0]=Qj;
	  qnums[1]=Sj;
	  int jj_corr=pAbasis->GetBlockFromQNumbers(qnums);
	  int NstjblBC=pAbasis->GetBlockSize(jj_corr);

	  if ( (!dEqual(Qj,pAbasis->GetQNumber(jj_corr,0)))||
	       (!dEqual(Sj,pAbasis->GetQNumber(jj_corr,1))) )
	    {
	      cout << "Ops, problems with GetBlockFromQNumbers" << endl;
	      pAbasis->PrintQNumbers();
	      cout << qnums[0] << " " << qnums[1] << endl;
	      cout << pAbasis->GetBlockFromQNumbers(qnums) << endl;
	    }



	  //Loop over channels: check if the matrix elements are non-zero
	  if (  dEqual(Qi+1.0,Qj)&&
		( dEqual(Si+0.5,Sj)||dEqual(Si-0.5,Sj) )  )
	    {
	      for (int ich=1;ich<=2;ich++)
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
			      // New approach

			      // Get SingleSite QNumbers
			      int pos=0;
			      int iblssi=pSingleSite->GetBlockFromSt(typei,pos);
			      int iblssj=pSingleSite->GetBlockFromSt(typej,pos);


			      double Qtildei=pSingleSite->GetQNumber(iblssi,0);
			      double Stildei=pSingleSite->GetQNumber(iblssi,1);
			      double Sztildei=pSingleSite->GetQNumber(iblssi,2);
			 
			      double Qtildej=pSingleSite->GetQNumber(iblssj,0);
			      double Stildej=pSingleSite->GetQNumber(iblssj,1);
			      double Sztildej=pSingleSite->GetQNumber(iblssj,2);

			      double siteqnumsi[]={Qtildei,Stildei,Sztildei};
			      double siteqnumsj[]={Qtildej,Stildej,Sztildej};

			      double Scfi=Si-Sztildei;
			      double Scfj=Sj-Sztildej;

			      double auxBasis[2]={0.0,0.0};
			      double auxCG[2]={0.0,0.0};

			      double CGnorm=1.0;

			      int sigma=0;

			      // Choose which sigma
			      auxCG[0]=CGordan(Si,Si,0.5,-0.5,Sj,Sj);
			      auxCG[1]=CGordan(Si,Si,0.5,0.5,Sj,Sj);
			      if ( dNEqual(auxCG[0],0.0) )
				{sigma=-1;CGnorm=auxCG[0];}
			      else
				{
				  if ( dNEqual(auxCG[1],0.0) )
				    {sigma=1;CGnorm=auxCG[1];}
				  else fnbasis(istbl,jstbl)=0.0;
				}
		      
			      double dSigma=0.5*sigma;
			      // Now Sj=Si+dSigma necessarily
			      // Calculate
			      if (sigma!=0)
				{
				  for (double Szcfj=Scfj;
				       Szcfj>=-Scfj;Szcfj-=1.0)
				    {
				      Sztildei=Si-Szcfj; //Szcfi=Szcfj
				      Sztildej=Sj-Szcfj; 
				      siteqnumsi[2]=Sztildei;
				      int blockstatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);
				      siteqnumsj[2]=Sztildej;
				      int blockstatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);

				      auxCG[0]=CGordan(Scfi,Szcfj,
						       Stildei,Sztildei,
						       Si,Si);
				      auxCG[1]=CGordan(Scfj,Szcfj,
						       Stildej,Sztildej,
						       Sj,Sj);

				      for (
			  int sitestatei=pSingleSite->GetBlockLimit(blockstatei,0);
			  sitestatei<=pSingleSite->GetBlockLimit(blockstatei,1);
			  sitestatei++)
					{
					  for (
			  int sitestatej=pSingleSite->GetBlockLimit(blockstatej,0);
			  sitestatej<=pSingleSite->GetBlockLimit(blockstatej,1);
			  sitestatej++)
					    {
					      auxBasis[0]+=TwoChQS_fd_table(ich,sigma,sitestatej,sitestatei)*auxCG[0]*auxCG[1]*TwoChQS_SpSm_table(sitestatei,typei)*TwoChQS_SpSm_table(sitestatej,typej);

					    }
					  // loop in blockstatej
					}
				      // loop in blockstatei  
				    }
				  // end Loop in Szcf

				  //Calculate fbasis
				  fnbasis(istbl,jstbl)=auxBasis[0]/CGnorm;

				}
			      // end if (sigma!=0)


			    }
			  // end if icf=jcf
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
	      // end channel loop

	    }
	  //end if Q=Q'+1 etc



	}
      // end jbl loop

    }
  // end ibl loop




}
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////


// 			      for (Sztildej=-Stildej;
// 				   Sztildej<=Stildej;Sztildej+=1.0)
// 				{
// 				  siteqnumsj[2]=Sztildej;
// 				  int sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);
// 				  int i1=0;
// 				  for (double dSigma=-0.5;dSigma<=0.5;
// 				       dSigma+=1.0)
// 				    {
// 				      siteqnumsi[2]=Sztildej-dSigma;
// 				      int sitestatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);

// 				      auxCG[0]=CGordan(Scfi,Si-Sztildej+dSigma,
// 						       Stildei,Sztildej-dSigma,
// 						       Si,Si);
// 				      auxCG[1]=CGordan(Scfj,Sj-Sztildej,
// 						       Stildej,Sztildej,
// 						       Sj,Sj);
// 	  auxBasis[i1]+=TwoChQS_fd_table(ich,(int)(dSigma/0.5),sitestatej,sitestatei)*auxCG[0]*auxCG[1];

// 				      i1++;
// 				    }
// 				  // END loop in dSigma

// 				}
// 			      //END loop in Sztildej			      

			      // Check to see which of them is non-zero
// 			      auxCG[0]=CGordan(Si,Si,0.5,-0.5,Sj,Sj);
// 			      auxCG[1]=CGordan(Si,Si,0.5,0.5,Sj,Sj);
// 			      if ( !(dEqual(auxCG[0],0.0)) ) 
// 				fnbasis(istbl,jstbl)=auxBasis[0]/auxCG[0];
// 			      else
// 				{
// 				  if ( !(dEqual(auxCG[1],0.0)) ) 
// 				    fnbasis(istbl,jstbl)=auxBasis[1]/auxCG[1];
// 				  else fnbasis(istbl,jstbl)=0.0;
// 				}




			      // Loop in Sztildej,sigma:
			      //Sztildei=Sztildej-sigma
// 			      for (Sztildej=-Stildej;
// 				   Sztildej<=Stildej;Sztildej+=1.0)
// 				{
// 				  siteqnumsj[2]=Sztildej;
// 				  int sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);
// 				  int i1=0;
// 				  for (double dSigma=-0.5;dSigma<=0.5;
// 				       dSigma+=1.0)
// 				    {
// 				      siteqnumsi[2]=Sztildej-dSigma;
// 				      int sitestatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);

// 				      auxCG[0]=CGordan(Scfi,Si-Sztildej+dSigma,
// 						       Stildei,Sztildej-dSigma,
// 						       Si,Si);
// 				      auxCG[1]=CGordan(Scfj,Sj-Sztildej,
// 						       Stildej,Sztildej,
// 						       Sj,Sj);
// 				      auxBasis[i1]+=TwoChQS_fd_table(ich,(int)(dSigma/0.5),sitestatej,sitestatei)*auxCG[0]*auxCG[1];

// 				      i1++;
// 				    }
// 				  // END loop in dSigma

// 				}
// 			      //END loop in Sztildej
// 			      // Check to see which of them is non-zero
// 			      auxCG[0]=CGordan(Si,Si,0.5,-0.5,Sj,Sj);
// 			      auxCG[1]=CGordan(Si,Si,0.5,0.5,Sj,Sj);
// 			      if ( !(dEqual(auxCG[0],0.0)) ) 
// 				fnbasis(istbl,jstbl)=auxBasis[0]/auxCG[0];
// 			      else
// 				{
// 				  if ( !(dEqual(auxCG[1],0.0)) ) 
// 				    fnbasis(istbl,jstbl)=auxBasis[1]/auxCG[1];
// 				  else fnbasis(istbl,jstbl)=0.0;
// 				}

// //fnbasis(istbl,jstbl)=TwoChQS_fd_table(ich,-1,typej,typei)+TwoChQS_fd_table(ich,1,typej,typei);
