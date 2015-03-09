#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>


#include <vector>
#include <math.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "OneChQS.hpp"


void OneChQS_UpdateQm1fQ(CNRGmatrix* pQm1fNQ,CNRGarray* pAeig, 
			 CNRGbasisarray* pAbasis, CNRGbasisarray* pSingleSite){


  boost::numeric::ublas::matrix<double> Mtemp;
  boost::numeric::ublas::matrix<double> fnw;

  //Try this
  boost::numeric::ublas::matrix<double> Zibl;
  boost::numeric::ublas::matrix<double> Zjbl;

  boost::numeric::ublas::matrix<double> fnbasis;


  int icount=0;

  pQm1fNQ->ClearAll();

  pQm1fNQ->SyncNRGarray(*pAeig);

  // Find blocks such that Q'=Q+1, S'=S+-1/2

  for (int ii=0;ii<pAeig->NumBlocks();ii++)
    {
      double Qi=pAeig->GetQNumber(ii,0);
      double Si=pAeig->GetQNumber(ii,1);

      int Nstibl=pAeig->GetBlockSize(ii);

      for (int jj=ii+1;jj<pAeig->NumBlocks();jj++)
	{
      
	  double Qj=pAeig->GetQNumber(jj,0);
	  double Sj=pAeig->GetQNumber(jj,1);

	  int Nstjbl=pAeig->GetBlockSize(jj);

	  int kkp=0;
	  //if (Sj==Si+0.5) kkp=23;
	  //if (Sj==Si-0.5) kkp=32;
	  //if (fabs(Sj-(Si+0.5))<1.0E-10) kkp=23;
	  //if (fabs(Sj-(Si-0.5))<1.0E-10) kkp=32;

	  if (dEqual(Sj,(Si+0.5))) kkp=23;
	  if (dEqual(Sj,(Si-0.5))) kkp=32;

	  //	  if ( (Qj==Qi+1.0)&&( kkp!=0 ) )
	  //if ( ( fabs(Qj-(Qi+1.0))<1.0E-10 )&&( kkp!=0 ) )
	  if ( ( dEqual(Qj,(Qi+1.0)) )&&( kkp!=0 ) )
	    {
	      cout << " Nshell = " << pAeig->Nshell << endl;
	      cout << "  Setting up Block i : " << ii 
	         << " (size " <<  Nstibl 
	      	   << ") x Block j " << jj
	         << " (size " << Nstjbl << ")" << endl;

 	      pQm1fNQ->MatBlockMap.push_back(ii);
 	      pQm1fNQ->MatBlockMap.push_back(jj);
	      pQm1fNQ->MatBlockBegEnd.push_back(icount);
	      icount+=Nstibl*Nstjbl;
	      pQm1fNQ->MatBlockBegEnd.push_back(icount-1);



	      // Sets the uBLAS matrix for Z_QS(iblock)
	      // Sets the uBLAS matrix for Z_QS(jblock)


	      cout << "    Setting up BLAS matrices..." << endl;
// 	      boost::numeric::ublas::matrix<double> Zibl=pAeig->EigVec2BLAS(ii);
// 	      boost::numeric::ublas::matrix<double> Zjbl=pAeig->EigVec2BLAS(jj);

	      //boost::numeric::ublas::matrix<double> Zibl(Nstibl,Nstibl);
	      //boost::numeric::ublas::matrix<double> Zjbl(Nstjbl,Nstjbl);
	      Zibl.resize(Nstibl,Nstibl);
	      Zjbl.resize(Nstjbl,Nstjbl);
	      Zibl=pAeig->EigVec2BLAS(ii);
	      Zjbl=pAeig->EigVec2BLAS(jj);


	      cout << "    ...Zs done ..." << endl;
	      // Sets the uBLAS matrix <type|fn|type^'> with a bunch of zeroes

	      //boost::numeric::ublas::matrix<double> fnbasis (Nstibl,Nstjbl);
	      fnbasis.resize(Nstibl,Nstjbl);

	      // Check time
	      boost::timer t;
	      double time_elapsed;


	      // Loop on each block
	      int istbl=0;
	      for (int ist=pAbasis->GetBlockLimit(ii,0);
		   ist<=pAbasis->GetBlockLimit(ii,1);ist++)
		{
		  int typei=pAbasis->iType[ist];
		  int stcfi=pAbasis->StCameFrom[ist];
		  //cout << "Type i : " << typei << endl;
		  int jstbl=0;
		  for (int jst=pAbasis->GetBlockLimit(jj,0);
		       jst<=pAbasis->GetBlockLimit(jj,1);jst++)
		    {
		      int typej=pAbasis->iType[jst];
		      int stcfj=pAbasis->StCameFrom[jst];
		      //cout << "Type j : " << typej << endl;
		      // setting MatBlockMap
  		      ////pQm1fNQ->MatBlockMap.push_back(ii);
  		      ////pQm1fNQ->MatBlockMap.push_back(jj);

		      fnbasis(istbl,jstbl)=0.0;

		      // Both basis sates have to originate from the same 
		      // previous state! Otherwise, zero.

		      if (stcfi==stcfj)
			{
// 			  if (  
// 			      ( (typei==1) && (typej==2) &&(kkp==23) )||
// 			      ( (typei==1) && (typej==3) &&(kkp==32) )
// 			      ) 
// 			    fnbasis(istbl,jstbl)=1.0;
// 			  if ( (typej==4)&&(typei==2)&&(kkp==32) )
// 			    fnbasis(istbl,jstbl)=-sqrt((2.0*Si+1.0)/(2.0*Si));
// 			  if ( (typej==4)&&(typei==3)&&(kkp==23) )
// 			    fnbasis(istbl,jstbl)=sqrt((2.0*Si+1.0)/(2.0*Si+2.0)); 
			  // New approach
			  double Qtildei=pSingleSite->GetQNumber(typei,0);
			  double Stildei=pSingleSite->GetQNumber(typei,1);
			  double Sztildei=pSingleSite->GetQNumber(typei,2);
			 
			  double Qtildej=pSingleSite->GetQNumber(typej,0);
			  double Stildej=pSingleSite->GetQNumber(typej,1);
			  double Sztildej=pSingleSite->GetQNumber(typej,2);

			  double siteqnumsi[]={Qtildei,Stildei,Sztildei};
			  double siteqnumsj[]={Qtildej,Stildej,Sztildej};

			  double Scfi=Si-Sztildei;
			  double Scfj=Sj-Sztildej;

			  double auxBasis[2]={0.0,0.0};
			  double auxCG[2]={0.0,0.0};


			  //cout << " Si = " << Si << " Scfi = " << Scfi << " S~i = "<< Stildei << endl;
			  //cout << " Sj = " << Sj << " Scfj = " << Scfj << " S~j = "<< Stildej << endl;


			  // Loop in Sztildej,sigma: sztildei=Sztildej-sigma

			  for (Sztildej=-Stildej;
			       Sztildej<=Stildej;Sztildej+=1.0)
			    {
			      siteqnumsj[2]=Sztildej;
			      int sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);
			      //cout << " Sz~j = "<<Sztildej << " ssj =  "<< sitestatej << endl; 
			      int i1=0;
			      for (double sigma=-0.5;sigma<=0.5;sigma+=1.0)
				{
			      siteqnumsi[2]=Sztildej-sigma;
			      int sitestatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);
			      //cout << " Sz~i = "<<Sztildej-sigma << " ssi =  "<< sitestatei << " sigma = " << sigma << endl; 

			      auxCG[0]=CGordan(Scfi,Si-Sztildej+sigma,Stildei,Sztildej-sigma,Si,Si);
			      auxCG[1]=CGordan(Scfj,Sj-Sztildej,Stildej,Sztildej,Sj,Sj);
			      auxBasis[i1]+=fd_table((int)(sigma/0.5),sitestatej,sitestatei)*auxCG[0]*auxCG[1];

			      //cout << " CG 1 = " << auxCG[0] << " CG 2 = " << auxCG[1] << " mat_el = " <<  fd_table((int)(sigma/0.5),sitestatej,sitestatei) << endl;
			      i1++;
				}
			    }

			  // Check to see which of them is non-zero
			  auxCG[0]=CGordan(Si,Si,0.5,-0.5,Sj,Sj);
			  auxCG[1]=CGordan(Si,Si,0.5,0.5,Sj,Sj);
			  if ( !(dEqual(auxCG[0],0.0)) ) 
			  fnbasis(istbl,jstbl)=auxBasis[0]/auxCG[0];
			  else
			    {
			    if ( !(dEqual(auxCG[1],0.0)) ) 
			      fnbasis(istbl,jstbl)=auxBasis[1]/auxCG[1];
			    else fnbasis(istbl,jstbl)=0.0;
			    }
			    

			}
		      // END if sci=scf

		      //cout << "fn(" << istbl << "," << jstbl << ") = ";
		      //cout << fnbasis(istbl,jstbl) << endl;

		      jstbl++;
		    }
		  istbl++;

		}
      
	      cout << "    ...fnbasis done." << endl;

	      // Multiply the 3 of them. Zi has eigvecs in ROWS!

	      cout << "    Multiplying BLAS matrices... " << endl;

//  	      boost::numeric::ublas::matrix<double> Mtemp=prod(fnbasis,trans(Zjbl));
//  	      boost::numeric::ublas::matrix<double> fnw=prod (Zibl,Mtemp);
// 	      fnw=prod (Zibl,Mtemp);

	      //Mtemp.resize(Nstibl,Nstjbl);
	      fnw.resize(Nstibl,Nstjbl);
 	      //noalias(Mtemp)=prod(fnbasis,trans(Zjbl));
	      //cout << "    ...one done..." << endl;
 	      //noalias(fnw)=prod (Zibl, Mtemp );

	      t.restart();
 	      noalias(fnw)=prod (Zibl, 
			boost::numeric::ublas::matrix<double>(prod(fnbasis,trans(Zjbl))) );
	      time_elapsed=t.elapsed();
	      cout << "    ...done. Elapsed time:" << time_elapsed << endl;



   	      //cout << "Z(i="<<ii<<")  : " <<  Zibl << endl;
   	      //cout << "Z(j="<<jj<<")  : " <<  Zjbl << endl;
    	      //cout << "fbasis   :" <<  fnbasis << endl;
    	      //cout << "Zi.fbasis.ZjT : " <<  fnw << endl;


	      // Add to Qm1fNQ
	      for (int ii=0;ii<fnw.size1();ii++)
		{
		  for (int jj=0;jj<fnw.size2();jj++)
		    {
		      pQm1fNQ->MatEl.push_back(fnw(ii,jj));
		    }
		}

	    }
	}
      //Loop in jblock
    }
  // Loop in iblock


}


///////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////

void OneChQS_UpdateMatrixAfterCutting(CNRGmatrix* pQm1fNQ,CNRGmatrix* pMQQp1,
				      CNRGbasisarray* pAeigCut,
			 CNRGbasisarray* pAbasis, CNRGbasisarray* pSingleSite){


  boost::numeric::ublas::matrix<double> Mtemp;
  boost::numeric::ublas::matrix<double> fnw;

  //Try this
  boost::numeric::ublas::matrix<double> Zibl;
  boost::numeric::ublas::matrix<double> Zjbl;

  boost::numeric::ublas::matrix<double> fnbasis;

  double qnums[2];

  int icount=0;

  pQm1fNQ->ClearAll();

  pQm1fNQ->SyncNRGarray(*pAeigCut);


  // Find blocks such that Q'=Q+1, S'=S+-1/2

  for (int ii=0;ii<pAeigCut->NumBlocks();ii++)
    {

      double Qi=pAeigCut->GetQNumber(ii,0);
      double Si=pAeigCut->GetQNumber(ii,1);

      int Nstibl=pAeigCut->GetBlockSize(ii);

      // Get the corresponding block in pAbasis!
      qnums[0]=Qi;
      qnums[1]=Si;
      int ii_corr=pAbasis->GetBlockFromQNumbers(qnums);
      int NstiblBC=pAbasis->GetBlockSize(ii_corr);

      //cout << "ii = " << ii << " ii_corr = " << ii_corr << endl;
      //cout << " Nstibl = " << Nstibl << " NstiblBC = " << NstiblBC << endl;

      for (int jj=ii+1;jj<pAeigCut->NumBlocks();jj++)
	{
      
	  double Qj=pAeigCut->GetQNumber(jj,0);
	  double Sj=pAeigCut->GetQNumber(jj,1);

	  int Nstjbl=pAeigCut->GetBlockSize(jj);

	  // Get the corresponding block in pAbasis

	  qnums[0]=Qj;
	  qnums[1]=Sj;
	  int jj_corr=pAbasis->GetBlockFromQNumbers(qnums);
	  int NstjblBC=pAbasis->GetBlockSize(jj_corr);

// 	  cout << " Qj = " << Qj << " Sj = " << Sj << endl;
// 	  cout << " QjBC = " << pAbasis->GetQNumber(jj_corr,0) << " SjBC = " <<  pAbasis->GetQNumber(jj_corr,1) << endl;
	  if ( !dEqual(Qj,pAbasis->GetQNumber(jj_corr,0)) )
	    {
	      cout << "Ops, problems" << endl;
	      pAbasis->PrintQNumbers();
	      cout << qnums[0] << " " << qnums[1] << endl;
	      cout << pAbasis->GetBlockFromQNumbers(qnums) << endl;
	      exit(0);
	    }
// 	  cout << "jj = " << jj << " jj_corr = " << jj_corr << endl;
// 	  cout << " Nstjbl = " << Nstjbl << " NstjblBC = " << NstjblBC << endl;


	  int kkp=0;
	  if (dEqual(Sj,(Si+0.5))) kkp=23;
	  if (dEqual(Sj,(Si-0.5))) kkp=32;

	  if ( ( dEqual(Qj,(Qi+1.0)) )&&( kkp!=0 ) )
	    {
	      cout << " Nshell = " << pAeigCut->Nshell << endl;
	      cout << "  Setting up Block i : " << ii 
		   << " (size " <<  Nstibl << " was " << NstiblBC
	      	   << ") x Block j " << jj
		   << " (size " << Nstjbl << " was " << NstjblBC
		   << ")" << endl;

 	      pQm1fNQ->MatBlockMap.push_back(ii);
 	      pQm1fNQ->MatBlockMap.push_back(jj);
	      pQm1fNQ->MatBlockBegEnd.push_back(icount);
	      icount+=Nstibl*Nstjbl;
	      pQm1fNQ->MatBlockBegEnd.push_back(icount-1);



	      // Sets the uBLAS matrix for Z_QS(iblock)
	      // Sets the uBLAS matrix for Z_QS(jblock)


	      cout << "    Setting up BLAS matrices..." << endl;

	      Zibl.resize(Nstibl,NstiblBC);
	      Zjbl.resize(Nstjbl,NstjblBC);
	      Zibl=pAeigCut->EigVecCut2BLAS(ii);
	      Zjbl=pAeigCut->EigVecCut2BLAS(jj);


	      cout << "    ...Zs done ..." << endl;
//    	      cout << "Z(i="<<ii<<")  : " <<  Zibl << endl;
//    	      cout << "Z(j="<<jj<<")  : " <<  Zjbl << endl;


	      // Sets the uBLAS matrix <type|fn|type^'> with a bunch of zeroes

	      //boost::numeric::ublas::matrix<double> fnbasis (Nstibl,Nstjbl);
	      fnbasis.resize(NstiblBC,NstjblBC);

	      // Check time
	      boost::timer t;
	      double time_elapsed;


	      // Loop on each block
	      int istbl=0;
	      for (int ist=pAbasis->GetBlockLimit(ii_corr,0);
		   ist<=pAbasis->GetBlockLimit(ii_corr,1);ist++)
		{
		  int typei=pAbasis->iType[ist];
		  int stcfi=pAbasis->StCameFrom[ist];
		  //cout << "Type i : " << typei << endl;
		  int jstbl=0;
		  for (int jst=pAbasis->GetBlockLimit(jj_corr,0);
		       jst<=pAbasis->GetBlockLimit(jj_corr,1);jst++)
		    {
		      int typej=pAbasis->iType[jst];
		      int stcfj=pAbasis->StCameFrom[jst];

		      fnbasis(istbl,jstbl)=0.0;

		      // Both basis sates have to originate from the same 
		      // previous state! Otherwise, zero.

		      if (stcfi==stcfj)
			{
			  double Qtildei=pSingleSite->GetQNumber(typei,0);
			  double Stildei=pSingleSite->GetQNumber(typei,1);
			  double Sztildei=pSingleSite->GetQNumber(typei,2);
			 
			  double Qtildej=pSingleSite->GetQNumber(typej,0);
			  double Stildej=pSingleSite->GetQNumber(typej,1);
			  double Sztildej=pSingleSite->GetQNumber(typej,2);

			  double siteqnumsi[]={Qtildei,Stildei,Sztildei};
			  double siteqnumsj[]={Qtildej,Stildej,Sztildej};

			  double Scfi=Si-Sztildei;
			  double Scfj=Sj-Sztildej;

			  double auxBasis[2]={0.0,0.0};
			  double auxCG[2]={0.0,0.0};

			  // Loop in Sztildej,sigma: sztildei=Sztildej-sigma

			  for (Sztildej=-Stildej;
			       Sztildej<=Stildej;Sztildej+=1.0)
			    {
			      siteqnumsj[2]=Sztildej;
			      int sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);
			      int i1=0;
			      for (double sigma=-0.5;sigma<=0.5;sigma+=1.0)
				{
			      siteqnumsi[2]=Sztildej-sigma;
			      int sitestatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);

			      auxCG[0]=CGordan(Scfi,Si-Sztildej+sigma,Stildei,Sztildej-sigma,Si,Si);
			      auxCG[1]=CGordan(Scfj,Sj-Sztildej,Stildej,Sztildej,Sj,Sj);
			      auxBasis[i1]+=fd_table((int)(sigma/0.5),sitestatej,sitestatei)*auxCG[0]*auxCG[1];

			      i1++;
				}
			    }

			  // Check to see which of them is non-zero
			  auxCG[0]=CGordan(Si,Si,0.5,-0.5,Sj,Sj);
			  auxCG[1]=CGordan(Si,Si,0.5,0.5,Sj,Sj);
			  if ( !(dEqual(auxCG[0],0.0)) ) 
			  fnbasis(istbl,jstbl)=auxBasis[0]/auxCG[0];
			  else
			    {
			    if ( !(dEqual(auxCG[1],0.0)) ) 
			      fnbasis(istbl,jstbl)=auxBasis[1]/auxCG[1];
			    else fnbasis(istbl,jstbl)=0.0;
			    }
			    

			}
		      // END if sci=scf

		      //cout << "fn(" << istbl << "," << jstbl << ") = ";
		      //cout << fnbasis(istbl,jstbl) << endl;

		      jstbl++;
		    }
		  istbl++;

		}
      
	      cout << "    ...fnbasis done." << endl;
//     	      cout << "fbasis   :" <<  fnbasis << endl;

	      // Multiply the 3 of them. Zi has eigvecs in ROWS!

	      cout << "    Multiplying BLAS matrices... " << endl;

	      fnw.resize(Nstibl,Nstjbl);
 	      //noalias(Mtemp)=prod(fnbasis,trans(Zjbl));
	      //cout << "    ...one done..." << endl;
 	      //noalias(fnw)=prod (Zibl, Mtemp );

	      t.restart();
 	      noalias(fnw)=prod (Zibl, 
			boost::numeric::ublas::matrix<double>(prod(fnbasis,trans(Zjbl))) );
	      time_elapsed=t.elapsed();
	      cout << "    ...done. Elapsed time:" << time_elapsed << endl;



//    	      cout << "Z(i="<<ii<<")  : " <<  Zibl << endl;
//    	      cout << "Z(j="<<jj<<")  : " <<  Zjbl << endl;
//     	      cout << "fbasis   :" <<  fnbasis << endl;
//     	      cout << "Zi.fbasis.ZjT : " <<  fnw << endl;


	      // Add to Qm1fNQ
	      for (int ii=0;ii<fnw.size1();ii++)
		{
		  for (int jj=0;jj<fnw.size2();jj++)
		    {
		      pQm1fNQ->MatEl.push_back(fnw(ii,jj));
		    }
		}

	    }
	}
      //Loop in jblock
    }
  // Loop in iblock


}
// END 
