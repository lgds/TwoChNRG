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


void OneChQS_UpdateQm1fQ(CNRGmatrix &Qm1fNQ,CNRGarray Aeig, 
			 CNRGbasisarray Abasis, CNRGbasisarray SingleSite){


  boost::numeric::ublas::matrix<double> Mtemp;
  boost::numeric::ublas::matrix<double> fnw;

  int icount=0;

  Qm1fNQ.ClearAll();

  Qm1fNQ.SyncNRGarray(Aeig);

  // Find blocks such that Q'=Q+1, S'=S+-1/2

  for (int ii=0;ii<Aeig.NumBlocks();ii++)
    {
      double Qi=Aeig.GetQNumber(ii,0);
      double Si=Aeig.GetQNumber(ii,1);

      int Nstibl=Aeig.GetBlockSize(ii);

      for (int jj=ii+1;jj<Aeig.NumBlocks();jj++)
	{
      
	  double Qj=Aeig.GetQNumber(jj,0);
	  double Sj=Aeig.GetQNumber(jj,1);

	  int Nstjbl=Aeig.GetBlockSize(jj);

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

	      cout << "  Setting up Block i : " << ii 
	         << " (size " <<  Nstibl 
	      	   << ") x Block j " << jj
	         << " (size " << Nstjbl << ")" << endl;

 	      Qm1fNQ.MatBlockMap.push_back(ii);
 	      Qm1fNQ.MatBlockMap.push_back(jj);
	      Qm1fNQ.MatBlockBegEnd.push_back(icount);
	      icount+=Nstibl*Nstjbl;
	      Qm1fNQ.MatBlockBegEnd.push_back(icount-1);



	      // Sets the uBLAS matrix for Z_QS(iblock)
	      // Sets the uBLAS matrix for Z_QS(jblock)


	      cout << "    Setting up BLAS matrices..." << endl;
// 	      boost::numeric::ublas::matrix<double> Zibl=Aeig.EigVec2BLAS(ii);
// 	      boost::numeric::ublas::matrix<double> Zjbl=Aeig.EigVec2BLAS(jj);

	      boost::numeric::ublas::matrix<double> Zibl(Nstibl,Nstibl);
	      boost::numeric::ublas::matrix<double> Zjbl(Nstjbl,Nstjbl);
	      Zibl=Aeig.EigVec2BLAS(ii);
	      Zjbl=Aeig.EigVec2BLAS(jj);


	      cout << "    ...Zs done ..." << endl;
	      // Sets the uBLAS matrix <type|fn|type^'> with a bunch of zeroes

	      boost::numeric::ublas::matrix<double> fnbasis (Nstibl,Nstjbl);

	      // Check time
	      boost::timer t;
	      double time_elapsed;


	      // Loop on each block
	      int istbl=0;
	      for (int ist=Abasis.GetBlockLimit(ii,0);
		   ist<=Abasis.GetBlockLimit(ii,1);ist++)
		{
		  int typei=Abasis.iType[ist];
		  int stcfi=Abasis.StCameFrom[ist];
		  //cout << "Type i : " << typei << endl;
		  int jstbl=0;
		  for (int jst=Abasis.GetBlockLimit(jj,0);
		       jst<=Abasis.GetBlockLimit(jj,1);jst++)
		    {
		      int typej=Abasis.iType[jst];
		      int stcfj=Abasis.StCameFrom[jst];
		      //cout << "Type j : " << typej << endl;
		      // setting MatBlockMap
  		      ////Qm1fNQ.MatBlockMap.push_back(ii);
  		      ////Qm1fNQ.MatBlockMap.push_back(jj);

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
			  double Qtildei=SingleSite.GetQNumber(typei,0);
			  double Stildei=SingleSite.GetQNumber(typei,1);
			  double Sztildei=SingleSite.GetQNumber(typei,2);
			 
			  double Qtildej=SingleSite.GetQNumber(typej,0);
			  double Stildej=SingleSite.GetQNumber(typej,1);
			  double Sztildej=SingleSite.GetQNumber(typej,2);

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
			      int sitestatej=SingleSite.GetBlockFromQNumbers(siteqnumsj);
			      //cout << " Sz~j = "<<Sztildej << " ssj =  "<< sitestatej << endl; 
			      int i1=0;
			      for (double sigma=-0.5;sigma<=0.5;sigma+=1.0)
				{
			      siteqnumsi[2]=Sztildej-sigma;
			      int sitestatei=SingleSite.GetBlockFromQNumbers(siteqnumsi);
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
		      Qm1fNQ.MatEl.push_back(fnw(ii,jj));
		    }
		}

	    }
	}
      //Loop in jblock
    }
  // Loop in iblock


}


