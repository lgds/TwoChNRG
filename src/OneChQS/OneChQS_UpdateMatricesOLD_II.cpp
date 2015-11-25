#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>


#include <vector>
#include <math.h>

#include "NRGclasses.hpp"


void OneChQS_UpdateQm1fQ(CNRGmatrix &Qm1fNQ,CNRGarray Aeig, 
			 CNRGbasisarray Abasis){


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
	  if (fabs(Sj-(Si+0.5))<1.0E-10) kkp=23;
	  if (fabs(Sj-(Si-0.5))<1.0E-10) kkp=32;

	  //	  if ( (Qj==Qi+1.0)&&( kkp!=0 ) )
	  if ( ( fabs(Qj-(Qi+1.0))<1.0E-10 )&&( kkp!=0 ) )
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
//  		      Qm1fNQ.MatBlockMap.push_back(ii);
//  		      Qm1fNQ.MatBlockMap.push_back(jj);

		      fnbasis(istbl,jstbl)=0.0;

		      // Both basis sates have to originate from the same 
		      // previous state! Otherwise, zero.

		      if (stcfi==stcfj)
			{
			  if (  
			      ( (typei==1) && (typej==2) &&(kkp==23) )||
			      ( (typei==1) && (typej==3) &&(kkp==32) )
			      ) 
			    fnbasis(istbl,jstbl)=1.0;
			  if ( (typej==4)&&(typei==2)&&(kkp==32) )
			    fnbasis(istbl,jstbl)=-sqrt((2.0*Si+1.0)/(2.0*Si));
			  if ( (typej==4)&&(typei==3)&&(kkp==23) )
			    fnbasis(istbl,jstbl)=sqrt((2.0*Si+1.0)/(2.0*Si+2.0)); 
			}
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



//   	      cout << "Z(i="<<ii<<")  : " <<  Zibl << endl;
//   	      cout << "Z(j="<<jj<<")  : " <<  Zjbl << endl;
//    	      cout << "fbasis   :" <<  fnbasis << endl;
//    	      cout << "Zi.fbasis.ZjT : " <<  fnw << endl;


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




//    cout << "Qm1fNQ MatBlockMap : ";
//    for (int ii=0;ii<Qm1fNQ.MatBlockMap.size();ii++)
//      cout << Qm1fNQ.MatBlockMap[ii] << " ";
//    cout << endl;

// Too slow.
//     Qm1fNQ.FilterMap_SetBegEnd();


//    cout << "Filtered Qm1fNQ MatBlockMap : ";
//    for (int ii=0;ii<Qm1fNQ.MatBlockMap.size();ii++)
//      cout << Qm1fNQ.MatBlockMap[ii] << " ";
//    cout << endl;
//    cout << "Qm1fNQ MatBlockBegEnd : ";
//    for (int ii=0;ii<Qm1fNQ.MatBlockBegEnd.size();ii++)
//      cout << Qm1fNQ.MatBlockBegEnd[ii] << " ";
//    cout << endl;


//   for (int ibl=0;ibl<Qm1fNQ.NumMatBlocks();ibl++)
//     {
//       Qm1fNQ.PrintMatBlock(ibl);
//     }

}


