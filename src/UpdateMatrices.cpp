#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>


#include <vector>
#include <math.h>

extern "C"{
  #include <cblas.h>
}

#define INVALID -1


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

void UpdateMatrices_uBLAS(CNRGbasisarray* pSingleSite,
			  CNRGbasisarray* pAeigCut, 
			  CNRGbasisarray* pAbasis,
			  CNRGmatrix* NRGMats, int NumNRGMats, bool display){

  // Boost matrices
  boost::numeric::ublas::matrix<double> Zibl;
  boost::numeric::ublas::matrix<double> Zjbl;
  boost::numeric::ublas::matrix<double> fnbasis;


  boost::numeric::ublas::matrix<complex<double> > cZibl;
  boost::numeric::ublas::matrix<complex<double> > cZjbl;
  boost::numeric::ublas::matrix<complex<double> > cfnbasis;


  boost::numeric::ublas::matrix<double> fnw;
  boost::numeric::ublas::matrix<complex<double> > cfnw;

  // Check that eigenvectors are properly normalized!

  boost::numeric::ublas::matrix<complex<double> > cOneTest;


  // Check time
  boost::timer t;
  double time_elapsed;

  boost::timer t2;
  double time_lap2;

  
  CNRGmatrix* MatOLD = new CNRGmatrix [NumNRGMats];

  //Clear all
  for (int imats=0;imats<NumNRGMats;imats++){
    if (NRGMats[imats].NeedOld) MatOLD[imats]=NRGMats[imats];
    NRGMats[imats].ClearAll();
    NRGMats[imats].SyncNRGarray(*pAeigCut);
  }
  // end for

  // Note: assumes  Abasis and Aeig have 
  // EXACTLY the same block structure... should check.

  // Loop over blocks (icount counts for each matrix)
  int icount[NumNRGMats];
  for (int ii=0;ii<NumNRGMats;ii++) icount[ii]=0;

  t2.restart();
  time_lap2=0.0;
  // Note Aeig can be a "cut" array
  for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++){
    int Nstibl=pAeigCut->GetBlockSize(ibl);

    cout << "Block i : " << ibl << "/" << pAeigCut->NumBlocks()-1
	 << " (sz=" <<  Nstibl << ") Nshell=" << pAeigCut->Nshell << endl;
    // Corresponding block size in pAbasis
    int iblBC=FindMatchBlock(pAeigCut,ibl,pAbasis);
    int NstiblBC=pAbasis->GetBlockSize(iblBC);


    for (int jbl=0;jbl<pAeigCut->NumBlocks();jbl++){
      int Nstjbl=pAeigCut->GetBlockSize(jbl);
	  
      cout << "  Block j : " << jbl << "/" << pAeigCut->NumBlocks()-1
	   << " (sz=" <<  Nstjbl 
	   << ") Time to here " << t2.elapsed() 
	   << "; last lap " << t2.elapsed()-time_lap2 << endl;
      time_lap2=t2.elapsed();

      int jblBC=FindMatchBlock(pAeigCut,jbl,pAbasis);
      int NstjblBC=pAbasis->GetBlockSize(jblBC);


      //Loop matrices: check if the matrix elements are non-zero
      for (int imats=0;imats<NumNRGMats;imats++){
	if ( NRGMats[imats].CheckForMatEl(pAeigCut,ibl,jbl)   ){
	  cout << "   Updating Op imats=: " << imats << " of " << NumNRGMats-1 << endl;
	  // Set Zi,Zj
	  if (display){
	    cout << " Nshell = " << pAeigCut->Nshell << endl;
	    //cout << " Updating Op imats=: " << imats << " of " << NumNRGMats-1 << endl;
	    cout << " icount = " << icount[imats] << " to " 
		 << icount[imats]+Nstibl*Nstjbl-1 << endl;
	    cout << "  Setting up Block i : " << ibl 
		 << " (size " <<  Nstibl 
		 << "  was bl " << iblBC << " sz " << NstiblBC
		 << ") x Block j " << jbl
		 << " (size " << Nstjbl 
		 << "  was bl " << jblBC << " sz " << NstjblBC
		 << ") of " << pAeigCut->NumBlocks() << endl;
	  }
	  NRGMats[imats].MatBlockMap.push_back(ibl);
	  NRGMats[imats].MatBlockMap.push_back(jbl);
	  NRGMats[imats].MatBlockBegEnd.push_back(icount[imats]);
	  icount[imats]+=Nstibl*Nstjbl;
	  NRGMats[imats].MatBlockBegEnd.push_back(icount[imats]-1);

	  // 		  cout << "    Setting up BLAS matrices..." << endl;
	  //boost::numeric::ublas::matrix<double> Zibl(Nstibl,Nstibl);
	  //boost::numeric::ublas::matrix<double> Zjbl(Nstjbl,Nstjbl);
	  if (NRGMats[imats].IsComplex){
	    cZibl.resize(Nstibl,NstiblBC);
	    cZjbl.resize(Nstjbl,NstjblBC);
	    cZibl=pAeigCut->cEigVecCut2BLAS(ibl);
	    cZjbl=pAeigCut->cEigVecCut2BLAS(jbl);
	    cfnbasis.resize(NstiblBC,NstjblBC);
	  }else{
	    Zibl.resize(Nstibl,NstiblBC);
	    Zjbl.resize(Nstjbl,NstjblBC);
	    Zibl=pAeigCut->EigVecCut2BLAS(ibl);
	    Zjbl=pAeigCut->EigVecCut2BLAS(jbl);
	    fnbasis.resize(NstiblBC,NstjblBC);
	  }

	  if (display) cout << "    ...Zs done ..." << endl;

// 	  fnbasis.resize(NstiblBC,NstjblBC);
	  // Set-up fn basis: Loop in the basis states
	  int istbl=0;
	  for (int ist=pAbasis->GetBlockLimit(iblBC,0);
	       ist<=pAbasis->GetBlockLimit(iblBC,1);ist++){
	    //int typei=pAbasis->iType[ist];
	    //int stcfi=pAbasis->StCameFrom[ist];

	    int jstbl=0;
	    for (int jst=pAbasis->GetBlockLimit(jblBC,0);
		 jst<=pAbasis->GetBlockLimit(jblBC,1);jst++){
	      //int typej=pAbasis->iType[jst];
	      //int stcfj=pAbasis->StCameFrom[jst];

	      	  
	      if (NRGMats[imats].IsComplex) 
		cfnbasis(istbl,jstbl)=ZeroC;
	      else fnbasis(istbl,jstbl)=0.0;

	      //Calculate fbasis
	      if (NRGMats[imats].NeedOld){
		int iold=pAbasis->StCameFrom[ist];
		int jold=pAbasis->StCameFrom[jst];
		if (NRGMats[imats].IsComplex){
		  cfnbasis(istbl,jstbl)=(NRGMats[imats].CalcMatElCplx(pAbasis,pSingleSite,ist,jst))*MatOLD[imats].cGetMatEl(iold,jold);
		}
		  else
		fnbasis(istbl,jstbl)=(NRGMats[imats].CalcMatEl(pAbasis,pSingleSite,ist,jst))*MatOLD[imats].GetMatEl(iold,jold);
	      }
	      else{
		if (NRGMats[imats].IsComplex)
		  cfnbasis(istbl,jstbl)=NRGMats[imats].CalcMatElCplx(pAbasis,pSingleSite,ist,jst);
 
		else
		  fnbasis(istbl,jstbl)=NRGMats[imats].CalcMatEl(pAbasis,pSingleSite,ist,jst);
	      }
	      jstbl++;
	    }
	    // end jst loop
	    istbl++;
	  }
	  // end ist loop
	  if (display) cout << "    ...fnbasis done." << endl;

	  if (display) cout << "    Multiplying BLAS matrices (final size:" 
			    << Nstibl << " x " << Nstjbl << ") ... " << endl;

	  t.restart();
	  if (NRGMats[imats].IsComplex){
	    cfnw.resize(Nstibl,Nstjbl);
// 	    noalias(cfnw)=prod (cZibl, 
// 		     boost::numeric::ublas::matrix<complex<double> >(prod(cfnbasis,trans(cZjbl))) );
	    // Wrong. It should be:
	    noalias(cfnw)=prod (cZibl, 
		     boost::numeric::ublas::matrix<complex<double> >(prod(cfnbasis,herm(cZjbl))) );

	  }
	  else{
	    fnw.resize(Nstibl,Nstjbl);
	    noalias(fnw)=prod (Zibl, 
		     boost::numeric::ublas::matrix<double>(prod(fnbasis,trans(Zjbl))) );
	  }
	  time_elapsed=t.elapsed();
	  if (display) cout << "    ...done. Elapsed time:" << time_elapsed << endl;
	  // Debug
	  // if ( (display)&&(pAeigCut->Nshell==1)&&
	  //      ( (ibl==1)||(jbl==1) )&&
	  //      ( (imats==1) )
	  //      ) {
	  //   cout << "Z(i="<<ibl<<")  : " <<  Zibl << endl;
	  //   cout << "Z(j="<<jbl<<")  : " <<  Zjbl << endl;
	  //   cout << "fbasis   :" <<  fnbasis << endl;
	  //   cout << "Zi.fbasis.ZjT : " <<  fnw << endl;
	  //   ///
	  //   // cout << "cZ(i="<<ibl<<")  : " <<  cZibl << endl;
	  //   // cout << "cZ(j="<<jbl<<")  : " <<  cZjbl << endl;

	  //   // cOneTest.resize(NstiblBC,NstiblBC);
	  //   // noalias(cOneTest)=prod (herm(cZibl), cZibl);    

	  //   // cout << "cZi+.cZi   : " <<  cOneTest << endl;

	  //   // cOneTest.resize(NstjblBC,NstjblBC);
	  //   // noalias(cOneTest)=prod (herm(cZjbl), cZjbl);    

	  //   // cout << "cZj+.cZj   : " <<  cOneTest << endl;

	  //   // cout << "cfbasis   :" <<  cfnbasis << endl;
	  //   // cout << "cZi.cfbasis.cZjT : " <<  cfnw << endl;
	  //   //cout << "faux   :" <<  faux << endl;
	  //  }
	  // end debug

	  // Add to NRGMats[imats-1]
	  if (NRGMats[imats].IsComplex){
	    for (int ii=0;ii<cfnw.size1();ii++)
	      for (int jj=0;jj<cfnw.size2();jj++)
		NRGMats[imats].MatElCplx.push_back(cfnw(ii,jj));
	  }
	  else{
	    for (int ii=0;ii<fnw.size1();ii++)
	      for (int jj=0;jj<fnw.size2();jj++)
		NRGMats[imats].MatEl.push_back(fnw(ii,jj));
	  }
	  // end if is complex 
	}
	//end if Q=Q'+1 etc

      }
      // end imats loop



    }
    // end jbl loop

  }
  // end ibl loop

  // Release MatOLD
  delete[] MatOLD;

  cout << " DONE updating operators in Nshell = " << pAeigCut->Nshell << endl; 

}
// END subroutine

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

void UpdateMatrices(CNRGbasisarray* pSingleSite,CNRGbasisarray* pAeigCut, 
		    CNRGbasisarray* pAbasis,
		    CNRGmatrix* NRGMats, int NumNRGMats, bool display){

  // STL vectors 
  vector<double> fnbasis;
  vector<complex<double> > cfnbasis;

  vector<double> fnw;
  vector<complex<double> > cfnw;

  vector<double> auxMat1;
  vector<complex<double> > cauxMat1;

  
  double auxEl;
  complex<double> cauxEl;

  int ibegZibl,ibegZjbl; // Point where the blocks begin

  int ldaBC,lda,ldaaux; //=max(SizeBl,SizeBl2);

  // cblas_dgemm variables
  double ALPHA=1.0, BETA=0.0;
  //   C <- ALPHA*A.B + BETA*C

  //typedef CBLAS_ORDER CBLAS_LAYOUT; // this for backward compatibility 

  // Check that eigenvectors are properly normalized!

  //boost::numeric::ublas::matrix<complex<double> > cOneTest;


  // Check time
  boost::timer t;
  double time_elapsed;

  boost::timer t2;
  double time_lap2;

  
  CNRGmatrix* MatOLD = new CNRGmatrix [NumNRGMats];

  //Clear all
  for (int imats=0;imats<NumNRGMats;imats++){
    if (NRGMats[imats].NeedOld) MatOLD[imats]=NRGMats[imats];
    NRGMats[imats].ClearAll();
    NRGMats[imats].SyncNRGarray(*pAeigCut);
  }
  // end for

  // Note: assumes  Abasis and Aeig have 
  // EXACTLY the same block structure... should check.

  // Loop over blocks (icount counts for each matrix)
  int icount[NumNRGMats];
  for (int ii=0;ii<NumNRGMats;ii++) icount[ii]=0;

  t2.restart();
  time_lap2=0.0;
  // Note Aeig can be a "cut" array
  for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++){
    int Nstibl=pAeigCut->GetBlockSize(ibl);

    cout << "Block i : " << ibl << "/" << pAeigCut->NumBlocks()-1
	 << " (sz=" <<  Nstibl << ") Nshell=" << pAeigCut->Nshell << endl;
    // Corresponding block size in pAbasis
    int iblBC=FindMatchBlock(pAeigCut,ibl,pAbasis);
    int NstiblBC=pAbasis->GetBlockSize(iblBC);


    for (int jbl=0;jbl<pAeigCut->NumBlocks();jbl++){
      int Nstjbl=pAeigCut->GetBlockSize(jbl);
	  
      cout << "  Block j : " << jbl << "/" << pAeigCut->NumBlocks()-1
	   << " (sz=" <<  Nstjbl 
	   << ") Time to here " << t2.elapsed() 
	   << "; last lap " << t2.elapsed()-time_lap2 << endl;
      time_lap2=t2.elapsed();

      int jblBC=FindMatchBlock(pAeigCut,jbl,pAbasis);
      int NstjblBC=pAbasis->GetBlockSize(jblBC);


      //Loop matrices: check if the matrix elements are non-zero
      for (int imats=0;imats<NumNRGMats;imats++){
	if ( NRGMats[imats].CheckForMatEl(pAeigCut,ibl,jbl)   ){
	  cout << "   Updating Op imats=: " << imats << " of " << NumNRGMats-1 << endl;
	  // Set Zi,Zj
	  if (display){
	    cout << " Nshell = " << pAeigCut->Nshell << endl;
	    //cout << " Updating Op imats=: " << imats << " of " << NumNRGMats-1 << endl;
	    cout << " icount = " << icount[imats] << " to " 
		 << icount[imats]+Nstibl*Nstjbl-1 << endl;
	    cout << "  Setting up Block i : " << ibl 
		 << " (size " <<  Nstibl 
		 << "  was bl " << iblBC << " sz " << NstiblBC
		 << ") x Block j " << jbl
		 << " (size " << Nstjbl 
		 << "  was bl " << jblBC << " sz " << NstjblBC
		 << ") of " << pAeigCut->NumBlocks() << endl;
	  }
	  NRGMats[imats].MatBlockMap.push_back(ibl);
	  NRGMats[imats].MatBlockMap.push_back(jbl);
	  NRGMats[imats].MatBlockBegEnd.push_back(icount[imats]);
	  icount[imats]+=Nstibl*Nstjbl;
	  NRGMats[imats].MatBlockBegEnd.push_back(icount[imats]-1);

	  ibegZibl=pAeigCut->GetBlockLimitEigv(ibl,0);
	  ibegZjbl=pAeigCut->GetBlockLimitEigv(jbl,0);

	  cfnbasis.clear();
	  fnbasis.clear();	  

	  if (display) cout << "    ...Zs done ..." << endl;

	  // Set-up fn basis: Loop in the basis states
	  int istbl=0;
	  for (int ist=pAbasis->GetBlockLimit(iblBC,0);
	       ist<=pAbasis->GetBlockLimit(iblBC,1);ist++){
	    //int typei=pAbasis->iType[ist];
	    //int stcfi=pAbasis->StCameFrom[ist];

	    int jstbl=0;
	    for (int jst=pAbasis->GetBlockLimit(jblBC,0);
		 jst<=pAbasis->GetBlockLimit(jblBC,1);jst++){
	      //int typej=pAbasis->iType[jst];
	      //int stcfj=pAbasis->StCameFrom[jst];

	      //Calculate fbasis
	      if (NRGMats[imats].NeedOld){
		int iold=pAbasis->StCameFrom[ist];
		int jold=pAbasis->StCameFrom[jst];
		if (NRGMats[imats].IsComplex){
		  cauxEl=(NRGMats[imats].CalcMatElCplx(pAbasis,pSingleSite,ist,jst))*MatOLD[imats].cGetMatEl(iold,jold);
		  cfnbasis.push_back(cauxEl);
		}
		else{
		  auxEl=(NRGMats[imats].CalcMatEl(pAbasis,pSingleSite,ist,jst))*MatOLD[imats].GetMatEl(iold,jold);
		  fnbasis.push_back(auxEl);
		}
	      }
	      else{
		if (NRGMats[imats].IsComplex){
		  cauxEl=NRGMats[imats].CalcMatElCplx(pAbasis,pSingleSite,ist,jst);
		  cfnbasis.push_back(cauxEl);
		}
		else{
		  auxEl=NRGMats[imats].CalcMatEl(pAbasis,pSingleSite,ist,jst);
		  fnbasis.push_back(auxEl);
		}
	      }// end if NeedOld else
	      jstbl++;
	    }
	    // end jst loop
	    istbl++;
	  }
	  // end ist loop
	  if (display) cout << "    ...fnbasis done." << endl;

	  if (display) cout << "    Multiplying BLAS matrices (final size:" 
			    << Nstibl << " x " << Nstjbl << ") ... " << endl;

	  t.restart();

	  // Let's see:
	  // LDA-> No. of cols of the ORIGINAL matrix
	  //
	  //   Zibl -> (Nstibl x NstiblBC) matrix -> LDA=NstiblBC
	  //   Zjbl -> (Nstjbl x NstjblBC) matrix -> LDA=NstjblBC
	  //   fnbasis -> (NstiblBC x NstjblBC) matrix -> LDA=NstjblBC
	  //   fnw -> (Nstibl x Nstjbl) matrix    -> LDA=Nstjbl
	  // auxMat1 -> (NstiblBC x Nstjbl) matrix -> LDA=Nstjbl
	  // auxMat1 = cfnbasis_(NstiblBC x NstjblBC) . herm(cZjbl)_(NstjblBC x Nstjbl)
	  // 
	  //  -> (NstiblBC x Nstjbl) matrix
	  // cfnw = cZibl_(Nstibl x NstiblBC) . auxMat1_(NstiblBC x Nstjbl) 
	  
	  if (NRGMats[imats].IsComplex){

	    cauxMat1.clear();
	    cauxMat1.resize(NstiblBC*Nstjbl,ZeroC);

	    // fnbasis_NstiblBC x NstjblBC) x [Zjbl_(Nstjbl x NstjblBC)]^+
	    // -> cMat1_(NstiblBC x Nstjbl)
	    cblas_zgemm(CblasRowMajor,  CblasNoTrans, CblasConjTrans,
			NstiblBC, Nstjbl, NstjblBC,
			&ALPHA, &cfnbasis[0], NstjblBC, &(pAeigCut->cEigVec[ibegZjbl]), NstjblBC,
			&BETA, &cauxMat1[0], Nstjbl);

	    
	    cfnw.clear();
	    cfnw.resize(Nstibl*Nstjbl,ZeroC);


	    //  Zibl_(Nstibl x NstiblBC) x cMat1_(NstiblBC x Nstjbl)
	    //  -> fnw_(Nstibl x Nstjbl)
	    cblas_zgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans,
			Nstibl, Nstjbl, NstiblBC,
			&ALPHA, &(pAeigCut->cEigVec[ibegZibl]), NstiblBC,
			&cauxMat1[0], Nstjbl, &BETA, &cfnw[0], Nstjbl );

	  }else{

	    auxMat1.clear();
	    auxMat1.resize(NstiblBC*Nstjbl,0.0);

	    // fnbasis_NstiblBC x NstjblBC) x [Zjbl_(Nstjbl x NstjblBC)]^+
	    // -> cMat1_(NstiblBC x Nstjbl)
	    cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasConjTrans,
			NstiblBC, Nstjbl, NstjblBC,
			ALPHA, &fnbasis[0], NstjblBC, &(pAeigCut->dEigVec[ibegZjbl]), NstjblBC,
			BETA, &auxMat1[0], Nstjbl);

	    fnw.clear();
	    fnw.resize(Nstibl*Nstjbl,0.0);

	    //  Zibl_(Nstibl x NstiblBC) x auxMat1_(NstiblBC x Nstjbl)
	    //  -> fnw_(Nstibl x Nstjbl)
	    cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans,
			Nstibl, Nstjbl, NstiblBC,
			ALPHA, &(pAeigCut->dEigVec[ibegZibl]), NstiblBC,
			&auxMat1[0], Nstjbl, BETA, &fnw[0], Nstjbl );


	  }
	  time_elapsed=t.elapsed();
	  if (display) cout << "    ...done. Elapsed time:" << time_elapsed << endl;
	  // Debug
	  // if ( (display)&&(pAeigCut->Nshell==1)&&
	  //      ( (ibl==1)||(jbl==1) )&&
	  //      ( (imats==1) )
	  //      ) {
	  //   cout << " ibl= " << ibl
	  // 	 << " jbl= " << jbl
	  // 	 << endl;
	  //   cout << " Nstibl=" << Nstibl
	  // 	 << " NstiblBC=" << NstiblBC
	  // 	 << " Nstjbl=" << Nstjbl
	  // 	 << " NstjblBC=" << NstjblBC
	  // 	 << endl;
	  //   cout << " fnbasis = ";
	  //   for (int ii=0; ii<fnbasis.size(); ii++){
	  //     if ( ii % NstjblBC ==0) cout << endl;
	  //     cout << fnbasis[ii] << " ";
	  //   }
	  //   cout<< endl;
	  //   cout << " fnw = ";
	  //   for (int ii=0; ii<fnw.size(); ii++){
	  //     if ( ii % Nstjbl ==0) cout << endl;
	  //     cout << fnw[ii] << " ";
	  //   }
	  //   cout<< endl;

	    // cout << " cfnbasis = ";
	    // for (int ii=0; ii<cfnbasis.size(); ii++){
	    //   if ( ii % NstjblBC ==0) cout << endl;
	    //   cout << cfnbasis[ii] << " ";
	    // }
	    // cout<< endl;
	    // cout << " cfnw = ";
	    // for (int ii=0; ii<cfnw.size(); ii++){
	    //   if ( ii % Nstjbl ==0) cout << endl;
	    //   cout << cfnw[ii] << " ";
	    // }
	    // cout<< endl;

	  // }
	// end debug

	  // Add to NRGMats[imats-1]
	  if (NRGMats[imats].IsComplex){
	    for (int ii=0;ii<cfnw.size();ii++)
	      NRGMats[imats].MatElCplx.push_back(cfnw[ii]);
	  }
	  else{
	    for (int ii=0;ii<fnw.size();ii++)
	      NRGMats[imats].MatEl.push_back(fnw[ii]);
	  }
	  // end if is complex 
	}
	//end if Q=Q'+1 etc

      }
      // end imats loop

    }
    // end jbl loop
  }
  // end ibl loop

  // Release MatOLD
  delete[] MatOLD;

  cout << " DONE updating operators in Nshell = " << pAeigCut->Nshell << endl; 

}
// END subroutine



//////////////////////////////////
//////////////////////////////////
//////////////////////////////////


