
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>


#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

using namespace std;


////////////////////////////////////
// Class CNRGmatrix functions //
////////////////////////////////////


//////////////
int CNRGmatrix::NumMatBlocks(){

  // Needs MatBlockMap updated
  
  return(MatBlockBegEnd.size()/2);

}


//////////////
int CNRGmatrix::FindMatBlock(int iblock1, int iblock2){

  int iMatBlock=0;

  int imatch[]={iblock1,iblock2};

  vector<int>::iterator it;

  // STL search algorithm

  it =search(MatBlockMap.begin(),MatBlockMap.end(),imatch,imatch+2 );

  int it_pos=int(it-MatBlockMap.begin());
  if ( (it!= MatBlockMap.end())&&((it_pos % 2)==0) )
    {
      iMatBlock=it_pos/2;
//    cout <<  "Found iblock1 = "<< iblock1 << " iblock2 = "<< iblock2 
//         << "at " << iMatBlock << endl;
    }
  else
    {
      iMatBlock=-1;
//       cout <<  "Cannot find iblock1 = "<< iblock1 << " iblock2 = "<< iblock2 << endl;
    }

  return(iMatBlock);


}

//////////////
int CNRGmatrix::GetMatBlockLimit(int imatblock, int whichlimit){


  return(MatBlockBegEnd[2*imatblock+whichlimit]);


}

//////////////
int CNRGmatrix::GetMatBlockLimit(int iblock1, int iblock2, int whichlimit){

  int ii=FindMatBlock(iblock1,iblock2);

  if (ii<0)
    return(0);
  else
  return(MatBlockBegEnd[2*ii+whichlimit]);


}

//////////////
int CNRGmatrix::GetMatBlockSize(int imatblock ){


  return(MatBlockBegEnd[2*imatblock+1]-MatBlockBegEnd[2*imatblock]+1);

}



//////////////
int CNRGmatrix::GetMatBlockSize(int iblock1, int iblock2){


  int ii=FindMatBlock(iblock1,iblock2);

  if (ii<0)
    return(0);
  else
  return(MatBlockBegEnd[2*ii+1]-MatBlockBegEnd[2*ii]+1);

}

//////////////
double CNRGmatrix::GetBlockMatEl(int iblock1, int iblock2, int iel, int jel){

  int iMatBlock=FindMatBlock(iblock1,iblock2);

  //if (iMatBlock<0)
  // Watch out: if BlockBegEnd has "holes", iel or jel can be negative (state is not there)
  if ( (iMatBlock<0)||(iel<0)||(jel<0) )
      return(0.0);
  else{
    int i0=GetMatBlockLimit(iMatBlock, 0);
    int Nbl1=CNRGarray::GetBlockSize(iblock1);
    int Nbl2=CNRGarray::GetBlockSize(iblock2);
    int rbl=ij2rNSq(Nbl1,Nbl2,iel,jel);
    
    // watch out for this one!! 
    // It was commented... WHY???? April 09.
    if ((UpperTriangular)&&(iblock1==iblock2)){
      if (iel>jel)
	rbl=ij2rUpTr(Nbl1,jel,iel);
      else
	rbl=ij2rUpTr(Nbl1,iel,jel);
    }

    return(MatEl[i0+rbl]);
  }
}
// complex
complex<double> CNRGmatrix::cGetBlockMatEl(int iblock1, int iblock2, int iel, int jel){

  int iMatBlock=FindMatBlock(iblock1,iblock2);

  //if (iMatBlock<0)
  // Watch out: if BlockBegEnd has "holes", iel or jel can be negative (state is not there)
  if ( (iMatBlock<0)||(iel<0)||(jel<0) )
      return(ZeroC);
  else{
    int i0=GetMatBlockLimit(iMatBlock, 0);
    int Nbl1=CNRGarray::GetBlockSize(iblock1);
    int Nbl2=CNRGarray::GetBlockSize(iblock2);
    int rbl=ij2rNSq(Nbl1,Nbl2,iel,jel);
    
    // watch out for this one!! 
    // It was commented... WHY???? April 09.
    if ((UpperTriangular)&&(iblock1==iblock2)){
      if (iel>jel)
	rbl=ij2rUpTr(Nbl1,jel,iel);
      else
	rbl=ij2rUpTr(Nbl1,iel,jel);
    }
    return(MatElCplx[i0+rbl]);
  }
}



//////////////
double CNRGmatrix::GetMatEl(int ist, int jst){

  int ibl,jbl;

  int iblock1=CNRGarray::GetBlockFromSt(ist,ibl);
  int iblock2=CNRGarray::GetBlockFromSt(jst,jbl);

  return(GetBlockMatEl(iblock1,iblock2,ibl,jbl));
}
// complex
complex<double> CNRGmatrix::cGetMatEl(int ist, int jst){

  int ibl,jbl;

  int iblock1=CNRGarray::GetBlockFromSt(ist,ibl);
  int iblock2=CNRGarray::GetBlockFromSt(jst,jbl);

  return(cGetBlockMatEl(iblock1,iblock2,ibl,jbl));

}

//////////////
// Replace existing matrix element by El

int CNRGmatrix::GetMatElPosition(int iblock1, int iblock2, int iel, int jel){

  int iMatBlock=FindMatBlock(iblock1,iblock2);

  // Watch out: if BlockBegEnd has "holes", iel or jel can be negative (state is not there)

  // iel,jel index states within the block (starting from 0)

  if ( (iMatBlock<0)||(iel<0)||(jel<0) ){
    // Message is harmless.
    //cout << " CNRGmatrix::GetMatElPosition: iMatBlock not found " << endl; 
    return(-1);
  }
  else{
    int i0=GetMatBlockLimit(iMatBlock, 0);
    int Nbl1=CNRGarray::GetBlockSize(iblock1);
    int Nbl2=CNRGarray::GetBlockSize(iblock2);
    int rbl=ij2rNSq(Nbl1,Nbl2,iel,jel);
    
    // watch out for this one!! 
    if ((UpperTriangular)&&(iblock1==iblock2)){
      if (iel>jel)
	rbl=ij2rUpTr(Nbl1,jel,iel);
      else
	rbl=ij2rUpTr(Nbl1,iel,jel);
    }
    return(i0+rbl);
  }

}
///////

void CNRGmatrix::PushBlockMatEl(double El, int iblock1, int iblock2, int iel, int jel){
  // Deprecated by GetMatElPosition(Oct 2013)
  int iMatBlock=FindMatBlock(iblock1,iblock2);

  // Watch out: if BlockBegEnd has "holes", iel or jel can be negative (state is not there)

  // iel,jel index states within the block (starting from 0)

  if ( (iMatBlock<0)||(iel<0)||(jel<0) ){
    // Message is harmless.
    //cout << " CNRGmatrix::PushBlockMatEl: iMatBlock not found " << endl; 
    return;
  }
  else{
    int i0=GetMatBlockLimit(iMatBlock, 0);
    int Nbl1=CNRGarray::GetBlockSize(iblock1);
    int Nbl2=CNRGarray::GetBlockSize(iblock2);
    int rbl=ij2rNSq(Nbl1,Nbl2,iel,jel);
    
    // watch out for this one!! 
    if ((UpperTriangular)&&(iblock1==iblock2)){
      if (iel>jel)
	rbl=ij2rUpTr(Nbl1,jel,iel);
      else
	rbl=ij2rUpTr(Nbl1,iel,jel);
    }
//     MatEl[i0+rbl]=El;
  //if (IsComplex) {MatElCplx[i0+rbl].real()=El;MatElCplx[i0+rbl].imag()=0.0;}
    if (IsComplex) MatElCplx[i0+rbl]=(complex<double>)El;
    else MatEl[i0+rbl]=El;
  }


}

//////////////
// Replace EXISTING matrix element ist,jst by El
void CNRGmatrix::PushMatEl(double El, int ist, int jst){
  int ibl,jbl;

  int iblock1=CNRGarray::GetBlockFromSt(ist,ibl);
  int iblock2=CNRGarray::GetBlockFromSt(jst,jbl);

  int iPos=GetMatElPosition(iblock1,iblock2,ibl,jbl);

//   if (iPos>-1) MatEl[iPos]=El;
  if (iPos>-1)
    if (IsComplex){MatElCplx[iPos]=(complex<double>)El;}
    else MatEl[iPos]=El;

  //PushBlockMatEl(El,iblock1,iblock2,ibl,jbl);

}
void CNRGmatrix::PushMatEl(complex<double> El, int ist, int jst){
  int ibl,jbl;

  int iblock1=CNRGarray::GetBlockFromSt(ist,ibl);
  int iblock2=CNRGarray::GetBlockFromSt(jst,jbl);

  int iPos=GetMatElPosition(iblock1,iblock2,ibl,jbl);

  if (iPos>-1) MatElCplx[iPos]=El;
}
///////////////

///////////// TEMPLATE TESTING /////////////
// Needs to be defined in header file. Otherwise won't work
//
///////////// END TEMPLATE TESTING ////////////////////



void CNRGmatrix::FilterMap_SetBegEnd(){

// Original
// [e1 e2 e3 e4 e5]
// [(b1 b2), (b1 b2), (b4 b5), (b4 b5), (b5 b6) ]
//
// End up with
//
// [e1 e2 e3 e4 e5]
//[(b1 b2) (b4 b5) (b5 b6) ]
//[(0 1),(2 3), (4 4)] 
// 

//
// Note: needs that MatBlockMap entering to have two indexes *per element*
//

  vector<int>::iterator map_iter,map_iter2;
  vector<int>::iterator it;
   int ipos=0;
   int nelblock=0;

   MatBlockBegEnd.clear();


   cout << "Setting up MatBlockBegEnd: " << endl;

   for (map_iter=MatBlockMap.begin(); map_iter<MatBlockMap.end(); 
       map_iter+=2)
     {
       MatBlockBegEnd.push_back(ipos);


       nelblock=1;
       map_iter2=map_iter+2;
       it=search(map_iter2,MatBlockMap.end(),
		 map_iter,map_iter+1);
       // Search for sequence map_iter,map_iter+1
//        while (it!=MatBlockMap.end())
// 	 {
// 	   cout << "Found repetition. MatBlockMap size =  " 
// 		<< MatBlockMap.size() << endl;
// 	   nelblock++;
// 	   MatBlockMap.erase(it,it+2);
// 	   it=search(map_iter2,MatBlockMap.end(),
// 		     map_iter,map_iter+1);
// 	 }

       while (map_iter2<MatBlockMap.end())
	 {
	   if ( (*map_iter2==*map_iter)&&(*(map_iter2+1)==*(map_iter+1))  )
	     {
	       nelblock++;
	       MatBlockMap.erase(map_iter2,map_iter2+2);
	       map_iter2-=2;

	     }
	   map_iter2+=2;
	 }


       ipos+=nelblock;

       MatBlockBegEnd.push_back(ipos-1);

     }

}

//////////////
// double* Block2Array(int imatblock){


//   return();
// }

//////////////
void CNRGmatrix::ClearAll(){

  CNRGarray::ClearAll();
  MatEl.clear();
  MatElCplx.clear();
  MatBlockMap.clear();
  MatBlockBegEnd.clear();
  

}
//////////////

//////////////
void CNRGmatrix::PrintMatBlockQNumbers(int imatblock){


  int ibl=MatBlockMap[2*imatblock]; 
  int jbl=MatBlockMap[2*imatblock+1]; 

  cout << " Block (" << ibl << ") x Block (" << jbl << ") : ";
  for (int ii=0;ii<NQNumbers;ii++) cout << fixed << GetQNumber(ibl,ii) << " ";  
  cout << "  x  ";
  for (int ii=0;ii<NQNumbers;ii++) cout << fixed << GetQNumber(jbl,ii) << " ";  
  cout << endl;

  cout << "Nsti = " << CNRGarray::GetBlockLimit(ibl,0) << " to " << CNRGarray::GetBlockLimit(ibl,1) << "   ;   "; 
  cout << "Nstj = " << CNRGarray::GetBlockLimit(jbl,0) << " to " << CNRGarray::GetBlockLimit(jbl,1) << endl; 

  cout << resetiosflags (ios_base::floatfield);
  
}

//////////////
void CNRGmatrix::PrintMatBlock(int imatblock){

  int iblock1=MatBlockMap[2*imatblock];
  int iblock2=MatBlockMap[2*imatblock+1];

  int SizeBl1=GetBlockSize(iblock1);
  int SizeBl2=GetBlockSize(iblock2);

  int i0=GetMatBlockLimit(imatblock,0);
  int i1=i0;

  double Mij;
  complex<double> Mij_cplx;
  int pos=0;

  if (QNumbers.size()>0)
    PrintMatBlockQNumbers(imatblock);

  for (int ii=0;ii<SizeBl1;ii++){
    for (int jj=0;jj<SizeBl2;jj++){
      if (UpperTriangular){
	if (ii>jj){
	  pos=i0+ij2rUpTr(SizeBl1,jj,ii);
	  if (IsComplex) Mij_cplx=conj(MatElCplx[pos]);
	  }
	else{
	  pos=i0+ij2rUpTr(SizeBl1,ii,jj);
	  if (IsComplex) Mij_cplx=MatElCplx[pos];
	}
	if (!IsComplex) Mij=MatEl[pos];
      }
      else{
	if (!IsComplex)Mij=MatEl[i1];
	else{
	  // Is this correct? I don't think so...
	  // This is NOT upper triangular!!
// 	  if (ii>jj) Mij_cplx=conj(MatElCplx[i1]);
// 	  else Mij_cplx=MatElCplx[i1];
 	  Mij_cplx=MatElCplx[i1];
	}
	i1++;
      }
      if (!IsComplex)
      cout << Mij  << " ";
      else
      cout << Mij_cplx  << " ";
    }
    cout << endl;
  }
  cout << endl;

  
}

//////////////
void CNRGmatrix::PrintMatBlock(int iblock1, int iblock2){

  int imatblock=FindMatBlock(iblock1,iblock2);

  int SizeBl1=GetBlockSize(iblock1);
  int SizeBl2=GetBlockSize(iblock2);

  int i0=GetMatBlockLimit(imatblock,0);
  int i1=i0;

  double Mij;
  complex<double> Mij_cplx;
  int pos=0;

  if (QNumbers.size()>0)
    PrintMatBlockQNumbers(imatblock);
  for (int ii=0;ii<SizeBl1;ii++){
    for (int jj=0;jj<SizeBl2;jj++){
      if (UpperTriangular){
	if (ii>jj){
	  pos=i0+ij2rUpTr(SizeBl1,jj,ii);
	  if (IsComplex) Mij_cplx=conj(MatElCplx[pos]);
	}
	else{
	  pos=i0+ij2rUpTr(SizeBl1,ii,jj);
	  if (IsComplex) Mij_cplx=MatElCplx[pos];
	}
	if (!IsComplex) Mij=MatEl[pos];
      }
      else{
	if (!IsComplex)Mij=MatEl[i1];
	else{
	  // Is this correct? I don't think so...
	  // This is NOT upper triangular!!
// 	  if (ii>jj) Mij_cplx=conj(MatElCplx[i1]);
// 	  else Mij_cplx=MatElCplx[i1];
 	  Mij_cplx=MatElCplx[i1];
	}
	i1++;
      }
      //cout << Mij  << " ";
      if (!IsComplex)
	cout << Mij  << " ";
      else
	cout << Mij_cplx  << " ";
    }
    cout << endl;
  }
  cout << endl;

}

//////////////
void CNRGmatrix::PrintAllBlocks(){

  for (int ibl=0;ibl<NumMatBlocks();ibl++)
    PrintMatBlock(ibl);

}


//////////////
void CNRGmatrix::DiagBlock(int iblock,
			   vector<double> &eigvalues, 
			   vector<double> &eigvectors){
// Regular
//
//   SUBROUTINE DSYEV( JOBZ, UPLO,	N, A, LDA, W, WORK, LWORK, INFO	)
//       CHARACTER	    JOBZ, UPLO
//       INTEGER	    INFO, LDA, LWORK, N
//       DOUBLE	    PRECISION A( LDA, *	), W( *	), WORK( * )
//
// UpperTriangular
//
//   SUBROUTINE DSPEV( JOBZ, UPLO,	N, AP, W, Z, LDZ, WORK,	INFO )
//       CHARACTER	    JOBZ, UPLO
//       INTEGER	    INFO, LDZ, N
//       DOUBLE	    PRECISION AP( * ), W( * ), WORK( * ), Z( LDZ, * )
//
//  AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//           On entry, the upper or lower triangle of the symmetric matrix
//           A, packed columnwise in a linear array.  The j-th column of A
//           is stored in the array AP as follows:
//           if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//           if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
//
//
// Z is DOUBLE PRECISION array, dimension (LDZ, N)
//           If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
//           eigenvectors of the matrix A, with the i-th column of Z
//           holding the eigenvector associated with W(i).
//           If JOBZ = 'N', then Z is not referenced.



  // LAPACK variables
  char JOBZ='V'; //'N'
  char UPLO='U';
  int Nst;
  int INFO=0,LDA,LWORK;
  int LDZ;

  // LDA=Nst
  // LWORK=3*Nst-1

  double *eigv;
  double *work;
  double *Zvec;

  // Get Block info
  int iMatBlock=FindMatBlock(iblock,iblock);

  if (iMatBlock<0){
    cout << "Cannot diagonalize block " << iblock << endl;
    return;
  }
  Nst=CNRGarray::GetBlockSize(iblock);
  //    cout << "Size block                 : " 
  //         << GetMatBlockSize(iMatBlock) 
  //         << " should be = " << CNRGarray::GetBlockSize(iblock) <<"^2"<< endl;
  LDA=Nst;
  LWORK=3*Nst-1;

  // Construct the STL vector. Remember: ().begin is an iterator

  vector<double> A1(MatEl.begin()+GetMatBlockLimit(iMatBlock,0),
		    MatEl.begin()+GetMatBlockLimit(iMatBlock,1)+1);

  

  eigv = new double [Nst];
  work = new double [LWORK];

  // Let me try this:
  if (UpperTriangular){
    cout << "Upper triangular" << endl;
    UPLO='L'; // LAPACK stores by cols, I store by lines: my 'U' is their 'L'
    LDZ=Nst;

    Zvec= new double [Nst*Nst];

    dspev_(&JOBZ,&UPLO,&Nst,&A1[0],eigv,Zvec,&LDZ,work,&INFO);

    int i1=0;
    for (int ii=0;ii<Nst;ii++){
      eigvalues.push_back(eigv[ii]);
      for (int jj=0;jj<Nst;jj++){
	eigvectors.push_back(Zvec[i1]);
	i1++;
      }
    }

    delete[] Zvec;

  }else{

    cout << "Regular block" << endl;

    // Diagonalization works on the STL vector directly!

    dsyev_(&JOBZ,&UPLO,&Nst,&A1[0],&LDA,eigv,work,&LWORK,&INFO);
    
    int i1=0;
    for (int ii=0;ii<Nst;ii++){
      eigvalues.push_back(eigv[ii]);
      for (int jj=0;jj<Nst;jj++){
	eigvectors.push_back(A1[i1]);
	i1++;
      }
    }

  }
  //else if UpTriang

  //cout << "Diag Block " << iblock << " Size: " << Nst << endl;
  //PrintMatBlockQNumbers(iMatBlock);
  //cout << "Lowest States :" << endl;
  //for (int ii=0;ii<( (20 < Nst)?20:Nst );ii++)
  //cout << eigv[ii] << endl;

  delete[] eigv;
  delete[] work;

  A1.clear();

  //Not needed
  //double *A; 
  //A = new double [Nst*Nst];
  //   for (int ii=0;ii < A1.size();ii++) A[ii]=A1[ii];
  //dsyev_(&JOBZ,&UPLO,&Nst,A,&LDA,eigv,work,&LWORK,&INFO);
  //delete[] A;

}

//////////////
//
//  Adding complex matrix! (Oct 2013)
//////////////
void CNRGmatrix::DiagBlock(int iblock,
			   vector<double> &eigvalues, 
			   vector<complex<double> > &eigvectors){
  // OLD

// Regular
//
// SUBROUTINE ZHEEV( JOBZ, UPLO,N, A, LDA, W, WORK, LWORK,RWORK, INFO)
//       CHARACTER	    JOBZ, UPLO
//       INTEGER	    INFO, LDA, LWORK, N
//       DOUBLE	    PRECISION RWORK( * ), W( * )
//       COMPLEX*16    A( LDA, * ), WORK( * )
//
// UpperTriangular
//
//   SUBROUTINE ZHPEV( JOBZ, UPLO,N, AP, W, Z, LDZ, WORK,RWORK, INFO )
//       CHARACTER	    JOBZ, UPLO
//       INTEGER	    INFO, LDZ, N
//       DOUBLE	    PRECISION RWORK( * ), W( * )
//       COMPLEX*16    AP(	* ), WORK( * ),	Z( LDZ,	* )
//
//AP is COMPLEX*16 array, dimension (N*(N+1)/2)
//           On entry, the upper or lower triangle of the Hermitian matrix
//           A, packed columnwise in a linear array.  The j-th column of A
//           is stored in the array AP as follows:
//           if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//           if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
// Z is COMPLEX*16 array, dimension (LDZ, N)
//           If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
//           eigenvectors of the matrix A, with the i-th column of Z
//           holding the eigenvector associated with W(i).


  // LAPACK variables
  char JOBZ='V'; //'N'
  char UPLO='U';
  int Nst;
  int INFO=0,LDA,LWORK,LRWORK;
  int LDZ;

  // LDA=Nst
  // LWORK=3*Nst-1

  double *eigv; //W
  complex<double> *work;
  double *rwork;
  complex<double> *Zvec;

  // Get Block info
  int iMatBlock=FindMatBlock(iblock,iblock);

  if (iMatBlock<0){
    cout << "Cannot diagonalize block " << iblock << endl;
    return;
  }
  Nst=CNRGarray::GetBlockSize(iblock);
  //    cout << "Size block                 : " 
  //         << GetMatBlockSize(iMatBlock) 
  //         << " should be = " << CNRGarray::GetBlockSize(iblock) <<"^2"<< endl;
  LDA=Nst;
  //LWORK=3*Nst-1;
  LWORK=std::max(1,2*Nst-1);
  LRWORK=std::max(1,3*Nst-2);

  // Construct the STL vector. Remember: ().begin is an iterator

  vector<complex<double> > A1(MatElCplx.begin()+GetMatBlockLimit(iMatBlock,0),
		    MatElCplx.begin()+GetMatBlockLimit(iMatBlock,1)+1);

  

  eigv = new double [Nst];
  rwork = new double [LRWORK];
  work = new complex<double> [LWORK];

  // Let me try this:
  if (UpperTriangular){
    cout << "Upper triangular" << endl;
    UPLO='L'; // LAPACK stores by cols, I store by lines: my 'U' is their 'L'
              // HOWEVER, I need to conjugte the off-diagonal elements first
              // Example: 4x4 matrix
              // Mine 'U': A11 A12 A13 A22 A23 A33
              // LAPACK 'L': A11 A12* A13* A22 A23* A33
              // LAPACK 'U' (wrong): A11 A12 A22 A13 A23 A33 
    LDZ=Nst;

    Zvec= new complex<double> [Nst*Nst];

    complex<double> caux=ZeroC;
    for (int ii=0; ii<A1.size();ii++){
      caux=std::conj(A1[ii]);
      A1[ii]=caux;
    }
  
//     cout << " Block : " << iblock;
//     for (int ii=0; ii<A1.size();ii++)
//       cout << " - A1(" << ii <<")= " << A1[ii] << endl;
    

    //dspev_(&JOBZ,&UPLO,&Nst,&A1[0],eigv,Zvec,&LDZ,work,&INFO);
    zhpev_(&JOBZ,&UPLO,&Nst,&A1[0],eigv,Zvec,&LDZ,work,rwork,&INFO);

    int i1=0;
    for (int ii=0;ii<Nst;ii++){
      eigvalues.push_back(eigv[ii]);
      for (int jj=0;jj<Nst;jj++){
 	eigvectors.push_back(Zvec[i1]);
///
// 	cout << " Block : " << iblock 
// 	     << " - Z("<< i1 <<")= " << Zvec[i1] << endl;
///
	i1++;
      }
    }

    delete[] Zvec;

  }else{

    cout << "Regular block" << endl;

    // Diagonalization works on the STL vector directly!

    zheev_(&JOBZ,&UPLO,&Nst,&A1[0],&LDA,eigv,work,&LWORK,rwork,&INFO);
    
    int i1=0;
    for (int ii=0;ii<Nst;ii++){
      eigvalues.push_back(eigv[ii]);
      for (int jj=0;jj<Nst;jj++){
	eigvectors.push_back(A1[i1]);
	i1++;
      }
    }

  }
  //else if UpTriang

  //cout << "Diag Block " << iblock << " Size: " << Nst << endl;
  //PrintMatBlockQNumbers(iMatBlock);
  //cout << "Lowest States :" << endl;
  //for (int ii=0;ii<( (20 < Nst)?20:Nst );ii++)
  //cout << eigv[ii] << endl;

  delete[] eigv;
  delete[] rwork;
  delete[] work;

  A1.clear();

}

//////////////




//////////////
//
// New April 08- But not going to use it...
void CNRGmatrix::PutInRegularForm(int iblock){

  // Get block info
  int iMatBlock=FindMatBlock(iblock,iblock);

  if (iMatBlock<0)
    {
      cout << "PutInDiagForm: Cannot find block " 
	   << iblock << endl;
      return;
    }
  int NstSq=GetMatBlockSize(iMatBlock);
  int Nst=GetBlockSize(iblock);
  int iblBeg=GetMatBlockLimit(iMatBlock,0);

  // Set "unscrambled" array

  double *BlArray = new double [NstSq];
  int i1,j1,ipos;

  for (int ii=0;ii<NstSq;ii++)
    {
      // Get corresponding i,j (new)
      r2ij(Nst,ii,i1,j1);
      // Get Pmn
      ipos=Pmn(Nst,i1,j1); //(old)
      BlArray[ii]=MatEl[iblBeg+ipos];
    }

  // Unscramble
   for (int ii=0;ii<NstSq;ii++)
     MatEl[iblBeg+ii]=BlArray[ii];

  delete[] BlArray;

}

//////////////
//
// New May 08- Save matrix elements in binary file

void CNRGmatrix::SaveBinary(char filename[]){

  //
  // Saves actually: < c+ + c > (Hermitian)
  //

  FILE *pFile;
  int i,j;
  double auxEl=0.0;

  pFile=fopen(filename,"wb");
  for (i=0;i<Nstates();i++)
    {
      for (j=0;j<Nstates();j++)
	{
	  auxEl=GetMatEl(i,j); 
	  fwrite (&auxEl, sizeof(double), 1 ,pFile); 
	}
    }
  fclose(pFile);



}

//////////////
//

void CNRGmatrix::SaveInOldFormat(){

  char arqname[32];
  char CNsites[8];
  
  sprintf(CNsites,"%d",Nshell);

  strcpy(arqname,"DataN");
  strcat(arqname,CNsites);
  strcat(arqname,".bin");

  cout << " Matrix file is " << arqname << endl;

  SaveBinary(arqname);


}


//////////////


/////////////
/// Save object (Mar. 09)
/////////////


void CNRGmatrix::SaveBin(char arqname[]){

  CNRGarray::SaveBin(arqname);

  vector<double>::iterator dit;
  vector<int>::iterator iit;

  vector<complex<double> >::iterator cit;

  ofstream OutFile(arqname, ios::out | ios::binary | ios::app);
  double daux;
  int iaux;
  complex<double> caux;

  if (!OutFile){cout << "SaveBin: Cannot save data in " << arqname << endl; return;}

  // Save MatEl
  iaux=MatEl.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(dit=MatEl.begin();dit<MatEl.end();dit++){
    daux=*dit;
    OutFile.write(reinterpret_cast<char*>(&daux), sizeof(double));
  }
  // Save MatElCplx
  iaux=MatElCplx.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(cit=MatElCplx.begin();cit<MatElCplx.end();cit++){
    caux=*cit;
    OutFile.write(reinterpret_cast<char*>(&caux), sizeof(complex<double>));
  }
  // Save MatBlockMap
  iaux=MatBlockMap.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=MatBlockMap.begin();iit<MatBlockMap.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }
  // Save MatBlockBegEnd
  iaux=MatBlockBegEnd.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=MatBlockBegEnd.begin();iit<MatBlockBegEnd.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }
  // Save NeedOld
  OutFile.write((char*)&NeedOld, sizeof(bool));
  // Save UpperTriangular
  OutFile.write((char*)&UpperTriangular, sizeof(bool));
  // Save IsComplex
  OutFile.write((char*)&IsComplex, sizeof(bool));

  OutFile.close();

  if (!OutFile.good()){cout << "SaveBin: error saving" << endl;}


}
///

void CNRGmatrix::ReadBin(char arqname[]){

  // ifstream
  ifstream InFile;

  InFile.open(arqname, ios::in | ios::binary);
  //InFile.seekg(0,ios::beg);

  if (!InFile){cout << "ReadBin: Cannot read data in " << arqname << endl; return;}

  CNRGmatrix::ReadBinInStream(InFile);

  InFile.close();

  if (!InFile.good()){cout << "ReadBin: overall error reading" << endl;}


}

/////
void CNRGmatrix::ReadBinInStream(ifstream &InFile){


  int isize,iaux, iaux2;
  double daux;
  complex<double> caux;

  CNRGarray::ReadBinInStream(InFile);
  
  // Read QNumbers_stcf
  MatEl.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&daux), sizeof(double));
    MatEl.push_back(daux);
  }
  // Read MatElCplx
  MatElCplx.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&caux), sizeof(complex<double>));
    MatElCplx.push_back(caux);
  }
  // Read MatBlockMap
  MatBlockMap.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    MatBlockMap.push_back(iaux2);
  }
  // Read MatBlockBegEnd
  MatBlockBegEnd.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    MatBlockBegEnd.push_back(iaux2);
  }
  // Read NeedOld
  InFile.read((char*)&NeedOld, sizeof(bool));
  // Read UpperTriangular
  InFile.read((char*)&UpperTriangular, sizeof(bool));
  // Read IsComplex
  InFile.read((char*)&IsComplex, sizeof(bool));


}
/////////////
void CNRGmatrix::GetConnectingBlocks(int ibl_base,
				     vector < vector<int> > &Connect_BlockMatBlock, 
				     bool ToBase){

  // Given a block ibl_base, and stores the connecting blocks
  // "From" (<ibl_base|Op|ibl>) or "To" (<ibl|Op|ibl_base>)
  // and respective BlockMat in the two-column array
  // Connect_BlockMatBlock

  Connect_BlockMatBlock.clear();
  int ifound=0;
  int iblFrom=ibl_base;
  int iblTo=0;
  for (int iblOp=0;iblOp<NumBlocks();iblOp++){
    iblFrom=ibl_base;
    iblTo=iblOp;
    if (ToBase){
      iblFrom=iblOp;
      iblTo=ibl_base;
    }
//     cout << "NRGmatrix::GetConnectingBlocks: ib_base = " << ibl_base
//  	 << " ibOp = " << iblOp
//  	 << " iblFrom = " << iblFrom
//  	 << " iblTo = " << iblTo 
//     	 << " -> <From|Op|To> ? " << FindMatBlock(iblFrom,iblTo)
//  	 << endl;

    if (FindMatBlock(iblFrom,iblTo)>=0){
      Connect_BlockMatBlock.push_back( vector<int>() );
      Connect_BlockMatBlock[ifound].push_back(iblOp);
      Connect_BlockMatBlock[ifound].push_back(FindMatBlock(iblFrom,iblTo));
      ifound++;
    }
  }
  // end loop in blocks

}
//////////////


void CNRGmatrix::GetBlocksFromMatBlock(int iMatBlock, int &ibl1, int &ibl2){


  vector<int>::iterator iit=MatBlockMap.begin()+2*iMatBlock;

  ibl1=(*iit);
  ibl2=(*(iit+1));


}



//////////////
boost::numeric::ublas::matrix<double> CNRGmatrix::MatBlock2BLAS(int iMatBlock){

  int ibl1=0;
  int ibl2=0;

  GetBlocksFromMatBlock(iMatBlock,ibl1,ibl2);

 
  return(MatBlock2BLAS(ibl1,ibl2));


}
//////////////////
boost::numeric::ublas::matrix<double> CNRGmatrix::MatBlock2BLAS(int ibl1, int ibl2){


  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);
 
  int pos=GetMatBlockLimit(ibl1,ibl2,0);

  boost::numeric::ublas::matrix<double> dMatAux(Nst1,Nst2);

  for (int ii=0;ii<Nst1;ii++){
    for (int jj=0;jj<Nst2;jj++){
      dMatAux(ii,jj)=GetBlockMatEl(ibl1,ibl2,ii,jj);
    }
  }

  return(dMatAux);

}

/// complex
boost::numeric::ublas::matrix<complex<double> > CNRGmatrix::cMatBlock2BLAS(int ibl1, int ibl2){


  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);
 
  int pos=GetMatBlockLimit(ibl1,ibl2,0);

  boost::numeric::ublas::matrix<complex<double> > cMatAux(Nst1,Nst2);

  for (int ii=0;ii<Nst1;ii++){
    for (int jj=0;jj<Nst2;jj++){
      cMatAux(ii,jj)=cGetBlockMatEl(ibl1,ibl2,ii,jj);
    }
  }

  return(cMatAux);

}



//////////////
boost::numeric::ublas::matrix<double> CNRGmatrix::MatBlock2BLAS(int iMatBlock,bool kpi, bool kpj){

  int ibl1=0;
  int ibl2=0;

  GetBlocksFromMatBlock(iMatBlock,ibl1,ibl2);

 
  return(MatBlock2BLAS(ibl1,ibl2,kpi,kpj));


}
//////////////////
boost::numeric::ublas::matrix<double> CNRGmatrix::MatBlock2BLAS(int ibl1, int ibl2,bool kp1, bool kp2){


  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);
 
  int Nstkp1=CNRGarray::GetBlockSize(ibl1,kp1);
  int Nstkp2=CNRGarray::GetBlockSize(ibl2,kp2);
  
  //int pos=GetMatBlockLimit(ibl1,ibl2,0);

  boost::numeric::ublas::matrix<double> dMatAux(Nstkp1,Nstkp2);

  int istkp=0;
  for (int ii=0;ii<Nst1;ii++){
    int jstkp=0;
    for (int jj=0;jj<Nst2;jj++){
      if ((CheckKept(ibl1,ii,kp1))&&(CheckKept(ibl2,jj,kp2))){
	dMatAux(istkp,jstkp)=GetBlockMatEl(ibl1,ibl2,ii,jj);
	jstkp++;
      }
    }
    //end jj
    if (CheckKept(ibl1,ii,kp1)) istkp++; 
  }
  //end ii

  return(dMatAux);

}
//////////////////
boost::numeric::ublas::matrix<complex<double> > CNRGmatrix::cMatBlock2BLAS(int ibl1, int ibl2,bool kp1, bool kp2){


  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);
 
  int Nstkp1=CNRGarray::GetBlockSize(ibl1,kp1);
  int Nstkp2=CNRGarray::GetBlockSize(ibl2,kp2);
  
  //int pos=GetMatBlockLimit(ibl1,ibl2,0);

  boost::numeric::ublas::matrix<complex<double> > cMatAux(Nstkp1,Nstkp2);

  int istkp=0;
  for (int ii=0;ii<Nst1;ii++){
    int jstkp=0;
    for (int jj=0;jj<Nst2;jj++){
      if ((CheckKept(ibl1,ii,kp1))&&(CheckKept(ibl2,jj,kp2))){
	cMatAux(istkp,jstkp)=cGetBlockMatEl(ibl1,ibl2,ii,jj);
	jstkp++;
      }
    }
    //end jj
    if (CheckKept(ibl1,ii,kp1)) istkp++; 
  }
  //end ii

  return(cMatAux);

}
///////////





//////////////
boost::numeric::ublas::matrix<double> CNRGmatrix::RotateBlock2BLAS(CNRGbasisarray* pAcut, int iMatBlock){

  int ibl1=0;
  int ibl2=0;

  GetBlocksFromMatBlock(iMatBlock,ibl1,ibl2);


  return(RotateBlock2BLAS(pAcut,ibl1,ibl2));

}

boost::numeric::ublas::matrix<double> CNRGmatrix::RotateBlock2BLAS(CNRGbasisarray* pAcut, int iMatBlock, bool kp1, bool kp2){

  int ibl1=0;
  int ibl2=0;

  GetBlocksFromMatBlock(iMatBlock,ibl1,ibl2);


  return(RotateBlock2BLAS(pAcut,ibl1,ibl2));

}

//////////////
boost::numeric::ublas::matrix<double> CNRGmatrix::RotateBlock2BLAS(CNRGbasisarray* pAcut, int ibl1, int ibl2){


  // This routine takes an op in the Aeig basis (Nst_bl1,Nst_bl2) and rotates it back
  // to the "before cut" basis (NstBC_bl1,NstBC_bl2)


  // Boost matrices
  
  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);

  int Nst_w1=pAcut->GetBlockSize(ibl1);
  int Nst_w2=pAcut->GetBlockSize(ibl2);


  if ( (Nst1!=Nst_w1)||(Nst2!=Nst_w2) ){
    cout << "RotateBlock: Error: Block sizes do not match." << endl;
    return(boost::numeric::ublas::matrix<double> ());
  }

  int Nst_basis1=pAcut->GetBlockSizeBC(ibl1);
  int Nst_basis2=pAcut->GetBlockSizeBC(ibl2);


  boost::numeric::ublas::matrix<double> Z1bl(Nst_w1,Nst_basis1);
  boost::numeric::ublas::matrix<double> Z2bl(Nst_w2,Nst_basis2);
  boost::numeric::ublas::matrix<double> MatInitial(Nst_w1,Nst_w2);
  boost::numeric::ublas::matrix<double> MatFinal(Nst_basis1,Nst_basis2);

  Z1bl=pAcut->EigVecCut2BLAS(ibl1);
  Z2bl=pAcut->EigVecCut2BLAS(ibl2);

  MatInitial=MatBlock2BLAS(ibl1,ibl2);

  noalias(MatFinal)=prod ( trans(Z1bl), 
	    boost::numeric::ublas::matrix<double>(prod(MatInitial,Z2bl)) );


  return(MatFinal);

}
/// complex
boost::numeric::ublas::matrix<complex<double> > CNRGmatrix::cRotateBlock2BLAS(CNRGbasisarray* pAcut, int ibl1, int ibl2){

  // This routine takes an op in the Aeig basis (Nst_bl1,Nst_bl2) and rotates it back
  // to the "before cut" basis (NstBC_bl1,NstBC_bl2)
  // Boost matrices
  
  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);

  int Nst_w1=pAcut->GetBlockSize(ibl1);
  int Nst_w2=pAcut->GetBlockSize(ibl2);


  if ( (Nst1!=Nst_w1)||(Nst2!=Nst_w2) ){
    cout << "RotateBlock: Error: Block sizes do not match." << endl;
    return(boost::numeric::ublas::matrix<complex<double> > ());
  }

  int Nst_basis1=pAcut->GetBlockSizeBC(ibl1);
  int Nst_basis2=pAcut->GetBlockSizeBC(ibl2);


  boost::numeric::ublas::matrix<complex<double> > Z1bl(Nst_w1,Nst_basis1);
  boost::numeric::ublas::matrix<complex<double> > Z2bl(Nst_w2,Nst_basis2);
  boost::numeric::ublas::matrix<complex<double> > MatInitial(Nst_w1,Nst_w2);
  boost::numeric::ublas::matrix<complex<double> > MatFinal(Nst_basis1,Nst_basis2);

  Z1bl=pAcut->cEigVecCut2BLAS(ibl1);
  Z2bl=pAcut->cEigVecCut2BLAS(ibl2);

  MatInitial=cMatBlock2BLAS(ibl1,ibl2);

//   noalias(MatFinal)=prod ( trans(Z1bl), 
// 	    boost::numeric::ublas::matrix<complex<double> >(prod(MatInitial,Z2bl)) );
// Ops, should be:

  noalias(MatFinal)=prod ( herm(Z1bl), 
	    boost::numeric::ublas::matrix<complex<double> >(prod(MatInitial,Z2bl)) );


  return(MatFinal);

}


////////

//////////////

// Oct 2011
boost::numeric::ublas::matrix<double> CNRGmatrix::RotateBlock2BLAS_NoCut(CNRGarray* pAeig, int ibl1, int ibl2, bool forward){

  // ONLY WORKS IF MATRIX IS IN SYNC WITH Aeig and NO CUTTING HAS BEEN DONE!!!

  // forward=0
  // Takes an op in the "Aeig" basis (Nst_bl1,Nst_bl2) and rotates it backward
  // to the "Abasis" basis.


  // forward=1
  // Takes an op in the "Abasis" basis  and rotates it forward
  // to the "Aeig" basis. 

  // Boost matrices
  
  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);

  int Nst_w1=pAeig->GetBlockSize(ibl1);
  int Nst_w2=pAeig->GetBlockSize(ibl2);


  if ( (Nst1!=Nst_w1)||(Nst2!=Nst_w2) ){
    cout << "RotateForwardBackwardBlock: Error: Block sizes do not match." << endl;
    return(boost::numeric::ublas::matrix<double> ());
  }


  boost::numeric::ublas::matrix<double> Z1bl(Nst_w1,Nst_w1);
  boost::numeric::ublas::matrix<double> Z2bl(Nst_w2,Nst_w2);
  boost::numeric::ublas::matrix<double> MatInitial(Nst_w1,Nst_w2);
  boost::numeric::ublas::matrix<double> MatFinal(Nst_w1,Nst_w2);

  Z1bl=pAeig->EigVec2BLAS(ibl1);
  Z2bl=pAeig->EigVec2BLAS(ibl2);

  MatInitial=MatBlock2BLAS(ibl1,ibl2);


  if (forward)
    noalias(MatFinal)=prod ( Z1bl, 
			   boost::numeric::ublas::matrix<double>(prod(MatInitial,trans(Z2bl))) );
  else
    noalias(MatFinal)=prod ( trans(Z1bl), 
			     boost::numeric::ublas::matrix<double>(prod(MatInitial,Z2bl)) );


  return(MatFinal);
}

//////////////

// Nov 2013
boost::numeric::ublas::matrix<complex<double> > CNRGmatrix::cRotateBlock2BLAS_NoCut(CNRGarray* pAeig, int ibl1, int ibl2, bool forward){

  // ONLY WORKS IF MATRIX IS IN SYNC WITH Aeig and NO CUTTING HAS BEEN DONE!!!

  // forward=0
  // Takes an op in the "Aeig" basis (Nst_bl1,Nst_bl2) and rotates it backward
  // to the "Abasis" basis.


  // forward=1
  // Takes an op in the "Abasis" basis  and rotates it forward
  // to the "Aeig" basis. 

  // Boost matrices
  
  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);

  int Nst_w1=pAeig->GetBlockSize(ibl1);
  int Nst_w2=pAeig->GetBlockSize(ibl2);


  if ( (Nst1!=Nst_w1)||(Nst2!=Nst_w2) ){
    cout << "RotateForwardBackwardBlock: Error: Block sizes do not match." << endl;
    return(boost::numeric::ublas::matrix<complex<double> > ());
  }


  boost::numeric::ublas::matrix<complex<double> > Z1bl(Nst_w1,Nst_w1);
  boost::numeric::ublas::matrix<complex<double> > Z2bl(Nst_w2,Nst_w2);
  boost::numeric::ublas::matrix<complex<double> > MatInitial(Nst_w1,Nst_w2);
  boost::numeric::ublas::matrix<complex<double> > MatFinal(Nst_w1,Nst_w2);

  Z1bl=pAeig->cEigVec2BLAS(ibl1);
  Z2bl=pAeig->cEigVec2BLAS(ibl2);

  MatInitial=cMatBlock2BLAS(ibl1,ibl2);

  //cout << "cRotateBlock2BLAS_NoCut: I'm here III." << endl; 


  if (forward)
//     noalias(MatFinal)=prod ( Z1bl, 
// 			   boost::numeric::ublas::matrix<complex<double> >(prod(MatInitial,trans(Z2bl))) );
    noalias(MatFinal)=prod ( Z1bl, 
			   boost::numeric::ublas::matrix<complex<double> >(prod(MatInitial,herm(Z2bl))) );

  else
//     noalias(MatFinal)=prod ( trans(Z1bl), 
// 			     boost::numeric::ublas::matrix<complex<double> >(prod(MatInitial,Z2bl)) );
    noalias(MatFinal)=prod ( herm(Z1bl), 
			     boost::numeric::ublas::matrix<complex<double> >(prod(MatInitial,Z2bl)) );

  return(MatFinal);

}

////////

////////
boost::numeric::ublas::matrix<double> CNRGmatrix::RotateBlock2BLAS(CNRGbasisarray* pAcut, int ibl1, int ibl2, bool kp1, bool kp2){


  // Boost matrices
  
  int Nst1=CNRGarray::GetBlockSize(ibl1);
  int Nst2=CNRGarray::GetBlockSize(ibl2);

  int Nstkp1=CNRGarray::GetBlockSize(ibl1,kp1);
  int Nstkp2=CNRGarray::GetBlockSize(ibl2,kp2);


  int Nst_w1=pAcut->GetBlockSize(ibl1);
  int Nst_w2=pAcut->GetBlockSize(ibl2);


  if ( (Nst1!=Nst_w1)||(Nst2!=Nst_w2) ){
    cout << "RotateBlock: Error: Block sizes do not match." << endl;
    return(boost::numeric::ublas::matrix<double> ());
  }

  int Nst_basis1=pAcut->GetBlockSizeBC(ibl1);
  int Nst_basis2=pAcut->GetBlockSizeBC(ibl2);

  if ( (Nst1!=Nst_basis1)||(Nst2!=Nst_basis2) ){
    cout << "RotateBlock kept: Error: states are cut! (size in basis does not match)" << endl;
    return(boost::numeric::ublas::matrix<double> ());
  }


  boost::numeric::ublas::matrix<double> Z1bl(Nstkp1,Nstkp1);
  boost::numeric::ublas::matrix<double> Z2bl(Nstkp2,Nstkp2);
  boost::numeric::ublas::matrix<double> MatInitial(Nstkp1,Nstkp2);
  boost::numeric::ublas::matrix<double> MatFinal(Nstkp1,Nstkp2);

  Z1bl=pAcut->EigVecCut2BLAS(ibl1,kp1);
  Z2bl=pAcut->EigVecCut2BLAS(ibl2,kp2);

  MatInitial=MatBlock2BLAS(ibl1,ibl2,kp1,kp2);

  noalias(MatFinal)=prod ( trans(Z1bl), 
	    boost::numeric::ublas::matrix<double>(prod(MatInitial,Z2bl)) );


  return(MatFinal);

}

//////////////////////////

void CNRGmatrix::SetDiagonalMatrix(double El){
  // Set matrix as the identity matrix (UpperTriangular==false only)
  if (UpperTriangular){
    cout << " SetIdentityMatrix for UpperTriangular not set up " << endl;
    return;
  }

  MatBlockMap.clear();
  MatBlockBegEnd.clear();
  MatEl.clear();
  MatElCplx.clear();
  
  complex<double> CEl (1.0,0.0);
  CEl=OneC*El; // complex version of El
  
  int i1=0;
  for (int ibl=0; ibl<NumBlocks();ibl++){
    // Diagonal matrix
    MatBlockMap.push_back(ibl);
    MatBlockMap.push_back(ibl);
    MatBlockBegEnd.push_back(i1);

    int blsize=GetBlockSize(ibl);

    for (int istbl=0;istbl<blsize;istbl++){
      for (int jstbl=0;jstbl<blsize;jstbl++){
	if (istbl==jstbl)
	  if (IsComplex) MatElCplx.push_back(CEl);
	  else MatEl.push_back(El);
	else
	  if (IsComplex) MatElCplx.push_back(ZeroC);
	  else MatEl.push_back(0.0);
	i1++;
      }
    }
    // end loop in block states
    MatBlockBegEnd.push_back(i1-1);
  }
  // end loop in block

}
//////////////////////////

void CNRGmatrix::SetIdentityMatrix(){

  SetDiagonalMatrix(1.0);

}
//////////////////////////

bool CNRGmatrix::SelfCheckForMatEl(int iblock1, int iblock2){

  CNRGbasisarray Aempty;

  ReverseSyncNRGarray(Aempty);

  return(CheckForMatEl(&Aempty,iblock1,iblock2));

}


//////////////////////////

void CNRGmatrix::SetZeroMatrix(){

  // Sets the matrix based on CheckForMatEl (Matrix elements are zero)

  MatBlockMap.clear();
  MatBlockBegEnd.clear();
  MatEl.clear();
  MatElCplx.clear();
  
  complex<double> C1 (0.0,0.0);

  int i1=0;
  for (int ibl1=0;ibl1<NumBlocks();ibl1++){
    for (int ibl2=0;ibl2<NumBlocks();ibl2++){
      if (SelfCheckForMatEl(ibl1,ibl2)){
	// Diagonal matrix
	MatBlockMap.push_back(ibl1);
	MatBlockMap.push_back(ibl2);
	MatBlockBegEnd.push_back(i1);

	int bl1size=GetBlockSize(ibl1);
	int bl2size=GetBlockSize(ibl2);

	for (int istbl=0;istbl<bl1size;istbl++){
	  for (int jstbl=0;jstbl<bl2size;jstbl++){
	    if (IsComplex)
	      MatElCplx.push_back(C1);
	    else
	      MatEl.push_back(0.0);
	    i1++;
	  }
	}
	// end loop in block states
	MatBlockBegEnd.push_back(i1-1);
      }
      // end if check for mat
    }
  }
  // end loops in blocks
  
}
//////////////////////////


//////////////////////////

void CNRGmatrix::SetMatrix(CNRGbasisarray* pAbasis,
			   CNRGbasisarray* pSingleSite){
  // Just to make sure...
  if (!ChkSync(pAbasis)){SyncNRGarray(*pAbasis);}

  MatEl.clear();
  MatElCplx.clear();
  MatBlockMap.clear();
  MatBlockBegEnd.clear();

  double aux=0.0;
  complex<double> caux=ZeroC;
  int icounter=0;
  // Routine Set matrix!!
  for (int ibl=0; ibl<pAbasis->NumBlocks();ibl++){
    for (int jbl=0; jbl<pAbasis->NumBlocks();jbl++){
      if ( CheckForMatEl(pAbasis,ibl,jbl) ){
	MatBlockMap.push_back(ibl);
	MatBlockMap.push_back(jbl);
	MatBlockBegEnd.push_back(icounter);
	for (int ist=pAbasis->GetBlockLimit(ibl,0);ist<=pAbasis->GetBlockLimit(ibl,1);ist++) {
	  for (int jst=pAbasis->GetBlockLimit(jbl,0);jst<=pAbasis->GetBlockLimit(jbl,1);jst++) {
	    if (!IsComplex){
	      aux=CalcMatEl(pAbasis,pSingleSite,ist,jst);
	      MatEl.push_back(aux);
	    }else{
	      caux=CalcMatElCplx(pAbasis,pSingleSite,ist,jst);
	      MatElCplx.push_back(caux);
	    }
	    icounter++;
	  }
	}
	// end loop in states
	MatBlockBegEnd.push_back(icounter-1);
      }
      // end if check
    }
  }
  // end loop in ibl,jbl


}
//////////////////////////
