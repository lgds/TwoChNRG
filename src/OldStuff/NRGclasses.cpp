
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include <math.h>

#include "NRGclasses.hpp"
using namespace std;



//////////////
double CNRGarray::EniBlockj(int i, int j){

  // returns the energy of the i-th state in block j

  double aux=0.0;
  
  if ( (2*j-2)<BlockBegEnd.size() )
    {
      if ( (BlockBegEnd[2*j-2]+(i-1))<dEn.size() )
	{
	  aux=dEn[BlockBegEnd[2*j-2]+(i-1)];
	}
    }
  
  return(aux);

}
//////////////
// void CNRGarray::SetsizeBlock(){

//   int ii,i1;
//   for (ii=0;ii<BlockBegEnd.size();ii+=2)
//     {
//       sizeBlock.push_back(BlockBegEnd[ii+1]-BlockBegEnd[ii]+1);
//     }
  

// }
//////////////
void CNRGarray::PrintQNumbers(){


  int ii;
  cout << "Qnumbers : " << endl;
  for (ii=0;ii<QNumbers.size();ii++)
    cout << QNumbers[ii] << "  ";
  cout << endl;
}

//////////////
void CNRGarray::PrintAll(){


  int ii;
  void CNRGarray::PrintQNumbers();

  PrintQNumbers();

  cout << "dEn : " << endl;
  for (ii=0;ii<dEn.size();ii++)
    cout << dEn[ii] << " ";
  cout << endl;

  cout << "BlockBegEnd : " << endl;
  for (ii=0;ii<BlockBegEnd.size();ii+=2)
    cout << BlockBegEnd[ii] << " "<< BlockBegEnd[ii+1] << ", ";
  cout << endl;

  cout << "Size Blocks : " << endl;
  for (ii=0;ii<NumBlocks();ii++)
    cout << GetBlockSize(ii) << " ";
  cout << endl;
}
//////////////

int CNRGarray::NumBlocks(){
  
  return((int)(QNumbers.size()/NQNumbers));

}
//////////////
int CNRGarray::Nstates(){

  return(BlockBegEnd[NQNumbers*(NumBlocks()-1)+1]-BlockBegEnd[0]+1);
  
}


//////////////
void CNRGarray::SetDegen(){

  int nst,ii,jj;

  // Needs work
  
//   iDegen.push_back(1);
//   ii=1;
//   while ( ii<dEn.size() )
//     {
//       if (fabs(dEn[ii]-dEn[ii-1])<1.0E-10) nst=iDegen[ii-1]+1;
//       else nst=1;
//       iDegen.push_back(nst);
//       ii++;
//     }

}


void CNRGarray::FilterQNumbers(){

  int equal=0;
  int ii,i1;

  vector<double>::iterator qnums_iter,qnums_iter2,qnums_ii1,qnums_ii2;

  for (qnums_iter=QNumbers.begin(); qnums_iter<QNumbers.end(); 
       qnums_iter+=NQNumbers)
    {

      //      cout << "Iter1 = " << (*qnums_iter) << endl;
      ii=1;
      qnums_iter2=qnums_iter+ii*NQNumbers;
      while (qnums_iter2<QNumbers.end())
	{
	  //  cout << "Iter2 = " << (*qnums_iter2) << endl;
	  equal=1;
	  i1=0;
	  for (qnums_ii2=qnums_iter2;
	       qnums_ii2<qnums_iter2+NQNumbers;qnums_ii2++)
	    {
	      qnums_ii1=qnums_iter+i1;
	      if ( (*qnums_ii2==*qnums_ii1)&&(equal==1) ) equal=1;
	      else equal=0;
	      i1++;
	    }
	  if (equal==1)
	    {
	      QNumbers.erase(qnums_iter2,qnums_iter2+NQNumbers);             
	    }
	  ii++;
	  qnums_iter2=qnums_iter+ii*NQNumbers;
	}
    }

}

//////////////
void CNRGarray::ClearAll(){

  QNumbers.clear();
  dEn.clear();
  dEigVec.clear();
  iDegen.clear(); 

  BlockBegEnd.clear();
  

}

//////////////
double CNRGarray::GetQNumber(int iblock, int whichqn){

  int iq=iblock*NQNumbers+whichqn;
      
  return(QNumbers[iq]);
}

//////////////
int CNRGarray::GetBlockLimit(int iblock, int whichlimit){

  return(BlockBegEnd[2*iblock+whichlimit]);
  
}

//////////////
int CNRGarray::GetBlockSize(int iblock){

  return(BlockBegEnd[2*iblock+1]-BlockBegEnd[2*iblock]+1);
  
}


////////////////////////////////////
// Class CNRGbasisarray functions //
////////////////////////////////////


//////////////
//void CNRGbasisarray::ClearVecsBasis(){
void CNRGbasisarray::ClearAll(){

  CNRGarray::ClearAll();

  iType.clear();
  StCameFrom.clear();
  

}
//////////////
void CNRGbasisarray::PrintAll(){

  CNRGarray::PrintAll();

  cout << "Type : " << endl;
  for (int ii=0;ii<iType.size();ii++)
    cout << iType[ii] << " ";
  cout << endl;

  cout << "StCameFrom : " << endl;
  for (int ii=0;ii<StCameFrom.size();ii++)
    cout << StCameFrom[ii] << " ";
  cout << endl;
  

}
///////////////////

////////////////////////////////////
// Class CNRGmatrix functions //
////////////////////////////////////

// vector<double> CNRGmatrix::GetFullMatrix(int Nst){

//   vector<double> Mat((Nst*Nst),0.0); //All zeros
  
//   int i1;
//   for (int ii=0;ii<vec.size();ii++)
// 	{
// 	  i1=ij2r(Nst,vec[ii].ist,vec[ii].jst);
// 	  if (i1<Mat.size()) Mat[i1]=vec[ii].val;
// 	  i1=ij2r(Nst,vec[ii].jst,vec[ii].ist);
// 	  if (i1<Mat.size()) Mat[i1]=vec[ii].val;
// 	}

//   return(Mat);

// }
//////////////
int CNRGmatrix::GetMatBlockLimit(int iblock1, int iblock2, int whichlimit){

  int ii=0;

  while ( ((iblock1!=MatBlockMap[ii])||(iblock2!=MatBlockMap[ii+1]))
	  &&(ii<MatBlockMap.size()) ) ii+=2;

  return(MatBlockBegEnd[2*ii+whichlimit]);


}

//////////////
int CNRGmatrix::GetMatBlockSize(int iblock1, int iblock2){

  int ii=0;

  while ( ((iblock1!=MatBlockMap[ii])||(iblock2!=MatBlockMap[ii+1]))
	  &&(ii<MatBlockMap.size()) ) ii+=2;

  return(MatBlockBegEnd[2*ii+1]-MatBlockBegEnd[2*ii]+1);

}

void CNRGmatrix::DiagBlock(int iblock,
			   vector<double> &eigvalues, 
			   vector<double> &eigvectors){

  // LAPACK variables
  char JOBZ='N',UPLO='U';
  int Nst;
  int INFO=0,LDA,LWORK;

  // LDA=Nst
  // LWORK=3*Nst-1

  double *A;
  double *eigv;
  double *work;

  // Get Block info

  
  Nst=GetMatBlockSize(iblock,iblock);
//   cout << "Size block                 : " << GetBlockSize(iblock) << endl;
//   cout << "No of independent elements : " << GetMatrixBlockSize(iblock) 
//        << endl;
  
  //Nst=2;
  LDA=Nst;
  LWORK=3*Nst-1;
  eigv = new double [Nst];
  work = new double [LWORK];
  A = new double[Nst*Nst];


  // I'll work on a routine to map upper triangular to r. 
  // For now, let's use this:
  double **Aaux;
  Aaux = new double*[Nst];
  for (int ii=0;ii<Nst;ii++) Aaux[ii]=new double[Nst];


  for (int ii=0;ii<MatEl.size();ii++)
    {
//       // This takes a tremendous amount of CPU time! Must do better!
//       if ( (vec[ii].Bli==iblock)&&(vec[ii].Blj==iblock) )
// 	{
// 	  // Awful construcion. An invitation to segmentation faults.
// 	  Aaux[vec[ii].ibl][vec[ii].jbl]=vec[ii].val;
// 	  Aaux[vec[ii].jbl][vec[ii].ibl]=vec[ii].val;
// 	}
    }
  int i1=0;
  for (int ii=0;ii<Nst;ii++)
    {
      for (int jj=0;jj<Nst;jj++)
	{
	  A[i1]=Aaux[ii][jj];
	  i1++;
	}
    }
//   int i1=0;
//   int r1=0;
//   for (int i=GetBlockLimit(iblock,0);
//            i<GetBlockLimit(iblock,1);i++)
//     {
//       for (int j=i;
//            j<GetBlockLimit(iblock,1);j++)
// 	{
// 	  r1=ij2r(Nstates(),i,j);
//       Aaux[]=vec[r1].val;
// 	}
//       i1++;
//     }

//   A[0]=1.0;
//   A[1]=-sqrt(2.0);
//   A[2]=-sqrt(2.0);
//   A[3]=1.0;

  dsyev_(&JOBZ,&UPLO,&Nst,&A[0],&LDA,eigv,work,&LWORK,&INFO);

  for (int ii=0;ii<Nst;ii++)
    {
      eigvalues.push_back(eigv[ii]);
      eigvectors.push_back(eigv[ii]);
      cout << " Block  : " << iblock << endl;
      std::cout << "E = " << eigvalues[ii] << std::endl;
    }


  delete[] eigv;
  delete[] work;
  delete[] A;

  for (int ii=0;ii<Nst;ii++) delete[] Aaux[ii];
  delete[] Aaux;


}
