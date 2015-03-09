#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>


#include <vector>
#include <cmath>
using namespace std;

#include "NRGclasses.hpp"
//#include "NRGfunctions.hpp"


void RotateMatrix(CNRGmatrix* pMat, CNRGbasisarray* pAcut,
		  CNRGmatrix* pMatRot){

  pMatRot->ClearAll();

  //pMatRot->SyncNRGarray(*pAcut);

  pMatRot->MatBlockMap=pMat->MatBlockMap;
  pMatRot->UpperTriangular=pMat->UpperTriangular;
  pMatRot->IsComplex=pMat->IsComplex;


  // Same Block strucuture as pAcut (THIS IS CRUCIAL!)
  pMatRot->NQNumbers=pAcut->NQNumbers;
  pMatRot->QNumbers=pAcut->QNumbers;
  if (pAcut->BlockBegEndBC.size()!=0)
    pMatRot->BlockBegEnd=pAcut->BlockBegEndBC;
  else
    pMatRot->BlockBegEnd=pAcut->BlockBegEnd;



  // Need to Check Sync!!
  bool ChkSync=pMat->ChkSync(pAcut);

  if (!ChkSync){
    cout << "RotateMatrix: pMat and pAcut not in sync" << endl;
    return;
  }

  // Writing the blocks of pMat to the basis in pAcut

  int i1=0;
  for (int iMatbl=0; iMatbl<pMat->NumMatBlocks(); iMatbl++){

    int ibl1=0;
    int ibl2=0;

    pMatRot->MatBlockBegEnd.push_back(i1);

    pMat->GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);

    
 
    int Nst_basis1=pAcut->GetBlockSizeBC(ibl1); // Should work: 
    int Nst_basis2=pAcut->GetBlockSizeBC(ibl2); // ibl1, ibl2 in the current 
                                                // block structure

//      cout << " Rotating block " << iMatbl << " of " << pMat->NumMatBlocks()-1
//  	 << " Final Size: " << Nst_basis1 << " x " << Nst_basis2 << endl; 


    // Boost matrices
    boost::numeric::ublas::matrix<double> MatFinal(Nst_basis1,Nst_basis2);
    boost::numeric::ublas::matrix<complex<double> > cMatFinal(Nst_basis1,Nst_basis2);


    if (pMat->IsComplex){
      //cout << "pMat is complex " << endl;
       cMatFinal=pMat->cRotateBlock2BLAS(pAcut,ibl1,ibl2);
    }
    else
      MatFinal=pMat->RotateBlock2BLAS(pAcut,ibl1,ibl2);

    int j0=0;
    for (int ii=0;ii<Nst_basis1;ii++){
      if ( (pMatRot->UpperTriangular)&&(ibl1==ibl2) )
	j0=ii;
      else 
	j0=0;
      for (int jj=j0;jj<Nst_basis2;jj++){
	if (pMatRot->IsComplex)
	  pMatRot->MatElCplx.push_back(cMatFinal(ii,jj));
	else
	  pMatRot->MatEl.push_back(MatFinal(ii,jj));
	i1++;
      }
    }
    // end loop in MatFinal
    pMatRot->MatBlockBegEnd.push_back(i1-1);
  }
  // end loop in MatBlocks
}
// end subroutine


// Rotate in the UNCUT eigenvector
void RotateMatrix_NoCut(CNRGmatrix* pMat, CNRGarray* pAeig,
		  CNRGmatrix* pMatRot, bool forward){

  pMatRot->ClearAll();

  //pMatRot->SyncNRGarray(*pAeig);

  pMatRot->MatBlockMap=pMat->MatBlockMap;
  pMatRot->UpperTriangular=pMat->UpperTriangular;
  pMatRot->IsComplex=pMat->IsComplex;

  // Same Block strucuture as pAeig (THIS IS CRUCIAL!)
  pMatRot->NQNumbers=pAeig->NQNumbers;
  pMatRot->QNumbers=pAeig->QNumbers;
  pMatRot->BlockBegEnd=pAeig->BlockBegEnd;



  // Need to Check Sync!!
  bool ChkSync=pMat->ChkSync(pAeig);

  if (!ChkSync){
    cout << "RotateMatrix: pMat and pAeig not in sync" << endl;
    return;
  }

  // Writing the blocks of pMat to the basis in pAeig

//   CNRGbasisarray Acut;
//   Acut.FalseCut(pAeig);
//   Acut.PrintAll();

  int i1=0;
  for (int iMatbl=0; iMatbl<pMat->NumMatBlocks(); iMatbl++){

    int ibl1=0;
    int ibl2=0;

    pMatRot->MatBlockBegEnd.push_back(i1);

    pMat->GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);

    
    int Nst_basis1=pAeig->GetBlockSize(ibl1); // Should work: 
    int Nst_basis2=pAeig->GetBlockSize(ibl2); // ibl1, ibl2 in the current 
                                                // block structure

//      cout << " Rotating block " << iMatbl << " of " << pMat->NumMatBlocks()-1
//  	 << " Final Size: " << Nst_basis1 << " x " << Nst_basis2 << endl; 


    // Boost matrices
    boost::numeric::ublas::matrix<double> MatFinal(Nst_basis1,Nst_basis2);
    boost::numeric::ublas::matrix<complex<double> > cMatFinal(Nst_basis1,Nst_basis2);

    //MatFinal=pMat->RotateBlock2BLAS(&Acut,ibl1,ibl2);

    if (pMat->IsComplex){
      cMatFinal=pMat->cRotateBlock2BLAS_NoCut(pAeig,ibl1,ibl2,forward);
    }
    else
      MatFinal=pMat->RotateBlock2BLAS_NoCut(pAeig,ibl1,ibl2,forward);

//     if ( (ibl1==3)&&(ibl2==6) ){
    //cout << " RotMat(ibl1="<<ibl1<<",ibl2="<<ibl2<<")="<< MatFinal << endl;
    //cout << " RotMat(ibl1="<<ibl1<<",ibl2="<<ibl2<<")="<< cMatFinal << endl;
//     }

    int j0=0;
    for (int ii=0;ii<Nst_basis1;ii++){
      if ( (pMatRot->UpperTriangular)&&(ibl1==ibl2) )
	j0=ii;
      else 
	j0=0;
      for (int jj=j0;jj<Nst_basis2;jj++){
	if (pMatRot->IsComplex){
	  pMatRot->MatElCplx.push_back(cMatFinal(ii,jj));
	}
	else{
	  pMatRot->MatEl.push_back(MatFinal(ii,jj));
	}
	i1++;
      }
    }
    // end loop in MatFinal

    pMatRot->MatBlockBegEnd.push_back(i1-1);
    
  }
  // end loop in MatBlocks


}
// end subroutine


///////////////////////////////////////////////////////
///// Moving all RotateBlock2BLAS to here  ////////////
///////////////////////////////////////////////////////

