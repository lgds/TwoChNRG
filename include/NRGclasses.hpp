
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
////using namespace boost::numeric::ublas;

#include <gsl/gsl_integration.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <complex>      // std::complex
using namespace std;

#ifndef _CPLXCONST_
#define _CPLXCONST_

const  complex<double> ZeroC (0.0,0.0);
const  complex<double> OneC (1.0,0.0);
const  complex<double> OneImC (0.0,1.0);

#endif


// Declares the classes
////////////////////////////////////////////////
////////////// Class CNRGarray \\\\\\\\\\\\\\\\\
////////////////////////////////////////////////


#ifndef _CNRGARRAY_
#define _CNRGARRAY_

class CNRGarray { 
/* Declares the class */

public:

  // Default constructor

  //CNRGarray():totalS(false){}
  CNRGarray(){}
  // Constructor with NQNumbers defined
  CNRGarray(int Nqns){
    NQNumbers=Nqns;
  }
  // Default  Destructor
  ~CNRGarray(){}


  // Member variables (regular)
  int Nshell;
  int NQNumbers;

  // Member variables (STL vectors)
  vector<double> QNumbers;


  // Vectors that spam the full basis
  vector<double> dEigVec;
  vector<double> dEn;
  vector<complex<double> > cEigVec;


  vector<int> iDegen; 

  vector<int> BlockBegEnd;

  // April 09: list of N+1 states generated
  // (to be used in DM-NRG)

  vector< vector<int> > ChildStates;

 // August 2011: mark kept states. (spans the full basis)
  vector<bool> Kept;

  // April 10: totalS qnumbers
  bool totalS; // tru -> there are totalS QNs 
  vector<int> Sqnumbers; // position(s) of total S in QNumbers
  // Idea: This can be used in an updated version of BuildBasis...
  //bool ComplexEVec; // tru -> complex eigenvectors


 
  // Member functions
  int NumBlocks();
  int Nstates();
  double EniBlockj(int i, int j);
  void SetDegen();
  void SetE0zero();


  void FilterQNumbers();

  void PrintQNumbers();

  void PrintEn();


  void PrintAll();

  void PrintBlockQNumbers(int iblock);

  void PrintBlockEn(int iblock);

  void PrintBlock(int iblock);


  void ClearAll();


  double GetQNumber(int iblock, int whichqn);

  int GetBlockLimit(int iblock, int whichlimit);
  int GetBlockSize(int iblock);

  int GetBlockSize(int iblock,bool kp);


  int GetBlockFromSt(int ist,int &pos_in_block);

  int GetBlockLimitEigv(int iblock, int whichlimit);

  boost::numeric::ublas::matrix<double> EigVec2BLAS(int iblock);

  // 2013
  boost::numeric::ublas::matrix<complex<double> > cEigVec2BLAS(int iblock);


  double GetEigVecComponent(int ist,int istbasis); // April 09

  double Ecut(int Ncut);

  double GetBlockEmin(int iblock);

  double GetBlockEmax(int iblock);

  int GetBlockEcutPos(int iblock,double Ecut);

  void SetKept(int Ncutoff);
  bool CheckKept(int ist, bool kp);
  bool CheckKept(int iblock, int istbl, bool kp);


  int GetBlockFromQNumbers(double *qnums);

  double GetQNumberFromSt(int ist, int whichqn);

  // May 08
  int GetiGS();
  void SaveQSParameters();

  // Sep 08
  void AddQNumber();

  // Feb 09
  void SaveBin(char arqname[]);
  void ReadBin(char arqname[]);

  void ReadBinInStream(ifstream &InFile);

  // April 09

  double PartitionFuncTeq0();

  double PartitionFunc(double betabar);

  // September 09
  void SetEigVecToOne();

  // Nov 2010
  boost::numeric::ublas::matrix<double> ExpEi2BLAS(int iblock, double betabar);

  // March 2014
  bool CheckComplex();



protected :



};

#endif

////////////////////////////////////////////////
/////////// Class CNRGbasisarray \\\\\\\\\\\\\\\
////////////////////////////////////////////////
//
//  Subclass of CNRGarray
//


#ifndef _CNRGBASISARRAY_
#define _CNRGBASISARRAY_

class CNRGbasisarray: public CNRGarray{

public :

  // Default constructor
  CNRGbasisarray(){}

  CNRGbasisarray(int Nqns){
    NQNumbers=Nqns;
  }
  // Default  Destructor

   // Constructor that copies members from a CNRGarray object
   CNRGbasisarray(CNRGarray A1){
     CNRGarray::Nshell=A1.Nshell;
     CNRGarray::NQNumbers=A1.NQNumbers;
     CNRGarray::QNumbers=A1.QNumbers;
     CNRGarray::dEn=A1.dEn;
     CNRGarray::dEigVec=A1.dEigVec;
     CNRGarray::cEigVec=A1.cEigVec;
     CNRGarray::Kept=A1.Kept;
     CNRGarray::BlockBegEnd=A1.BlockBegEnd;
     CNRGarray::totalS=A1.totalS; // tru -> there are totalS QNs 
     CNRGarray::Sqnumbers=A1.Sqnumbers; // position(s) of total S in QNumbers
   }

  ~CNRGbasisarray(){}


  vector<int> iType;
  vector<int> StCameFrom;

  // moved here



  // NEW
  vector<int> BlockBegEndBC;
  vector<int> BlockBegEndEigVec;

  // Newer: if the symmetry is reduced at a certain point, need to retain 
  // old Qnumbers

  int NQNumbers_stcf;
  vector<double>  StCameFromQNumbers;

  // Member functions

  void ClearAll();
  void PrintAll();
  void PrintBlockBasis(int iblock);

  void PrintBasisAll(){
    for (int ibl=0;ibl<NumBlocks();ibl++) 
      PrintBlockBasis(ibl);
  }


  void SetSCFfromdEn();

  int  GetBlockSizeBC(int iblock);

  int GetBlockFromStBC(int ist,int &pos_in_block);

  void SyncNRGarray(CNRGarray &A1){
     CNRGarray::Nshell=A1.Nshell;
     CNRGarray::NQNumbers=A1.NQNumbers;
     CNRGarray::QNumbers=A1.QNumbers;
     CNRGarray::dEn=A1.dEn;
     CNRGarray::dEigVec=A1.dEigVec;
     CNRGarray::cEigVec=A1.cEigVec;
     CNRGarray::Kept=A1.Kept;
     CNRGarray::BlockBegEnd=A1.BlockBegEnd;
     CNRGarray::totalS=A1.totalS; // tru -> there are totalS QNs 
     CNRGarray::Sqnumbers=A1.Sqnumbers; // position(s) of total S in QNumbers
  }


  // New
  void FalseCut(CNRGarray* pA1){
    CNRGbasisarray::ClearAll();
    CNRGbasisarray::SyncNRGarray(*pA1);
    CNRGbasisarray::SetSCFfromdEn();
    CNRGbasisarray::BlockBegEndBC=pA1->BlockBegEnd;
    CNRGbasisarray::SetBlockBegEndEigVec();
  }


  void SetBlockBegEndEigVec();

  void RemoveBlock(int iblock); // StCameFrom keeps the pre-cut state indexes


  void RemoveState(int ist); // StCameFrom keeps the pre-cut state indexes

  void RemoveStatesFromBlock(int ibl, int ist1,int ist2);

  int GetBlockLimitEigv(int iblock, int whichlimit); // April 09, overloads

  double GetEigVecComponent(int ist,int istbasis); // April 09 overloads

  // Ubuntu's g++ does not like this... 
  boost::numeric::ublas::matrix<double> EigVecCut2BLAS(int iblock); 
  boost::numeric::ublas::matrix<complex<double> > cEigVecCut2BLAS(int iblock); 
  // converts the non-square eigenvector matrix to BLAS format

  boost::numeric::ublas::matrix<double> EigVecCut2BLAS(int iblock, bool kp); 
  // converts the non-square eigenvector matrix to BLAS format (kp states)


  // April 08
  double GetStCameFromQNumberFromSt(int ist, int whichqn);

  // Sets dEn as the energies of the previous states in StCameFrom
  // SetdEnPrevious();


  void CopyBlock(int iblock, int NoCopies); 
  // Create copies of block states and set iDegen  


  // Sep 08
  void SetLastQNumber(vector<double> InputVec);

  // Feb 09

  void SaveBin(char arqname[]);
  void ReadBin(char arqname[]);

  void ReadBinInStream(ifstream &InFile);




};

#endif



////////////////////////////////////////////////
/////////// Class CNRGmatrix     \\\\\\\\\\\\\\\
////////////////////////////////////////////////
//
//  Subclass of CNRGarray
//

#ifndef _CNRGMATRIX_
#define _CNRGMATRIX_

class CNRGmatrix: public CNRGarray{

public :

  // Default constructor (with initialization list)
  CNRGmatrix():UpperTriangular(false),NeedOld(false),SaveMatYN(false),CalcAvg(false),IsComplex(false),WignerEckartL(0.5)
  {}
  //CNRGmatrix(){}

   // Constructor that copies members from a CNRGarray object
   CNRGmatrix(CNRGbasisarray A1){
     UpperTriangular=false;
     SaveMatYN=false;
     CalcAvg=false;
     IsComplex=false;
     CNRGarray::Nshell=A1.Nshell;
     CNRGarray::NQNumbers=A1.NQNumbers;
     CNRGarray::QNumbers=A1.QNumbers;
     CNRGarray::BlockBegEnd=A1.BlockBegEnd;
     CNRGarray::Kept=A1.Kept;
     //     CNRGarray::sizeBlock=A1.sizeBlock;
     iType=A1.iType;
     CNRGarray::totalS=A1.totalS; // tru -> there are totalS QNs 
     CNRGarray::Sqnumbers=A1.Sqnumbers; // position(s) of total S in QNumbers
   } 
   // Default  Destructor
   ~CNRGmatrix(){}

  //////////////////////
  //   Static members //
  //////////////////////

  vector<double> MatEl;
  // new
  //template <typename AnyType>
  vector< complex<double> > MatElCplx;
  // Best solution: look for MatEl everywhere and 
  // add a if (IsComplex)-> use MatElCmplx instead


  vector<int> MatBlockMap;
  vector<int> MatBlockBegEnd;

  // New: same as in CNRGbasisarray
  vector<int> iType;

  // NEW: indicate kept-kept,disc-disc, disc-kept,...
  vector<bool> MatKept;

  // Does this matrix needs the "old" one on updating?
  bool NeedOld;

  // Stores upper triangular part of the blocks only?
  bool UpperTriangular;

  // Save matrix?
  bool SaveMatYN;

  // CalcAvg?
  bool CalcAvg;
  char MatName[20];

  // Is complex?
  bool IsComplex;

  // Wigner-Eckart rank (usually 1/2 but it can be 1)
  // Needed when using SU(2) symmetry.
  double WignerEckartL;

  
  //////////////////////
  // Function members //
  //////////////////////

  void SyncNRGarray(CNRGarray A1){
     CNRGarray::Nshell=A1.Nshell;
     CNRGarray::NQNumbers=A1.NQNumbers;
     CNRGarray::QNumbers=A1.QNumbers;
     CNRGarray::BlockBegEnd=A1.BlockBegEnd;
     CNRGarray::Kept=A1.Kept;
     CNRGarray::totalS=A1.totalS; // tru -> there are totalS QNs 
     CNRGarray::Sqnumbers=A1.Sqnumbers; // position(s) of total S in QNumbers
  }
  /// Same for basis array
  void SyncNRGarray(CNRGbasisarray A1){
     CNRGarray::Nshell=A1.Nshell;
     CNRGarray::NQNumbers=A1.NQNumbers;
     CNRGarray::QNumbers=A1.QNumbers;
     CNRGarray::BlockBegEnd=A1.BlockBegEnd;
     CNRGarray::Kept=A1.Kept;
     CNRGarray::totalS=A1.totalS; // tru -> there are totalS QNs 
     CNRGarray::Sqnumbers=A1.Sqnumbers; // position(s) of total S in QNumbers
     iType=A1.iType;
  }
//   ////
//   void ReverseSyncNRGarray(CNRGarray &A1){
//      A1.Nshell=CNRGarray::Nshell;
//      A1.NQNumbers=CNRGarray::NQNumbers;
//      A1.QNumbers=CNRGarray::QNumbers;
//      A1.BlockBegEnd=CNRGarray::BlockBegEnd;
//   }
  /// Same for basis array
  void ReverseSyncNRGarray(CNRGbasisarray &A1){
     A1.Nshell=CNRGarray::Nshell;
     A1.NQNumbers=CNRGarray::NQNumbers;
     A1.QNumbers=CNRGarray::QNumbers;
     A1.BlockBegEnd=CNRGarray::BlockBegEnd;
     A1.Kept=CNRGarray::Kept;
     A1.totalS=CNRGarray::totalS; // tru -> there are totalS QNs 
     A1.Sqnumbers=CNRGarray::Sqnumbers; // position(s) of total S in QNumbers
     A1.iType=iType;
  }



  bool ChkSync(CNRGarray* pA1){
    if ( (CNRGarray::Nshell==pA1->Nshell)&&
	 (CNRGarray::NQNumbers==pA1->NQNumbers)&&
	 (CNRGarray::QNumbers==pA1->QNumbers)&&
	 (CNRGarray::BlockBegEnd==pA1->BlockBegEnd)&&
	 (CNRGarray::Kept==pA1->Kept)&&
	 (CNRGarray::totalS==pA1->totalS)&&
	 (CNRGarray::Sqnumbers==pA1->Sqnumbers)
	 ){
      return(true);
    }
    else{return(false);}
  }
  /// Same for basis array
  bool ChkSync(CNRGbasisarray* pA1){
    if ( (CNRGarray::Nshell==pA1->Nshell)&&
	 (CNRGarray::NQNumbers==pA1->NQNumbers)&&
	 (CNRGarray::QNumbers==pA1->QNumbers)&&
	 (CNRGarray::BlockBegEnd==pA1->BlockBegEnd)&&
	 (CNRGarray::Kept==pA1->Kept)&&
	 (CNRGarray::totalS==pA1->totalS)&&
	 (CNRGarray::Sqnumbers==pA1->Sqnumbers)&&
	 (iType==pA1->iType) // new (Sep 09, checking)
	 ){
      return(true);
    }
    else{return(false);}
  }

  ////

  void CopyData(CNRGmatrix* pMat1){
    // Array data
    Nshell=pMat1->Nshell;
    NQNumbers=pMat1->NQNumbers;
    QNumbers=pMat1->QNumbers;
    BlockBegEnd=pMat1->BlockBegEnd;
    Kept=pMat1->Kept;
    iType=pMat1->iType;
    totalS=pMat1->totalS;
    Sqnumbers=pMat1->Sqnumbers;
    // Matrix data
    MatEl=pMat1->MatEl;
    MatElCplx=pMat1->MatElCplx;
    MatBlockMap=pMat1->MatBlockMap;
    MatBlockBegEnd=pMat1->MatBlockBegEnd;
    MatKept=pMat1->MatKept;
    UpperTriangular=pMat1->UpperTriangular; 
    IsComplex=pMat1->IsComplex;
    WignerEckartL=pMat1->WignerEckartL;
  }


  int NumMatBlocks();

  int FindMatBlock(int iblock1, int iblock2);

  int GetMatBlockLimit(int iblock1, int iblock2, int whichlimit);
  int GetMatBlockLimit(int imatblock, int whichlimit);

  int GetMatBlockSize(int iblock1, int iblock2);
  int GetMatBlockSize(int imatblock );

  double GetBlockMatEl(int iblock1, int iblock2, int iel, int jel);
  double GetMatEl(int ist, int jst);
  complex<double> cGetBlockMatEl(int iblock1, int iblock2, int iel, int jel);
  complex<double> cGetMatEl(int ist, int jst);



  int GetMatElPosition(int iblock1, int iblock2, int iel, int jel);

  // deprecated
  void PushBlockMatEl(double El, int iblock1, int iblock2, int iel, int jel);
  ///
  void PushMatEl(double El, int ist, int jst);
  void PushMatEl(complex<double> El, int ist, int jst);

  ////// TEMPLATES! Testing. Not easy.
  template<class T> void TemplateTest(T El){

    cout << " This is the element: " << El << endl;
      return;
  } 

// Replace EXISTING matrix element ist,jst by El
  template <class T> void TPushMatEl(T El, int ist, int jst){
    int ibl,jbl;
    
    int iblock1=CNRGarray::GetBlockFromSt(ist,ibl);
    int iblock2=CNRGarray::GetBlockFromSt(jst,jbl);

    int iPos=CNRGmatrix::GetMatElPosition(iblock1,iblock2,ibl,jbl);

    // This does not compile!
//     if (iPos>-1){
//       if (IsComplex) MatElCplx[iPos]=El;
//       else MatEl[iPos]=El;
//     }
    
  }
  //// END TEMPLATE TESTING

  void FilterMap_SetBegEnd();

  //  double* Block2Array(int imatblock);

  void ClearAll();

  void PrintMatBlockQNumbers(int imatblock);

  void PrintMatBlock(int imatblock);
  void PrintMatBlock(int iblock1, int iblock2);

  void PrintAllBlocks();

  void DiagBlock(int iblock,
		 vector<double> &eigvalues, 
		 vector<double> &eigvectors);
  void DiagBlock(int iblock,
		 vector<double> &eigvalues, 
		 vector<complex<double> > &eigvectors);

  


  // Pointers to functions to be defined


  // Checks whether there are non-zero matrix elements 
  bool (*CheckForMatEl)(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);


  // Calculates matrix elements between iblock1 and iblock2
  double (*CalcMatEl)(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
		      int isti, int istj);

  complex<double> (*CalcMatElCplx)(CNRGbasisarray *pAbasis, 
			       CNRGbasisarray *pSingleSite, 
			       int isti, int istj);


  double (*CalcHNMatEl)(vector<double> Params,
		      CNRGbasisarray* pAbasis, 
		      CNRGbasisarray* pSingleSite,
		      CNRGmatrix* MatArray,
		      int ist, int jst);

  complex<double> (*CalcHNMatElCplx)(vector<double> Params,
		      CNRGbasisarray* pAbasis, 
		      CNRGbasisarray* pSingleSite,
		      CNRGmatrix* MatArray,
		      int ist, int jst);


  //
  // THIS one is looking good!!!
  // 


  void DiagHN(vector<double> ParamsHN,
	      CNRGbasisarray* pAbasis, 
	      CNRGbasisarray* pSingleSite,
	      CNRGmatrix* MatArray,
	      CNRGarray* pAeig,bool display=false);

  //
  // This one too!
  // 

  void PutInRegularForm(int iblock);

  //
  // Save in Binary
  //

  void SaveBinary(char arqname[]);

  void SaveInOldFormat();


  // Mar 09

  void SaveBin(char arqname[]);
  void ReadBin(char arqname[]);

  void ReadBinInStream(ifstream &InFile);

  // April 09

  // Given a block ibl_base, and stores the connecting blocks
  // and respective BlockMat in the two-column array
  // Connect_BlockMatBlock
  void GetConnectingBlocks(int ibl_base,
			   vector < vector<int> > &Connect_BlockMatBlock,
			   bool ToBase=false); 
  // May 09


  void GetBlocksFromMatBlock(int iMatBlock, int &ibl1, int &ibl2);

  boost::numeric::ublas::matrix<double> MatBlock2BLAS(int iMatBlock);


  boost::numeric::ublas::matrix<double> MatBlock2BLAS(int ibl1, int ibl2);
  boost::numeric::ublas::matrix<complex<double> > cMatBlock2BLAS(int ibl1, int ibl2);
  
  boost::numeric::ublas::matrix<double> RotateBlock2BLAS(CNRGbasisarray* pAcut, int iMatBlock);

  boost::numeric::ublas::matrix<double> RotateBlock2BLAS(CNRGbasisarray* pAcut, int ibl1, int ibl2);
  boost::numeric::ublas::matrix<complex<double> > cRotateBlock2BLAS(CNRGbasisarray* pAcut, int ibl1, int ibl2);



  boost::numeric::ublas::matrix<double> RotateBlock2BLAS_NoCut(CNRGarray* pAeig, int ibl1, int ibl2, bool forward);
  boost::numeric::ublas::matrix<complex<double> > cRotateBlock2BLAS_NoCut(CNRGarray* pAeig, int ibl1, int ibl2, bool forward);



  // August 2011

  boost::numeric::ublas::matrix<double> MatBlock2BLAS(int iMatBlock, 
						      bool kp1, bool kp2);


  boost::numeric::ublas::matrix<double> MatBlock2BLAS(int ibl1, int ibl2, bool kp1, bool kp2);
  boost::numeric::ublas::matrix<complex<double> > cMatBlock2BLAS(int ibl1, int ibl2, bool kp1, bool kp2);

  
  boost::numeric::ublas::matrix<double> RotateBlock2BLAS(CNRGbasisarray* pAcut, int iMatBlock, bool kp1, bool kp2);

  boost::numeric::ublas::matrix<double> RotateBlock2BLAS(CNRGbasisarray* pAcut, int ibl1, int ibl2, bool kp1, bool kp2);



  // Sep 09
  // Set as Diagonal matrix (diagonal, UpperTriang=false)
  void SetDiagonalMatrix(double El);
  void SetIdentityMatrix();

  // Uses it's own QNumbers and iType to check for mat els. 
  bool SelfCheckForMatEl(int iblock1, int iblock2);
  // Sets Matrix as zero, based on CheckForMatEl (SelfCheckForMatEl)
  void SetZeroMatrix();

  // May 2011
  // IF (*CalcMatEl) and (*CheckForMatEl) are defined, then it sets the matrix
  void SetMatrix(CNRGbasisarray* pAbasis,
		 CNRGbasisarray* pSingleSite);

};

#endif





#ifndef _CNRGCHAIN_
#define _CNRGCHAIN_


class CNRGchain{


public :

  // Default constructor (with initialization list)

  CNRGchain():Lambda(2.5),z_twist(1.0),BandOffset(0.0),HybRef(1.0),HybRefIsSet(false),HybFuncIsSumDeltas(false),Fsq(2.0),DiscScheme(0)
  {}

  CNRGchain(double ThisLambda):z_twist(1.0),BandOffset(0.0),HybRef(1.0),HybRefIsSet(false),HybFuncIsSumDeltas(false),Fsq(2.0),DiscScheme(0){
    Lambda=ThisLambda;
  }

  CNRGchain(double ThisLambda, int Nsitesmax):z_twist(1.0),BandOffset(0.0),HybRef(1.0),HybRefIsSet(false),HybFuncIsSumDeltas(false),Fsq(2.0),DiscScheme(0)
  {
    Lambda=ThisLambda;
    //z_twist=1.0;
    SetChainWilson(Nsitesmax);
  }

  // Default  Destructor
  ~CNRGchain(){}

  // Static members
  //int Nshell;
  //vector<double> eps_chi;

  double Lambda;
  double z_twist;
  double BandOffset;

  vector<double> en;
  vector<double> chin;

  vector<double> *HybFuncParams; // pointer to an STL vector
  // It SHOULD NOT be a POINTER!! Should allocate this WITHIN THE OBJECT!
  // Leave it like that for now...

  //vector< vector<double> > HybFuncVecParams; // STL array

  double HybRef;
  bool HybRefIsSet;
  bool HybFuncIsSumDeltas;
  double Fsq;
  int DiscScheme; // 0- usual; 1- CampoOliveira


  // Function members

  //void SetChain(double Lambda);


  void PrintChain(int Nsites);

  void PrintAll();

  void SetChainWilson(int Nsitesmax);

  void ReadParams(char arqname[],int Nparams);

  void ReadParams(char arqname[]); // Reads the whole thing into the STL vector
                                   // LAST NUMBER: Nlines.

  //void ReadVecParams(char arqname[], int Ncols);

  
  // Hybridization function in a GSL format (need to be defined externally!)
  double (*HybFunction)(double omega,
			void *params);

  double (*HybFuncWithEn)(double omega, 
  			  void *params);

  // Ok, this is a bummer but... 
  // hybridizations (Delta(en)) AND en*Delta(en) (or Delta(en)/en) 
  // have to be defined OUTSIDE as static members
  // in order to be passed in here.
  // So, for every static HybFunction (=Delta(e))
  // there HAS TO BE a HybFuncWithEn defined as well  

  // Tried hard to pass non-static function to GSL... 
  // ...doesn't work. Ditto (may 2010)
  //   double HFtimesEn(double omega, 
  // 		   void *params);
  
  //   static double HFtimesEnCallBack(double omega,
  // 				  void* params);


  ////////////

  double GetEn(int Nsites);
  double GetChin(int Nsites);

  // GSL integration
  double GSL_Integration(gsl_function *pIntegrandGSL,
			 gsl_integration_workspace *w_gsl, 
			 double x0, double x1);


  // Lanczos procedure
  void SetChainLanczos(int Nsitesmax);


  double GetHyb_w(double omega);

  double ScaleFactor(int n);


private:

  // GSL stuff
//   gsl_function HybFunGSL;
//   gsl_function EnHybF;
  gsl_integration_workspace *wsp_gsl;

  void SetHybRef();
  double CalcF2(gsl_integration_workspace *w_gsl);

  // Adds elements of *pFm and *pEm
  void CalcFabmEabm(int m,
		    gsl_integration_workspace *w_gsl,
		    char which,
		    vector<double> *pFm,
		    vector<double> *pEm);


  double RHS(int n, double eam, double um, double um1, 
	     double en, double chinm1, bool scale=true);


};

// struct NRGchainGSLCallbackHolder
// { 
//   //CNRGchain* cls;
//   //void* data;
//   CNRGchain cls;
//   void* data;

// };


#endif



/////////////////////////////////////////////////
////////////// Class CNRGthermo \\\\\\\\\\\\\\\\\
/////////////////////////////////////////////////


#ifndef _CNRGTHERMO_
#define _CNRGTHERMO_

class CNRGthermo{

public :

  //   // Static members

  // Input/Output
  char ThermoName[64];
  char ChainArqName[64];
  char ArqName[64];

  double dImpValue;
  bool CalcChain;
  int Nsite0Chain;

  double betabar;

   // 
  vector<double> dChainValues;
  vector<double> dValues;
  vector<double> dTempValues;

   // pointer to function:
  double (*Calc)(vector<double> Params, 
		 CNRGarray* pAeig, 
		 int sqnumber,bool totalS);


//   // Function members

  void AddValue(vector<double> Params, 
		CNRGarray* pAeig, 
		int sqnumber,bool totalS,
		double Temp);

  int ReadChain();

  void ReadNChainValue(int Nsites,int Nsites0);

  void SaveNValue(int Nsites, int Nsites0);

  void SaveToFile();

  bool CheckForFile();

};
// end CNRGthermo
#endif


//////////////////////////////////////////////////////
////////////// Class CNRGCodeHandler \\\\\\\\\\\\\\\\\
//////////////////////////////////////////////////////

#ifndef _CNRGCODEHANDLER_
#define _CNRGCODEHANDLER_

class CNRGCodeHandler{

public :

  // Static members
  char ModelOption[32];
  char Symmetry[32];
  char BandType[32];

  int ModelNo;
  int SymNo;
  int BandNo;


  // Files

  char ParamFileName[64];
  char BeginFileName[64];
  char EndFileName[64];
  char SaveArraysFileName[64];

  char LancInFileName[64];


  int NoInputParamsDouble;
  vector <double> dInitParams;
  //  vector <string> cParamsName;

  int calcdens;

  int Ncutoff;
  int Nsitesmax;
  int Nsites0;
  int UpdateBefCut;
  int Nsites;

  int NumNRGmats;
  int NumThermoMats;
  int NumChannels;
  int NopsSaved;

  double Lambda;
//  double Dband;
  double ChemPot;
  double DN;
  double betabar;

  double HalfLambdaFactor;

  // totalS QNumbers
  bool totalS;
  int Sqnumber;

  // instream (perhaps not such a good idea)
  std::ifstream CodeInFile;

  // CNRGarray objects (pointers)

  CNRGbasisarray* pAcut;
  CNRGbasisarray* pAbasis;
  CNRGmatrix* MatArray;


  vector<CNRGthermo> ThermoSTLArray; // New Sep 09



  bool SaveData;

  CNRGchain chain;

  // Function members

  //void InitialSetUp();
  // Trying this
  void InitialSetUp(bool ReadParamsOnly=false);


  void ReadParams();

  void ReadParams(char FileName[], int NoAdditionalParams);

  void SetChain();

  void SaveBegFile();

  void WrapUp();

  void SetSingleSite(CNRGbasisarray* pSingleSite);

  void CalcThermo(CNRGthermo* ThermoArray, 
		  CNRGarray* pAeig);

  void CalcStuff(CNRGarray* pAeig, CNRGbasisarray* pAeigCut);

  void DoesNothing();

  void ZeroParams();

  void SetTotS();

  //void SaveGenPars(char idname[]);
  void SaveGenPars();

  void SaveArrays();
  void SaveArrays(char idname[]);


  //void ReadGenPars(char idname[]);
  //void ReadGenPars();
  // Keep this until all ThisCodePars.dat have been updated :)
  void ReadGenPars(bool ReadModelBand=false);


  void ReadArrays();
  void ReadArrays(char idname[]);

  void SetCurrentDN();

  void SetNopsSaved();

  bool CheckFileExists(char arqname[]);

  // New: May 2010

  void ModelSwitch(  vector<int> &CommonQNs,
		     vector<int> &totSpos,   
		     CNRGarray* pAeig,
		     CNRGbasisarray* pSingleSite,
		     CNRGmatrix* pHN, 
		     vector<CNRGmatrix> &STLMatArray,
		     CNRGthermo* ThermoArray,
		     double &chi_m1);

  double chiN(int Nsites, double Lambda);


  void PrintSettings();

};

#endif
