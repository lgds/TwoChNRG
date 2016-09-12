///////////////////////////////
//                           //
//                           //
// Useful little functions   //
//                           //
///////////////////////////////


#include "NRGclasses.hpp"


#ifndef _R2IJ_
#define _R2IJ_
int ij2r(int Nel, int i, int j);
int ij2rNSq(int Neli, int Nelj, int i, int j);
int r2ij(int Nel, int r, int &i, int &j);

int Pmn(int Nel, int i, int j);
int ij2rUpTr(int Nel, int i, int j);
#endif


// LAPACK routines

#ifndef _DSYEV_
#define _DSYEV_
extern "C" void  dsyev_(char *JOBZ, char *UPLO, int *N, double *A, 
			int *LDA, double *W, double *WORK, int *LWORK, int *INFO);

#endif

#ifndef _DSPEV_
#define _DSPEV_
extern "C" void  dspev_(char *JOBZ, char *UPLO, int *N, double *AP,  
			 double *W, double *Z, int *LDZ, double *WORK, int *INFO);

#endif

#ifndef _DGESV_
#define _DGESV_
extern "C" void  dgesv_(int *N, int *NRHS, double *A, int *LDA, 
			double *IPIV, double *B, int *LDB, int *INFO );

// DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
//
// *  N       (input) INTEGER
// *          The number of linear equations, i.e., the order of the
// *          matrix A.  N >= 0.
// *
// *  NRHS    (input) INTEGER
// *          The number of right hand sides, i.e., the number of columns
// *          of the matrix B.  NRHS >= 0.
// *
// *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
// *          On entry, the N-by-N coefficient matrix A.
// *          On exit, the factors L and U from the factorization
// *          A = P*L*U; the unit diagonal elements of L are not stored.
// *
// *  LDA     (input) INTEGER
// *          The leading dimension of the array A.  LDA >= max(1,N).
// *
// *  IPIV    (output) INTEGER array, dimension (N)
// *          The pivot indices that define the permutation matrix P;
// *          row i of the matrix was interchanged with row IPIV(i).
// *
// *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
// *          On entry, the N-by-NRHS matrix of right hand side matrix B.
// *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
// *
// *  LDB     (input) INTEGER
// *          The leading dimension of the array B.  LDB >= max(1,N).
// *
// *  INFO    (output) INTEGER
// *          = 0:  successful exit
// *          < 0:  if INFO = -i, the i-th argument had an illegal value
// *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
// *               has been completed, but the factor U is exactly
// *                singular, so the solution could not be computed.
// *
// *  =====================================================================



#endif


#ifndef _ZHEEV_
#define _ZHEEV_
// SUBROUTINE ZHEEV( JOBZ, UPLO,N, A, LDA, W, WORK, LWORK,RWORK, INFO)
//       CHARACTER	    JOBZ, UPLO
//       INTEGER	    INFO, LDA, LWORK, N
//       DOUBLE	    PRECISION RWORK( * ), W( * )
//       COMPLEX*16    A( LDA, * ), WORK( * )

extern "C" void  zheev_(char *JOBZ, char *UPLO, int *N, complex<double> *A, 
			int *LDA, double *W, complex<double> *WORK, int *LWORK, 
			double *RWORK, int *INFO);

#endif

#ifndef _ZHPEV_
#define _ZHPEV_
//   SUBROUTINE ZHPEV( JOBZ, UPLO,N, AP, W, Z, LDZ, WORK,RWORK, INFO )
//       CHARACTER	    JOBZ, UPLO
//       INTEGER	    INFO, LDZ, N
//       DOUBLE	    PRECISION RWORK( * ), W( * )
//       COMPLEX*16    AP(	* ), WORK( * ),	Z( LDZ,	* )

extern "C" void  zhpev_(char *JOBZ, char *UPLO, int *N, complex<double> *A, 
			double *W, complex<double> *Zvec, int *LDZ,
			complex<double> *WORK,double *RWORK, int *INFO);

#endif




#ifndef _CUTSTATES_
#define _CUTSTATES_

CNRGbasisarray CutStates(CNRGarray* Ain, int Ncutoff);

void CutStates(CNRGarray &Ain, CNRGbasisarray &Aout, int Ncutoff);


#endif

#ifndef _QS_BUILDBASIS_
#define _QS_BUILDBASIS_


 int QS_BuildBasis(CNRGbasisarray* pAeigCut,  
 		  CNRGbasisarray* pAbasis, 
                   CNRGbasisarray* pSingleSite,
 		  int BefCut, int KeepSz=0);



// int QS_BuildBasis(CNRGbasisarray &AeigCut,  
// 		  CNRGbasisarray &Abasis, 
//                   CNRGbasisarray &SingleSite,
// 		  int BefCut);

//		  int Ncutoff, int BefCut);


#endif


#ifndef _BUILDBASIS_
#define _BUILDBASIS_

//void BuildBasis_RedSymmetry(vector<int> CommonQNs, 
//
// Should work even if symmetry is reduced!
// Allows for totS
//
void BuildBasis(vector<int> CommonQNs, vector<int> totSpos,
		CNRGbasisarray* pAeigCut,  
		CNRGbasisarray* pAbasis, 
		CNRGbasisarray* pSingleSite, int BefCut);

#endif

#ifndef _Q1Q2Sz_BUILDBASIS_
#define _Q1Q2Sz_BUILDBASIS_

int Q1Q2Sz_BuildBasis(CNRGbasisarray* pAeigCut,  CNRGbasisarray* pAbasis, 
		      CNRGbasisarray* pSingleSite, int BefCut);

#endif

#ifndef _QSz_BUILDBASIS_
#define _QSz_BUILDBASIS_

int QSz_BuildBasis(CNRGbasisarray* pAeigCut,  CNRGbasisarray* pAbasis, 
		   CNRGbasisarray* pSingleSite, int BefCut);


#endif

#ifndef _CGordan_
#define _CGordan_

double CGordan(double J1, double M1, double J2, double M2, double J, double M);

bool TriangIneq(double J1, double J2, double J);

#endif


#ifndef _CGordan_SStot_
#define _CGordan_SStot_

double CGordan_SStot(double S, double Stilde, double Sztilde, double Stot);

#endif

#ifndef _DEQUAL_
#define _DEQUAL_
bool dEqual(double A, double B);
bool dNEqual(double A, double B);
bool dLEqual(double A, double B);
bool dGEqual(double A, double B);

bool dEqualPrec(double A, double B, double prec);
bool dNEqualPrec(double A, double B, double prec);
bool dLEqualPrec(double A, double B, double prec);
bool dGEqualPrec(double A, double B, double prec);

bool dGT(double A, double B);
bool dLT(double A, double B);

bool dGTPrec(double A, double B, double prec);
bool dLTPrec(double A, double B, double prec);

#endif

#ifndef _QS_CALCSUSCEP_
#define _QS_CALCSUSCEP_

double QS_CalcSuscep(vector<double> Params, CNRGarray* pAeig);


#endif

#ifndef _CALCTHERMO_
#define _CALCTHERMO_

double CalcSuscep(vector<double> Params, CNRGarray* pAeig, 
		  int sqnumber,bool totalS);

double CalcEntropy(vector<double> Params, CNRGarray* pAeig, 
		  int sqnumber,bool totalS);
#endif


#ifndef _SDOTS_
#define _SDOTS_

double Sdots_totalS(double ST, double STz, double Sold, double Stilde);

double Sdots_Sz(double Sold, double Stilde,
		double Szoldp,  double Sztildep,
		double Szold,  double Sztilde);

double Splus(double Sp, double Szp,
	      double S,  double Sz);

double Sminus(double Sp, double Szp,
	      double S,  double Sz);


#endif

#ifndef _OPTABLES_
#define _OPTABLES_

double OneCh_fd_table(int sigma, int type_i, int type_j);

double TwoCh_fd_table(int channel, int sigma, int type_i, int type_j);

double TwoDotQS_fd_reduced_table(int channel, int type_i, int type_j);

double OneChS_fd_table(int sigma, int type_i, int type_j);

double OneChS_fdupfddn_table(int type_i, int type_j);


#endif

#ifndef _TWODOTQSSZ_FD_TABLE_
#define _TWODOTQSSZ_FD_TABLE_

double TwoDotQSSz_fd_table(int channel, int sigma, int type_i, int type_j);

#endif

#ifndef _QFDQM1TOTS_
#define _QFDQM1TOTS_

double QfdQm1_totS(vector<double> Params, 
		   vector<int> Indexes, 
		   CNRGbasisarray* pSingleSite);

#endif


#ifndef _UPDATEMATRICES_
#define _UPDATEMATRICES_

void UpdateMatrices(CNRGbasisarray* pSingleSite,CNRGbasisarray* pAeigCut, 
		    CNRGbasisarray* pAbasis,
		    CNRGmatrix* NRGMats, int NumNRGMats, bool display=false);

#endif

#ifndef _FINDMATCHBLOCK_
#define _FINDMATCHBLOCK_


int FindMatchBlock(CNRGarray* pAeigCut, int iblock, 
		   CNRGbasisarray* pAbasis);

#endif

#ifndef _CALCOPAVG_
#define _CALCOPAVG_

double CalcOpAvg(vector<double> Params, 
		 CNRGbasisarray* pAeig, CNRGmatrix* pAop, 
		 bool totalS,int sqnumber);

#endif


#ifndef _DIAGCHECK_
#define _DIAGCHECK_
bool Diag_check(CNRGbasisarray *pAeigCut, 
		      int iblock1, 
		int iblock2);

#endif

// should not be here...
//#ifndef _IMPONLYMATEL_
//#define _IMPONLYMATEL_
//double ImpOnly_MatEl(CNRGbasisarray *pAbasis,
//		     CNRGbasisarray *pSingleSite,
//		     int ist, int jst);
//#endif



#ifndef _CALCPHMATEL_
#define _CALCPHMATEL_

double Calc_phMatEl(double mi, double mj, double lambda, double w0);

double Calc_apad(double mi, double mj);


#endif

#ifndef _FILEEXISTS_
#define _FILEEXISTS_

bool FileExists(char arqname[]);

#endif
////////////////////////


#ifndef _IONRGARRAY_
#define _IONRGARRAY_

void SaveNRGarrayBin(CNRGarray *pArray,char arqname[]);

void ReadNRGarrayBin(CNRGarray *pArray,char arqname[]);

#endif


#ifndef _BROADDELTA_
#define _BROADDELTA_

double BroadDelta(double omega, double Ep, double b);

double BroadDelta2GSL(double omega, void *pars);

double LorentzDelta(double omega, double Ep, double b);


double GaussDelta(double omega, double Ep, double b);


double LorentzDeltaAnders(double omega, double Ep, double b);


#endif


#ifndef _CALCDN_
#define _CALCDN_

double CalcDN(double Lambda,int Nshell, double z_twist=1.0);

#endif


#ifndef _ROTATEMATRIX_
#define _ROTATEMATRIX_

void RotateMatrix(CNRGmatrix* pMat, CNRGbasisarray* pAcut,
		  CNRGmatrix *pMatRot);

// Rotate in the UNCUT eigenvector
void RotateMatrix_NoCut(CNRGmatrix* pMat, CNRGarray* pAeig,
			CNRGmatrix* pMatRot, bool forward);

#endif


#ifndef _GSLINTEGRATOR_
#define _GSLINTEGRATOR_

double GSL_Integrator(gsl_function *pIntegrandGSL,
		      gsl_integration_workspace *w_gsl, 
		      double x0, double x1);

#endif


////////////////////////
