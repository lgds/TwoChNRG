// The goal is to put ALL checks and rules here... eventually.

#include "NRGclasses.hpp"

////////////////////////////////////////////
///////////     MatElChecks    /////////////
////////////////////////////////////////////


#ifndef _DIAGCHECK_
#define _DIAGCHECK_
bool Diag_check(CNRGbasisarray *pAeigCut, 
		      int iblock1, 
		int iblock2);

#endif



#ifndef _ONECHQSZ_CD_CHECK_
#define _ONECHQSZ_CD_CHECK_

bool OneChQSz_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

#endif

#ifndef _ONECHQSZ_CDUP_CHECK_
#define _ONECHQSZ_CDUP_CHECK_

bool OneChQSz_cdup_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

#endif

#ifndef _ONECHQSZ_CDDN_CHECK_
#define _ONECHQSZ_CDDN_CHECK_

bool OneChQSz_cddn_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

#endif

#ifndef _ONECHQ_CD_CHECK_
#define _ONECHQ_CD_CHECK_

bool OneChQ_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

#endif

#ifndef _ONECHQSOPMATRULES_ // moved from OneChQS.hpp
#define _ONECHQSOPMATRULES_

bool OneChQS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

// Not included in OpMatRules.cpp
//bool OneChQS_nd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);


double OneChQS_fN_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);

double OneChQS_cd_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);


double OneChQS_Sz_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst);


// Deprecated: Not included in OpMatRules.cpp (it's in OneChQS_OpMatRules.cpp)
double OneChQS_nd_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 		int ist, int jst);



#endif

#ifndef _ONECHNUPPDN_CD_CHECK_ 
#define _ONECHNUPPDN_CD_CHECK_

bool OneChNupPdn_cdup_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

bool OneChNupPdn_cddn_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);


#endif



#ifndef _ONECHSOPMATRULES_ 
#define _ONECHSOPMATRULES_

bool OneChS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

// Not included in OpMatRules.cpp
//bool OneChQS_nd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);


double OneChS_fN_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);

double OneChS_cd_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);

// Not included in OpMatRules.cpp
//double OneChS_nd_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 		int ist, int jst);


#endif


#ifndef _ONECHSZOPMATRULES_ 
#define _ONECHSZOPMATRULES_

bool OneChSz_cdup_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

bool OneChSz_cddn_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

double OneChSz_fN_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst, int isigma);

double OneChSz_fNup_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);

double OneChSz_fNdn_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);


double OneChSz_cd_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst,int isigma);

double OneChSz_cdup_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);

double OneChSz_cddn_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);



#endif




////////////////////////////////////////////
////////////////////////////////////////////
///////////   Mat ElCalculations   /////////
////////////////////////////////////////////
////////////////////////////////////////////


//////////////////////////////
///////  Operators   /////////
//////////////////////////////



#ifndef _IMPONLYMATEL_
#define _IMPONLYMATEL_
double ImpOnly_MatEl(CNRGbasisarray *pAbasis,
		     CNRGbasisarray *pSingleSite,
		     int ist, int jst);

complex<double> ImpOnly_MatElCplx(CNRGbasisarray *pAbasis,
				  CNRGbasisarray *pSingleSite,
				  int ist, int jst);
#endif

#ifndef _ALWAYSONEMATEL_
#define _ALWAYSONEMATEL_
double AlwaysOne_MatEl(CNRGbasisarray *pAbasis,
		       CNRGbasisarray *pSingleSite,
		       int ist, int jst);
#endif




#ifndef _ONECHQSZ_CD_MATEL_
#define _ONECHQSZ_CD_MATEL_

double OneChQSz_cd_MatEl(CNRGbasisarray* pAbasis, 
			 CNRGbasisarray* pSingleSite,
			 int ist, int jst, int isigma);

double OneChQSz_cdup_MatEl(CNRGbasisarray* pAbasis, 
			   CNRGbasisarray* pSingleSite,
			   int ist, int jst);


double OneChQSz_cddn_MatEl(CNRGbasisarray* pAbasis, 
			   CNRGbasisarray* pSingleSite,
			   int ist, int jst);


#endif

#ifndef _ONECHQSZ_FN_MATEL_
#define _ONECHQSZ_FN_MATEL_


double OneChQSz_fN_MatEl(CNRGbasisarray* pAbasis, 
			 CNRGbasisarray* pSingleSite,
			 int ist, int jst, int isigma);

double OneChQSz_fNup_MatEl(CNRGbasisarray* pAbasis, 
			   CNRGbasisarray* pSingleSite,
			   int ist, int jst);

double OneChQSz_fNdn_MatEl(CNRGbasisarray* pAbasis, 
			   CNRGbasisarray* pSingleSite,
			   int ist, int jst);

#endif


#ifndef _ONECHQ_CD_MATEL_
#define _ONECHQ_CD_MATEL_

double OneChQ_cd_MatEl(CNRGbasisarray* pAbasis, 
			 CNRGbasisarray* pSingleSite,
			 int ist, int jst);

#endif


#ifndef _ONECHQ_FN_MATEL_
#define _ONECHQ_FN_MATEL_

double OneChQ_fN_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst, int isigma);
double OneChQ_fNup_MatEl(CNRGbasisarray* pAbasis, 
			 CNRGbasisarray* pSingleSite,
			 int ist, int jst);
double OneChQ_fNdn_MatEl(CNRGbasisarray* pAbasis, 
			 CNRGbasisarray* pSingleSite,
			 int ist, int jst);

#endif


#ifndef _ONECHNUPPDN_CD_MATEL_ 
#define _ONECHNUPPDN_CD_MATEL_

double OneChNupPdn_cd_MatEl(CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    int ist, int jst, int isigma);

double OneChNupPdn_cdup_MatEl(CNRGbasisarray* pAbasis, 
			      CNRGbasisarray* pSingleSite,
			      int ist, int jst);

double OneChNupPdn_cddn_MatEl(CNRGbasisarray* pAbasis, 
			      CNRGbasisarray* pSingleSite,
			      int ist, int jst);

complex<double> OneChNupPdn_cdup_MatElCplx(CNRGbasisarray* pAbasis, 
					   CNRGbasisarray* pSingleSite,
					   int ist, int jst);

complex<double> OneChNupPdn_cddn_MatElCplx(CNRGbasisarray* pAbasis, 
					   CNRGbasisarray* pSingleSite,
					   int ist, int jst);


#endif

#ifndef _ONECHNUPPDN_FN_MATEL_ 
#define _ONECHNUPPDN_FN_MATEL_

double OneChNupPdn_fN_MatEl(CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    int ist, int jst, int isigma);

double OneChNupPdn_fNup_MatEl (CNRGbasisarray* pAbasis, 
			       CNRGbasisarray* pSingleSite,
			       int ist, int jst);

double OneChNupPdn_fNdn_MatEl (CNRGbasisarray* pAbasis, 
			       CNRGbasisarray* pSingleSite,
			       int ist, int jst);


complex<double> OneChNupPdn_fNup_MatElCplx(CNRGbasisarray* pAbasis, 
					   CNRGbasisarray* pSingleSite,
					   int ist, int jst);

complex<double> OneChNupPdn_fNdn_MatElCplx(CNRGbasisarray* pAbasis, 
					   CNRGbasisarray* pSingleSite,
					   int ist, int jst);



#endif


////////////////////////////////////
///////  Hamiltonians HN   /////////
////////////////////////////////////
// H0 Hamiltonians in NRG_main.hpp

#ifndef _ONECHQ_HN_MATEL_
#define _ONECHQ_HN_MATEL_ 

double OneChQ_HN_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
		       int ist, int jst);


#endif

#ifndef _ONECHQSZ_HN_MATEL_
#define _ONECHQSZ_HN_MATEL_ 
double OneChQSz_HN_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
		       int ist, int jst);
#endif

#ifndef _ONECHQS_HN_MATEL_
#define _ONECHQS_HN_MATEL_ 
double OneChQS_HN_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
		       int ist, int jst);
#endif


#ifndef _ONECHNUPPDN_HN_MATEL_
#define _ONECHNUPPDN_HN_MATEL_

complex<double> OneChNupPdn_HN_MatElCplx(vector<double> Params,
					 CNRGbasisarray* pAbasis, 
					 CNRGbasisarray* pSingleSite,
					 CNRGmatrix* MatArray,
					 int ist, int jst);


double OneChNupPdn_HN_MatEl(vector<double> Params,
			    CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    CNRGmatrix* MatArray,
			    int ist, int jst);



#endif

#ifndef _ONECHS_HN_MATEL_
#define _ONECHS_HN_MATEL_ 
double OneChS_HNsc_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
		       int ist, int jst);
#endif


#ifndef _ONECHSZ_HN_MATEL_
#define _ONECHSZ_HN_MATEL_ 
double OneChSz_HNsc_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst);
#endif
