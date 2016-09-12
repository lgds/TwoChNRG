
//#include "NRGclasses.hpp"
//using namespace std;


#ifndef _ONECHQSBUILDBASIS_
#define _ONECHQSBUILDBASIS_

int OneChQS_BuildBasis(CNRGarray Aeig,  CNRGbasisarray &Abasis, int Ncutoff);

#endif


#ifndef _ONECHQSDIAGHN_
#define _ONECHQSDIAGHN_
int OneChQS_DiagHN( CNRGmatrix Qm1fNQ, 
		    CNRGbasisarray Abasis,  CNRGbasisarray SingleSite,
		    CNRGarray &Aeig,
		    double eps_N, double chi_N, double Lambda,bool display=false);

#endif

#ifndef _ONECHQSPREFACTORHN_
#define _ONECHQSPREFACTORHN_

double OneChQS_prefactorHN(int type_i, int type_j, double S);
#endif

#ifndef _ONECHQSUPDATEQM1FQ_
#define _ONECHQSUPDATEQM1FQ_

//void OneChQS_UpdateQm1fQ(CNRGmatrix &Qm1fNQ,CNRGarray Aeig,
//			 CNRGbasisarray Abasis, CNRGbasisarray SingleSite);

void OneChQS_UpdateQm1fQ(CNRGmatrix* pQm1fNQ,CNRGarray* pAeig,
			 CNRGbasisarray* pAbasis, CNRGbasisarray* pSingleSite);

void OneChQS_UpdateMatrixAfterCutting(CNRGmatrix* pQm1fNQ,CNRGmatrix* pMQQp1,
				      CNRGbasisarray* pAeigCut,
				      CNRGbasisarray* pAbasis, CNRGbasisarray* pSingleSite);


#endif 

#ifndef _ONECHQSSETH0_
#define _ONECHQSSETH0_

void OneChQS_SetAndersonHm1(vector<double> Params,
 			    CNRGarray* pAeig, CNRGmatrix* pQm1fNQ, 
			    CNRGmatrix* NRGMats);

//			    CNRGmatrix* pMQQp1);


void OneChQS_SetKondoH0(vector<double> Params,
			CNRGarray* pAeig, CNRGmatrix* pQm1fNQ, 
			CNRGbasisarray* pSingleSite );


void OneChQS_SetH0Chain(CNRGarray* pAeig, 
			CNRGmatrix* pQm1fNQ);


void OneChQS_SetSMM_Hm1(vector<double> Params,
		       CNRGarray* pAeig, 
		       CNRGmatrix* pQm1fNQ,
			CNRGmatrix* NRGMats);


#endif

#ifndef _ONECHQSSETSINGLESITE_
#define _ONECHQSSETSINGLESITE_
int OneChQS_SetSingleSite(CNRGbasisarray &SingleSite);
#endif

// Operator Tables

#ifndef _FD_TABLE_
#define _FD_TABLE_
double fd_table(int sigma, int typep, int type);
#endif


#ifndef _ONECHQSOPMATRULES_
#define _ONECHQSOPMATRULES_

bool OneChQS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

bool OneChQS_nd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);


double OneChQS_fN_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);

double OneChQS_cd_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);

double OneChQS_nd_MatEl(CNRGbasisarray *pAbasis, CNRGbasisarray *pSingleSite, 
			int ist, int jst);


#endif
