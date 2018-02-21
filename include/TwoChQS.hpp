#ifndef _TWOCHQSSETSINGLESITE_
#define _TWOCHQSSETSINGLESITE_


int TwoChQS_SetSingleSite(CNRGbasisarray* pSingleSite);

#endif

#ifndef _TWOCHQSSETH0ANDERSON_
#define _TWOCHQSSETH0ANDERSON_

int TwoChQS_SetH0Anderson(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasis);


void TwoChQS_SetH0CMphonon(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0,CNRGmatrix* NRGMats);

void TwoChQS_SetH0Kondo(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0,CNRGmatrix* MatArray);


void TwoChQS_SetH0Chain(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasis);

void TwoChQS_SetH0CMphononwTransf(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasisH0,CNRGmatrix* NRGMats);



#endif

#ifndef _TWOCHQS_UPDATEMATRICES_
#define _TWOCHQS_UPDATEMATRICES_

void TwoChQS_UpdateQm1fQ(CNRGbasisarray* pSingleSite, CNRGarray* pAeig, 
			  CNRGbasisarray* pAbasis,CNRGmatrix* Qm1fNQ);


void TwoChQS_UpdateMatrixAfterCutting(CNRGbasisarray* pSingleSite,
				       CNRGbasisarray* pAeigCut, 
				       CNRGbasisarray* pAbasis,
				       CNRGmatrix* Qm1fNQ, 
				       CNRGmatrix* pMQQp1);

#endif

#ifndef _TWOCHQSDIAGHN_
#define _TWOCHQSDIAGHN_

void TwoChQS_DiagHN(vector<double> Params, 
		     CNRGbasisarray* pAbasis,CNRGbasisarray* pSingleSite,
		     CNRGmatrix* Qm1fNQ, CNRGarray* pAeig);
#endif

#ifndef _TWOCHQS_OPMATRULES_
#define _TWOCHQS_OPMATRULES_


bool TwoChQS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);


bool TwoChQSP_fA_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);

double TwoChQS_d_fd_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst);


double TwoChQS_H0ph_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst);

double TwoChQS_H0phwTransf_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst);



double TwoChQS_H0Kondo_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			     int ist, int jst);

double TwoChQS_HN_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			int ist, int jst);


double TwoChQS_cdPhonon_MatEl(CNRGbasisarray* pAbasis, 
			      CNRGbasisarray* pSingleSite,
			      int ist, int jst, int ich);

double TwoChQS_cd_ich1_Phonon_MatEl(CNRGbasisarray* pAbasis, 
				    CNRGbasisarray* pSingleSite,
				    int ist, int jst);

double TwoChQS_cd_ich2_Phonon_MatEl(CNRGbasisarray* pAbasis, 
				    CNRGbasisarray* pSingleSite,
				    int ist, int jst);


double TwoChQS_fN_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
			int ist, int jst, int ich);

double TwoChQS_fNch1_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
			   int ist, int jst);

double TwoChQS_fNch2_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
			   int ist, int jst);



#endif

#ifndef _TWOCHQS_OPTABLES_
#define _TWOCHQS_OPTABLES_

double TwoChQS_fd_table(int channel, int sigma, int type_i, int type_j);


double TwoChQS_SpSm_table(int type_i, int type_j);


double TwoChQS_fdupfdn_table(int channel, 
			     int type_i, int type_j);

double TwoChQS_Szch_table(int channel, int type_i);



#endif


#ifndef _TWOCHQSPSETSINGLESITE_
#define _TWOCHQSPSETSINGLESITE_


int TwoChQSP_SetSingleSite(CNRGbasisarray* pSingleSite);

#endif
