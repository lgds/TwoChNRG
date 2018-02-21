
#ifndef _TWOCHQSZSETSINGLESITE_
#define _TWOCHQSZSETSINGLESITE_


int TwoChQSz_SetSingleSite(CNRGbasisarray* pSingleSite);

#endif


#ifndef _TWOCHQSZSETH0ANDERSON_
#define _TWOCHQSZSETH0ANDERSON_

//int TwoChQSz_SetH0Anderson(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGmatrix* Qm1fNQ);
int TwoChQSz_SetH0Anderson(vector<double> Params,CNRGbasisarray* pSingleSite,CNRGarray* pAeig,CNRGbasisarray* pAbasis);

#endif

#ifndef _TWOCHQSZOPTABLES_
#define _TWOCHQSZOPTABLES_

double fd_table(int channel, int sigma, int type_i, int type_j);

#endif

#ifndef _TWOCHQSZ_UPDATEMATRICES_
#define _TWOCHQSZ_UPDATEMATRICES_

void TwoChQSz_UpdateQm1fQ(CNRGbasisarray* pSingleSite, CNRGarray* pAeig, 
			  CNRGbasisarray* pAbasis,CNRGmatrix* Qm1fNQ);


void TwoChQSz_UpdateMatrixAfterCutting(CNRGbasisarray* pSingleSite,
				       CNRGbasisarray* pAeigCut, 
				       CNRGbasisarray* pAbasis,
				       CNRGmatrix* Qm1fNQ, 
				       CNRGmatrix* pMQQp1);

#endif


#ifndef _TWOCHQSDIAGHN_
#define _TWOCHQSDIAGHN_

void TwoChQSz_DiagHN(vector<double> Params, 
		     CNRGbasisarray* pAbasis,CNRGbasisarray* pSingleSite,
		     CNRGmatrix* Qm1fNQ, CNRGarray* pAeig);
#endif
