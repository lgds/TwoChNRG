
#ifndef _MAIN_SETH0_
#define _MAIN_SETH0_

void OneChQ_SetAnderson_Hm1_old(vector<double> Params,
 			    CNRGarray* pAeig, 
 			    CNRGmatrix* NRGMats);

void OneChQ_SetAnderson_Hm1(vector<double> Params,
			    CNRGarray* pAeig, 
			    vector<CNRGmatrix> &STLNRGMats);



// void OneChQ_SetSMM_H0(vector<double> Params,
// 		       CNRGarray* pAeig, 
// 		       CNRGbasisarray* pSingleSite,
// 		       CNRGmatrix* NRGMats);
void OneChQ_SetSMM_H0(vector<double> Params,
		      CNRGarray* pAeig,
		      CNRGbasisarray* pSingleSite,
		      vector<CNRGmatrix> &STLNRGMats);


void TwoChQSP_SetH0CMphonon(vector<double> Params, 
			    CNRGarray* pAeig,
			    CNRGbasisarray* pSingleSite,
			    CNRGmatrix* NRGMats);



void OneChQSz_SetAnderson_Hm1(vector<double> Params,
			    CNRGarray* pAeig,
			    vector<CNRGmatrix> &STLNRGMats);


void OneChQSz_SetH0_DQD(vector<double> Params, 
			CNRGarray* pAeig,
			CNRGbasisarray* pSingleSite,
			vector<CNRGmatrix> &STLNRGMats);
// Includes Zeeman in BOTH dots

void OneChQSz_SetChainH0(CNRGarray* pAeig, 
			vector<CNRGmatrix> &STLNRGMats);




void OneChQS_SetH0_DQD(vector<double> Params, 
		       CNRGarray* pAeig,
		       CNRGbasisarray* pSingleSite,
		       vector<CNRGmatrix> &STLNRGMats);



void OneChQS_SetChainH0(CNRGarray* pAeig, 
			vector<CNRGmatrix> &STLNRGMats);



void OneChQS_SetAnderson_Hm1(vector<double> Params, 
			     CNRGarray* pAeig,
			     vector<CNRGmatrix> &STLNRGMats);

void OneChQS_SetKondoH0(vector<double> Params, 
			CNRGarray* pAeig,
			CNRGbasisarray* pSingleSite,
			vector<CNRGmatrix> &STLNRGMats);

void TwoChQS_SetKondoH0(vector<double> Params,
			CNRGbasisarray* pSingleSite,
			CNRGarray* pAeig,
			vector<CNRGmatrix> &STLNRGMats);

void TwoChQS_SetH0Chain(CNRGarray* pAeig, 
			CNRGbasisarray* pSingleSite, 
			vector<CNRGmatrix> &STLNRGMats);


void OneChNupPdn_SetH0_AndersonMajorana(vector<double> Params,
					CNRGbasisarray* pSingleSite,
					CNRGarray* pAeig, 
					vector<CNRGmatrix> &STLNRGMats);

// void OneChNupPdn_SetH0_AndersonMajorana(vector<double> Params,
//  					CNRGarray* pAeig, 
//  					vector<CNRGmatrix> &STLNRGMats);

void OneChNupPdn_SetChainH0(CNRGarray* pAeig, 
			    vector<CNRGmatrix> &STLNRGMats);


void OneChS_SetAnderson_Hm1(vector<double> Params,
			    CNRGarray* pAeig, 
			    vector<CNRGmatrix> &STLNRGMats);


void OneChSz_SetAnderson_Hm1(vector<double> Params,
			     CNRGarray* pAeig, 
			     vector<CNRGmatrix> &STLNRGMats);



#endif


#ifndef _MAIN_OPMATRULES_
#define _MAIN_OPMATRULES_

// bool OneChQ_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2);


/// Operators


double TwoChQSNoSz_fm1_dot1up_MatEl(CNRGbasisarray* pAbasis, 
				  CNRGbasisarray* pSingleSite,
				    int ist, int jst);

double TwoChQSNoSz_fm1_dot1dn_MatEl(CNRGbasisarray* pAbasis, 
				  CNRGbasisarray* pSingleSite,
				    int ist, int jst);

double TwoChQSNoSz_fm1_dot2up_MatEl(CNRGbasisarray* pAbasis, 
				  CNRGbasisarray* pSingleSite,
				    int ist, int jst);


double TwoChQSNoSz_fm1_dot2dn_MatEl(CNRGbasisarray* pAbasis, 
				  CNRGbasisarray* pSingleSite,
				    int ist, int jst);

///



/// Hamiltonians

// Moved to NRGOpMatRules.hpp
// double OneChQ_HN_MatEl(vector<double> Params,
// 		       CNRGbasisarray* pAbasis, 
// 		       CNRGbasisarray* pSingleSite,
// 		       CNRGmatrix* MatArray,
// 		       int ist, int jst);



double TwoChQSNoSz_Hm1SMM_MatEl(vector<double> Params,
			       CNRGbasisarray* pAbasis, 
			       CNRGbasisarray* pSingleSite,
			       CNRGmatrix* MatArray,
			       int ist, int jst);


double OneChQ_H0SMM_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst);

double TwoChQSP_H0ph_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
			   int ist, int jst);


// Moved to NRGOpMatRules.hpp
// double OneChQSz_HN_MatEl(vector<double> Params,
// 		       CNRGbasisarray* pAbasis, 
// 		       CNRGbasisarray* pSingleSite,
// 		       CNRGmatrix* MatArray,
// 		       int ist, int jst);


double OneChQS_H0DQD_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			   int ist, int jst);


double TwoDotQSz_Hm1_MatEl(vector<double> Params,
			   CNRGbasisarray* pAbasis, 
			   CNRGbasisarray* pSingleSite,
			   CNRGmatrix* MatArray,
			   int ist, int jst);


double OneChQSz_H0DQD_MatEl(vector<double> Params,
			    CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    CNRGmatrix* MatArray,
			    int ist, int jst);


complex<double> OneChNupPdn_Hm1_Majorana_MatEl(vector<double> Params,
					       CNRGbasisarray* pAbasis, 
					       CNRGbasisarray* pSingleSite,
					       CNRGmatrix* MatArray,
					       int ist, int jst);

#endif

#ifndef _MAIN_SETSINGLESITE_
#define _MAIN_SETSINGLESITE_

int OneChQSz_SetSingleSite(CNRGbasisarray &SingleSite);


int TwoChQSNoSz_SetSingleSite(CNRGbasisarray* pSingleSite);

// Not needed now...
void TwoChQSP_SetSingleSite(CNRGbasisarray &SingleSite);

// New
void TwoChQSz_SetSingleSite(CNRGbasisarray* pSingleSite);

#endif

#ifndef _TWODOTQS_SETINITIALSITE_
#define _TWODOTQS_SETINITIALSITE_

void TwoDotQS_SetInitialSite(CNRGbasisarray* pSingleSite);

#endif
