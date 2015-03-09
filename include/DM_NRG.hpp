#include "SpecFuncClass.hpp"

#ifndef _DM_NRG_COMMANDLINEREAD_
#define _DM_NRG_COMMANDLINEREAD_


void DM_NRG_CommandLineRead(int argc, char* argv[], int &Mtemp, 
			    double &betabar, double &twindow,
			    double &broadtemp,
			    double &bbroad,
			    int &UseCFS,
			    int &Nw);

#endif


#ifndef _DM_NRG_SETCHILDST_
#define _DM_NRG_SETCHILDST_


void DM_NRG_SetChildSt(CNRGbasisarray* pAcutN,
		       CNRGbasisarray* pAbasisNp1);


void DM_NRG_SetChildSameType(int ist_N, int jst_N,
			     CNRGbasisarray* pAcutN, 
			     CNRGbasisarray* pAcutNp1,
			     vector< vector<int> > &ChildSt_SameType);


#endif

#ifndef _DM_NRG_CALCRHON_
#define _DM_NRG_CALCRHON_


void DM_NRG_CalcRhoN_old(CNRGbasisarray* pAcutN,
		     CNRGbasisarray* pAcutNp1,
		     CNRGbasisarray* pAbasisNp1,
		     CNRGmatrix* pRhoN,
		     CNRGmatrix* pRhoNp1);


void DM_NRG_CalcRhoN(CNRGbasisarray* pAcutN,
		     CNRGbasisarray* pAcutNp1,
		     CNRGbasisarray* pAbasisNp1,
		     CNRGmatrix* pRhoN,
		     CNRGmatrix* pRhoNp1);

void DM_NRG_CalcRhoN_withSU2(CNRGbasisarray* pAcutN,
			     CNRGbasisarray* pAcutNp1,
			     CNRGbasisarray* pAbasisNp1,
			     CNRGbasisarray* pSingleSite,
			     CNRGmatrix* pRhoN,
			     CNRGmatrix* pRhoNp1);


#endif

#ifndef _DM_NRG_SETRHONMAX_
#define _DM_NRG_SETRHONMAX_

void DM_NRG_SetRhoNmax(vector<double> ParamsTemp,
		       CNRGbasisarray* pAcutNp1,
		       CNRGmatrix* pRhoNmax);


#endif

#ifndef _DM_NRG_CALCSPECFUNCS_
#define _DM_NRG_CALCSPECFUNCS_

void DM_NRG_CalcSpecFuncs(CNRGCodeHandler* pThisCode,
			  CNRGbasisarray* AcutN,
			  CNRGmatrix* RhoN,
			  CNRGmatrix** OpArrayN,
			  int iop, int jop);


// void DM_NRG_CalcSpecFunc_ij(CNRGCodeHandler* pThisCode,
// 			    CNRGbasisarray* AcutN,
// 			    CNRGmatrix* RhoN,
// 			    CNRGmatrix** OpArrayN,
// 			    int iop, int jop);

void DM_NRG_CalcSpecFunc_ij(CSpecFunction* pSpec,
			    CNRGmatrix** OpArrayN,
			    int iop, int jop, int UseCFS=0, int Nw=1);


#endif
