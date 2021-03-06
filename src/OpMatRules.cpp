
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"


///////////////////
//  Checks      /// 
///////////////////

////////////////////////

bool Diag_check(CNRGbasisarray *pAeigCut, 
		      int iblock1, 
		      int iblock2)
{

  if ( iblock1==iblock2 )
    return(true);
  else
    return(false);
}


////////////////////////////////////////////

////////////////////////
////////////////////////
// OneChQSz check


bool OneChQSz_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

  double Qi=pAeigCut->GetQNumber(iblock1,0);
  double Szi=pAeigCut->GetQNumber(iblock1,1);

  double Qj=pAeigCut->GetQNumber(iblock2,0);
  double Szj=pAeigCut->GetQNumber(iblock2,1);

  
  if (  (dEqual(Qj,(Qi+1.0)))&&
	( (dEqual(Szj,(Szi+0.5)))||(dEqual(Szj,(Szi-0.5))) )
	) 
    return(true);
  else
    return(false);

}

////////////////////////////////////////////


////////////////////////
bool OneChQSz_cdup_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

  double Qi=pAeigCut->GetQNumber(iblock1,0);
  double Szi=pAeigCut->GetQNumber(iblock1,1);

  double Qj=pAeigCut->GetQNumber(iblock2,0);
  double Szj=pAeigCut->GetQNumber(iblock2,1);

  
  if (  (dEqual(Qj,(Qi+1.0)))&&(dEqual(Szj,(Szi+0.5))) )
    return(true);
  else
    return(false);

}

////////////////////////////////////////////
////////////////////////
bool OneChQSz_cddn_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

  double Qi=pAeigCut->GetQNumber(iblock1,0);
  double Szi=pAeigCut->GetQNumber(iblock1,1);

  double Qj=pAeigCut->GetQNumber(iblock2,0);
  double Szj=pAeigCut->GetQNumber(iblock2,1);

  
  if (  (dEqual(Qj,(Qi+1.0)))&&(dEqual(Szj,(Szi-0.5))) )
    return(true);
  else
    return(false);

}


////////////////////////////////////////////
// OneChQ check

bool OneChQ_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

  double Qi=pAeigCut->GetQNumber(iblock1,0);
  double Qj=pAeigCut->GetQNumber(iblock2,0);
  
  if ( dEqual(Qj,(Qi+1.0)) ) 
    return(true);
  else
    return(false);

}

////////////////////////
// OneChQS check

bool OneChQS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

  double Qi=pAeigCut->GetQNumber(iblock1,0);
  double Si=pAeigCut->GetQNumber(iblock1,1);

  double Qj=pAeigCut->GetQNumber(iblock2,0);
  double Sj=pAeigCut->GetQNumber(iblock2,1);

  
  if (  (dEqual(Qj,(Qi+1.0)))&&
        ( (dEqual(Sj,(Si+0.5)))||(dEqual(Sj,(Si-0.5))) )
        ) 
    return(true);
  else
    return(false);

}

////////////////////////////////////////////

////////////////////////
// OneChNupPdn check

////////////////////////
bool OneChNupPdn_cdup_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){
  // cdup|up>=|0> (not dagger)
  // cdup|Nup Pdn>=+/-|Nup-1 Pdn>

  double Nupi=pAeigCut->GetQNumber(iblock1,0);
  double Pdni=pAeigCut->GetQNumber(iblock1,1);

  double Nupj=pAeigCut->GetQNumber(iblock2,0);
  double Pdnj=pAeigCut->GetQNumber(iblock2,1);

  
  if (  (dEqual(Nupj,(Nupi+1.0)))&&(dEqual(Pdni,Pdnj)) )
    return(true);
  else
    return(false);

}

////////////////////////
bool OneChNupPdn_cddn_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){
  // cddn|Nup Pdn>=+/-|Nup -Pdn>

  double Nupi=pAeigCut->GetQNumber(iblock1,0);
  double Pdni=pAeigCut->GetQNumber(iblock1,1);

  double Nupj=pAeigCut->GetQNumber(iblock2,0);
  double Pdnj=pAeigCut->GetQNumber(iblock2,1);

  
  if (  (dEqual(Nupj,Nupi))&&(dEqual(Pdni,-1.0*Pdnj)) )
    return(true);
  else
    return(false);

}


////////////////////////////////////////////

////////////////////////
// OneChS check

bool OneChS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

  double Si=pAeigCut->GetQNumber(iblock1,0);
  double Sj=pAeigCut->GetQNumber(iblock2,0);

  if ( (dEqual(Sj,(Si+0.5)))||(dEqual(Sj,(Si-0.5))) ) 
    return(true);
  else
    return(false);

}

////////////////////////////////////////////

////////////////////////
// OneChSz check
bool OneChSz_cdup_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

  double Szi=pAeigCut->GetQNumber(iblock1,0);
  double Szj=pAeigCut->GetQNumber(iblock2,0);

  if ( dEqual(Szj,(Szi+0.5)) ) 
    return(true);
  else
    return(false);

}
bool OneChSz_cddn_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

  double Szi=pAeigCut->GetQNumber(iblock1,0);
  double Szj=pAeigCut->GetQNumber(iblock2,0);

  if ( dEqual(Szj,(Szi-0.5)) ) 
    return(true);
  else
    return(false);

}

////////////////////////
// OneChPupPdn check

////////////////////////
bool OneChPupPdn_cdup_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){
  // cdup|up>=|0> (not dagger)
  // cdup|Pup Pdn>=+/-|-Pup Pdn>

  double Pupi=pAeigCut->GetQNumber(iblock1,0);
  double Pdni=pAeigCut->GetQNumber(iblock1,1);

  double Pupj=pAeigCut->GetQNumber(iblock2,0);
  double Pdnj=pAeigCut->GetQNumber(iblock2,1);

  
  if (  (dEqual(Pupj,-1.0*Pupi))&&(dEqual(Pdni,Pdnj)) )
    return(true);
  else
    return(false);

}

////////////////////////
bool OneChPupPdn_cddn_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){
  // cddn|Pup Pdn>=+/-|Pup -Pdn>

  double Pupi=pAeigCut->GetQNumber(iblock1,0);
  double Pdni=pAeigCut->GetQNumber(iblock1,1);

  double Pupj=pAeigCut->GetQNumber(iblock2,0);
  double Pdnj=pAeigCut->GetQNumber(iblock2,1);

  
  if (  (dEqual(Pupj,Pupi))&&(dEqual(Pdni,-1.0*Pdnj)) )
    return(true);
  else
    return(false);

}


////////////////////////////////////////////




////////////////////////////////
//////                    //////
////// Operator Rules     //////
//////                    //////
////////////////////////////////

////////////////////////

double ImpOnly_MatEl(CNRGbasisarray *pAbasis,
		     CNRGbasisarray *pSingleSite,
		     int ist, int jst)
{
                                                                   
  double fbasis=0.0;
                                                                               
  int typei=pAbasis->iType[ist];
  int typej=pAbasis->iType[jst];
                                                                               
  // Not diagonal necessarily
  if (typei==typej) fbasis=1.0;
 
  return(fbasis);
 
}

complex<double> ImpOnly_MatElCplx(CNRGbasisarray *pAbasis,
				  CNRGbasisarray *pSingleSite,
				  int ist, int jst)
{
                                                                   
  double fbasis=0.0;
                                                                               
  int typei=pAbasis->iType[ist];
  int typej=pAbasis->iType[jst];
                                                                               
  // Not diagonal necessarily
  if (typei==typej) fbasis=1.0;
 
  return(complex<double>(fbasis,0.0));
 
}

 
////////////////////////

////////////////////////

double AlwaysOne_MatEl(CNRGbasisarray *pAbasis,
		       CNRGbasisarray *pSingleSite,
		       int ist, int jst)
{
  // returns fbasis=1.0 always. 
  // Use with empty single site and NeedOld=true
  return(1.0);
}
 
////////////////////////





//////////////////////////////////////////
/////////  OneChQSz symmetry   ///////////
//////////////////////////////////////////


////////////////////////

double OneChQSz_cd_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst, int isigma){

  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Szi=pAbasis->GetQNumberFromSt(ist,1);
  //int stcfi=pAbasis->StCameFrom[ist];
  int typei=pAbasis->iType[ist];

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Szj=pAbasis->GetQNumberFromSt(jst,1);
  //int stcfj=pAbasis->StCameFrom[jst];
  int typej=pAbasis->iType[jst];

  double Qtildej=pSingleSite->GetQNumber(typej,0);

//   cout << " isigma = " << isigma
//        << " stcfi = " << stcfi
//        << " stcfj = " << stcfj
//        << " typei = " << typei
//        << " typej = " << typej
//        << endl
//        << " <ti|f+|tj> = " << 
//     OneCh_fd_table(isigma,typej,typei)
//        << endl;

  double dSigma=((double)isigma)*0.5;

  // Have to make this more clear but let's see:
  double FermiSign=1.0;
  if (dEqual(Qtildej,0.0)) FermiSign=-1.0;

  if (typei==typej){
    if ( (dEqual(Qj,Qi+1.0))&&(dEqual(Szj,Szi+dSigma)) )
      return(FermiSign);
    else return(0.0);
  }
  else
    return(0.0);

}

//////


double OneChQSz_cdup_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChQSz_cd_MatEl(pAbasis,pSingleSite,
			 ist, jst,1));
}

double OneChQSz_cddn_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChQSz_cd_MatEl(pAbasis,pSingleSite,
			 ist, jst,-1));
}
//////


////////////////////////

double OneChQSz_fN_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst, int isigma){

  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Szi=pAbasis->GetQNumberFromSt(ist,1);
  int stcfi=pAbasis->StCameFrom[ist];
  int typei=pAbasis->iType[ist];

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Szj=pAbasis->GetQNumberFromSt(jst,1);
  int stcfj=pAbasis->StCameFrom[jst];
  int typej=pAbasis->iType[jst];

  if ( (stcfi==stcfj) ){
    return(OneCh_fd_table(isigma,typej,typei));
  }else
    return(0.0);
}

//////

double OneChQSz_fNup_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChQSz_fN_MatEl(pAbasis,pSingleSite,
			 ist, jst,1));
}

double OneChQSz_fNdn_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChQSz_fN_MatEl(pAbasis,pSingleSite,
			 ist, jst,-1));
}
//////


/////////////////////////////////
///
/// <Q' Sz'|HN|Q Sz> :
///
/////////////////////////////////
/////////////////////////////////
double OneChQSz_HN_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
		       int ist, int jst){

  // Feb 09
  // Calculates matrix elements for HN when Sz is explicitly included
  //
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  eN            - Params[2] (added)
  //  <Qm1 Sz+/-1/2|fN_up|Q Sz> - MatArray[0] (not reduced)
  //  <Qm1Sz+/-1/2|fN_dn|Q Sz> - MatArray[1] (not reduced)

  double MatEl=0.0;

  double Lambda=Params[0];
  double chi_N=Params[1];
  double eN=Params[2];


  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Szi=pAbasis->GetQNumberFromSt(ist,1);

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Szj=pAbasis->GetQNumberFromSt(jst,1);

  if ( (dNEqual(Qi,Qj))||(dNEqual(Szi,Szj)) ) return(0.0);

  if (ist==jst){
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
  // add eN term here
  int type_i=pAbasis->iType[ist];
  double Qtilde=pSingleSite->GetQNumberFromSt(type_i,0);
  MatEl+=eN*(Qtilde+1.0);
  }
  else{
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    // Sum over spins
    int imat=0;
    for (int isigma=1;isigma>=-1;isigma-=2){
      int typep=type_i;
      int type=type_j;

      int stcf_i=pAbasis->StCameFrom[ist];
      int stcf_j=pAbasis->StCameFrom[jst];

      // Matrix element from <istcf| fN_sigma | jstcf>
      double OldEl=MatArray[imat].GetMatEl(stcf_i,stcf_j);

      // if zero, try the h.c term, meaning f+N f_N+1 in original case
      if (dEqual(OldEl,0.0)){
	OldEl=MatArray[imat].GetMatEl(stcf_j,stcf_i);
	typep=type_j;
	type=type_i;
      }

      // Get SingleSite QNumbers

      double Qtilde=pSingleSite->GetQNumberFromSt(type,0);
      double Sztilde=pSingleSite->GetQNumberFromSt(type,1);


      // Fermi Sign
      double FermiSign=1.0;
      if (dEqual(Qtilde,0.0)) FermiSign=-1.0;

      // Debug
      // 	if (ist==1){
      // 	  cout << " ist = " << ist 
      // 	       << " jst = " << jst
      // 	       << " stcf_i = " << stcf_i
      // 	       << " stcf_j = " << stcf_j
      // 	       << " type_i = " << type_i
      // 	       << " type_j = " << type_j
      // 	       << " sigma = " << isigma << endl;
      // 	  cout << " typep = " << typep
      // 	       << " type = " << type 
      // 	       << " OldEl = " << OldEl
      // 	       << endl;
      // 	  cout << " fd_table = "
      // 	       << OneCh_fd_table(isigma,typep,type) 
      // 	       << endl;
      // 	  cout << " Qtilde = " << Qtilde << endl;
      // 	  cout << " FermiSign = " << FermiSign << endl;
      // 	  cout << " chi_N = " << chi_N << endl;
      // 	}

      MatEl+=FermiSign*OldEl*OneCh_fd_table(isigma,typep,type);
      imat++;

    }
    // End sum over isigma
      
    MatEl*=chi_N;
  }
  // end if i=j


  return(MatEl);
}

/////////////////////////////////


/////////////////////////////////////////////
/////////  OneChNupPdn symmetry   ///////////
/////////////////////////////////////////////

////////////////////////

double OneChNupPdn_cd_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst, int isigma){

  double Nupi=pAbasis->GetQNumberFromSt(ist,0);
  double Pdni=pAbasis->GetQNumberFromSt(ist,1);
  int stcfi=pAbasis->StCameFrom[ist];
  int typei=pAbasis->iType[ist];

  double Nupj=pAbasis->GetQNumberFromSt(jst,0);
  double Pdnj=pAbasis->GetQNumberFromSt(jst,1);
  int stcfj=pAbasis->StCameFrom[jst];
  int typej=pAbasis->iType[jst];


//   cout << " isigma = " << isigma
//        << " stcfi = " << stcfi
//        << " stcfj = " << stcfj
//        << " typei = " << typei
//        << " typej = " << typej
//        << endl
//        << " <ti|f+|tj> = " << 
//     OneCh_fd_table(isigma,typej,typei)
//        << endl;

  double dSigma=((double)isigma)*0.5;

  // Have to make this more clear but let's see:
  // if |typej> = |up> or |dn>, we sign=-1
  double FermiSign=1.0;
  if ((typej==1)||(typej==2)) FermiSign=-1.0;

  if (typei==typej){
    if (   (  (isigma==1)&&(dEqual(Nupj,(Nupi+1.0)))&&(dEqual(Pdni,Pdnj)) )||
	   (  (isigma==-1)&&(dEqual(Nupj,Nupi))&&(dEqual(Pdni,-1.0*Pdnj)) )   )
      return(FermiSign);
    else return(0.0);
  }
  else
    return(0.0);

}

//////
double OneChNupPdn_cdup_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChNupPdn_cd_MatEl(pAbasis,pSingleSite,
			 ist, jst,1));
}

double OneChNupPdn_cddn_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChNupPdn_cd_MatEl(pAbasis,pSingleSite,
			 ist, jst,-1));
}

complex<double> OneChNupPdn_cdup_MatElCplx(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  double aux=OneChNupPdn_cd_MatEl(pAbasis,pSingleSite,
				  ist, jst,1);
  return(complex<double>(aux,0.0));
}

complex<double> OneChNupPdn_cddn_MatElCplx(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){

  double aux=OneChNupPdn_cd_MatEl(pAbasis,pSingleSite,
				  ist, jst,-1);
  return(complex<double>(aux,0.0));
}



//////




double OneChNupPdn_fN_MatEl(CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    int ist, int jst, int isigma){

//   double Nupi=pAbasis->GetQNumberFromSt(ist,0);
//   double Pdni=pAbasis->GetQNumberFromSt(ist,1);
  int stcfi=pAbasis->StCameFrom[ist];
  int typei=pAbasis->iType[ist];

//   double Nupj=pAbasis->GetQNumberFromSt(jst,0);
//   double Pdnj=pAbasis->GetQNumberFromSt(jst,1);
  int stcfj=pAbasis->StCameFrom[jst];
  int typej=pAbasis->iType[jst];

  if ( (stcfi==stcfj) ){
      return(OneCh_fd_table(isigma,typej,typei));
  }else return(0.0);



}

//////

double OneChNupPdn_fNup_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChNupPdn_fN_MatEl(pAbasis,pSingleSite,
			 ist, jst,1));
}

double OneChNupPdn_fNdn_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChNupPdn_fN_MatEl(pAbasis,pSingleSite,
			 ist, jst,-1));
}

complex<double> OneChNupPdn_fNup_MatElCplx(CNRGbasisarray* pAbasis, 
					   CNRGbasisarray* pSingleSite,
					   int ist, int jst){
  double aux=OneChNupPdn_fN_MatEl(pAbasis,pSingleSite,
				  ist, jst,1);
  return(complex<double>(aux,0.0));
}

complex<double> OneChNupPdn_fNdn_MatElCplx(CNRGbasisarray* pAbasis, 
					   CNRGbasisarray* pSingleSite,
					   int ist, int jst){
  double aux=OneChNupPdn_fN_MatEl(pAbasis,pSingleSite,
				  ist, jst,-1);
  return(complex<double>(aux,0.0));
}
//////

///
/// <Nup' Pdn'|HN|Nup Pdn> :
///
/////////////////////////////////
/////////////////////////////////
complex<double> OneChNupPdn_HN_MatElCplx(vector<double> Params,
			    CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    CNRGmatrix* MatArray,
			    int ist, int jst){

  // Feb 05/2014
  // Calculates matrix elements for HN in the Nup Pdn basis
  //
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  eN            - Params[2] (added)
  //  <Nup-1 Pdn|fN_up|Nup Pdn> - MatArray[0] (not reduced)
  //  <Nup -Pdn|fN_dn|Nup Pdn> - MatArray[1] (not reduced)

  complex<double> cMatEl=(0.0,0.0);

  double Lambda=Params[0];
  double chi_N=Params[1];
  double eN=Params[2];


  double Nupi=pAbasis->GetQNumberFromSt(ist,0);
  double Pdni=pAbasis->GetQNumberFromSt(ist,1);

  double Nupj=pAbasis->GetQNumberFromSt(jst,0);
  double Pdnj=pAbasis->GetQNumberFromSt(jst,1);



  if ( (dNEqual(Nupi,Nupj))||(dNEqual(Pdni,Pdnj)) ) return(0.0);

  if (ist==jst){
    double MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
    // add eN term here
    int type_i=pAbasis->iType[ist];
    double Nuptilde=pSingleSite->GetQNumberFromSt(type_i,0);
    double Pdntilde=pSingleSite->GetQNumberFromSt(type_i,1);
    // Calculate Qtilde
    double Qtilde=0.0;
    // WRONG!
//     if ( (dEqual(Nuptilde,0.0))&&(dEqual(Pdntilde,-1.0)) ) Qtilde=-1.0;
//     if ( (dEqual(Nuptilde,1.0))&&(dEqual(Pdntilde,1.0)) ) Qtilde=1.0;
    if ( (dEqual(Nuptilde,0.0))&&(dEqual(Pdntilde,1.0)) ) Qtilde=-1.0; //|0>
    if ( (dEqual(Nuptilde,1.0))&&(dEqual(Pdntilde,-1.0)) ) Qtilde=1.0; // | up dn>
    MatEl+=eN*(Qtilde+1.0);
    cMatEl=complex<double>(MatEl,0.0);
  }else{
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    // Sum over spins
    int imat=0;
    for (int isigma=1;isigma>=-1;isigma-=2){
      int typep=type_i;
      int type=type_j;

      int stcf_i=pAbasis->StCameFrom[ist];
      int stcf_j=pAbasis->StCameFrom[jst];

      // Matrix element from <istcf| fN_sigma | jstcf>
      complex<double> cOldEl=MatArray[imat].cGetMatEl(stcf_i,stcf_j);

      complex<double> cOldEl_old=cOldEl;

      // if zero, try the h.c term, meaning f+N f_N+1 in original case
      // Ops, it can be non-zero but the wrong one!!
//       if ( (dEqualPrec(cOldEl.real(),0.0,1e-10))&&(dEqualPrec(cOldEl.imag(),0.0,1e-10)) ){
// 	// debug
// 	cOldEl=MatArray[imat].cGetMatEl(stcf_j,stcf_i);
// 	typep=type_j;
// 	type=type_i;
//       }
      // end if
      double OpTable=OneCh_fd_table(isigma,typep,type);
      if ( dEqual(OpTable,0.0) ){
	// debug
	cOldEl=MatArray[imat].cGetMatEl(stcf_j,stcf_i);
	typep=type_j;
	type=type_i;
      }
      // end if




      // Get SingleSite QNumbers

      double Nuptilde=pSingleSite->GetQNumberFromSt(type,0);
      double Pdntilde=pSingleSite->GetQNumberFromSt(type,1);

      // Fermi Sign
      double FermiSign=1.0;
      // 	if (dEqual(Qtilde,0.0)) FermiSign=-1.0;
      // |tilde> is |up>=|Nup=1 Pdn=1> or |dn>=|Nup=0 Pdn=-1>
      if (  ( (dEqual(Nuptilde,1.0))&&(dEqual(Pdntilde,1.0)) )||
	    ( (dEqual(Nuptilde,0.0))&&(dEqual(Pdntilde,-1.0)) )
	   )  
	FermiSign=-1.0;

      cMatEl+=FermiSign*cOldEl*OneCh_fd_table(isigma,typep,type);

      // Debug 
//       if (  (pAbasis->Nshell==2)&&
//  	    (imat==1)&&
// 	    ( ( (ist==160)&&(jst>=187)&&(jst<=187) )
// 	      )  ){
// 	cout << " imat = " << imat
// 	     << " ist= " << ist
// 	     << " jst= " << jst
// 	     << " stcfi= " << stcf_i
// 	     << " stcfj= " << stcf_j
// 	     << " typei= " << type_i
// 	     << " typej= " << type_j 
//  	     << " typep= " << typep
// 	     << " type= " << type 
// 	     << endl
// 	     << " <stcfi|f_dn|scfj>= " << MatArray[imat].cGetMatEl(stcf_i,stcf_j)
// 	     << " <stcfj|f_dn|scfi>= " << MatArray[imat].cGetMatEl(stcf_j,stcf_i)
// 	     << " cOldEl_old= " << cOldEl_old
// 	     << " cOldEl= " << cOldEl 
// 	     << endl
// 	     << " Fermi sign= " << FermiSign
// 	     << " |Re(cOldEl_old)| < 1e-10? " << dEqualPrec(cOldEl_old.real(),0.0,1e-10)
// 	     << " |Im(cOldEl_old)| < 1e-10? " << dEqualPrec(cOldEl_old.imag(),0.0,1e-10)
// 	     << " |Re(cOldEl_old)| < 1e-15? " << dEqualPrec(cOldEl_old.real(),0.0,1e-15)
// 	     << " |Im(cOldEl_old)| < 1e-15? " << dEqualPrec(cOldEl_old.imag(),0.0,1e-15)
// 	     << " (f_dn)^d_table= " << OneCh_fd_table(isigma,typep,type)
// 	     << " Mat El/chi_N= " << cMatEl
// 	     << " Mat El= " << cMatEl*chi_N 
// 	     << endl;
//       }
      ///

      imat++;
    }
    // End sum over isigma
      
    cMatEl*=chi_N;
  }
  // end if i=j


  return(cMatEl);
}

/////////////////////////////////

double OneChNupPdn_HN_MatEl(vector<double> Params,
			    CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    CNRGmatrix* MatArray,
			    int ist, int jst){

  // Feb 05/2014
  // Calculates matrix elements for HN in the Nup Pdn basis
  //
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  eN            - Params[2] (added)
  //  <Nup-1 Pdn|fN_up|Nup Pdn> - MatArray[0] (not reduced)
  //  <Nup -Pdn|fN_dn|Nup Pdn> - MatArray[1] (not reduced)

  double MatEl=0.0;

  double Lambda=Params[0];
  double chi_N=Params[1];
  double eN=Params[2];


  double Nupi=pAbasis->GetQNumberFromSt(ist,0);
  double Pdni=pAbasis->GetQNumberFromSt(ist,1);

  double Nupj=pAbasis->GetQNumberFromSt(jst,0);
  double Pdnj=pAbasis->GetQNumberFromSt(jst,1);



  if ( (dNEqual(Nupi,Nupj))||(dNEqual(Pdni,Pdnj)) ) return(0.0);

  if (ist==jst){
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
    // add eN term here
    int type_i=pAbasis->iType[ist];
    double Nuptilde=pSingleSite->GetQNumberFromSt(type_i,0);
    double Pdntilde=pSingleSite->GetQNumberFromSt(type_i,1);
    // Calculate Qtilde
    double Qtilde=0.0;
    if ( (dEqual(Nuptilde,0.0))&&(dEqual(Pdntilde,-1.0)) ) Qtilde=-1.0;
    if ( (dEqual(Nuptilde,1.0))&&(dEqual(Pdntilde,1.0)) ) Qtilde=1.0;
    MatEl+=eN*(Qtilde+1.0);
  }else{
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    // Sum over spins
    int imat=0;
    for (int isigma=1;isigma>=-1;isigma-=2){
      int typep=type_i;
      int type=type_j;

      int stcf_i=pAbasis->StCameFrom[ist];
      int stcf_j=pAbasis->StCameFrom[jst];

      // Matrix element from <istcf| fN_sigma | jstcf>
      double OldEl=MatArray[imat].GetMatEl(stcf_i,stcf_j);

      double OldEl_old=OldEl;

      // if zero, try the h.c term, meaning f+N f_N+1 in original case
      // Ops, it can be non-zero but the wrong one!!
//       if ( (dEqualPrec(cOldEl.real(),0.0,1e-10))&&(dEqualPrec(cOldEl.imag(),0.0,1e-10)) ){
// 	// debug
// 	cOldEl=MatArray[imat].cGetMatEl(stcf_j,stcf_i);
// 	typep=type_j;
// 	type=type_i;
//       }
      // end if
      double OpTable=OneCh_fd_table(isigma,typep,type);
      if ( dEqual(OpTable,0.0) ){
	// debug
	OldEl=MatArray[imat].GetMatEl(stcf_j,stcf_i);
	typep=type_j;
	type=type_i;
      }
      // end if


      // Get SingleSite QNumbers

      double Nuptilde=pSingleSite->GetQNumberFromSt(type,0);
      double Pdntilde=pSingleSite->GetQNumberFromSt(type,1);

      // Fermi Sign
      double FermiSign=1.0;
      // 	if (dEqual(Qtilde,0.0)) FermiSign=-1.0;
      // |tilde> is |up>=|Nup=1 Pdn=1> or |dn>=|Nup=0 Pdn=-1>
      if (  ( (dEqual(Nuptilde,1.0))&&(dEqual(Pdntilde,1.0)) )||
	    ( (dEqual(Nuptilde,0.0))&&(dEqual(Pdntilde,-1.0)) )
	   )  
	FermiSign=-1.0;

      MatEl+=FermiSign*OldEl*OneCh_fd_table(isigma,typep,type);
      imat++;
    }
    // End sum over isigma
      
    MatEl*=chi_N;
  }
  // end if i=j

  return(MatEl);
}





/////////////////////////////////////////
/////////   OneChQ symmetry   ///////////
/////////////////////////////////////////


////////////////////////
////////////////////////

double OneChQ_cd_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){

  // For both spin up and down. 
  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  int typei=pAbasis->iType[ist];

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  int typej=pAbasis->iType[jst];

  // WRONG! typej marks the state not the block (not one state per block!!)
  // 
  //double Qtildej=pSingleSite->GetQNumber(typej,0);
  double Qtildej=pSingleSite->GetQNumberFromSt(typej,0);


  // Have to make this more clear but let's see:
  // site states | up > or | dn > give an extra Fermi sign. 
  double FermiSign=1.0;
  if (dEqual(Qtildej,0.0)) FermiSign=-1.0;

  if ( (ist==2)&&(jst==7) )
    cout << " Qtildej= "<< Qtildej
	 << " typej= " << typej
	 << " FermiSign = " << FermiSign
	 << endl;

  if (typei==typej){
    if (dEqual(Qj,Qi+1.0))
      return(FermiSign);
    else return(0.0);
  }
  else
    return(0.0);

}
////////////////////////

/////////////////////////////////
///
/// <Q'|fN_sigma |Q> : 
/// <Q'Sz'|fN_sigma |Q Sz> : 
///

double OneChQ_fN_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst, int isigma){

  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  int stcfi=pAbasis->StCameFrom[ist];
  int typei=pAbasis->iType[ist];

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  int stcfj=pAbasis->StCameFrom[jst];
  int typej=pAbasis->iType[jst];

  double DQ=fabs(Qi-Qj);
//   cout << " isigma = " << isigma
//        << " stcfi = " << stcfi
//        << " stcfj = " << stcfj
//        << " typei = " << typei
//        << " typej = " << typej
//        << endl
//        << " <ti|f+|tj> = " << 
//     OneCh_fd_table(isigma,typej,typei)
//        << endl;

  if ( (stcfi==stcfj) )
    {
      return(OneCh_fd_table(isigma,typej,typei));
    }
  else
    return(0.0);

}

//////

double OneChQ_fNup_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChQ_fN_MatEl(pAbasis,pSingleSite,
			 ist, jst,1));

}
//////

double OneChQ_fNdn_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChQ_fN_MatEl(pAbasis,pSingleSite,
			 ist, jst,-1));

}

//////////////////////////////



////////////////////////

///
/// <Q'|HN|Q> :
///
/////////////////////////////////
/////////////////////////////////
double OneChQ_HN_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
		       int ist, int jst){

  // July 08
  // Calculates matrix elements for HN when Sz is explicitly included
  //
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  <Qm1|fN_up|Q> - MatArray[0] (not reduced)
  //  <Qm1|fN_dn|Q> - MatArray[1] (not reduced)

  double MatEl=0.0;

  double Lambda=Params[0];
  double chi_N=Params[1];


  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Qj=pAbasis->GetQNumberFromSt(jst,0);

  if (dNEqual(Qi,Qj)) return(0.0);

  if (ist==jst)
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
  else{
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    // Sum over spins
    int imat=0;
    for (int isigma=1;isigma>=-1;isigma-=2){
      int typep=type_i;
      int type=type_j;

      int stcf_i=pAbasis->StCameFrom[ist];
      int stcf_j=pAbasis->StCameFrom[jst];

      // Matrix element from <istcf| fN_sigma | jstcf>
      double OldEl=MatArray[imat].GetMatEl(stcf_i,stcf_j);

      // if zero, try the h.c term, meaning f+N f_N+1 in original case
      if (dEqual(OldEl,0.0)){
	OldEl=MatArray[imat].GetMatEl(stcf_j,stcf_i);
	typep=type_j;
	type=type_i;
      }

      // Get SingleSite QNumbers

      double Qtilde=pSingleSite->GetQNumberFromSt(type,0);

      // Fermi Sign
      double FermiSign=1.0;
      if (dEqual(Qtilde,0.0)) FermiSign=-1.0;

      // 	  cout << " ist = " << ist 
      // 	       << " jst = " << jst
      // 	       << " stcf_i = " << stcf_i
      // 	       << " stcf_j = " << stcf_j
      // 	       << " type_i = " << type_i
      // 	       << " type_j = " << type_j
      // 	       << " sigma = " << isigma << endl;
      // 	  cout << " typep = " << typep
      // 	       << " type = " << type 
      // 	       << " OldEl = " << OldEl
      // 	       << endl;
      // 	  cout << " fd_table = "
      // 	       << OneCh_fd_table(isigma,typep,type) 
      // 	       << endl;
      // 	  cout << " Qtilde = " << Qtilde << endl;
      // 	  cout << " FermiSign = " << FermiSign << endl;
      // 	  cout << " chi_N = " << chi_N << endl;

      MatEl+=FermiSign*OldEl*OneCh_fd_table(isigma,typep,type);
      imat++;

    }
    // End sum over isigma
      
    MatEl*=chi_N;
  }
  // end if i=j

  return(MatEl);

}

/////////////////////////////////



/////////////////////////////////////////
/////////  OneChQS symmetry   ///////////
/////////////////////////////////////////

////////////////////////
////////////////////////

double OneChQS_fN_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){

  double fbasis=0.0;

  double Si=pAbasis->GetQNumberFromSt(ist,1);
  int typei=pAbasis->iType[ist];
  int stcfi=pAbasis->StCameFrom[ist];

  double Sj=pAbasis->GetQNumberFromSt(jst,1);
  int typej=pAbasis->iType[jst];
  int stcfj=pAbasis->StCameFrom[jst];

  if (stcfi==stcfj)
    {
  
      // New approach
      double Qtildei=pSingleSite->GetQNumber(typei,0);
      double Stildei=pSingleSite->GetQNumber(typei,1);
      double Sztildei=pSingleSite->GetQNumber(typei,2);
                         
      double Qtildej=pSingleSite->GetQNumber(typej,0);
      double Stildej=pSingleSite->GetQNumber(typej,1);
      double Sztildej=pSingleSite->GetQNumber(typej,2);

      double siteqnumsi[]={Qtildei,Stildei,Sztildei};
      double siteqnumsj[]={Qtildej,Stildej,Sztildej};
  
      double Scfi=Si-Sztildei;
      double Scfj=Sj-Sztildej;
  
      double auxBasis[2]={0.0,0.0};
      double auxCG[2]={0.0,0.0};

      // Let me try this... Works!
      double sigma=Sj-Si;



      for (Sztildej=-Stildej;
           Sztildej<=Stildej;Sztildej+=1.0)
        {
          siteqnumsj[2]=Sztildej;
          int sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);
          //cout << " Sz~j = "<<Sztildej << " ssj =  "<< sitestatej << endl; 
          int i1=0;


          siteqnumsi[2]=Sztildej-sigma;
          int sitestatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);
          
          auxCG[0]=CGordan(Scfi,Si-Sztildej+sigma,Stildei,Sztildej-sigma,Si,Si);
          auxCG[1]=CGordan(Scfj,Sj-Sztildej,Stildej,Sztildej,Sj,Sj);
          auxBasis[i1]+=OneCh_fd_table((int)(sigma/0.5),sitestatej,sitestatei)*auxCG[0]*auxCG[1];

          //cout << " CG 1 = " << auxCG[0] << " CG 2 = " << auxCG[1] << " mat_el = " <<  OneCh_fd_table((int)(sigma/0.5),sitestatej,sitestatei) << endl;
          i1++;
        }
      // end loop in Sztildej
      // Check to see which of them is non-zero
      auxCG[0]=CGordan(Si,Si,0.5,sigma,Sj,Sj);
      if ( !(dEqual(auxCG[0],0.0)) ) 
        fbasis=auxBasis[0]/auxCG[0];
      else
        {
          if ( !(dEqual(auxCG[1],0.0)) ) 
            fbasis=auxBasis[1]/auxCG[1];
          else fbasis=0.0;
        }

    }
  else fbasis=0.0;
  // end if stcfi==stcfj

  return(fbasis);

}

////////////////////////
////////////////////////

double OneChQS_cd_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){

  double fbasis=0.0;

  double Si=pAbasis->GetQNumberFromSt(ist,1);
  int typei=pAbasis->iType[ist];
  int stcfi=pAbasis->StCameFrom[ist];

  double Sj=pAbasis->GetQNumberFromSt(jst,1);
  int typej=pAbasis->iType[jst];
  int stcfj=pAbasis->StCameFrom[jst];

  double Qtildei=pSingleSite->GetQNumber(typei,0);
  double Stildei=pSingleSite->GetQNumber(typei,1);
  double Sztildei=pSingleSite->GetQNumber(typei,2);
                         
  double Qtildej=pSingleSite->GetQNumber(typej,0);
  double Stildej=pSingleSite->GetQNumber(typej,1);
  double Sztildej=pSingleSite->GetQNumber(typej,2);


  double siteqnumsi[]={Qtildei,Stildei,Sztildei};
  double siteqnumsj[]={Qtildej,Stildej,Sztildej};
  
  double Scfi=Si-Sztildei;
  double Scfj=Sj-Sztildej;
  
  double auxCG[4]={0.0,0.0,0.0,0.0};

  // Let me try this... Works!
  double sigma=Sj-Si;

  auxCG[3]=CGordan(Si,Si,0.5,sigma,Sj,Sj);

  // Check sigma=Si-Sj
  if (dEqual(auxCG[3],0.0)){
    sigma=Si-Sj;
    auxCG[3]=CGordan(Si,Si,0.5,sigma,Sj,Sj);
    if (dEqual(auxCG[3],0.0)) return(0.0);
  }


  for (double Szcfj=-Scfj;
       Szcfj<=Scfj;Szcfj+=1.0){
    double Szcfi=Szcfj-sigma;
    Sztildej=Sj-Szcfj;
    Sztildei=Si-Szcfi;
    
    siteqnumsj[2]=Sztildej;
    int sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);

    siteqnumsi[2]=Sztildei;
    int sitestatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);

    auxCG[0]=CGordan(Scfi,Szcfi,Stildei,Sztildei,Si,Si);
    auxCG[1]=CGordan(Scfj,Szcfj,Stildej,Sztildej,Sj,Sj);
    auxCG[2]=CGordan(Scfi,Szcfi,0.5,sigma,Scfj,Szcfj);

    // Have to make this more clear but let's see:
    double FermiSign=1.0;
    // Jul 09: Does FermiSign enter here????
    // Yes, it does.
    if (dEqual(Qtildej,0.0)) FermiSign=-1.0;


    if (sitestatei==sitestatej)
      //       fbasis+=auxCG[0]*auxCG[1]*auxCG[2];
      fbasis+=auxCG[0]*auxCG[1]*auxCG[2]*FermiSign;

    //cout << " CG 1 = " << auxCG[0] << " CG 2 = " << auxCG[1] << " mat_el = " <<  OneCh_fd_table((int)(sigma/0.5),sitestatej,sitestatei) << endl;
  }
  // end loop in Szcfj

  if (dNEqual(auxCG[3],0.0)) fbasis=fbasis/auxCG[3];
  else fbasis=0.0;

  return(fbasis);

}

////////////////////////
////////////////////////

double OneChQS_Sz_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){
  // Calculate <Q S  a|| Sz ||Q S b>:
  // 
  // The trick:
  //  
  // <Q S Sz a| Szimp |Q S Sz b> = <S Sz 1 0| S Sz><Q S a||Szimp||Q S b>
  //
  // Keep only reduced matrix elements
  //
  //

  double fbasis=0.0;

  double Si=pAbasis->GetQNumberFromSt(ist,1);
  int typei=pAbasis->iType[ist];
  int stcfi=pAbasis->StCameFrom[ist];

  double Sj=pAbasis->GetQNumberFromSt(jst,1);
  int typej=pAbasis->iType[jst];
  int stcfj=pAbasis->StCameFrom[jst];

  double auxCG[3]={0.0,0.0,0.0};

  auxCG[2]=CGordan(Si,Si,1.0,0.0,Si,Si);

  if ( (typei==typej)&&(dNEqual(auxCG[2],0.0)) ){
    double Qtildei=pSingleSite->GetQNumber(typei,0);
    double Stildei=pSingleSite->GetQNumber(typei,1);
    double Sztildei=pSingleSite->GetQNumber(typei,2);
                         
    double siteqnumsi[]={Qtildei,Stildei,Sztildei};
  
    double Scfi=Si-Sztildei; // Ok, I will buy that
  

    // Loop in Szold
    for (double Szcfi=-Scfi; Szcfi<=Scfi;Szcfi+=1.0){
      
      Sztildei=Si-Szcfi;
          
      auxCG[0]=CGordan(Scfi,Szcfi,Stildei,Sztildei,Si,Si);
      auxCG[1]=CGordan(Scfi,Szcfi,1.0,0.0,Scfi,Szcfi);
  
      fbasis+=auxCG[0]*auxCG[0]*auxCG[1];

//       cout << " ist = " << ist << " jst = " << jst
// 	   << " typei = " << typei << "  typej = " <<  typej << endl
// 	   << " Scfi = " << Scfi << " Szcfi = " << Szcfi 
// 	   << " Stildei = " << Stildei << " Sztildei = " << Sztildei << endl 
// 	   << " CG 0 = " << auxCG[0] << " CG 1 = " << auxCG[1] << " CG 2 = " 
// 	   << auxCG[2] << " mat_el = " <<  fbasis << endl;
    }
    // end loop in Szcfi

    fbasis=fbasis/auxCG[2]; // reduced element

  }
  else fbasis=0.0;
  // end if stcfi==stcfj

  return(fbasis);

}


////////////////////////
////////////////////////

///////////////////////////////////////////
/////////  <Q' S' | HN | Q S>   ///////////
///////////////////////////////////////////



double OneChQS_HN_MatEl(vector<double> Params,
                          CNRGbasisarray* pAbasis,
                          CNRGbasisarray* pSingleSite,
                          CNRGmatrix* MatArray,
                          int ist, int jst)
{
  // Sept 09
  //
  // Finally implementing OneChQS_HN_MatEl (it's embedded in OneChQS_DiagHN in earlier versions
  // Calculates matrix elements for HN in the QS basis
  //
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  eN            - Params[2] (NEW! 2012)
  //  <Qm1 S'||fN||Q S> - MatArray[0] (reduced)

  double MatEl=0.0;

  double Lambda=Params[0];
  double chi_N=Params[1];
  double eN=Params[2]; // NEW! (Jan 2012)


  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Si=pAbasis->GetQNumberFromSt(ist,1);

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Sj=pAbasis->GetQNumberFromSt(jst,1);

  if (dNEqual(Qi,Qj)||dNEqual(Si,Sj)) return(0.0);

  if (ist==jst){
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
    // add eN term here
    int type_i=pAbasis->iType[ist];
    double Qtilde=pSingleSite->GetQNumberFromSt(type_i,0);
    MatEl+=eN*(Qtilde+1.0);
  }
  else{
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    // Matrix element from <istcf|| fN ||jstcf>
    double OldEl=MatArray[0].GetMatEl(pAbasis->StCameFrom[ist],
				      pAbasis->StCameFrom[jst]);
    int typep=type_i;
    int type=type_j;

    // if zero, try the h.c term
    if (dEqual(OldEl,0.0)){
      OldEl=MatArray[0].GetMatEl(pAbasis->StCameFrom[jst],
				 pAbasis->StCameFrom[ist]);
      typep=type_j;
      type=type_i;
    }

    // Get SingleSite QNumbers

    double Qtildep=pSingleSite->GetQNumber(typep,0);
    double Stildep=pSingleSite->GetQNumber(typep,1);
    double Sztildep=pSingleSite->GetQNumber(typep,2);

    double Qtilde=pSingleSite->GetQNumber(type,0);
    double Stilde=pSingleSite->GetQNumber(type,1);
    double Sztilde=pSingleSite->GetQNumber(type,2);

    double Scfp=Si-Sztildep;
    double Scf=Si-Sztilde;

    // Check Fermi Sign
    double FermiSign=1.0;
    if (dEqual(Qtilde,0.0)) FermiSign=-1.0;
    // Fermi sign appears IF there is an
    // odd number of fermionic ops in site state

    double CoefMatEl=0.0;
    double auxCG[]={0.0,0.0,0.0,0.0};
    double siteqnumsp[]={Qtildep,Stildep,Sztildep};
    double siteqnums[]={Qtilde,Stilde,Sztilde};
    // Sum in Sztilde
    for (Sztilde=-Stilde;Sztilde<=Stilde;Sztilde+=1.0){
      // Sum in sigma: Note that Sztildep=Sztilde+sigma
      siteqnums[2]=Sztilde;
      // Only one site per block in 1chQS
      int sitestate=pSingleSite->GetBlockFromQNumbers(siteqnums); 
      for (int sigma=-1;sigma<=1;sigma+=2){
	double dSigma=0.5*(double)sigma;
	double Szcf=Si-Sztilde;
	double Szcfp=Szcf-dSigma;
	Sztildep=Sztilde+dSigma;

	siteqnumsp[2]=Sztilde+dSigma;
	// Test Find site state
	int sitestatep=pSingleSite->GetBlockFromQNumbers(siteqnumsp);

	//cout << "Sz~ = " << Sztilde 
	//   << " dSigma = " << dSigma << endl;
	
	auxCG[0]=CGordan(Scf,Szcf,
			 Stilde,Sztilde,Si,Si);

	auxCG[1]=CGordan(Scfp,Szcfp,
			 Stildep,Sztildep,Si,Si);

	auxCG[2]=CGordan(Scfp,Szcfp,0.5,dSigma,Scf,Szcf);
 
	auxCG[3]=OneCh_fd_table(sigma,sitestatep,sitestate);
	//cout << "<sitestp|f|sitest> = " << auxCG[3] << endl;
	CoefMatEl+=auxCG[0]*auxCG[1]*auxCG[2]*auxCG[3];
				
      }
      // end loop in sigma
    }
    // end loop in Sztilde

    //auxEl*=CoefMatEl;
    MatEl=chi_N*OldEl*CoefMatEl*FermiSign;
    //cout << " MatEl = " << auxEl << endl;
  }
  // end if i=j

  
  return(MatEl);

}

////////////////////////
////////////////////////








////////////////////////
////////////////////////////////////////
/////////  OneChS symmetry   ///////////
////////////////////////////////////////

////////////////////////
////////////////////////

double OneChS_fN_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){

  double fbasis=0.0;

  double Si=pAbasis->GetQNumberFromSt(ist,0);
  int typei=pAbasis->iType[ist];
  int stcfi=pAbasis->StCameFrom[ist];

  double Sj=pAbasis->GetQNumberFromSt(jst,0);
  int typej=pAbasis->iType[jst];
  int stcfj=pAbasis->StCameFrom[jst];

  // Debug
  bool disp=false;
//   if (   ( (ist==5)&&(jst==1) )||
// 	 ( (ist==6)&&(jst==3) )||
// 	 ( (ist==7)&&(jst==4) )
// 	 )
//     disp=true;

  if (stcfi==stcfj){
  
    if (disp) 
      cout << "I'm in: " 
	   << " ist=" << ist << " jst= " << jst 
	   << " typei=" << typei << " typej= " << typej
	   << " Si= " << Si << " Sj= " << Sj
	   << " stcfi=" << stcfi << " stcfj= " << stcfj 
	   << endl; 
    // New approach
    double Stildei=pSingleSite->GetQNumberFromSt(typei,0);
    double Sztildei=pSingleSite->GetQNumberFromSt(typei,1);
                         
    double Stildej=pSingleSite->GetQNumberFromSt(typej,0);
    double Sztildej=pSingleSite->GetQNumberFromSt(typej,1);

    double siteqnumsi[]={Stildei,Sztildei};
    double siteqnumsj[]={Stildej,Sztildej};
  
    double Scfi=Si-Sztildei;
    double Scfj=Sj-Sztildej;
  
    double auxBasis[2]={0.0,0.0};
    double auxCG[2]={0.0,0.0};

    // Let me try this... Works!
    double sigma=Sj-Si;

    if (disp) 
      cout << " Stildei=" << Stildei << " Sztildei= " << Sztildei
	   << " Stildej=" << Stildej << " Sztildej= " << Sztildej
	   << endl 
	   << " Scfi=" << Scfi << " Scfj= " << Scfj << " sigma= " << sigma
	   << endl; 



    for (Sztildej=-Stildej;
	 Sztildej<=Stildej;Sztildej+=1.0){
      siteqnumsj[1]=Sztildej;
      int sitestatej=typej; // Ok if Stildej=0
      if ((typej==2)||(typej==3)) // If not...
	sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj)+1;
      if (disp) cout << " Sz~j = "<<Sztildej << " ssj =  "<< sitestatej << endl; 
      int i1=0;


      siteqnumsi[1]=Sztildej-sigma;
      int sitestatei=typei; // Ok if Stildei=0
      if ((typei==2)||(typei==3))
	sitestatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi)+1;
          
      auxCG[0]=CGordan(Scfi,Si-Sztildej+sigma,Stildei,Sztildej-sigma,Si,Si);
      auxCG[1]=CGordan(Scfj,Sj-Sztildej,Stildej,Sztildej,Sj,Sj);
      // Check this below!!
      auxBasis[i1]+=OneChS_fd_table((int)(sigma/0.5),sitestatej,sitestatei)*auxCG[0]*auxCG[1];

      if (disp) cout << " CG 1 = " << auxCG[0] << " CG 2 = " << auxCG[1] << " mat_el = " <<  OneCh_fd_table((int)(sigma/0.5),sitestatej,sitestatei) << endl;
      i1++;
    }
    // end loop in Sztildej
    // Check to see which of them is non-zero
    auxCG[0]=CGordan(Si,Si,0.5,sigma,Sj,Sj);
    if (disp) cout << " auxBasis = " << auxBasis[0] << " auxGC0 = " << auxCG[0] << endl;
    if ( !(dEqual(auxCG[0],0.0)) ) 
      fbasis=auxBasis[0]/auxCG[0];
    else{
      if ( !(dEqual(auxCG[1],0.0)) ) 
	fbasis=auxBasis[1]/auxCG[1];
      else fbasis=0.0;
    }
  }
  else fbasis=0.0;
  // end if stcfi==stcfj

  return(fbasis);

}

////////////////////////
////////////////////////

double OneChS_cd_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){

  double fbasis=0.0;

  double Si=pAbasis->GetQNumberFromSt(ist,0);
  int typei=pAbasis->iType[ist];
  int stcfi=pAbasis->StCameFrom[ist];

  double Sj=pAbasis->GetQNumberFromSt(jst,0);
  int typej=pAbasis->iType[jst];
  int stcfj=pAbasis->StCameFrom[jst];

  // Not one block per state!!
  double Stildei=pSingleSite->GetQNumberFromSt(typei,0);
  double Sztildei=pSingleSite->GetQNumberFromSt(typei,1);
                         
  double Stildej=pSingleSite->GetQNumberFromSt(typej,0);
  double Sztildej=pSingleSite->GetQNumberFromSt(typej,1);


  double siteqnumsi[]={Stildei,Sztildei};
  double siteqnumsj[]={Stildej,Sztildej};
  
  double Scfi=Si-Sztildei;
  double Scfj=Sj-Sztildej;
  
  double auxCG[4]={0.0,0.0,0.0,0.0};

  // Let me try this... Works!
  double sigma=Sj-Si;


  // Debug
  bool disp=false;
//    if (   ( (ist==5)&&(jst==9) ) )
//      disp=true;


  auxCG[3]=CGordan(Si,Si,0.5,sigma,Sj,Sj);

  // Check sigma=Si-Sj
  if (dEqual(auxCG[3],0.0)){
    sigma=Si-Sj;
    auxCG[3]=CGordan(Si,Si,0.5,sigma,Sj,Sj);
    if (dEqual(auxCG[3],0.0)) return(0.0);
  }


  for (double Szcfj=-Scfj;
       Szcfj<=Scfj;Szcfj+=1.0){
    double Szcfi=Szcfj-sigma;
    Sztildej=Sj-Szcfj;
    Sztildei=Si-Szcfi;
      
    siteqnumsj[1]=Sztildej;
    //int sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj); // Does not work
    int sitestatej=typej; // Ok if Stildej=0
    if ((typej==2)||(typej==3)) // If not...
      sitestatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj)+1;


    siteqnumsi[1]=Sztildei;
    int sitestatei=typei; // Ok if Stildei=0
    if ((typei==2)||(typei==3))
      sitestatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi)+1;


   if (disp) 
      cout << "I'm in: " 
	   << " ist=" << ist << " jst= " << jst 
	   << " typei=" << typei << " typej= " << typej
	   << " Si= " << Si << " Sj= " << Sj 
	   << " Scfi= " << Scfi << " Scfj= " << Scfj 
	   << " sigma= " << sigma << endl 
	   << " Szcfi=" << Szcfi << " Szcfj= " << Szcfj 
	   << " Sztildei=" << Sztildei << " Sztildej= " << Sztildej 
	   << " sitestatei=" << sitestatei << " sitestatei= " << sitestatej 
	   << endl; 


    auxCG[0]=CGordan(Scfi,Szcfi,Stildei,Sztildei,Si,Si);
    auxCG[1]=CGordan(Scfj,Szcfj,Stildej,Sztildej,Sj,Sj);
    auxCG[2]=CGordan(Scfi,Szcfi,0.5,sigma,Scfj,Szcfj);

    // Have to make this more clear but let's see:
    double FermiSign=1.0;
    // Jul 09: Does FermiSign enter here????
    // Yes, it does
    // OK, need to check this!!!
    //
    if ((typej==2)||(typej==3)) FermiSign=-1.0;
    // Fermi sign appears IF there is an
    // odd number of fermionic ops in site state
    // if (dEqual(Stildej,0.5)) FermiSign=-1.0;


    if (sitestatei==sitestatej)
      fbasis+=auxCG[0]*auxCG[1]*auxCG[2]*FermiSign;


    if (disp) cout << " CG 0= " << auxCG[0] 
		   << " CG 1= " << auxCG[1] 
		   << " CG 2= " << auxCG[2] 
		   << " CG 3= " << auxCG[3] 
		   << " Fermi Sign= " <<  FermiSign
		   << " fbasis= " <<  fbasis
		   << endl;

  }
  // end loop in Szcfj

  

  if (dNEqual(auxCG[3],0.0)) fbasis=fbasis/auxCG[3];
  else fbasis=0.0;

  return(fbasis);

}

////////////////////////
////////////////////////


////////////////////////
////////////////////////

///////////////////////////////////////////
/////////  <S' | HN | S>   ///////////
///////////////////////////////////////////



double OneChS_HNsc_MatEl(vector<double> Params,
                          CNRGbasisarray* pAbasis,
                          CNRGbasisarray* pSingleSite,
                          CNRGmatrix* MatArray,
                          int ist, int jst){

  // Calculates matrix elements for HN in the S basis
  //
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  eN            - Params[2]
  //  DeltaSC       - Params[3]
  //  <S'||fN||S>   - MatArray[0] (reduced)


  // Need to add the superconducting term HERE!!

  double MatEl=0.0;

  double Lambda=Params[0];
  double chi_N=Params[1];
  double eN=Params[2]; 
  double DeltaSC=Params[3]; 


  double Si=pAbasis->GetQNumberFromSt(ist,0);

  double Sj=pAbasis->GetQNumberFromSt(jst,0);

  if (dNEqual(Si,Sj)) return(0.0);

  if (ist==jst){
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
    // add eN term here
    int type_i=pAbasis->iType[ist];
    // Charge in basis |0>, |up dn>, |up>, |dn>
    double auxN[4]={0.0,2.0,1.0,1.0};
    MatEl+=eN*auxN[type_i];

  }else{
    
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    int stcf_i=pAbasis->StCameFrom[ist];
    int stcf_j=pAbasis->StCameFrom[jst];

    //Adding the SC term.
    if(stcf_i==stcf_j){
      MatEl-=DeltaSC*OneChS_fdupfddn_table(type_i,type_j);
      MatEl-=DeltaSC*OneChS_fdupfddn_table(type_j,type_i); // h.c.
    }else{
      // Matrix element from <istcf|| fN ||jstcf>
      double OldEl=MatArray[0].GetMatEl(stcf_i,stcf_j);
      int typep=type_i;
      int type=type_j;

      // if zero, try the h.c term
//       if (dEqual(OldEl,0.0)){
// 	OldEl=MatArray[0].GetMatEl(stcf_j,stcf_i);
// 	typep=type_j;
// 	type=type_i;
//       }
      // ALSO: if both basis are zero, try the h.c.
      if (   (dEqual(OldEl,0.0))||
	     (  (dEqual(OneChS_fd_table(1,type_i,type_j),0.0))&&
		(dEqual(OneChS_fd_table(-1,type_i,type_j),0.0)) )   ){
	OldEl=MatArray[0].GetMatEl(stcf_j,stcf_i);
	typep=type_j;
	type=type_i;
      }
      //

      // Get SingleSite QNumbers (watch out! NOT a single state per block!)

      double Stildep=pSingleSite->GetQNumberFromSt(typep,0);
      double Sztildep=pSingleSite->GetQNumberFromSt(typep,1);

      double Stilde=pSingleSite->GetQNumberFromSt(type,0);
      double Sztilde=pSingleSite->GetQNumberFromSt(type,1);

      double Scfp=Si-Sztildep;
      double Scf=Si-Sztilde;

      // Check Fermi Sign
      double FermiSign=1.0;
      if ((type==2)||(type==3)) FermiSign=-1.0;
      // Fermi sign appears IF there is an
      // odd number of fermionic ops in site state

      double CoefMatEl=0.0;
      double auxCG[]={0.0,0.0,0.0,0.0};
      double siteqnumsp[]={Stildep,Sztildep};
      double siteqnums[]={Stilde,Sztilde};
      // Sum in Sztilde
      for (Sztilde=-Stilde;Sztilde<=Stilde;Sztilde+=1.0){
	// Sum in sigma: Note that Sztildep=Sztilde+sigma
	siteqnums[1]=Sztilde;
	// Problem: more that one state per block in 1chS single site
	int sitestate=type; // Ok if Stilde=0
	if ((type==2)||(type==3)) // But not if Stilde=1/2...
	  sitestate=pSingleSite->GetBlockFromQNumbers(siteqnums)+1; 

	for (int sigma=-1;sigma<=1;sigma+=2){
	  double dSigma=0.5*(double)sigma;
	  double Szcf=Si-Sztilde;
	  double Szcfp=Szcf-dSigma;
	  Sztildep=Sztilde+dSigma;

	  siteqnumsp[1]=Sztildep;
	  // Test Find site state
	  int sitestatep=pSingleSite->GetBlockFromQNumbers(siteqnumsp)+1;
	  if (dEqual(Sztildep,0.0))
	    sitestatep=typep; // Ok if Sztildep=0

	  //cout << "Sz~ = " << Sztilde 
	  //   << " dSigma = " << dSigma << endl;
	
	  auxCG[0]=CGordan(Scf,Szcf,
			   Stilde,Sztilde,Si,Si);

	  auxCG[1]=CGordan(Scfp,Szcfp,
			   Stildep,Sztildep,Si,Si);

	  auxCG[2]=CGordan(Scfp,Szcfp,0.5,dSigma,Scf,Szcf);
 
	  auxCG[3]=OneChS_fd_table(sigma,sitestatep,sitestate);
	  //cout << "<sitestp|f|sitest> = " << auxCG[3] << endl;
	  CoefMatEl+=auxCG[0]*auxCG[1]*auxCG[2]*auxCG[3];

	  // Debugging
// 	  if ((jst==0)&&(ist==1)){
//        	  cout << " ist = " << ist 
//        	       << " jst = " << jst
//        	       << " stcf_i = " << stcf_i
//        	       << " stcf_j = " << stcf_j
//        	       << " type_i = " << type_i
//        	       << " type_j = " << type_j
// 	       << endl;
//        	  cout << " typep = " << typep
//        	       << " type = " << type 
//        	       << " OldEl = " << OldEl
//        	       << endl;
//        	  cout << " sigma, Si = " << sigma << " " << Si
//        	       << " sitestatep = " << sitestatep
//        	       << " sitestate = " << sitestate
// 	       << " fd_table = "
//        	       << OneChS_fd_table(sigma,sitestatep,sitestate) 
//        	       << endl;
// 	  cout << " Scf, Szcf = " << Scf << "  " << Szcf
// 	       << " Stilde,Sztilde = " << Stilde << "  " << Sztilde
// 	       << " Scfp, Szcfp = " << Scfp << "  " << Szcfp
// 	       << " Stildep,Sztildep = " << Stildep << "  " << Sztildep
// 	       << endl;
//       	  cout << " CG0 = " << auxCG[0]
// 	       << " CG1 = " << auxCG[1] 
// 	       << " CG2 = " << auxCG[2] 
// 	       << endl;
// 	  cout << " chi_N = " << chi_N 
// 	       << " OldEl = " << OldEl 
// 	       << " CoefMatEl = " << CoefMatEl
// 	       << " FermiSign = " << FermiSign 
// 	       << endl;
// 	  }
	  //////////////////////
				
	}
	// end loop in sigma
      }
      // end loop in Sztilde

      //auxEl*=CoefMatEl;
      MatEl+=chi_N*OldEl*CoefMatEl*FermiSign;
    }
    // end if stcf_i=stcf_j
  }
  // end if i=j

  return(MatEl);
}

////////////////////////
////////////////////////

////////////////////////////////////////
/////////  OneChSz symmetry   ///////////
////////////////////////////////////////

////////////////////////
////////////////////////

double OneChSz_fN_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst, int isigma){
  double fbasis=0.0;

  //double Szi=pAbasis->GetQNumberFromSt(ist,0);
  int typei=pAbasis->iType[ist];
  int stcfi=pAbasis->StCameFrom[ist];

  //double Szj=pAbasis->GetQNumberFromSt(jst,0);
  int typej=pAbasis->iType[jst];
  int stcfj=pAbasis->StCameFrom[jst];

  // Debug
  bool disp=false;
//   if (   ( (ist==5)&&(jst==1) )||
// 	 ( (ist==6)&&(jst==3) )||
// 	 ( (ist==7)&&(jst==4) )
// 	 )
//     disp=true;

  if ( (stcfi==stcfj) ){
    return(OneChS_fd_table(isigma,typej,typei));
  }else
    return(0.0);

}

////////////////////////

double OneChSz_fNup_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){


  return(OneChSz_fN_MatEl(pAbasis,pSingleSite,
			  ist, jst,1));


}

double OneChSz_fNdn_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){

  return(OneChSz_fN_MatEl(pAbasis,pSingleSite,
			  ist, jst,-1));

}

////////////////////////

////////////////////////
////////////////////////

double OneChSz_cd_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
		       int ist, int jst, int isigma){

  double Szi=pAbasis->GetQNumberFromSt(ist,0);
  int typei=pAbasis->iType[ist];
  int stcfi=pAbasis->StCameFrom[ist];

  double Szj=pAbasis->GetQNumberFromSt(jst,0);
  int typej=pAbasis->iType[jst];
  int stcfj=pAbasis->StCameFrom[jst];

  // Have to make this more clear but let's see:
  double FermiSign=1.0;
  if ((typej==2)||(typej==3)) FermiSign=-1.0;
  // Fermi sign appears IF there is an
  // odd number of fermionic ops in site state

  double dSigma=((double)isigma)*0.5;

  if (typei==typej){
    if ( dEqual(Szj,Szi+dSigma) ) return(FermiSign);
    else return(0.0);
  }
  else
    return(0.0);

}


////////////////////////
////////////////////////

////////////////////////

double OneChSz_cdup_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){
  return(OneChSz_cd_MatEl(pAbasis,pSingleSite,
			  ist, jst,1));
}

double OneChSz_cddn_MatEl(CNRGbasisarray *pAbasis, 
                        CNRGbasisarray *pSingleSite, 
                        int ist, int jst){
  return(OneChSz_cd_MatEl(pAbasis,pSingleSite,
			  ist, jst,-1));
}

////////////////////////

///////////////////////////////////////////
/////////  <Sz' | HN | Sz>   ///////////
///////////////////////////////////////////



double OneChSz_HNsc_MatEl(vector<double> Params,
                          CNRGbasisarray* pAbasis,
                          CNRGbasisarray* pSingleSite,
                          CNRGmatrix* MatArray,
                          int ist, int jst){

  // Calculates matrix elements for HN in the S basis
  //
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  eN            - Params[2]
  //  DeltaSC       - Params[3]
  //  <Sz-1/2|fN_up|Sz>   - MatArray[0] (not reduced)
  //  <Sz+1/2|fN_dn|Sz>   - MatArray[1] (not reduced)


  // Need to add the superconducting term HERE!!

  double MatEl=0.0;

  double Lambda=Params[0];
  double chi_N=Params[1];
  double eN=Params[2]; 
  double DeltaSC=Params[3]; 


  double Szi=pAbasis->GetQNumberFromSt(ist,0);

  double Szj=pAbasis->GetQNumberFromSt(jst,0);

  if (dNEqual(Szi,Szj)) return(0.0);

  if (ist==jst){
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
    // add eN term here
    int type_i=pAbasis->iType[ist];
    // Charge in basis |0>, |up dn>, |up>, |dn>
    double auxN[4]={0.0,2.0,1.0,1.0};
    MatEl+=eN*auxN[type_i];

  }else{
    
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    int stcf_i=pAbasis->StCameFrom[ist];
    int stcf_j=pAbasis->StCameFrom[jst];

    //Adding the SC term. Works for OneChSz as well
    if(stcf_i==stcf_j){
      MatEl-=DeltaSC*OneChS_fdupfddn_table(type_i,type_j);
      MatEl-=DeltaSC*OneChS_fdupfddn_table(type_j,type_i); // h.c.
    }else{
      // Matrix element from <istcf| fN_up |jstcf> and <istcf| fN_dn |jstcf>
      // Sum over spins
      int imat=0;
      for (int isigma=1;isigma>=-1;isigma-=2){
	double OldEl=MatArray[imat].GetMatEl(stcf_i,stcf_j);
	int typep=type_i;
	int type=type_j;

	// if zero, try the h.c term, meaning f+N f_N+1 in original case
 	// ALSO try the h.c. if the basis element is zero.
	// This because we can have a non-zero OldEl but not the correct one
	if (   (dEqual(OldEl,0.0))||
	       (dEqual(OneChS_fd_table(isigma,type_i,type_j),0.0))   ){
	  OldEl=MatArray[imat].GetMatEl(stcf_j,stcf_i);
	  typep=type_j;
	  type=type_i;
	}
	//

	// Check Fermi Sign
	double FermiSign=1.0;
	if ((type==2)||(type==3)) FermiSign=-1.0;
	// Fermi sign appears IF there is an
	// odd number of fermionic ops in site state

       	// fbasis(typep,type) should be ok in the Sz case!
	MatEl+=FermiSign*OldEl*OneChS_fd_table(isigma,typep,type);

	// Debugging
// 	if (   ( ((jst==29)&&(ist==20))||
// 		 ((jst==44)&&(ist==35)) )&&(pAbasis->Nshell<=1)   ){
// 	  cout << " ist = " << ist 
// 	       << " jst = " << jst
// 	       << " stcf_i = " << stcf_i
// 	       << " stcf_j = " << stcf_j
// 	       << " type_i = " << type_i
// 	       << " type_j = " << type_j
// 	       << endl;
// 	  cout << " typep = " << typep
// 	       << " type = " << type 
// 	       << " isigma = " << isigma
// 	       << " imat = " << imat
// 	       << endl;
// 	  cout << " fd_table = "
// 	       << OneChS_fd_table(isigma,typep,type) 
// 	       << " OldEl = " << OldEl 
// 	       << " FermiSign = " << FermiSign 
// 	       << " MatEl = " << MatEl 
// 	       << endl;
// 	  cout << " chi_N*MatEl = " << chi_N*MatEl 
// 	       << endl;

// 	}
	//////////////////////
	imat++;
      }
      // End sum over isigma

      MatEl*=chi_N;
    }
    // end if stcf_i=stcf_j
  }
  // end if i=j

  return(MatEl);
}

////////////////////////
////////////////////////


/////////////////////////////////////////////
/////////  OneChPupPdn symmetry   ///////////
/////////////////////////////////////////////

////////////////////////

double OneChPupPdn_cd_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst, int isigma){

  double Pupi=pAbasis->GetQNumberFromSt(ist,0);
  double Pdni=pAbasis->GetQNumberFromSt(ist,1);
  int stcfi=pAbasis->StCameFrom[ist];
  int typei=pAbasis->iType[ist];

  double Pupj=pAbasis->GetQNumberFromSt(jst,0);
  double Pdnj=pAbasis->GetQNumberFromSt(jst,1);
  int stcfj=pAbasis->StCameFrom[jst];
  int typej=pAbasis->iType[jst];


//   cout << " isigma = " << isigma
//        << " stcfi = " << stcfi
//        << " stcfj = " << stcfj
//        << " typei = " << typei
//        << " typej = " << typej
//        << endl
//        << " <ti|f+|tj> = " << 
//     OneCh_fd_table(isigma,typej,typei)
//        << endl;

  double dSigma=((double)isigma)*0.5;

  // Have to make this more clear but let's see:
  // if |typej> = |up> or |dn>, we sign=-1
  double FermiSign=1.0;
  if ((typej==1)||(typej==2)) FermiSign=-1.0;

  if (typei==typej){
    if (   (  (isigma==1)&&(dEqual(Pupj,-1.0*Pupi))&&(dEqual(Pdni,Pdnj)) )||
	   (  (isigma==-1)&&(dEqual(Pupj,Pupi))&&(dEqual(Pdni,-1.0*Pdnj)) )   )
      return(FermiSign);
    else return(0.0);
  }
  else
    return(0.0);

}

//////
double OneChPupPdn_cdup_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChPupPdn_cd_MatEl(pAbasis,pSingleSite,
			 ist, jst,1));
}

double OneChPupPdn_cddn_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChPupPdn_cd_MatEl(pAbasis,pSingleSite,
			 ist, jst,-1));
}

complex<double> OneChPupPdn_cdup_MatElCplx(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  double aux=OneChPupPdn_cd_MatEl(pAbasis,pSingleSite,
				  ist, jst,1);
  return(complex<double>(aux,0.0));
}

complex<double> OneChPupPdn_cddn_MatElCplx(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){

  double aux=OneChPupPdn_cd_MatEl(pAbasis,pSingleSite,
				  ist, jst,-1);
  return(complex<double>(aux,0.0));
}



//////




double OneChPupPdn_fN_MatEl(CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    int ist, int jst, int isigma){

  int stcfi=pAbasis->StCameFrom[ist];
  int typei=pAbasis->iType[ist];

  int stcfj=pAbasis->StCameFrom[jst];
  int typej=pAbasis->iType[jst];

  if ( (stcfi==stcfj) ){
      return(OneCh_fd_table(isigma,typej,typei));
  }else return(0.0);



}

//////

double OneChPupPdn_fNup_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChPupPdn_fN_MatEl(pAbasis,pSingleSite,
			 ist, jst,1));
}

double OneChPupPdn_fNdn_MatEl(CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       int ist, int jst){
  return(OneChPupPdn_fN_MatEl(pAbasis,pSingleSite,
			 ist, jst,-1));
}

complex<double> OneChPupPdn_fNup_MatElCplx(CNRGbasisarray* pAbasis, 
					   CNRGbasisarray* pSingleSite,
					   int ist, int jst){
  double aux=OneChPupPdn_fN_MatEl(pAbasis,pSingleSite,
				  ist, jst,1);
  return(complex<double>(aux,0.0));
}

complex<double> OneChPupPdn_fNdn_MatElCplx(CNRGbasisarray* pAbasis, 
					   CNRGbasisarray* pSingleSite,
					   int ist, int jst){
  double aux=OneChPupPdn_fN_MatEl(pAbasis,pSingleSite,
				  ist, jst,-1);
  return(complex<double>(aux,0.0));
}
//////

///
/// <Pup' Pdn'|HN|Pup Pdn> :
///
/////////////////////////////////
/////////////////////////////////
complex<double> OneChPupPdn_HN_MatElCplx(vector<double> Params,
			    CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    CNRGmatrix* MatArray,
			    int ist, int jst){

  // June 20/2017
  // Calculates matrix elements for HN in the Pup Pdn basis
  //
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  eN            - Params[2] (added)
  //  <-Pup Pdn|fN_up|Pup Pdn> - MatArray[0] (not reduced)
  //  <Pup -Pdn|fN_dn|Pup Pdn> - MatArray[1] (not reduced)

  complex<double> cMatEl=(0.0,0.0);

  double Lambda=Params[0];
  double chi_N=Params[1];
  double eN=Params[2];


  double Pupi=pAbasis->GetQNumberFromSt(ist,0);
  double Pdni=pAbasis->GetQNumberFromSt(ist,1);

  double Pupj=pAbasis->GetQNumberFromSt(jst,0);
  double Pdnj=pAbasis->GetQNumberFromSt(jst,1);



  if ( (dNEqual(Pupi,Pupj))||(dNEqual(Pdni,Pdnj)) ) return(0.0);

  if (ist==jst){
    double MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
    // add eN term here
    int type_i=pAbasis->iType[ist];
    double Puptilde=pSingleSite->GetQNumberFromSt(type_i,0);
    double Pdntilde=pSingleSite->GetQNumberFromSt(type_i,1);
    // Calculate Qtilde
    double Qtilde=0.0;
    if ( (dEqual(Puptilde,1.0))&&(dEqual(Pdntilde,1.0)) ) Qtilde=-1.0;  // |0>
    if ( (dEqual(Puptilde,-1.0))&&(dEqual(Pdntilde,-1.0)) ) Qtilde=1.0; // |up dn>
    MatEl+=eN*(Qtilde+1.0);
    cMatEl=complex<double>(MatEl,0.0);
  }else{
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    // Sum over spins
    int imat=0;
    for (int isigma=1;isigma>=-1;isigma-=2){
      int typep=type_i;
      int type=type_j;

      int stcf_i=pAbasis->StCameFrom[ist];
      int stcf_j=pAbasis->StCameFrom[jst];

      // Matrix element from <istcf| fN_sigma | jstcf>
      complex<double> cOldEl=MatArray[imat].cGetMatEl(stcf_i,stcf_j);

      complex<double> cOldEl_old=cOldEl;

      // if zero, try the h.c term, meaning f+N f_N+1 in original case
      // Ops, it can be non-zero but the wrong one!!
//       if ( (dEqualPrec(cOldEl.real(),0.0,1e-10))&&(dEqualPrec(cOldEl.imag(),0.0,1e-10)) ){
// 	// debug
// 	cOldEl=MatArray[imat].cGetMatEl(stcf_j,stcf_i);
// 	typep=type_j;
// 	type=type_i;
//       }
      // end if
      double OpTable=OneCh_fd_table(isigma,typep,type);
      if ( dEqual(OpTable,0.0) ){
	// debug
	cOldEl=MatArray[imat].cGetMatEl(stcf_j,stcf_i);
	typep=type_j;
	type=type_i;
      }
      // end if


      // PROBLEM: There might be a non-zero Hermitian term!!

      // Get SingleSite QNumbers

      double Puptilde=pSingleSite->GetQNumberFromSt(type,0);
      double Pdntilde=pSingleSite->GetQNumberFromSt(type,1);

      // Fermi Sign
      double FermiSign=1.0;
      // 	if (dEqual(Qtilde,0.0)) FermiSign=-1.0;
      // |tilde> is |up>=|Pup=-1 Pdn=1> or |dn>=|Pup=1 Pdn=-1>
      if (  ( (dEqual(Puptilde,-1.0))&&(dEqual(Pdntilde,1.0)) )||
	    ( (dEqual(Puptilde,1.0))&&(dEqual(Pdntilde,-1.0)) )
	   )  
	FermiSign=-1.0;

      cMatEl+=FermiSign*cOldEl*OneCh_fd_table(isigma,typep,type);

      // Debug 
//       if (  (pAbasis->Nshell==1)&&
// 	    ( ( (ist==0)&&(jst>=45)&&(jst<=45) )
// 	      )  ){
// 	cout << " imat = " << imat
// 	     << " ist= " << ist
// 	     << " jst= " << jst
// 	     << " stcfi= " << stcf_i
// 	     << " stcfj= " << stcf_j
// 	     << " typei= " << type_i
// 	     << " typej= " << type_j 
//  	     << " typep= " << typep
// 	     << " type= " << type 
// 	     << endl;
// 	if (imat==0)
// 	  cout << " <stcfi|f_up|scfj>= " << MatArray[imat].cGetMatEl(stcf_i,stcf_j)
// 	       << " <stcfj|f_up|scfi>= " << MatArray[imat].cGetMatEl(stcf_j,stcf_i)
// 	       << " cOldEl_old= " << cOldEl_old
// 	       << " cOldEl= " << cOldEl 
// 	       << endl;
// 	else
// 	  cout << " <stcfi|f_dn|scfj>= " << MatArray[imat].cGetMatEl(stcf_i,stcf_j)
// 	       << " <stcfj|f_dn|scfi>= " << MatArray[imat].cGetMatEl(stcf_j,stcf_i)
// 	       << " cOldEl_old= " << cOldEl_old
// 	       << " cOldEl= " << cOldEl 
// 	       << endl;
// 	cout << " Fermi sign= " << FermiSign
// 	     << " |Re(cOldEl_old)| < 1e-10? " << dEqualPrec(cOldEl_old.real(),0.0,1e-10)
// 	     << " |Im(cOldEl_old)| < 1e-10? " << dEqualPrec(cOldEl_old.imag(),0.0,1e-10)
// 	     << " |Re(cOldEl_old)| < 1e-15? " << dEqualPrec(cOldEl_old.real(),0.0,1e-15)
// 	     << " |Im(cOldEl_old)| < 1e-15? " << dEqualPrec(cOldEl_old.imag(),0.0,1e-15)
// 	     << " (f_dn)^d_table= " << OneCh_fd_table(isigma,typep,type)
// 	     << " Mat El/chi_N= " << cMatEl
// 	     << " Mat El= " << cMatEl*chi_N 
// 	     << endl;
//       }
      ///

      imat++;
    }
    // End sum over isigma
      
    cMatEl*=chi_N;
  }
  // end if i=j


  return(cMatEl);
}

/////////////////////////////////
/////////////////////////////////
double OneChPupPdn_HN_MatEl(vector<double> Params,
			    CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite,
			    CNRGmatrix* MatArray,
			    int ist, int jst){

  // July 21/2017
  // Calculates matrix elements for HN in the Pup Pdn basis
  // REAL for chain calculations!
  // Parameters: 
  //  Lambda        - Params[0]
  //  chi_N         - Params[1]
  //  eN            - Params[2] (added)
  //  <-Pup Pdn|fN_up|Pup Pdn> - MatArray[0] (not reduced)
  //  <Pup -Pdn|fN_dn|Pup Pdn> - MatArray[1] (not reduced)

  double MatEl=0.0;

  double Lambda=Params[0];
  double chi_N=Params[1];
  double eN=Params[2];


  double Pupi=pAbasis->GetQNumberFromSt(ist,0);
  double Pdni=pAbasis->GetQNumberFromSt(ist,1);

  double Pupj=pAbasis->GetQNumberFromSt(jst,0);
  double Pdnj=pAbasis->GetQNumberFromSt(jst,1);



  if ( (dNEqual(Pupi,Pupj))||(dNEqual(Pdni,Pdnj)) ) return(0.0);

  if (ist==jst){
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
    // add eN term here
    int type_i=pAbasis->iType[ist];
    double Puptilde=pSingleSite->GetQNumberFromSt(type_i,0);
    double Pdntilde=pSingleSite->GetQNumberFromSt(type_i,1);
    // Calculate Qtilde
    double Qtilde=0.0;
    if ( (dEqual(Puptilde,1.0))&&(dEqual(Pdntilde,1.0)) ) Qtilde=-1.0;  // |0>
    if ( (dEqual(Puptilde,-1.0))&&(dEqual(Pdntilde,-1.0)) ) Qtilde=1.0; // |up dn>
    MatEl+=eN*(Qtilde+1.0);

  }else{
    int type_i=pAbasis->iType[ist];
    int type_j=pAbasis->iType[jst];

    // Sum over spins
    int imat=0;
    for (int isigma=1;isigma>=-1;isigma-=2){
      int typep=type_i;
      int type=type_j;

      int stcf_i=pAbasis->StCameFrom[ist];
      int stcf_j=pAbasis->StCameFrom[jst];

      // Matrix element from <istcf| fN_sigma | jstcf>
      double OldEl=MatArray[imat].GetMatEl(stcf_i,stcf_j);

      double OldEl_old=OldEl;

      // if zero, try the h.c term, meaning f+N f_N+1 in original case
      // Ops, it can be non-zero but the wrong one!!
//       if ( (dEqualPrec(cOldEl.real(),0.0,1e-10))&&(dEqualPrec(cOldEl.imag(),0.0,1e-10)) ){
// 	// debug
// 	cOldEl=MatArray[imat].cGetMatEl(stcf_j,stcf_i);
// 	typep=type_j;
// 	type=type_i;
//       }
      // end if
      double OpTable=OneCh_fd_table(isigma,typep,type);
      if ( dEqual(OpTable,0.0) ){
	// debug
	OldEl=MatArray[imat].GetMatEl(stcf_j,stcf_i);
	typep=type_j;
	type=type_i;
      }
      // end if

      // Get SingleSite QNumbers

      double Puptilde=pSingleSite->GetQNumberFromSt(type,0);
      double Pdntilde=pSingleSite->GetQNumberFromSt(type,1);

      // Fermi Sign
      double FermiSign=1.0;
      // 	if (dEqual(Qtilde,0.0)) FermiSign=-1.0;
      // |tilde> is |up>=|Pup=-1 Pdn=1> or |dn>=|Pup=1 Pdn=-1>
      if (  ( (dEqual(Puptilde,-1.0))&&(dEqual(Pdntilde,1.0)) )||
	    ( (dEqual(Puptilde,1.0))&&(dEqual(Pdntilde,-1.0)) )
	   )  
	FermiSign=-1.0;

      MatEl+=FermiSign*OldEl*OneCh_fd_table(isigma,typep,type);

      // Debug 
//       if (  (pAbasis->Nshell==0)&&
// 	    ( ( (ist==5)&&(jst>=10)&&(jst<=10) )
// 	      )  ){
// 	cout << " imat = " << imat
// 	     << " ist= " << ist
// 	     << " jst= " << jst
// 	     << " stcfi= " << stcf_i
// 	     << " stcfj= " << stcf_j
// 	     << " typei= " << type_i
// 	     << " typej= " << type_j 
//  	     << " typep= " << typep
// 	     << " type= " << type 
// 	     << endl;
// 	if (imat==0)
// 	  cout << " <stcfi|f_up|scfj>= " << MatArray[imat].cGetMatEl(stcf_i,stcf_j)
// 	       << " <stcfj|f_up|scfi>= " << MatArray[imat].cGetMatEl(stcf_j,stcf_i)
// 	       << " cOldEl_old= " << cOldEl_old
// 	       << " cOldEl= " << cOldEl 
// 	       << endl;
// 	else
// 	  cout << " <stcfi|f_dn|scfj>= " << MatArray[imat].cGetMatEl(stcf_i,stcf_j)
// 	       << " <stcfj|f_dn|scfi>= " << MatArray[imat].cGetMatEl(stcf_j,stcf_i)
// 	       << " cOldEl_old= " << cOldEl_old
// 	       << " cOldEl= " << cOldEl 
// 	       << endl;
// 	cout << " Fermi sign= " << FermiSign
// 	     << " |Re(cOldEl_old)| < 1e-10? " << dEqualPrec(cOldEl_old.real(),0.0,1e-10)
// 	     << " |Im(cOldEl_old)| < 1e-10? " << dEqualPrec(cOldEl_old.imag(),0.0,1e-10)
// 	     << " |Re(cOldEl_old)| < 1e-15? " << dEqualPrec(cOldEl_old.real(),0.0,1e-15)
// 	     << " |Im(cOldEl_old)| < 1e-15? " << dEqualPrec(cOldEl_old.imag(),0.0,1e-15)
// 	     << " (f_dn)^d_table= " << OneCh_fd_table(isigma,typep,type)
// 	     << " Mat El/chi_N= " << cMatEl
// 	     << " Mat El= " << cMatEl*chi_N 
// 	     << endl;
//       }
//       ///

      imat++;
    }
    // End sum over isigma
      
    MatEl*=chi_N;
  }
  // end if i=j


  return(MatEl);
}
////////////////////////
////////////////////////
