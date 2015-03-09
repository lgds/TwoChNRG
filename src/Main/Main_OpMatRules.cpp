
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "TwoChQS.hpp"

// Note (Oct 2013)
// This file should contains only drivers for setting up
// specific matrices such as Hamiltonians for N=-1 or N=0
// and NOT generic NRG matrices (fN, cd, etc.)
//


////////////////////////////////////////////
///////////     MatElChecks    /////////////
////////////////////////////////////////////

// Moved to OpMatRules.cpp

// bool OneChQ_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

//   double Qi=pAeigCut->GetQNumber(iblock1,0);
//   double Qj=pAeigCut->GetQNumber(iblock2,0);
  
//   if ( dEqual(Qj,(Qi+1.0)) ) 
//     return(true);
//   else
//     return(false);

// }


////////////////////////////////////////////

// Moved to OpMatRules.cpp
//
// bool OneChQSz_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2){

//   double Qi=pAeigCut->GetQNumber(iblock1,0);
//   double Szi=pAeigCut->GetQNumber(iblock1,1);

//   double Qj=pAeigCut->GetQNumber(iblock2,0);
//   double Szj=pAeigCut->GetQNumber(iblock2,1);

  
//   if (  (dEqual(Qj,(Qi+1.0)))&&
// 	( (dEqual(Szj,(Szi+0.5)))||(dEqual(Szj,(Szi-0.5))) )
// 	) 
//     return(true);
//   else
//     return(false);

// }


////////////////////////////////////////////
////////////////////////////////////////////
///////////   Mat ElCalculations   /////////
////////////////////////////////////////////
////////////////////////////////////////////


//////////////////////////////
///////  Operators   /////////
//////////////////////////////



///////////////////////////
///
/// Basis for f_dot1up, f_dot1dn, f_dot2up, f_dot2dn
///

// Use TwoChQS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2)
// for check

double TwoChQSNoSz_fm1_dot1up_MatEl(CNRGbasisarray* pAbasis, 
				  CNRGbasisarray* pSingleSite,
				  int ist, int jst)
{
  // Basis has the SAME structure as TwoChQS_SingleSite. 
  /// Let's use that
  // channel 1 <-> dot 1
  return(TwoChQS_fd_table(1,1,jst,ist));

}

///
double TwoChQSNoSz_fm1_dot1dn_MatEl(CNRGbasisarray* pAbasis, 
				  CNRGbasisarray* pSingleSite,
				  int ist, int jst)
  {return(TwoChQS_fd_table(1,-1,jst,ist));}

///
double TwoChQSNoSz_fm1_dot2up_MatEl(CNRGbasisarray* pAbasis, 
				  CNRGbasisarray* pSingleSite,
				  int ist, int jst)
  {return(TwoChQS_fd_table(2,1,jst,ist));}

///
double TwoChQSNoSz_fm1_dot2dn_MatEl(CNRGbasisarray* pAbasis, 
				  CNRGbasisarray* pSingleSite,
				  int ist, int jst)
  {return(TwoChQS_fd_table(2,-1,jst,ist));}


//////////////////////////






/////////////////////////////////
///////  Hamiltonians   /////////
/////////////////////////////////




/////////////////////////////////
/////////////////////////////////
///
/// <Q'S' | H_SMM |Q S> :
///


double TwoChQSNoSz_Hm1SMM_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
		       int ist, int jst)
{

  // August 08
  // Calculates matrix elements for H0_SMM 
  // Sz enters in pAbasis->iDegen[i]/10

  double MatEl=0.0;


  // Get Params
  double Lambda=Params[0];
  double HalfLambdaFactor=Params[1];


  double U=Params[2];
  double ed=Params[3];
  double U2=Params[4];
  double ed2=Params[5];

  double J12=Params[6];

  double BmagPar=Params[7];
  double BmagPerp=Params[8];
 
  double Dz=Params[9];
  double B2=Params[10];

//   cout << " U1 = " << U
//        << " ed1 = " << ed
//        << " U2 = " << U2
//        << " ed2 = " << ed2
//        << " J12 = " << J12
//        << " BmagPar = " << BmagPar
//        << " BmagPerp = " << BmagPerp
//        << " Dz = " << Dz
//        << " B2 = " << B2
//        << endl;
  // Dot energies
  double en[3][2]={{0.5*U,0.5*U2},
		   {ed+0.5*U,ed2+0.5*U2},
		   {2.0*ed+1.5*U,2.0*ed2+1.5*U2}};

  // divide by Lambda*HalfLambdaFactor // check to see if is sqrt(Lambda)
  for (int ii=0;ii<3;ii++)
    for (int jj=0;jj<2;jj++)
      en[ii][jj]/=(Lambda*HalfLambdaFactor);


  // Hd1 in the basis

  double Hd1[16]={en[0][0],en[0][0],en[1][0],en[0][0],en[1][0],en[1][0],en[0][0],en[2][0],
		  en[1][0],en[1][0],en[1][0],en[1][0],en[2][0],en[1][0],en[2][0],en[2][0]};

  double Hd2[16]={en[0][1],en[1][1],en[0][1],en[1][1],en[0][1],en[1][1],en[2][1],en[0][1],
		  en[1][1],en[1][1],en[1][1],en[2][1],en[1][1],en[2][1],en[1][1],en[2][1]};

  double HJ[16]={0.0,0.0,0.0,0.0,0.0,-3.0*J12/4.0,0.0,0.0,
		 J12/4.0,J12/4.0,J12/4.0,0.0,0.0,0.0,0.0,0.0};

//   cout << " Ed1 = "; 
//   for (int ii=0;ii<16;ii++) cout << Hd1[ii] << "  ";
//   cout << endl;

//   cout << " Ed2 = "; 
//   for (int ii=0;ii<16;ii++) cout << Hd2[ii] << "  ";
//   cout << endl;

//   cout << " EJ = "; 
//   for (int ii=0;ii<16;ii++) cout << HJ[ii] << "  ";
//   cout << endl;


  // Get Szi, Szj


  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Si=pAbasis->GetQNumberFromSt(ist,1);


  int type_i=pAbasis->iType[ist];
  int type_j=pAbasis->iType[jst];

  int iSzi=pAbasis->iDegen[ist];
  int iSzj=pAbasis->iDegen[jst];

  double Szi=(double)(iSzi/10);
  double Szj=(double)(iSzj/10);

  
//   cout << " iSzi = " << iSzi
//        << " iSzj = " << iSzj 
//        << " Si = " << Si 
//        << " Szi = " << Szi 
//        << " Szj = " << Szj
//        << endl; 

  if (ist==jst)
    {
      MatEl=Hd1[ist]+Hd2[ist]+HJ[ist]+Dz*Szi*Szi+BmagPar*Szi;
    }
  else
    {
      MatEl=0.0;
//       cout << " Si*(Si+1.0)-Szj*(Szj+1.0) = " << Si*(Si+1.0)-Szj*(Szj+1.0) << endl;
//       cout << " Si*(Si+1.0)-Szj*(Szj-1.0) = " << Si*(Si+1.0)-Szj*(Szj-1.0) << endl;

      if (abs(iSzi-iSzj)==10) MatEl=BmagPerp;
      if (iSzi==(iSzj+20)) 
	MatEl=B2*sqrt(Si*(Si+1.0)-Szj*(Szj+1.0))*sqrt(Si*(Si+1.0)-(Szj+2.0)*(Szj+1.0));
      if (iSzi==(iSzj-20)) 
	MatEl=B2*sqrt(Si*(Si+1.0)-Szj*(Szj-1.0))*sqrt(Si*(Si+1.0)-(Szj-2.0)*(Szj-1.0));
    }



  return(MatEl);


}

/////////////////////////////////
/////////////////////////////////
///
/// <Q' | H_SMM + site | Q> :
///


double OneChQ_H0SMM_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* fHm1Array,
			  int ist, int jst){

  double MatEl=0.0;

  double Lambda=Params[0];

  double chi_N[2]={Params[11],Params[12]};


  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Qj=pAbasis->GetQNumberFromSt(jst,0);

  int stcf_i=pAbasis->StCameFrom[ist];
  int stcf_j=pAbasis->StCameFrom[jst];


  if (dNEqual(Qi,Qj)) return(0.0);

  if (ist==jst)
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist]; // Why Lambda???
  else
    {
      int type_i=pAbasis->iType[ist];
      int type_j=pAbasis->iType[jst];

      // Sum over dots

      int imat=0;
      for (int idot=1;idot<=2;idot++)
	{
	  // Sum over spins
	  for (int isigma=1;isigma>=-1;isigma-=2)
	    {
	      int typep=type_i;
	      int type=type_j;

	      // Matrix element from <istcf| fN_sigma | jstcf>
	      double OldEl=fHm1Array[imat].GetMatEl(stcf_i,stcf_j);

	      // if zero, try the h.c term, meaning f+N f_N+1 in original case
	      if (dEqual(OldEl,0.0))
		{
		  OldEl=fHm1Array[imat].GetMatEl(stcf_j,stcf_i);
		  typep=type_j;
		  type=type_i;
		}

	      // Get SingleSite QNumbers

	      double Qtilde=pSingleSite->GetQNumberFromSt(type,0);

	      // Fermi Sign
	      double FermiSign=1.0;
	      if (dEqual(Qtilde,0.0)) FermiSign=-1.0;

// 	      	  cout << " ist = " << ist 
// 	      	       << " jst = " << jst
// 	      	       << " stcf_i = " << stcf_i
// 	      	       << " stcf_j = " << stcf_j
// 	      	       << " type_i = " << type_i
// 	      	       << " type_j = " << type_j
// 	      	       << " sigma = " << isigma << endl;
// 	      	  cout << " typep = " << typep
// 	      	       << " type = " << type 
// 	      	       << " OldEl = " << OldEl
// 	      	       << endl;
// 	      	  cout << " fd_table = "
// 	      	       << OneCh_fd_table(isigma,typep,type) 
// 	      	       << endl;
// 	      	  cout << " Qtilde = " << Qtilde << endl;
// 	      	  cout << " FermiSign = " << FermiSign << endl;
// 	      	  cout << " chi_N = " << chi_N[0] << " " << chi_N[1] << endl;



	      MatEl+=chi_N[idot-1]*FermiSign*OldEl*OneCh_fd_table(isigma,typep,type);
	      imat++;

	    }
	  // End sum over isigma     
	}
      // End sum over dots

    }
  // end if i=j





  return(MatEl);

}

/////////////////////////////////
/////////////////////////////////
///
/// <Q S P (m) | He_ph + site | Q S P (m')> :
///


double TwoChQSP_H0ph_MatEl(vector<double> Params,
		       CNRGbasisarray* pAbasis, 
		       CNRGbasisarray* pSingleSite,
		       CNRGmatrix* MatArray,
		       int ist, int jst)
{

  double MatEl=0.0;
  // <mi| a + a^+ |mj>

  double chi_N[2]={Params[0],Params[1]};
  double w0=Params[2];
  double lambda=Params[3];
  double alpha=Params[4];

////// Types

  int stcf_i=pAbasis->StCameFrom[ist];
  int stcf_j=pAbasis->StCameFrom[jst];

  int type_i=pAbasis->iType[ist];
  int type_j=pAbasis->iType[jst];

  int mi=pAbasis->iDegen[ist];
  int mj=pAbasis->iDegen[jst];

  double delta_mimj=(mi==mj)?1.0:0.0;
  
  // Phonon Mat el

  double phonon_matel=Calc_apad((double)mi, (double)mj);

  // MatArray[0]: Hm1 in the OLD basis
  double Hm1MatEl=MatArray[0].GetMatEl(stcf_i,stcf_j);

  // MatArray[1]: reduced <||d||> elements in the OLD basis
  CNRGmatrix Mat_d[1]; // Need for d_fd_MatEl
  Mat_d[0]=MatArray[1];


  // Off-diagonal terms
  // TwoChQS_d_fd should work for QSP symmetry too.
  vector<double> ParamsCh;
  double Hd_ch[2]={0.0,0.0};
  ParamsCh.push_back(1.0); //ch=1



  Hd_ch[0]=delta_mimj*chi_N[0]*TwoChQS_d_fd_MatEl(ParamsCh,pAbasis,
						  pSingleSite,Mat_d,
						  ist,jst);


  ParamsCh[0]=2.0; //ch=2
  Hd_ch[1]=alpha*chi_N[1]*phonon_matel*TwoChQS_d_fd_MatEl(ParamsCh,
							  pAbasis,
							  pSingleSite,
							  Mat_d,
							  ist,jst);
//   cout << " alpha = " << alpha
//        << " chi1  = " << chi_N[1]
//        << " ph_matel = " << phonon_matel
//        << " d_fd  = " << TwoChQS_d_fd_MatEl(ParamsCh,pAbasis,pSingleSite,
// 					    Mat_d,ist,jst)
//        << endl;


  if (ist==jst){
    MatEl=delta_mimj*Hm1MatEl; // Hubbard
    MatEl+=delta_mimj*mj*w0; // Phonon
  }
  else{
    MatEl=Hd_ch[0]+Hd_ch[1];
  }
  // end if diagonal


  return(MatEl);
}
///////////////////


/////////////////////////////////
/////////////////////////////////
///
/// <Q' S' | H_DQD + site | Q S> :
///

double OneChQS_H0DQD_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst){

  // Sep 09
  // Calculates matrix elements for H0 (DQD+site) in the QS basis.
  //
  // Parameters: 
  //  chi1        - Params[0]
  //  chi2        - Params[1]
  //  <Qm1 S +/-1/2|| f1 |Q S> - MatArray[0] (reduced)
  //  <Qm1 S+/-1/2 || f2 ||Q S> - MatArray[1] (reduced)
  // Abasis       - DQD + site
  // pSingleSite  - OneChQS single site

  double MatEl=0.0;

  double chi_N[2]={Params[0],Params[1]};

  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Si=pAbasis->GetQNumberFromSt(ist,1);

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Sj=pAbasis->GetQNumberFromSt(jst,1);

  if ( (dNEqual(Qi,Qj))||(dNEqual(Si,Sj)) ) return(0.0);

  //////////**///////////////


  double OldEl[2]={0.0,0.0};
  
  // Includes the possibility of partiy/Other QNs
  int NqnsSite=pSingleSite->NQNumbers;
  double* siteqnumsp=new double [NqnsSite];
  double* siteqnums=new double [NqnsSite];

  if (ist==jst){ // Diagonal terms
  // Already multiplied U1~ and ed~ by sqrt(Lambda)
  // (i.e. divide by sqrt(Lambda) instead of Lambda)
  // Thus, no lmabda^{1/2} factor in setting up H0.
    MatEl=pAbasis->dEn[ist];
  }
  else{
    int typei=pAbasis->iType[ist];
    int stcfi=pAbasis->StCameFrom[ist];
    
    int typej=pAbasis->iType[jst];
    int stcfj=pAbasis->StCameFrom[jst];

    // Get c1 matrix element from MatArray[0]
    // Get c2 matrix element from MatArray[1]
    // Get f0 matrix element using pSingleSite


    // Sum in dots

    for (int idot=1;idot<=2;idot++){
      OldEl[idot-1]=MatArray[idot-1].GetMatEl(stcfi,stcfj); // This is a reduced MatEl

      int typep=typei;
      int type=typej;

      if (dEqual(fabs(OldEl[idot-1]),0.0)){
	OldEl[idot-1]=MatArray[idot-1].GetMatEl(stcfj,stcfi);
	typep=typej;
	type=typei;
      }

      // Get SingleSite QNumbers
      int pos=0;
      int iblssp=pSingleSite->GetBlockFromSt(typep,pos);
      int iblss=pSingleSite->GetBlockFromSt(type,pos);

      double Qtildep=pSingleSite->GetQNumber(iblssp,0);
      double Stildep=pSingleSite->GetQNumber(iblssp,1);
      double Sztildep=pSingleSite->GetQNumber(iblssp,2);

      double Qtilde=pSingleSite->GetQNumber(iblss,0);
      double Stilde=pSingleSite->GetQNumber(iblss,1);
      double Sztilde=pSingleSite->GetQNumber(iblss,2);

      double Soldp=Si-Sztildep;
      double Sold=Si-Sztilde;


      // Check Fermi Sign (single channel)
      double FermiSign=1.0;
      if (dEqual(Qtilde,0.0)) FermiSign=-1.0;

      for (int iqn=0;iqn<NqnsSite;iqn++){
	siteqnumsp[iqn]=pSingleSite->GetQNumber(iblssp,iqn);
	siteqnums[iqn]=pSingleSite->GetQNumber(iblss,iqn);
      }

      //Loop in spins
      for (int sigma=-1;sigma<=1;sigma+=2){
	double dSigma=0.5*(double)sigma;
	double FullMatEl=0.0;
	double auxCG[]={0.0,0.0,0.0,0.0};

	// Loop in Szold
	for (double Szold=Sold;Szold>=-Sold;Szold-=1.0){
	  double Szoldp=Szold-dSigma;
	  Sztilde=Si-Szold;
	  Sztildep=Sztilde+dSigma;
	  // Changing dSigma changes Sztildep

	  // Site matrix element: finds the bleock
	  siteqnumsp[2]=Sztildep;
	  siteqnums[2]=Sztilde;

	  // OneChQS: one block per state only
	  int sitestatep=pSingleSite->GetBlockFromQNumbers(siteqnumsp);
	  int sitestate=pSingleSite->GetBlockFromQNumbers(siteqnums);

	  if ( (dLEqual(fabs(Sztildep),Stildep))&&
	       (dLEqual(fabs(Szold),Sold))&&
	       (dLEqual(fabs(Szoldp),Soldp)) ){
	    // CG coefs
	    auxCG[0]=CGordan(Sold,Szold,
			     Stilde,Sztilde,
			     Si,Si);				 

	    auxCG[1]=CGordan(Soldp,Szoldp,
			     Stildep,Sztildep,
			     Si,Si);

	    auxCG[2]=CGordan(Soldp,Szoldp,
			     0.5,dSigma,
			     Sold,Szold);
				  
	    auxCG[3]=OneCh_fd_table(sigma,sitestatep,sitestate); 

	    FullMatEl+=auxCG[0]*auxCG[1]*auxCG[2]*auxCG[3]*FermiSign;


	  }
	  // END Calc coefs safeguard
	}
	// End loop in Szold

	MatEl+=chi_N[idot-1]*OldEl[idot-1]*FullMatEl;

// 	if (( (ist==5)||(ist==6) )&&( (jst==9)||(jst==9) ))
// 	  cout << " idot = " << idot
// 	       << " FermiSign = " << FermiSign
// 	       << " sigma = " << sigma
// 	    //<< " Szold = " << Szold
// 	       << " OldEl = " << OldEl[idot-1]
// 	       << " FullMatEl = " << FullMatEl
// 	       << " MatEl = " << MatEl
// 	       << endl;

      }
      // end loop in sigma

    }
    // end sum in dots

  }
  // end if ist==jst


  delete[] siteqnumsp;
  delete[] siteqnums;


  return(MatEl);

}

/////////////////////////

/////////////////////////////////
/////////////////////////////////
///
/// <Q' Sz' | Hm1_DQD| Q Sz> :
///

double TwoDotQSz_Hm1_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst){

  // Oct 11
  // Calculates matrix elements for Hm1 (DQD) in the QSz basis.
  //
  // Parameters: 
  //  small_lambda  - Params[0]
  // Abasis       - DQD with small_lambda=0

  double MatEl=0.0;

  double small_lambda=Params[0];

  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Szi=pAbasis->GetQNumberFromSt(ist,1);

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Szj=pAbasis->GetQNumberFromSt(jst,1);

  if ( (dNEqual(Qi,Qj))||(dNEqual(Szi,Szj)) ) return(0.0);


  //////////**///////////////

  if (ist==jst){ // Diagonal terms
    MatEl=pAbasis->dEn[ist];
  }
  else{

//     if ( (dEqual(Szi,1.0))||(dEqual(Szi,-1.0)) ) MatEl=0.0;
    
//     if ( (dEqual(Qi,2.0))||(dEqual(Qi,-2.0)) ) MatEl=0.0;


    if (dEqual(Qi,-1.0)) MatEl=small_lambda;

    if (dEqual(Qi,1.0)) MatEl=-small_lambda;


    if ( (dEqual(Qi,0.0))&&(dEqual(Szi,0.0)) ) {
      //if ( (ist==5) MatEl=small_lambda; WRONG!!
      if ( (ist==5)&&
	   ((jst==6)||(jst==7)) ) MatEl=small_lambda;
      //if (ist==8) MatEl=-small_lambda; WRONG!!
      if ( (ist==8)&&
	   ((jst==6)||(jst==7)) ) MatEl=-small_lambda;


      if ( ((ist==6)||(ist==7))&&
	   (jst==5) ) MatEl=small_lambda;


      if ( ((ist==6)||(ist==7))&&
	   (jst==8) ) MatEl=-small_lambda;

    }
    // end if Qi=0, Szi=0


  }
  // end if ist==jst


  return(MatEl);

}

/////////////////////////


/////////////////////////////////
/////////////////////////////////
///
/// <Q' Sz' | H_DQD + site | Q Sz> :
///

double OneChQSz_H0DQD_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst){

  // October 2011
  // Calculates matrix elements for H0 (DQD+site) in the QSz basis.
  //
  // Parameters: 
  //  chi1        - Params[0]
  //  chi2        - Params[1]
  //  <Qold-1 Szold-1/2|| f1_up |Qold Szold> - MatArray[0] 
  //  <Qold-1 Szold+1/2|| f1_dn |Qold Szold> - MatArray[1] 
  //  <Qold-1 Szold-1/2|| f2_up |Qold Szold> - MatArray[2] 
  //  <Qold-1 Szold+1/2|| f2_dn |Qold Szold> - MatArray[3] 
  // Abasis       - DQD + site
  // pSingleSite  - OneChQSz single site

  double MatEl=0.0;

  double chi_N[2]={Params[0],Params[1]};

  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Szi=pAbasis->GetQNumberFromSt(ist,1);

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Szj=pAbasis->GetQNumberFromSt(jst,1);

  if ( (dNEqual(Qi,Qj))||(dNEqual(Szi,Szj)) ) return(0.0);

  //

  double OldEl[4]={0.0,0.0,0.0,0.0};
  

  if (ist==jst){ // Diagonal terms
    MatEl=pAbasis->dEn[ist];
  }
  else{
    int typei=pAbasis->iType[ist];
    int stcfi=pAbasis->StCameFrom[ist];
    
    int typej=pAbasis->iType[jst];
    int stcfj=pAbasis->StCameFrom[jst];

    // Get c1_up matrix element from MatArray[0]
    // Get c1_dn matrix element from MatArray[1]
    // Get c2_up matrix element from MatArray[2]
    // Get c2_dn matrix element from MatArray[3]
    // Get f0_sigma matrix element using pSingleSite


    // Sum in dots

    int icounter=0;
    for (int idot=1;idot<=2;idot++){

      // Add a sum in spins.
      for (int sigma=1;sigma>=-1;sigma-=2){
	double dSigma=0.5*(double)sigma;

	OldEl[icounter]=MatArray[icounter].GetMatEl(stcfi,stcfj); //

	int typep=typei;
	int type=typej;

	// if zero, try h.c.
	if (dEqual(fabs(OldEl[icounter]),0.0)){
	  OldEl[icounter]=MatArray[icounter].GetMatEl(stcfj,stcfi);
	  typep=typej;
	  type=typei;
	}

	// Get SingleSite QNumbers

 	double Qtilde=pSingleSite->GetQNumber(type,0);

	// Check Fermi Sign (single channel)
	double FermiSign=1.0;
	if (dEqual(Qtilde,0.0)) FermiSign=-1.0;

	double FullMatEl=OneCh_fd_table(sigma,typep,type)*FermiSign;

	MatEl+=chi_N[idot-1]*OldEl[icounter]*FullMatEl;

	// 	if (( (ist==5)||(ist==6) )&&( (jst==9)||(jst==9) ))
	// 	  cout << " idot = " << idot
	// 	       << " FermiSign = " << FermiSign
	// 	       << " sigma = " << sigma
	// 	    //<< " Szold = " << Szold
	// 	       << " OldEl = " << OldEl[icounter]
	// 	       << " FullMatEl = " << FullMatEl
	// 	       << " MatEl = " << MatEl
	// 	       << endl;


	icounter++;
      }
      // end loop in sigma
 
    }
    // end sum in dots

  }
  // end if ist==jst


  return(MatEl);

}
/////////////////////////


/////////////////////////////////
/////////////////////////////////
///
/// <Nup Pdn | H_Majorana | Nup  Pdn>_N=-1> :
///
/// Complex BUT Params is double!

complex<double> OneChNupPdn_Hm1_Majorana_MatEl(vector<double> Params,
					       CNRGbasisarray* pAbasis, 
					       CNRGbasisarray* pSingleSite,
					       CNRGmatrix* MatArray,
					       int ist, int jst){

  complex<double> cMatEl=ZeroC;

  // define t+ and t- (complex variables)

  complex<double> tplus,tminus;

  double t1=Params[0];
  double t2=Params[1];
  double phi_mag=Params[2];
  double en=Params[3];

//   tplus.real()=t1+t2*cos(phi_mag);
//   tminus.real()=t1-t2*cos(phi_mag);

//   tplus.imag()=t2*sin(phi_mag);
//   tminus.imag()=-t2*sin(phi_mag);



  double re_tplus  = (t1+t2*cos(phi_mag))/sqrt(2.0);
  double re_tminus = (t1-t2*cos(phi_mag))/sqrt(2.0);
  double im_tplus  = (t2*sin(phi_mag))/sqrt(2.0);
  double im_tminus = (-1.0*t2*sin(phi_mag))/sqrt(2.0);

  tplus  = complex<double>(re_tplus , im_tplus );
  tminus = complex<double>(re_tminus, im_tminus);

  // Matrix elements

  double Nupi=pAbasis->GetQNumberFromSt(ist,0);
  double Pdni=pAbasis->GetQNumberFromSt(ist,1);

  double Nupj=pAbasis->GetQNumberFromSt(jst,0);
  double Pdnj=pAbasis->GetQNumberFromSt(jst,1);

  if ( (dNEqual(Nupi,Nupj))||(dNEqual(Pdni,Pdnj)) ) return(ZeroC);

  //
  if (ist==jst){ // Diagonal terms
//     cMatEl.real()=pAbasis->dEn[ist];
//     cMatEl.imag()=0.0;
    cMatEl=complex<double>(pAbasis->dEn[ist],0.0);
  }
  else{
    if ( (ist==1)&&(jst==0) ) cMatEl=tminus;
    if ( (ist==0)&&(jst==1) ) cMatEl=std::conj(tminus);

    if ( (ist==3)&&(jst==2) ) cMatEl=tplus;
    if ( (ist==2)&&(jst==3) ) cMatEl=std::conj(tplus);

//     if ( (ist==5)&&(jst==4) ) cMatEl=(-1.0)*tminus;
//     if ( (ist==4)&&(jst==5) ) cMatEl=std::conj((-1.0)*tminus);
    if ( (ist==5)&&(jst==4) ) cMatEl=tminus;
    if ( (ist==4)&&(jst==5) ) cMatEl=std::conj(tminus);

//     if ( (ist==7)&&(jst==6) ) cMatEl=(-1.0)*tplus;
//     if ( (ist==6)&&(jst==7) ) cMatEl=std::conj((-1.0)*tplus);
    if ( (ist==7)&&(jst==6) ) cMatEl=tplus;
    if ( (ist==6)&&(jst==7) ) cMatEl=std::conj(tplus);

  }
  // end if ist==jst

  return(cMatEl);
}
// end OneChNupPdn_Hm1_Majorana_MatEl
/////////////////////////
