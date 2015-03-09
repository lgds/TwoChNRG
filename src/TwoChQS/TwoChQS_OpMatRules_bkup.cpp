

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

#include "TwoChQS.hpp"

//////
////// Two Channel QS, <Q' S'|| cd || Q S>
//////
////////////////////////////////////////////
///////////     MatElChecks    /////////////
////////////////////////////////////////////


bool TwoChQS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2)
{

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

////////////////////////
////////////////////////
////////////////////////




////////////////////////////////////////////
///////////   Mat ElCalculations   /////////
////////////////////////////////////////////


/////////////////////////////////////////
////                                 ////
////  Starting with a basic one:     ////
////  dot + first site               ////
////                                 ////
/////////////////////////////////////////

double TwoChQS_d_fd_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst)
{

  //
  // Calculates sum_{sigma} <ist| d_{sigma} f+_{0 ich sigma} + h.c |jst>
  // 
  // Two channels with Q,S symmetry
  //
  //  - ich=1 or 2 is the first input in Params 
  //  - ist!=jst (returns zero otherwise)
  //  - MatArray[0] carries the REDUCED <istcf || d || jstcf> elements

  int ich=(int)Params[0];

  if ( (ich!=1)&&(ich!=2) ) 
    {
      cout << " d_fd_MatEl: Invalid channel choice" << endl;
      return(0.0);
    }

  if (ist==jst) return(0.0); //Off-diagonal terms only


  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Si=pAbasis->GetQNumberFromSt(ist,1);   

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Sj=pAbasis->GetQNumberFromSt(jst,1);   

  if (dNEqual(Qi,Qj)||dNEqual(Si,Sj)) return(0.0);

  double HmatEl=0.0;
  
  int type_i=pAbasis->iType[ist];
  int type_j=pAbasis->iType[jst];

  // Reduced matrix element from <istcf| d | jstcf>
  double OldEl=MatArray[0].GetMatEl(pAbasis->StCameFrom[ist],
				    pAbasis->StCameFrom[jst]);
  int typep=type_i;
  int type=type_j;

  // if zero, try the h.c term
  if (dEqual(OldEl,0.0))
    {
      OldEl=MatArray[0].GetMatEl(pAbasis->StCameFrom[jst],
				 pAbasis->StCameFrom[ist]);
      typep=type_j;
      type=type_i;
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


  // Check Fermi Sign
  double FermiSign=1.0;
  // will be -1 if only one of them is zero (A XOR B)
  // 

  if ( (dEqual(Qtilde,1.0)||dEqual(Qtilde,-1.0))&&
       (dEqual(Sztilde,0.5)||dEqual(Sztilde,-0.5)) ) FermiSign=-1.0;
  // Fermi sign appears IF there is an 
  // odd number of fermionic ops in site state 
  // Only one electron "alone" in one channel, 
  // e.g. |0 up> or |dn 0> (not |up dn>

  // now get
  // auxEl=sum_ch <fd_{ch 0 sigma} f_{-1 sigma}>
  // REMEMBER:
  // |Q S>_N+1 = sum <So Szo Stil Sztil|SS>|QoSo>|Qtil Stil>

  // NO Loop in channels


  double auxCG[]={0.0,0.0,0.0,0.0};
  double siteqnumsp[]={Qtildep,Stildep,Sztildep};
  double siteqnums[]={Qtilde,Stilde,Sztilde};
  // Loop in sigma
  for (int sigma=-1;sigma<=1;sigma+=2)
    {
      double dSigma=0.5*(double)sigma;

      // Loop in Szold
      for (double Szold=Sold;Szold>=-Sold;Szold-=1.0)
	{
	  double Szoldp=Szold-dSigma;
	  Sztilde=Si-Szold;
	  Sztildep=Sztilde+dSigma;
	  // Changing dSigma changes Sztildep
	  // Need an extra loop in the "p" block

	  // Site matrix element: finds the block
	  siteqnums[2]=Sztilde;
	  siteqnumsp[2]=Sztildep;
	  int siteblock=pSingleSite->GetBlockFromQNumbers(siteqnums);
	  int siteblockp=pSingleSite->GetBlockFromQNumbers(siteqnumsp);

  
	  if ( (siteblock<pSingleSite->NumBlocks())&&
	       (siteblockp<pSingleSite->NumBlocks())&&
	       (dLEqual(fabs(Sztilde),Stilde))&&
	       (dLEqual(fabs(Sztildep),Stildep))&&
	       (dLEqual(fabs(Szold),Sold))&&
	       (dLEqual(fabs(Szoldp),Soldp)) )
	    {
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

// 	      cout << " Stilde = " << Stilde 
// 		   << " Sztilde = " << Sztilde
// 		   << " Stildep = " << Stildep 
// 		   << " Sztildep = " << Sztildep 
// 		   << endl;
				  
	      //Loop in site blocks!! Necessary!
	      for (int sitestatep=pSingleSite->GetBlockLimit(siteblockp,0);sitestatep<=pSingleSite->GetBlockLimit(siteblockp,1);sitestatep++)
		{
		  for (int sitestate=pSingleSite->GetBlockLimit(siteblock,0);sitestate<=pSingleSite->GetBlockLimit(siteblock,1);sitestate++)
		    {
		      double SpSm[2];
		      SpSm[0]=TwoChQS_SpSm_table(sitestate,type);
		      SpSm[1]=TwoChQS_SpSm_table(sitestatep,typep);
		      auxCG[3]=TwoChQS_fd_table(ich,sigma,
						sitestatep,sitestate);
 
		      HmatEl+=auxCG[0]*auxCG[1]*auxCG[2]*OldEl*auxCG[3]*FermiSign*SpSm[0]*SpSm[1];
		    }
		  // end loop in site block
		}
	      // end loop in site blockp
	    }
	  // END Calc coefs safeguard
	}
      // End loop in Szold
    }
  //END Sum in sigma


  return(HmatEl);

}



////////////////////////////////////////////
///////////       CM phonons       /////////
////////////////////////////////////////////

double TwoChQS_H0ph_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst)
{

  double HmatEl=0.0;

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

  // MatArray[1]: (nd_sigma-1) in the OLD basis
  double ndm1MatEl=MatArray[1].GetMatEl(stcf_i,stcf_j);
  // MatArray[2]: reduced <||d||> elements in the OLD basis
  CNRGmatrix Mat_d[1];
  Mat_d[0]=MatArray[2];

//   cout << " mi = " << mi 
//        << " mj = " << mj  
//        << endl;
//    cout << " Hm1 = " << Hm1MatEl
//         << " nd-1 = " << ndm1MatEl
//         << " phonon = " << phonon_matel 
//         << endl;

  vector<double> ParamsCh;
  double Hd_ch[2]={0.0,0.0};
  ParamsCh.push_back(1.0);
  Hd_ch[0]=delta_mimj*chi_N[0]*TwoChQS_d_fd_MatEl(ParamsCh,pAbasis,
				       pSingleSite,Mat_d,
				       ist,jst);


  ParamsCh[0]=2.0;
  Hd_ch[1]=alpha*chi_N[1]*phonon_matel*TwoChQS_d_fd_MatEl(ParamsCh,pAbasis,
							  pSingleSite,Mat_d,
							  ist,jst);
  //////////////////

  HmatEl=Hd_ch[0]+Hd_ch[1]; // always. It will be zero if not diagonal

//   cout << "HmatEl 1 = " << HmatEl << endl;
  
  if (type_i==type_j)   // Diagonal terms
    {
      HmatEl-=lambda*ndm1MatEl*phonon_matel; // Electron-phonon
//       cout << "HmatEl 2 = " << HmatEl << endl;
      HmatEl+=delta_mimj*Hm1MatEl; // Hubbard
//  	  cout << "HmatEl 3 = " << HmatEl << endl;
      if (stcf_i==stcf_j)
	HmatEl+=delta_mimj*mj*w0; // Phonon
      // 	  cout << "HmatEl 4 = " << HmatEl << endl;
    }


  
  return(HmatEl);
}
////////////////////////////////////////////
////////////////////////////////////////////



////////////////////////////////////////////
///////////   CM phonons w/ Transf /////////
////////////////////////////////////////////



double TwoChQS_H0phwTransf_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst)
{

  double HmatEl=0.0;

  double MatEl[3]={0.0,0.0,0.0}; // A bit confusing...
  // These are the 3 matrix elements:
  // <|cddag f1+h.c.|> 
  // <|cddag f2+h.c.|> 
  // <|(nd-1)cddag f2+h.c.|> 
  double OldEl[2]={0.0,0.0};
  // Two "old" matrix elements
  // <old'|cd |old>
  // <old'|(nd-1)cd|old>

  double phonon_matel[3]={0.0,0.0,0.0};
  // <mi|exp()|mj>
  // <mi|exp()|mj-1>
  // <mi|exp()|mj+1>

  double chi_N[2]={Params[0],Params[1]};
  double w0=Params[2];
  double lambda=Params[3];
  double alpha=Params[4];

  // should be block diagonal. Getting qnumbers
  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Si=pAbasis->GetQNumberFromSt(ist,1);   

  int type_i=pAbasis->iType[ist];
  double mphi=pAbasis->GetStCameFromQNumberFromSt(ist,2);

  int type_j=pAbasis->iType[jst];
  double mphj=pAbasis->GetStCameFromQNumberFromSt(jst,2);


//   cout << " ist = " << ist 
//        << " mi = " << mphi
//        << " jst = " << jst
//        << " mj = " << mphj
//        << endl;


  if (ist==jst)   // Diagonal terms
    {
      HmatEl=pAbasis->dEn[ist];
    }
  else 		  //Off-diagonal terms
    {
      for (int ich=0;ich<=1;ich++) // Loop in channels
	{
	  OldEl[ich]=MatArray[ich].GetMatEl(pAbasis->StCameFrom[ist],pAbasis->StCameFrom[jst]);
//  	  cout << " ich = " << ich
// 	       << " OldEl = " << OldEl[ich]
// 	       << endl;

	  double CheckHc=MatArray[0].GetMatEl(pAbasis->StCameFrom[ist],pAbasis->StCameFrom[jst]);


	  int typep=type_i;
	  int type=type_j;
	  double mphp=mphi;
	  double mph=mphj;

	  // if <Qoldp Soldp||cd||Qold Sold> zero, use the h.c term
	  // if (dEqual(fabs(OldEl[0]),0.0))
	  if (dEqual(fabs(CheckHc),0.0))
	    {
	      OldEl[ich]=MatArray[ich].GetMatEl(pAbasis->StCameFrom[jst],
						pAbasis->StCameFrom[ist]);

// 	      cout << " ich = " << ich
// 		   << " New OldEl = " << OldEl[ich]
// 		   << endl;

	      typep=type_j;
	      type=type_i;
	      mphp=mphj;
	      mph=mphi;
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


	  // Check Fermi Sign
	  double FermiSign=1.0;
	  // will be -1 if only one of them is zero (A XOR B)
	  // Check this one. 
	  // Fermi sign now depends on which type the thing is!!!
	  // In other words, it will be -1 if the state has only one 
	  // electron in only one of the channels
	  // Single site states +-1,+-1/2!!
	  // 
	  if ( (dEqual(Sztilde,0.5)||dEqual(Sztilde,-0.5)) ) FermiSign=-1.0;

	  double siteqnumsp[]={Qtildep,Stildep,Sztildep};
	  double siteqnums[]={Qtilde,Stilde,Sztilde};

	  double FullMatEl=0.0;
	  //Loop in spins
	  for (int sigma=-1;sigma<=1;sigma+=2)
	    {
	      double dSigma=0.5*(double)sigma;
	      double auxCG[]={0.0,0.0,0.0,0.0};

	      // Loop in Szold
	      for (double Szold=Sold;Szold>=-Sold;Szold-=1.0)
		{
		  double Szoldp=Szold-dSigma;
		  Sztilde=Si-Szold;
		  Sztildep=Sztilde+dSigma;
		  // Changing dSigma changes Sztildep
		  // Need an extra loop in the "p" block

		  // Site matrix element: finds the block
		  siteqnums[2]=Sztilde;
		  siteqnumsp[2]=Sztildep;
		  int siteblock=pSingleSite->GetBlockFromQNumbers(siteqnums);
		  int siteblockp=pSingleSite->GetBlockFromQNumbers(siteqnumsp);
			      
		  if ( (siteblock<pSingleSite->NumBlocks())&&
		       (siteblockp<pSingleSite->NumBlocks())&&
		       (dLEqual(fabs(Sztildep),Stildep))&&
		       (dLEqual(fabs(Szold),Sold))&&
		       (dLEqual(fabs(Szoldp),Soldp)) )
		    {
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
				  
		      //Loop in site blocks!! Necessary!
		      for (int sitestatep=pSingleSite->GetBlockLimit(siteblockp,0);sitestatep<=pSingleSite->GetBlockLimit(siteblockp,1);sitestatep++)
			{
			  for (int sitestate=pSingleSite->GetBlockLimit(siteblock,0);sitestate<=pSingleSite->GetBlockLimit(siteblock,1);sitestate++)
			    {
			      double SpSm[2];
			      SpSm[0]=TwoChQS_SpSm_table(sitestate,type);
			      SpSm[1]=TwoChQS_SpSm_table(sitestatep,typep);
			      auxCG[3]=TwoChQS_fd_table(ich+1,sigma,
							sitestatep,sitestate);
				  
			      FullMatEl+=auxCG[0]*auxCG[1]*auxCG[2]*auxCG[3]*FermiSign*SpSm[0]*SpSm[1];

			    }
			  // end loop in site block
			}
		      // end loop in site blockp
		    }
		  // END Calc coefs safeguard
		}
	      // End loop in Szold
	    }
	  // end loop in sigma

	  // Notice
	  // MatEl[0]=FullMat(ch=0)*OldMat[0]
	  // MatEl[1]=FullMat(ch=1)*OldMat[0]
	  // MatEl[2]=FullMat(ch=1)*OldMat[1]

	  if (ich==0)
	    {
	      MatEl[ich]=OldEl[ich]*FullMatEl; // <|ddag f1+h.c.|>
	      // Phonon part
	      phonon_matel[ich]=Calc_phMatEl(mphp,mph,lambda,w0);
	    }
	  else
	    {
	      MatEl[ich]=OldEl[ich-1]*FullMatEl; // <|ddag f2+h.c.|> 
	      MatEl[ich+1]=OldEl[ich]*FullMatEl; //<|(nd-1)ddag f2+h.c.|>
	      // Phonon part
	      phonon_matel[ich]=sqrt(mph)*Calc_phMatEl(mphp,mph-1.0,lambda,w0);
	      phonon_matel[ich+1]=sqrt(mph+1.0)*Calc_phMatEl(mphp,mph+1.0,lambda,w0);
	    }

// 	  cout << " ich = " << ich
// 	       << " FullMatEl = " << FullMatEl
// 	       << endl;
        }
      //end loop in channels

//       cout << " Gamma1 = " << chi_N[0]
// 	   << " Gamma2 = " << chi_N[1]
// 	   << " lambda = " << lambda
// 	   << " w0 = " << w0
// 	   << endl
// 	   << "  MatEl[0] = " << MatEl[0]
// 	   << "  MatEl[1] = " << MatEl[1]
// 	   << "  MatEl[2] = " << MatEl[2]
// 	   << endl
// 	   << "  OldEl[0] = " << OldEl[0]
// 	   << "  OldEl[1] = " << OldEl[1]
// 	   << endl
// 	   << " phonon_matel[0] = " << phonon_matel[0]
// 	   << " phonon_matel[1] = " << phonon_matel[1]
// 	   << " phonon_matel[2] = " << phonon_matel[2]
// 	   << endl;

      if (dNEqual(w0,0.0))
	HmatEl=chi_N[0]*phonon_matel[0]*MatEl[0]
	  +chi_N[1]*(  ( phonon_matel[1]+phonon_matel[2] )*MatEl[1]
		       +2.0*(lambda/w0)*phonon_matel[0]*MatEl[2]  );
      else
	HmatEl=0.0;


    }
  // END if ist=jst

  

  return(HmatEl);

}
////////////////////////////////////////////
////////////////////////////////////////////



////////////////////////////////////////////
///////////     Kondo model    /////////////
////////////////////////////////////////////

double TwoChQS_H0Kondo_MatEl(vector<double> Params,
			  CNRGbasisarray* pAbasis, 
			  CNRGbasisarray* pSingleSite,
			  CNRGmatrix* MatArray,
			  int ist, int jst)
{

  // Diagonal in QS basis!!

  if (ist!=jst) return(0.0);

  double HmatEl=0.0;
  double auxterm[]={0.0,0.0,0.0};

  double J[2];
  J[0]=Params[0];
  J[1]=Params[1];

  // should be block diagonal. Getting qnumbers
  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Si=pAbasis->GetQNumberFromSt(ist,1);   

  int type_i=pAbasis->iType[ist];
  int type_j=pAbasis->iType[jst];

  double auxCG[]={0.0,0.0,0.0,0.0};

  //       // Get SingleSite QNumbers
  int pos=0;
  int iblss=pSingleSite->GetBlockFromSt(type_i,pos);

  double Qtilde=pSingleSite->GetQNumber(iblss,0);
  double Stilde=pSingleSite->GetQNumber(iblss,1);
  double Sztilde=pSingleSite->GetQNumber(iblss,2);

  double Sold=Si-Sztilde;

  double siteqnums[]={Qtilde,Stilde,Sztilde};
  double siteqnumsp1[]={Qtilde,Stilde,Sztilde};
  double siteqnumsm1[]={Qtilde,Stilde,Sztilde};


  // Loop in Szoldp, Szold

  for (double Szold=-Sold;Szold<=Sold;Szold+=1.0)
    {
      Sztilde=Si-Szold;
      double Sztildep1=Sztilde+1.0;
      double Sztildem1=Sztilde-1.0;

      siteqnums[2]=Sztilde;
      if (dLEqual(fabs(Sztildep1),Stilde))
	siteqnumsp1[2]=Sztildep1;
      if (dLEqual(fabs(Sztildem1),Stilde))
	siteqnumsm1[2]=Sztildem1;

      int siteblock=pSingleSite->GetBlockFromQNumbers(siteqnums);
      int siteblockp1=pSingleSite->GetBlockFromQNumbers(siteqnumsp1);
      int siteblockm1=pSingleSite->GetBlockFromQNumbers(siteqnumsm1);

      // CG coefs
      auxCG[0]=CGordan(Sold,Szold,
		       Stilde,Sztilde,
		       Si,Si);


      auxCG[1]=CGordan(Sold,Szold+1.0,
		       Stilde,Sztildem1,
		       Si,Si);
      auxCG[2]=CGordan(Sold,Szold-1.0,
		       Stilde,Sztildep1,
		       Si,Si);

      //Loop in site blocks!! Necessary!

      for (int sitestate=pSingleSite->GetBlockLimit(siteblock,0);sitestate<=pSingleSite->GetBlockLimit(siteblock,1);sitestate++)
	{
	  
	  double SpSm[2];
	  SpSm[0]=TwoChQS_SpSm_table(sitestate,type_j);


	  // Term S+ fd_dn f_up
	  if (siteblock!=siteblockm1)
	    for (int sitestatem1=pSingleSite->GetBlockLimit(siteblockm1,0);
		 sitestatem1<=pSingleSite->GetBlockLimit(siteblockm1,1);sitestatem1++)
	      {
		SpSm[1]=TwoChQS_SpSm_table(sitestatem1,type_i);

		// Loop in channels
		auxterm[0]=0.0;
		for (int ich=1;ich<=2;ich++)
		  {
		    // Term S+ fd_dn f_up
		    double ChFlipTerm=TwoChQS_fdupfdn_table(ich,sitestate,sitestatem1)
		      *SpSm[0]*SpSm[1];
		    double SpinTerm=Splus(Sold,Szold+1.0,Sold,Szold);
		    HmatEl+=2.0*J[ich-1]*auxCG[0]*auxCG[1]*SpinTerm*ChFlipTerm;
		  }
		// END Loop in channels

	      }
	  // end loop in site blockm1


	  // Term S- fd_up f_dn
	  if (siteblock!=siteblockp1)
	    for (int sitestatep1=pSingleSite->GetBlockLimit(siteblockp1,0);
		 sitestatep1<=pSingleSite->GetBlockLimit(siteblockp1,1);sitestatep1++)
	      {

		SpSm[1]=TwoChQS_SpSm_table(sitestatep1,type_i);

		// Loop in channels
		auxterm[0]=0.0;
		for (int ich=1;ich<=2;ich++)
		  {
		    // Term S- fd_up f_dn
		    double ChFlipTerm=TwoChQS_fdupfdn_table(ich,sitestatep1,sitestate)
		      *SpSm[0]*SpSm[1];
		    double SpinTerm=Sminus(Sold,Szold-1.0,Sold,Szold);
		    HmatEl+=2.0*J[ich-1]*auxCG[0]*auxCG[2]*SpinTerm*ChFlipTerm;
		  }
		// END Loop in channels

	      }
	  // end loop in site blockp1

	  // Term Sz ni_up - ni_dn

	  // Loop in channels
	  auxterm[0]=0.0;
	  for (int ich=1;ich<=2;ich++)
	    {
	      HmatEl+=J[ich-1]*Szold*(auxCG[0]*auxCG[0])*TwoChQS_Szch_table(ich,sitestate)*SpSm[0];

	    }
	  // END Loop in channels

	}
      // end loop in site block
    }
  // end loop in Szold

  return(HmatEl);

}
//
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////


double TwoChQS_cdPhonon_MatEl(CNRGbasisarray* pAbasis, 
			       CNRGbasisarray* pSingleSite,
			       int ist, int jst, int ich)
{


  double Qi=pAbasis->GetQNumberFromSt(ist,0);
  double Si=pAbasis->GetQNumberFromSt(ist,1);   

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Sj=pAbasis->GetQNumberFromSt(jst,1);   


  int typei=pAbasis->iType[ist];
  int stcfi=pAbasis->StCameFrom[ist];

  int typej=pAbasis->iType[jst];
  int stcfj=pAbasis->StCameFrom[jst];

  int mi=pAbasis->iDegen[ist];
  int mj=pAbasis->iDegen[jst];

  double auxEl=0.0;

  // delta mi mj
  if ( (stcfi==stcfj)&&(mi==mj) )
    {
      // New approach

      // Get SingleSite QNumbers
      int pos=0;
      int iblssi=pSingleSite->GetBlockFromSt(typei,pos);
      int iblssj=pSingleSite->GetBlockFromSt(typej,pos);


      double Qtildei=pSingleSite->GetQNumber(iblssi,0);
      double Stildei=pSingleSite->GetQNumber(iblssi,1);
      double Sztildei=pSingleSite->GetQNumber(iblssi,2);
			 
      double Qtildej=pSingleSite->GetQNumber(iblssj,0);
      double Stildej=pSingleSite->GetQNumber(iblssj,1);
      double Sztildej=pSingleSite->GetQNumber(iblssj,2);

      double siteqnumsi[]={Qtildei,Stildei,Sztildei};
      double siteqnumsj[]={Qtildej,Stildej,Sztildej};

      double Scfi=Si-Sztildei;
      double Scfj=Sj-Sztildej; // should be the same!!

      double auxBasis[2]={0.0,0.0};
      double auxCG[2]={0.0,0.0};
      double CGnorm=1.0;

      int sigma=0;

      // Choose which sigma
      auxCG[0]=CGordan(Si,Si,0.5,-0.5,Sj,Sj);
      auxCG[1]=CGordan(Si,Si,0.5,0.5,Sj,Sj);
      if ( dNEqual(auxCG[0],0.0) )
	{sigma=-1;CGnorm=auxCG[0];}
      else
	{
	  if ( dNEqual(auxCG[1],0.0) )
	    {sigma=1;CGnorm=auxCG[1];}
	  else auxEl=0.0;
	}

//       cout << " UpMat: "
// 	   << " ist = " << ist 
// 	   << " jst = " << jst 
// 	   << " CG0 = " << auxCG[0]
// 	   << " CG1 = " << auxCG[1]
// 	   << " sigma = " << sigma
// 	   << endl;

      double dSigma=0.5*sigma;
      // Now Sj=Si+dSigma necessarily
      // Calculate
      if (sigma!=0)
	{
	  for (double Szcfj=Scfj;
	       Szcfj>=-Scfj;Szcfj-=1.0)
	    {
	      Sztildei=Si-Szcfj; //Szcfi=Szcfj
	      Sztildej=Sj-Szcfj; 
	      siteqnumsi[2]=Sztildei;
	      int blockstatei=pSingleSite->GetBlockFromQNumbers(siteqnumsi);
	      siteqnumsj[2]=Sztildej;
	      int blockstatej=pSingleSite->GetBlockFromQNumbers(siteqnumsj);

	      auxCG[0]=CGordan(Scfi,Szcfj,
			       Stildei,Sztildei,
			       Si,Si);
	      auxCG[1]=CGordan(Scfj,Szcfj,
			       Stildej,Sztildej,
			       Sj,Sj);

	      for (
		   int sitestatei=pSingleSite->GetBlockLimit(blockstatei,0);
		   sitestatei<=pSingleSite->GetBlockLimit(blockstatei,1);
		   sitestatei++)
		{
		  for (
		       int sitestatej=pSingleSite->GetBlockLimit(blockstatej,0);
		       sitestatej<=pSingleSite->GetBlockLimit(blockstatej,1);
		       sitestatej++)
		    {
		      auxBasis[0]+=TwoChQS_fd_table(ich,sigma,sitestatej,sitestatei)*auxCG[0]*auxCG[1]*TwoChQS_SpSm_table(sitestatei,typei)*TwoChQS_SpSm_table(sitestatej,typej);

// 		      cout << " Szcfj = " << Szcfj
// 			   << " Stildei = " << Stildei 
// 			   << " Sztildei = " << Sztildei 
// 			   << " Stildej = " << Stildej 
// 			   << " Sztildej = " << Sztildej 
// 			   << " typei = " << typei
// 			   << " sitestatei = " << sitestatei
// 			   << " typej = " << typej
// 			   << " sitestatej = " << sitestatej
// 			   << " CG0 = " << auxCG[0]
// 			   << " CG1 = " << auxCG[1]
// 			   << endl;


		    }
		  // loop in blockstatej
		}
	      // loop in blockstatei  
	    }
	  // end Loop in Szcf

	  //Calculate fbasis
	  auxEl=auxBasis[0]/CGnorm;

	  //  				  cout << " auxBasis = " << auxBasis[0]
	  // 				       << " CGnorm = " << CGnorm
	  //  				       << endl;


	}
      // end if (sigma!=0)


    }
  // end if scfi=scfj


  return(auxEl);
}

////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////

double TwoChQS_cd_ich1_Phonon_MatEl(CNRGbasisarray* pAbasis, 
			      CNRGbasisarray* pSingleSite,
			      int ist, int jst){


  return( TwoChQS_cdPhonon_MatEl(pAbasis, pSingleSite,ist,jst,1) );

}

double TwoChQS_cd_ich2_Phonon_MatEl(CNRGbasisarray* pAbasis, 
			      CNRGbasisarray* pSingleSite,
			      int ist, int jst){


  return( TwoChQS_cdPhonon_MatEl(pAbasis, pSingleSite,ist,jst,2) );

}

////////////////////////////////////////////
////////////////////////////////////////////
