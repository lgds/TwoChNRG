
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

#include "OneChQS.hpp"


//////
////// One Channel QS, <Q' S'|| cd || Q S>
//////
////////////////////////////////////////////
///////////     MatElChecks    /////////////
////////////////////////////////////////////

bool OneChQS_cd_check(CNRGbasisarray *pAeigCut, int iblock1, int iblock2)
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

bool OneChQS_nd_check(CNRGbasisarray *pAeigCut, 
		      int iblock1, 
		      int iblock2)
{

//   double Qi=pAeigCut->GetQNumber(iblock1,0);
//   double Si=pAeigCut->GetQNumber(iblock1,1);

//   double Qj=pAeigCut->GetQNumber(iblock2,0);
//   double Sj=pAeigCut->GetQNumber(iblock2,1);

  
//   if (  (dEqual(Qj,Qi))&&(dEqual(Sj,Si)) )
//     return(true);
//   else
//     return(false);

  if ( iblock1==iblock2 )
    return(true);
  else
    return(false);



}

////////////////////////





////////////////////////////////////////////
///////////     Calculations   /////////////
////////////////////////////////////////////

double OneChQS_fN_MatEl(CNRGbasisarray *pAbasis, 
		        CNRGbasisarray *pSingleSite, 
			int ist, int jst)
{

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
	  auxBasis[i1]+=fd_table((int)(sigma/0.5),sitestatej,sitestatei)*auxCG[0]*auxCG[1];

	  //cout << " CG 1 = " << auxCG[0] << " CG 2 = " << auxCG[1] << " mat_el = " <<  fd_table((int)(sigma/0.5),sitestatej,sitestatei) << endl;
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
			int ist, int jst)
{

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
  if (dEqual(auxCG[3],0.0))
    {
      sigma=Si-Sj;
      auxCG[3]=CGordan(Si,Si,0.5,sigma,Sj,Sj);
      if (dEqual(auxCG[3],0.0)) return(0.0);
    }


  for (double Szcfj=-Scfj;
       Szcfj<=Scfj;Szcfj+=1.0)
    {
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

      //cout << " CG 1 = " << auxCG[0] << " CG 2 = " << auxCG[1] << " mat_el = " <<  fd_table((int)(sigma/0.5),sitestatej,sitestatei) << endl;
    }
  // end loop in Szcfj

  if (dNEqual(auxCG[3],0.0)) fbasis=fbasis/auxCG[3];
  else fbasis=0.0;

  return(fbasis);

}

////////////////////////
////////////////////////
////////////////////////


double OneChQS_nd_MatEl(CNRGbasisarray *pAbasis, 
			CNRGbasisarray *pSingleSite, 
			int ist, int jst)
{

  double fbasis=0.0;

  int typei=pAbasis->iType[ist];
  int typej=pAbasis->iType[jst];

  // Not diagnonal necessarily
  if (typei==typej) fbasis=1.0;


  return(fbasis);


}

////////////////////////
////////////////////////
////////////////////////

double OneChQS_HN_MatEl_KeepSz(vector<double> Params,
                          CNRGbasisarray* pAbasis,
                          CNRGbasisarray* pSingleSite,
                          CNRGmatrix* MatArray,
                          int ist, int jst)
{
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
  double Si=pAbasis->GetQNumberFromSt(ist,1);

  double Qj=pAbasis->GetQNumberFromSt(jst,0);
  double Sj=pAbasis->GetQNumberFromSt(jst,1);

  if (dNEqual(Qi,Qj)||dNEqual(Si,Sj)) return(0.0);

  if (ist==jst)
    MatEl=sqrt(Lambda)*pAbasis->dEn[ist];
  else
    {
      int type_i=pAbasis->iType[ist];
      int type_j=pAbasis->iType[jst];

      // Sum over spins
      int imat=0;
      for (int isigma=-1;isigma<=1;isigma+=2)
	{
	  // Matrix element from <istcf| fN_sigma | jstcf>
	  double OldEl=MatArray[imat].GetMatEl(pAbasis->StCameFrom[ist],
					       pAbasis->StCameFrom[jst]);
	  int typep=type_i;
	  int type=type_j;

	  // if zero, try the h.c term
	  if (dEqual(OldEl,0.0))
	    {
	      OldEl=MatArray[imat].GetMatEl(pAbasis->StCameFrom[jst],
					 pAbasis->StCameFrom[ist]);
	      typep=type_j;
	      type=type_i;
	    }

	  // Get SingleSite QNumbers

	  double Qtilde=pSingleSite->GetQNumber(type,0);
	  double Stilde=pSingleSite->GetQNumber(type,1);
	  double Sztilde=pSingleSite->GetQNumber(type,2);


	  // Check Fermi Sign
	  double FermiSign=1.0;
	  if (dEqual(Qtilde,0.0)) FermiSign=-1.0;
	  // Fermi sign appears IF there is an
	  // odd number of fermionic ops in site state

	  MatEl+=FermiSign*OldEl*fd_table(isigma,typep,type);
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
////////////////////////
