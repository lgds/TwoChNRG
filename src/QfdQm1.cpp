#include <iostream>
#include <vector>
//#include <algorithm>
#include <cmath>
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

using namespace std;

// Operator table one channel

//double fd_table(int isigma, int type_i, int type_j);

//double fd_table(int channel, int sigma, int type_i, int type_j);

double QfdQm1_totS(vector<double> Params, vector<int> Indexes, 
		   CNRGbasisarray* pSingleSite){

  //
  // Calculates the reduced matrix element
  //  <ST'|| f+n || ST >
  // given:
  //
  // Structure of SingleSite
  //
  // Params: in this order
  //     0 - STp
  //     1 - Scfp
  //     2 - ST
  //     3 - Scf
  //
  // Important!Assumes 
  //       Stildep is at position 3+NQNumbers-1 and
  //
  //       Stilde is at position 3+2*NQNumbers-1
  //
  // Indexes:
  //    1 - sitestatep
  //    2 - sitestate
  //    3 - call one ch (1) or two ch (2) table
  //    4 - if two channel, channel index (1 or 2)
  //
  // SingleSite: the last two qnumbers should be S and Sz.
  //

  int Nqns=pSingleSite->NQNumbers;


  double Sip=Params[0];
  double Scfp=Params[1];
  double Si=Params[2];
  double Scf=Params[3];

  int sitestp=Indexes[0];
  int sitest=Indexes[1];
  int Nochannels=Indexes[2];
  int ichannel=1;
  if (Nochannels>1) ichannel=Indexes[3];
    
  if (!dEqual(Scfp,Scf)) return(0.0);



  // needs to define arrays siteqnums and siteqnumsp

  double* SitepQnums = new double [Nqns];
  double* SiteQnums = new double [Nqns];

  for (int ii=0;ii<Nqns;ii++)
    {
      SitepQnums[ii]=pSingleSite->GetQNumber(sitestp,ii);
      SiteQnums[ii]=pSingleSite->GetQNumber(sitest,ii);
    }

  double Stildep=SitepQnums[Nqns-2];
  double Stilde=SiteQnums[Nqns-2];
  double MatEl=0.0;
  double auxEl=0.0;
  double auxCG[3]={0.0,0.0,0.0};

  // Calculating REDUCED element: need an extra loop in isigma 
  //    isigma =+1/2 (1) or -1/2 (-1)


  for (double Sztilde=-Stilde;Sztilde<=Stilde;Sztilde+=1.0)
    {

      SiteQnums[Nqns-1]=Sztilde;
      double Sztildep=Sip-Si+Sztilde; // fixed
      SitepQnums[Nqns-1]=Sztildep;
      int sitestate=pSingleSite->GetBlockFromQNumbers(SiteQnums);
      int sitestatep=pSingleSite->GetBlockFromQNumbers(SitepQnums);
      double Szcf=Si-Sztilde; // fixed. Szcfp=Szcf

//       cout << "Stildep = "  << Stildep << " Sztildep = " << Sztildep << endl; 
//       cout << "SitepQnums = " << SitepQnums[0] << " " 
// 	   << SitepQnums[1] << " " 
// 	   << SitepQnums[2];
//       cout << " istp = " << sitestatep << endl;
//       cout << "Scfp = "  << Scfp << " Szcfp = " << Szcf << endl;

//       cout << "Stilde = "  << Stilde << " Sztilde = " << Sztilde << endl;
//       cout << "SiteQnums = " << SiteQnums[0] << " " 
// 	   << SiteQnums[1] << " " 
// 	   << SiteQnums[2];
//       cout << " ist = " << sitestate << endl;
//       cout << "Scf = "  << Scf << " Szcf = " << Szcf << endl;
//       cout << " Sip = " << Sip << " Si = " << Si << endl;


      auxCG[0]=CGordan(Scf,Szcf,
		       Stilde,Sztilde,Si,Si);
      auxCG[1]=CGordan(Scfp,Szcf,
		       Stildep,Sztildep,Sip,Sip); // Szcfp=Szcf


      // Check which sigma matrix element is non-zero and stick with that
      int isigma=-1;
      if (Nochannels==1) MatEl=OneCh_fd_table(isigma,sitestatep,sitestate);
      if (Nochannels==2) MatEl=TwoCh_fd_table(ichannel,isigma,sitestatep,sitestate);
      if (dEqual(MatEl,0.0))
	{
	isigma=1;
	if (Nochannels==1) MatEl=OneCh_fd_table(isigma,sitestatep,sitestate);
	if (Nochannels==2) MatEl=TwoCh_fd_table(ichannel,isigma,sitestatep,sitestate);
      }

      double dsigma=0.5*(double)(isigma);
      auxCG[2]=CGordan(Si,Si,0.5,dsigma,Sip,Sip);
//       cout << " CG(" << Scf <<","<< Szcf <<","
// 	   << 0.5 <<","<< dsigma <<","
// 	   << Sip <<","<< Sip <<") = " << auxCG[2] << endl; 
      if (dEqual(auxCG[2],0.0)) {MatEl=0.0;auxCG[2]=1.234;}
//       cout << " matEl( " << isigma <<" , " << sitestatep << " , " << sitestate << ") = " << MatEl << endl;

//       cout << " CG 1 = " << auxCG[0] 
// 	   << " CG 2 = " << auxCG[1] 
// 	   << " CG 3 = " << auxCG[2] << endl; 

      auxEl+=auxCG[0]*auxCG[1]*MatEl/auxCG[2];

 
    }
  // end Sum in Stildez




  delete[] SitepQnums;
  delete[] SiteQnums;


  return(auxEl);

}
