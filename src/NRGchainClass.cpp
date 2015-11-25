
#include <vector>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#define pi 3.141592653589793238462643383279502884197169
#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
using namespace std;

//////////////

//////////////

void CNRGchain::PrintChain(int Nsites){

  cout << " Lambda = " << Lambda
       << " z_twist = " << z_twist
       << " BandOffset = " << BandOffset 
       << endl;
  cout << " ======================== " << endl;
  for (int ii=0;ii<Nsites;ii++){
    if (ii<en.size())
      cout << " en("<< ii<< ")= " << en[ii];
    if (ii<chin.size()){
      cout << "  -- chin("<< ii<< ")= " << chin[ii];
      cout << " - tau("<< ii<< ")= " << chin[ii]/ScaleFactor(ii);
    }
    cout << endl;
  }


}


//////////////

void CNRGchain::PrintAll(){

  int MaxSize=en.size()>chin.size()? en.size():chin.size();

  PrintChain(MaxSize);

}

//////////////

void CNRGchain::ReadParams(char arqname[],  
			   int Nparams){
  // instream
  ifstream InFile;
  double aux;

  vector<double>::iterator dit;

  InFile.open(arqname,ifstream::in);
  if (InFile.is_open()){
    for (int ii=0;ii<Nparams;ii++){
      InFile >> aux;
      HybFuncParams->push_back(aux);
    }
  }
  else{
    cout << "CNRGchain::ReadParams - error: can't open " << arqname << endl;
    exit(0);
  }
  InFile.close();

  int ic=0;
  cout << " CNRGchain: Params read from " << arqname << " ;"<< endl;
  for (dit=HybFuncParams->begin();dit<HybFuncParams->end();dit++){
    cout << " param" << ic << " = "<< *dit << endl;
    ic++;
  }
  cout << "..." << arqname << " closed."<< endl;
  // end print params
}


//////////////

void CNRGchain::ReadParams(char arqname[]){

  // 
  // Reads ALL parameters from an arbitrary file
  //
  // Example: a file x,Hyb(x)
  //

  // instream
  ifstream InFile;
  double daux;
  int iline=0;


  HybFuncParams->clear();

  InFile.open(arqname,ifstream::in);
  if (InFile.is_open()){
    while (!InFile.eof()){
      InFile >> daux;
      HybFuncParams->push_back(daux);
      iline++;
    }
    InFile.close();
    HybFuncParams->pop_back(); // Last element is garbage
  }
  else{
    cout << "CNRGchain::ReadParams - error: can't open " << arqname << endl;
    exit(0);
  }
  InFile.close();

  cout << "CNRGchain::ReadParams - No elements : " << iline-1 << endl;

  cout << " Read data (in 2 cols): " << endl;
  for (int il=0;il<HybFuncParams->size()-1;il+=2){
    cout << (*HybFuncParams)[il] << " -- " << (*HybFuncParams)[il+1] << endl; 
  }
  // end print data


}
// end ReadVecParams

//////////////

void CNRGchain::SetChainWilson(int Nsitesmax){

  en.clear();
  chin.clear();

  for (int Nsites=0;Nsites<=Nsitesmax;Nsites++){

    double daux[4];
    daux[0]=1.0-pow( Lambda,-((double)(Nsites)+1.0) );
    daux[1]=sqrt( 1.0-pow(Lambda,-(2.0*(double)(Nsites)+1.0)) );
    daux[2]=sqrt( 1.0-pow(Lambda,-(2.0*(double)(Nsites)+3.0)) );  
    
    daux[3]=daux[0]/(daux[1]*daux[2]);

    en.push_back(0.0);
    chin.push_back(daux[3]);

  }
  // end for


  // set Fsq also
  Fsq=2.0; 
}


//////////////

double CNRGchain::GetChin(int Nsites){

  vector<double>::iterator dit;
  dit=chin.begin()+Nsites;

  if (dit<chin.end())
    return(*dit);
  else
    return(-1.0);
}


//////////////

double CNRGchain::GetEn(int Nsites){

  vector<double>::iterator dit;
  dit=en.begin()+Nsites;

  if (dit<en.end())
    return(*dit);
  else
    return(-1.0);
}

//////////////

double CNRGchain::ScaleFactor(int n){

  //
  // Why n-1 instead of n?
  //
  // good question. Actually, it's n-1 for en and n for chi! (see notes).
  //
  // z-trick: 

  //double aux=pow(Lambda,0.5*(double)(n));
  double aux=pow(Lambda,0.5*(double)(n)+z_twist-1.0);
  //double aux=pow(Lambda,0.5*(double)(n-1));
  aux/=0.5*(1.0+(1.0/Lambda));
  

  return(aux);

}


//////////////

double CNRGchain::RHS(int n, double eam, double um, double um1, 
		      double en, double chinm1, bool scale){

//   if (scale){
//     eam*=ScaleFactor(n);
//     chinm1*=sqrt(Lambda);
//   }
// Note: en scales with SF(n-1) while chin scales with SF(n)
// sqrt(Lambda) factor multiplies the whole RHS if scaled
//

  if (scale){
    eam*=ScaleFactor(n-1);
    um*=sqrt(Lambda);
    um1*=sqrt(Lambda);
  }

  return((eam-en)*um-chinm1*um1);

}


//////////////


void CNRGchain::SetChainLanczos(int Nsitesmax){
  // Lanczos procedure
  en.clear();
  chin.clear();

  vector <double> ea;
  vector <double> F2a;

  vector <double> eb;
  vector <double> F2b;

  vector < vector<double> > umn; // actually umn[1..Nmax][1..Mmax]
  vector < vector<double> > vmn; // i.e., "n" is in ROWS

  vector<double>::iterator dit;


  // LAPACK variables
  vector<double> MatA;
  vector<double> B;
  vector<double> pivot;
  int NVARS=0;
  int INFO=0;
  int NRHS=1;


  //double LambdaFactor=0.5*(1.0+(1.0/Lambda))*sqrt(Lambda);

  int Mmax=Nsitesmax%2==0?Nsitesmax/2:(Nsitesmax+1)/2;
  ////////////
  en.clear();  // en0
  chin.clear();   // chin0=0.0

  // Ok, now GSL integration

  double aux[4]={0.0};
  //double Fsq=0.0;


  /* Allocate integration workspace */
  gsl_integration_workspace *w_gsl;
  w_gsl = gsl_integration_workspace_alloc (10000);


  //
  //  Calculate F
  //

  // w^2=Gamma/GammaRef
  // Set GammaRef
  SetHybRef();

  Fsq=CalcF2(w_gsl);
  cout << " F^2 = " << Fsq << endl;
  cout << " GammaRef=" << HybRef << endl;

  cout << " Gamma(-1.0)=" << GetHyb_w(-1.0) << endl;
  cout << " Gamma(0.0)=" << GetHyb_w(0.0) << endl;
  cout << " Gamma(1.0)=" << GetHyb_w(1.0) << endl;

  double IntegGammaNorm=Fsq*HybRef; // gives int^{D}_{-D} Gamma(e) de

  //
  //  Calculate ea, eb, Fas, Fbs
  //

  for (int m=0;m<Mmax;m++){
    CalcFabmEabm(m,w_gsl,'a',&F2a,&ea);
//     cout << " ea_"<<m<<" = " << ea.back()
//  	 << "  (Fa)^2_"<<m<<" = " << F2a.back()/HybRef
//  	 << endl;
    CalcFabmEabm(m,w_gsl,'b',&F2b,&eb);
//     cout << " eb_"<<m<<" = " << eb.back()
//  	 << "  (Fb)^2_"<<m<<" = " << F2b.back()/HybRef
//  	 << endl;
  }


  gsl_integration_workspace_free(w_gsl);


  // Need to redefine Mmax so that if either F2a or F2b is zero, it stops.


  //
  // Calculate u_{0m},v_{0m}
  //

  umn.push_back( vector<double>() );
  vmn.push_back( vector<double>() );
  aux[2]=0.0; //norm
  for (int m=0;m<Mmax;m++){
    aux[0]=sqrt(F2a[m]/IntegGammaNorm);
    aux[1]=sqrt(F2b[m]/IntegGammaNorm);
//     cout << " u_"<<m<<"0 = " << aux[0];
//     cout << " v_"<<m<<"0 = " << aux[1] << endl;

    aux[2]+=(aux[0])*(aux[0])+(aux[1])*(aux[1]);
    umn[0].push_back(aux[0]);
    vmn[0].push_back(aux[1]);
  }

  // Normalize, calculate, en0
  aux[2]=sqrt(aux[2]);
  aux[0]=0.0;

//   cout << " Norm: " << aux[2] << endl; 

  for (int m=0;m<Mmax;m++){
    if (aux[2]>0.0){
       umn[0][m]/=aux[2];
       vmn[0][m]/=aux[2];
//        cout << " unorm_"<<m<<"0 = " << umn[0][m];
//        cout << " vnorm_"<<m<<"0 = " << vmn[0][m] << endl;
    } 
    aux[0]+=(umn[0][m])*(umn[0][m])*ea[m]+(vmn[0][m])*(vmn[0][m])*eb[m];
  }
  //aux[0]*=ScaleFactor(0);
//   cout << " e(n=0) = " << aux[0] << " x "<<ScaleFactor(-1)<<" = " 
//        << aux[0]*ScaleFactor(-1) << endl; 
  aux[0]*=ScaleFactor(-1); // this is the way to do it
  en.push_back(aux[0]);  // en0

  // Calculate chi_0

  aux[0]=0.0;
  for (int m=0;m<Mmax;m++){
//      aux[0]+=pow(RHS(0,ea[m],umn[0][m],0.0,en[0],0.0,false),2.0)
//        +pow(RHS(0,eb[m],vmn[0][m],0.0,en[0],0.0,false),2.0);
    aux[0]+=pow(RHS(0,ea[m],umn[0][m],0.0,en[0],0.0),2.0)
      +pow(RHS(0,eb[m],vmn[0][m],0.0,en[0],0.0),2.0);
  }
  aux[0]=sqrt(aux[0]);
  chin.push_back(aux[0]);
//   cout << " chi(n=0) = " << aux[0] 
//        << "  and tau0 = " << aux[0]/ScaleFactor(0) <<endl;

  //
  // Use recurrence relations to calculate e(n>0), chi(n>0)
  //
  aux[2]=0.0; // norm
  for (int n=1;n<Nsitesmax;n++){
  //for (int n=1;n<=3;n++){
    // 
    // Calculate u_m(n),v_mn using the Chen-Jayaprakash regularization proc.
    //
    umn.push_back( vector<double>() );
    vmn.push_back( vector<double>() );
    // Allocate umn[n], vmn[n]
    umn[n].insert(umn[n].begin(),Mmax,0.0);
    vmn[n].insert(vmn[n].begin(),Mmax,0.0);
    //
    int mmin=n%2==0?n/2:(n-1)/2;
    //cout << " mmin = " << mmin << endl;
    aux[2]=0.0; // en
    for (int m=mmin;m<Mmax;m++){
      if (n==1){
	aux[0]=(RHS(n-1,ea[m],umn[n-1][m],0.0,en[n-1],0.0))/chin[n-1];
	aux[1]=(RHS(n-1,eb[m],vmn[n-1][m],0.0,en[n-1],0.0))/chin[n-1];
      }else{
	aux[0]=(RHS(n-1,ea[m],umn[n-1][m],umn[n-2][m],en[n-1],chin[n-2]))/chin[n-1];
	aux[1]=(RHS(n-1,eb[m],vmn[n-1][m],vmn[n-2][m],en[n-1],chin[n-2]))/chin[n-1];
      }
//       cout << " u_"<<m<<"(n="<<n<<") = " << aux[0];
//       cout << " v_"<<m<<"(n="<<n<<") = " << aux[1];
//       cout << endl;
      umn[n][m]=aux[0];
      vmn[n][m]=aux[1];
      aux[2]+=(aux[0])*(aux[0])+(aux[1])*(aux[1]);
    }
    // loop in m

    //
    // Chen-Jayaprakash procedure
    //
    if (mmin>0){
      //
      // Set system of equations for other umn:
      //
      MatA.clear();
      B.clear();
      pivot.clear();

      // Lookout:
      // for a00 x + a01 y = B0
      //     a10 x + a11 y = B1
      // MatA is in form (a_00,a_10,a_01,a11)!!
      //
      // cols
      for (int ivar2=0;ivar2<2*mmin;ivar2++){
	// rows
	for (int ivar1=0;ivar1<2*mmin;ivar1++){
	  if (ivar2<mmin) 
	    MatA.push_back(umn[ivar1][ivar2]);
	  else
	    MatA.push_back(vmn[ivar1][ivar2-mmin]);
	}
	// Calc B
	aux[0]=0.0;
	for (int m=mmin;m<Mmax;m++){
	  aux[0]+=umn[ivar2][m]*umn[n][m]+vmn[ivar2][m]*vmn[n][m];
	}
	B.push_back(-aux[0]);
      }
      // Solve the system
      NVARS=B.size();
      pivot.insert(pivot.begin(),NVARS,0.0);

//       for (int ii=0;ii<MatA.size();ii++){
// 	cout << "MatA(" << ii <<")="<< MatA[ii] << endl;
//       }
//       for (int ii=0;ii<NVARS;ii++){
// 	cout << "B(" << ii <<")="<< B[ii] << endl;
//       }
      dgesv_(&NVARS, &NRHS, &MatA[0], &NVARS, &pivot[0], &B[0], &NVARS, &INFO);
//       for (int ii=0;ii<NVARS;ii++){
// 	cout << "X(" << ii <<")="<< B[ii] << endl;
//       }
      // Get the other umn, vmn, complete the norm.
      int ii=0;
      for (int m=0;m<mmin;m++){
	umn[n][m]=B[ii];
	vmn[n][m]=B[ii+mmin];
	ii++;
	aux[2]+=(umn[n][m])*(umn[n][m])+(vmn[n][m])*(vmn[n][m]);
      }
    }
    // end Chen-Jayaprakash
    // Renormalize and calculate en
    aux[2]=sqrt(aux[2]); // Norm
    //cout << " Renorm by " << aux[2] << endl; 

    aux[0]=0.0;
    aux[3]=0.0;
    for (int m=0;m<Mmax;m++){
      if (aux[2]>0.0){
	umn[n][m]/=aux[2];
	vmn[n][m]/=aux[2];
// 	cout << " u_"<<m<<"(n="<<n<<") = " << umn[n][m];
// 	cout << " v_"<<m<<"(n="<<n<<") = " << vmn[n][m];
// 	cout << endl;
      } 
      aux[0]+=(umn[n][m])*(umn[n][m])*ea[m]+(vmn[n][m])*(vmn[n][m])*eb[m];
      aux[3]+=(umn[n][m])*(umn[n][m])+(vmn[n][m])*(vmn[n][m]);
    }
    //cout << " New Norm: " << aux[3] << endl; 
    //aux[0]*=ScaleFactor(n);
    aux[0]*=ScaleFactor(n-1);
//     cout << " e(n="<<n<<") = " << aux[0] << endl;
    en.push_back(aux[0]);
    //
    // Calculate tau_n using the new umn
    //
    aux[2]=0.0;
    for (int m=0;m<Mmax;m++){
      aux[0]=(RHS(n,ea[m],umn[n][m],umn[n-1][m],en[n],chin[n-1]));
      aux[1]=(RHS(n,eb[m],vmn[n][m],vmn[n-1][m],en[n],chin[n-1]));
      aux[2]+=(aux[0])*(aux[0])+(aux[1])*(aux[1]);
    }
    aux[2]=sqrt(aux[2]);
//     cout << " chi(n="<<n<<") = " << aux[2]
// 	 << "  and taun = " << aux[2]/ScaleFactor(n) <<endl;
    chin.push_back(aux[2]);
  }
  // end loop in n


  // How to solve a linear system using LAPACK dgesv_ routine
  // Testing: 2 x 2 system: x +3y = 2; x - 3y=4 (x=3, y=-1/3)
  //
  // MatA is TRANSPOSED in relation to the usual: [ 1 1 ; 3 -3 ]
  // rather than "1 3; 1 -3" 
//   MatA.clear();
//   B.clear();
//   MatA.push_back(1.0);
//   MatA.push_back(2.0);
//   MatA.push_back(3.0);
//   MatA.push_back(-3.0);

//   B.push_back(2.0);
//   B.push_back(4.0);
//   NVARS=B.size();

//   for (int ii=0;ii<NVARS;ii++){
//     pivot.push_back(0.0);
//   }

//   dgesv_(&NVARS, &NRHS, &MatA[0], &NVARS, &pivot[0], &B[0], &NVARS, &INFO);

//   for (int ii=0;ii<NVARS;ii++){
//     cout << "B(" << ii <<")="<< B[ii] << endl;
//   }




}
 
//////////////

double CNRGchain::GSL_Integration(gsl_function *pIntegrandGSL,
				  gsl_integration_workspace *w_gsl, 
				  double x0, double x1){

  double Integ=0.0;
  
  /* GSL integration variables */
  double error,epsabs=1.0E-6,epsrel=1.0E-6;
  /* Error handling */
  gsl_error_handler_t *old_handler;
  int status=0;

  old_handler=gsl_set_error_handler_off();
  // Integral int_{-1}^{1} Delta(e))
  status=gsl_integration_qag(pIntegrandGSL,x0,x1,epsabs,epsrel,10000,3,
			     w_gsl,&Integ,&error);
  if (status){
    //      while (status!=0)
    while (status==GSL_EROUND){
      cout << " GSL_Integration: Reducing error tolerance. " << endl;
      epsabs*=10.0;
      epsrel*=10.0;
      status=gsl_integration_qag(pIntegrandGSL,x0,x1,epsabs,epsrel,10000,3,
				 w_gsl,&Integ,&error);
    }
  }
  // end if status

  return(Integ);



}

//////////////

void CNRGchain::SetHybRef(){

  // Normalization factor in w^2=Gamma(e)/Gamma0.
  // Pseudogap cases: HybRef is set in the definition 
  // (see  CNRGCodeHandler::SetChain())

  if(HybRefIsSet){return;}

  gsl_function HybFunGSL;
 

  HybFunGSL.function=HybFunction;
  HybFunGSL.params=HybFuncParams;

  double aux=GSL_FN_EVAL(&HybFunGSL,0.0);
  if (dEqual(aux,0.0)){HybRef=1.0;}else{HybRef=aux;}


}


//////////////

double CNRGchain::CalcF2(gsl_integration_workspace *w_gsl){

  double Integ=1.0;

  gsl_function HybFunGSL;
 
  HybFunGSL.function=HybFunction;
  HybFunGSL.params=HybFuncParams;

  SetHybRef();

  // Consider the case of a sum of deltas
  //
  if (HybFuncIsSumDeltas){
     vector<double>::iterator dit;
     int ic=0;
     double ei=0.0;
     double Ai=0.0;
     cout << " CNRGchain:CalcF2 Params : "<< endl;
     //for (dit=HybFuncParams->begin();dit<HybFuncParams->end();dit++){
     for (dit=HybFuncParams->begin();dit<HybFuncParams->end();dit+=2){
       ei=*dit;
       Ai=*(dit+1);
//        cout << " e_" << ic << " = "<< ei << 
//                "   t_" << ic << " = "<< Ai << endl;
       Integ+=fabs(Ai);
       ic++;
     }
     //exit(0);
   }
   else{
//      Integ=GSL_Integration(&HybFunGSL,w_gsl,-1.0,1.0);
     Integ=GSL_Integration(&HybFunGSL,w_gsl,-(1.0+BandOffset),(1.0-BandOffset));
   }

//   double aux=GSL_FN_EVAL(&HybFunGSL,0.0);
//   cout << "Gamma(0.0) = " << aux << endl;

   //Integ=GSL_Integration(&HybFunGSL,w_gsl,-1.0,1.0);

  //return(Integ);
  return(Integ/HybRef);

}


//////////////

void CNRGchain::CalcFabmEabm(int m,
			     gsl_integration_workspace *w_gsl,
			     char which,
			     vector<double> *pF2m,
			     vector<double> *pEm){
  //
  // Calculates BOTH 
  //  F2a(b)_m = int^{em}_{em+1} Delta(e) de
  //  ea(b)_m = (1/Fa(b)_m)^2 int^{em}_{em+1} e*Delta(e) de 
  //

  double dM=(double)m;
//   double x1=pow(Lambda,-(dM+z_twist-1.0));  
//   double x0=pow(Lambda,-(dM+z_twist));
  // New (June 2013)
  double x1=(1.0-BandOffset)*pow(Lambda,-(dM+z_twist-1.0));  
  double x0=(1.0-BandOffset)*pow(Lambda,-(dM+z_twist));

  // Watch out for dM=0!
  if (m==0) x1=(1.0-BandOffset);


  double Integ=0.0;
  double Integ2=0.0;

  if (which=='b'){
//     Integ=x0; // aux
//     x0=-x1;
//     x1=-Integ;
    x0=-(1.0+BandOffset)*pow(Lambda,-(dM+z_twist-1.0));  
    x1=-(1.0+BandOffset)*pow(Lambda,-(dM+z_twist));
    // Watch out for dM=0!
    if (m==0) x0=-(1.0+BandOffset);
  }
  // end if "b"

//   cout << " dM = " << dM
//        << " z_twist = " << z_twist
//        << " x0 = " << x0 
//        << " x1 = " << x1 
//        << endl;

  gsl_function HybFunGSL;
  gsl_function EnHybF;

  HybFunGSL.function=HybFunction;
  HybFunGSL.params=HybFuncParams;

  // en*HybFun : static, user supplied function!!
  EnHybF.function=HybFuncWithEn;
  EnHybF.params=HybFuncParams;

  Integ=0.0;
  Integ2=0.0;

  if (HybFuncIsSumDeltas){
     vector<double>::iterator dit;
     int ic=0;
     double ei=0.0;
     double Ai=0.0;
     // Ok, let's do the whole thing (can't assume energies are ordered)
     for (dit=HybFuncParams->begin();dit<HybFuncParams->end();dit+=2){
       ei=*dit;
       Ai=*(dit+1);
       
       if ((ei>=x0)&&(ei<=x1) ){

	 Integ+=fabs(Ai);
	 Integ2+=fabs(Ai)*ei;
	 cout << " e_" << ic << " = "<< ei << 
	   "   t_" << ic << " = "<< Ai << endl;

       }
       ic++;
     }
     //exit(0);
   }
   else{
     Integ=GSL_Integration(&HybFunGSL,w_gsl,x0,x1);
     
     Integ2=GSL_Integration(&EnHybF,w_gsl,x0,x1);
   }


//   cout << " Integ = " << Integ <<
//     " Integ2 = " << Integ2 << endl;


  pF2m->push_back(Integ);

// Note: for the Campo-Oliveira, this will be different!!
//     pEm->push_back(Integ/Integ2); 

  switch (DiscScheme){
  case 0:
    if (Integ!=0.0)
      pEm->push_back(Integ2/Integ);
    else
      pEm->push_back(-1.0);
    break;
  case 1: // Campo-Oliveira discretization scheme
    if (Integ2!=0.0)
      pEm->push_back(Integ/Integ2);
    else
      pEm->push_back(-1.0);
    break;
  default:
    if (Integ!=0.0)
      pEm->push_back(Integ2/Integ);
    else
      pEm->push_back(-1.0);
    break;
  }
  // end DiscScheme switch


}


//////////////

double CNRGchain::GetHyb_w(double omega){

  // Returns the GSL-generated Hybridization function Hyb(omega)

  gsl_function HybFunGSL;
 
  HybFunGSL.function=HybFunction;
  HybFunGSL.params=HybFuncParams;

  return(GSL_FN_EVAL(&HybFunGSL,omega));

}



//////////////









// Ok, now I quit:
// No way to pass a non-static method that as a pointer to GSL
// AND use dynamic method HybFunction at the same time !!!
// code below DOES NOT WORK:
//
// double CNRGchain::HFtimesEn(double omega, 
// 			     void *params){

//   vector<double> VecParam=*(vector<double> *)(params);

//   cout << " First param: " << VecParam[0] << endl;
//   double aux=HybFunction(omega,params);

//   return(aux*omega);

// }

// double CNRGchain::HFtimesEnCallBack(double omega,
// 				    void* params){
//   //NRGchainGSLCallbackHolder* h =
//   //  static_cast<NRGchainGSLCallbackHolder*>(params); 
//   //vector<double> VecParam=*(vector<double> *)(h->data);


//   NRGchainGSLCallbackHolder h;
//   h.data=params;
//   vector<double> VecParam=*(vector<double> *)(h.data);
//   cout << " First param HERE: " << VecParam[0] << endl;
//   h.HybFunction=HybFunction;
//   double aux=h.cls.HFtimesEn(omega,h.data);

//   return(aux);
//   //return h->cls->HFtimesEn(omega,h->data);
// }
//
// As long as we are using GSL, we'll write static functions.
//

////////////////////
