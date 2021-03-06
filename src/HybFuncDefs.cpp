
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)


#include "NRGclasses.hpp"
#include "SpecFuncClass.hpp"
using namespace std;

//
// Static definitions for Hynridization functions to be used with 
// CNRGchain objects
//
//  For each HybFunc, a second function ("HybFuncWithEn")
//  need to be defined. It can be either Delta(en)*en 
//  (usual) or something a la Campo-Oliveira 
//  
//

// VecParam[0] is ALWAYS whichbandtype!!
// for the "usual" bandtypes, the order in VecParam is
//
//   0- whichbandtype
//   1- e2
//   2 - Gamma2
//   3 - Gamma1
//   4 - small_lambda
//   5 - MagFlux

////////////////////////////////////////////////
double HybFunc_Lorentzian(double omega,
			  void *params){
  // A*omega^2

  vector<double> VecParam=*(vector<double> *)(params);

  //Gamma(e)=pi*lambda^2*(Gamma2/Pi)/[(en-e2)^2+Gamma2^2]
  // pure Lorentzian shape (peak)
  // e2 in VecParam[1]
  // Gamma2 in VecParam[2]
  // lambda in VecParam[4]

  double e2=VecParam[1];
  double Gamma2=VecParam[2];
  double small_lambda=VecParam[4];

//    cout << " e2= " << e2 
//         << " Gamma2= " << Gamma2 
//         << " lambda = " << small_lambda << endl;

  if (VecParam.size()<6)
    return(1.0);
  else{
    double aux=pow((omega-e2),2.0)+pow(Gamma2,2.0);
    if (aux!=0.0)
      aux=small_lambda*small_lambda*Gamma2/aux;
    return(aux);
  }

}

///////

double HybFunc_Lorentzian_timesEn(double omega,
				  void *params){

  return(omega*HybFunc_Lorentzian(omega,params));

}

////////////////////////////////////////////////

double HybFunc_Const(double omega,
		     void *params){
  // Constant Densty of States (usual case)
  // Delta(e)=Gamma0
  // 
  vector<double> VecParam=*(vector<double> *)(params);

  double Gamma0=1.0;
  if (VecParam.size()==1)
    Gamma0=VecParam[0];
  
  return(Gamma0);


}

///////

double HybFunc_Const_timesEn(double omega,
			     void *params){

  return(omega*HybFunc_Const(omega,params));

}

///////

double HybFunc_Const_divEn(double omega,
			     void *params){

  if (omega!=0.0)
    return(HybFunc_Const(omega,params)/omega);
  else
    return(-1.0);


}

//////////////////

////////////////////////////////////////////////
double HybFunc_PowerLaw(double omega,
			void *params){
  // Pseudogapped : Delta(w)=Gamma0*|w-w0|^r 
  // or Gapped (Delta(w)=Gamma0*|w-w0|^r for |w0|>Delta and zero otherwise)
  // Actually, generalizing to (2014):
  // Delta(w)=Gamma0*|w-w0|^r + small_gamma

  vector<double> VecParam=*(vector<double> *)(params);


// VecParam[0] is ALWAYS whichbandtype!!
// Here, we have
//   0- whichbandtype
//   1- r_exp
//   2 - Gamma0
//   3 - Delta (small_gamma)
//   4 - w0
//   5 - 


  double r_exp=VecParam[1];
  double Gamma0=VecParam[2];
//   double Delta=VecParam[3];
  double small_gamma=VecParam[3];
  double w0=VecParam[4];

  double aux=1.0;

//    cout << " r_exp= " << r_exp 
//         << " Gamma0= " << Gamma0 
//         << " Delta = " << Delta << endl;

  if (VecParam.size()<5)
    return(1.0);
  else{
//     if (abs(omega-w0)>=Delta)
//       aux=Gamma0*pow(abs(omega-w0),r_exp);
//     else
//       aux=0.0;
    aux=Gamma0*pow(abs(omega-w0),r_exp)+small_gamma;
    
    return(aux);
  }

}

///////

double HybFunc_PowerLaw_timesEn(double omega,
				  void *params){

  return(omega*HybFunc_PowerLaw(omega,params));

}


///////

double HybFunc_PowerLaw_divEn(double omega,
				  void *params){

  if (omega!=0.0)
    return(HybFunc_PowerLaw(omega,params)/omega);
  else
    return(-1.0);



}


////////////////////////////////////////////////

////////////////////////////////////////////////
double HybFunc_FromFile(double omega,
			void *params){
  //
  // Takes the data from file HybFunc.dat 
  // and interpolates it
  // 


  // Ok, params needs to be a pointer to a double STL vector with the data!!

  // I think it's doable. Modify NRGcodehandlerClass.cpp -> SetChain

  // Need to write chain.ReadParams_new to accomodate this change.
  // Meaning: reads into a vector < vector<double> >
  // Then it should be ok!

  // DONE! Need to pass 

  // vector<double> VecParam=*(vector<double> *)(params);
  // Needs to be a double array

  // This contains the WHOLE THING!!
  // It SHOULD NOT be a POINTER!! Should allocate this WITHIN THE OBJECT!
  vector<double> VecParam= *(vector<double> *)(params);

  // First 2N numbers : w,Dos(w)
  // LAST NUMBER: Nlines.

  double aux=1.234;

  CSpecFunction spec_aux;
  // To use RhoInterpol
  // Assuming 2 cols.

  // Read data file: omega_data, HybFunc_data

  spec_aux.Omega_Rhow.push_back( vector<double> () ); // omegas
  spec_aux.Omega_Rhow.push_back( vector<double> () ); // HybFunc

  for (int ii=0;ii<VecParam.size()-1;ii+=2){
    spec_aux.Omega_Rhow[0].push_back(VecParam[ii]);
    spec_aux.Omega_Rhow[1].push_back(VecParam[ii+1]); 
  }
  // end print files

  // Debug
  //spec_aux.PrintOmegaRhow();

  // Sort data (if not sorted).

  // Interpolated it. Need a function similar to RhoInterpol to do it.
  // Using RhoInterpol... WORKS!

  aux=spec_aux.RhoInterpol(omega);


  return(aux);
}

///////

double HybFunc_FromFile_timesEn(double omega,
				  void *params){

  return(omega*HybFunc_FromFile(omega,params));

}

///////

double HybFunc_FromFile_divEn(double omega,
			      void *params){

  if (omega!=0.0)
    return(HybFunc_FromFile(omega,params)/omega);
  else
    return(-1.0);

}


////////////////////////////////////////////////
double HybFunc_Cavity(double omega,
			  void *params){
 

  vector<double> VecParam=*(vector<double> *)(params);

  // Parameters
  //
  // dot-leads: \Gamma_{dS}=\Gamma{dR}=0.125U=0.0625D$.
  // dot-cavity: \small_lambda = 0.06U--0.1U=0.03D--0.05D$.
  // Cavity-reservoir: \Gamma_{cR}= 0.06U = 0.03D$.
  // Level spacing: \delta = 0.32 U = 0.016D$.
  // cavity level e2
  // Nlevels 

  // Gamma(e)=-IMAG(Sigma(e))
  // where
  //
  // Sigma(e)=-i(Gamma_{dS}+Gamma_{dR})
  //          +(\small_lambda-i(sqrt{Gamma_{dR}*Gamma_{cR}}))^2
  //          x S(e)/(1+iS(e)Gamma_{cR})
  // 
  // S(e) =sum_{j=0}^{Nlevels-1} [e-(e2+j*delta)+I*small_gamma]^-1
  // 
  // let's set small_gamma=1e-10.
  // 
  //

  // cavity level e2 in VecParam[1]
  // Cavity-reservoir: \Gamma_{cR} in VecParam[2]
  // dot-leads: \Gamma_{dS}=\Gamma{dR} in VecParam[3]
  // dot-cavity: \small_lambda in VecParam[4]
  // Level spacing: \delta in VecParam[5]
  // Nlevels in VecParam[6]


  double e2=VecParam[1];
  double GammacR=VecParam[2];
  double GammadR=VecParam[3];
//  double GammadS=GammadR; // not anymore
  double small_lambda=VecParam[4];
  double delta=VecParam[5];
  int Nlevels=(int)VecParam[6];
  double GammadS=VecParam[7];

  double small_gamma=1e-10;
 

//    cout << " e2= " << e2 
//         << " GammacR= " << GammacR 
//         << " GammadS(R)= " << GammadS 
//         << " small_lambda = " << small_lambda 
//         << " delta= " << delta
//         << " Nlevels= " << Nlevels
// 	<< endl;

  if (VecParam.size()<7)
    return(1.0);
  else{
    // Set S(w)
    gsl_complex Sw;
    gsl_complex caux[3];
    double aux=0.0;
    GSL_SET_COMPLEX(&Sw,0.0,0.0);

    for (int ii=0;ii<Nlevels;ii++){
      aux=e2+ii*delta;
      GSL_SET_COMPLEX(&caux[0],0.0,small_gamma);
      caux[0]=gsl_complex_add_real(caux[0],(omega-aux));
      caux[0]=gsl_complex_inverse(caux[0]);
      Sw=gsl_complex_add(Sw,caux[0]);
    }
    // Set Sigma(w)
  // Sigma(e)=-i(Gamma_{dS}+Gamma_{dR})
  //          +(\small_lambda-i(sqrt{Gamma_{dR}*Gamma_{cR}}))^2
  //          x S(e)/(1+iS(e)Gamma_{cR})

    GSL_SET_COMPLEX(&caux[0],0.0,-1.0);
    caux[0]=gsl_complex_mul_real(caux[0],(GammadS+GammadR));

    GSL_SET_COMPLEX(&caux[1],small_lambda,-sqrt(GammadR*GammacR));
    caux[1]=gsl_complex_mul(caux[1],caux[1]);

    GSL_SET_COMPLEX(&caux[2],0.0,1.0);
    caux[2]=gsl_complex_mul(caux[2],Sw);
    caux[2]=gsl_complex_mul_real(caux[2],GammacR);
    caux[2]=gsl_complex_add_real(caux[2],1.0);
    caux[2]=gsl_complex_div(Sw,caux[2]);
    
    caux[1]=gsl_complex_mul(caux[1],caux[2]);
    
    //Sigma
    caux[0]=gsl_complex_add(caux[0],caux[1]);

    // Set Delta(w)    
    aux=-GSL_IMAG(caux[0]);

    return(aux);
  }

}

///////

double HybFunc_Cavity_timesEn(double omega,
			      void *params){

  return(omega*HybFunc_Cavity(omega,params));

}

////////////////////////////////////////////////

////////////////////////////////////////////////

// Sign function
double SignF(double xin){
  double Y = 1.0;
  if (xin < 0.0){
    Y = -1.0;
  }
  return Y;
}

// Heaviside theta function
double HeavF(double xin){
  double Y = 0.0;
    if (xin >= 0.0){
    Y = 1.0;
  }
  
  return Y; 
}


double HybFunc_CoSilicene(double omega,
			  void *params){
 
  vector<double> VecParam=*(vector<double> *)(params);

  // Parameters
  //
  // params contains alpha, beta, vF,
  // lambda_SO, lEz, D, in that order

  // VecParam[0] is ALWAYS whichbandtype!!
  
  // params[1] = alpha =  0.28*22.8 // units of eV*Angstrom
  // params[2] = beta = -0.489     // units of eV
  // params[3] = vF =  22.8      // units of eV*Angstrom
  // params[4] = lambda_SO = 0.0039    // units of eV
  // params[5] = lEz =  0.00      // units of eV, TUNABLE,
  //                      // input as fracc. of D=22.5eV
  // params[6] = 2*ich - 1 // ich=0->dchi=-1 ich=1->dchi=+1

  // All energy parameters in units of D =  22.5 eV
  
  double Gamma=0.0;

  
  if (VecParam.size()<7)
    return(1.0);
  
  

  double alpha = VecParam[1];
  double beta  = VecParam[2];
  double vF    = VecParam[3];
  double lSO   = VecParam[4];
  double lEz   = VecParam[5];

  double dchi  = VecParam[6]; // equals +/- 1.0 --> depends on channel


  double V0Sq= (alpha/vF)*(beta)*pow((lEz - dchi*lSO),2.0);
  double V1Sq = pow(beta,2.0) - (pow(alpha/vF,2.0))*(pow(lEz - dchi*lSO,2.0));
  double V2Sq = -(alpha/vF)*beta;
  double V3Sq = pow(alpha/vF,2.0);

  
  double Gamma0=SignF(omega)*V0Sq;
  double Gamma1=fabs(omega)*V1Sq;
  double Gamma2=SignF(omega)*(pow(omega,2.0))*V2Sq;
  double Gamma3=fabs(pow(omega,3.0))*V3Sq;
 

  Gamma  = (2.0*M_PI)*(Gamma0 + Gamma1 + Gamma2 + Gamma3)*HeavF(1.0-fabs(omega))*HeavF(fabs(omega)-fabs(lEz - dchi*lSO));

  return(Gamma);
  
  // enf 
 
}

///////

double HybFunc_CoSilicene_timesEn(double omega,
				  void *params){

  return(omega*HybFunc_CoSilicene(omega,params));

}

////////////////////////////////////////////////
