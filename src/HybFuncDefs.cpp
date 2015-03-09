
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
