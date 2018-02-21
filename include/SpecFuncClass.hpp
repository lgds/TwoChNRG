
#include "NRGclasses.hpp"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


#ifndef _CSPECFUNCTION_
#define _CSPECFUNCTION_

class CSpecFunction{

public:

  // Default constructor (with an initialization list!)
  CSpecFunction(): Lambda(1.0),NshellMax(1),NshellMin(0),NonDiagGF(false),dBroad(0.3),Temp(0.0),Mtemp(1000),Betabar(0.727),dBroadTemp(0.1),TwindowFac(1.0),IsBosonic(false),z_twist(1.0),Gap(-1.0)
  {
  }

  // Default destructor
  ~CSpecFunction(){}

  // Member variables

  double Lambda;
  double z_twist;
  int NshellMax;
  int NshellMin;

  bool NonDiagGF;
  bool IsBosonic;

  // Temperature (Feb 2011)
  double Temp; // ="Real Temp" (in units of D)=DN(Mtemp)/betabar
  int Mtemp;
  double Betabar;

  // General broadening parameter

  double dBroad;

  // Finite-temperature parameters


  double dBroadTemp;
  double TwindowFac;

  // Adding Gap (no broadening)

  double Gap;

  // Arrays of NshellMax objects
  CNRGbasisarray* AcutN;
  CNRGmatrix* RhoN;
  CNRGmatrix* Op1N;
  CNRGmatrix* Op2N;

  // Pairs of (omega, rho) at fixed WN values 

  vector < vector<double> > Omega_Rhow_Even; // two ROWS (omega, rhow)
  vector < vector<double> > Omega_Rhow_Odd;  // two ROWS
  vector < vector<double> > Omega_Rhow;      // two ROWS

  // Name

  char Name[48];


  // Member functions

  double HalfLambdaFactor();

  double CalcDN(int Nshell);

  int CalcNfromOmega(double omega);

  double SumDeltas(double bbroad,
		   vector<double> Msq, 
		   vector<double> En,
		   double omega);

  vector < vector<double> > CalcSpecCosti_Teq0();

  vector < vector<double> >  CalcSpecCosti_Teq0(double bbroad,
			    double factorWN);


  vector < vector<double> >  CalcSpecCosti_T(int Mtemp);


  // Bulla's N,N+2 method (TO BE DONE)
  void SpecBulla_SetMsqEdiff(int Mtemp);

  double CalcSpecBulla_T(int Mtemp, double omega);

  void CalcSpecBulla_T_FixedOmegas(int Mtemp, double factorWN);


  //double CalcSpecCosti_T(int Mtemp);

  double CalcSpecCosti_T_N(int Nshell, double betabar,
			   double factorWN, bool CalcNorm=false);


  double DMNRG_SpecDens_M(double omega,int Nshell, bool CalcNorm=false);

  // Doesn't work yet... Double counting. Need FULL basis.
  double CalcSpecDM_NRG(double omega, int UseCFS=0);
  


  void CalcSpecDM_NRG_FixedOmegas(double factorWN, int UseCFS=0,int NwEachN=1);

  void CalcSpec_ManyOmegas(int NomegasEachN, double factorWN, int UseCFS);

  
  double CalcSpecDM_NRG_InterpolOddEven(double omega);

  void DMNRG_SpecDens_ChkPHS(int Nshell);


  // Complete Fock Space

  double CFS_SpecDens_M(double omegabar,int Nshell, bool CalcNorm=false);

  // BroadDeltaFunction

  double (*BDelta)(double omega, double Ei, double bbroad);

  double (*BDeltaTemp)(double omega, double Ei, double bbroad);


  // Aij*BDelta(Ei-Ej);
  boost::numeric::ublas::matrix<double> MijxBDeltaEij(boost::numeric::ublas::matrix<double>Mat, CNRGarray* pAeig, int iblock1, int iblock2, double omega);
  boost::numeric::ublas::matrix<double> MijxBDeltaEij(CNRGmatrix* pMat, CNRGarray* pAeig, int iblock1, int iblock2, double omega);

  // complex
  boost::numeric::ublas::matrix<complex<double> > cMijxBDeltaEij(boost::numeric::ublas::matrix<complex<double> >Mat, CNRGarray* pAeig, int iblock1, int iblock2, double omega);
  boost::numeric::ublas::matrix<complex<double> > cMijxBDeltaEij(CNRGmatrix* pMat, CNRGarray* pAeig, int iblock1, int iblock2, double omega);


  // Aij*BDelta(Ei-Ej);
  boost::numeric::ublas::matrix<double> MijxBDeltaEij(boost::numeric::ublas::matrix<double> Mat, CNRGarray* pAeig, int iblock1, int iblock2, bool kp1, bool kp2, double omega);
  boost::numeric::ublas::matrix<double> MijxBDeltaEij(CNRGmatrix* pMat, CNRGarray* pAeig, int iblock1, int iblock2, bool kp1, bool kp2, double omega);
  // complex
  boost::numeric::ublas::matrix<complex<double> > cMijxBDeltaEij(boost::numeric::ublas::matrix<complex<double> > Mat, CNRGarray* pAeig, int iblock1, int iblock2, bool kp1, bool kp2, double omega);
  boost::numeric::ublas::matrix<complex<double> > cMijxBDeltaEij(CNRGmatrix* pMat, CNRGarray* pAeig, int iblock1, int iblock2, bool kp1, bool kp2, double omega);



  // KK transformation

  double RhoInterpol( double omega );

  static double RhoInterpol2GSL( double omega, void *params);

  double KKtransf( gsl_function Rhow, double x0, double x1, double omega );

  double KKRho(double omega );

  // Check Normalization

  double CalcNorm(int UseCFS=1);

  double CalcNormInteg();


  // Green's function

  gsl_complex GreensFunction (double omega);

  // Save/Read Interpol

  void SaveOmegaRhow();

//   void ReadOmegaRhow();
  void ReadOmegaRhow(char arqextension[]=(char*)"_OmegaRhow.dat");

  void PrintOmegaRhow();

  void ClearOmegaRhow(){
    Omega_Rhow_Even.clear();
    Omega_Rhow_Odd.clear();
    Omega_Rhow.clear();
  }

  void PrintOmegaEven();

  void GetSubGapData(int Nshell, 
		     vector< double > &Eb, 
		     vector< double > &wb,
		     bool NegOmega=false);

};
// Never forget the ; at the end of a class declaration!

#endif
