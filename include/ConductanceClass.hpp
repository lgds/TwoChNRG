
#include "NRGclasses.hpp"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
//#include <gsl/gsl_complex_math.h>
#include "SpecFuncClass.hpp"
#include "NRGfunctions.hpp"

#ifndef _CCONDUCTANCE_
#define _CCONDUCTANCE_

class CConductance{

public:

  // Default constructor (with an initialization list!)
  CConductance(): ModelNo(0),SymNo(0),BandNo(0),Nops(1),Temp(0.0),Mtemp(1001),UseSpec(0),CalcResistivity(false),SpinResolved(false),IsGraphene(false)
  {
  }

  // Default destructor
  ~CConductance(){}

  // Member variables


  vector <CSpecFunction> SpecVec;

  //Set from ThisCode
  int ModelNo;
  int SymNo;
  int BandNo;
  int Nops;
  vector <double> CondParams;

  double Temp;
  int Mtemp;

  // To use one of the spec functions only

  int UseSpec;
  // 0 - Use all spectral functions in the calculation
  // n>0 - Use spectral functions SpecVec[n-1].

  bool CalcResistivity;

  bool SpinResolved;

  bool IsGraphene;


  // Member functions

  void Initialize ();

  void ReorderParams();

  void SetIntegrand();

  void PrintParams();

  double CalcIntegral(int WhichIntegrand=0);

  // Pointers to functions (need to think about this...)
  // Integrand in GSL format (needs to be defined externally)
  double (*Integrand)(double omega,
		       void *params);

  // If spin resolved
  double (*Integrand2)(double omega,
		       void *params);


  // TMatrix (needs to be defined externally)
  gsl_complex (*TMatrix)(double omega,
			 void *params);


  double FermiFunction(double En, double Temp);


};
// Never forget the ; at the end of a class declaration!

#endif


////////////////////////////////////////////////
///         Tmatrix definitions             ////
///                                         ////
///   Static functions (model dependenent)  ////
///                                         ////
////////////////////////////////////////////////

#ifndef _TMATRIXANDERSON_
#define _TMATRIXANDERSON_

gsl_complex TMatrix_Anderson(double omega,
			     void *params);

double Integrand_Anderson(double omega,
			  void *params);

double Integrand_Anderson_Resistivity(double omega,
				      void *params);

double Integrand_Anderson_Resistivity_Graphene(double omega,
					       void *params);


#endif

#ifndef _TMATRIXDQD_
#define _TMATRIXDQD_

gsl_complex TMatrix_DQD(double omega,
			void *params, int UpDn=0);

double Integrand_DQD(double omega,
		     void *params);

double Integrand_DQD_up(double omega,
		     void *params);

double Integrand_DQD_dn(double omega,
		     void *params);


#endif


#ifndef _TMATRIXSIDEDOT_
#define _TMATRIXSIDEDOT_

gsl_complex TMatrix_SideDot(double omega,
			    void *params);

double Integrand_SideDot(double omega,
			 void *params);


#endif

#ifndef _TMATRIXCAVITY_
#define _TMATRIXCAVITY_

gsl_complex TMatrix_Cavity(double omega,
			   void *params);

double Integrand_Cavity(double omega,
			 void *params);


#endif
