#include <iostream>
#include <vector>
#include <math.h>

#include <gsl/gsl_sf_gamma.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"


// Calculates the Clebsh-Gordan coefficient:
//
//
// <S S; S~ Sz~=S-ST) | ST ST >  
//
//
// Uses gsl factorial function

double CGordan_SStot(double S, double Stilde, double Sztilde, double Stot){

  double aux1,aux2;

  if ( (fabs(S+Sztilde-Stot)>1.0e-10)||(fabs(Sztilde)>Stilde) ) return(0.0);
  if ( (S<0.0)||(Stilde<0.0)||(Stot<0.0) ) return(0.0);

  cout << " S = " << S << " S~ = " << Stilde << " Sztilde = " << Sztilde << " Stot = " << Stot << endl;
  // Numerator
  aux1=gsl_sf_fact((int)(2*S-1))*gsl_sf_fact((int)(2*Stot-1));
  aux1*=Stot*(2*Stot+1)*S;
  aux1=2*sqrt(aux1);

  // Denominator
  aux2=gsl_sf_fact((int)(S+Stot-Stilde-1))*gsl_sf_fact((int)(S+Stot+Stilde-1));
  aux2*=(S+Stot-Stilde)*(S+Stot+Stilde)*(S+Stot+Stilde+1);
  aux2=sqrt(aux2);

  cout << "aux1 = " << aux1 << " aux2 = " << aux2 << endl;
  cout << "CG = " << aux1/aux2 << endl;

  return(aux1/aux2);

}
