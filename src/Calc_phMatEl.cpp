
#include <iostream>
#include <math.h>

#include <gsl/gsl_sf_gamma.h>
#include "NRGfunctions.hpp"

using namespace std;



double Calc_phMatEl(double mi, double mj, double lambda, double w0){

  //
  // Calculates the following element:
  //
  // <mi| Exp(-(lambda/w0)(a^+ - a) |mj>
  //
  //

  int imi=(int)mi;
  int imj=(int)mj;

  if (dEqual(w0,0.0)) return(0.0);

  if ( (mi<0.0)||(mj<0.0) ) return(0.0);


  double loverw0=lambda/w0;
  int n=0;
  int Dm=imi-imj;

  int nmin=(Dm>0?0:abs(Dm));

  double Expfac=exp(-0.5*loverw0*loverw0);
  double termo=0.0;
  double aux[2];

//   cout << " nmin = " << nmin
//        << " imj = " << imj
//        << " Dm  = " << Dm
//        << " loverw0 = " << loverw0
//        << " Expfac = " << Expfac
//        << endl;

  for (n=nmin;n<=imj;n++)
    {
      aux[0]=gsl_sf_fact(n)*gsl_sf_fact(n+Dm)*gsl_sf_fact(imj-n);
      aux[1]=sqrt(gsl_sf_fact(imj)*gsl_sf_fact(imj+Dm));
      termo+=pow(-loverw0,2*n+Dm)*pow(-1.0,n)*aux[1]/aux[0];
    }
  //end loop in n

  termo*=Expfac;


  return(termo);

}
///////////////////////////////////////

double Calc_apad(double mi, double mj){

  //
  // Calculates the following element:
  //
  // <mi| (a^+ + a) |mj> = sqrt(mj) delta_{mi+1 mj} + sqrt(mj+1)delta_{mi-1 mj}
  //
  //

  int imi=(int)mi;
  int imj=(int)mj;

  if (imi==imj-1) return(sqrt(mj));
  else
    if (imi==imj+1) return(sqrt(mj+1.0));
    else
      return(0.0);


}
//////////////////////////
