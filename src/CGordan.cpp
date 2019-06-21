#include <iostream>
//#include <vector>
#include <algorithm>
#include <math.h>

#include <gsl/gsl_sf_gamma.h>
using namespace std;

//
// Triangular Inequality
//

bool TriangIneq(double J1, double J2, double J){

  // J1+J2>= J >= |J1-J2| -> true
  // J>J1+J2 or J<|J1-J2| -> false

  if ( (J-(J1+J2)>1.0e-10)||(J-fabs(J1-J2)<-1.0e-10) ) return(false);
  else return(true);


}



// Calculates the Clebsh-Gordan coefficient:
//
//
// <S S; S~ Sz~=S-ST) | ST ST >  
//
//
// Uses gsl factorial function

double CGordan(double J1, double M1, double J2, double M2, double J, double M){

  // As defined in Landau-Lifshitz' "Quantum Mechanics"

  double aux[8];
  
  double rho,sigma,tau;
  int rmax,rmin;
  
  // Safeguards

  if ( (J1<0.0)||(J2<0.0)||(J<0.0) ) return(0.0);
  if ( (fabs(M1)>J1)||(fabs(M2)>J2)||(fabs(M)>J) ) return(0.0);
  if ( fabs(M-(M1+M2))>1.0e-10 ) return(0.0);

  // Triangular inequality:

  if ( (J>(J1+J2))||(J<fabs(J1-J2)) ) return(0.0);

  //cout << scientific;
  //cout << "J1 = " << J1 << " M1 = " << M1 << endl;
  //cout << "J2 = " << J2 << " M2 = " << M2 << endl;
  //cout << "J  = " << J  << " M  = " << M  << endl;

  // Calculate rho

  aux[0]=gsl_sf_fact((int)(J1+J2-J));
  aux[1]=gsl_sf_fact((int)(J1-J2+J));
  aux[2]=gsl_sf_fact((int)(-J1+J2+J));
  aux[3]=2*J+1;
  aux[4]=gsl_sf_fact((int)(J1+J2+J+1));

  aux[5]=aux[0]*aux[1]*aux[2]*aux[3]/aux[4];

  if (aux[5]>0.0) rho=sqrt(aux[5]);
  else rho=0.0;

  //cout << "rho = " << rho << endl; 

  // Calculate sigma

  aux[0]=gsl_sf_fact((int)(J+M));
  aux[1]=gsl_sf_fact((int)(J-M));
  aux[2]=gsl_sf_fact((int)(J1+M1));
  aux[3]=gsl_sf_fact((int)(J1-M1));
  aux[4]=gsl_sf_fact((int)(J2+M2));
  aux[5]=gsl_sf_fact((int)(J2-M2));

  aux[6]=aux[0]*aux[1]*aux[2]*aux[3]*aux[4]*aux[5];

  if (aux[6]>0.0) sigma=sqrt(aux[6]);
  else sigma=0.0;

  //cout << "sigma = " << sigma << endl; 


  // Calculate tau
  // rmin = Max -(j-j2+m1) -(j-j2+m1)
  // rmax = Min (j1-m1), (j2+m2), j1+j2-j
  int a1[]={(int)(J1-M1),(int)(J2+M2),(int)(J1+J2-J)};
  rmax=*min_element(a1,a1+3);
  int a2[]={-(int)(J-J2+M1),-(int)(J-J1-M2),0};
  rmin=*max_element(a2,a2+3);

  //cout << "rmin = " << rmin << " "; 
  //cout << "rmax = " << rmax << endl; 


  tau=0.0;
  for (int r=rmin;r<=rmax;r++)
    {
      aux[0]=gsl_sf_fact((int)(J1-M1)-r);
      aux[1]=gsl_sf_fact((int)(J2+M2)-r);
      aux[2]=gsl_sf_fact((int)(J-J2+M1)+r);
      aux[3]=gsl_sf_fact((int)(J-J1-M2)+r);
      aux[4]=gsl_sf_fact((int)(J1+J2-J)-r);
      aux[5]=gsl_sf_fact(r);
      aux[6]=1.0/(aux[0]*aux[1]*aux[2]*aux[3]*aux[4]*aux[5]);
      tau+=pow(-1.0,r)*aux[6];
    }

  //cout << "tau = " << tau << endl; 


  return(rho*sigma*tau);

}
