#include <iostream>
//#include <vector>
//#include <algorithm>
#include <cmath>
#include "NRGfunctions.hpp"

using namespace std;

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////


double Sdots_totalS(double ST, double STz, double Sold, double Stilde){

// Calculates the matrix element
//
//
// <ST STz| S.s| ST STz >  
//
// where 
//
// | ST STz > = \sum_{Sz Sztilde} CG(...) |Sold >x|Stilde>
//
// Uses CGordan.cpp

  double auxCG[3]={0.0,0.0,0.0};
  double fac[3]={0.0,0.0,0.0};
  double aux=0.0;

  for (double Sztilde=(-Stilde);Sztilde<=Stilde;Sztilde+=1.0)
    {
      double Szold=STz-Sztilde;
//       cout << "ST = " << ST << " STz = " << STz << endl;
//       cout << "Stilde = " << Stilde << " Sztilde = " << Sztilde << endl;
//       cout << " Sold = " << Sold << " Szold = " << Szold << endl;
      if (Szold<=fabs(Sold))
	{
	  auxCG[0]=CGordan(Sold, Szold, Stilde, Sztilde, ST, STz);
	  auxCG[1]=CGordan(Sold, Szold+1.0, Stilde, Sztilde-1.0, ST, STz);
	  auxCG[2]=CGordan(Sold, Szold-1.0, Stilde, Sztilde+1.0, ST, STz);

	  double Soldp1=Sold*(Sold+1.0);
	  double Stildep1=Stilde*(Stilde+1.0);
	  fac[0]=0.5*sqrt(Soldp1-Szold*(Szold+1.0))*
	    sqrt(Stildep1-Sztilde*(Sztilde-1.0));
	  fac[1]=0.5*sqrt(Soldp1-Szold*(Szold-1.0))*
	    sqrt(Stildep1-Sztilde*(Sztilde+1.0));
	  fac[2]=Szold*Sztilde;

	  // First term
	  aux+=auxCG[0]*auxCG[1]*fac[0];
	  // Second term
	  aux+=auxCG[0]*auxCG[2]*fac[1];
	  // Third term
	  aux+=auxCG[0]*auxCG[0]*fac[2];
	}
    }


  return(aux);
}


//////////////////////////////////
//////////////////////////////////
//////////////////////////////////


double Sdots_Sz(double Sold, double Stilde,
		double Szoldp,  double Sztildep,
		double Szold,  double Sztilde){

// Calculates the matrix element
//
//
// <Sold Szoldp|x<Stilde Sztilde| S.s | Sold Szold>x| Stilde Sztilde >  
//

  if ( !dEqual((Szoldp+Sztildep),(Szold+Sztilde)) ) return(0.0);
  if (  (fabs(Szoldp)>Sold)||
        (fabs(Szold)>Sold)||
        (fabs(Sztildep)>Stilde)||
        (fabs(Sztilde)>Stilde)
	) return(0.0);

  double aux=0.0;
  
  double Soldp1=Sold*(Sold+1.0);
  double Stildep1=Stilde*(Stilde+1.0);
  double fac[3]={0.0,0.0,0.0};

  fac[0]=0.5*sqrt(Soldp1-Szold*(Szold+1.0))*
    sqrt(Stildep1-Sztilde*(Sztilde-1.0));
  fac[1]=0.5*sqrt(Soldp1-Szold*(Szold-1.0))*
    sqrt(Stildep1-Sztilde*(Sztilde+1.0));
  fac[2]=Szold*Sztilde;


  // Term S+ s-

  if (dEqual(Szoldp,Szold+1.0)&&dEqual(Sztildep,Sztilde-1.0) ) aux+=fac[0];

  // Term S- s+

  if (dEqual(Szoldp,Szold-1.0)&&dEqual(Sztildep,Sztilde+1.0) ) aux+=fac[1];

  // Term Sz sz

  if (dEqual(Szoldp,Szold)&&dEqual(Sztildep,Sztilde) ) aux+=fac[2];

  return(aux);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

double Splus(double Sp, double Szp,
	     double S,  double Sz){

  // Calculates the matrix element
  // <Sp Szp|S+|S Sz> = sqrt((S-Sz)(S+Sz+1))
  // for S=Sp Szp=Sz+1

  if ( (dLT(Sp,0.0))||
       (dLT(S,0.0))||
       (dGT(fabs(Szp),Sp))||
       (dGT(fabs(Sz),S)) ) return(0.0);


  if ( (dEqual(S,Sp))&&(dEqual(Szp,Sz+1.0)) )
    return( sqrt((S-Sz)*(S+Sz+1.0)) );
  else
    return(0.0);


}


//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

double Sminus(double Sp, double Szp,
	      double S,  double Sz){

  // Calculates the matrix element
  // <Sp Szp|S-|S Sz> = sqrt((S-Sz)(S+Sz+1))
  // for S=Sp Szp=Sz+1

  if ( (dLT(Sp,0.0))||
       (dLT(S,0.0))||
       (dGT(fabs(Szp),Sp))||
       (dGT(fabs(Sz),S)) ) return(0.0);


  if ( (dEqual(S,Sp))&&(dEqual(Szp,Sz-1.0)) )
    return( sqrt((S+Sz)*(S-Sz+1.0)) );
  else
    return(0.0);


}
