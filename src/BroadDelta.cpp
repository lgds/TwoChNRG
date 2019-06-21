
#include <math.h>
#include "MyStructures.h"
#define pi 3.141592653589793238462643383279502884197169

/****************************************************************************/
/**********************     Broadened Delta Function      *******************/
/****************************************************************************/



double BroadDelta(double omega, double Ep, double b){


double aux[3];
 if (Ep==0.0) Ep=1.0E-100;
////////////////

// omega and Ep need to have the same sign other wise, returns zero

 if ( (omega*Ep)<0.0 ){return(0.0);}

 /* aux[0]=exp(-b*b/4.0)/(b*Ep*sqrt(pi)); */
 aux[0]=exp(-b*b/4.0)/(b*fabs(Ep)*sqrt(pi));
 aux[1]=log(fabs(omega/Ep))/b;
 aux[2]=aux[0]*exp(-aux[1]*aux[1]);

 return(aux[2]);
}
// END BroadDelta

/*********************************************************************/
/******************      BroadDelta2GSL          *********************/
/*********************************************************************/
/*********************************************************************/
//
/*  Same thing as before, but in a format to be used */
/*  with Gnu Scientific Library routines */
//
double BroadDelta2GSL(double omega, void *pars){

  double Ep; 
  double b;
  double aux;


  struct TwoDouble  *params1 =(struct TwoDouble *)pars;

  Ep = (params1->dvar1);
  b = (params1->dvar2);


  aux=BroadDelta(omega, Ep, b);
  return aux;
}
// END of BroadDelta2GSL

/****************************************************************************/
/**************     Lorentzian-broadened Delta Function      *****************/
/****************************************************************************/


double LorentzDelta(double omega, double Ep, double b){

double aux[3];

////////////////
    aux[0]=1.0/(2.0*pi);
    aux[1]=b/( (omega-Ep)*(omega-Ep)+b*b );
    aux[2]=aux[0]*aux[1];

    return(aux[2]);
}
// END LorentzDelta

double LorentzDeltaAnders(double omega, double Ep, double b){

double aux[3];

////////////////
    aux[0]=1.0/(pi);
    aux[1]=b/( (omega-Ep)*(omega-Ep)+b*b );
    aux[2]=aux[0]*aux[1];

    return(aux[2]);
}
// END LorentzDeltaAnders



/****************************************************************************/
/**************     Gaussian-broadened Delta Function      *****************/
/****************************************************************************/


double GaussDelta(double omega, double Ep, double b){

double aux[3];

////////////////
    aux[0]=1.0/(b*sqrt(pi));
    aux[1]=(omega-Ep)/b;
    aux[2]=aux[0]*exp(-aux[1]*aux[1]);

    return(aux[2]);
}
// END GaussDelta


/****************************************************************************/
/************     Log-gaussian  Delta Function (Weichselbaum) ***************/
/****************************************************************************/



double LogGaussDelta(double omega, double Ep, double b){


double aux[3];
 if (Ep==0.0) Ep=1.0E-100;
////////////////

// omega and Ep need to have the same sign other wise, returns zero

 if ( (omega*Ep)<0.0 ){return(0.0);}

 aux[0]=1.0/(b*fabs(Ep)*sqrt(pi));
 aux[1]=log(fabs(omega/Ep))/b;
 aux[2]=aux[0]*exp(-(aux[1]-b/4.0)*(aux[1]-b/4.0));
 
 return(aux[2]);
}
// END BroadDelta
