
#include <math.h>
//
// Compare two doubles
//
bool dEqual(double A, double B){
  if (fabs(A-B)<1.e-20) return(true);
  else return(false);
}

bool dNEqual(double A, double B){
  if (fabs(A-B)>1.e-20) return(true);
  else return(false);
} 


bool dLEqual(double A, double B){
  if ((A<B)||(fabs(A-B)<1.e-20)) return(true);
  else return(false);
}

bool dGEqual(double A, double B){
  if ((A>B)||(fabs(A-B)<1.e-20)) return(true);
  else return(false);
}

//
// Compare two doubles setting the precision
//

bool dEqualPrec(double A, double B, double prec){
  if (fabs(A-B)<prec) return(true);
  else return(false);
}

bool dNEqualPrec(double A, double B, double prec){
  if (fabs(A-B)>prec) return(true);
  else return(false);
} 


bool dLEqualPrec(double A, double B, double prec){
  if ((A<B)||(fabs(A-B)<prec)) return(true);
  else return(false);
}

bool dGEqualPrec(double A, double B, double prec){
  if ((A>B)||(fabs(A-B)<prec)) return(true);
  else return(false);
}




//
// See if A>B
//

bool dGT(double A, double B){
  if (A>B) return(true);
  else return(false);
}

bool dLT(double A, double B){
  if (A<B) return(true);
  else return(false);
}

bool dGTPrec(double A, double B, double prec){
  if ((A>B)&&(fabs(A-B)>prec)) return(true);
  else return(false);
}

bool dLTPrec(double A, double B, double prec){
  if ((A<B)&&(fabs(A-B)>prec)) return(true);
  else return(false);

}


