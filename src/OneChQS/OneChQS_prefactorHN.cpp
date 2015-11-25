
#include <math.h>
using namespace std;

double OneChQS_prefactorHN(int type_i, int type_j, double Si)
{

  double aux=0.0;
  double TwoS=2.0*Si;

  if ( ((type_i==1)&&(type_j==2))||
       ((type_i==2)&&(type_j==1))||
       ((type_i==1)&&(type_j==3))||
       ((type_i==3)&&(type_j==1)) )
    {
      aux=1.0;
    }

  if ( ((type_i==2)&&(type_j==4))||
       ((type_i==4)&&(type_j==2)) )
    {
      aux=sqrt(TwoS/(TwoS+1.0));
    }

  if ( ((type_i==3)&&(type_j==4))||
       ((type_i==4)&&(type_j==3)) )
    {
      aux=-sqrt((TwoS+2.0)/(TwoS+1.0));
    }

  return(aux);

}
