
#include <cmath>
using namespace std;

double CalcDN(double Lambda, int Nshell, double z_twist){

  double HalfLambdaFactor=1.0;

  if (Lambda>0.0)
    HalfLambdaFactor=0.5*(1.0+(1.0/Lambda));

//   return(HalfLambdaFactor*pow(Lambda,(-(Nshell-1)/2.0) ));
  return(HalfLambdaFactor*pow(Lambda,(-(Nshell-1)/2.0) - z_twist + 1.0));


}
