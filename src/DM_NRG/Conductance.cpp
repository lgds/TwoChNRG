
#include <iostream>
#include <cstring>
#include <unistd.h> // contains getopt

#include <cstdlib>  // contains atoi

//#include "NRGclasses.hpp"
//#include "NRGfunctions.hpp"
//#include "SpecFuncClass.hpp"
#include "ConductanceClass.hpp"
using namespace std;


void CommandLineRead(int argc, char* argv[], CConductance* pG){

  int c;
  opterr = 0;

  char cmtemp[8],cbetabar[8];

  while ((c = getopt (argc, argv, "h R G ")) != -1)
    switch (c){
      case 'G':
	//	strcpy(cmtemp,optarg);
	//Mtemp = atoi(cmtemp);
	cout << " Graphene calculation " << endl;
	pG->IsGraphene=true;
	break;
      case 'R':
	//	strcpy(cmtemp,optarg);
	//Mtemp = atoi(cmtemp);
	cout << "Resistivity calculation " << endl;
	pG->CalcResistivity=true;
	break;
      case 'h':
	cout << "Usage: " << argv[0] << " (-R)" << endl;
	cout << "  -R : resistivity calculation " << endl;
	cout << "  -G : uses graphene model " << endl;
	exit(0);
	return;
       case '?':
	cout << "Usage: " << argv[0] << " (-R)" << endl;
	cout << "  -R : resistivity calculation " << endl;
	cout << "  -G : uses graphene model " << endl;
	exit(0);
 	return;
      default:
	return;
      }


}
////////////////////


////////////////////
//    Main code   //  
////////////////////




int main (int argc, char* argv[]){


  //CNRGCodeHandler ThisCode;

 
  CConductance G1;

  CommandLineRead(argc, argv, &G1);

  // Read Params, Temp
  G1.Initialize();
  G1.UseSpec=0; // 0 - Full calculation n>0: will use SpecVec[n-1]

  double omega1=0.0;

  // test
  cout << "G11(omega=" << omega1 <<") = " 
       << GSL_REAL(G1.SpecVec[0].GreensFunction(omega1)) << " + i "
       << GSL_IMAG(G1.SpecVec[0].GreensFunction(omega1))
       << endl; 	

  double Norm=G1.SpecVec[0].CalcNormInteg();

  cout << " Norm_spec = " << Norm << endl;  

  // Integral: Calculate T matrix (model) dependent and integrate

  double Cond=G1.CalcIntegral();


  G1.PrintParams();
  cout << " GoverG0= " << Cond ;

  if (G1.SpinResolved){
    double Cond2=G1.CalcIntegral(1);
    cout << " GoverG0_dn= " << Cond2 
	 << " (G_up+G_dn)/2= " << 0.5*(Cond+Cond2);
  } 
  // end if spin-resolved

  cout << endl;


//   double Delta1=G1.CondParams[1]; // Check this.

//   cout << "pi*Delta1*rho(0) = " <<  M_PI*Delta1*G1.SpecVec[0].RhoInterpol(0.0) << endl;


}
// end of main
