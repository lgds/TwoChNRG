
#include <iostream>
#include <cstring>
#include <unistd.h> // contains getopt

#include <cstdlib>  // contains atoi
using namespace std;


void DM_NRG_CommandLineRead(int argc, char* argv[], int &Mtemp, 
			    double &betabar, double &twindow,
			    double &broadtemp,
			    double &bbroad,
			    int &UseCFS,
			    bool &UseGap,
			    int &Nw){

  int c;
  opterr = 0;

  char cmtemp[8],cbetabar[8];

  char ctwindow[8],cbroadtemp[8],cbbroad[8];

  //All initialized in DMNRG.cpp
  //Mtemp=1000;
  //betabar=0.727;
  //UseCFS=0;
  //broadtemp=0.8;
  //NwEachShell=1

  while ((c = getopt (argc, argv, "h M: b: w: B: e: n: C G")) != -1)
    switch (c){
      case 'M':
	strcpy(cmtemp,optarg);
	Mtemp = atoi(cmtemp);
	break;
      case 'b':
	strcpy(cbetabar,optarg);
	betabar=atof(cbetabar);
	break;
      case 'w':
	strcpy(ctwindow,optarg);
	twindow=atof(ctwindow);
	break;
      case 'B':
	strcpy(cbroadtemp,optarg);
	broadtemp=atof(cbroadtemp);
	break;
      case 'e':
	strcpy(cbbroad,optarg);
	bbroad=atof(cbbroad);
	break;
      case 'n':
	strcpy(cmtemp,optarg);
	Nw = atoi(cmtemp);
	break;
      case 'G':
	UseGap=true;
	cout << " Gapped host. Adding peaks below the gap " << endl;
	break;
      case 'C':
	UseCFS=1;
	cout << " Using CFS to calculate spectral functions " << endl;
	break;
      case 'h':
	cout << "Usage: " << argv[0] << " (-e bbroad -M Mtemp -b betabar -w twindow -B broadtemp)" << endl;
	cout << " Default values: bbroad=0.5*log(Lambda); betabar=0.727 ; twindow=2 ; broadtemp=0.8 ; " << endl;
// 	cout << " Additional option: -C (uses CFS method) " << endl;
	cout << " Additional options: " << endl;
	cout << " -C : uses CFS method " << endl;
	cout << " -G : gapped host (uses gap in input_nrg.dat) " << endl;
	cout << " -n NomegasEachN : Calculates a denser mesh in omega (default=1) " << endl; 
	exit(0);
	return;
       case '?':
	cout << "Usage: " << argv[0] << " (-M Mtemp -b betabar -w twindow -B broadtemp)" << endl;
	cout << " Additional options: " << endl;
	cout << " -C : uses CFS method " << endl;
	cout << " -G : gapped host (uses gap in input_nrg.dat) " << endl;
	cout << " -n NomegasEachN : Calculates a denser mesh in omega (default=1) " << endl; 
	exit(0);
 	return;
      default:
	return;
      }


//   cout << " cMtemp = " << cmtemp 
//        << " cbetabar = " << cbetabar 
//        << endl; 


  cout << " Inside Mtemp = " << Mtemp 
       << " betabar = " << betabar 
       << endl; 




}
////////////////////
