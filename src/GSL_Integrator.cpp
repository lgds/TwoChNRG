
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>
using namespace std;


double GSL_Integrator(gsl_function *pIntegrandGSL,
		      gsl_integration_workspace *w_gsl, 
		      double x0, double x1){

  double Integ=0.0;
  
  /* GSL integration variables */
  double error,epsabs=1.0E-6,epsrel=1.0E-6;
  /* Error handling */
  gsl_error_handler_t *old_handler;
  int status=0;

  old_handler=gsl_set_error_handler_off();
  // Integral int_{-1}^{1} Delta(e))
  status=gsl_integration_qag(pIntegrandGSL,x0,x1,epsabs,epsrel,10000,3,
			     w_gsl,&Integ,&error);
  if (status){
    //      while (status!=0)
    while (status==GSL_EROUND){
      cout << " GSL_Integrator: Reducing error tolerance. " << endl;
      epsabs*=10.0;
      epsrel*=10.0;
      status=gsl_integration_qag(pIntegrandGSL,x0,x1,epsabs,epsrel,10000,3,
				 w_gsl,&Integ,&error);
    }
  }
  // end if status

  return(Integ);



}

//////////////
