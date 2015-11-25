

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)


#include <iostream>
#include <cstring>

//#include "SpecFuncClass.hpp"
#include "ConductanceClass.hpp"

using namespace std;

///////////////

void CConductance::Initialize(){

  CNRGCodeHandler ThisCode;
  CSpecFunction EachSpec;
  ifstream InFile;

  // Read code parameters (including ModelNo and BandNo)
  //ThisCode.ReadGenPars();
  ThisCode.ReadGenPars(true);
  //ThisCode.InitialSetUp();
  ThisCode.InitialSetUp(true);
  ThisCode.PrintSettings();

  Nops=ThisCode.NopsSaved;
  ModelNo=ThisCode.ModelNo;
  SymNo=ThisCode.SymNo;
  BandNo=ThisCode.BandNo;

  CondParams=ThisCode.dInitParams;

  cout << " Old Parameters : " << endl;
  PrintParams();
  cout << endl;


  // Adding chain params to conductance pars, if needed
  if (BandNo!=0){
    cout << "Non-flat band: BandNo= " << BandNo << endl;
//     cout << " Original Cond Parameters : " << endl;
//     PrintParams();
//     cout << endl;
    // I looove classes!
    ThisCode.chain.HybFuncParams=&CondParams;
    ThisCode.chain.ReadParams((char *)"lanc.in",6);
    // Ok, need to re-arrange the lanc.in parameters
    ReorderParams();
    //(char *) needed to avoid the warning message here
    cout << " New Cond Parameters : " << endl;
    PrintParams();
    cout << endl;
//     exit(0);
  }
  // end if special band.

  // Read spectral functions (method)
  int iop,jop;
  char Ciop[4],Cjop[4];

  int ispec=0;
  for (int iop=0;iop<Nops;iop++){
    for (int jop=0;jop<Nops;jop++){
      // Name
      sprintf(Ciop,"%d",iop);
      sprintf(Cjop,"%d",jop);
      strcpy(EachSpec.Name,"rho_");
      strcat(EachSpec.Name,Ciop);
      strcat(EachSpec.Name,"_");
      strcat(EachSpec.Name,Cjop);
      EachSpec.ClearOmegaRhow();
      EachSpec.ReadOmegaRhow();
      //EachSpec.PrintOmegaRhow();
      // Test if Omega_Rhow is there!
      if (EachSpec.Omega_Rhow[0].size()!=0){
	cout << " Accepting " << EachSpec.Name << "  Size = " 
	     << EachSpec.Omega_Rhow[0].size() 
	     << endl;
	SpecVec.push_back(EachSpec);
	ispec++;
      }
    }
  }
  // end loop in iop,jop

  cout << "No of Spectral functions read: " << SpecVec.size() << endl;

  // Read Temp
  //strcpy(arqname,"Temp_Mtemp.dat");
  InFile.open("Temp_Mtemp.dat");
  if (InFile.is_open()){
    InFile >> Temp >> Mtemp;
    cout << "(Initialize): Temp_Mtemp.dat found: Temp = " << Temp 
	 << " Mtemp = " << Mtemp << endl; 
  }
  else{
    cout << "(Initialize) Error: Cant open Temp_Mtemp.dat : using Temp=0.0 " << endl;
    Temp=0.0;
  }
  InFile.close();

  // Spin resolved?

  if ( (SymNo==2)||(SymNo==5) ) SpinResolved=true;


  // Set T(w) and Integrand

  SetIntegrand();

}
// end initialize

///////////////
void CConductance::ReorderParams(){

  // input_nrg.dat: DQD order is 
  // U1 Delta1 e1 U2 Delta2 e2 lambda hz_dot1 hz_dot2
  // 0.5 0 -0.25 0.0 0.02 -0.05 0.07 0 0

  double daux=0.0;

  if ( (BandNo>=4)&&(BandNo<=9) ){
  // with lanc.in, wrong order is 
  // U1 Delta1 e1 whichbandtype e2 Delta2 Delta1 lambda MagField
  // 0.5 0 -0.25 4 -0.05 0.02 0 0.07 0
    if (CondParams.size()>8){
      CondParams[1]=CondParams[6]; // Sets Delta1 ok
      CondParams[3]=0.0; // sets U2=0
      daux=CondParams[4];CondParams[4]=CondParams[5];
      CondParams[5]=daux;  // Corrects e2 and Delta2
      CondParams[6]=CondParams[7]; //Sets lambda
      CondParams[7]=CondParams[8]; //Sets hz1
      CondParams[8]=0.0; // hz_dot2=0.0
    }
  } 
  // end if 4<=BandNo<=9
  else if (BandNo==12){
    // Power-Law band (=Graphene)
  // with lanc.in, order is 
  // U1 Delta1 e1 whichbandtype  r  Gamma0 small_gamma w0
  // 0.5 0   -0.25    1         1.0 0.016      0.0    0.0
    CondParams[1]=CondParams[5]; // Sets Gamma as Gamma0
  }
  // end if BandNo==1 
}


///////////////
void CConductance::PrintParams(){

  cout << "Params: ";
  for (int ip=0; ip<CondParams.size(); ip++){
    cout << CondParams[ip] << " ";
  }
  cout << " Temp: " << Temp << " Mtemp " << Mtemp;

}
// end PrintParams

///////////////
void CConductance::SetIntegrand(){

  switch (ModelNo){
  case (0): // Anderson model
    switch (BandNo){
    case(0): // Flat band
      Integrand=Integrand_Anderson;
      if (CalcResistivity) Integrand=Integrand_Anderson_Resistivity;
      break;
    case(12): // Power-Law (graphene)
      Integrand=Integrand_Anderson; // Need to set this correctly!!
      if ( (IsGraphene)&&(CalcResistivity) ) Integrand=Integrand_Anderson_Resistivity_Graphene;
      else cout << "Warning: conductance is not set to Graphene..." << endl;
      break;
    case(4): // Side dot
      Integrand=Integrand_SideDot;
      break;
    default:
      Integrand=Integrand_Anderson;
    } // end switch BandNo
    break;
  case (6): // DQD
    //TMatrix=TMatrix_DQD; // HOW TO CALCULATE THIS???
    if (SpinResolved){
      cout << " Spin-resolved conductance " << endl;
      Integrand=Integrand_DQD_up;
      Integrand2=Integrand_DQD_dn;
    }
    else
      Integrand=Integrand_DQD;
    break;
  default :
    cout << " Conductance for ModelNo = " << ModelNo << " not implemented. Exiting..." << endl;
    exit(0);
  }

}


///////////////
double CConductance::CalcIntegral(int WhichIntegrand){

  double Cond=0.0;

  // Ok, need to think about this...
  //Integrand=Integrand_DQD;

  gsl_function IntegrandGSL;
 
  switch (WhichIntegrand){
  case 0:
    IntegrandGSL.function=Integrand;
    break;
  case 1:
    IntegrandGSL.function=Integrand2;
    break;
  default:
    IntegrandGSL.function=Integrand;
  }    
  //IntegrandGSL.function=Integrand;
  // Will THIS work??? Yes!!! Holy cow!
  CConductance Gaux; 
  Gaux.SpecVec=SpecVec;
  Gaux.CondParams=CondParams;
  Gaux.Temp=Temp;
  Gaux.Mtemp=Mtemp;
  Gaux.UseSpec=UseSpec;

  IntegrandGSL.params=&Gaux;


  // Debug
//   for (double omega=-1.0;omega<=1.0; omega+=0.05){
//     cout << " omega: " << omega 
// 	 << " Integrand: " << GSL_FN_EVAL(&IntegrandGSL,omega) 
// 	 << endl;
//   }
  //

  if (dEqual(Temp,0.0)){

    Cond=GSL_FN_EVAL(&IntegrandGSL,0.0);
    
  } else{
    //double w0=-1.0;
    //double w1=1.0;
    double w0=-10*Temp;
    double w1=10*Temp;

    //Allocate integration workspace
    gsl_integration_workspace *w_gsl;
    w_gsl = gsl_integration_workspace_alloc (10000);

    Cond=GSL_Integrator(&IntegrandGSL,w_gsl,w0, w1);

  // Debug
//     cout << "CondA = " << Cond << endl;
//     double CondPlus=0.0;
//     double CondMinus=0.0;
//     w0=-10.0*Temp;
//     w1=10*Temp;
//     double L1=2.5;
//     double w1p=w1;
//     double w0p=w1/L1;
//     double w0m=w0;  //assuming w0<0
//     double w1m=w0/L1;
//     double auxIntegrand=0.0;
//     for (int N_int=0; N_int<31; N_int++){

//       CondPlus+=GSL_Integrator(&IntegrandGSL,w_gsl,w0p, w1p);
//       CondMinus+=GSL_Integrator(&IntegrandGSL,w_gsl,w0m, w1m);

//       auxIntegrand=GSL_FN_EVAL(&IntegrandGSL,0.5*(w0p+w1p));
//       cout << " Integrand: w= " << 0.5*(w0p+w1p) 
//       << " tau= " <<  auxIntegrand << endl;
//       auxIntegrand=GSL_FN_EVAL(&IntegrandGSL,0.5*(w0m+w1m));
//       cout << " Integrand: w= " << 0.5*(w0m+w1m) 
//       << " tau= " <<  auxIntegrand << endl;

// //       cout << " i=" << N_int << " w0p= " << w0p << " w1p= " << w1p << endl;
// //       cout << " i=" << N_int << " w0m= " << w0m << " w1m= " << w1m << endl;
  
//       w1p=w0p;
//       w0p/=L1;

//       w0m=w1m;
//       w1m/=L1;
//     }
//     //
//     cout << "CondB = " << CondPlus+CondMinus << endl;


    gsl_integration_workspace_free(w_gsl);

  }


  if (CalcResistivity) Cond=1.0/Cond;

  return(Cond);

}

///////////////

double CConductance::FermiFunction(double En, double Temp)
{
  // Calculates e^(-En/Temp) //
  // If Temp==0.0 : En==0.0, returns 1 //
  //              : En!=0.0, returns 0 //

  double aux=0.0;

  if (dEqual(Temp,0.0)){
    if (En>=0.0)
      aux=0.0;
    else aux=1.0;    
  }
  else
    aux=1.0/(exp(En/Temp)+1.0);

  return (aux);
}

///////////////
////////////////////////////////////////////////
///                                         ////
///      END OF CConductance functions      ////
///                                         ////
////////////////////////////////////////////////


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
///                                         ////
///                                         ////
///         Tmatrix definitions             ////
///                                         ////
///   Static functions (model dependent)    ////
///                                         ////
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

gsl_complex TMatrix_Anderson(double omega,
			void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);

  CConductance Gaux=*(CConductance *)params;


  gsl_complex caux;
  //GSL_SET_COMPLEX(&caux,0.0,3.0);

  // Delta1*G11(omega)

  double Delta1=Gaux.CondParams[1]; // Check this.

  // Debugging
//   cout << "Delta 1 = " << Delta1 
//        << " rho(w=" << omega << ") = " << Gaux.SpecVec[0].RhoInterpol(omega) 
//        << endl;

  //cout << "pi*Delta1*rho(0) = " <<  M_PI*Delta1*Gaux.SpecVec[0].RhoInterpol(0.0) << endl;
  caux=Gaux.SpecVec[0].GreensFunction(omega);
  caux=gsl_complex_mul_real (caux,Delta1);

  return(caux);

}

//////


double Integrand_Anderson(double omega,
		     void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);
  //double A=*(double *)params;

  // Would THIS work??? Yes!!! Holy cow!

  CConductance Gaux=*(CConductance *)params;


  gsl_complex Tmat;

  // Why this??
  //Tmat=TMatrix_DQD(omega,params);
  Tmat=TMatrix_Anderson(omega,params);

  // testing!

  double aux=-GSL_IMAG(Tmat); // /M_PI

  // Put temperature in here, my friend!

  double Temp=Gaux.Temp;
  double FD=Gaux.FermiFunction(omega,Temp);

  // add -df/dw =(1.0/Temp)*FD*(1.0-FD)


  if (dNEqual(Temp,0.0)){
    aux*=(1/Temp)*FD*(1.0-FD);    
  } 

  return(aux);

}

////////////////

double Integrand_Anderson_Resistivity(double omega,
		     void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);
  //double A=*(double *)params;

 
  CConductance Gaux=*(CConductance *)params;

  gsl_complex Tmat;

  Tmat=TMatrix_Anderson(omega,params);

  // Tmat=Gamma*G_d(omega) AND
  // tau(omega)=(Gamma*rho(omega))^-1=[-Gamma*Im(G_d(omega))/pi]^-1
  // Therefore: tau(omega)=-pi/Gamma*Im(G_d(omega)=-pi/Im(Tmat)


  // tau

  double aux=-M_PI/GSL_IMAG(Tmat); // Is the pi factor ok???

  // Put temperature in here, my friend!

  double Temp=Gaux.Temp;
  double FD=Gaux.FermiFunction(omega,Temp);

  // add -df/dw =(1.0/Temp)*FD*(1.0-FD)


  if (dNEqual(Temp,0.0)){
    aux*=(1/Temp)*FD*(1.0-FD);    
  } 

  return(aux);

}

////////////////


double Integrand_Anderson_Resistivity_Graphene(double omega,
		     void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);
  //double A=*(double *)params;

 
  CConductance Gaux=*(CConductance *)params;

  gsl_complex Tmat;

  Tmat=TMatrix_Anderson(omega,params);

  // Tmat=Gamma*G_d(omega) AND
  // tau(omega)=(Gamma*rho(omega))^-1=[-Gamma*Im(G_d(omega))/pi]^-1
  // Therefore: tau(omega)=-pi/Gamma*Im(G_d(omega)=-pi/Im(Tmat)

  // tau


  double r=Gaux.CondParams[4];
  double Gamma0=Gaux.CondParams[5];
  double w0=Gaux.CondParams[6];

  double aux=-M_PI/GSL_IMAG(Tmat); // Is the pi factor ok???

  // Graphene: Integrand is tau(omega)*|omega-w0| (not simply tau)
  // In our language mu=-w0: Gamma(w)=Gamma0*|w-w0|


  aux*=Gamma0*fabs(omega-w0);

  // Put temperature in here, my friend!

  double Temp=Gaux.Temp;
  double FD=Gaux.FermiFunction(omega,Temp);

  // add -df/dw =(1.0/Temp)*FD*(1.0-FD)


  if (dNEqual(Temp,0.0)){
    aux*=(1/Temp)*FD*(1.0-FD);    
  } 

  return(aux);

}

////////////////


gsl_complex TMatrix_DQD(double omega,
			void *params,int UpDn){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);

  // UpDn=0 : Not spin resolved (QS symmetry)
  // UpDn=1 : Spin up
  // UpDn=2 : Spin dn

  int iSpec[4]={0,1,2,3};

  switch (UpDn){
  case 1: // spin up: iSpec={0,1,4,5};
    iSpec[0]=0;
    iSpec[1]=1;
    iSpec[2]=4;
    iSpec[3]=5;
    break;
  case 2: // spin dn: iSpec={2,3,6,7};
    iSpec[0]=2;
    iSpec[1]=3;
    iSpec[2]=6;
    iSpec[3]=7;
    break;
   default:
    iSpec[0]=0;
    iSpec[1]=1;
    iSpec[2]=2;
    iSpec[3]=3;
  }
  // end Switch UpDn
 
  CConductance Gaux=*(CConductance *)params;

  gsl_complex caux[3];
  //GSL_SET_COMPLEX(&caux,0.0,3.0);

  // Delta1*G11(omega)

  double Delta1=Gaux.CondParams[1]; // Check this: input_params order
  double Delta2=Gaux.CondParams[4]; 

  double Delta12=sqrt(Delta1*Delta2);


  // Debugging
//   cout << " Delta 1 = " << Delta1 
//        << " Delta 2 = " << Delta2 << endl
//        << " rho11(w=" << omega << ") = " << Gaux.SpecVec[iSpec[0]].RhoInterpol(omega) 
//        << endl
//        << " rho12(w=" << omega << ") = " << Gaux.SpecVec[iSpec[1]].RhoInterpol(omega) 
//        << endl
//        << " rho21(w=" << omega << ") = " << Gaux.SpecVec[iSpec[2]].RhoInterpol(omega) 
//        << endl
//        << " rho22(w=" << omega << ") = " << Gaux.SpecVec[iSpec[3]].RhoInterpol(omega) 
//        << endl;

//   cout << "pi*Delta1*rho(0) = " << M_PI*Delta1*Gaux.SpecVec[iSpec[0]].RhoInterpol(0.0) << endl;

  ////////////

  // lambda term is missing!!

  // Delta1 G_11
  caux[0]=Gaux.SpecVec[iSpec[0]].GreensFunction(omega);
  caux[0]=gsl_complex_mul_real (caux[0],Delta1);
  //Delta_12 ( G_12 + G_21 )
  caux[1]=Gaux.SpecVec[iSpec[1]].GreensFunction(omega);
  caux[1]=gsl_complex_add(caux[1],Gaux.SpecVec[iSpec[2]].GreensFunction(omega));
  caux[1]=gsl_complex_mul_real (caux[1],Delta12);
  //Delta2*G22
  caux[2]=Gaux.SpecVec[iSpec[3]].GreensFunction(omega);
  caux[2]=gsl_complex_mul_real (caux[2],Delta2);




  // T_DQD=Delta1 G_11 + Delta_12 ( G_12 + G_21 ) + Delta2*G22
  caux[0]=gsl_complex_add(caux[0],caux[1]);
  caux[0]=gsl_complex_add(caux[0],caux[2]);


//   cout << " Delta2*G22(0.0) = " 
//        << GSL_REAL(caux[2]) << " + i "
//        << GSL_IMAG(caux[2])
//        << endl; 	
//   cout << " T(0.0) = " 
//        << GSL_REAL(caux[0]) << " + i "
//        << GSL_IMAG(caux[0])
//        << endl; 	


  return(caux[0]);

}

//////

double Integrand_DQD(double omega,
		     void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);
  //double A=*(double *)params;

  // Would THIS work??? Yes!!! Holy cow!

  CConductance Gaux=*(CConductance *)params;
  gsl_complex Tmat;

  Tmat=TMatrix_DQD(omega,params);

  double aux=-GSL_IMAG(Tmat);
  double Temp=Gaux.Temp;
  double FD=Gaux.FermiFunction(omega,Temp);

  // add -df/dw =(1.0/Temp)*FD*(1.0-FD)

  if (dNEqual(Temp,0.0))
    aux*=(1/Temp)*FD*(1.0-FD);
  
  return(aux);

}

////////////////


double Integrand_DQD_up(double omega,
			void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);
  //double A=*(double *)params;

  // Would THIS work??? Yes!!! Holy cow!

  CConductance Gaux=*(CConductance *)params;
  gsl_complex Tmat;

  Tmat=TMatrix_DQD(omega,params,1);

  double aux=-GSL_IMAG(Tmat);
  double Temp=Gaux.Temp;
  double FD=Gaux.FermiFunction(omega,Temp);

  // add -df/dw =(1.0/Temp)*FD*(1.0-FD)

  if (dNEqual(Temp,0.0))
    aux*=(1/Temp)*FD*(1.0-FD);
  
  return(aux);

}
/// dn

double Integrand_DQD_dn(double omega,
			void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);
  //double A=*(double *)params;

  // Would THIS work??? Yes!!! Holy cow!

  CConductance Gaux=*(CConductance *)params;
  gsl_complex Tmat;

  Tmat=TMatrix_DQD(omega,params,2);

  double aux=-GSL_IMAG(Tmat);
  double Temp=Gaux.Temp;
  double FD=Gaux.FermiFunction(omega,Temp);

  // add -df/dw =(1.0/Temp)*FD*(1.0-FD)

  if (dNEqual(Temp,0.0))
    aux*=(1/Temp)*FD*(1.0-FD);
  
  return(aux);

}


////////////////





////////////////

gsl_complex TMatrix_SideDot(double omega,
			void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);

  CConductance Gaux=*(CConductance *)params;

  gsl_complex caux[3];
  gsl_complex Tmat;


  gsl_complex zaux[5],AI;
  gsl_complex G022,G11;

  //gsl_complex termo[3];


  // Delta1*G11(omega)
  
  // input_nrg params order 
  double Delta1=Gaux.CondParams[1]; 
  double Delta2=Gaux.CondParams[4]; 
  double e2=Gaux.CondParams[5];
  double lambda=Gaux.CondParams[6];

  double Delta12=sqrt(Delta1*Delta2);


//   cout << " SideDot: e2 = " << e2
//        << " Delta2 = " << Delta2 
//        << " Delta1 = " << Delta1
//        << " lambda = " << lambda
//        << endl;

  // Debugging

  // T(w) =  Delta_1 G_11 + 2 Delta_12 [G^0_22 (lambda - iDelta_12)G_11] + 
  //            
  //          Delta_2 [G^0_22 + (G^0_22)^2 (lambda - iDelta_12)^2 G_11] =
  //

  // CHECK 0.5 factors!!

  GSL_SET_COMPLEX(&caux[0],0.0,0.0);
  GSL_SET_COMPLEX(&caux[1],0.0,0.0);
  GSL_SET_COMPLEX(&caux[2],0.0,0.0);
  GSL_SET_COMPLEX(&AI,0.0,1.0);
  GSL_SET_COMPLEX(&Tmat,0.0,0.0);
  GSL_SET_COMPLEX(&G11,0.0,0.0);

  // Calculate G0_22
  if ( (fabs(Delta2)>1.0e-10) ){   
      /*Calculate G022=(omega-e2+iDelta2)^-1*/
      zaux[1]=gsl_complex_mul_real(AI,Delta2);
      zaux[1]=gsl_complex_add_real(zaux[1],omega-e2); 
      /* omega-en0 corrected in Nov 16/06 */
      G022=gsl_complex_inverse(zaux[1]);
  } 
  else 
    GSL_SET_COMPLEX(&G022,0.0,0.0); /*Disregard G022 */

  // Calculate term WITHOUT G11!

  //Second term: Delta12*[G22*(lambda- i Delta12)*G11] | 
  // CHECK: is there a factor of 2 here? Yes

  // -i Delta12
  zaux[0]=gsl_complex_mul_real(AI,-Delta12);

  // (lambda- i Delta12) 
  zaux[0]=gsl_complex_add_real(zaux[0],lambda); 

  // (lambda- i Delta12)^2 
  zaux[1]=gsl_complex_mul(zaux[0],zaux[0]); 
  
  // 2*Delta12*G022*(lambda- i Delta12)
  zaux[2]=gsl_complex_mul(G022,zaux[0]);
  zaux[2]=gsl_complex_mul_real(zaux[2],2.0*Delta12);

  //Third term: termo[2]=Delta2*[G22 + G22^2*(lambda - i Delta12)^2*G11] */

  //  G022^2*(lambda - i Delta12)^2
  zaux[3]=gsl_complex_mul(G022,G022);
  zaux[3]=gsl_complex_mul(zaux[3],zaux[1]);

  // Delta1 G_11 (spin up)
  if (Gaux.UseSpec<1){
    G11=Gaux.SpecVec[0].GreensFunction(omega);
    caux[0]=gsl_complex_mul_real (G11,Delta1);
    //Second term: 2*Delta12*[G022*(lambda-i Delta12)*G11] | 
    caux[1]=gsl_complex_mul(G11,zaux[2]); /* G022 (lambda- i Delta12) G11 */
    // Third term: termo[2]=Delta2*[G22 + G22^2*(lambda- i Delta12)^2*G11] */
    caux[2]=gsl_complex_mul(zaux[3],G11);
    caux[2]=gsl_complex_add(G022,caux[2]);
    caux[2]=gsl_complex_mul_real(caux[2],Delta2);

    Tmat=gsl_complex_add(caux[0],caux[1]);
    Tmat=gsl_complex_add(Tmat,caux[2]); //spin up

    if (Gaux.SpecVec.size()>1){
      Tmat=gsl_complex_mul_real(Tmat,0.5);

      // Use last (should be diagonal)
      G11=Gaux.SpecVec[Gaux.SpecVec.size()-1].GreensFunction(omega);
      caux[0]=gsl_complex_mul_real (G11,Delta1);
      //Second term: 2*Delta12*[G022*(lambda-i Delta12)*G11] | 
      caux[1]=gsl_complex_mul(G11,zaux[2]); /* G022 (lambda- i Delta12) G11 */
      // Third term: termo[2]=Delta2*[G22 + G22^2*(lambda- i Delta12)^2*G11] */
      caux[2]=gsl_complex_mul(zaux[3],G11);
      caux[2]=gsl_complex_add(G022,caux[2]);
      caux[2]=gsl_complex_mul_real(caux[2],Delta2);

      caux[0]=gsl_complex_add(caux[0],caux[1]);
      caux[0]=gsl_complex_add(caux[0],caux[2]); //spin down

      caux[0]=gsl_complex_mul_real(caux[0],0.5);
     
      Tmat=gsl_complex_add(Tmat,caux[0]); // 1/2*(spin up + spin down)


    }// if we have spin down as well

  }
  else{

    G11=Gaux.SpecVec[Gaux.UseSpec-1].GreensFunction(omega);
    caux[0]=gsl_complex_mul_real (G11,Delta1);
    //Second term: 2*Delta12*[G022*(lambda-i Delta12)*G11] | 
    caux[1]=gsl_complex_mul(G11,zaux[2]); /* G022 (lambda- i Delta12) G11 */
    // Third term: termo[2]=Delta2*[G22 + G22^2*(lambda- i Delta12)^2*G11] */
    caux[2]=gsl_complex_mul(zaux[3],G11);
    caux[2]=gsl_complex_add(G022,caux[2]);
    caux[2]=gsl_complex_mul_real(caux[2],Delta2);

    Tmat=gsl_complex_add(caux[0],caux[1]);
    Tmat=gsl_complex_add(Tmat,caux[2]);


  }
  // end if Gaux.UseSpec==0 or <1


//   cout << " Tmat(0.0) = " 
//        << GSL_REAL(Tmat) << " + i "
//        << GSL_IMAG(Tmat)
//        << endl; 	



  return(Tmat);

}


///////


///////

double Integrand_SideDot(double omega,
			 void *params){
  // Get parameters from params. Example
  //vector<double> VecParam=*(vector<double> *)(params);
  //double A=*(double *)params;

  // Would THIS work??? Yes!!! Holy cow!

  CConductance Gaux=*(CConductance *)params;


  gsl_complex Tmat;

  Tmat=TMatrix_SideDot(omega,params);

  // testing!

  double aux=-GSL_IMAG(Tmat);
  double Temp=Gaux.Temp;
  double FD=Gaux.FermiFunction(omega,Temp);

  // add -df/dw =(1.0/Temp)*FD*(1.0-FD)


  if (dNEqual(Temp,0.0)){

    aux*=(1/Temp)*FD*(1.0-FD);
    
  } 

  return(aux);

}

////////////////

