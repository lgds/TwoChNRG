

#include <cmath>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "SpecFuncClass.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)


#include "MyStructures.h"



using namespace std;



//////////////////
double CSpecFunction::CFS_SpecDens_M(double omegabar,int Nshell, bool CalcNorm){

  // Need to get the Connecting matrix elements and energies
  // to and from a given block (say, blGS)
  // Interesting procedure...

  // Loop in blocks of Operator


  double rho_wM=0.0;

  if ((Nshell>Mtemp)||(Nshell>NshellMax-1)){
    cout << "CFS_SpecDens_M: Nshell > M . Return 0" << endl;
    return(0.0);
  }

  if ((Nshell==Mtemp)||(Nshell==NshellMax-1)){
    if (dEqual(Temp,0.0))
      Betabar=1.0e10;
//     cout << "CFS_SpecDens_M: Last Nshell: Costi. CalcNorm= " << CalcNorm << endl;
//     cout << "CFS_SpecDens_M: T=0 calculation (Costi)." << endl;
//     cout << " Betabar = " << Betabar << endl;
    rho_wM=CalcSpecCosti_T_N(Nshell,Betabar,omegabar,CalcNorm);
    return(rho_wM);
  }
  //

  double Positive_w=0.0;
  double Negative_w=0.0;

  bool printstuff=false;

  // Debug
 //  if ( (Nshell==7)&&(fabs(fabs(omegabar)-0.0457)<0.001) ) printstuff=true; else printstuff=false;

  // Need to Check Sync!!
  bool ChkSyncRhoAcut=RhoN[Nshell].ChkSync((&AcutN[Nshell]));
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));

  if ( (!ChkSyncRhoAcut)||(!ChkSyncOp1Acut) ){
    cout << "CFS_SpecDens_M: RhoN or Op1N not in sync with AcutN " << endl;
    return(0.0);
  }

  if (printstuff)
    cout << " CFS_SpecDens_M: Nshell = " << Nshell
	 << " omegabar = " << omegabar
	 << " CalcNorm = " << CalcNorm
	 << endl;


  for (int iMatbl=0;iMatbl<Op1N[Nshell].NumMatBlocks();iMatbl++){
  // Checking a few blocks first
  //for (int iMatbl=0;iMatbl<3;iMatbl++){

    int ibl1=0;
    int ibl2=0;

    Op1N[Nshell].GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);


    // Ok, remember the input Op2N is the COMPLEX CONJUGATE of what we want!
    //  
    //  So (B)_{ibl2 ibl1} = ((B+)_{ibl1 ibl2})* <- THIS is Op2

    //if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl2,ibl1)<0) ){
    if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl1,ibl2)<0) ){
      cout << "CFS_SpecDens_M: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return(0.0);
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 
 
    // Check block sizes
    if ( (Nst_bl1!=RhoN[Nshell].GetBlockSize(ibl1))||
	 (Nst_bl2!=RhoN[Nshell].GetBlockSize(ibl2)) ){
      cout << "CFS_SpecDens_M: Block in Op1 not the same size as block in RhoN " << endl;
      return(0.0);
    }
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "CFS_SpecDens_M: Block in Op1 not the same size as block in Op2 " << endl;
      return(0.0);
    }

    int Nst_bl1_kept=Op1N[Nshell].GetBlockSize(ibl1,true); 
    int Nst_bl1_disc=Op1N[Nshell].GetBlockSize(ibl1,false);  

    int Nst_bl2_kept=Op1N[Nshell].GetBlockSize(ibl2,true); 
    int Nst_bl2_disc=Op1N[Nshell].GetBlockSize(ibl2,false); 

    if (printstuff){
      cout << " ibl1 = " << ibl1
	   << " ibl2 = " << ibl2
	   << endl;

      cout << " Nst_bl1= " << Nst_bl1 
	   << " Nst_bl1_kept= " << Nst_bl1_kept 
	   << " Nst_bl1_disc=" << Nst_bl1_disc	
	   << endl
	   << " Nst_bl2= " << Nst_bl2 
	   << " Nst_bl2_kept= " << Nst_bl2_kept 
	   << " Nst_bl2_disc=" << Nst_bl2_disc
	   << endl;

    }
    // end printstuff


    if ( (Nst_bl1!=Nst_bl1_kept+Nst_bl1_disc)||
	 (Nst_bl2!=Nst_bl2_kept+Nst_bl2_disc) ) {
      cout << "CFS_SpecDens_M: Error: Kept+Dics NOT equal to Nstblock! " << endl;
      return(0.0);
    }

    // Gotta check if ANY of these are zero!

    if (  ( (Nst_bl1_kept==0)&&(Nst_bl2_kept==0) )||
	  ( (Nst_bl1_disc==0)&&(Nst_bl2_disc==0) )   ){
      if (printstuff) cout << " Blocks " << ibl1 << " and " << ibl2 
	   << " are either both kept or both discarded. Moving on..." << endl;
    }else{

      // Blas matrices

      double Termo1=0.0;
//       if ( (Nst_bl1_kept!=0)&&(Nst_bl2_disc!=0)&&(omegabar>=0.0) ){
      if ( (Nst_bl1_kept!=0)&&(Nst_bl2_disc!=0)&&((omegabar>=0.0)||(CalcNorm)) ){
	// caculate term 1

	// REAL operators
	if ( ( !(Op1N[Nshell].IsComplex) )&&( !(RhoN[Nshell].IsComplex) ) ) {

	  boost::numeric::ublas::matrix<double> Rho1(Nst_bl1_kept,Nst_bl1_kept);

	  boost::numeric::ublas::matrix<double> opA1w(Nst_bl1_kept,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> opA1(Nst_bl1_kept,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> opB1(Nst_bl2_disc,Nst_bl1_kept);
	  // aux
	  boost::numeric::ublas::matrix<double> Aux2disc(Nst_bl2_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> Prod1(Nst_bl2_disc,Nst_bl1_kept);

	  Rho1=RhoN[Nshell].MatBlock2BLAS(ibl1,ibl1,true,true);

	  opA1=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2,true,false);
	  // Nst_bl1_kept x Nst_bl2_disc matrix (assuming real!)


	  if (NonDiagGF)
	    opB1=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2,true,false)); // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
	  else
	    opB1=trans(opA1);

	  // Positive and negative energies matrix elements
	  // Note for the future: sign is "-" if bosonic operators
	  //noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix

	  noalias(Prod1)=prod(opB1,Rho1); // Nb2_disc x Nb1_kept

	  if (CalcNorm) opA1w=opA1; else
	    noalias(opA1w)=MijxBDeltaEij(opA1,(&AcutN[Nshell]),ibl1,ibl2,true,false,omegabar); // Nb1_kept x Nb2_disc

	  noalias(Aux2disc)=prod(Prod1,opA1w);  // Nb2_disc x Nb2_disc

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	    diag1(Aux2disc,
		  boost::numeric::ublas::range(0,Aux2disc.size1()),
		  boost::numeric::ublas::range(0,Aux2disc.size2())); 


	  //  	if ( (ibl1==2)&&(ibl2==3) ){

	  // 	  Op1N[Nshell].PrintBlock(ibl1);
	  // 	  Op1N[Nshell].PrintBlock(ibl2);

	  // 	  //Op1N[Nshell].PrintMatBlock(ibl1,ibl2);
	  
	  // 	  cout << "Rho1 : " << endl;
	  // 	  for (unsigned ii = 0; ii < Rho1.size1 (); ++ii){
	  // 	    for (unsigned jj = 0; jj < Rho1.size2 (); ++jj)
	  // 	      if ( dEqualPrec(Rho1(ii,jj),0.0,1.e-15) ) cout << "0,";
	  // 	      else cout << Rho1(ii,jj) << ",";
	  // 	    cout << endl;
	  // 	  }

	  // 	  cout << "opA1 : " << endl;
	  // 	  for (unsigned ii = 0; ii < opA1.size1 (); ++ii){
	  // 	    for (unsigned jj = 0; jj < opA1.size2 (); ++jj)
	  // 	      if ( dEqualPrec(opA1(ii,jj),0.0,1.e-15) ) cout << "0,";
	  // 	      else cout << opA1(ii,jj) << ",";
	  // 	    cout << endl;
	  // 	  }

	  // 	  cout << "opA1w : " << endl;
	  // 	  for (unsigned ii = 0; ii < opA1w.size1 (); ++ii){
	  // 	    for (unsigned jj = 0; jj < opA1w.size2 (); ++jj)
	  // 	      if ( dEqualPrec(opA1(ii,jj),0.0,1.e-15) ) cout << "0,";
	  // 	      else cout << opA1w(ii,jj) << ",";
	  // 	    cout << endl;
	  // 	  }

	  // // 	  cout << "opB1 : " << endl;
	  // // 	  for (unsigned ii = 0; ii < opB1.size1 (); ii++){
	  // // 	    for (unsigned jj = 0; jj < opB1.size2 (); jj++)
	  // // 	      if ( dEqualPrec(opB1(ii,jj),0.0,1.e-15) ) cout << "0,";
	  // // 	      else cout << opB1(ii,jj) << ",";
	  // // 	    cout << endl;
	  // // 	  }


	  // 	  cout << "Prod1 : " << endl;
	  // 	  for (unsigned ii = 0; ii < Prod1.size1 (); ++ii){
	  // 	    for (unsigned jj = 0; jj < Prod1.size2 (); ++jj)
	  // 	      if ( dEqualPrec(Prod1(ii,jj),0.0,1.e-15) ) cout << "0,";
	  // 	      else cout << Prod1(ii,jj) << ",";
	  // 	    cout << endl;
	  // 	  }

	  // 	  cout << "diag : " << endl;
	  // 	  for (unsigned ii = 0; ii < diag1.size(); ++ii){
	  // 	    if ( dEqualPrec(diag1(ii),0.0,1.e-15) ) cout << "0,";
	  // 	    else cout << diag1(ii) << ",";
	  // 	  }

	  // 	}
	  // end debug

	  Termo1=sum(diag1);


	// COMPLEX operators
	} else {

	  if (printstuff)
	    cout << " ComplexOperator. Termo1 " 
		 << endl;
	  // if printstuff

	  boost::numeric::ublas::matrix<complex<double> > cRho1(Nst_bl1_kept,Nst_bl1_kept);

	  boost::numeric::ublas::matrix<complex<double> > copA1w(Nst_bl1_kept,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > copA1(Nst_bl1_kept,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > copB1(Nst_bl2_disc,Nst_bl1_kept);
	  // aux
	  boost::numeric::ublas::matrix<complex<double> > cAux2disc(Nst_bl2_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > cProd1(Nst_bl2_disc,Nst_bl1_kept);


	  cRho1=RhoN[Nshell].cMatBlock2BLAS(ibl1,ibl1,true,true);

	  copA1=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2,true,false);
	  // Nst_bl1_kept x Nst_bl2_disc matrix (assuming real!)

	  if (NonDiagGF)
	    //opB=Op2N[Nshell].MatBlock2BLAS(ibl2,ibl1);
	    copB1=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2,true,false)); // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
	  else
	    copB1=trans(copA1);

	  // Positive and negative energies matrix elements
	  // Note for the future: sign is "-" if bosonic operators

	  noalias(cProd1)=prod(copB1,cRho1); // Nb2_disc x Nb1_kept

	  if (CalcNorm) copA1w=copA1; else
	    noalias(copA1w)=cMijxBDeltaEij(copA1,(&AcutN[Nshell]),ibl1,ibl2,true,false,omegabar); // Nb1_kept x Nb2_disc

	  noalias(cAux2disc)=prod(cProd1,copA1w);  // Nb2_disc x Nb2_disc

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	    cdiag1(cAux2disc,
		  boost::numeric::ublas::range(0,cAux2disc.size1()),
		  boost::numeric::ublas::range(0,cAux2disc.size2())); 

	  Termo1=sum(cdiag1).real();


	  if (printstuff)
	    cout << " ComplexOperator. Termo1= " << Termo1 
		 << endl;
	  // if printstuff


	}
	// end if Op1 is REAL else Comples
 
      }
      // end calculate termo1

      // Termo2
      double Termo2=0.0;
//       if ( (Nst_bl1_disc!=0)&&(Nst_bl2_kept!=0)&&(omegabar<0.0) ){
      if ( (Nst_bl1_disc!=0)&&(Nst_bl2_kept!=0)&&((omegabar<0.0)||(CalcNorm)) ){
	// REAL operators
	if ( !(Op1N[Nshell].IsComplex) ){

	  boost::numeric::ublas::matrix<double> Rho2(Nst_bl2_kept,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<double> opA2(Nst_bl1_disc,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<double> opA2w(Nst_bl1_disc,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<double> opB2(Nst_bl2_kept,Nst_bl1_disc);

	  // Can I use opB1 and opB2 here instead? (also: do I need two opA matrices??) Check...
	  boost::numeric::ublas::matrix<double> Prod2(Nst_bl2_kept,Nst_bl1_disc);

	  Rho2=RhoN[Nshell].MatBlock2BLAS(ibl2,ibl2,true,true);

	  opA2=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,true);
	  if (NonDiagGF)
	    //opB=Op2N[Nshell].MatBlock2BLAS(ibl2,ibl1);
	    opB2=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,true)); 
	  // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
	  else
	    opB2=trans(opA2);

	  noalias(Prod2)=prod(Rho2,opB2); // Nb2_kept x Nb1_disc

	  // Calculate A_ij * Delta(w-(Ej-Ei))

	  if (CalcNorm) opA2w=opA2; else
	    noalias(opA2w)=MijxBDeltaEij(opA2,(&AcutN[Nshell]),ibl1,ibl2,false,true,omegabar); // Nb1_disc x Nb2_kept

	  // Take the trace: (use preallocated Rho2 as temp...)

	  noalias(Rho2)=prod(Prod2,opA2w);  // Nb2_kept x Nb2_kept

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	    diag2(Rho2,
		  boost::numeric::ublas::range(0,Rho2.size1()),
		  boost::numeric::ublas::range(0,Rho2.size2())); 

	  Termo2=sum(diag2);

	  // COMPLEX operators
	} else {

	  boost::numeric::ublas::matrix<complex<double> > cRho2(Nst_bl2_kept,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<complex<double> > copA2(Nst_bl1_disc,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<complex<double> > copA2w(Nst_bl1_disc,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<complex<double> > copB2(Nst_bl2_kept,Nst_bl1_disc);

	  // Can I use opB1 and opB2 here instead? (also: do I need two opA matrices??) Check...
	  boost::numeric::ublas::matrix<complex<double> > cProd2(Nst_bl2_kept,Nst_bl1_disc);

	  cRho2=RhoN[Nshell].cMatBlock2BLAS(ibl2,ibl2,true,true);

	  copA2=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,true);
	  if (NonDiagGF)
	    copB2=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,true)); 
	  // Nst_bl2_disc x Nst_bl1_kept matrix
	  else
	    copB2=trans(copA2);

	  noalias(cProd2)=prod(cRho2,copB2); // Nb2_kept x Nb1_disc

	  // Calculate A_ij * Delta(w-(Ej-Ei))

	  if (CalcNorm) copA2w=copA2; else
	    noalias(copA2w)=cMijxBDeltaEij(copA2,(&AcutN[Nshell]),ibl1,ibl2,false,true,omegabar); // Nb1_disc x Nb2_kept


	  // Take the trace Use cRho2 as aux // Nb2_kept x Nb2_kept:
	  // aux
// 	  boost::numeric::ublas::matrix<complex<double> > cAux2kept(Nst_bl2_kept,Nst_bl2_kept);
	  noalias(cRho2)=prod(cProd2,copA2w);  // Nb2_kept x Nb2_kept

// 	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
// 	    cdiag2(cAux2kept,
// 		   boost::numeric::ublas::range(0,cAux2kept.size1()),
// 		   boost::numeric::ublas::range(0,cAux2kept.size2())); 
	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	    cdiag2(cRho2,
		   boost::numeric::ublas::range(0,cRho2.size1()),
		   boost::numeric::ublas::range(0,cRho2.size2())); 

	  Termo2=sum(cdiag2).real();

	}
	// end if Op1 is REAL else Comples

      }
      // end calc Termo2

      // SU(2) symmetry: <ibl1|A|ibl2> (or <ibl2|B|ibl1> )
      // Given ibl1, ibl2, calculate the following scalars
      // Get Sibl1,Sibl2
      // 
      // OpA: auxCG1[0]= sum_Szilb1 CGordan(Sibl1, Szibl1, 0.5, sigma, Sibl2, Szibl2=sigma+Szibl1)

      // OpB: auxCG1[1]= sum_Szilb1 CGordan(Sibl1, Szibl1, 0.5, sigma, Sibl2, Szibl2=sigma+Szibl1)
      //  They are the same!!
      double CGfactor=1.0;
      // Get S from ibl
      double Sbl1=0.0;
      double Sbl2=0.0;
      if (AcutN[Nshell].totalS){
	CGfactor=0.0;
	Sbl1=AcutN[Nshell].GetQNumber(ibl1,AcutN[Nshell].Sqnumbers[0]); 
	Sbl2=AcutN[Nshell].GetQNumber(ibl2,AcutN[Nshell].Sqnumbers[0]); 
	// only a single SU(2) for now
	for (double Szbl1=-Sbl1;Szbl1<=Sbl1;Szbl1+=1.0){
	  double dSigma=0.5; 
	  double Szbl2=Szbl1+dSigma;
	  double auxCG=CGordan(Sbl1,Szbl1, 0.5, dSigma, Sbl2, Szbl2);
	  CGfactor+=auxCG*auxCG;
	}
	if (dEqual(CGfactor,0.0)){cout << "Ops. CGfactor = 0.0 " << endl;}
      }
      // end if totalS

//       rho_wM+=CGfactor*Termo1;

//       rho_wM+=CGfactor*Termo2;

      Positive_w+=CGfactor*Termo1;

      Negative_w+=CGfactor*Termo2;

    }
    // end if there are kept states in blocks

  }
  // end loop in Op1N blocks

  if (printstuff)
    cout << " Positive_w = " << Positive_w
	 << " Negative_w = " << Negative_w
	 << endl;

  double BosonSign=1.0;
  // if (IsBosonic) BosonSign=-1.0;
  if (CalcNorm){rho_wM=Positive_w+(BosonSign)*Negative_w;}else{
    if (omegabar>=0.0){rho_wM=Positive_w;}else{rho_wM=Negative_w;}
  } 
  // end if CalcNorm

  return(rho_wM);

}

//////////////////
///////////////////
///////////////////
///////////////////

// REAL 

boost::numeric::ublas::matrix<double> CSpecFunction::MijxBDeltaEij(CNRGmatrix* pMat, CNRGarray* pAeig, int ibl1, int ibl2, bool kp1, bool kp2, double omega){

  bool ChkSync=pMat->ChkSync(pAeig);
  int iMatBl=pMat->FindMatBlock(ibl1,ibl2);

  if ( (!ChkSync)||(iMatBl<0) ){
    cout << " Cannot calculate Mij*BD(Ei-Ej) " << endl;
    return(boost::numeric::ublas::matrix<double> () );
  }

  boost::numeric::ublas::matrix<double> Maux=pMat->MatBlock2BLAS(ibl1,ibl2,kp1,kp2);

  return( MijxBDeltaEij(Maux,pAeig,ibl1,ibl2,kp1,kp2,omega) );

}
//

boost::numeric::ublas::matrix<double> CSpecFunction::MijxBDeltaEij(boost::numeric::ublas::matrix<double> Mat, CNRGarray* pAeig, int ibl1, int ibl2,  bool kp1, bool kp2, double omega){


  int Nst1=pAeig->GetBlockSize(ibl1);
  int Nst2=pAeig->GetBlockSize(ibl2);


  int Nst1_kp1=pAeig->GetBlockSize(ibl1,kp1);
  int Nst2_kp2=pAeig->GetBlockSize(ibl2,kp2);


  if ( (Nst1_kp1!=Mat.size1())||(Nst2_kp2!=Mat.size2())
       ||(Nst1_kp1<=0)||(Nst2_kp2<=0) ){
    cout << " Cannot calculate Mij_kpi,kpi*BD(Ei-Ej) " << endl;
    return(boost::numeric::ublas::matrix<double> () );
  }

  boost::numeric::ublas::matrix<double> Maux(Nst1_kp1,Nst2_kp2);

  int ist0=pAeig->GetBlockLimit(ibl1,0);
  int jst0=pAeig->GetBlockLimit(ibl2,0);

  double Ei=0.0;
  double Ej=0.0;

  bool print=false;

  //if ( (ibl1==2)&&(ibl2==3) ) print=true;

  if (print) cout << " omega = " << omega 
		  << " dBroad = " << dBroad <<endl;

  double BD=0.0;
//   double TempDN=Temp/CalcDN(pAeig->Nshell);
  double ThisDN=CalcDN(pAeig->Nshell);
  double TempDN=Temp/ThisDN;
  double GapDN=Gap/ThisDN;


  int istkp=0;
  for (int ii=0;ii<Nst1;ii++){
    if (pAeig->CheckKept(ist0+ii,kp1)){
      Ei=pAeig->dEn[ist0+ii];
      if (print) cout << " Ei = " << Ei;
      int jstkp=0;
      for (int jj=0;jj<Nst2;jj++){
	if (pAeig->CheckKept(jst0+jj,kp2)){
	Ej=pAeig->dEn[jst0+jj];
	if (print) cout << " Ej = " << Ej;
	// NEED TO INCLUDE TEMPERATURE HERE!
	//BD=BDelta(omega,(Ej-Ei),dBroad);
	// No double counting. We can do this:
	BD=BDelta(omega,Ej-Ei,dBroad);
	// NEW (Jan 2016): No broadening if it is inside the Gap
	if (fabs(Ej-Ei)<GapDN){
// 	  cout << " MijxBDeltaEij: Inside the gap. Ei= " 
// 	       << Ei << " Ej= " << Ej << endl;
	  //BD=0.0; // Will give zero... Is this correct??
	}
	if (fabs(Ej-Ei)<TwindowFac*TempDN){
	  //BD=BDeltaTemp(omega,Ej-Ei,dBroadTemp*TempDN); //
	  BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp); // Anders
	  //  	cout << " |Ej-Ei|<TempDN: |Ej-Ei|=" << fabs(Ej-Ei) 
	  //  	     << " TempDN = " << TempDN << endl;
	}
	// end if |Ej-Ei|<TempDN
	if (print) cout << " BD =  " << BD;
	Maux(istkp,jstkp)=Mat(istkp,jstkp)*BD; // Uncomment this to work
	jstkp++;
	// END TEMP
	}
	// end check kept in state j
      }
      if (pAeig->CheckKept(ist0+ii,kp1)){istkp++;}
      if (print) cout << endl;
    }
    // end check kept in state i
  }
  // end loop in i,j

  return(Maux);

}
//



///// COMPLEX

boost::numeric::ublas::matrix<complex<double> > CSpecFunction::cMijxBDeltaEij(CNRGmatrix* pMat, CNRGarray* pAeig, int ibl1, int ibl2, bool kp1, bool kp2, double omega){

  bool ChkSync=pMat->ChkSync(pAeig);
  int iMatBl=pMat->FindMatBlock(ibl1,ibl2);

  if ( (!ChkSync)||(iMatBl<0) ){
    cout << " Cannot calculate Mij*BD(Ei-Ej) " << endl;
    return(boost::numeric::ublas::matrix<double> () );
  }

  boost::numeric::ublas::matrix<complex<double> > cMaux=pMat->cMatBlock2BLAS(ibl1,ibl2,kp1,kp2);

  return( cMijxBDeltaEij(cMaux,pAeig,ibl1,ibl2,kp1,kp2,omega) );

}
//

boost::numeric::ublas::matrix<complex<double> > CSpecFunction::cMijxBDeltaEij(boost::numeric::ublas::matrix<complex<double> > cMat, CNRGarray* pAeig, int ibl1, int ibl2,  bool kp1, bool kp2, double omega){


  int Nst1=pAeig->GetBlockSize(ibl1);
  int Nst2=pAeig->GetBlockSize(ibl2);


  int Nst1_kp1=pAeig->GetBlockSize(ibl1,kp1);
  int Nst2_kp2=pAeig->GetBlockSize(ibl2,kp2);


  if ( (Nst1_kp1!=cMat.size1())||(Nst2_kp2!=cMat.size2())
       ||(Nst1_kp1<=0)||(Nst2_kp2<=0) ){
    cout << " Cannot calculate Mij_kpi,kpi*BD(Ei-Ej) " << endl;
    return(boost::numeric::ublas::matrix<double> () );
  }

  boost::numeric::ublas::matrix<complex<double> > cMaux(Nst1_kp1,Nst2_kp2);

  int ist0=pAeig->GetBlockLimit(ibl1,0);
  int jst0=pAeig->GetBlockLimit(ibl2,0);

  double Ei=0.0;
  double Ej=0.0;

  bool print=false;

  //if ( (ibl1==2)&&(ibl2==3) ) print=true;

  if (print) cout << " omega = " << omega 
		  << " dBroad = " << dBroad <<endl;

  double BD=0.0;
//   double TempDN=Temp/CalcDN(pAeig->Nshell);
  double ThisDN=CalcDN(pAeig->Nshell);
  double TempDN=Temp/ThisDN;
  double GapDN=Gap/ThisDN;



  int istkp=0;
  for (int ii=0;ii<Nst1;ii++){
    if (pAeig->CheckKept(ist0+ii,kp1)){
      Ei=pAeig->dEn[ist0+ii];
      if (print) cout << " Ei = " << Ei;
      int jstkp=0;
      for (int jj=0;jj<Nst2;jj++){
	if (pAeig->CheckKept(jst0+jj,kp2)){
	Ej=pAeig->dEn[jst0+jj];
	if (print) cout << " Ej = " << Ej;
	// NEED TO INCLUDE TEMPERATURE HERE!
	//BD=BDelta(omega,(Ej-Ei),dBroad);
	// No double counting. We can do this:
	BD=BDelta(omega,Ej-Ei,dBroad);
	// NEW (Jan 2016): No broadening if it is inside the Gap
	if (fabs(Ej-Ei)<GapDN) BD=0.0; // Will give zero...
	if (fabs(Ej-Ei)<TwindowFac*TempDN){
	  //BD=BDeltaTemp(omega,Ej-Ei,dBroadTemp*TempDN); //
	  BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp); // Anders
	  //  	cout << " |Ej-Ei|<TempDN: |Ej-Ei|=" << fabs(Ej-Ei) 
	  //  	     << " TempDN = " << TempDN << endl;
	}
	// end if |Ej-Ei|<TempDN
	if (print) cout << " BD =  " << BD;
	cMaux(istkp,jstkp)=cMat(istkp,jstkp)*BD; // Uncomment this to work
	jstkp++;
	// END TEMP
	}
	// end check kept in state j
      }
      if (pAeig->CheckKept(ist0+ii,kp1)){istkp++;}
      if (print) cout << endl;
    }
    // end check kept in state i
  }
  // end loop in i,j

  return(cMaux);

}
//



// Full Density-Matrix calculations


//////////////////
double CSpecFunction::FDM_SpecDens_M(double omegabar,int Nshell, bool CalcNorm){

  // Need to get the Connecting matrix elements and energies
  // to and from a given block (say, blGS)
  // Interesting procedure...

  // Loop in blocks of Operator


  double rho_wM=0.0;

  if (Nshell>NshellMax-1){
    cout << "FDM_SpecDens_M: Nshell > NshellMax . Return 0" << endl;
    return(0.0);
  }

  if (Nshell==NshellMax-1){
    if (dEqual(Temp,0.0))
      Betabar=1.0e10;
//     cout << "FDM_SpecDens_M: Last Nshell: Costi. CalcNorm= " << CalcNorm << endl;
//     cout << "FDM_SpecDens_M: T=0 calculation (Costi)." << endl;
//     cout << " Betabar = " << Betabar << endl;
    rho_wM=CalcSpecCosti_T_N(Nshell,Betabar,omegabar,CalcNorm);
    return(rho_wM);
  }
  //

  double Positive_w=0.0;
  double Negative_w=0.0;

  bool printstuff=false;

  // Need to Check Sync!!
  bool ChkSyncRhoAcut=RhoN[Nshell].ChkSync((&AcutN[Nshell]));
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));

  if ( (!ChkSyncRhoAcut)||(!ChkSyncOp1Acut) ){
    cout << "FDM_SpecDens_M: RhoN or Op1N not in sync with AcutN " << endl;
    return(0.0);
  }

  if (printstuff)
    cout << " FDM_SpecDens_M: Nshell = " << Nshell
	 << " omegabar = " << omegabar
	 << " CalcNorm = " << CalcNorm
	 << endl;


  for (int iMatbl=0;iMatbl<Op1N[Nshell].NumMatBlocks();iMatbl++){
  // Checking a few blocks first
  //for (int iMatbl=0;iMatbl<3;iMatbl++){

    int ibl1=0;
    int ibl2=0;

    Op1N[Nshell].GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);


    // Ok, remember the input Op2N is the COMPLEX CONJUGATE of what we want!
    //  
    //  So (B)_{ibl2 ibl1} = ((B+)_{ibl1 ibl2})* <- THIS is Op2

    //if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl2,ibl1)<0) ){
    if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl1,ibl2)<0) ){
      cout << "FDM_SpecDens_M: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return(0.0);
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 
 
    // Check block sizes
    if ( (Nst_bl1!=RhoN[Nshell].GetBlockSize(ibl1))||
	 (Nst_bl2!=RhoN[Nshell].GetBlockSize(ibl2)) ){
      cout << "FDM_SpecDens_M: Block in Op1 not the same size as block in RhoN " << endl;
      return(0.0);
    }
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "FDM_SpecDens_M: Block in Op1 not the same size as block in Op2 " << endl;
      return(0.0);
    }

    int Nst_bl1_kept=Op1N[Nshell].GetBlockSize(ibl1,true); 
    int Nst_bl1_disc=Op1N[Nshell].GetBlockSize(ibl1,false);  

    int Nst_bl2_kept=Op1N[Nshell].GetBlockSize(ibl2,true); 
    int Nst_bl2_disc=Op1N[Nshell].GetBlockSize(ibl2,false); 

    if (printstuff){
      cout << " ibl1 = " << ibl1
	   << " ibl2 = " << ibl2
	   << endl;

      cout << " Nst_bl1= " << Nst_bl1 
	   << " Nst_bl1_kept= " << Nst_bl1_kept 
	   << " Nst_bl1_disc=" << Nst_bl1_disc	
	   << endl
	   << " Nst_bl2= " << Nst_bl2 
	   << " Nst_bl2_kept= " << Nst_bl2_kept 
	   << " Nst_bl2_disc=" << Nst_bl2_disc
	   << endl;

    }
    // end printstuff


    if ( (Nst_bl1!=Nst_bl1_kept+Nst_bl1_disc)||
	 (Nst_bl2!=Nst_bl2_kept+Nst_bl2_disc) ) {
      cout << "FDM_SpecDens_M: Error: Kept+Dics NOT equal to Nstblock! " << endl;
      return(0.0);
    }

    // Gotta check if ANY of these are zero!

    if (  ( (Nst_bl1_kept==0)&&(Nst_bl2_kept==0) )||
	  ( (Nst_bl1_disc==0)&&(Nst_bl2_disc==0) )   ){
      if (printstuff) cout << " Blocks " << ibl1 << " and " << ibl2 
	   << " are either both kept or both discarded. Moving on..." << endl;
    }else{

      // Blas matrices

      // Termo1KK
      double Termo1KK=0.0;
      if ( (Nst_bl1_kept!=0)&&(Nst_bl2_disc!=0)&&((omegabar>=0.0)||(CalcNorm)) ){
	// caculate term 1
	  if (printstuff)
	    cout << "  Calculating Termo1KK " << endl;
	  // if printstuff


	// REAL operators
	if ( ( !(Op1N[Nshell].IsComplex) )&&( !(RhoN[Nshell].IsComplex) ) ) {

	  boost::numeric::ublas::matrix<double> Rho1(Nst_bl1_kept,Nst_bl1_kept);

	  boost::numeric::ublas::matrix<double> opA1w(Nst_bl1_kept,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> opA1(Nst_bl1_kept,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> opB1(Nst_bl2_disc,Nst_bl1_kept);
	  // aux
	  boost::numeric::ublas::matrix<double> Aux2disc(Nst_bl2_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> Prod1(Nst_bl2_disc,Nst_bl1_kept);

	  Rho1=RhoN[Nshell].MatBlock2BLAS(ibl1,ibl1,true,true);

	  opA1=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2,true,false);
	  // Nst_bl1_kept x Nst_bl2_disc matrix (assuming real!)


	  if (NonDiagGF)
	    opB1=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2,true,false)); // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
	  else
	    opB1=trans(opA1);

	  // Positive and negative energies matrix elements
	  // Note for the future: sign is "-" if bosonic operators
	  //noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix

	  noalias(Prod1)=prod(opB1,Rho1); // Nb2_disc x Nb1_kept

	  if (CalcNorm) opA1w=opA1; else
	    noalias(opA1w)=MijxBDeltaEij(opA1,(&AcutN[Nshell]),ibl1,ibl2,true,false,omegabar); // Nb1_kept x Nb2_disc

	  noalias(Aux2disc)=prod(Prod1,opA1w);  // Nb2_disc x Nb2_disc

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	    diag1(Aux2disc,
		  boost::numeric::ublas::range(0,Aux2disc.size1()),
		  boost::numeric::ublas::range(0,Aux2disc.size2())); 

	  Termo1KK=sum(diag1);


	// COMPLEX operators
	} else {

	  if (printstuff)
	    cout << " ComplexOperator. Termo1KK " << endl;
	  // if printstuff

	  boost::numeric::ublas::matrix<complex<double> > cRho1(Nst_bl1_kept,Nst_bl1_kept);

	  boost::numeric::ublas::matrix<complex<double> > copA1w(Nst_bl1_kept,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > copA1(Nst_bl1_kept,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > copB1(Nst_bl2_disc,Nst_bl1_kept);
	  // aux
	  boost::numeric::ublas::matrix<complex<double> > cAux2disc(Nst_bl2_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > cProd1(Nst_bl2_disc,Nst_bl1_kept);


	  cRho1=RhoN[Nshell].cMatBlock2BLAS(ibl1,ibl1,true,true);

	  copA1=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2,true,false);
	  // Nst_bl1_kept x Nst_bl2_disc matrix (assuming real!)

	  if (NonDiagGF)
	    //opB=Op2N[Nshell].MatBlock2BLAS(ibl2,ibl1);
	    copB1=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2,true,false)); // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
	  else
	    copB1=trans(copA1);

	  // Positive and negative energies matrix elements
	  // Note for the future: sign is "-" if bosonic operators

	  noalias(cProd1)=prod(copB1,cRho1); // Nb2_disc x Nb1_kept

	  if (CalcNorm) copA1w=copA1; else
	    noalias(copA1w)=cMijxBDeltaEij(copA1,(&AcutN[Nshell]),ibl1,ibl2,true,false,omegabar); // Nb1_kept x Nb2_disc

	  noalias(cAux2disc)=prod(cProd1,copA1w);  // Nb2_disc x Nb2_disc

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	    cdiag1(cAux2disc,
		  boost::numeric::ublas::range(0,cAux2disc.size1()),
		  boost::numeric::ublas::range(0,cAux2disc.size2())); 

	  Termo1KK=sum(cdiag1).real();


	  if (printstuff)
	    cout << " ComplexOperator. Termo1= " << Termo1KK 
		 << endl;
	  // if printstuff


	}
	// end if Op1 is REAL else Comples
 
      }
      // end calculate termo1KK


      // Termo1DDADD and Termo1DDADK
      double Termo1DDADD=0.0;
      double Termo1DDADK=0.0;
      if ( (Nst_bl1_disc!=0)&&(Nst_bl2_disc!=0) ){
	// caculate term 1

	if (printstuff)
	  cout << "  Calculating Termo1DDADD e Termo1DDADK " << endl;
	// if printstuff


	// REAL operators
	if ( ( !(Op1N[Nshell].IsComplex) )&&( !(RhoN[Nshell].IsComplex) ) ) {

	  boost::numeric::ublas::matrix<double> Rho1dd(Nst_bl1_disc,Nst_bl1_disc);

	  boost::numeric::ublas::matrix<double> opA1ddw(Nst_bl1_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> opA1dd(Nst_bl1_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> opB1dd(Nst_bl2_disc,Nst_bl1_disc);
	  // aux
	  boost::numeric::ublas::matrix<double> Aux2dd(Nst_bl2_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> Prod1dd(Nst_bl2_disc,Nst_bl1_disc);


	  Rho1dd=RhoN[Nshell].MatBlock2BLAS(ibl1,ibl1,false,false);
	  // Nst_bl1_disc x Nst_bl2_disc matrix (assuming real!)

	  opA1dd=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,false);

	  if (NonDiagGF)
	    opB1dd=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,false));
	  // Nst_bl2_disc x Nst_bl1_disc matrix (assuming real!)
	  else
	    opB1dd=trans(opA1dd);

	  noalias(Prod1dd)=prod(opB1dd,Rho1dd); // Nb2_disc x Nb1_disc

	  if (CalcNorm)
	    opA1ddw=opA1dd; 
	  else
	    noalias(opA1ddw)=MijxBDeltaEij(opA1dd,(&AcutN[Nshell]),ibl1,ibl2,false,false,omegabar); // Nb1_disc x Nb2_disc
	  
	  noalias(Aux2dd)=prod(Prod1dd,opA1ddw);  // Nb2_disc x Nb2_disc

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	    diagDD(Aux2dd,
		   boost::numeric::ublas::range(0,Aux2dd.size1()),
		   boost::numeric::ublas::range(0,Aux2dd.size2())); 

	  
	  Termo1DDADD=sum(diagDD);
	  
	  
	  if (Nst_bl2_kept!=0){
	    // Added condition if Nst_bl2_kept!=0 !
	    boost::numeric::ublas::matrix<double> opA1dkw(Nst_bl1_disc,Nst_bl2_kept);
	    boost::numeric::ublas::matrix<double> opA1dk(Nst_bl1_disc,Nst_bl2_kept);
	    boost::numeric::ublas::matrix<double> opB1kd(Nst_bl2_kept,Nst_bl1_disc);
	    // aux
	    boost::numeric::ublas::matrix<double> Aux2kk(Nst_bl2_kept,Nst_bl2_kept);
	    boost::numeric::ublas::matrix<double> Prod1kd(Nst_bl2_kept,Nst_bl1_disc);

	    opA1dk=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,true);

	    if (NonDiagGF)
	      opB1kd=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2,true,false));
	    // Nst_bl2_kept x Nst_bl1_disc matrix (assuming real!)
	    else
	      opB1kd=trans(opA1dk);
	  

	    // Positive and negative energies matrix elements
	    // Note for the future: sign is "-" if bosonic operators

	    noalias(Prod1kd)=prod(opB1kd,Rho1dd); // Nb2_kept x Nb1_disc

	    if (CalcNorm)
	      opA1dkw=opA1dk; 
	    else
	      noalias(opA1dkw)=MijxBDeltaEij(opA1dk,(&AcutN[Nshell]),ibl1,ibl2,false,true,omegabar); // Nb1_disc x Nb2_kept

	    noalias(Aux2kk)=prod(Prod1kd,opA1dkw);  // Nb2_kept x Nb2_kept


	    boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	      diagKK(Aux2kk,
		     boost::numeric::ublas::range(0,Aux2kk.size1()),
		     boost::numeric::ublas::range(0,Aux2kk.size2())); 

	    Termo1DDADK=sum(diagKK);
  
	  }
	  // end if Nst_ketp_bl2 neq 0


	// COMPLEX operators
	} else {

	  if (printstuff)
	    cout << " ComplexOperator. Termo1DDADD and Termo1DDADK " 
		 << endl;
	  // if printstuff

	  boost::numeric::ublas::matrix<complex<double> > cRho1dd(Nst_bl1_disc,Nst_bl1_disc);

	  boost::numeric::ublas::matrix<complex<double> > copA1ddw(Nst_bl1_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > copA1dd(Nst_bl1_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > copB1dd(Nst_bl2_disc,Nst_bl1_disc);
	  // aux
	  boost::numeric::ublas::matrix<complex<double> > cAux2dd(Nst_bl2_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > cProd1dd(Nst_bl2_disc,Nst_bl1_disc);

	  cRho1dd=RhoN[Nshell].cMatBlock2BLAS(ibl1,ibl1,false,false);
	  // Nst_bl1_disc x Nst_bl2_disc matrix (assuming real!)

	  copA1dd=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,false);


	  if (NonDiagGF)
	    copB1dd=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,false));
	  // Nst_bl2_disc x Nst_bl1_disc matrix (assuming real!)
	  else
	    copB1dd=trans(copA1dd);
	    

	  // Positive and negative energies matrix elements
	  // Note for the future: sign is "-" if bosonic operators
	  //noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix

	  noalias(cProd1dd)=prod(copB1dd,cRho1dd); // Nb2_disc x Nb1_disc

	  if (CalcNorm){
	    copA1ddw=copA1dd; 
	  }else{
	    noalias(copA1ddw)=cMijxBDeltaEij(copA1dd,(&AcutN[Nshell]),ibl1,ibl2,false,false,omegabar); // Nb1_disc x Nb2_disc
	  }
	  noalias(cAux2dd)=prod(cProd1dd,copA1ddw);  // Nb2_disc x Nb2_disc

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	    cdiagDD(cAux2dd,
		  boost::numeric::ublas::range(0,cAux2dd.size1()),
		  boost::numeric::ublas::range(0,cAux2dd.size2())); 


	  Termo1DDADD=sum(cdiagDD).real();

	  if (Nst_bl2_kept!=0){
	    // Added condition if Nst_bl2_kept!=0 !
	    boost::numeric::ublas::matrix<complex<double> > copA1dkw(Nst_bl1_disc,Nst_bl2_kept);
	    boost::numeric::ublas::matrix<complex<double> > copA1dk(Nst_bl1_disc,Nst_bl2_kept);
	    boost::numeric::ublas::matrix<complex<double> > copB1kd(Nst_bl2_kept,Nst_bl1_disc);

	    // aux
	    boost::numeric::ublas::matrix<complex<double> > cAux2kk(Nst_bl2_kept,Nst_bl2_kept);
	    boost::numeric::ublas::matrix<complex<double> > cProd1kd(Nst_bl2_kept,Nst_bl1_disc);

	    copA1dk=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,true);

	    if (NonDiagGF)
	      copB1kd=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2,true,false));
	    // Nst_bl2_kept x Nst_bl1_disc matrix (assuming real!)
	    else
	      copB1kd=trans(copA1dk);

	    noalias(cProd1kd)=prod(copB1kd,cRho1dd); // Nb2_kept x Nb1_disc

	    if (CalcNorm){
	      copA1dkw=copA1dk; 
	    }else{
	      noalias(copA1dkw)=cMijxBDeltaEij(copA1dk,(&AcutN[Nshell]),ibl1,ibl2,false,true,omegabar); // Nb1_disc x Nb2_kept
	    }
	    noalias(cAux2kk)=prod(cProd1kd,copA1dkw);  // Nb2_kept x Nb2_kept

	    boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	      cdiagKK(cAux2kk,
		      boost::numeric::ublas::range(0,cAux2kk.size1()),
		      boost::numeric::ublas::range(0,cAux2kk.size2())); 

	    Termo1DDADK=sum(cdiagKK).real();
	    
	  }
	  // end if Nst_bl2_kept!=0
	  
	  if (printstuff)
	    cout << " ComplexOperator. Termo1DDADD= " << Termo1DDADD
		 << endl;
	  // if printstuff
	}
	// end if Op1 is REAL else Comples
      }
      // end calculate termo1DD


      // Termo2KK 
      double Termo2KK=0.0;

      if ( (Nst_bl1_disc!=0)&&(Nst_bl2_kept!=0)&&((omegabar<0.0)||(CalcNorm)) ){

	if (printstuff)
	  cout << "  Calculating Termo2KK " << endl;
	// if printstuff


	// REAL operators
	if ( !(Op1N[Nshell].IsComplex) ){

	  boost::numeric::ublas::matrix<double> Rho2(Nst_bl2_kept,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<double> opA2(Nst_bl1_disc,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<double> opA2w(Nst_bl1_disc,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<double> opB2(Nst_bl2_kept,Nst_bl1_disc);

	  // Can I use opB1 and opB2 here instead? (also: do I need two opA matrices??) Check...
	  boost::numeric::ublas::matrix<double> Prod2(Nst_bl2_kept,Nst_bl1_disc);

	  Rho2=RhoN[Nshell].MatBlock2BLAS(ibl2,ibl2,true,true);

	  opA2=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,true);
	  if (NonDiagGF)
	    //opB=Op2N[Nshell].MatBlock2BLAS(ibl2,ibl1);
	    opB2=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,true)); 
	  // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
	  else
	    opB2=trans(opA2);

	  noalias(Prod2)=prod(Rho2,opB2); // Nb2_kept x Nb1_disc

	  // Calculate A_ij * Delta(w-(Ej-Ei))

	  if (CalcNorm) opA2w=opA2; else
	    noalias(opA2w)=MijxBDeltaEij(opA2,(&AcutN[Nshell]),ibl1,ibl2,false,true,omegabar); // Nb1_disc x Nb2_kept

	  // Take the trace: (use preallocated Rho2 as temp...)

	  noalias(Rho2)=prod(Prod2,opA2w);  // Nb2_kept x Nb2_kept

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	    diag2(Rho2,
		  boost::numeric::ublas::range(0,Rho2.size1()),
		  boost::numeric::ublas::range(0,Rho2.size2())); 

	  Termo2KK=sum(diag2);

	  // COMPLEX operators
	} else {

	  boost::numeric::ublas::matrix<complex<double> > cRho2(Nst_bl2_kept,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<complex<double> > copA2(Nst_bl1_disc,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<complex<double> > copA2w(Nst_bl1_disc,Nst_bl2_kept);
	  boost::numeric::ublas::matrix<complex<double> > copB2(Nst_bl2_kept,Nst_bl1_disc);

	  // Can I use opB1 and opB2 here instead? (also: do I need two opA matrices??) Check...
	  boost::numeric::ublas::matrix<complex<double> > cProd2(Nst_bl2_kept,Nst_bl1_disc);

	  cRho2=RhoN[Nshell].cMatBlock2BLAS(ibl2,ibl2,true,true);

	  copA2=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,true);
	  if (NonDiagGF)
	    copB2=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,true)); 
	  // Nst_bl2_disc x Nst_bl1_kept matrix
	  else
	    copB2=trans(copA2);

	  noalias(cProd2)=prod(cRho2,copB2); // Nb2_kept x Nb1_disc

	  // Calculate A_ij * Delta(w-(Ej-Ei))

	  if (CalcNorm) copA2w=copA2; else
	    noalias(copA2w)=cMijxBDeltaEij(copA2,(&AcutN[Nshell]),ibl1,ibl2,false,true,omegabar); // Nb1_disc x Nb2_kept


	  // Take the trace Use cRho2 as aux // Nb2_kept x Nb2_kept:
	  // aux
// 	  boost::numeric::ublas::matrix<complex<double> > cAux2kept(Nst_bl2_kept,Nst_bl2_kept);
	  noalias(cRho2)=prod(cProd2,copA2w);  // Nb2_kept x Nb2_kept

// 	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
// 	    cdiag2(cAux2kept,
// 		   boost::numeric::ublas::range(0,cAux2kept.size1()),
// 		   boost::numeric::ublas::range(0,cAux2kept.size2())); 
	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	    cdiag2(cRho2,
		   boost::numeric::ublas::range(0,cRho2.size1()),
		   boost::numeric::ublas::range(0,cRho2.size2())); 

	  Termo2KK=sum(cdiag2).real();

	}
	// end if Op1 is REAL else Comples

      }
      // end calc Termo2


      // Termo2DDADD and Termo2DDAKD
      double Termo2DDADD=0.0;
      double Termo2DDAKD=0.0;
      if ( (Nst_bl1_disc!=0)&&(Nst_bl2_disc!=0) ){
	// caculate term 2

	if (printstuff)
	  cout << "  Calculating Termo2DDADD e Termo2DDAKD " << endl;
	// if printstuff


	// REAL operators
	if ( ( !(Op1N[Nshell].IsComplex) )&&( !(RhoN[Nshell].IsComplex) ) ) {

	  // Check sizes! Just put it here...
	  
	  boost::numeric::ublas::matrix<double> Rho2dd(Nst_bl2_disc,Nst_bl2_disc);

	  boost::numeric::ublas::matrix<double> opA2ddw(Nst_bl1_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> opA2dd(Nst_bl1_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> opB2dd(Nst_bl2_disc,Nst_bl1_disc);
	  // aux
	  boost::numeric::ublas::matrix<double> Aux2dd(Nst_bl2_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<double> Prod2dd(Nst_bl2_disc,Nst_bl1_disc);


	  Rho2dd=RhoN[Nshell].MatBlock2BLAS(ibl2,ibl2,false,false);
	  // Nst_bl2_disc x Nst_bl2_disc matrix (assuming real!)

	  opA2dd=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,false);

	  if (NonDiagGF)
	    opB2dd=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2,false,false));
	  // Nst_bl2_disc x Nst_bl1_disc matrix (assuming real!)
	  else
	    opB2dd=trans(opA2dd);

	  noalias(Prod2dd)=prod(Rho2dd,opB2dd); // Nb2_disc x Nb1_disc

	  if (CalcNorm)
	    opA2ddw=opA2dd; 
	  else
	    noalias(opA2ddw)=MijxBDeltaEij(opA2dd,(&AcutN[Nshell]),ibl1,ibl2,false,false,omegabar); // Nb1_disc x Nb2_disc
	  
	  noalias(Aux2dd)=prod(Prod2dd,opA2ddw);  // Nb2_disc x Nb2_disc

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	    diagDD(Aux2dd,
		   boost::numeric::ublas::range(0,Aux2dd.size1()),
		   boost::numeric::ublas::range(0,Aux2dd.size2())); 

	  
	  Termo2DDADD=sum(diagDD);
	  
	  
	  if (Nst_bl1_kept!=0){
	    // Added condition if Nst_bl2_kept!=0 !
	    boost::numeric::ublas::matrix<double> opA2kdw(Nst_bl1_kept,Nst_bl2_disc);
	    boost::numeric::ublas::matrix<double> opA2kd(Nst_bl1_kept,Nst_bl2_disc);
	    boost::numeric::ublas::matrix<double> opB2dk(Nst_bl2_disc,Nst_bl1_kept);
	    // aux
	    boost::numeric::ublas::matrix<double> Prod2dk(Nst_bl2_disc,Nst_bl1_kept);
	    boost::numeric::ublas::matrix<double> Aux2dd(Nst_bl2_disc,Nst_bl2_disc);

	    opA2kd=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2,true,false);

	    if (NonDiagGF)
	      opB2dk=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2,true,false));
	    // Nst_bl2_kept x Nst_bl1_disc matrix (assuming real!)
	    else
	      opB2dk=trans(opA2kd);
	  

	    // Positive and negative energies matrix elements
	    // Note for the future: sign is "-" if bosonic operators

	    noalias(Prod2dk)=prod(Rho2dd,opB2dk); // Nb2_disc x Nb1_kept

	    if (CalcNorm)
	      opA2kdw=opA2kd; 
	    else
	      noalias(opA2kdw)=MijxBDeltaEij(opA2kd,(&AcutN[Nshell]),ibl1,ibl2,true,false,omegabar); // Nb1_kept x Nb2_disc

	    // HERE
	    noalias(Aux2dd)=prod(Prod2dk,opA2kdw);  // Nb2_disc x Nb2_disc


	    boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	      diagDD(Aux2dd,
		     boost::numeric::ublas::range(0,Aux2dd.size1()),
		     boost::numeric::ublas::range(0,Aux2dd.size2())); 

	    Termo2DDAKD=sum(diagDD);
  
	  }
	  // end if Nst_ketp_bl2 neq 0


	// COMPLEX operators
	} else {

	  if (printstuff)
	    cout << " ComplexOperator. Termo2ADD and Termo2AKD" 
		 << endl;
	  // if printstuff


	  boost::numeric::ublas::matrix<complex<double> > cRho2dd(Nst_bl2_disc,Nst_bl2_disc);

	  boost::numeric::ublas::matrix<complex<double> > copA2ddw(Nst_bl1_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > copA2dd(Nst_bl1_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > copB2dd(Nst_bl2_disc,Nst_bl1_disc);
	  // aux
	  boost::numeric::ublas::matrix<complex<double> > cAux2dd(Nst_bl2_disc,Nst_bl2_disc);
	  boost::numeric::ublas::matrix<complex<double> > cProd2dd(Nst_bl2_disc,Nst_bl1_disc);


	  cRho2dd=RhoN[Nshell].cMatBlock2BLAS(ibl2,ibl2,false,false);
	  // Nst_bl2_disc x Nst_bl2_disc matrix (assuming real!)

	  copA2dd=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,false);

	  if (NonDiagGF)
	    copB2dd=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2,false,false));
	  // Nst_bl2_disc x Nst_bl1_disc matrix (assuming real!)
	  else
	    copB2dd=trans(copA2dd);

	  noalias(cProd2dd)=prod(cRho2dd,copB2dd); // Nb2_disc x Nb1_disc

	  if (CalcNorm)
	    copA2ddw=copA2dd; 
	  else
	    noalias(copA2ddw)=cMijxBDeltaEij(copA2dd,(&AcutN[Nshell]),ibl1,ibl2,false,false,omegabar); // Nb1_disc x Nb2_disc
	  
	  noalias(cAux2dd)=prod(cProd2dd,copA2ddw);  // Nb2_disc x Nb2_disc

	  boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	    cdiagDD(cAux2dd,
		    boost::numeric::ublas::range(0,cAux2dd.size1()),
		    boost::numeric::ublas::range(0,cAux2dd.size2())); 

	  
	  Termo2DDADD=sum(cdiagDD).real();
	  
	  
	  if (Nst_bl1_kept!=0){
	    // Added condition if Nst_bl2_kept!=0 !
	    boost::numeric::ublas::matrix<complex<double> > copA2kdw(Nst_bl1_kept,Nst_bl2_disc);
	    boost::numeric::ublas::matrix<complex<double> > copA2kd(Nst_bl1_kept,Nst_bl2_disc);
	    boost::numeric::ublas::matrix<complex<double> > copB2dk(Nst_bl2_disc,Nst_bl1_kept);
	    // aux
	    boost::numeric::ublas::matrix<complex<double> > cProd2dk(Nst_bl2_disc,Nst_bl1_kept);
	    boost::numeric::ublas::matrix<complex<double> > cAux2dd(Nst_bl2_disc,Nst_bl2_disc);

	    copA2kd=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2,true,false);

	    if (NonDiagGF)
	      copB2dk=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2,true,false));
	    // Nst_bl2_kept x Nst_bl1_disc matrix (assuming real!)
	    else
	      copB2dk=trans(copA2kd);
	  

	    // Positive and negative energies matrix elements
	    // Note for the future: sign is "-" if bosonic operators

	    noalias(cProd2dk)=prod(cRho2dd,copB2dk); // Nb2_disc x Nb1_kept

	    if (CalcNorm)
	      copA2kdw=copA2kd; 
	    else
	      noalias(copA2kdw)=cMijxBDeltaEij(copA2kd,(&AcutN[Nshell]),ibl1,ibl2,true,false,omegabar); // Nb1_kept x Nb2_disc

	    // HERE
	    noalias(cAux2dd)=prod(cProd2dk,copA2kdw);  // Nb2_disc x Nb2_disc


	    boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	      cdiagDD(cAux2dd,
		     boost::numeric::ublas::range(0,cAux2dd.size1()),
		     boost::numeric::ublas::range(0,cAux2dd.size2())); 

	    Termo2DDAKD=sum(cdiagDD).real();
  
	  }
	  // end if Nst_ketp_bl2 neq 0


	  if (printstuff)
	    cout << " ComplexOperator: Termo2DDADD= " << Termo2DDADD
		 << " Termo2DDAKD= " << Termo2DDAKD
		 << endl;
	  // if printstuff
	}
	// end if Op1 is REAL else Comples
      }
      // end calculate termo1DD



      // SU(2) symmetry: <ibl1|A|ibl2> (or <ibl2|B|ibl1> )
      // Given ibl1, ibl2, calculate the following scalars
      // Get Sibl1,Sibl2
      // 
      // OpA: auxCG1[0]= sum_Szilb1 CGordan(Sibl1, Szibl1, 0.5, sigma, Sibl2, Szibl2=sigma+Szibl1)

      // OpB: auxCG1[1]= sum_Szilb1 CGordan(Sibl1, Szibl1, 0.5, sigma, Sibl2, Szibl2=sigma+Szibl1)
      //  They are the same!!
      double CGfactor=1.0;
      // Get S from ibl
      double Sbl1=0.0;
      double Sbl2=0.0;
      if (AcutN[Nshell].totalS){
	CGfactor=0.0;
	Sbl1=AcutN[Nshell].GetQNumber(ibl1,AcutN[Nshell].Sqnumbers[0]); 
	Sbl2=AcutN[Nshell].GetQNumber(ibl2,AcutN[Nshell].Sqnumbers[0]); 
	// only a single SU(2) for now
	for (double Szbl1=-Sbl1;Szbl1<=Sbl1;Szbl1+=1.0){
	  double dSigma=0.5; 
	  double Szbl2=Szbl1+dSigma;
	  double auxCG=CGordan(Sbl1,Szbl1, 0.5, dSigma, Sbl2, Szbl2);
	  CGfactor+=auxCG*auxCG;
	}
	if (dEqual(CGfactor,0.0)){cout << "Ops. CGfactor = 0.0 " << endl;}
      }
      // end if totalS

      Positive_w+=CGfactor*(Termo1KK+Termo1DDADD+Termo2DDADD+Termo2DDAKD);

      Negative_w+=CGfactor*(Termo2KK+Termo1DDADD+Termo1DDADK+Termo2DDADD);

    }
    // end if there are kept states in blocks

  }
  // end loop in Op1N blocks

      
  // STOPPED HERE! Need to continue...


  if (printstuff)
    cout << " Positive_w = " << Positive_w
	 << " Negative_w = " << Negative_w
	 << endl;

  double BosonSign=1.0;
  // if (IsBosonic) BosonSign=-1.0;
  if (CalcNorm){rho_wM=Positive_w+(BosonSign)*Negative_w;}else{
    if (omegabar>=0.0){rho_wM=Positive_w;}else{rho_wM=Negative_w;}
  } 
  // end if CalcNorm

  return(rho_wM);

}

//////////////////
///////////////////



