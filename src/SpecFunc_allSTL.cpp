

#include <cmath>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "SpecFuncClass.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)


#include "MyStructures.h"



using namespace std;




/// REAL USING STL VECTORS

void CSpecFunction::MijxBDeltaEij(vector<double>& Mout, CNRGmatrix* pMat, CNRGarray* pAeig, int ibl1, int ibl2, double omega){

  vector<double> Maux;
  // Neet so
  pMat->MatBlock2STL(Maux,ibl1,ibl2);

  MijxBDeltaEij(Mout, Maux, pAeig, ibl1, ibl2, omega);

}
///

void CSpecFunction::MijxBDeltaEij(vector<double>& Mout, vector<double>& Min, CNRGarray* pAeig, int ibl1, int ibl2, double omega){


  int Nst1=pAeig->GetBlockSize(ibl1);
  int Nst2=pAeig->GetBlockSize(ibl2);


  // Check sizea
  if ( (Nst1<=0)||(Nst2<=0) ){
    cout << " MijxBDeltaEij: Missing blocks. Cannot calculate Mij_kpi,kpi*BD(Ei-Ej) " << endl;
    return;
  }

  // Check sizes
  if ( (Min.size()!=Nst1*Nst2) ){
    cout << " MijxBDeltaEij: Matrix size mismatch. Cannot calculate Mij_kpi,kpi*BD(Ei-Ej) " << endl;
    return;
  }

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


  int ist=0;
  int iout=0;
  for (int ii=0;ii<Nst1;ii++){
    Ei=pAeig->dEn[ist0+ii];
    if (print) cout << " Ei = " << Ei;
    int jst=0;
    for (int jj=0;jj<Nst2;jj++){
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
	double PreFac=1.0-exp(log(fabs((Ej-Ei)/dBroadTemp*TempDN)));
	// 	    if (UseCFS==2) BD=PreFac*BDeltaTemp(omega,Ej-Ei,dBroadTemp*TempDN); // FDM-NRG
	// 	    else BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp); // CFS Anders
	BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp); // CFS Anders
	//  	cout << " |Ej-Ei|<TempDN: |Ej-Ei|=" << fabs(Ej-Ei) 
	//  	     << " TempDN = " << TempDN << endl;
      }
      // end if |Ej-Ei|<TempDN
      if (print) cout << " BD =  " << BD;
      Mout.push_back(Min[iout]*BD); // Uncomment this to work
      iout++;
      jst++;
      // END TEMP
    }
    ist++;
    if (print) cout << endl;
  }
  // end loop in i,j

}
//


// Using keep/disc states


void CSpecFunction::MijxBDeltaEij(vector<double>& Mout, CNRGmatrix* pMat, CNRGarray* pAeig, int ibl1, int ibl2, bool kp1, bool kp2, double omega){

  vector<double> Maux;
  // Neet so
  pMat->MatBlock2STL(Maux,ibl1,ibl2,kp1,kp2);

  MijxBDeltaEij(Mout, Maux, pAeig, ibl1, ibl2, kp1, kp2,  omega);

}
///

void CSpecFunction::MijxBDeltaEij(vector<double>& Mout, vector<double>& Min, CNRGarray* pAeig, int ibl1, int ibl2, bool kp1, bool kp2, double omega){


  int Nst1=pAeig->GetBlockSize(ibl1);
  int Nst2=pAeig->GetBlockSize(ibl2);


  int Nst1_kp1=pAeig->GetBlockSize(ibl1,kp1);
  int Nst2_kp2=pAeig->GetBlockSize(ibl2,kp2);

  // Check sizea
  if ( (Nst1_kp1<=0)||(Nst2_kp2<=0) ){
    cout << " MijxBDeltaEij: Missing blocks. Cannot calculate Mij_kpi,kpi*BD(Ei-Ej) " << endl;
    return;
  }

  // Check sizes
  if ( (Min.size()!=Nst1_kp1*Nst2_kp2) ){
    cout << " MijxBDeltaEij: Matrix size mismatch. Cannot calculate Mij_kpi,kpi*BD(Ei-Ej) " << endl;
    return;
  }

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
  int iout=0;
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
	    double PreFac=1.0-exp(log(fabs((Ej-Ei)/dBroadTemp*TempDN)));
// 	    if (UseCFS==2) BD=PreFac*BDeltaTemp(omega,Ej-Ei,dBroadTemp*TempDN); // FDM-NRG
// 	    else BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp); // CFS Anders
	    BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp); // CFS Anders
	    //  	cout << " |Ej-Ei|<TempDN: |Ej-Ei|=" << fabs(Ej-Ei) 
	    //  	     << " TempDN = " << TempDN << endl;
	  }
	  // end if |Ej-Ei|<TempDN
	  if (print) cout << " BD =  " << BD;
	  Mout.push_back(Min[iout]*BD); // Uncomment this to work
	  iout++;
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

}
//
/////
//////////////////
double CSpecFunction::DMNRG_SpecDens_M_STL(double omegabar,int Nshell,
					   bool CalcNorm){

  // Need to get the Connecting matrix elements and energies
  // to and from a given block (say, blGS)
  // Interesting procedure...

  // Loop in blocks of Operator


  double rho_wM=0.0;

  if (Nshell>NshellMax-1){
    cout << "DMNRG_SpecDens_M_STL: Nshell > NshellMax . Return 0" << endl;
    return(0.0);
  }

  if (Nshell==NshellMax-1){
    if (dEqual(Temp,0.0))
      Betabar=1.0e10;
    rho_wM=CalcSpecCosti_T_N(Nshell,Betabar,omegabar,CalcNorm);
    return(rho_wM);
  }
  //

  double Positive_w=0.0;
  double Negative_w=0.0;


  // cblas_dgemm variables
  double ALPHA=1.0, BETA=0.0;
//   C <- ALPHA*A.B + BETA*C

  complex<double> cALPHA=(1.0,0.0), cBETA=(0.0,0.0);

  bool printstuff=false;

  // Need to Check Sync!!
  bool ChkSyncRhoAcut=RhoN[Nshell].ChkSync((&AcutN[Nshell]));
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));

  if ( (!ChkSyncRhoAcut)||(!ChkSyncOp1Acut) ){
    cout << "DMNRG_SpecDens_M_STL: RhoN or Op1N not in sync with AcutN " << endl;
    return(0.0);
  }

  if (printstuff)
    cout << " DMNRG_SpecDens_M_STL: Nshell = " << Nshell
	 << " omegabar = " << omegabar
	 << " CalcNorm = " << CalcNorm
	 << endl;


  for (int iMatbl=0;iMatbl<Op1N[Nshell].NumMatBlocks();iMatbl++){
    int ibl1=0;
    int ibl2=0;

    Op1N[Nshell].GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);


    // Ok, remember the input Op2N is the COMPLEX CONJUGATE of what we want!
    //  
    //  So (B)_{ibl2 ibl1} = ((B+)_{ibl1 ibl2})* <- THIS is Op2

    //if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl2,ibl1)<0) ){
    if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl1,ibl2)<0) ){
      cout << "DMNRG_SpecDens_M_STL: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return(0.0);
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 
 
    // Check block sizes
    if ( (Nst_bl1!=RhoN[Nshell].GetBlockSize(ibl1))||
	 (Nst_bl2!=RhoN[Nshell].GetBlockSize(ibl2)) ){
      cout << "DMNRG_SpecDens_M_STL: Block in Op1 not the same size as block in RhoN " << endl;
      return(0.0);
    }
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "DMNRG_SpecDens_M_STL: Block in Op1 not the same size as block in Op2 " << endl;
      return(0.0);
    }


    if (printstuff){
      cout << " ibl1 = " << ibl1
	   << " ibl2 = " << ibl2
	   << endl;

      cout << " Nst_bl1= " << Nst_bl1 
	   << " Nst_bl2= " << Nst_bl2 
	   << endl;

    }
    // end printstuff


    // Termos
    double Termo1=0.0;
    double Termo2=0.0;

    // REAL operators
    if ( ( !(Op1N[Nshell].IsComplex) )&&( !(RhoN[Nshell].IsComplex) ) ) {

      vector<double> Rho1,opA1,opB1t,opA1w;

      RhoN[Nshell].MatBlock2STL(Rho1,ibl1,ibl1);
      Op1N[Nshell].MatBlock2STL(opA1,ibl1,ibl2);
      // Nst_bl1 x Nst_bl2 matrix (assuming real!)


      if (NonDiagGF)
	Op2N[Nshell].MatBlock2STL(opB1t,ibl1,ibl2);
      // Nst_bl2 x Nst_bl1 matrix (assuming real!)
      else
	opB1t=opA1;

      // Positive and negative energies matrix elements
      // Note for the future: sign is "-" if bosonic operators
      //noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix

      vector<double> Prod1(Nst_bl2*Nst_bl1);
      vector<double> Aux2disc(Nst_bl2*Nst_bl2);

      // Here's the first test

      // (opB1t_{Nst_bl1 x Nst_bl2})^T. Rho1_{Nst_bl1 x Nst_bl1} = Prod1_{Nst_bl2 x Nst_bl1}

      if (printstuff) cout << " BLAS: Prod1 = opB1.Rho1 ..." << endl;

      cblas_dgemm(CblasRowMajor,  CblasTrans, CblasNoTrans, 
		  Nst_bl2, Nst_bl1, Nst_bl1,
		  ALPHA,
		  &opB1t[0],Nst_bl2,
		  &Rho1[0], Nst_bl1, 
		  BETA, 
		  &Prod1[0], Nst_bl1);

      //(A(n1 x n0))^T . B_(n1x n2) = C_(n0 x n2)

      //M=n0 -> No of Rows in C (=rows in A^T): n0
      //N=n2 -> No of Cols in B and C: n2
      //K=n1 -> No of rows in ORIGINAL A (= no rows in B).

      //lda A = n0 no cols in ORIGINAL A.
      //lda B = n2 no cols in B.
      //lda C = n2 no cols in C (or no. of rows in C^T)

      if (CalcNorm) opA1w=opA1; else
	MijxBDeltaEij(opA1w,opA1,(&AcutN[Nshell]),
		      ibl1,ibl2,omegabar); // Nbl1 x Nbl2
	  
	  
      if (printstuff) cout 
			<< " ... done. Now Aux2 = Prod1.opA1w opA1w.size() = " 
			<< opA1w.size()<< endl;

      cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
		  Nst_bl2, Nst_bl2, Nst_bl1,
		  ALPHA,
		  &Prod1[0],Nst_bl1,
		  &opA1w[0], Nst_bl2, 
		  BETA, 
		  &Aux2disc[0], Nst_bl2);

      // A = Prod1_{Nst_bl2 x Nst_bl1}
      // B = opA1w_{Nst_bl1 x Nst_bl2}
      // C = Aux2disc_{Nst_bl2 x Nst_bl2}
      //
      // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

      // [DEFAULT]
      // M=n0 -> No of Rows in A and C: 
      // N=n2 -> No of Cols in B and C: 
      // K=n1 -> No of Cols in A and rows in B): 
	  
      // lda A = n1 no. cols in A.
      // lda B = n2 no. cols in B.
      // lda C = n2 no cols in C. 


      // Taking the trace
      int il=0;
      for (int ii=0;ii<Nst_bl2;ii++){
	il=ij2r(Nst_bl2,ii,ii);
	Termo1+=Aux2disc[il];
      }

      if (printstuff){ 
	cout 
	  << " ... done. Termo1 = " 
	  << Termo1 << endl;
	cout << " (opB1)^T : " << endl;
	il=0; 
	for (int ii=0;ii<Nst_bl1;ii++){
	  for (int jj=0;jj<Nst_bl2;jj++){
	    cout << opB1t[il] << " ";
	    il++;
	  }
	  cout << endl;
	}
	cout << " Rho1 : " << endl;
	il=0;
	for (int ii=0;ii<Nst_bl1;ii++){
	  for (int jj=0;jj<Nst_bl1;jj++){
	    cout << Rho1[il] << " ";
	    il++;
	  }
	  cout << endl;
	}
	cout << " Prod1=opB1.Rho1 : " << endl;
	il=0;
	for (int ii=0;ii<Nst_bl2;ii++){
	  for (int jj=0;jj<Nst_bl1;jj++){
	    cout << Prod1[il] << " ";
	    il++;
	  }
	  cout << endl;
	}
	cout << " opA1w : " << endl;
	il=0;
	for (int ii=0;ii<Nst_bl1;ii++){
	  for (int jj=0;jj<Nst_bl2;jj++){
	    cout << opA1w[il] << " ";
	    il++;
	  }
	  cout << endl;
	}
	cout << "Aux2 = Prod1.opA1w : " << endl;
	il=0;
	for (int ii=0;ii<Nst_bl2;ii++){
	  for (int jj=0;jj<Nst_bl2;jj++){
	    cout << Aux2disc[il] << " ";
	    il++;
	  }
	  cout << endl;
	}
      }
      // end printstuff

      // Termo2
      
      vector<double> Rho2,opA2,opB2t,opA2w;

      vector<double> Prod2(Nst_bl2*Nst_bl1);
      vector<double> Aux2(Nst_bl2*Nst_bl2);

      RhoN[Nshell].MatBlock2STL(Rho2,ibl2,ibl2);
      // Nst_bl2 x Nst_bl2 matrix (assuming real!)

      Op1N[Nshell].MatBlock2STL(opA2,ibl1,ibl2);
      // Nst_bl1 x Nst_bl2 matrix (assuming real!)

      if (NonDiagGF)
	Op2N[Nshell].MatBlock2STL(opB2t,ibl1,ibl2);
      // Nst_bl2 x Nst_bl1 matrix (assuming real!)
      else
	opB2t=opA2;

      cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasTrans, 
		  Nst_bl2, Nst_bl1, Nst_bl2,
		  ALPHA,
		  &Rho2[0], Nst_bl2, 
		  &opB2t[0],Nst_bl2,
		  BETA, 
		  &Prod2[0], Nst_bl1);
	  
      // A = Rho2 (Nst_bl2 x Nst_bl2)
      // B = opB2t(Nst_bl1 x Nst_bl2)
      // C = Prod2(Nst_bl2 x Nst_bl1)

      // A(n0 x n1) . (B_(n2 x n1))^T = C_(n0 x n2)

      //M=n0 -> No of Rows in A and C 
      //N=n2 -> No of Cols in C (=no rows in ORIGINAL B).
      //K=n1 -> No of Cols in A (= no cols in ORIGINAL B).

      //lda A = n1 no cols in  A.
      //lda B = n1 no cols in ORIGINAL B.
      //lda C = n2 no cols in C

	  

      // Calculate A_ij * Delta(w-(Ej-Ei))

      if (CalcNorm) opA2w=opA2; 
      else
	MijxBDeltaEij(opA2w,opA2,(&AcutN[Nshell]),
		      ibl1,ibl2,omegabar); 
      // Nbl1 x Nbl2

      cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
		  Nst_bl2, Nst_bl2, Nst_bl1,
		  ALPHA,
		  &Prod2[0],Nst_bl1,
		  &opA2w[0], Nst_bl2, 
		  BETA, 
		  &Aux2[0], Nst_bl2);

      // A = Prod2(Nst_bl2 x Nst_bl1)
      // B = opA2w_{Nst_bl1 x Nst_bl2}
      // C = Aux2kk_{Nst_bl2 x Nst_bl2}
      //
      // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)
      
      // [DEFAULT]
      // M=n0 -> No of Rows in A and C: 
      // N=n2 -> No of Cols in B and C: 
      // K=n1 -> No of Cols in A and rows in B): 
	  
      // lda A = n1 no. cols in A.
      // lda B = n2 no. cols in B.
      // lda C = n2 no cols in C. 


      // Taking the trace
      il=0;
      for (int ii=0;ii<Nst_bl2;ii++){
	il=ij2r(Nst_bl2,ii,ii);
	Termo2+=Aux2[il];
      }

      
      // COMPLEX operators
    } else {

      if (printstuff)
	cout << " ComplexOperator. Termo1 " << endl;
      // if printstuff

      boost::numeric::ublas::matrix<complex<double> > cRho1(Nst_bl1,Nst_bl1);

      boost::numeric::ublas::matrix<complex<double> > copA1w(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > copA1(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > copB1(Nst_bl2,Nst_bl1);
      // aux
      boost::numeric::ublas::matrix<complex<double> > cAux2disc(Nst_bl2,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > cProd1(Nst_bl2,Nst_bl1);


      cRho1=RhoN[Nshell].cMatBlock2BLAS(ibl1,ibl1);

      copA1=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2);
      // Nst_bl1 x Nst_bl2 matrix (assuming real!)

      if (NonDiagGF)
	//opB=Op2N[Nshell].MatBlock2BLAS(ibl2,ibl1);
	copB1=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2)); // Nst_bl2 x Nst_bl1 matrix (assuming real!)
      else
	copB1=trans(copA1);

      // Positive and negative energies matrix elements
      // Note for the future: sign is "-" if bosonic operators

      noalias(cProd1)=prod(copB1,cRho1); // Nb2_disc x Nb1_kept

      if (CalcNorm) copA1w=copA1; else
	noalias(copA1w)=cMijxBDeltaEij(copA1,(&AcutN[Nshell]),ibl1,ibl2,omegabar); // Nbl1 x Nbl2

      noalias(cAux2disc)=prod(cProd1,copA1w);  // Nbl2 x Nbl2

      boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	cdiag1(cAux2disc,
	       boost::numeric::ublas::range(0,cAux2disc.size1()),
	       boost::numeric::ublas::range(0,cAux2disc.size2())); 

      Termo1=sum(cdiag1).real();


      if (printstuff)
	cout << " ComplexOperator. Termo1= " << Termo1 
	     << endl;
      // if printstuff


      // Termo2
      
      boost::numeric::ublas::matrix<complex<double> > cRho2(Nst_bl2,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > copA2(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > copA2w(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > copB2(Nst_bl2,Nst_bl1);

      // Can I use opB1 and opB2 here instead? (also: do I need two opA matrices??) Check...
      boost::numeric::ublas::matrix<complex<double> > cProd2(Nst_bl2,Nst_bl1);

      cRho2=RhoN[Nshell].cMatBlock2BLAS(ibl2,ibl2);

      copA2=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2);
      if (NonDiagGF)
	copB2=trans(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2)); 
      // Nst_bl2 x Nst_bl1 matrix
      else
	copB2=trans(copA2);

      noalias(cProd2)=prod(cRho2,copB2); // Nbl2 x Nbl1

      // Calculate A_ij * Delta(w-(Ej-Ei))

      if (CalcNorm) copA2w=copA2; else
	noalias(copA2w)=cMijxBDeltaEij(copA2,(&AcutN[Nshell]),ibl1,ibl2,omegabar); // Nbl1 x Nbl2


      noalias(cRho2)=prod(cProd2,copA2w);  // Nbl2 x Nbl2

      boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	cdiag2(cRho2,
	       boost::numeric::ublas::range(0,cRho2.size1()),
	       boost::numeric::ublas::range(0,cRho2.size2())); 

      Termo2=sum(cdiag2).real();

      
    }
    // end if Op1 is REAL else Complex


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

    double BosonSign=1.0;
    if (IsBosonic) BosonSign=-1.0; // To be implemented
    // No need for special attention to CalcNorm...
    // if (CalcNorm){rho_wM=Positive_w+(BosonSign)*Negative_w;}else{
    //   if (omegabar>=0.0){rho_wM=Positive_w;}else{rho_wM=Negative_w;}
    // } 
    // // end if CalcNorm

    
    rho_wM+=CGfactor*(Termo1+(BosonSign)*Termo2);

    if (printstuff){
      cout << " DMNRG_Spec_M_STL: trace = " << Termo1+Termo2
	   << " rho_wM acc. = " << rho_wM 
	   << endl;
    }

  }
  // end loop in Op1N blocks


  return(rho_wM);

}

//////////////////



//////


//////////////////
double CSpecFunction::CFS_SpecDens_M_STL(double omegabar,int Nshell,
					 bool CalcNorm){

  // Need to get the Connecting matrix elements and energies
  // to and from a given block (say, blGS)
  // Interesting procedure...

  // Loop in blocks of Operator


  double rho_wM=0.0;

  if (Nshell>NshellMax-1){
    cout << "CFS_SpecDens_M: Nshell > NshellMax . Return 0" << endl;
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


  // cblas_dgemm variables
  double ALPHA=1.0, BETA=0.0;
//   C <- ALPHA*A.B + BETA*C

  complex<double> cALPHA=(1.0,0.0), cBETA=(0.0,0.0);

  bool printstuff=false;

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
      cout << "CFS_SpecDens_M_STL: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return(0.0);
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 
 
    // Check block sizes
    if ( (Nst_bl1!=RhoN[Nshell].GetBlockSize(ibl1))||
	 (Nst_bl2!=RhoN[Nshell].GetBlockSize(ibl2)) ){
      cout << "CFS_SpecDens_M_STL: Block in Op1 not the same size as block in RhoN " << endl;
      return(0.0);
    }
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "CFS_SpecDens_M_STL: Block in Op1 not the same size as block in Op2 " << endl;
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
      cout << "CFS_SpecDens_M_STL: Error: Kept+Dics NOT equal to Nstblock! " << endl;
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

	  vector<double> Rho1,opA1,opB1t,opA1w;

	  RhoN[Nshell].MatBlock2STL(Rho1,ibl1,ibl1,true,true);

	  Op1N[Nshell].MatBlock2STL(opA1,ibl1,ibl2,true,false);
	  // Nst_bl1_kept x Nst_bl2_disc matrix (assuming real!)


 	  if (NonDiagGF)
 	    Op2N[Nshell].MatBlock2STL(opB1t,ibl1,ibl2,true,false);
	  // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
 	  else
 	    opB1t=opA1;

	  // Positive and negative energies matrix elements
	  // Note for the future: sign is "-" if bosonic operators
	  //noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix

	  vector<double> Prod1(Nst_bl2_disc*Nst_bl1_kept);
	  vector<double> Aux2disc(Nst_bl2_disc*Nst_bl2_disc);

	  // Here's the first test

	  // (opB1t_{Nst_bl1_kept x Nst_bl2_disc})^T. Rho1_{Nst_bl1_kept x Nst_bl1_kept} = Prod1_{Nst_bl2_disc x Nst_bl1_kept}

	  if (printstuff) cout << " BLAS: Prod1 = opB1.Rho1 ..." << endl;

	  cblas_dgemm(CblasRowMajor,  CblasTrans, CblasNoTrans, 
		      Nst_bl2_disc, Nst_bl1_kept, Nst_bl1_kept,
		      ALPHA,
		      &opB1t[0],Nst_bl2_disc,
		      &Rho1[0], Nst_bl1_kept, 
		      BETA, 
		      &Prod1[0], Nst_bl1_kept);

	  //(A(n1 x n0))^T . B_(n1x n2) = C_(n0 x n2)

	  //M=n0 -> No of Rows in C (=rows in A^T): n0
	  //N=n2 -> No of Cols in B and C: n2
	  //K=n1 -> No of rows in ORIGINAL A (= no rows in B).

	  //lda A = n0 no cols in ORIGINAL A.
	  //lda B = n2 no cols in B.
	  //lda C = n2 no cols in C (or no. of rows in C^T)

	  if (CalcNorm) opA1w=opA1; else
	    MijxBDeltaEij(opA1w,opA1,(&AcutN[Nshell]),
			  ibl1,ibl2,true,false,omegabar); // Nb1_kept x Nb2_disc
	  
	  
	  if (printstuff) cout 
	    << " ... done. Now Aux2 = Prod1.opA1w opA1w.size() = " 
	    << opA1w.size()<< endl;

	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
		      Nst_bl2_disc, Nst_bl2_disc, Nst_bl1_kept,
		      ALPHA,
		      &Prod1[0],Nst_bl1_kept,
		      &opA1w[0], Nst_bl2_disc, 
		      BETA, 
		      &Aux2disc[0], Nst_bl2_disc);

	  // A = Prod1_{Nst_bl2_disc x Nst_bl1_kept}
	  // B = opA1w_{Nst_bl1_kept x Nst_bl2_disc}
	  // C = Aux2disc_{Nst_bl2_disc x Nst_bl2_disc}
	  //
	  // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

	  // [DEFAULT]
	  // M=n0 -> No of Rows in A and C: 
	  // N=n2 -> No of Cols in B and C: 
	  // K=n1 -> No of Cols in A and rows in B): 
	  
	  // lda A = n1 no. cols in A.
	  // lda B = n2 no. cols in B.
	  // lda C = n2 no cols in C. 


	  // Taking the trace
	  int il=0;
	  for (int ii=0;ii<Nst_bl2_disc;ii++){
	    il=ij2r(Nst_bl2_disc,ii,ii);
	    Termo1KK+=Aux2disc[il];
	  }

	  if (printstuff){ 
	    cout 
	    << " ... done. Termo1KK = " 
	    << Termo1KK << endl;
	    cout << " (opB1)^T : " << endl;
	    il=0; 
	    for (int ii=0;ii<Nst_bl1_kept;ii++){
	      for (int jj=0;jj<Nst_bl2_disc;jj++){
		cout << opB1t[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	    cout << " Rho1 : " << endl;
	    il=0;
	    for (int ii=0;ii<Nst_bl1_kept;ii++){
	      for (int jj=0;jj<Nst_bl1_kept;jj++){
		cout << Rho1[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	    cout << " Prod1=opB1.Rho1 : " << endl;
	    il=0;
	    for (int ii=0;ii<Nst_bl2_disc;ii++){
	      for (int jj=0;jj<Nst_bl1_kept;jj++){
		cout << Prod1[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	    cout << " opA1w : " << endl;
	    il=0;
	    for (int ii=0;ii<Nst_bl1_kept;ii++){
	      for (int jj=0;jj<Nst_bl2_disc;jj++){
		cout << opA1w[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	    cout << "Aux2 = Prod1.opA1w : " << endl;
	    il=0;
	    for (int ii=0;ii<Nst_bl2_disc;ii++){
	      for (int jj=0;jj<Nst_bl2_disc;jj++){
		cout << Aux2disc[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	  }
	  // end printstuff

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


      // Termo2KK 
      double Termo2KK=0.0;

      if ( (Nst_bl1_disc!=0)&&(Nst_bl2_kept!=0)&&((omegabar<0.0)||(CalcNorm)) ){

	if (printstuff)
	  cout << "  Calculating Termo2KK " << endl;
	// if printstuff


	// REAL operators
	if ( !(Op1N[Nshell].IsComplex) ){

	  vector<double> Rho2,opA2,opB2t,opA2w;

	  vector<double> Prod2(Nst_bl2_kept*Nst_bl1_disc);
	  vector<double> Aux2kk(Nst_bl2_kept*Nst_bl2_kept);

	  RhoN[Nshell].MatBlock2STL(Rho2,ibl2,ibl2,true,true);
	  // Nst_bl2_kept x Nst_bl2_kept matrix (assuming real!)

	  Op1N[Nshell].MatBlock2STL(opA2,ibl1,ibl2,false,true);
	  // Nst_bl1_disc x Nst_bl2_kept matrix (assuming real!)

 	  if (NonDiagGF)
 	    Op2N[Nshell].MatBlock2STL(opB2t,ibl1,ibl2,false,true);
	  // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
 	  else
 	    opB2t=opA2;

	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasTrans, 
		      Nst_bl2_kept, Nst_bl1_disc, Nst_bl2_kept,
		      ALPHA,
		      &Rho2[0], Nst_bl2_kept, 
		      &opB2t[0],Nst_bl2_kept,
		      BETA, 
		      &Prod2[0], Nst_bl1_disc);
	  
	  // A = Rho2 (Nst_bl2_kept x Nst_bl2_kept)
	  // B = opB2t(Nst_bl1_disc x Nst_bl2_kept)
	  // C = Prod2(Nst_bl2_kept x Nst_bl1_disc)

	  // A(n0 x n1) . (B_(n2 x n1))^T = C_(n0 x n2)

	  //M=n0 -> No of Rows in A and C 
	  //N=n2 -> No of Cols in C (=no rows in ORIGINAL B).
	  //K=n1 -> No of Cols in A (= no cols in ORIGINAL B).

	  //lda A = n1 no cols in  A.
	  //lda B = n1 no cols in ORIGINAL B.
	  //lda C = n2 no cols in C

	  

	  // Calculate A_ij * Delta(w-(Ej-Ei))

	  if (CalcNorm) opA2w=opA2; 
	  else
	    MijxBDeltaEij(opA2w,opA2,(&AcutN[Nshell]),
			  ibl1,ibl2,false,true,omegabar); 
	  // Nb1_disc x Nb2_kept

	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
		      Nst_bl2_kept, Nst_bl2_kept, Nst_bl1_disc,
		      ALPHA,
		      &Prod2[0],Nst_bl1_disc,
		      &opA2w[0], Nst_bl2_kept, 
		      BETA, 
		      &Aux2kk[0], Nst_bl2_kept);

	  // A = Prod2(Nst_bl2_kept x Nst_bl1_disc)
	  // B = opA2w_{Nst_bl1_disc x Nst_bl2_kept}
	  // C = Aux2kk_{Nst_bl2_kept x Nst_bl2_kept}
	  //
	  // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

	  // [DEFAULT]
	  // M=n0 -> No of Rows in A and C: 
	  // N=n2 -> No of Cols in B and C: 
	  // K=n1 -> No of Cols in A and rows in B): 
	  
	  // lda A = n1 no. cols in A.
	  // lda B = n2 no. cols in B.
	  // lda C = n2 no cols in C. 


	  // Taking the trace
	  int il=0;
	  for (int ii=0;ii<Nst_bl2_kept;ii++){
	    il=ij2r(Nst_bl2_kept,ii,ii);
	    Termo2KK+=Aux2kk[il];
	  }
	  
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

      Positive_w+=CGfactor*Termo1KK;

      Negative_w+=CGfactor*Termo2KK;

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



// Full Density-Matrix calculations


///////////////////


//////////////////
double CSpecFunction::FDM_SpecDens_M_STL(double omegabar,int Nshell,
					 bool CalcNorm){

  // Need to get the Connecting matrix elements and energies
  // to and from a given block (say, blGS)
  // Interesting procedure...

  // Loop in blocks of Operator


  double rho_wM=0.0;

  if (Nshell>NshellMax-1){
    cout << "FDM_SpecDens_M_STL: Nshell > NshellMax . Return 0" << endl;
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


  // cblas_dgemm variables
  double ALPHA=1.0, BETA=0.0;
//   C <- ALPHA*A.B + BETA*C

  complex<double> cALPHA=(1.0,0.0), cBETA=(0.0,0.0);

  bool printstuff=false;

  // Need to Check Sync!!
  bool ChkSyncRhoAcut=RhoN[Nshell].ChkSync((&AcutN[Nshell]));
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));

  if ( (!ChkSyncRhoAcut)||(!ChkSyncOp1Acut) ){
    cout << "FDM_SpecDens_M_STL: RhoN or Op1N not in sync with AcutN " << endl;
    return(0.0);
  }

  if (printstuff)
    cout << " FDM_SpecDens_M_STL: Nshell = " << Nshell
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
      cout << "FDM_SpecDens_M_STL: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return(0.0);
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 
 
    // Check block sizes
    if ( (Nst_bl1!=RhoN[Nshell].GetBlockSize(ibl1))||
	 (Nst_bl2!=RhoN[Nshell].GetBlockSize(ibl2)) ){
      cout << "FDM_SpecDens_M_STL: Block in Op1 not the same size as block in RhoN " << endl;
      return(0.0);
    }
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "FDM_SpecDens_M_STL: Block in Op1 not the same size as block in Op2 " << endl;
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
      cout << "FDM_SpecDens_M_STL: Error: Kept+Dics NOT equal to Nstblock! " << endl;
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

	  vector<double> Rho1,opA1,opB1t,opA1w;

	  RhoN[Nshell].MatBlock2STL(Rho1,ibl1,ibl1,true,true);

	  Op1N[Nshell].MatBlock2STL(opA1,ibl1,ibl2,true,false);
	  // Nst_bl1_kept x Nst_bl2_disc matrix (assuming real!)


 	  if (NonDiagGF)
 	    Op2N[Nshell].MatBlock2STL(opB1t,ibl1,ibl2,true,false);
	  // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
 	  else
 	    opB1t=opA1;

	  // Positive and negative energies matrix elements
	  // Note for the future: sign is "-" if bosonic operators
	  //noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix

	  vector<double> Prod1(Nst_bl2_disc*Nst_bl1_kept);
	  vector<double> Aux2disc(Nst_bl2_disc*Nst_bl2_disc);

	  // Here's the first test

	  // (opB1t_{Nst_bl1_kept x Nst_bl2_disc})^T. Rho1_{Nst_bl1_kept x Nst_bl1_kept} = Prod1_{Nst_bl2_disc x Nst_bl1_kept}

	  if (printstuff) cout << " BLAS: Prod1 = opB1.Rho1 ..." << endl;

	  cblas_dgemm(CblasRowMajor,  CblasTrans, CblasNoTrans, 
		      Nst_bl2_disc, Nst_bl1_kept, Nst_bl1_kept,
		      ALPHA,
		      &opB1t[0],Nst_bl2_disc,
		      &Rho1[0], Nst_bl1_kept, 
		      BETA, 
		      &Prod1[0], Nst_bl1_kept);

	  //(A(n1 x n0))^T . B_(n1x n2) = C_(n0 x n2)

	  //M=n0 -> No of Rows in C (=rows in A^T): n0
	  //N=n2 -> No of Cols in B and C: n2
	  //K=n1 -> No of rows in ORIGINAL A (= no rows in B).

	  //lda A = n0 no cols in ORIGINAL A.
	  //lda B = n2 no cols in B.
	  //lda C = n2 no cols in C (or no. of rows in C^T)

	  if (CalcNorm) opA1w=opA1; else
	    MijxBDeltaEij(opA1w,opA1,(&AcutN[Nshell]),
			  ibl1,ibl2,true,false,omegabar); // Nb1_kept x Nb2_disc
	  
	  
	  if (printstuff) cout 
	    << " ... done. Now Aux2 = Prod1.opA1w opA1w.size() = " 
	    << opA1w.size()<< endl;

	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
		      Nst_bl2_disc, Nst_bl2_disc, Nst_bl1_kept,
		      ALPHA,
		      &Prod1[0],Nst_bl1_kept,
		      &opA1w[0], Nst_bl2_disc, 
		      BETA, 
		      &Aux2disc[0], Nst_bl2_disc);

	  // A = Prod1_{Nst_bl2_disc x Nst_bl1_kept}
	  // B = opA1w_{Nst_bl1_kept x Nst_bl2_disc}
	  // C = Aux2disc_{Nst_bl2_disc x Nst_bl2_disc}
	  //
	  // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

	  // [DEFAULT]
	  // M=n0 -> No of Rows in A and C: 
	  // N=n2 -> No of Cols in B and C: 
	  // K=n1 -> No of Cols in A and rows in B): 
	  
	  // lda A = n1 no. cols in A.
	  // lda B = n2 no. cols in B.
	  // lda C = n2 no cols in C. 


	  // Taking the trace
	  int il=0;
	  for (int ii=0;ii<Nst_bl2_disc;ii++){
	    il=ij2r(Nst_bl2_disc,ii,ii);
	    Termo1KK+=Aux2disc[il];
	  }

	  if (printstuff){ 
	    cout 
	    << " ... done. Termo1KK = " 
	    << Termo1KK << endl;
	    cout << " (opB1)^T : " << endl;
	    il=0; 
	    for (int ii=0;ii<Nst_bl1_kept;ii++){
	      for (int jj=0;jj<Nst_bl2_disc;jj++){
		cout << opB1t[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	    cout << " Rho1 : " << endl;
	    il=0;
	    for (int ii=0;ii<Nst_bl1_kept;ii++){
	      for (int jj=0;jj<Nst_bl1_kept;jj++){
		cout << Rho1[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	    cout << " Prod1=opB1.Rho1 : " << endl;
	    il=0;
	    for (int ii=0;ii<Nst_bl2_disc;ii++){
	      for (int jj=0;jj<Nst_bl1_kept;jj++){
		cout << Prod1[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	    cout << " opA1w : " << endl;
	    il=0;
	    for (int ii=0;ii<Nst_bl1_kept;ii++){
	      for (int jj=0;jj<Nst_bl2_disc;jj++){
		cout << opA1w[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	    cout << "Aux2 = Prod1.opA1w : " << endl;
	    il=0;
	    for (int ii=0;ii<Nst_bl2_disc;ii++){
	      for (int jj=0;jj<Nst_bl2_disc;jj++){
		cout << Aux2disc[il] << " ";
		il++;
	      }
	      cout << endl;
	    }
	  }
	  // end printstuff

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

	  vector<double> Rho1dd,opA1dd,opB1ddt,opA1ddw;

	  vector<double> Prod1dd(Nst_bl2_disc*Nst_bl1_disc);
	  vector<double> Aux2dd(Nst_bl2_disc*Nst_bl2_disc);

	  RhoN[Nshell].MatBlock2STL(Rho1dd,ibl1,ibl1,false,false);
	  // Nst_bl1_disc x Nst_bl1_disc matrix (assuming real!)


	  Op1N[Nshell].MatBlock2STL(opA1dd,ibl1,ibl2,false,false);
	  // Nst_bl1_disc x Nst_bl2_disc matrix (assuming real!)


 	  if (NonDiagGF)
 	    Op2N[Nshell].MatBlock2STL(opB1ddt,ibl1,ibl2,false,false);
	  // Nst_bl2_disc x Nst_bl1_disc matrix (assuming real!)
 	  else
 	    opB1ddt=opA1dd;


	  cblas_dgemm(CblasRowMajor,  CblasTrans, CblasNoTrans, 
		      Nst_bl2_disc, Nst_bl1_disc, Nst_bl1_disc,
		      ALPHA,
		      &opB1ddt[0],Nst_bl2_disc,
		      &Rho1dd[0], Nst_bl1_disc, 
		      BETA, 
		      &Prod1dd[0], Nst_bl1_disc);
	  
	  // A = opB1ddt(Nst_bl1_disc x Nst_bl2_disc)
	  // B = Rho1dd (Nst_bl1_disc x Nst_bl1_disc)
	  // C = Prod1dd(Nst_bl2_disc x Nst_bl1_disc)

	  //(A(n1 x n0))^T . B_(n1x n2) = C_(n0 x n2)

	  //M=n0 -> No of Rows in C (=rows in A^T): n0
	  //N=n2 -> No of Cols in B and C: n2
	  //K=n1 -> No of rows in ORIGINAL A (= no rows in B).

	  //lda A = n0 no cols in ORIGINAL A.
	  //lda B = n2 no cols in B.
	  //lda C = n2 no cols in C (or no. of rows in C^T)


	  if (CalcNorm) opA1ddw=opA1dd; 
	  else
	    MijxBDeltaEij(opA1ddw,opA1dd,(&AcutN[Nshell]),
			  ibl1,ibl2,false,false,omegabar); 
	  // Nb1_disc x Nb2_disc

	  
	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
		      Nst_bl2_disc, Nst_bl2_disc, Nst_bl1_disc,
		      ALPHA,
		      &Prod1dd[0],Nst_bl1_disc,
		      &opA1ddw[0], Nst_bl2_disc, 
		      BETA, 
		      &Aux2dd[0], Nst_bl2_disc);

	  // A = Prod1dd_{Nst_bl2_disc x Nst_bl1_disc}
	  // B = opA1ddw_{Nst_bl1_disc x Nst_bl2_disc}
	  // C = Aux2dd_{Nst_bl2_disc x Nst_bl2_disc}
	  //
	  // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

	  // [DEFAULT]
	  // M=n0 -> No of Rows in A and C: 
	  // N=n2 -> No of Cols in B and C: 
	  // K=n1 -> No of Cols in A and rows in B): 
	  
	  // lda A = n1 no. cols in A.
	  // lda B = n2 no. cols in B.
	  // lda C = n2 no cols in C. 


	  // Taking the trace
	  int il=0;
	  for (int ii=0;ii<Nst_bl2_disc;ii++){
	    il=ij2r(Nst_bl2_disc,ii,ii);
	    Termo1DDADD+=Aux2dd[il];
	  }

	  
	  if (Nst_bl2_kept!=0){
	    // Added condition if Nst_bl2_kept!=0 !

	    vector<double> opA1dk,opB1kdt,opA1dkw;

	    vector<double> Prod1kd(Nst_bl2_kept*Nst_bl1_disc);
	    vector<double> Aux2kk(Nst_bl2_kept*Nst_bl2_kept);


	    Op1N[Nshell].MatBlock2STL(opA1dk,ibl1,ibl2,false,true);
	    // Nst_bl1_disc x Nst_bl2_kept matrix (assuming real!)


	    if (NonDiagGF)
	      Op2N[Nshell].MatBlock2STL(opB1kdt,ibl1,ibl2,true,false);
	    // Nst_bl2_kept x Nst_bl1_disc matrix (assuming real!)
	    else
	      opB1kdt=opA1dk;


	    cblas_dgemm(CblasRowMajor,  CblasTrans, CblasNoTrans, 
			Nst_bl2_kept, Nst_bl1_disc, Nst_bl1_disc,
			ALPHA,
			&opB1kdt[0],Nst_bl2_kept,
			&Rho1dd[0], Nst_bl1_disc, 
			BETA, 
			&Prod1kd[0], Nst_bl1_disc);
	    // A = opB1kdt(Nst_bl1_disc x Nst_bl2_kept)
	    // B = Rho1dd (Nst_bl1_disc x Nst_bl1_disc)
	    // C = Prod1kd(Nst_bl2_kept x Nst_bl1_disc)


	    //(A(n1 x n0))^T . B_(n1x n2) = C_(n0 x n2)

	    //M=n0 -> No of Rows in C (=rows in A^T): n0
	    //N=n2 -> No of Cols in B and C: n2
	    //K=n1 -> No of rows in ORIGINAL A (= no rows in B).

	    //lda A = n0 no cols in ORIGINAL A.
	    //lda B = n2 no cols in B.
	    //lda C = n2 no cols in C (or no. of rows in C^T)


	    if (CalcNorm) opA1dkw=opA1dk; 
	    else
	      MijxBDeltaEij(opA1dkw,opA1dk,(&AcutN[Nshell]),
			    ibl1,ibl2,false,true,omegabar); 
	    // Nb1_disc x Nb2_kept


	    // Positive and negative energies matrix elements
	    // Note for the future: sign is "-" if bosonic operators


	    cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
			Nst_bl2_kept, Nst_bl2_kept, Nst_bl1_disc,
			ALPHA,
			&Prod1kd[0],Nst_bl1_disc,
			&opA1dkw[0], Nst_bl2_kept, 
			BETA, 
			&Aux2kk[0], Nst_bl2_kept);

	    // A = Prod1kd_{Nst_bl2_kept x Nst_bl1_disc}
	    // B = opA1dkw_{Nst_bl1_disc x Nst_bl2_kept}
	    // C = Aux2kk_{Nst_bl2_kept x Nst_bl2_kept}
	    //
	    // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

	    // [DEFAULT]
	    // M=n0 -> No of Rows in A and C: 
	    // N=n2 -> No of Cols in B and C: 
	    // K=n1 -> No of Cols in A and rows in B): 
	  
	    // lda A = n1 no. cols in A.
	    // lda B = n2 no. cols in B.
	    // lda C = n2 no cols in C. 


	    // Taking the trace
	    int il=0;
	    for (int ii=0;ii<Nst_bl2_kept;ii++){
	      il=ij2r(Nst_bl2_kept,ii,ii);
	      Termo1DDADK+=Aux2kk[il];
	    }

  
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

	  vector<double> Rho2,opA2,opB2t,opA2w;

	  vector<double> Prod2(Nst_bl2_kept*Nst_bl1_disc);
	  vector<double> Aux2kk(Nst_bl2_kept*Nst_bl2_kept);

	  RhoN[Nshell].MatBlock2STL(Rho2,ibl2,ibl2,true,true);
	  // Nst_bl2_kept x Nst_bl2_kept matrix (assuming real!)

	  Op1N[Nshell].MatBlock2STL(opA2,ibl1,ibl2,false,true);
	  // Nst_bl1_disc x Nst_bl2_kept matrix (assuming real!)

 	  if (NonDiagGF)
 	    Op2N[Nshell].MatBlock2STL(opB2t,ibl1,ibl2,false,true);
	  // Nst_bl2_disc x Nst_bl1_kept matrix (assuming real!)
 	  else
 	    opB2t=opA2;

	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasTrans, 
		      Nst_bl2_kept, Nst_bl1_disc, Nst_bl2_kept,
		      ALPHA,
		      &Rho2[0], Nst_bl2_kept, 
		      &opB2t[0],Nst_bl2_kept,
		      BETA, 
		      &Prod2[0], Nst_bl1_disc);
	  
	  // A = Rho2 (Nst_bl2_kept x Nst_bl2_kept)
	  // B = opB2t(Nst_bl1_disc x Nst_bl2_kept)
	  // C = Prod2(Nst_bl2_kept x Nst_bl1_disc)

	  // A(n0 x n1) . (B_(n2 x n1))^T = C_(n0 x n2)

	  //M=n0 -> No of Rows in A and C 
	  //N=n2 -> No of Cols in C (=no rows in ORIGINAL B).
	  //K=n1 -> No of Cols in A (= no cols in ORIGINAL B).

	  //lda A = n1 no cols in  A.
	  //lda B = n1 no cols in ORIGINAL B.
	  //lda C = n2 no cols in C

	  

	  // Calculate A_ij * Delta(w-(Ej-Ei))

	  if (CalcNorm) opA2w=opA2; 
	  else
	    MijxBDeltaEij(opA2w,opA2,(&AcutN[Nshell]),
			  ibl1,ibl2,false,true,omegabar); 
	  // Nb1_disc x Nb2_kept

	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
		      Nst_bl2_kept, Nst_bl2_kept, Nst_bl1_disc,
		      ALPHA,
		      &Prod2[0],Nst_bl1_disc,
		      &opA2w[0], Nst_bl2_kept, 
		      BETA, 
		      &Aux2kk[0], Nst_bl2_kept);

	  // A = Prod2(Nst_bl2_kept x Nst_bl1_disc)
	  // B = opA2w_{Nst_bl1_disc x Nst_bl2_kept}
	  // C = Aux2kk_{Nst_bl2_kept x Nst_bl2_kept}
	  //
	  // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

	  // [DEFAULT]
	  // M=n0 -> No of Rows in A and C: 
	  // N=n2 -> No of Cols in B and C: 
	  // K=n1 -> No of Cols in A and rows in B): 
	  
	  // lda A = n1 no. cols in A.
	  // lda B = n2 no. cols in B.
	  // lda C = n2 no cols in C. 


	  // Taking the trace
	  int il=0;
	  for (int ii=0;ii<Nst_bl2_kept;ii++){
	    il=ij2r(Nst_bl2_kept,ii,ii);
	    Termo2KK+=Aux2kk[il];
	  }
	  
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
	  vector<double> Rho2dd,opA2dd,opB2ddt,opA2ddw;

	  vector<double> Prod2dd(Nst_bl2_disc*Nst_bl1_disc);
	  vector<double> Aux2dd(Nst_bl2_disc*Nst_bl2_disc);

	  RhoN[Nshell].MatBlock2STL(Rho2dd,ibl2,ibl2,false,false);
	  // Nst_bl2_disc x Nst_bl2_disc matrix (assuming real!)

	  Op1N[Nshell].MatBlock2STL(opA2dd,ibl1,ibl2,false,false);
	  // Nst_bl1_disc x Nst_bl2_disc matrix (assuming real!)

 	  if (NonDiagGF)
 	    Op2N[Nshell].MatBlock2STL(opB2ddt,ibl1,ibl2,false,false);
	  // Nst_bl2_disc x Nst_bl1_disc matrix (assuming real!)
 	  else
 	    opB2ddt=opA2dd;

	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasTrans, 
		      Nst_bl2_disc, Nst_bl1_disc, Nst_bl2_disc,
		      ALPHA,
		      &Rho2dd[0], Nst_bl2_disc, 
		      &opB2ddt[0],Nst_bl2_disc,
		      BETA, 
		      &Prod2dd[0], Nst_bl1_disc);
	  
	  // A = Rho2dd (Nst_bl2_disc x Nst_bl2_disc)
	  // B = opB2t(Nst_bl1_disc x Nst_bl2_disc)
	  // C = Prod2dd(Nst_bl2_disc x Nst_bl1_disc)

	  // A(n0 x n1) . (B_(n2 x n1))^T = C_(n0 x n2)

	  //M=n0 -> No of Rows in A and C 
	  //N=n2 -> No of Cols in C (=no rows in ORIGINAL B).
	  //K=n1 -> No of Cols in A (= no cols in ORIGINAL B).

	  //lda A = n1 no cols in  A.
	  //lda B = n1 no cols in ORIGINAL B.
	  //lda C = n2 no cols in C

	  if (CalcNorm) opA2ddw=opA2dd; 
	  else
	    MijxBDeltaEij(opA2ddw,opA2dd,(&AcutN[Nshell]),
			  ibl1,ibl2,false,false,omegabar); 
	  // Nb1_disc x Nb2_disc

	  cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
		      Nst_bl2_disc, Nst_bl2_disc, Nst_bl1_disc,
		      ALPHA,
		      &Prod2dd[0],Nst_bl1_disc,
		      &opA2ddw[0], Nst_bl2_disc, 
		      BETA, 
		      &Aux2dd[0], Nst_bl2_disc);

	  // A = Prod2dd_{Nst_bl2_disc x Nst_bl1_disc}
	  // B = opA2ddw_{Nst_bl1_disc x Nst_bl2_disc}
	  // C = Aux2dd_{Nst_bl2_disc x Nst_bl2_disc}
	  //
	  // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

	  // [DEFAULT]
	  // M=n0 -> No of Rows in A and C: 
	  // N=n2 -> No of Cols in B and C: 
	  // K=n1 -> No of Cols in A and rows in B): 
	  
	  // lda A = n1 no. cols in A.
	  // lda B = n2 no. cols in B.
	  // lda C = n2 no cols in C. 

	  // Taking the trace
	  int il=0;
	  for (int ii=0;ii<Nst_bl2_disc;ii++){
	    il=ij2r(Nst_bl2_disc,ii,ii);
	    Termo2DDADD+=Aux2dd[il];
	  }
	  
	  
	  if (Nst_bl1_kept!=0){
	    // Added condition if Nst_bl2_kept!=0 !

	    vector<double> opA2kd,opB2dkt,opA2kdw;
	    vector<double> Prod2dk(Nst_bl2_disc*Nst_bl1_kept);


	    Op1N[Nshell].MatBlock2STL(opA2kd,ibl1,ibl2,true,false);
	    // Nst_bl1_kept x Nst_bl2_disc matrix (assuming real!)
	    
	    if (NonDiagGF)
	      Op2N[Nshell].MatBlock2STL(opB2dkt,ibl1,ibl2,true,false);
	    // Nst_bl1_kept x Nst_bl2_disc matrix (assuming real!)
	    else
	      opB2dkt=opA2kd;

	  

	    // Positive and negative energies matrix elements
	    // Note for the future: sign is "-" if bosonic operators

	    cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasTrans, 
			Nst_bl2_disc, Nst_bl1_kept, Nst_bl2_disc,
			ALPHA,
			&Rho2dd[0], Nst_bl2_disc, 
			&opB2dkt[0],Nst_bl2_disc,
			BETA, 
			&Prod2dk[0], Nst_bl1_kept);
	    
	    // A = Rho2dd (Nst_bl2_disc x Nst_bl2_disc)
	    // B = opB2dkt(Nst_bl1_kept x Nst_bl2_disc)
	    // C = Prod2dk(Nst_bl2_disc x Nst_bl1_kept)

	    // A(n0 x n1) . (B_(n2 x n1))^T = C_(n0 x n2)

	    //M=n0 -> No of Rows in A and C 
	    //N=n2 -> No of Cols in C (=no rows in ORIGINAL B).
	    //K=n1 -> No of Cols in A (= no cols in ORIGINAL B).

	    //lda A = n1 no cols in  A.
	    //lda B = n1 no cols in ORIGINAL B.
	    //lda C = n2 no cols in C

	    if (CalcNorm) opA2kdw=opA2kd; 
	    else
	      MijxBDeltaEij(opA2kdw,opA2kd,(&AcutN[Nshell]),
			    ibl1,ibl2,true,false,omegabar); 
	    // Nb1_kept x Nb2_disc

	    cblas_dgemm(CblasRowMajor,  CblasNoTrans, CblasNoTrans, 
			Nst_bl2_disc, Nst_bl2_disc, Nst_bl1_kept,
			ALPHA,
			&Prod2dk[0],Nst_bl1_kept,
			&opA2kdw[0], Nst_bl2_disc, 
			BETA, 
			&Aux2dd[0], Nst_bl2_disc);

	    // A = Prod2dk_{Nst_bl2_disc x Nst_bl1_kept}
	    // B = opA2kdw_{Nst_bl1_kept x Nst_bl2_disc}
	    // C = Aux2dd_{Nst_bl2_disc x Nst_bl2_disc}
	    //
	    // A(n0 x n1).B( n1 x n2 ) = C(n0 x n2)

	    // [DEFAULT]
	    // M=n0 -> No of Rows in A and C: 
	    // N=n2 -> No of Cols in B and C: 
	    // K=n1 -> No of Cols in A and rows in B): 
	  
	    // lda A = n1 no. cols in A.
	    // lda B = n2 no. cols in B.
	    // lda C = n2 no cols in C. 

	    // Taking the trace
	    int il=0;
	    for (int ii=0;ii<Nst_bl2_disc;ii++){
	      il=ij2r(Nst_bl2_disc,ii,ii);
	      Termo2DDAKD+=Aux2dd[il];
	    }

  
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

