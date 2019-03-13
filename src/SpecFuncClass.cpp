
#include <cmath>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "SpecFuncClass.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)


#include "MyStructures.h"



using namespace std;


///////////////////
double CSpecFunction::HalfLambdaFactor(){

  if (Lambda!=0.0)
    return(0.5*(1.0+(1.0/Lambda)));
  else
    return(0.0);
}

///////////////////
double CSpecFunction::CalcDN(int Nsites){

//   double aux=HalfLambdaFactor()*pow(Lambda,(-(Nsites-1)/2.0) );
  double aux=HalfLambdaFactor()*pow(Lambda,(-(Nsites-1)/2.0)-z_twist+1.0 );

  return(aux);


}

///////////////////
int CSpecFunction::CalcNfromOmega(double omega){

  int iaux=(int)(2.0*log(2.0*abs(omega)/HalfLambdaFactor())/log(Lambda));
  iaux=1-iaux;
  if (iaux<0) iaux=0;

  return(iaux);

}


///////////////////
double  CSpecFunction::SumDeltas(double bbroad,
				 vector<double> Msq, 
				 vector<double> En,
				 double omega){

  double sum=0.0;

  if (Msq.size()!=En.size()){return(-1.0);}
  int Nterms=En.size();
  for (int ii=0;ii<Nterms;ii++){
    sum+=Msq[ii]*BDelta(omega,En[ii],bbroad);
  }

  return(sum);


}
///////////////////
vector< vector<double> >  CSpecFunction::CalcSpecCosti_Teq0(){

  double bbroad=0.5*log(Lambda);
  double factorWN=pow(Lambda,1.25);

  return(CalcSpecCosti_Teq0(bbroad,factorWN));

}


///////////////////
vector< vector<double> >  CSpecFunction::CalcSpecCosti_Teq0(double bbroad,
					 double factorWN){

  vector< vector<double> > Rhow;


  bool printstuff=false;

  //AcutN[0].PrintEn();
  
  //Op1N[0].PrintQNumbers();

  // Works only with EVEN sites!!
  int Nshell0=(NshellMin%2==0?NshellMin:NshellMin+1);

  int iNsh=0;
  //for (int Nshell=0;Nshell<NshellMax;Nshell+=2){
  //for (int Nshell=NshellMin;Nshell<NshellMax;Nshell+=2){
  for (int Nshell=Nshell0;Nshell<NshellMax;Nshell+=2){
    double DN=CalcDN(Nshell);
    double WN=factorWN*DN;

    // Get iGS and blGS

    int iGS=AcutN[Nshell].GetiGS();
    int posInBl=0;
    int iblGS=AcutN[Nshell].GetBlockFromSt(iGS,posInBl);

    // Debug
 //    if (Nshell==0) printstuff=true; else printstuff=false;


    if (printstuff){
      //AcutN[Nshell].PrintEn();
      cout << " Costi_Teq0: Nshell = " << Nshell << endl;
      cout << " iGS = " << iGS 
	   << " iblGS = " << iblGS 
	   << endl;
    }

    // SU(2) symmetry
    double CGfactor=1.0;
    double SblGS=0.0;
    if (AcutN[Nshell].totalS)
      SblGS=AcutN[Nshell].GetQNumber(iblGS,AcutN[Nshell].Sqnumbers[0]); 


    vector< vector<int> > BlFromGS_Op1;
    vector< vector<int> > BlToGS_Op1;

    vector< vector<int> >::iterator iitarray;

    // check blocks in Op1 that connect to iblGS;
    Op1N[Nshell].GetConnectingBlocks(iblGS,BlFromGS_Op1);
    Op1N[Nshell].GetConnectingBlocks(iblGS,BlToGS_Op1,true);
//     if (NonDiagGF){
//       Op2N[Nshell].GetConnectingBlocks(iblGS,BlFromGS_Op2);
//       Op2N[Nshell].GetConnectingBlocks(iblGS,BlToGS_Op2,true);
//     }
    // end if NonDiagGF

      if (printstuff) {
	cout << "BlFromGS size : " << BlFromGS_Op1.size() << endl;
        cout << "BlToGS size : " << BlToGS_Op1.size() << endl;
	cout << "Op1 MatBlockMap :" << endl;
	for (int iit=0;iit<Op1N[Nshell].MatBlockMap.size();iit+=2){
	  cout << Op1N[Nshell].MatBlockMap[iit] << " " <<Op1N[Nshell].MatBlockMap[iit+1] << endl;
	}
      }

    vector<double> En_0i;
    vector<double> En_i0;
    vector<double> Msq_0i;
    vector<double> Msq_i0;

    double MaxMsq=0.0;
    double EnMax=0.0;

    double auxM=0.0;
    complex<double> cauxM=ZeroC;

    // Loop in connecting blocks (positive omega) - WATCH OUT! This vector can be null
    if (printstuff)   cout << " Loop in connecting blocks: <iGS|Op1|ibl>  " << endl;
    for (iitarray=BlFromGS_Op1.begin();
	 iitarray<BlFromGS_Op1.end();iitarray++){
      int ibl=(*iitarray)[0];
      int iMatBl=(*iitarray)[1];

      if (printstuff)
           cout << "Positive w: Block : " <<  ibl
 	   << " MatBlock: " << iMatBl
 	   << endl;

      // SU(2) symmetry
      double Sbl=0.0;
      if (AcutN[Nshell].totalS){
	CGfactor=0.0;
	Sbl=AcutN[Nshell].GetQNumber(ibl,AcutN[Nshell].Sqnumbers[0]); 
	// only a single SU(2) for now
	for (double SzblGS=-SblGS;SzblGS<=SblGS;SzblGS+=1.0){
	  double dSigma=0.5; 
	  double Szbl=SzblGS+dSigma;
	  double auxCG=CGordan(SblGS,SzblGS, 0.5, dSigma, Sbl, Szbl);
	  CGfactor+=auxCG*auxCG;
	}
	if (dEqual(CGfactor,0.0)){cout << "Ops. CGfactor = 0.0 " << endl;}
      }
      // end if totalS

      if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(iGS,ibl)<0) ){
	cout << " Op2 does not connect blocks " 
	     << iGS << "  " << ibl << endl;
      } 
      else { // Loop in block states


	for (int ist=AcutN[Nshell].GetBlockLimit(ibl,0);
	     ist<=AcutN[Nshell].GetBlockLimit(ibl,1);
	     ist++){

// 	  double auxM=Op1N[Nshell].GetMatEl(iGS,ist);
// 	  if (NonDiagGF){auxM*=Op2N[Nshell].GetMatEl(ist,iGS);}
// 	  else{auxM*=auxM;}

	  if (Op1N[Nshell].IsComplex){
	    cauxM=Op1N[Nshell].cGetMatEl(iGS,ist);
	    if (NonDiagGF){cauxM*=conj(Op2N[Nshell].cGetMatEl(ist,iGS));}
	    else{cauxM*=conj(cauxM);}
	    auxM=cauxM.real(); // Only takes the real part! 
	                       // Watch out for NonDiagGF AND IsComplex
	  }
	  else {
	    auxM=Op1N[Nshell].GetMatEl(iGS,ist);
	    if (NonDiagGF){auxM*=Op2N[Nshell].GetMatEl(ist,iGS);}
	    else{auxM*=auxM;}
	  }

	  if (printstuff)
	    cout << " ist = " << ist 
		 << " En_0i = " << AcutN[Nshell].dEn[ist]
		 << " Op1_0i = " << auxM
		 << " Op1_0i(complex) = " << cauxM
		 << endl;


	  En_0i.push_back(AcutN[Nshell].dEn[ist]);
	  // Don't forget the 1/DN factor
	  //Msq_0i.push_back(auxM/DN);
	  // SU(2) symmetry
	  Msq_0i.push_back(CGfactor*auxM/DN);
	  
	  // debug
	  if ((auxM/DN)>MaxMsq){
	    MaxMsq=auxM/DN;
	    EnMax=AcutN[Nshell].dEn[ist];
	  }

	}
	// end loop in block ibl
      }
      // end if NonDiagGF and Op2 does not connect blocks
    }
    // end loop in connecting blocks (positive omega)

//     cout << " Nshell = " << Nshell 
// 	 << " MaxMsq = " << MaxMsq 
// 	 << " En = " << EnMax 
// 	 << endl;

    // Loop in connecting blocks (negative omega)
    if (printstuff)  cout << " Loop in connecting blocks: <ibl|Op1|iGS>  " << endl;
    for (iitarray=BlToGS_Op1.begin();
	 iitarray<BlToGS_Op1.end();iitarray++){
      int ibl=(*iitarray)[0];
      int iMatBl=(*iitarray)[1];

      if (printstuff)
           cout << "Negative w: Block : " <<  ibl
 	   << " MatBlock: " << iMatBl
 	   << endl;

      double Sbl=0.0;
      if (AcutN[Nshell].totalS){
	CGfactor=0.0;
	Sbl=AcutN[Nshell].GetQNumber(ibl,AcutN[Nshell].Sqnumbers[0]); 
	// only a single SU(2) for now
	for (double Szbl=-Sbl;Szbl<=Sbl;Szbl+=1.0){
	  double dSigma=0.5; 
	  double SzblGS=Szbl+dSigma;
	  double auxCG=CGordan(Sbl, Szbl, 0.5, dSigma, SblGS,SzblGS);
	  CGfactor+=auxCG*auxCG;
	}
	if (dEqual(CGfactor,0.0)){cout << "Ops. CGfactor = 0.0 " << endl;}
      }
      // end if totalS


      if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl,iGS)<0) ){
	cout << " Op2 does not connect blocks " 
	     << ibl << "  " << iGS << endl;
      } 
      else { // Loop in block states
	for (int ist=AcutN[Nshell].GetBlockLimit(ibl,0);
	     ist<=AcutN[Nshell].GetBlockLimit(ibl,1);
	     ist++){
//	  double auxM=Op1N[Nshell].GetMatEl(ist,iGS);
//	  if (NonDiagGF){auxM*=Op2N[Nshell].GetMatEl(ist,iGS);}
//	  else{auxM*=auxM;}

	  if (Op1N[Nshell].IsComplex){
	    cauxM=Op1N[Nshell].cGetMatEl(ist,iGS);
	    if (NonDiagGF){cauxM*=conj(Op2N[Nshell].cGetMatEl(iGS,ist));}
	    else{cauxM*=conj(cauxM);}
	    auxM=cauxM.real(); // Only takes the real part! 
	                       // Watch out for NonDiagGF AND IsComplex
	  }
	  else {
	    auxM=Op1N[Nshell].GetMatEl(ist,iGS);
	    if (NonDiagGF){auxM*=Op2N[Nshell].GetMatEl(iGS,ist);}
	    else{auxM*=auxM;}
	  }

	  if (printstuff)
	    cout << " ist = " << ist 
		 << " En_0i = " << AcutN[Nshell].dEn[ist]
		 << " Op1_0i = " << auxM
		 << " Op1_0i(complex) = " << cauxM
		 << endl;


	  En_i0.push_back(AcutN[Nshell].dEn[ist]);
	  // Don't forget the 1/DN factor
	  //Msq_i0.push_back(auxM/DN);
	  // SU(2) symmetry
	  Msq_i0.push_back(CGfactor*auxM/DN);

	}
	// end loop in block ibl
      }
      // end if NonDiagGF and Op2 does not connect blocks
    }
    // end loop in connecting blocks (negative omega)

    // valid only for Q,Sz. Need to improve on this... (e.g., totS).
    //double ZN=AcutN[Nshell].PartitionFuncTeq0(); // Not so good...
    double ZN=AcutN[Nshell].PartitionFunc(1.0e3);
    double auxRho=0.0;
    Rhow.push_back( vector<double> () );
    // positive omega
    Rhow[iNsh].push_back(WN);
    auxRho=SumDeltas(bbroad,Msq_0i,En_0i,factorWN);
    Rhow[iNsh].push_back(auxRho/ZN);
    if (printstuff){
      cout << "w+ : auxRho = " << auxRho << " ZN = " << ZN 
	   << " auxRho/ZN = " << auxRho/ZN 
	   << " ZN_Teq0 = " << AcutN[Nshell].PartitionFuncTeq0() << endl;
    }
    // negative omega
    Rhow[iNsh].push_back(-WN);
    auxRho=SumDeltas(bbroad,Msq_i0,En_i0,factorWN);
    Rhow[iNsh].push_back(auxRho/ZN);
    if (printstuff){
      cout << "w- : auxRho = " << auxRho << " ZN = " << ZN 
	   << " auxRho/ZN = " << auxRho/ZN << endl;
    }

    iNsh++;
  }
  // end loop in Nshell

  return(Rhow);


}

///////////////////
vector < vector<double> >  CSpecFunction::CalcSpecCosti_T(int Mtemp){


  vector < vector<double> > RhoCostiT;
  double auxRho=0.0;
  cout << "CalcSpecCosti_T: Implementing it " << endl;

  // Calculate Temp

  double DM=1.0;
  double Temp=0.0;
  double betabar=0.727; // OR should I add as a parameter?
  double factorWN=pow(Lambda,1.25);  // OR should I add as a parameter?

  bool ZeroTemp=false;
  if (Mtemp>200){
    cout << "CalcSpecCosti_T: Zero Temp " << endl;
    ZeroTemp=true;
    Temp=0.0;
    betabar=1.0e10;
  }
  else{
    DM=CalcDN(Mtemp);
    Temp=DM/betabar;
  }


  // Working with EVEN sites (can do better)!!
  int Nshell0=(NshellMin%2==0?NshellMin:NshellMin+1);

  int iNsh=0;
  //for (int Nshell=0;Nshell<NshellMax;Nshell+=2){
  //for (int Nshell=NshellMin;Nshell<NshellMax;Nshell+=2){
  for (int Nshell=Nshell0;Nshell<(Mtemp>NshellMax?NshellMax:Mtemp);Nshell+=2){
    double DN=CalcDN(Nshell);
    double WN=factorWN*DN;

    if (ZeroTemp) betabar=1.0e10; else betabar=DN/Temp;

    RhoCostiT.push_back( vector<double> () );
    // positive omega
    RhoCostiT[iNsh].push_back(WN);
    auxRho=CalcSpecCosti_T_N(Nshell,betabar,factorWN);
    RhoCostiT[iNsh].push_back(auxRho/DN);
    // negative omega
    RhoCostiT[iNsh].push_back(-WN);
    auxRho=CalcSpecCosti_T_N(Nshell,betabar,-factorWN);
    RhoCostiT[iNsh].push_back(auxRho/DN);

    iNsh++;
  }

  return(RhoCostiT);


}

/////////////////////

double CSpecFunction::CalcSpecCosti_T_N(int Nshell, double betabar,
					double factorWN,bool CalcNorm){
  //
  // Finite-T Costi calculation for Nshell (<Mtemp)
  //
  // Already includes ZN!

  double rhoCosti=0.0;

  bool printstuff=false;
  // Debug
  //if ( (Nshell==0) ) printstuff=true; else printstuff=false;



  // Need to Check Sync!!
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));

  if ( (!ChkSyncOp1Acut) ){
    cout << "CalcSpecCosti_T_N: Op1N not in sync with AcutN " << endl;
    return(0.0);
  }


  for (int iMatbl=0;iMatbl<Op1N[Nshell].NumMatBlocks();iMatbl++) {
    int ibl1=0;
    int ibl2=0;

    Op1N[Nshell].GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);


    // Ok, remember the input Op2N is the COMPLEX CONJUGATE of what we want!
    //  
    //  So (B)_{ibl2 ibl1} = ((B+)_{ibl1 ibl2})* <- THIS is Op2

    //if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl2,ibl1)<0) ){
    if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl1,ibl2)<0) ){
      cout << "CalcSpecCosti_T_N: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return(0.0);
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);  
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 

    // Check block sizes
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "CalcSpecCosti_T_N: Block in Op1 not the same size as block in Op2 " << endl;
      return(0.0);
    }


    //////////

    // Blas matrices
    boost::numeric::ublas::matrix<double> Rho1(Nst_bl1,Nst_bl1);
    boost::numeric::ublas::matrix<double> Rho2(Nst_bl2,Nst_bl2);

    // In the Costi scheme, these are DIAGONAL matrices e^{-betabar*Ei}!!
    // New subroutine
    
    Rho1=AcutN[Nshell].ExpEi2BLAS(ibl1,betabar); // Nst_bl1 x Nst_bl1
    Rho2=AcutN[Nshell].ExpEi2BLAS(ibl2,betabar); // Nst_bl2 x Nst_bl2


    double Trace=0.0; // Calculating the final trace.

    if ( !(Op1N[Nshell].IsComplex) ){

      ///////// Real operators ////////////


      boost::numeric::ublas::matrix<double> opAw(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<double> opA(Nst_bl1,Nst_bl2);

      boost::numeric::ublas::matrix<double> opB(Nst_bl2,Nst_bl1);

      boost::numeric::ublas::matrix<double> Sum12(Nst_bl2,Nst_bl1);


      opA=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2);

      if (NonDiagGF)
	//opB=Op2N[Nshell].MatBlock2BLAS(ibl2,ibl1);
	opB=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2)); // Nst_bl2 x Nst_bl1 matrix (assuming real!)
      else
	opB=trans(opA);
    
      // Positive and negative energies matrix elements
      // Note for the future: sign is "-" if bosonic operators
      noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix

      // Calculate A_ij * Delta(w-(Ej-Ei))

      if (CalcNorm) opAw=opA; else
      noalias(opAw)=MijxBDeltaEij(opA,(&AcutN[Nshell]),ibl1,ibl2,factorWN);

      // Take the trace: (use preallocated Rho2 as temp...)

      noalias(Rho2)=prod(Sum12,opAw);
    

      boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	diag(Rho2,
	     boost::numeric::ublas::range(0,Rho2.size1()),
	     boost::numeric::ublas::range(0,Rho2.size2())); 

      Trace=sum(diag);

    } else {
    //////
    /// COMPLEX operators
    /////

//       // Blas matrices
//       boost::numeric::ublas::matrix<complex<double> > cRho1(Nst_bl1,Nst_bl1);
// Just checking...
      boost::numeric::ublas::matrix<complex<double> > cRho2(Nst_bl2,Nst_bl2);

      boost::numeric::ublas::matrix<complex<double> > copAw(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > copA(Nst_bl1,Nst_bl2);

      boost::numeric::ublas::matrix<complex<double> > copB(Nst_bl2,Nst_bl1);

      boost::numeric::ublas::matrix<complex<double> > cSum12(Nst_bl2,Nst_bl1);

//       cRho1=RhoN[Nshell].cMatBlock2BLAS(ibl1,ibl1);
//       cRho2=RhoN[Nshell].cMatBlock2BLAS(ibl2,ibl2);

      copA=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2);


      if (NonDiagGF)
	copB=herm(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2)); // Nst_bl2 x Nst_bl1 matrix 
      else
	copB=herm(copA);

      if (printstuff){
	cout << " Costi_T: Nshell= " << Nshell  
	     << " omegabar = " << factorWN
	     << " betabar = " << betabar
	     << " ibl1 = " << ibl1
	     << " ibl2 = " << ibl2
	     << endl;
	Op1N[Nshell].PrintBlockQNumbers(ibl1);
	Op1N[Nshell].PrintBlockQNumbers(ibl2);
	if ( (ibl1==0)&&(ibl2==2) ){
	  cout << "Rho1 : " << endl << Rho1 << endl;
	  cout << "Rho2 : " << endl << Rho2 << endl;
	  cout << "opA : " << endl << copA << endl;
	  cout << "opB : " << endl << copB << endl;
	}
      }
      // end if printstuff


      // Positive and negative energies matrix elements
      // Note for the future: sign is "-" if bosonic operators
      noalias(cSum12)=prod(copB,Rho1)+prod(Rho2,copB); // Nst2 x Nst1 matrix

      //noalias(Sum12)=prod(opB,Rho1); // stronger asymmetry (negative higher)
      //noalias(Sum12)=prod(Rho2,opB); // weaker asymmetry (positive higher)

      // Calculate A_ij * Delta(w-(Ej-Ei))

      if (CalcNorm) copAw=copA; else
      noalias(copAw)=cMijxBDeltaEij(copA,(&AcutN[Nshell]),ibl1,ibl2,factorWN);

      if (printstuff){
	if ( (ibl1==0)&&(ibl2==2) ){
	  cout << "Sum12 : " << endl << cSum12 << endl;

	  cout << "opAw : " << endl << copAw << endl;

	}
      }
      // end printstuff

      // Take the trace: (use preallocated Rho2 as temp...)
      
      noalias(cRho2)=prod(cSum12,copAw);

      boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	cdiag(cRho2,
	     boost::numeric::ublas::range(0,cRho2.size1()),
	     boost::numeric::ublas::range(0,cRho2.size2())); 


      Trace=(sum(cdiag)).real();

    }
    // end if Op1 is real else.

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

//     rhoCosti+=CGfactor*sum(diag);
    rhoCosti+=CGfactor*Trace;
 
    if (printstuff){
      cout << " Costi_T: trace= " << Trace
	   << " Accumulated RhoCosti= " << rhoCosti
	   << endl;
    }


  }
  // end loop in Op1N blocks

  double ZN=AcutN[Nshell].PartitionFunc(betabar);

  if (printstuff)
    cout << " CalcSpecCosti_T_N: Nshell= " << Nshell
	 << " factorWN=" << factorWN 
	 << " rhoCosti = " << rhoCosti/ZN
	 << " ZN = " << ZN  << endl;

  return(rhoCosti/ZN);

}
///////////////////
///////////////////
///////////////////
///////////////////
///////////////////
///////////////////


///////////////////
void  CSpecFunction::CalcCondSIAM_ManyTs(int NwEachN){

  // Use Omega_Rhow vector to store

  Omega_Rhow_Even.clear();
  Omega_Rhow_Odd.clear();
  Omega_Rhow.clear();

  Omega_Rhow.push_back( vector<double> () ); // Temps
  Omega_Rhow.push_back( vector<double> () ); // CondN (Raw!)


  double auxCond=0.0;
  cout << "CalcCondSIAM_ManyTs: Implementing it " << endl;

  // Test
  int Nmax=NshellMax-1;

  int Nsh0=NshellMin;
  // if (UseCFS==1){Nsh0+=2;}


  // First shell

  //   double DN=CalcDN(NshellMin);
  double DN=CalcDN(Nsh0);
  double TempN=DN/Betabar;

  // "Positive" temperatures (decreasing)
  cout << "CalcCondSIAM_ManyTs : Loop in Nshells " << endl;
  for (int Nsh=Nmax;Nsh>=Nsh0;Nsh--){
    DN=CalcDN(Nsh);
    TempN=DN/Betabar;
    cout << " N = " << Nsh 
	 << " DN = " << DN 
	 << " TempN = " << TempN 
	 << endl;

    for (int itemp=NwEachN-1;itemp>=0;itemp--){
      double RedFactor=(double)itemp/(2.0*NwEachN);
      double betatemp=Betabar*pow(Lambda,RedFactor); // betabar increases a bit
      double ThisTemp=DN/betatemp;
      double gT=CalcCondSIAM_T_N(Nsh,betatemp);

      Omega_Rhow[0].push_back(ThisTemp);
      Omega_Rhow[1].push_back(gT*M_PI/ThisTemp);
      //Omega_Rhow[1].push_back(gT*M_1_PI*betatemp/(DN*DN));
    }
    // end loop in shells (temperature)

  }
  // end loop in Nshells

}
// end function
/////////////////////

double CSpecFunction::CalcCondSIAM_T_N(int Nshell, double betabar, bool CalcNorm){
  //
  // Finite-T conductance calculation for Nshell (<Mtemp)
  //
  // Already includes ZN and beta*pi!

  double CondT=0.0;

  bool printstuff=true;
  // Debug
  //if ( (Nshell==0) ) printstuff=true; else printstuff=false;



  // Need to Check Sync!!
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));


  if ( (!ChkSyncOp1Acut) ){
    cout << "CalcCondSIAM_T_N: Op1N not in sync with AcutN " << endl;
    return(0.0);
  }


  for (int iMatbl=0;iMatbl<Op1N[Nshell].NumMatBlocks();iMatbl++) {
    int ibl1=0;
    int ibl2=0;

    Op1N[Nshell].GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);



    // Ok, remember the input Op2N is the COMPLEX CONJUGATE of what we want!
    //  
    //  So (B)_{ibl2 ibl1} = ((B+)_{ibl1 ibl2})* <- THIS is Op2

    //if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl2,ibl1)<0) ){
    if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl1,ibl2)<0) ){
      cout << "CalcCondSIAM_T_N: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return(0.0);
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);  
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 

    // Check block sizes
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "CalcCondSIAM_T_N: Block in Op1 not the same size as block in Op2 " << endl;
      return(0.0);
    }


    //////////

    // Blas matrices
    boost::numeric::ublas::matrix<double> Diag1(Nst_bl1,Nst_bl1);

    double Trace=0.0; // Calculating the final trace.

    if ( !(Op1N[Nshell].IsComplex) ){

      ///////// Real operators ////////////

      boost::numeric::ublas::matrix<double> opAw(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<double> opA(Nst_bl1,Nst_bl2);

      boost::numeric::ublas::matrix<double> opB(Nst_bl2,Nst_bl1);

      opA=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2); // Nst_bl1 x Nst_bl2

      if (NonDiagGF)
	opB=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2)); // Nst_bl2 x Nst_bl1 matrix (assuming real!)
      else
	opB=trans(opA);
    
      // Calculate A_ij/(exp(betabar*Ei) + exp(betabar*Ej) )

      if (CalcNorm) opAw=opA; else
      noalias(opAw)=AijOverexpEij(opA,(&AcutN[Nshell]),ibl1,ibl2,betabar); 
      // Nst_bl1 x Nst_bl2 

      // Take the trace: (use preallocated Diag1 as temp...)

      noalias(Diag1)=prod(opAw,opB); // Nst_bl1 x Nst_bl1 

    
      boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	diag(Diag1,
	     boost::numeric::ublas::range(0,Diag1.size1()),
	     boost::numeric::ublas::range(0,Diag1.size2())); 

      Trace=sum(diag);

    } else {
    //////
    /// COMPLEX operators
    /////

//       // Blas matrices
//       boost::numeric::ublas::matrix<complex<double> > cRho1(Nst_bl1,Nst_bl1);
// Just checking...
      boost::numeric::ublas::matrix<complex<double> > cDiag1(Nst_bl1,Nst_bl1);

      boost::numeric::ublas::matrix<complex<double> > copAw(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > copA(Nst_bl1,Nst_bl2);

      boost::numeric::ublas::matrix<complex<double> > copB(Nst_bl2,Nst_bl1);


      copA=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2);


      if (NonDiagGF)
	copB=herm(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2)); // Nst_bl2 x Nst_bl1 matrix 
      else
	copB=herm(copA);

      // Calculate A_ij/(exp(betabar*Ei) + exp(betabar*Ej) )

      if (CalcNorm) copAw=copA; else
      noalias(copAw)=cAijOverexpEij(copA,(&AcutN[Nshell]),ibl1,ibl2,betabar);


      // Take the trace: (use preallocated Rho2 as temp...)
      
      noalias(cDiag1)=prod(copAw,copB);

      boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	cdiag(cDiag1,
	     boost::numeric::ublas::range(0,cDiag1.size1()),
	     boost::numeric::ublas::range(0,cDiag1.size2())); 


      Trace=(sum(cdiag)).real();

    }
    // end if Op1 is real else.

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

    CondT+=CGfactor*Trace;
 
    if (printstuff){
      cout << "Block " << iMatbl << " of " << Op1N[Nshell].NumMatBlocks()-1 << endl
	   << " Cond_T: trace_bl= " << Trace << " x CGfactor=" << CGfactor
	   << " Accumulated CondT= " << CondT
	   << endl;
    }
    //end print

  }
  // end loop in Op1N blocks


  double ZN=AcutN[Nshell].PartitionFunc(betabar);

  if (printstuff)
    cout << " CalcCondSIAM_T_N: Nshell= " << Nshell
	 << " betabar=" << betabar 
	 << " CondT = " << CondT
	 << " ZN = " << ZN  
	 << " CondT/ZN = " << CondT/ZN
	 << endl;

  return(CondT/ZN);

}




///////////////////
///////////////////
///////////////////
///////////////////
///////////////////
///////////////////

boost::numeric::ublas::matrix<double> CSpecFunction::MijxBDeltaEij(CNRGmatrix* pMat, CNRGarray* pAeig, int ibl1, int ibl2, double omega){

  bool ChkSync=pMat->ChkSync(pAeig);
  int iMatBl=pMat->FindMatBlock(ibl1,ibl2);

  if ( (!ChkSync)||(iMatBl<0) ){
    cout << " Cannot calculate Mij*BD(Ei-Ej) " << endl;
    return(boost::numeric::ublas::matrix<double> () );
  }

  boost::numeric::ublas::matrix<double> Maux=pMat->MatBlock2BLAS(ibl1,ibl2);

  return( MijxBDeltaEij(Maux,pAeig,ibl1,ibl2,omega) );

}
//

boost::numeric::ublas::matrix<double> CSpecFunction::MijxBDeltaEij(boost::numeric::ublas::matrix<double> Mat, CNRGarray* pAeig, int ibl1, int ibl2, double omega){

  int Nst1=pAeig->GetBlockSize(ibl1);
  int Nst2=pAeig->GetBlockSize(ibl2);


  if ( (Nst1!=Mat.size1())||(Nst2!=Mat.size2()) ){
    cout << " Cannot calculate Mij*BD(Ei-Ej) " << endl;
    return(boost::numeric::ublas::matrix<double> () );
  }

  boost::numeric::ublas::matrix<double> Maux(Nst1,Nst2);

  int ist0=pAeig->GetBlockLimit(ibl1,0);
  int jst0=pAeig->GetBlockLimit(ibl2,0);

  double Ei=0.0;
  double Ej=0.0;

  bool print=false;

  double BD=0.0;
  double ThisDN=CalcDN(pAeig->Nshell);
//   double TempDN=Temp/CalcDN(pAeig->Nshell);
  double TempDN=Temp/ThisDN;
  double GapDN=Gap/ThisDN;

  for (int ii=0;ii<Nst1;ii++){
    Ei=pAeig->dEn[ist0+ii];
    if (print) cout << " Ei = " << Ei;
    for (int jj=0;jj<Nst2;jj++){
      Ej=pAeig->dEn[jst0+jj];
      if (print) cout << " Ej = " << Ej;
      //Maux(ii,jj)=Mat(ii,jj)*BDelta(omega,(Ej-Ei),dBroad); // correct for T=0
      // NEED TO INCLUDE TEMPERATURE HERE!
      BD=BDelta(omega,(Ej-Ei),dBroad);
      //if (fabs(Ej-Ei)<TempDN){
      //if (fabs(Ej-Ei)<4.0*TempDN){
      // NEW (Jan 2016): No broadening if it is inside the Gap
      if (fabs(Ej-Ei)<GapDN) BD=0.0; // Will give zero...
      if (fabs(Ej-Ei)<TwindowFac*TempDN){
 	//BD=BDeltaTemp(omega,(Ej-Ei),0.1*TempDN); //
 	//BD=BDeltaTemp(omega,(Ej-Ei),TempDN); //
 	//BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp*TempDN); //
 	BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp); // Anders
	//double btemp=(TempDN<dBroadTemp?sqrt(dBroadTemp*TempDN):dBroadTemp);
 	//BD=BDeltaTemp(omega,(Ej-Ei),btemp); // Works if T is small!
      }
      // end if |Ej-Ei|<TempDN
      Maux(ii,jj)=Mat(ii,jj)*BD; // Uncomment this to work

      // END TEMP
    }
    if (print) cout << endl;
  }
  // end loop in i,j


  return(Maux);

}
//

//// COMPLEX
boost::numeric::ublas::matrix<complex<double> > CSpecFunction::cMijxBDeltaEij(CNRGmatrix* pMat, CNRGarray* pAeig, int ibl1, int ibl2, double omega){

  bool ChkSync=pMat->ChkSync(pAeig);
  int iMatBl=pMat->FindMatBlock(ibl1,ibl2);

  if ( (!ChkSync)||(iMatBl<0) ){
    cout << " Cannot calculate Mij*BD(Ei-Ej) " << endl;
    return(boost::numeric::ublas::matrix<complex<double> > () );
  }

  boost::numeric::ublas::matrix<complex<double> > cMaux=pMat->cMatBlock2BLAS(ibl1,ibl2);

  return( cMijxBDeltaEij(cMaux,pAeig,ibl1,ibl2,omega) );

}
//


///
boost::numeric::ublas::matrix<complex<double> > CSpecFunction::cMijxBDeltaEij(boost::numeric::ublas::matrix<complex<double> > cMat, CNRGarray* pAeig, int ibl1, int ibl2, double omega){

  int Nst1=pAeig->GetBlockSize(ibl1);
  int Nst2=pAeig->GetBlockSize(ibl2);


  if ( (Nst1!=cMat.size1())||(Nst2!=cMat.size2()) ){
    cout << " Cannot calculate Mij*BD(Ei-Ej) " << endl;
    return(boost::numeric::ublas::matrix<complex<double> > () );
  }

  boost::numeric::ublas::matrix<complex<double> > cMaux(Nst1,Nst2);

  int ist0=pAeig->GetBlockLimit(ibl1,0);
  int jst0=pAeig->GetBlockLimit(ibl2,0);

  double Ei=0.0;
  double Ej=0.0;

  bool print=false;

  // Debug
  //if ((pAeig->Nshell==0)&&(ibl1==0)&&(ibl2==2)) print=true;

  double BD=0.0;
//  double TempDN=Temp/CalcDN(pAeig->Nshell);
  double ThisDN=CalcDN(pAeig->Nshell);
  double TempDN=Temp/ThisDN;
  double GapDN=Gap/ThisDN;



  for (int ii=0;ii<Nst1;ii++){
    Ei=pAeig->dEn[ist0+ii];
    if (print) cout << " Ei = " << Ei;
    for (int jj=0;jj<Nst2;jj++){
      Ej=pAeig->dEn[jst0+jj];
      if (print) cout << " Ej = " << Ej;
      //Maux(ii,jj)=Mat(ii,jj)*BDelta(omega,(Ej-Ei),dBroad); // correct for T=0
      // NEED TO INCLUDE TEMPERATURE HERE!
      BD=BDelta(omega,(Ej-Ei),dBroad);
      // NEW (Jan 2016): No broadening if it is inside the Gap
      if (fabs(Ej-Ei)<GapDN) BD=0.0; // Will give zero...
      //if (fabs(Ej-Ei)<TempDN){
      //if (fabs(Ej-Ei)<4.0*TempDN){
      if (print) cout << endl << " omega = " << omega << " BD_T0= " << BD << " Temp= " << Temp << " TempDN = " << TempDN << endl;
      if (fabs(Ej-Ei)<TwindowFac*TempDN){
 	//BD=BDeltaTemp(omega,(Ej-Ei),0.1*TempDN); //
 	//BD=BDeltaTemp(omega,(Ej-Ei),TempDN); //
 	//BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp*TempDN); // 
 	//BD=BDeltaTemp(omega,(Ej-Ei),dBroadTemp); // Anders
	double btemp=(TempDN<dBroadTemp?sqrt(dBroadTemp*TempDN):dBroadTemp);
 	BD=BDeltaTemp(omega,(Ej-Ei),btemp); // Works if T is small!
      }
      // end if |Ej-Ei|<TempDN
      cMaux(ii,jj)=cMat(ii,jj)*BD; // Uncomment this to work
      if (print) cout << "Mat("<<ii<<","<<jj<<")= "<<cMat(ii,jj)<< " - BD= " << BD << endl;

      // END TEMP
    }
    if (print) cout << endl;
  }
  // end loop in i,j

  return(cMaux);

}
//



///////////////////

boost::numeric::ublas::matrix<double> CSpecFunction::AijOverexpEij(CNRGmatrix* pMat, CNRGarray* pAeig, int ibl1, int ibl2, double betabar){

  bool ChkSync=pMat->ChkSync(pAeig);
  int iMatBl=pMat->FindMatBlock(ibl1,ibl2);

  if ( (!ChkSync)||(iMatBl<0) ){
    cout << " Cannot calculate Aij/(expEi + expEj) " << endl;
    return(boost::numeric::ublas::matrix<double> () );
  }

  boost::numeric::ublas::matrix<double> Maux=pMat->MatBlock2BLAS(ibl1,ibl2);

  return( AijOverexpEij(Maux,pAeig,ibl1,ibl2,betabar) );

}
//

boost::numeric::ublas::matrix<double> CSpecFunction::AijOverexpEij(boost::numeric::ublas::matrix<double> Mat, CNRGarray* pAeig, int ibl1, int ibl2, double betabar){

  int Nst1=pAeig->GetBlockSize(ibl1);
  int Nst2=pAeig->GetBlockSize(ibl2);


  if ( (Nst1!=Mat.size1())||(Nst2!=Mat.size2()) ){
    cout << " Cannot calculate Aij/(expEi + expEj) " << endl;
    return(boost::numeric::ublas::matrix<double> () );
  }

  boost::numeric::ublas::matrix<double> Maux(Nst1,Nst2);

  int ist0=pAeig->GetBlockLimit(ibl1,0);
  int jst0=pAeig->GetBlockLimit(ibl2,0);

  double Ei=0.0;
  double Ej=0.0;

  bool print=false;

  double Expi=0.0;
  double Expj=0.0;
  double ThisDN=CalcDN(pAeig->Nshell);

  for (int ii=0;ii<Nst1;ii++){
    Ei=pAeig->dEn[ist0+ii];
    if (print) cout << " Ei = " << Ei;
    for (int jj=0;jj<Nst2;jj++){
      Ej=pAeig->dEn[jst0+jj];
      if (print) cout << " Ej = " << Ej;
      Expi=exp(betabar*Ei);
      Expj=exp(betabar*Ej);
      Maux(ii,jj)=Mat(ii,jj)/(Expi+Expj); // Uncomment this to work
    }
    if (print) cout << endl;
  }
  // end loop in i,j


  return(Maux);

}
//

// COMPLEX (To do!! Nov 2018)

boost::numeric::ublas::matrix<complex<double> > CSpecFunction::cAijOverexpEij(CNRGmatrix* pMat, CNRGarray* pAeig, int ibl1, int ibl2, double betabar){

  bool ChkSync=pMat->ChkSync(pAeig);
  int iMatBl=pMat->FindMatBlock(ibl1,ibl2);

  if ( (!ChkSync)||(iMatBl<0) ){
    cout << " Cannot calculate Aij/(expEi + expEj) " << endl;
    return(boost::numeric::ublas::matrix<complex<double> > () );
  }

  boost::numeric::ublas::matrix<complex<double> > cMaux=pMat->cMatBlock2BLAS(ibl1,ibl2);

  return( cAijOverexpEij(cMaux,pAeig,ibl1,ibl2,betabar) );

}
//

boost::numeric::ublas::matrix<complex<double> > CSpecFunction::cAijOverexpEij(boost::numeric::ublas::matrix<complex<double> > cMat, CNRGarray* pAeig, int ibl1, int ibl2, double betabar){

  int Nst1=pAeig->GetBlockSize(ibl1);
  int Nst2=pAeig->GetBlockSize(ibl2);


  if ( (Nst1!=cMat.size1())||(Nst2!=cMat.size2()) ){
    cout << " Cannot calculate Aij/(expEi + expEj) " << endl;
    return(boost::numeric::ublas::matrix<complex<double> > () );
  }

  boost::numeric::ublas::matrix<complex<double> > cMaux(Nst1,Nst2);

  int ist0=pAeig->GetBlockLimit(ibl1,0);
  int jst0=pAeig->GetBlockLimit(ibl2,0);

  double Ei=0.0;
  double Ej=0.0;

  bool print=false;

  double Expi=0.0;
  double Expj=0.0;
  double ThisDN=CalcDN(pAeig->Nshell);

  for (int ii=0;ii<Nst1;ii++){
    Ei=pAeig->dEn[ist0+ii];
    if (print) cout << " Ei = " << Ei;
    for (int jj=0;jj<Nst2;jj++){
      Ej=pAeig->dEn[jst0+jj];
      if (print) cout << " Ej = " << Ej;
      Expi=exp(betabar*Ei);
      Expj=exp(betabar*Ej);
      cMaux(ii,jj)=cMat(ii,jj)/(Expi+Expj);
    }
    if (print) cout << endl;
  }
  // end loop in i,j


  return(cMaux);

}
//






///////////////////
void CSpecFunction::SpecBulla_SetMsqEdiff(int Mtemp){

  // TO BE DONE

}
//

///////////////////
double CSpecFunction::CalcSpecBulla_T(int Mtemp, double omega){


  // TO BE DONE


  double Rhow=0.0;

  return(Rhow);


}


///////////////////
void CSpecFunction::CalcSpecBulla_T_FixedOmegas(int Mtemp, double factorWN){


  // TO BE DONE



}



//////////////////
double CSpecFunction::DMNRG_SpecDens_M(double omegabar,int Nshell, bool CalcNorm){

  // Need to get the Connecting matrix elements and energies
  // to and from a given block (say, blGS)
  // Interesting procedure...

  // Loop in blocks of Operator


  double rho_wM=0.0;

  bool printstuff=false;
  // Debug
  //if ( (Nshell==0) ) printstuff=true; else printstuff=false;



  // Need to Check Sync!!
  bool ChkSyncRhoAcut=RhoN[Nshell].ChkSync((&AcutN[Nshell]));
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));

  if ( (!ChkSyncRhoAcut)||(!ChkSyncOp1Acut) ){
    cout << "DMNRG_SpecDens_M: RhoN or Op1N not in sync with AcutN " << endl;
    return(0.0);
  }

  if (printstuff)
    cout << " Nshell = " << Nshell
	 << " omegabar = " << omegabar;


  for (int iMatbl=0;iMatbl<Op1N[Nshell].NumMatBlocks();iMatbl++) {
    int ibl1=0;
    int ibl2=0;

    Op1N[Nshell].GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);


    // Ok, remember the input Op2N is the COMPLEX CONJUGATE of what we want!
    //  
    //  So (B)_{ibl2 ibl1} = ((B+)_{ibl1 ibl2})* <- THIS is Op2

    //if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl2,ibl1)<0) ){
    if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl1,ibl2)<0) ){
      cout << "DMNRG_SpecDens_M: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return(0.0);
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);  
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 

    // Check block sizes
    if ( (Nst_bl1!=RhoN[Nshell].GetBlockSize(ibl1))||
	 (Nst_bl2!=RhoN[Nshell].GetBlockSize(ibl2)) ){
      cout << "DMNRG_SpecDens_M: Block in Op1 not the same size as block in RhoN " << endl;
      return(0.0);
    }
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "DMNRG_SpecDens_M: Block in Op1 not the same size as block in Op2 " << endl;
      return(0.0);
    }



    double Trace=0.0; // Calculating the final trace.

    if ( (RhoN[Nshell].IsComplex)&&(Op1N[Nshell].IsComplex) ){

    //////
    /// COMPLEX operators
    /////


      // Blas matrices
      boost::numeric::ublas::matrix<complex<double> > cRho1(Nst_bl1,Nst_bl1);
      boost::numeric::ublas::matrix<complex<double> > cRho2(Nst_bl2,Nst_bl2);

      boost::numeric::ublas::matrix<complex<double> > copAw(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<complex<double> > copA(Nst_bl1,Nst_bl2);

      boost::numeric::ublas::matrix<complex<double> > copB(Nst_bl2,Nst_bl1);

      boost::numeric::ublas::matrix<complex<double> > cSum12(Nst_bl2,Nst_bl1);



      cRho1=RhoN[Nshell].cMatBlock2BLAS(ibl1,ibl1);
      cRho2=RhoN[Nshell].cMatBlock2BLAS(ibl2,ibl2);

      copA=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2);

      if (NonDiagGF)
	copB=herm(Op2N[Nshell].cMatBlock2BLAS(ibl1,ibl2)); // Nst_bl2 x Nst_bl1 matrix 
      else
	copB=herm(copA);

      if (printstuff){
	cout << " DMNRG_SpecDens_M: Nshell= " << Nshell
	     << " ibl1 = " << ibl1
	     << " ibl2 = " << ibl2
	     << endl;
	Op1N[Nshell].PrintBlockQNumbers(ibl1);
	Op1N[Nshell].PrintBlockQNumbers(ibl2);
	if ( (ibl1==0)&&(ibl2==2) ){
	  cout << "Rho1 : " << endl << cRho1 << endl;
	  cout << "Rho2 : " << endl << cRho2 << endl;
	  cout << "opA : " << endl << copA << endl;
	  cout << "opB : " << endl << copB << endl;
	}
      }
      //end if printstuff


      // Positive and negative energies matrix elements
      // Note for the future: sign is "-" if bosonic operators
      noalias(cSum12)=prod(copB,cRho1)+prod(cRho2,copB); // Nst2 x Nst1 matrix

      //noalias(Sum12)=prod(opB,Rho1); // stronger asymmetry (negative higher)
      //noalias(Sum12)=prod(Rho2,opB); // weaker asymmetry (positive higher)

      // Calculate A_ij * Delta(w-(Ej-Ei))

      if (CalcNorm) copAw=copA; else
      noalias(copAw)=cMijxBDeltaEij(copA,(&AcutN[Nshell]),ibl1,ibl2,omegabar);


      if (printstuff){
	if ( (ibl1==0)&&(ibl2==2) ){
	  cout << "Sum12 : " << endl << cSum12 << endl;

	//       cout << " Energies: " << endl;
	//       AcutN[Nshell].PrintBlockEn(ibl1);
	//       AcutN[Nshell].PrintBlockEn(ibl2);

	//       cout << " Delta functions Delta(omegabar-(Ej-Ei)) " 
	// 	   << " bbroad = " << dBroad << endl;
	//       double dEi=0.0;
	//       double dEj=0.0;
	//       int ist0=AcutN[Nshell].GetBlockLimit(ibl1,0);
	//       int jst0=AcutN[Nshell].GetBlockLimit(ibl2,0);

	//       for (int ii=0;ii<AcutN[Nshell].GetBlockSize(ibl1);ii++){
	// 	dEi=AcutN[Nshell].dEn[ist0+ii];
	// 	cout << " E_(i:" << ii << ") = " << dEi << endl;
	// 	for (int jj=0;jj<AcutN[Nshell].GetBlockSize(ibl2);jj++){
	// 	  dEj=AcutN[Nshell].dEn[jst0+jj];
	// 	  cout << " E_(j:" << jj << ") = " << dEj << endl;
	// 	  cout << " Delta(Ej-Ei) = " << BDelta(omegabar,(dEj-dEi),dBroad) 
	// 	       << endl;
	// 	}      
	//       }

	  cout << "opAw : " << endl << copAw << endl;
	}
      }


      // Take the trace: (use preallocated Rho2 as temp...)
      
      noalias(cRho2)=prod(cSum12,copAw);

      boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<complex<double> > > 
	cdiag(cRho2,
	     boost::numeric::ublas::range(0,cRho2.size1()),
	     boost::numeric::ublas::range(0,cRho2.size2())); 

      Trace=(sum(cdiag)).real();

    } else {

      ///////// Real operators ////////////


      // Blas matrices
      boost::numeric::ublas::matrix<double> Rho1(Nst_bl1,Nst_bl1);
      boost::numeric::ublas::matrix<double> Rho2(Nst_bl2,Nst_bl2);

      boost::numeric::ublas::matrix<double> opAw(Nst_bl1,Nst_bl2);
      boost::numeric::ublas::matrix<double> opA(Nst_bl1,Nst_bl2);

      boost::numeric::ublas::matrix<double> opB(Nst_bl2,Nst_bl1);

      boost::numeric::ublas::matrix<double> Sum12(Nst_bl2,Nst_bl1);

      Rho1=RhoN[Nshell].MatBlock2BLAS(ibl1,ibl1);
      Rho2=RhoN[Nshell].MatBlock2BLAS(ibl2,ibl2);

      opA=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2);

      if (NonDiagGF)
	//opB=Op2N[Nshell].MatBlock2BLAS(ibl2,ibl1);
	opB=trans(Op2N[Nshell].MatBlock2BLAS(ibl1,ibl2)); // Nst_bl2 x Nst_bl1 matrix (assuming real!)
      else
	opB=trans(opA);


 
      // Positive and negative energies matrix elements
      // Note for the future: sign is "-" if bosonic operators
      noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix

      //noalias(Sum12)=prod(opB,Rho1); // stronger asymmetry (negative higher)
      //noalias(Sum12)=prod(Rho2,opB); // weaker asymmetry (positive higher)

      // Calculate A_ij * Delta(w-(Ej-Ei))
      if (CalcNorm) opAw=opA; else
      noalias(opAw)=MijxBDeltaEij(opA,(&AcutN[Nshell]),ibl1,ibl2,omegabar);

      // Take the trace: (use preallocated Rho2 as temp...)


      noalias(Rho2)=prod(Sum12,opAw);

      boost::numeric::ublas::matrix_vector_range<boost::numeric::ublas::matrix<double> > 
	diag(Rho2,
	     boost::numeric::ublas::range(0,Rho2.size1()),
	     boost::numeric::ublas::range(0,Rho2.size2())); 

      Trace=sum(diag);
    }
    /////// END IF REAL OPERATORS


//     double auxtrace=0.0;
//     for(int ii=0;ii<Rho2.size1();ii++){
//       auxtrace+=Rho2(ii,ii);
//     }

    // SU(2) symmetry: <ibl1|A|ibl2> (or <ibl2|B|ibl1> )
    // Given ibl1, ibl2, calculate the following scalars
    // Get Sibl1,Sibl2
    // 
    // OpA: auxCG1[0]= sum_Szilb1 CGordan(Sibl1, Szibl1, 0.5, sigma, Sibl2, Szibl2=sigma+Szibl1)

    // OpB: auxCG1[1]= sum_Szilb1 CGordan(Sibl1, Szibl1, 0.5, sigma, Sibl2, Szibl2=sigma+Szibl1)
    //  They are the same!!
    // 2015: Need to add the possibility of < ||Sz|| > reduced!
    // In this case,
    // auxCG1[0]= sum_Szilb1 CGordan(Sibl1, Szibl1, 1.0, 0.0, Sibl2, Szibl2=Szibl1)
    // Since it is DIAGONAL in the blocks!
      // So, I need a variable in the matrix WignerEckartL=1.0 (for instance)

    double CGfactor=1.0;
    // Get S from ibl
    double Sbl1=0.0;
    double Sbl2=0.0;
    if (AcutN[Nshell].totalS){
      CGfactor=0.0;
      Sbl1=AcutN[Nshell].GetQNumber(ibl1,AcutN[Nshell].Sqnumbers[0]);
      double L1=Op1N[Nshell].WignerEckartL; // Usually 0.5
      Sbl2=AcutN[Nshell].GetQNumber(ibl2,AcutN[Nshell].Sqnumbers[0]);
      double L2=L1;
      if (NonDiagGF) L2=Op2N[Nshell].WignerEckartL;
      // only a single SU(2) for now
      for (double Szbl1=-Sbl1;Szbl1<=Sbl1;Szbl1+=1.0){
// 	double dSigma=0.5; // NOT NECESSARILY!!!!
// 	double Szbl2=Szbl1+dSigma;
// 	double auxCG=CGordan(Sbl1,Szbl1, 0.5, dSigma, Sbl2, Szbl2);
// 	CGfactor+=auxCG*auxCG;
	double dSigmaOp1=0.5; 
	if (dEqual(L1,1.0)) dSigmaOp1=0.0; // Sz for instance
	double dSigmaOp2=0.5; 
	if (dEqual(L2,1.0)) dSigmaOp2=0.0; // Sz for instance
	double Szbl2=Szbl1+dSigmaOp1;
	double auxCG1=CGordan(Sbl1,Szbl1, L1, dSigmaOp1, Sbl2, Szbl2);
	Szbl2=Szbl1+dSigmaOp2;
	double auxCG2=CGordan(Sbl1,Szbl1, L2, dSigmaOp2, Sbl2, Szbl2);
	CGfactor+=auxCG1*auxCG2;
      }
      if (dEqual(CGfactor,0.0)){cout << "Ops. CGfactor = 0.0 L1= " 
				     << L1 << " L2 = " << L2 
				     << " Sbl1 = " << Sbl1
				     << " Sbl2 = " << Sbl2
				     << endl;}
    }
    // end if totalS

    //rho_wM+=sum(diag);
    //rho_wM+=CGfactor*sum(diag);
    rho_wM+=CGfactor*Trace;
    
    if (printstuff){
      cout << " DMNRG_Spec_M: trace = " << Trace
	   << " rho_wM acc. = " << rho_wM 
	   << endl;
    }
  }
  // end loop in Op1N blocks

  return(rho_wM);

}

//////////////////

//////////////////
void CSpecFunction::DMNRG_SpecDens_ChkPHS(int Nshell){

  // Check whether particle-hole symmetry is ok
  // in the DM-NRG spectral function 

  // Loop in blocks of Operator

  bool printstuff=false;

  // Need to Check Sync!!
  bool ChkSyncRhoAcut=RhoN[Nshell].ChkSync((&AcutN[Nshell]));
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));

  if ( (!ChkSyncRhoAcut)||(!ChkSyncOp1Acut) ){
    cout << "DMNRG_SpecDens_ChkPHS: RhoN or Op1N not in sync with AcutN " << endl;
    return;
  }

  bool UseComplex=false;
  complex<double> auxC;
  if ( (RhoN[Nshell].IsComplex)&&(Op1N[Nshell].IsComplex) ){
    UseComplex=true;
  }
  // end check if complex

  double auxtracePos=0.0;
  double auxtraceNeg=0.0;

  for (int ibl1=0;ibl1<Op1N[Nshell].NumBlocks();ibl1++) {
    for (int ibl2=0;ibl2<Op1N[Nshell].NumBlocks();ibl2++) {

      int iMatbl=Op1N[Nshell].FindMatBlock(ibl1,ibl2);

      if (iMatbl>=0){ // Only <ibl1|A|ibl2>
	
	if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl2,ibl1)<0) ){
	  cout << "DMNRG_SpecDens_ChkPHS: Op2 does not connect blocks " 
	       << ibl1 << "  " << ibl2 << endl;
	  return;
	} 

	int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);  
	int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 

	// Check block sizes
	if ( (Nst_bl1!=RhoN[Nshell].GetBlockSize(ibl1))||
	     (Nst_bl2!=RhoN[Nshell].GetBlockSize(ibl2)) ){
	  cout << "DMNRG_SpecDens_ChkPHS: Block in Op1 not the same size as block in RhoN " << endl;
	  return;
	}
	if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			    (Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
	  cout << "DMNRG_SpecDens_M: Block in Op1 not the same size as block in Op2 " << endl;
	  return;
	}


	// Ok, now separating the contributions

	// Calculating the CGfactor

	// SU(2) symmetry: <ibl1|A|ibl2> (or <ibl2|B|ibl1> )
	// auxCG1= sum_Szilb1 CGordan(Sibl1, Szibl1, 0.5, sigma, Sibl2, Szibl2=sigma+Szibl1)
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

	// Calculating the traces


	double tracePartialPos=0.0;
	double tracePartialNeg=0.0;

	int nn0=AcutN[Nshell].GetBlockLimit(ibl1,0);
	int np0=AcutN[Nshell].GetBlockLimit(ibl2,0);

	if ( UseComplex ){
	  
	  //////
	  /// COMPLEX operators
	  /////


	  // Blas matrices
	  boost::numeric::ublas::matrix<complex<double> > cRho1(Nst_bl1,Nst_bl1);
	  boost::numeric::ublas::matrix<complex<double> > cRho2(Nst_bl2,Nst_bl2);

	  boost::numeric::ublas::matrix<complex<double> > copAw(Nst_bl1,Nst_bl2);
	  boost::numeric::ublas::matrix<complex<double> > copA(Nst_bl1,Nst_bl2);

	  boost::numeric::ublas::matrix<complex<double> > copB(Nst_bl2,Nst_bl1);

	  boost::numeric::ublas::matrix<complex<double> > cSum12(Nst_bl2,Nst_bl1);


	  cRho1=RhoN[Nshell].cMatBlock2BLAS(ibl1,ibl1);
	  cRho2=RhoN[Nshell].cMatBlock2BLAS(ibl2,ibl2);

	  copA=Op1N[Nshell].cMatBlock2BLAS(ibl1,ibl2);

	  if (NonDiagGF)
	    copB=Op2N[Nshell].cMatBlock2BLAS(ibl2,ibl1);
	  else
	    copB=trans(copA); // B=A^\dagger
 
	  // Positive and negative energies matrix elements
	  noalias(cSum12)=prod(copB,cRho1)+prod(cRho2,copB); // Nst2 x Nst1 matrix

	  for(int np=0;np<Nst_bl2;np++){ // sum in np (Nst_bl2)
	    double Enp=AcutN[Nshell].dEn[np0+np];
	    for(int nn=0;nn<Nst_bl1;nn++){ // sum in n (Nst_bl1)
	      double En=AcutN[Nshell].dEn[nn0+nn];
	      auxC=cSum12(np,nn)*copA(nn,np);
	      if (Enp>En)
		//tracePartialPos+=Sum12(np,nn)*opA(nn,np);
		tracePartialPos+=CGfactor*auxC.real();
	      else
		tracePartialNeg+=CGfactor*auxC.real();
	    }
	    // end loop in ibl1 states (nn)
	  }
	  // end loop in ibl2 states (np)





	} else {

	  ///////// Real operators ////////////

	  // Blas matrices
	  boost::numeric::ublas::matrix<double> Rho1(Nst_bl1,Nst_bl1);
	  boost::numeric::ublas::matrix<double> Rho2(Nst_bl2,Nst_bl2);

	  boost::numeric::ublas::matrix<double> opAw(Nst_bl1,Nst_bl2);
	  boost::numeric::ublas::matrix<double> opA(Nst_bl1,Nst_bl2);

	  boost::numeric::ublas::matrix<double> opB(Nst_bl2,Nst_bl1);

	  boost::numeric::ublas::matrix<double> Sum12(Nst_bl2,Nst_bl1);

	  // Try something different

	  Rho1=RhoN[Nshell].MatBlock2BLAS(ibl1,ibl1);
	  Rho2=RhoN[Nshell].MatBlock2BLAS(ibl2,ibl2);

	// Debuging
// 	if ( ((ibl1==3)&&(ibl2==4))|| //  N=5 1000 and 1200 st
// 	     ((ibl1==4)&&(ibl2==5)) ){
// 	if ( ((ibl1==2)&&(ibl2==3))|| // 1200 st N=4
// 	     ((ibl1==7)&&(ibl2==8)) ){
// 	if ( ((ibl1==2)&&(ibl2==3))|| // 1000 st N =4
// 	     ((ibl1==5)&&(ibl2==6)) ){
// 	  Rho1=RhoN[Nshell].MatBlock2BLAS(ibl1,ibl1);
// 	  Rho2=RhoN[Nshell].MatBlock2BLAS(ibl2,ibl2);
//  	  if (ibl1==3){
// 	  for(int ii=0;ii<Rho1.size1();ii++){
// 	    for(int jj=0;jj<Rho1.size2();jj++){
// 	      if (ii==jj) Rho1(ii,jj)=1.0;
// 	      else Rho1(ii,jj)=0.0;
// 	    }
// 	  }
// 	  // set Rho1
// 	  }
//  	  if (ibl2==5){
// 	  for(int ii=0;ii<Rho2.size1();ii++){
// 	    for(int jj=0;jj<Rho2.size2();jj++){
// 	      if (ii==jj) Rho2(ii,jj)=1.0;
// 	      else Rho2(ii,jj)=0.0;
// 	    }
// 	  }
// 	  // set Rho2
// 	  }
// 	  if (ibl1==3){
// 	    for(int ii=0;ii<Rho1.size1();ii++){
// 	      int jj=ii;
// 		if (fabs(Rho1(ii,jj))>1.E-3)
// 		  cout << "Rho1("<<ii<<","<<jj<<") = " <<Rho1(ii,jj)<<endl;
// 	    }
// 	    // end print Rho1
// 	  }
// 	  if (ibl2==5){
// 	    for(int ii=0;ii<Rho2.size1();ii++){
// 	      int jj=ii;
// 		if (fabs(Rho2(ii,jj))>1.E-3){
// 		  cout << "Rho2("<<ii<<","<<jj<<") = " <<Rho2(ii,jj)<<endl;
// 		}
// 	    }
// 	    // end print Rho2
// 	  }
// 	}
	/////////// end Debugging


	    opA=Op1N[Nshell].MatBlock2BLAS(ibl1,ibl2);

	    if (NonDiagGF)
	      opB=Op2N[Nshell].MatBlock2BLAS(ibl2,ibl1);
	    else
	      opB=trans(opA); // B=A^\dagger
 
	    // Positive and negative energies matrix elements
	    noalias(Sum12)=prod(opB,Rho1)+prod(Rho2,opB); // Nst2 x Nst1 matrix


	    for(int np=0;np<Nst_bl2;np++){ // sum in np (Nst_bl2)
	      double Enp=AcutN[Nshell].dEn[np0+np];
	      for(int nn=0;nn<Nst_bl1;nn++){ // sum in n (Nst_bl1)
		double En=AcutN[Nshell].dEn[nn0+nn];
		if (Enp>En)
		  //tracePartialPos+=Sum12(np,nn)*opA(nn,np);
		  tracePartialPos+=CGfactor*Sum12(np,nn)*opA(nn,np);
		else
		  //tracePartialNeg+=Sum12(np,nn)*opA(nn,np);
		  tracePartialNeg+=CGfactor*Sum12(np,nn)*opA(nn,np);
	      }
	      // end loop in ibl1 states (nn)
	    }
	    // end loop in ibl2 states (np)

	}
	/////// END IF REAL OPERATORS




	cout << " ibl1 = " << ibl1 
	     << " ibl2 = " << ibl2 
	     << " (Nstb1:" << Nst_bl1 << " Nstb2:" << Nst_bl2 << ")"
	     << " partial Pos: " << tracePartialPos
	     << " partial Neg: " << tracePartialNeg
	     << endl;

	auxtracePos+=tracePartialPos;
	auxtraceNeg+=tracePartialNeg;
      }
      // end if iMatBl >= 0 
    }
    // end loop in ibl2
  }
  // end loop in ibl1

  cout << " Nshell = " << Nshell
       << " Pos contribution: " << auxtracePos
       << " Neg contribution: " << auxtraceNeg
       << endl;


}
// end subroutine check phs




//////////////////
double CSpecFunction::CalcSpecDM_NRG(double omega, int UseCFS){

  // Need to get the Connecting matrix elements and energies
  // to and from a given block (say, blGS)
  // Interesting procedure...

  // Loop in blocks of Operator


  int NshOmega=CalcNfromOmega(omega);

  int Nsh0=(NshOmega-8>NshellMin?NshOmega-8:NshellMin);
  int Nsh1=(NshOmega+14<NshellMax?NshOmega+14:NshellMax);

  double SpecM=0.0;
  double AM_w=0.0;
  //for (int Nshell=NshellMin;Nshell<NshellMax;Nshell+=1){
  for (int Nshell=Nsh0;Nshell<Nsh1;Nshell+=1){

    double DN=CalcDN(Nshell);
    double sqrtDN=sqrt(DN);

    //SpecM+=DMNRG_SpecDens_M(omega,Nshell);
    //AM_w=DMNRG_SpecDens_M(omega/DN,Nshell);

    if (UseCFS==1) AM_w=CFS_SpecDens_M(omega/DN,Nshell);
    else AM_w=DMNRG_SpecDens_M(omega/DN,Nshell); // does NOT work.

    SpecM+=AM_w/DN;


    // Debug (March 2016)
//     if (fabs(fabs(omega)-0.006)<0.001)
//       cout << " NShell=" << Nshell
// 	   << " DN = " << DN
// 	   << " - A_M(w=" << omega <<") = " << AM_w 
// 	   << " - A_M(w)/DN = " << AM_w/DN 
// 	   << " - Accumulated: " << SpecM << endl;

    AM_w/=DN;
  }
  // end loop in Nshell

  cout << " omega=" << omega
       << " NshOmega=" << NshOmega
       << "  Nsh0=" << Nsh0
       << "  Nsh1=" << Nsh1
       << "  LastAw=" << AM_w
       << "  Accumulated=" << SpecM
       << endl;

  return(SpecM);

}

//////////////////
double CSpecFunction::CalcNorm(int UseCFS){

  // Need to get the Connecting matrix elements and energies
  // to and from a given block (say, blGS)
  // Interesting procedure...

  // Loop in blocks of Operator


  int Nsh0=NshellMin;
  int Nsh1=NshellMax;

  double NormN=0.0;
  double Norm=0.0;
  //for (int Nshell=NshellMin;Nshell<NshellMax;Nshell+=1){
  for (int Nshell=Nsh0;Nshell<Nsh1;Nshell+=1){

    double DN=CalcDN(Nshell);
    double sqrtDN=sqrt(DN);

    if (UseCFS==1) NormN=CFS_SpecDens_M(0.0,Nshell,true);
    else NormN=DMNRG_SpecDens_M(0.0,Nshell,true); // does NOT work.

    Norm+=NormN;

//     cout << "  Last=" << NormN
// 	 << "  Accumulated=" << Norm
// 	 << endl;

  }
  // end loop in Nshell

  return(Norm);

}

//////////////////




void CSpecFunction::CalcSpecDM_NRG_FixedOmegas(double factorWN, int UseCFS, int NwEachN){

  Omega_Rhow_Even.clear();
  Omega_Rhow_Odd.clear();
  Omega_Rhow.clear();

  // First row: omega
  Omega_Rhow_Even.push_back( vector<double> () );
  // Second row: rhow
  Omega_Rhow_Even.push_back( vector<double> () );

  Omega_Rhow_Odd.push_back( vector<double> () ); // omega
  Omega_Rhow_Odd.push_back( vector<double> () ); // rhow

  Omega_Rhow.push_back( vector<double> () ); // omegas
  Omega_Rhow.push_back( vector<double> () ); // rhows

  int Nsh0Even=(NshellMin % 2==0?NshellMin:NshellMin+1);
  int Nsh0Odd=(NshellMin % 2==0?NshellMin+1:NshellMin);

  // The right thing to do is to define Mmindisc.
  // Until then, this will do.
  if (UseCFS==1){Nsh0Even+=2;Nsh0Odd+=2;}

  cout << " CalcSpecDM_NRG_FixedOmegas: Even negative " << endl;
  cout << "  Nsh0Even = " <<  Nsh0Even << endl;
  //int iNsh=0;
  int LastNsh=0;
  // Even: negative omega
  //for (int Nsh=0;Nsh<NshellMax;Nsh+=2){
  // Again: why rho_w/DN???
  for (int Nsh=Nsh0Even;Nsh<NshellMax;Nsh+=2){
    double DN=CalcDN(Nsh);
    double WN=factorWN*DN;
//     double rho_w=0.0;
//     if (UseCFS==1) rho_w=CalcSpecDM_NRG(-WN,1);
//     else {rho_w=DMNRG_SpecDens_M(-factorWN,Nsh);rho_w/=DN;}
    // DN(Nsh+1)=DN(Nsh)*Lambda^-1/2
    // Thus, we want omega=DN(Nsh)*Lambda^-(iw/2Nw) iw=0,1,Nw-1
    for (int iomega=0;iomega<NwEachN;iomega++){
      double RedFactor=((double)iomega/(2.0*NwEachN));
      double omega=pow(Lambda,-RedFactor)*WN;
      double rho_w=0.0;
      if (omega<Gap){
	cout << " Even Neg: omega < Gap -> rho_w=0 " << endl;
	rho_w=0.0;
      } else {
	if (UseCFS==1) rho_w=CalcSpecDM_NRG(-omega,1);
	else {
	rho_w=DMNRG_SpecDens_M(-factorWN*pow(Lambda,RedFactor),Nsh);rho_w/=DN;
	}
      }
      // end if omega<Gap
//       Omega_Rhow_Even[0].push_back(-WN);
      Omega_Rhow_Even[0].push_back(-omega);
      Omega_Rhow_Even[1].push_back(rho_w);
//       Omega_Rhow[0].push_back(-WN);
      Omega_Rhow[0].push_back(-omega);
    }
    // end loop in iomega
    LastNsh=Nsh;
  }


  // Even: positive omega
  cout << " CalcSpecDM_NRG_FixedOmegas: Even positive " << endl;
  //for (int Nsh=LastNsh;Nsh>=0;Nsh-=2){
  for (int Nsh=LastNsh;Nsh>=Nsh0Even;Nsh-=2){
    double DN=CalcDN(Nsh);
    double WN=factorWN*DN;
//     //double rho_w=DMNRG_SpecDens_M(factorWN,Nsh);
//     double rho_w=0.0;
//     if (UseCFS==1) rho_w=CalcSpecDM_NRG(WN,1);
//     else {rho_w=DMNRG_SpecDens_M(factorWN,Nsh);rho_w/=DN;}
//     Omega_Rhow_Even[0].push_back(WN);
// //     Omega_Rhow_Even[1].push_back(rho_w/DN);
//     Omega_Rhow_Even[1].push_back(rho_w);

//     Omega_Rhow[0].push_back(WN);
    // DN(Nsh+1)=DN(Nsh)*Lambda^-1/2
    // Thus, we want omega=DN(Nsh)*Lambda^-(iw/2Nw) iw=0,1,Nw-1
    for (int iomega=NwEachN-1;iomega>=0;iomega--){
      double RedFactor=((double)iomega/(2.0*NwEachN));
      double omega=pow(Lambda,-RedFactor)*WN;
      double rho_w=0.0;
      if (omega<Gap){
	cout << " Even Pos: omega < Gap -> rho_w=0" << endl;
	rho_w=0.0;
      }else{
	if (UseCFS==1) rho_w=CalcSpecDM_NRG(omega,1);
	else {
	rho_w=DMNRG_SpecDens_M(factorWN*pow(Lambda,RedFactor),Nsh);rho_w/=DN;
	}
      }
      // end if omega<Gap
      Omega_Rhow_Even[0].push_back(omega);
      Omega_Rhow_Even[1].push_back(rho_w);
      Omega_Rhow[0].push_back(omega);
    }
    // end loop in iomega
    LastNsh=Nsh;
  }
  // end loops in even Nshell

  LastNsh=0;
  // Odd: negative omega
  cout << " CalcSpecDM_NRG_FixedOmegas: Odd negative " << endl;
  cout << "  Nsh0Odd = " <<  Nsh0Odd << endl;
  //for (int Nsh=1;Nsh<NshellMax;Nsh+=2){
  for (int Nsh=Nsh0Odd;Nsh<NshellMax;Nsh+=2){
    double DN=CalcDN(Nsh);
    double WN=factorWN*DN;
//     double rho_w=0.0;
//     if (UseCFS==1) rho_w=CalcSpecDM_NRG(-WN,1);
//     else {rho_w=DMNRG_SpecDens_M(-factorWN,Nsh);rho_w/=DN;}
//     Omega_Rhow_Odd[0].push_back(-WN);
//     Omega_Rhow_Odd[1].push_back(rho_w);
//     Omega_Rhow[0].push_back(-WN);
    // DN(Nsh+1)=DN(Nsh)*Lambda^-1/2
    // Thus, we want omega=DN(Nsh)*Lambda^-(iw/2Nw) iw=0,1,Nw-1
    for (int iomega=0;iomega<NwEachN;iomega++){
      double RedFactor=((double)iomega/(2.0*NwEachN));
      double omega=pow(Lambda,-RedFactor)*WN;
      double rho_w=0.0;
      if (omega<Gap){
	cout << " Odd Neg: omega < Gap -> rho_w=0" << endl;
	rho_w=0.0;
      }else{
	if (UseCFS==1) rho_w=CalcSpecDM_NRG(-omega,1);
	else {
	  rho_w=DMNRG_SpecDens_M(-factorWN*pow(Lambda,RedFactor),Nsh);rho_w/=DN;
	}
      }
      // end if omega<Gap
      Omega_Rhow_Odd[0].push_back(-omega);
      Omega_Rhow_Odd[1].push_back(rho_w);
      Omega_Rhow[0].push_back(-omega);
    }
    // end loop in iomega
    LastNsh=Nsh;
  }
  // end loop in N odd

  // Odd: positive omega
  cout << " CalcSpecDM_NRG_FixedOmegas: Odd positive " << endl;
  //for (int Nsh=LastNsh;Nsh>=1;Nsh-=2){
  for (int Nsh=LastNsh;Nsh>=Nsh0Odd;Nsh-=2){
    double DN=CalcDN(Nsh);
    double WN=factorWN*DN;
//     double rho_w=0.0;
//     if (UseCFS==1) rho_w=CalcSpecDM_NRG(WN,1);
//     else {rho_w=DMNRG_SpecDens_M(factorWN,Nsh);rho_w/=DN;}
//     Omega_Rhow_Odd[0].push_back(WN);
//     Omega_Rhow_Odd[1].push_back(rho_w);
//     Omega_Rhow[0].push_back(WN);
    // DN(Nsh+1)=DN(Nsh)*Lambda^-1/2
    // Thus, we want omega=DN(Nsh)*Lambda^-(iw/2Nw) iw=0,1,Nw-1
    for (int iomega=NwEachN-1;iomega>=0;iomega--){
      double RedFactor=((double)iomega/(2.0*NwEachN));
      double omega=pow(Lambda,-RedFactor)*WN;
      double rho_w=0.0;
      if (omega<Gap){
	cout << " Odd Pos: omega < Gap -> rho_w=0 " << endl;
	rho_w=0.0;
      }else{
	if (UseCFS==1) rho_w=CalcSpecDM_NRG(omega,1);
	else {
	rho_w=DMNRG_SpecDens_M(factorWN*pow(Lambda,RedFactor),Nsh);rho_w/=DN;
	}
      }
      // end if omega<Gap
      Omega_Rhow_Odd[0].push_back(omega);
      Omega_Rhow_Odd[1].push_back(rho_w);
      Omega_Rhow[0].push_back(omega);
    }
    // end loop in iomega
  }
  // end loops in odd Nshell

  sort(Omega_Rhow[0].begin(),Omega_Rhow[0].end());

  // Interpolate even and odd 

  int NomegasEven=Omega_Rhow_Even[0].size();
  int NomegasOdd=Omega_Rhow_Odd[0].size();

  cout << " ... done. NomegasEven= " <<  NomegasEven 
       << " NomegasOdd= " << NomegasOdd << endl;


  //Alloc accelerator
  gsl_interp_accel *accEven = gsl_interp_accel_alloc ();
  gsl_interp_accel *accOdd = gsl_interp_accel_alloc ();

  // Create Omega_Rhow: start with even interpolated 
  gsl_interp *linearinterpEven=gsl_interp_alloc(gsl_interp_linear,
					    NomegasEven);
  gsl_interp *linearinterpOdd=gsl_interp_alloc(gsl_interp_linear,
					    NomegasOdd);


  //gsl_interp_init(linearinterp,&omega_x[0],&rhow[0],Nomegas);

  cout << " Interpolating... " << endl;
  gsl_interp_init(linearinterpEven,&Omega_Rhow_Even[0][0],
		  &Omega_Rhow_Even[1][0],NomegasEven);

  gsl_interp_init(linearinterpOdd,&Omega_Rhow_Odd[0][0],
		  &Omega_Rhow_Odd[1][0],NomegasOdd);


  int iomega=0;
  int iomega2=0;
  double omega=0.0;
  double rhoauxE=0.0;
  double rhoauxO=0.0;

  // GSL error handler
  gsl_error_handler_t *old_handler;

  vector<double>::iterator dit;
  for (dit=Omega_Rhow[0].begin();
       dit<Omega_Rhow[0].end();dit++){

    omega=(*dit);

    // turns off GSL error handler
    gsl_set_error_handler_off ();
    
    rhoauxE= gsl_interp_eval (linearinterpEven, 
			      &Omega_Rhow_Even[0][0],
			      &Omega_Rhow_Even[1][0], omega, accEven);

    rhoauxO= gsl_interp_eval (linearinterpOdd, 
			      &Omega_Rhow_Odd[0][0],
			      &Omega_Rhow_Odd[1][0], omega, accOdd);

    // turns it back on
    old_handler = gsl_set_error_handler (NULL);

    // Test
//     cout << "RhoE(w="<< omega <<")= " << rhoauxE << endl;
//     cout << "RhoO(w="<< omega <<")= " << rhoauxO << endl;

    // Error handling: Check if omega is outside the Even/Odd ranges

    if ( (omega<Omega_Rhow_Even[0][0])||
	 (omega>Omega_Rhow_Even[0][NomegasEven-1]) ){

      if ( (omega<Omega_Rhow_Odd[0][0])||
	   (omega>Omega_Rhow_Odd[0][NomegasOdd-1]) ){

	// outside BOTH ranges (impossible)	
	rhoauxE=-100.0;
	rhoauxO=-100.0;

      }else{rhoauxE=rhoauxO;}
      // end if outside odd range

    }else{
      if ( (omega<Omega_Rhow_Odd[0][0])||
	   (omega>Omega_Rhow_Odd[0][NomegasOdd-1]) ){
	
	rhoauxO=rhoauxE;

      }
      // end if outside odd range
    }
    // end if outside even range

    // Test
//     cout << "RhoE(w="<< omega <<")= " << rhoauxE << endl;
//     cout << "RhoO(w="<< omega <<")= " << rhoauxO << endl;
      
    Omega_Rhow[1].push_back(0.5*(rhoauxE+rhoauxO));

  }
  // end loop in Omega_Rhow


  // Free Interpolator
  gsl_interp_free(linearinterpOdd);
  gsl_interp_free(linearinterpEven);

  // Free accellerator
  gsl_interp_accel_free (accOdd);
  gsl_interp_accel_free (accEven);


  // Discard first/last Interpolation points
//   Omega_Rhow[0].erase(Omega_Rhow[0].begin());
//   Omega_Rhow[0].pull_back();

//   Omega_Rhow[1].erase(Omega_Rhow[1].begin());
//   Omega_Rhow[1].pull_back();


}


//////////////////

void CSpecFunction::CalcSpec_ManyOmegas(int NwEachN,  double factorWN, int UseCFS){

  // Nomegas per shell!

  Omega_Rhow_Even.clear();
  Omega_Rhow_Odd.clear();
  Omega_Rhow.clear();

  Omega_Rhow.push_back( vector<double> () ); // omegas
  Omega_Rhow.push_back( vector<double> () ); // rhows

  // Even for finite-T, we are going to energies well below Mtemp.
  //int Nmax=(NshellMax>51?NshellMax-1:50);
  // Test
  int Nmax=NshellMax-1;

  int Nsh0=NshellMin;
  if (UseCFS==1){Nsh0+=2;}


  // First omega

  //   double DN=CalcDN(NshellMin);
  double DN=CalcDN(Nsh0);
  double WN=factorWN*DN;

  // For each N
  // iomega=1,NwEachN
  //double RedFactor=(double)(iomega-HalfNomegas)/(2.0*NwEachN));
  // omega=pow(Lambda,RedFactor) 

  // 
  // Negative omegas
  cout << "CalcSpec_ManyOmegas : Negative w " << endl;
  for (int Nsh=Nsh0;Nsh<=Nmax;Nsh++){
    DN=CalcDN(Nsh);
    WN=factorWN*DN;
    cout << " N = " << Nsh 
	 << " DN = " << DN 
	 << " WN = " << WN 
	 << endl;
//     for (int iomega=NwEachN;iomega>=1;iomega--){
//       double RedFactor=(double)(iomega-HalfNomegas)/(2.0*NwEachN);
    // DN(Nsh+1)=DN(Nsh)*Lambda^-1/2
    // we want DN(Nsh)*Lambda^-(iw/2*Nw) iw=0,1,Nw-1
    for (int iomega=0;iomega<NwEachN;iomega++){
      double RedFactor=(double)iomega/(2.0*NwEachN);
      double omega=pow(Lambda,-RedFactor)*WN;
      double rho_w=0.0;

      // HERE I need to add a test to see if the energy is below the gap. If it is, we add the peaks!
      if (omega<Gap){
	cout << " Neg omega: omega < Gap -> rho_w=0 " << endl;
	rho_w=0.0;
      }else{
	if (UseCFS==1) rho_w=CalcSpecDM_NRG(-omega,1);
	else{
	rho_w=DMNRG_SpecDens_M(-factorWN*pow(Lambda,RedFactor),Nsh);rho_w/=DN;
	}
      } // end if omega<Gap
      Omega_Rhow[0].push_back(-omega);
      Omega_Rhow[1].push_back(rho_w);
    }
    // end loop in iomega
  }
  // end loop in Nshells

  // Positive omegas
  cout << "CalcSpec_ManyOmegas : Positive w " << endl;
  for (int Nsh=Nmax;Nsh>=Nsh0;Nsh--){
    DN=CalcDN(Nsh);
    WN=factorWN*DN;
    cout << " N = " << Nsh 
	 << " DN = " << DN 
	 << " WN = " << WN 
	 << endl;

    for (int iomega=NwEachN-1;iomega>=0;iomega--){
      double RedFactor=(double)iomega/(2.0*NwEachN);
      double omega=pow(Lambda,-RedFactor)*WN;
      double rho_w=0.0;

      if (omega<Gap){
	cout << " Pos omega: omega < Gap -> rho_w=0 " << endl;
	rho_w=0.0;
      }else{
	if (UseCFS==1) rho_w=CalcSpecDM_NRG(omega,1);
	else{ 
	rho_w=DMNRG_SpecDens_M(factorWN*pow(Lambda,RedFactor),Nsh);rho_w/=DN;  
	}
      } // end if omega<Gap

      Omega_Rhow[0].push_back(omega);
      Omega_Rhow[1].push_back(rho_w);
    }
    // end loop in iomega

  }
  // end loop in Nshells

}
//end function
//////////////////

double CSpecFunction::CalcSpecDM_NRG_InterpolOddEven(double omega){

  double rhow=0.0;

  //Alloc accelerator
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();

  // Free accellerator
  gsl_interp_accel_free (acc);

  return(rhow);

}


//////////////////

double CSpecFunction::RhoInterpol( double omega ){

  double rhow=0.0;

  if (Omega_Rhow.size()==0) {
    cout << " Omega_Rhow not set. Returning 0" << endl;
    return(0.0);
  }
  if (Omega_Rhow[0].size()==0) {
    cout << " Omega_Rhow[0] not set. Returning 0" << endl;
    return(0.0);
  }

  int Nomegas=Omega_Rhow[0].size();

  // Alloc accelerator
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();

  // Alloc Interpolator
  gsl_interp *linearinterp=gsl_interp_alloc(gsl_interp_linear,
					    Nomegas);

  // Intialize Interpolator
  gsl_interp_init(linearinterp,&Omega_Rhow[0][0],
		  &Omega_Rhow[1][0],Nomegas);


  // Evaluate
  rhow=gsl_interp_eval (linearinterp, 
			 &Omega_Rhow[0][0],
			 &Omega_Rhow[1][0], omega, acc);


  // Free Interpolator
  gsl_interp_free(linearinterp);
  // Free accellerator
  gsl_interp_accel_free (acc);


  return(rhow);

}

//////////////////


double CSpecFunction::RhoInterpol2GSL( double omega, void *params){
  // Needs to be a static function (cannot call RhoInterpol!)
  // "static" keyword is declared in hpp only.


  double aux=0.0;
  vector < vector<double> > *pOmega_Rhow=(vector < vector<double> > *)params;

  double rhow=0.0;

  // This does NOT work!
  //rhow=RhoInterpol(omega);
  // Let's redo things...

  if (pOmega_Rhow->size()==0) {
    cout << " Omega_Rhow not set. Returning 0" << endl;
    return(0.0);
  }
  if ((*pOmega_Rhow)[0].size()==0) {
    cout << " Omega_Rhow[0] not set. Returning 0" << endl;
    return(0.0);
  }

  int Nomegas=(*pOmega_Rhow)[0].size();

  // Alloc accelerator
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();

  // Alloc Interpolator
  gsl_interp *linearinterp=gsl_interp_alloc(gsl_interp_linear,
					    Nomegas);

  // Intialize Interpolator
  gsl_interp_init(linearinterp,&(*pOmega_Rhow)[0][0],
		  &(*pOmega_Rhow)[1][0],Nomegas);


  // Evaluate
  rhow=gsl_interp_eval (linearinterp, 
			 &(*pOmega_Rhow)[0][0],
			 &(*pOmega_Rhow)[1][0], omega, acc);


  // Free Interpolator
  gsl_interp_free(linearinterp);
  // Free accellerator
  gsl_interp_accel_free (acc);


  return(rhow);

 
}

//////////////////

double CSpecFunction::KKRho( double omega ){
  // Kramers-Kronig Transformation

  if (Omega_Rhow.size()==0) {
    cout << " Omega_Rhow not set. Returning 0" << endl;
    return(0.0);
  }
  if (Omega_Rhow[0].size()==0) {
    cout << " Omega_Rhow[0] not set. Returning 0" << endl;
    return(0.0);
  }

  int Nomegas=Omega_Rhow[0].size();

  gsl_function RhowGSL;
  double *EmptyParams;

  // Define GSL function
  RhowGSL.function = &RhoInterpol2GSL;
  RhowGSL.params = &Omega_Rhow;


  // Debug
//   double result[2]={-1.0,-2.0};
//   result[0]=GSL_FN_EVAL(&RhowGSL,omega);
//   result[1]=RhoInterpol(omega);

//   cout << " Rho(" << omega <<")= " 
//        << result[0] <<" (GSL) = "
//        << result[1] <<" (Interpol) " << endl;



  double auxSum=0.0;
  for (int ii=0;ii<Nomegas-1;ii++){
    auxSum+=KKtransf(RhowGSL,Omega_Rhow[0][ii],Omega_Rhow[0][ii+1],omega);
  }

  return(auxSum);

}
//////////////////


double CSpecFunction::KKtransf( gsl_function Rhow, double x0, double x1, double omega ){


  gsl_integration_workspace *w_gsl;

  int spacesize=10000;
  double result,error,epsabs=1.0E-10,epsrel=1.0E-15;

  // Error handling
  int status=1;

  // Excluding end-point singularities
  if (omega==x0) x0+=sqrt(epsabs);
  if (omega==x1) x1-=sqrt(epsabs);

  w_gsl = gsl_integration_workspace_alloc (spacesize);

  /* result=GSL_FN_EVAL(&Rhow,omega); */

/*    gsl_integration_qawc(&Rhow,x0,x1,omega,epsabs,epsrel,spacesize,w_gsl,&result,&error); */

   gsl_set_error_handler_off();
   status = gsl_integration_qawc(&Rhow,x0,x1,omega,epsabs,epsrel,spacesize,w_gsl,&result,&error);
   
   if (status){
     while (status!=0){
       epsabs*=10.0;
       epsrel*=10.0;
       printf(" KKtransf : QWAC error status=%d. Trying epsabs=%3.1e, epsrel = %3.1e \n",status,epsabs,epsrel);  
       status = gsl_integration_qawc(&Rhow,x0,x1,omega,epsabs,epsrel,spacesize,w_gsl,&result,&error);
     }
   }
   // end if status

   //gsl_set_error_handler(NULL);
   gsl_integration_workspace_free(w_gsl);

   //   result/=pi; 
   result*=M_1_PI;
   
   return(result);

}

///////////////////
double  CSpecFunction::CalcNormInteg(){

  // Integrate the spectral function (should be close to one)


  double Norm=0.0;

  gsl_function RhowGSL;

  // Define GSL function
  RhowGSL.function = &RhoInterpol2GSL;
  RhowGSL.params = &Omega_Rhow;

  Norm=GSL_FN_EVAL(&RhowGSL,0.0);
    
  double w0=-1.0;
  double w1=1.0;

  //Allocate integration workspace
  gsl_integration_workspace *w_gsl;
  w_gsl = gsl_integration_workspace_alloc (10000);

  Norm=GSL_Integrator(&RhowGSL,w_gsl,w0, w1);

  gsl_integration_workspace_free(w_gsl);

  return(Norm);


}


//////////////////

gsl_complex CSpecFunction::GreensFunction (double omega){

  gsl_complex gf;
  double ReG=-M_PI*KKRho(omega);
  double ImG=-M_PI*RhoInterpol(omega);

  GSL_SET_COMPLEX(&gf,ReG,ImG);

  return(gf);

}

//////////////////


void CSpecFunction::SaveOmegaRhow(){


  if (Omega_Rhow.size()==0) {
    cout << " SaveOmegaRhow Error: Omega_Rhow not set." << endl;
    return;
  }
  if (Omega_Rhow[0].size()==0) {
    cout << " SaveOmegaRhow Error: Omega_Rhow[0] not set." << endl;
    return;
  }

  char Zvalue[10];
  sprintf(Zvalue,"%4.2f",z_twist);

  char filename[50];
  strcpy(filename,Name);
  //  strcat(filename,"_OmegaRhow.dat");
  strcat(filename,"_OmegaRhow");
  if (dNEqual(z_twist,1.0)){
     strcat(filename,"_zEQ");
     strcat(filename,Zvalue);
  }
  strcat(filename,".dat");

  ofstream OutFile(filename, ios::out);
  double daux;
  int Nomegas;
  vector<double>::iterator dit;

  if (!OutFile){cout << "SaveOmegaRhow: Cannot save data in " << filename << endl; return;}

  Nomegas=Omega_Rhow[0].size();

  OutFile.setf(ios::scientific,ios::floatfield);
  OutFile.precision(15);

  for (int ii=0; ii<Nomegas; ii++){
    // Format: omega rho rho rho rho 0
     OutFile << Omega_Rhow[0][ii] << "  ";
     if (ii<Omega_Rhow[1].size()){
       for (int jj=0;jj<4;jj++)
	 OutFile << Omega_Rhow[1][ii] << "  ";
     }
     OutFile << "0" << endl;
  }

  OutFile.close();


  return;

}

//////////////////

void CSpecFunction::ReadOmegaRhow(char arqextension[]){



  Omega_Rhow.clear();
  Omega_Rhow.push_back( vector<double> () ); // omegas
  Omega_Rhow.push_back( vector<double> () ); // rhows

  char filename[50];
  strcpy(filename,Name);
  //  strcat(filename,"_OmegaRhow.dat");
  // arqextension is "_OmegaRhow.dat" by default
  strcat(filename,arqextension);


  // instream
  ifstream InFile;

  double aux[4];
  double omega_val;
  int iaux;

  InFile.open(filename);
  if (InFile.is_open()){
    while (!InFile.eof()){
      InFile >> omega_val >> aux[0] >> aux[1] >> aux[2] >> aux[3] >> iaux;
      Omega_Rhow[0].push_back(omega_val);
      Omega_Rhow[1].push_back(aux[3]);
    }
    Omega_Rhow[0].pop_back();
    Omega_Rhow[1].pop_back();
  }
  else cout << "Cant open " << filename << endl;
  InFile.close();

  return;

}

//////////////////

void CSpecFunction::PrintOmegaRhow(){

  if (Omega_Rhow.size()==0) {
    cout << " PrintOmegaRhow Error: Omega_Rhow not set." << endl;
    return;
  }
  if (Omega_Rhow[0].size()==0) {
    cout << " PrintOmegaRhow Error: Omega_Rhow[0] not set." << endl;
    return;
  }

  int NomegasTot=Omega_Rhow[0].size();
  for (int ii=0;ii<NomegasTot;ii++){
    cout << Name << "(omega= "<< Omega_Rhow[0][ii] << " )= ";
    if (ii<Omega_Rhow[1].size()){
      cout << Omega_Rhow[1][ii];
    }
    cout << endl;
  }
  // end loop in Omega_Rhow

  return;
}

//////////////////


void CSpecFunction::PrintOmegaEven(){

  if (Omega_Rhow_Even.size()==0) {
    cout << " PrintOmegaEven Error: Omega_Rhow_Even not set." << endl;
    return;
  }
  if (Omega_Rhow_Even[0].size()==0) {
    cout << " PrintOmegaEven Error: Omega_Rhow_Even[0] not set." << endl;
    return;
  }

  int NomegasTot=Omega_Rhow_Even[0].size();
  for (int ii=0;ii<NomegasTot;ii++){
    cout << Name << "_even(omega= "<< Omega_Rhow_Even[0][ii] << " )= ";
    if (ii<Omega_Rhow_Even[1].size()){
      cout << Omega_Rhow_Even[1][ii];
    }
    cout << endl;
  }
  // end loop in Omega_Rhow_Even

  return;
}

//////////////////


//////////////////
void CSpecFunction::GetSubGapData(int Nshell, 
		     vector< double > &Eb, 
		     vector< double > &wb,
		     bool NegOmega){

//   bool printstuff=true;
  bool printstuff=false;

  double DN=CalcDN(Nshell);
  double GapDN=Gap/DN;


  // Need to Check Sync!!
  bool ChkSyncOp1Acut=Op1N[Nshell].ChkSync((&AcutN[Nshell]));

  if ( (!ChkSyncOp1Acut) ){
    cout << "GetSubGapData: Op1N not in sync with AcutN " << endl;
    return;
  }

  for (int iMatbl=0;iMatbl<Op1N[Nshell].NumMatBlocks();iMatbl++) {
    int ibl1=0;
    int ibl2=0;

    Op1N[Nshell].GetBlocksFromMatBlock(iMatbl,ibl1,ibl2);

    if ( (NonDiagGF)&&(Op2N[Nshell].FindMatBlock(ibl1,ibl2)<0) ){
      cout << "GetSubGapData: Op2 does not connect blocks " 
	   << ibl2 << "  " << ibl1 << endl;
      return;
    } 

    int Nst_bl1=Op1N[Nshell].GetBlockSize(ibl1);  
    int Nst_bl2=Op1N[Nshell].GetBlockSize(ibl2); 

    // Check block sizes
    if ( (NonDiagGF)&&( (Nst_bl1!=Op2N[Nshell].GetBlockSize(ibl1))||
			(Nst_bl2!=Op2N[Nshell].GetBlockSize(ibl2)) ) ){
      cout << "GetSubGapData: Block in Op1 not the same size as block in Op2 " << endl;
      return;
    }

    int ist0=AcutN[Nshell].GetBlockLimit(ibl1,0);
    int jst0=AcutN[Nshell].GetBlockLimit(ibl2,0);

    double Ei=0.0;
    double Ej=0.0;

    double auxM=0.0;
    complex<double> cauxM=ZeroC;


    double ExpEi=1.0;
    double ExpEj=1.0;
    double betabar=1.0e10;

    // SU(2) symmetry
    double CGfactor=1.0;
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


    for (int ii=0;ii<Nst_bl1;ii++){
      Ei=AcutN[Nshell].dEn[ist0+ii];
      ExpEi=exp(-betabar*Ei);
      for (int jj=0;jj<Nst_bl2;jj++){
	Ej=AcutN[Nshell].dEn[jst0+jj];
	ExpEj=exp(-betabar*Ej);

	if (Op1N[Nshell].IsComplex){
	  cauxM=Op1N[Nshell].cGetMatEl(ii,jst0+jj);
	  if (NonDiagGF)
	    cauxM*=conj(Op2N[Nshell].cGetMatEl(ist0+ii,jst0+jj));
	  else cauxM*=conj(cauxM);
	  auxM=cauxM.real(); // Only takes the real part! 
	  // Watch out for NonDiagGF AND IsComplex
	} else {
	  auxM=Op1N[Nshell].GetMatEl(ist0+ii,jst0+jj);
	  if (NonDiagGF)
	    auxM*=Op2N[Nshell].GetMatEl(ist0+ii,jst0+jj);
	  // Is this the right thing to do?
	  // Yes, because I always have B and not B^dagger
	  else{auxM*=auxM;}
	}
	// end if Op1 is complex

	if ( (ExpEi>0.5)||(ExpEj>0.5) ){
	  // found a GS
	  if (fabs(Ej-Ei)<GapDN){
	    if (ExpEi>0.5){
	      // Ei is the GS 
	      // |<GS|Op1|Ej>|^2 delta(w-(Ej-0)) -> Positive omega
	      // It is ALWAYS <Ei|Op1|Ej>
	      // since we are looping in the connecting blocks 
	      // in series of 3 to make a sharp peak
	      if ( (!NegOmega)&&(dNEqual(auxM,0.0)) ){
		Eb.push_back((Ej-Ei)*DN);Eb.push_back((Ej-Ei)*DN);Eb.push_back((Ej-Ei)*DN);
		// Don't forget the 1/DN factor
		//Msq_0i.push_back(auxM/DN);
		// SU(2) symmetry
		wb.push_back(0.0);wb.push_back(CGfactor*auxM);wb.push_back(0.0);
	      }
	      // end if positive omega
	      if (printstuff){
		if (!NegOmega){
		  cout << " ibl1 (GS) = " << ibl1
		       << " ibl2  = " << ibl2
		       << " ist (GS) = " << ii 
		       << " jst = " << jj 
		       << " E_i (GS) = " << Ei
		       << " E_j = " << Ej
		       << " (E_j-E_i)*DN = " << (Ej-Ei)*DN
		       << " Op1_0j = " << auxM << " similar to " << Op1N[Nshell].GetMatEl(ii,jj)
		       << " Op1_0j(complex) = " << cauxM
		       << endl;
		}else{
		  cout << " ibl1 (GS) = " << ibl1
		       << " ibl2  = " << ibl2
		       << " ist (GS) = " << ii 
		       << " jst = " << jj 
		       << " E_i (GS) = " << Ei
		       << " E_j = " << Ej
		       << " (E_j-E_i)*DN = " << (Ej-Ei)*DN
		       << " is POSITIVE so does not enter Op1_0j " 
		       << endl;
		}
		// end if NegOmega
	      }
	      // end if printstuff
	    }
	    // end if ExpEI>0.5  
	    else {
	      // Ej is the GS 
	      // |<Ei|Op1|GS>|^2 delta(w+Ei)) -> Negative omega
	      // It is ALWAYS <Ei|Op1|Ej>
	      // since we are looping in the connecting blocks 
	      // in series of 3 to make a sharp peak
	      if ( (NegOmega)&&(dNEqual(auxM,0.0)) ){
		Eb.push_back((Ej-Ei)*DN);Eb.push_back((Ej-Ei)*DN);Eb.push_back((Ej-Ei)*DN);
		// Don't forget the 1/DN factor
		//Msq_0i.push_back(auxM/DN);
		// SU(2) symmetry
		wb.push_back(0.0);wb.push_back(CGfactor*auxM);wb.push_back(0.0);
	      }
	      // end if negative omega
	      if (printstuff){
		if (!NegOmega){
		  cout << " ibl1  = " << ibl1
		       << " ibl2 (GS) = " << ibl2
		       << " ist = " << ii 
		       << " jst (GS) = " << jj 
		       << " E_i = " << Ei
		       << " E_j (GS) = " << Ej
		       << " (E_j-E_i)*DN = " << (Ej-Ei)*DN
		       << " is NEGATIVE so does not enter here. " 
		       << endl;

		}else{
		  cout << " ibl1  = " << ibl1
		       << " ibl2 (GS) = " << ibl2
		       << " ist = " << ii 
		       << " jst (GS) = " << jj 
		       << " E_i = " << Ei
		       << " E_j (GS) = " << Ej
		       << " (E_j-E_i)*DN = " << (Ej-Ei)*DN
		       << " Op1_i0 = " << auxM << " similar to " << Op1N[Nshell].GetMatEl(ii,jj)
		       << " Op1_i0(complex) = " << cauxM
		       << endl;
		}
		// end if NegOmega
	      }
	      // end if printstuff

	    }
	    // end if ExpEJ>0.5
	  }
	  // end if inside the gap 
	}
	// end if found a GS
      }
      // end loop in ibl2
    }
    // end loop in ibl1

  }
  // end loop in MatBlocks
  /// Loop in blocks


  return;
}
//////////////////
