#include <iostream>
#include <fstream>


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>
#include <cmath>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "TwoChQSz.hpp"

#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif



bool AeqB(CNRGbasisarray *pAeigCut, int ibl1,int ibl2){

  return (ibl1==ibl2);

}

double AtimesB(CNRGbasisarray *pAeig, CNRGbasisarray *pSingleSite, 
	       int ibl1,int ibl2){

  return((double)ibl1*ibl2);

}




int main (){

  CNRGarray Aeig(2);

  CNRGbasisarray AeigCut(2);

  CNRGbasisarray Abasis(2);

  CNRGbasisarray SingleSite(2);

  // STL vector

  CNRGmatrix Qm1fNQ[2];

  CNRGmatrix MQQp1;

  double U,ed,Gamma1,Gamma2;
  vector<double> Params;
  double Lambda;
  double HalfLambdaFactor;
  double Dband=1.0;
  int calcdens;

  double DN=0.0;
  double TM=0.0;
  double betabar=0.727;
  double Temp=0.0;
  double Sus=0.0;
  //vector<double> SuscepChain;
  char arqSus[32],arqname[32];



  double chi_m1,chi_N[2];
  double daux[4];

  int Nsites,Nsitesmax=2;

  int Ncutoff=700;
  int UpdateBefCut=0;
  int auxIn;

  // outstream
  ofstream OutFile;
  // instream
  ifstream InFile;


  int ii,jj,i1,i2;

  // STL iterator:

  vector<double>::iterator diter;

  // Test
  CNRGmatrix HN;

  HN.CheckForMatEl=AeqB;
  HN.CalcMatEl=AtimesB;
  if (HN.CheckForMatEl(&Abasis,2,2)) cout << "TRUE!" << endl;
  else cout << "FALSE!" << endl;
  cout << HN.CalcMatEl(&Abasis,&SingleSite,3,4) << endl;
  exit(0);


  ///            ///
  /// Begin code ///
  ///            ///
  U=0.5;
  ed=-0.5*U;
  Gamma1=0.0282691;
  Gamma2=0.0;
  Lambda=2.5;
  
  // Steps: Input parameters
  InFile.open("nrg_input_TwoChAnderson.dat");
  if (InFile.is_open())
    {
      InFile >> Nsitesmax;
      InFile >> Ncutoff;
      InFile >> U;
      InFile >> Gamma1;
      InFile >> Gamma2;
      InFile >> ed;
      InFile >> Lambda;
      InFile >> Dband;
      InFile >> auxIn;
      InFile >> UpdateBefCut;
      InFile >> calcdens;
    }
  else cout << "can't open nrg_input_TwoChAnderson.dat" << endl;

  InFile.close();


  ///////////////////////////
  ///////////////////////////

  if (calcdens==2)
    {
      //SuscepChain.clear();
      strcpy(arqSus,"SuscepImp_25_726.dat");
  
      strcpy(arqname,"SuscepChain.dat");
      InFile.open(arqname);
      if (InFile.fail())
	{
	  cout << "Can't find " << arqname << endl;
	  strcpy(arqSus,"SuscepChain.dat");
	  calcdens=3;
	}
///////////////////////
//       else
// 	{
// 	  while (!InFile.eof())
// 	    {
// 	      InFile >> Temp >> daux[0] >> daux[0] >> Sus;
// 	      SuscepChain.push_back(Sus);
// 	    }
// 	  SuscepChain.pop_back();
// 	  if (SuscepChain.size()<Nsitesmax-1)
// 	    {
// 	      cout << "Nsuscep = " << SuscepChain.size() << " Nsitesmax = " << Nsitesmax << endl;
// 	      cout << " Longer SuscepChain needed! Doing it all over..." << endl;
// 	      SuscepChain.clear();
// 	      strcpy(arqSus,"SuscepChain.dat");
// 	      calcdens=3;
// 	    } 
// 	}
//       InFile.close();
//////////////////////////
     }
  if (calcdens==3)
    {
      strcpy(arqSus,"SuscepChain.dat");
      U=0.0;ed=0.0;Gamma1=0.0;Gamma2=0.0;
    }

  ///////////////////////////
  ///////////////////////////

  HalfLambdaFactor=0.5*(1.0+(1.0/Lambda));
  chi_m1=sqrt(2.0*Gamma1/pi)/(sqrt(Lambda)*HalfLambdaFactor);

  double U_tilde=0.5*U/HalfLambdaFactor;
  double ed_tilde=ed/HalfLambdaFactor;
  double Gamma1_tilde=Gamma1*(2.0/pi)/(HalfLambdaFactor*HalfLambdaFactor);
  double Gamma2_tilde=Gamma2*(2.0/pi)/(HalfLambdaFactor*HalfLambdaFactor);


  OutFile.open("NRG_in.txt");
  OutFile << "Begin NRG 2-ch calculation" << endl;
  OutFile << " Nsitesmax    = " << Nsitesmax-1 << endl;
  OutFile << " Ncutoff      = " << Ncutoff << endl;
  OutFile << " U            = " << U << endl;
  OutFile << " Gamma1        = " << Gamma1 << endl;
  OutFile << " Gamma2        = " << Gamma2 << endl;
  OutFile << " ed           = " << ed << endl;
  OutFile << " Lambda       = " << Lambda << endl;
  OutFile << " Dband        = " << Dband << endl;
  OutFile << " UpdateBefCut = " << UpdateBefCut << endl;
  OutFile << " calcdens     = " << calcdens << endl;
  OutFile << "=================================" << endl;
  OutFile << "U~ = " << U_tilde << endl;
  OutFile << "ed~ = " <<  ed_tilde << endl;
  OutFile << "Gamma1~ = " <<  Gamma1_tilde<< endl;
  OutFile << "Gamma2~ = " <<  Gamma2_tilde<< endl;
  OutFile << "=================================" << endl;
  OutFile.close();



  // Define H0 (impurity + 1s site)
  // Output:
  //       Aeig,
  //       fd_{1sigma} and fd_{2sigma} matrix elements
  //

  // Set single site (use pointers in the subroutines!!)

  TwoChQSz_SetSingleSite(&SingleSite);


  Params.push_back(U_tilde/sqrt(Lambda));
  Params.push_back(ed_tilde/sqrt(Lambda));
  Params.push_back(sqrt(Gamma1_tilde/Lambda));
  Params.push_back(sqrt(Gamma2_tilde/Lambda));

  Nsites=0;
  DN=HalfLambdaFactor*pow(Lambda,(-(Nsites-1)/2.0) );
  TM=DN/betabar;
  cout << "DN = " << DN << "TM = " << TM << endl;

  //TwoChQSz_SetH0Anderson(Params,&SingleSite,&Aeig,Qm1fNQ);
  TwoChQSz_SetH0Anderson(Params,&SingleSite,&Aeig,&Abasis);

  Aeig.PrintEn();

  // Calculate Susceptibility

  if ( (calcdens==2)||(calcdens==3) )
    {
      TM=DN/betabar;
      Params.clear();
      Params.push_back(betabar);
      Sus=CalcSuscep(Params,&Aeig,1,false);
      double Sus0=0.0;
      int nlines=0;
      if (calcdens==2)
	{
	  InFile.open("SuscepChain.dat");
	  InFile.clear(); 
	  InFile.seekg(0, ios::beg); // rewind
	  if (InFile.fail())
	    {
	      cout << "Can't open  SuscepChain.dat"<< endl;
	    }
	  else
	    {
	      while ( (!InFile.eof())&&(nlines<=Nsites) ) 
		{
		  InFile >> daux[0] >> daux[0] >> daux[0] >> Sus0;
		  nlines++;
		}
	    }
	  InFile.close();
	}
      //if (calcdens==2) Sus0=SuscepChain[Nsites];
      else Sus-=1.0/8.0; // exclude dot site in the chain.
      if (Nsites==0) OutFile.open(arqSus);
      else OutFile.open(arqSus,ofstream::app);
      OutFile.precision(20);
      OutFile << scientific << TM << " " << Sus << " " << Sus0 << " " << Sus-Sus0 << endl;
      OutFile.close();

    }

  TwoChQSz_UpdateQm1fQ(&SingleSite,&Aeig,&Abasis,Qm1fNQ);

  // Check Matrix
  for (int ich=1;ich<=2;ich++)
    {
      cout << " Printing Qm1fQ channel: " << ich << endl;
      for (int ibl=0; ibl<Qm1fNQ[ich-1].NumMatBlocks();ibl++)
	{
	  Qm1fNQ[ich-1].PrintMatBlock(ibl);
	}
    }

  // Loop on Nsites: start from Nsites=0 (imp+1)
  //

  Nsites=1;
  while (Nsites<=Nsitesmax)
    {

      DN=HalfLambdaFactor*pow(Lambda,(-(Nsites-1)/2.0) );
      TM=DN/betabar;
      cout << "DN = " << DN << "TM = " << TM << endl;

      // 0 - Update chi_N, eps_N

      daux[0]=(double)( 1.0-pow(Lambda,(-Nsites)) );
      daux[1]=(double)sqrt( 1.0-pow(Lambda,-(2*Nsites-1)) );
      daux[2]=(double)sqrt( 1.0-pow(Lambda,-(2*Nsites+1)) );  
      daux[3]=0.5*(1.0+(1.0/Lambda))*(double)sqrt(Lambda);

      chi_N[0]=daux[0]/(daux[1]*daux[2]);
      chi_N[1]=daux[0]/(daux[1]*daux[2]);

      cout << "chi_N = " << chi_N[0] << "  " << chi_N[1] << endl;


      cout << "Nsites = " << Nsites << endl;
      cout << "BEG Eig Nshell = " << Aeig.Nshell << endl;

      // 1 - Eliminate states and Build Abasis

      cout << "Cutting states..." << endl;
      AeigCut.CNRGbasisarray::ClearAll();
      AeigCut=CutStates(&Aeig, Ncutoff);

      cout << "... done cutting states." << endl;

      // Calculate new matrix elements using the CUT basis:
      // update Qm1fNQ, Qm1cdQ, etc.

      if ( (UpdateBefCut==0)&&(Nsites>1) )
	{
	  cout << "Updating matrices after cutting... " << endl;    
	  TwoChQSz_UpdateMatrixAfterCutting(&SingleSite,
					    &AeigCut, &Abasis, Qm1fNQ, &MQQp1);
	  cout << "... done updating matrices. " << endl;
	}



      //Q1Q2Sz_BuildBasis(&AeigCut,&Abasis,&SingleSite,UpdateBefCut); 
      // doesnt work
      QSz_BuildBasis(&AeigCut,&Abasis,&SingleSite,UpdateBefCut);


      cout << "No blocks = " << Abasis.NumBlocks() << endl;
      
      cout << "No states = " << Abasis.Nstates() << endl;

      //Abasis.PrintAll();


      // 2 - Build and diagonalize H_N+1

      Params.clear();
      Params.push_back(chi_N[0]);
      Params.push_back(chi_N[1]);
      Params.push_back(Lambda);

      cout << "Diagonalizing HN... " << endl;    

      TwoChQSz_DiagHN(Params,&Abasis,&SingleSite,Qm1fNQ,&Aeig);

      Aeig.PrintEn();

      cout << "..done diagonalizing HN. " << endl;    

      // 3 - Update Qm1f1NQ, Qm1f2Q, Qm1cdQ, etc.

      if ( (UpdateBefCut==1)&&(Nsites<Nsitesmax) )
	{
	  cout << "Updating matrices before cutting... " << endl;
	  TwoChQSz_UpdateQm1fQ(&SingleSite,&Aeig,&Abasis,Qm1fNQ);
	  cout << "... done updating matrices. " << endl;

	}


      // Calculate Susceptibility

      if ( (calcdens==2)||(calcdens==3) )
	{
	  TM=DN/betabar;
	  Params.clear();
	  Params.push_back(betabar);
	  Sus=CalcSuscep(Params,&Aeig,1,false);
	  double Sus0=0.0;
	  int nlines=0;
	  // Gets SuscepChain iteration by iteration
	  if (calcdens==2)
	    {
	      InFile.open("SuscepChain.dat");
	      InFile.clear(); 
	      InFile.seekg(0, ios::beg); // rewind
	      if (InFile.fail())
		{
		  cout << "Can't find " << arqname << endl;
		}
	      else
		{
		  while ( (!InFile.eof())&&(nlines<=Nsites) ) 
		    {
		      InFile >> daux[0] >> daux[0] >> daux[0] >> Sus0;
		      nlines++;
		    }
		}
	      InFile.close();
	    }
	  //if (calcdens==2) Sus0=SuscepChain[Nsites];
	  else Sus-=1.0/8.0; // exclude dot site in the chain.

	  if (Nsites==0) OutFile.open(arqSus);
	  else OutFile.open(arqSus,ofstream::app);
	  OutFile.precision(20);
	  OutFile << scientific << TM << " " << Sus << " " << Sus0 << " " << Sus-Sus0 << endl;
	  OutFile.close();
	}

      // 4 - Update Nsites

      Nsites++;

    }
 
  cout << "=== Calculation Finished! ==== "<< endl;
  OutFile.open("NRG_end.txt");
  OutFile << "END NRG calculation" << endl;
  OutFile.close();
  
  cout << "Calling destructors " << endl;

}
//END code


