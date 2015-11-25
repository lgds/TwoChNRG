#include <iostream>
#include <fstream>


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>
#include <cmath>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "TwoChQS.hpp"

#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif




int main (int argc, char* argv[]){

  //char ModelOption[]="Anderson";
  //int ModelNo=0;

  // command-line model options
#include"ModelOpt.cpp"

  CNRGarray Aeig(2);

  CNRGbasisarray AeigCut(2);

  CNRGbasisarray Abasis(2);

  CNRGbasisarray SingleSite(3);

  CNRGmatrix* MatArray;
  int NumNRGarrays=2;
  // MatArray 0 is nd
  // MatArray 1 is cd


  // Thermodynamics
  CNRGthermo Suscep;
  CNRGthermo Entropy;


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

  int Nsites=0,Nsitesmax=2;

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

  // Phonon
  // Test!! Remove later
  double Nph=1.0;
  double w0=0.2;
  double lambda=0.4;
  double alpha=0.4;
  double tp=0.2;


  ///            ///
  /// Begin code ///
  ///            ///
  U=0.5;
  ed=-0.5*U;
  Gamma1=0.0282691;
  Gamma2=0.0;
  Lambda=2.5;
  
  // Steps: Input parameters
  InFile.open("nrg_input_TwoCh.dat");
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
  else cout << "can't open nrg_input_TwoCh.dat" << endl;

  InFile.close();

  if (ModelNo==4)
    {
      InFile.open("Input_Phonon.dat");
      if (InFile.is_open())
	{
	  InFile >> Nph;
	  InFile >> w0;
	  InFile >> lambda;
	  InFile >> tp;
	  InFile >> alpha;
	}
      else cout << "can't open Input_Phonon.dat" << endl;
      InFile.close();
      Gamma1=pi*tp*tp;
      Gamma2=alpha*alpha*Gamma1;
    }
  // end if phonon


  ///////////////////////////
  ///////////////////////////

  // New stuff
   strcpy(Suscep.ArqName,"SuscepImp2Ch_25_726.dat");
   strcpy(Suscep.ChainArqName,"SuscepChain2Ch.dat");
   Suscep.Calc=CalcSuscep;
  
   strcpy(Entropy.ArqName,"EntropyImp2Ch_25_726.dat");
   strcpy(Entropy.ChainArqName,"EntropyChain2Ch.dat");
   Entropy.Calc=CalcEntropy;

   // Chain model
   if(ModelNo==3) calcdens=3;

   if (calcdens==2)
     {
       Suscep.CalcChain=false;
       Entropy.CalcChain=false;
       Suscep.Nsite0Chain=0; // Anderson chains
       Entropy.Nsite0Chain=0; // Anderson chains
       
       if ( (!FileExists(Suscep.ChainArqName))||
	    (!FileExists(Entropy.ChainArqName)) )
	 {
	   cout << " Can't find chain files: " << endl 
		<< "   " << Suscep.ChainArqName << endl
		<< "   " << Entropy.ChainArqName << endl;
	   cout << " Exiting... " << endl;
	   exit(0);
	 }
       else
	 cout << " Found files " 
	      << Suscep.ChainArqName << ", " 
	      << Entropy.ChainArqName 
	      << endl;
     }
   if (calcdens==3)
     {
       //strcpy(arqSus,"SuscepChain.dat");
       Suscep.CalcChain=true;
       Entropy.CalcChain=true;
       if ( (FileExists(Suscep.ChainArqName))||
	    (FileExists(Entropy.ChainArqName)) )
	 {
	   cout << " Error: Chain files already exist: " << endl
		<< "   " << Suscep.ChainArqName << endl
		<< "   " << Entropy.ChainArqName << endl;
	   cout << " Not rewriting. Exiting... " << endl;
	   exit(0);
	 }
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
  OutFile << " Model : " << ModelOption << endl;
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
  if (ModelNo==4)
    {
      OutFile << "Nph = " <<  Nph << endl;
      OutFile << "w0 = " <<  w0 << endl;
      OutFile << "lambda = " <<  lambda << endl;
      OutFile << "tp = " <<  tp << endl;
      OutFile << "alpha = " <<  alpha << endl;
      OutFile << "chi_S_tilde = " <<  sqrt(Gamma1_tilde) << endl;
      OutFile << "chi_A_tilde = " <<  alpha*sqrt(Gamma1_tilde) << endl;
      OutFile << "=================================" << endl;
    }
  OutFile.close();



  // Define H0 (impurity + 1s site)
  // Output:
  //       Aeig,
  //       fd_{1sigma} and fd_{2sigma} matrix elements
  //

  // Set single site (use pointers in the subroutines!!)

  TwoChQS_SetSingleSite(&SingleSite);

  Params.push_back(U_tilde/sqrt(Lambda));
  Params.push_back(ed_tilde/sqrt(Lambda));
  Params.push_back(sqrt(Gamma1_tilde/Lambda));
  Params.push_back(sqrt(Gamma2_tilde/Lambda));

  Nsites=0;
  DN=HalfLambdaFactor*pow(Lambda,(-(Nsites-1)/2.0) );
  TM=DN/betabar;
  cout << "DN = " << DN << "  TM = " << TM << endl;



  switch (ModelNo)
    {
    case 0 :
      TwoChQS_SetH0Anderson(Params,&SingleSite,&Aeig,&Abasis);
      Suscep.dImpValue=1.0/8.0;
      Entropy.dImpValue=2.0*log(2.0);
      break;

    case 1 :
      cout << " Testing 2ch Kondo model. " << endl;
      Params.clear();
      Params.push_back(Gamma1);
      Params.push_back(Gamma2);
      TwoChQS_SetH0Kondo(Params,&SingleSite,&Aeig,&Abasis,MatArray);
      Suscep.dImpValue=1.0/4.0;
      Entropy.dImpValue=log(2.0);
      break;

    case 3 :
      TwoChQS_SetH0Chain(Params,&SingleSite,&Aeig,&Abasis);
      Suscep.dImpValue=0.0;
      Entropy.dImpValue=0.0; // Check this!!
      Suscep.CalcChain=true;
      Entropy.CalcChain=true;
      break;
    case 4 :
      Params.push_back(w0); //
      Params.push_back(lambda); //
      Params.push_back(alpha); //
      Params.push_back(Nph); //
      // Using the untransformed Hamiltonian
      //TwoChQS_SetH0CMphononwTransf(Params,&SingleSite,&Aeig,&Abasis,MatArray);
      TwoChQS_SetH0CMphonon(Params,&SingleSite,&Aeig,&Abasis,MatArray);
      // Qm1fQ needs special updates
      for (int ich=0;ich<=1;ich++)
	{
	  Qm1fNQ[ich].NeedOld=false;
	  Qm1fNQ[ich].CheckForMatEl=TwoChQS_cd_check;
	}

      Qm1fNQ[0].CalcMatEl=TwoChQS_cd_ich1_Phonon_MatEl;
      Qm1fNQ[1].CalcMatEl=TwoChQS_cd_ich2_Phonon_MatEl;
      /////////
      Suscep.dImpValue=0.0;
      Entropy.dImpValue=0.0; // Check this!!
      strcpy(Suscep.ArqName,"SuscepImp2Ch_Phonon.dat");
      strcpy(Entropy.ArqName,"EntropyImp2Ch_Phonon.dat");
      // Update with a new scheme
      cout << " Updating! " << endl;
      AeigCut.CNRGbasisarray::ClearAll();
      AeigCut=CutStates(&Aeig, Ncutoff);
      UpdateMatrices(&SingleSite,&AeigCut, 
		     &Abasis,&Qm1fNQ[0], 2);

      break;
    default :
      cout << " Model not implemented. Exiting... " << endl;
      exit(0);
      break;
    }
  //end switch models

  Aeig.PrintEn();

  if (ModelNo!=4){TwoChQS_UpdateQm1fQ(&SingleSite,&Aeig,&Abasis,Qm1fNQ);}

  cout << " f0_ch1 : " << endl;
  //Qm1fNQ[0].PrintAllBlocks();
  Qm1fNQ[0].PrintMatBlock(1,3);
  Qm1fNQ[0].PrintMatBlock(2,3);
  cout << " f0_ch1 : " << endl;
  //Qm1fNQ[1].PrintAllBlocks();
  Qm1fNQ[1].PrintMatBlock(1,3);
  Qm1fNQ[1].PrintMatBlock(2,3);


  // Calculate Susceptibility at Nsites=0!

  if ( (calcdens==2)||(calcdens==3) )
    {
      TM=DN/betabar;
      Params.clear();
      Params.push_back(betabar);

      Suscep.ReadNChainValue(Nsites,0);
      Suscep.AddValue(Params,&Aeig,1,true,TM);
      Suscep.SaveNValue(Nsites,0);
      
      Entropy.ReadNChainValue(Nsites,0);
      Entropy.AddValue(Params,&Aeig,1,true,TM);
      Entropy.SaveNValue(Nsites,0);

    }
  //end if calcdens=2 or 3

  //TwoChQS_UpdateQm1fQ(&SingleSite,&Aeig,&Abasis,Qm1fNQ);



  // Check Matrix
//    for (int ich=1;ich<=2;ich++)
//      {
//        cout << " Printing Qm1fQ channel: " << ich << endl;
//        for (int ibl=0; ibl<Qm1fNQ[ich-1].NumMatBlocks();ibl++)
//  	{
//  	  Qm1fNQ[ich-1].PrintMatBlock(ibl);
//  	}
//      }

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

      cout << " Nstates: ";
      cout <<  Aeig.Nstates() << endl;


      cout << "Cutting states..." << endl;
      AeigCut.CNRGbasisarray::ClearAll();
      AeigCut=CutStates(&Aeig, Ncutoff);

      cout << "... done cutting states." << endl;

      // Calculate new matrix elements using the CUT basis:
      // update Qm1fNQ, Qm1cdQ, etc.

      if ( (UpdateBefCut==0)&&(Nsites>1) )
	{
	  cout << "Updating matrices after cutting... " << endl;    
	  TwoChQS_UpdateMatrixAfterCutting(&SingleSite,
	  				    &AeigCut, &Abasis, Qm1fNQ, &MQQp1);
	  cout << "... done updating matrices. " << endl;
	}

      // Spec density calculation, other similar things woudl go HERE

      // This is where the loop really begins...

      QS_BuildBasis(&AeigCut,&Abasis,&SingleSite,UpdateBefCut);


      cout << "No blocks = " << Abasis.NumBlocks() << endl;
      
      cout << "No states = " << Abasis.Nstates() << endl;

      //Abasis.PrintAll();


      // 2 - Build and diagonalize H_N+1

      Params.clear();
      Params.push_back(chi_N[0]);
      Params.push_back(chi_N[1]);
      Params.push_back(Lambda);

      cout << "Diagonalizing HN... " << endl;    

      TwoChQS_DiagHN(Params,&Abasis,&SingleSite,Qm1fNQ,&Aeig);

      Aeig.PrintEn();

      cout << "..done diagonalizing HN. " << endl;    

      // 3 - Update Qm1f1NQ, Qm1f2Q, Qm1cdQ, etc.

      if ( (UpdateBefCut==1)&&(Nsites<Nsitesmax) )
	{
	  cout << "Updating matrices before cutting... " << endl;
	  TwoChQS_UpdateQm1fQ(&SingleSite,&Aeig,&Abasis,Qm1fNQ);
	  cout << "... done updating matrices. " << endl;

	}


      // Calculate Susceptibility/Entropy

      if ( (calcdens==2)||(calcdens==3) )
	{
	  TM=DN/betabar;
	  Params.clear();
	  Params.push_back(betabar);

	  Suscep.ReadNChainValue(Nsites,0);
	  Suscep.AddValue(Params,&Aeig,1,true,TM);
	  Suscep.SaveNValue(Nsites,0);
      
	  Entropy.ReadNChainValue(Nsites,0);
	  Entropy.AddValue(Params,&Aeig,1,true,TM);
	  Entropy.SaveNValue(Nsites,0);
	}
  //end if calcdens=2 or 3


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

      //SuscepChain.clear();
//       strcpy(arqSus,"SuscepImp_25_726.dat");
  
//       strcpy(arqname,"SuscepChain.dat");
//       InFile.open(arqname);
//       if (InFile.fail())
// 	{
// 	  cout << "Can't find " << arqname << endl;
// 	  strcpy(arqSus,"SuscepChain.dat");
// 	  calcdens=3;
// 	}

//       Sus=CalcSuscep(Params,&Aeig,1,true);
//       double Sus0=0.0;
//       int nlines=0;
//       if (calcdens==2)
// 	{
// 	  InFile.open("SuscepChain.dat");
// 	  InFile.clear(); 
// 	  InFile.seekg(0, ios::beg); // rewind
// 	  if (InFile.fail())
// 	    {
// 	      cout << "Can't open  SuscepChain.dat"<< endl;
// 	    }
// 	  else
// 	    {
// 	      while ( (!InFile.eof())&&(nlines<=Nsites) ) 
// 		{
// 		  InFile >> daux[0] >> daux[0] >> daux[0] >> Sus0;
// 		  nlines++;
// 		}
// 	    }
// 	  InFile.close();
// 	}
//        else 
// 	Sus-=1.0/8.0; // exclude dot site in the chain.
//       if (Nsites==0) OutFile.open(arqSus);
//       else OutFile.open(arqSus,ofstream::app);
//       OutFile.precision(20);
//       OutFile << scientific << TM << " " << Sus << " " << Sus0 << " " << Sus-Sus0 << endl;
//       OutFile.close();


// 	  Sus=CalcSuscep(Params,&Aeig,1,true);
// 	  double Sus0=0.0;
// 	  int nlines=0;
// 	  // Gets SuscepChain iteration by iteration
// 	  if (calcdens==2)
// 	    {
// 	      InFile.open("SuscepChain.dat");
// 	      InFile.clear(); 
// 	      InFile.seekg(0, ios::beg); // rewind
// 	      if (InFile.fail())
// 		{
// 		  cout << "Can't find " << arqname << endl;
// 		}
// 	      else
// 		{
// 		  while ( (!InFile.eof())&&(nlines<=Nsites) ) 
// 		    {
// 		      InFile >> daux[0] >> daux[0] >> daux[0] >> Sus0;
// 		      nlines++;
// 		    }
// 		}
// 	      InFile.close();
// 	    }
// 	  //if (calcdens==2) Sus0=SuscepChain[Nsites];
// 	  else Sus-=1.0/8.0; // exclude dot site in the chain.

// 	  if (Nsites==0) OutFile.open(arqSus);
// 	  else OutFile.open(arqSus,ofstream::app);
// 	  OutFile.precision(20);
// 	  OutFile << scientific << TM << " " << Sus << " " << Sus0 << " " << Sus-Sus0 << endl;
// 	  OutFile.close();
