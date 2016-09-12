
#include <iostream>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "SpecFuncClass.hpp"

#include <vector>
using namespace std;



void DM_NRG_CalcSpecFunc_ij(CSpecFunction* pSpec,
			    CNRGmatrix** OpArrayN,
			    int iop, int jop,int UseCFS, bool UseGap, int Nw){

  // Nw : number of omegas in each Shell

  //double omega=0.1;

  //double rho_w=0.0;
   

  //CSpecFunction spec1;

  char Ciop[4],Cjop[4];

   if (iop==jop)
     pSpec->NonDiagGF=false;
   else{
     pSpec->NonDiagGF=true;
     int ibl1=0;
     int ibl2=0;
     //int Mtemp=pSpec->Mtemp-2;
     int Mtemp=1;

     // Ok... But OpArrayN[iop] is STILL a vector size NshellMax! 
     // For which Nshell? How about... Mtemp? 

//      OpArrayN[iop][Mtemp].PrintAllBlocks();
//      OpArrayN[iop][Mtemp].PrintMatBlock(0);

     //old
     //OpArrayN[iop]->GetBlocksFromMatBlock(0,ibl1,ibl2);
     OpArrayN[iop][Mtemp].GetBlocksFromMatBlock(0,ibl1,ibl2);
 
     // old
     //if (OpArrayN[jop]->FindMatBlock(ibl1,ibl2)<0){
     if (OpArrayN[jop][Mtemp].FindMatBlock(ibl1,ibl2)<0){
       cout << " Spectral function rho_"<< iop << "_" << jop << " apparently is not defined. Check operators." << endl;
       return;
     }
     // if blocks do not match
   }
   // end if diagonal

   pSpec->Op1N=OpArrayN[iop];
   if (pSpec->NonDiagGF)
     pSpec->Op2N=OpArrayN[jop];


   // Name
   sprintf(Ciop,"%d",iop);
   sprintf(Cjop,"%d",jop);
   strcpy(pSpec->Name,"rho_");
   strcat(pSpec->Name,Ciop);
   strcat(pSpec->Name,"_");
   strcat(pSpec->Name,Cjop);


   int Nomegas=( pSpec->NshellMax % 2 == 0  ? (pSpec->NshellMax / 2)
 		 : (pSpec->NshellMax-1)/2 );
   
   double factorWN=pow(pSpec->Lambda,1.25);

   if (UseCFS) {
     cout << " Calculating Norm of rho_"<<iop<< "_" <<jop << endl;
     cout << " NormCFS= " << pSpec->CalcNorm(UseCFS) << endl;}

   cout << " Calculating Spectral function rho_"<< iop << "_" << jop << "... " << endl;


   if (Nw>0)
     pSpec->CalcSpec_ManyOmegas(Nw,factorWN,UseCFS);
   else
     pSpec->CalcSpecDM_NRG_FixedOmegas(factorWN,UseCFS);

   cout << " ... done." << endl;

   pSpec->SaveOmegaRhow();

   cout << "NormInteg= " << pSpec->CalcNormInteg() << endl;

   if (UseGap){
     // Sub-gap spectrum. Testing.
     int LastNsh=pSpec->NshellMax-1;

     vector < vector<double> > Rhow_SubGap;

     Rhow_SubGap.push_back( vector<double> () ); // omegas
     Rhow_SubGap.push_back( vector<double> () ); // rhows

     // subGap data for negative omega
     pSpec->GetSubGapData(LastNsh,Rhow_SubGap[0],Rhow_SubGap[1],true);
     // subGap data for positive omega
     pSpec->GetSubGapData(LastNsh,Rhow_SubGap[0],Rhow_SubGap[1],false);

     cout << " DM_NRG_CalcSpec: SubGap data: " << endl;

     for (int ii=0; ii<Rhow_SubGap[0].size(); ii++){
       cout << Rhow_SubGap[0][ii] << "   " << Rhow_SubGap[1][ii] << endl;
     }

     char Zvalue[10];
     sprintf(Zvalue,"%4.2f",pSpec->z_twist);

     char filename[50];
     strcpy(filename,"rho_");
     strcat(filename,Ciop);
     strcat(filename,"_");
     strcat(filename,Cjop);
     strcat(filename,"_SubGap");
     if (dNEqual(pSpec->z_twist,1.0)){
       strcat(filename,"_zEQ");
       strcat(filename,Zvalue);
     }
     strcat(filename,".dat");

     ofstream OutFile(filename, ios::out);

     if (!OutFile){cout << "Cannot save SubGap data data in " << filename << endl; return;}

     OutFile.setf(ios::scientific,ios::floatfield);
     OutFile.precision(15);

     for (int ii=0; ii<Rhow_SubGap[0].size(); ii++){
       OutFile << Rhow_SubGap[0][ii] << "   " << Rhow_SubGap[1][ii] << endl;
     }
     OutFile.close();
   }
   // end if UseGap
}
// end subroutine




///////////////////
///  OLD STUFF  ///
///////////////////

///////////////////////////////////////////////////////

void DM_NRG_CalcSpecFuncs(CNRGCodeHandler* pThisCode,
			  CNRGbasisarray* AcutN,
			  CNRGmatrix* RhoN,
			  CNRGmatrix** OpArrayN,
			  int iop, int jop){


  // OLD STUFF!! Missing some details in finite-T...
  cout << " Spectral functions (Working on this)... "<< endl;

//   cout << " Mat 0 Nshell=0 " << endl;
//   OpArrayN[0][0].PrintAllBlocks();
//   cout << " Mat 0 Nshell=2 " << endl;
//   OpArrayN[0][2].PrintAllBlocks();

   double omega=0.1;

   double rho_w=0.0;
   
   vector < vector<double> > Rhow_Costi;

   CSpecFunction spec1;

   char Ciop[4],Cjop[4];

   spec1.NonDiagGF=false;
   spec1.Lambda=pThisCode->Lambda;
   spec1.NshellMax=pThisCode->Nsitesmax;
   spec1.NshellMin=pThisCode->Nsites0;
   spec1.AcutN=AcutN;
   spec1.RhoN=RhoN;
   spec1.Op1N=OpArrayN[iop];
   //spec1.Op2N=OpArrayN[jop];

   // Name
   sprintf(Ciop,"%d",iop);
   sprintf(Cjop,"%d",jop);
   strcpy(spec1.Name,"rhoCostiTeq0_");
   strcat(spec1.Name,Ciop);
   strcat(spec1.Name,"_");
   strcat(spec1.Name,Cjop);


   // Will this work?
   spec1.BDelta=BroadDelta;
   // BDeltaTemp??
   spec1.dBroad=0.5*log(pThisCode->Lambda);

   // Debugging: Check particle-hole symmetry
//    for (int Nsh=spec1.NshellMax-2; Nsh>=spec1.NshellMin;Nsh-=2){
//      cout << " Checking phs for Nshell = : " << Nsh << endl;
//      spec1.DMNRG_SpecDens_ChkPHS(Nsh);

//    }
   // end check phs


   //Rhow_Costi=spec1.CalcSpecCosti_Teq0();
   // OK replacint it with T=0 CalcSpecCosti_T
   Rhow_Costi=spec1.CalcSpecCosti_T(1000);


   //Rhow_Costi=spec1.CalcSpecCosti_Teq0(0.53,1.0);
   int Nomegas=( spec1.NshellMax % 2 == 0  ? (spec1.NshellMax / 2)
  		 : (spec1.NshellMax-1)/2 );
   cout << " Nomegas = " << Nomegas 
	<< " Size of RhowCosti = " << Rhow_Costi.size()
	<< endl;
   Nomegas=Rhow_Costi.size();

   //for (int inshell=0;inshell<0;inshell++){
   for (int iomega=0;iomega<Nomegas;iomega++){
     double omega=Rhow_Costi[iomega][0];
     cout << "RhoCosti" << Ciop<< "_" << Cjop <<"(omega= "<< Rhow_Costi[iomega][0] << " )= " 
	  << Rhow_Costi[iomega][1] << endl;
     cout << "RhoCosti" << Ciop<< "_" << Cjop <<"(omega= "<< Rhow_Costi[iomega][2] << " )= " 
	  << Rhow_Costi[iomega][3] << endl;


   }
   /// end Costi calculation

   double factorWN=pow(spec1.Lambda,1.25);


   // Debug
//    int Nsh_aux=spec1.NshellMax-1;
//    //int Nsh_aux=2;
//    double rhoaux=spec1.DMNRG_SpecDens_M(factorWN,Nsh_aux);
//    double DN=spec1.CalcDN(Nsh_aux);
//    cout << "RhoDMNRG_(Nshell=" << Nsh_aux  << ")= " 
// 	<< rhoaux/DN << " = " << rhoaux << " / " << DN << endl;
   // Costi finite-T
//    rhoaux=spec1.CalcSpecCosti_T_N(Nsh_aux, 10000000,factorWN);
//    cout << "RhoCosti_(Nshell=" << Nsh_aux  << ")= " 
// 	<< rhoaux/DN << " = " << rhoaux << " / " << DN << endl;
   //return;
   //////

//    cout << " Now: Calculating with DM-NRG: " << endl;

//    cout << " Additional tests " << endl;
//    spec1.CalcSpecDM_NRG_FixedOmegas(factorWN);

//    //   vector< vector<double> >::iterator iitRhow;

//    cout << " Shell Even: " << endl;
//    int NomegasEven=spec1.Omega_Rhow_Even[0].size();
//    for (int ii=0;ii<NomegasEven;ii++){
//      cout << "RhoEvenDM-NRG(omega= "<< spec1.Omega_Rhow_Even[0][ii] << " )= ";
//      if (ii<spec1.Omega_Rhow_Even[1].size()){
//        cout << spec1.Omega_Rhow_Even[1][ii];
//      }
//      cout << endl;
//    }
//    // end loop in Omega_Rhow_Even

//    cout << " Shell Odd: " << endl;
//    int NomegasOdd=spec1.Omega_Rhow_Odd[0].size();
//    for (int ii=0;ii<NomegasOdd;ii++){
//      cout << "RhoOddDM-NRG(omega= "<< spec1.Omega_Rhow_Odd[0][ii] << " )= ";
//      if (ii<spec1.Omega_Rhow_Odd[1].size()){
//        cout << spec1.Omega_Rhow_Odd[1][ii];
//      }
//      cout << endl;
//    }
//    // end loop in Omega_Rhow_Odd


//    cout << " Interpolate: " << endl;
//    int NomegasTot=spec1.Omega_Rhow[0].size();
//    for (int ii=0;ii<NomegasTot;ii++){
//      cout << "RhoInterpolDM-NRG(omega= "<< spec1.Omega_Rhow[0][ii] << " )= ";
//      if (ii<spec1.Omega_Rhow[1].size()){
//        cout << spec1.Omega_Rhow[1][ii];
//      }
//      cout << endl;
//    }
//    // end loop in Omega_Rhow


//    cout << " Green's function at omega=0.0: " << endl;

//    double omega1=0.0;

//    cout << "G11(omega=" << omega1 <<") = " 
// 	<< GSL_REAL(spec1.GreensFunction(omega1)) << " + i "
// 	<< GSL_IMAG(spec1.GreensFunction(omega1))
// 	<< endl; 
   

   spec1.SaveOmegaRhow();

// CalcSpec(Lambda,NshellMax,
// 		  AcutN,RhoN,
// 		  OpArrayN[0],OpArrayN[1],
// 		  omega);

}
// end subroutine


