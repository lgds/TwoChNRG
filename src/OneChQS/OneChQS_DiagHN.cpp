
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"

#include <iostream>
#include <vector>
#include <math.h>

#include <boost/timer.hpp>


#include "OneChQS.hpp"



int OneChQS_DiagHN( CNRGmatrix Qm1fNQ,
		    CNRGbasisarray Abasis, CNRGbasisarray SingleSite, 
		    CNRGarray &Aeig, 
		    double eps_N, double chi_N, double Lambda,bool display){




  int NstOld=Aeig.dEn.size();

  double prefactor,auxEl;

  CNRGmatrix HN(Abasis);

  double MyZero=1.0e-15;

  // Check time
  boost::timer t;




  // Setup HN

  // NEW: Special storage

  HN.UpperTriangular=true;

  int icount=0;

  for (int ii=0;ii<Abasis.NumBlocks();ii++)
    {
      double Qi=Abasis.GetQNumber(ii,0);
      double Si=Abasis.GetQNumber(ii,1);

      int sizebl=Abasis.GetBlockSize(ii);

      HN.MatBlockMap.push_back(ii);
      HN.MatBlockMap.push_back(ii);
      HN.MatBlockBegEnd.push_back(icount);
      // New
      if (HN.UpperTriangular)
 	icount+=(sizebl*(sizebl+1))/2;
      else
	icount+=sizebl*sizebl;
      
      HN.MatBlockBegEnd.push_back(icount-1);

      //cout << "  Setting up HN Block : " << ii ;
      //cout << "  size: " << Abasis.GetBlockSize(ii) << endl;

      int ibl=0;
      // Start timer
      t.restart();
      for (int ist=Abasis.GetBlockLimit(ii,0);
	       ist<=Abasis.GetBlockLimit(ii,1);ist++)
	{
	  int type_i = Abasis.iType[ist];

	  int jbl=0;
// 	  for (int jst=Abasis.GetBlockLimit(ii,0);
// 	           jst<=Abasis.GetBlockLimit(ii,1);jst++)
	  int j0;
	  // NEW THING! only calculates half of matrix els!
	  if (HN.UpperTriangular) j0=ist; 
	  else j0=Abasis.GetBlockLimit(ii,0);
	  for (int jst=j0;
	           jst<=Abasis.GetBlockLimit(ii,1);jst++)
	    {
	      int type_j = Abasis.iType[jst];
	      int rold=ij2r(NstOld,Abasis.StCameFrom[ist],
			    Abasis.StCameFrom[jst]);

 	      //cout << "  Limit in j : " << Abasis.BlockBegEnd[ii+1];
 	      //cout << "  ist : " << ist << "  jst : " << jst;
 	      //cout << "  type i : " << type_i << " type j : " << type_j <<endl;
	    
	      if (ist==jst)
		{
		  // Diagonal terms
		  auxEl=sqrt(Lambda)*Abasis.dEn[ist];
		  HN.MatEl.push_back(auxEl);
		  //cout << " H( = " << ist <<  "," << jst << ") =" 
		  //     << auxEl << endl;

		}
	      else
		{
		  // Off-diagonal terms: watch out for the h.c. term
		  // Nov 07: New off-diag approach, 
		  // valid for two-channel as well
		  auxEl=Qm1fNQ.GetMatEl(Abasis.StCameFrom[ist],
					Abasis.StCameFrom[jst]);
		  int typep=type_i;
		  int type=type_j;

      		  // if Zero, calculate <i|h.c|j>=<j|fd|i>
		  // DIAGONAL BLOCKS: Si=Sj
		  if (fabs(auxEl)<MyZero)
		    {
		      auxEl=Qm1fNQ.GetMatEl(Abasis.StCameFrom[jst],
					    Abasis.StCameFrom[ist]);
		      typep=type_j;
		      type=type_i;
		    }

		  // Find Qcfp, Scfp, Scf, Stildep, Stilde

		  // Qtildep, Qtilde

		  double Qtildep=SingleSite.GetQNumber(typep,0);
		  double Stildep=SingleSite.GetQNumber(typep,1);
		  double Sztildep=SingleSite.GetQNumber(typep,2);

		  double Qtilde=SingleSite.GetQNumber(type,0);
		  double Stilde=SingleSite.GetQNumber(type,1);
		  double Sztilde=SingleSite.GetQNumber(type,2);

		  // Watch out: fermionic sign: (one for fd_{N+1} f_{N} )
		  // if using fd_{N}f_{N+1}=-f_{N+1}fd_{N} put a minus sign in front

		  double FermiSign=1.0;
		  if (dEqual(Qtilde,0.0)) FermiSign=-1.0;

		  double Scfp=Si-Sztildep;
		  double Scf=Si-Sztilde;

		  // Round-off things
		  //Stilde=0.5*((int)(2.0*Stilde));

		  

		  //cout << " fN-1 = " << auxEl << endl;
		  //cout << "Si = " << Si 
		  //    << " typep = " << typep 
		  //    << " type = " << type << endl; 
 		  //cout << "Scfp = " << Scfp 
 		  //     << " Scf = " << Scf << endl;
		  //
 		  //cout << "Stildep = " << Stildep 
 		  //     << " Stilde = " << Stilde
		  //     << " FermiSign = " << FermiSign << endl;

		  double CoefMatEl=0.0;
		  double auxCG[]={0.0,0.0,0.0,0.0};
		  double siteqnumsp[]={Qtildep,Stildep,Sztildep};
		  double siteqnums[]={Qtilde,Stilde,Sztilde};
		  // Sum in Sztilde
		  for (Sztilde=-Stilde;Sztilde<=Stilde;Sztilde+=1.0)
		    {
		      // Sum in sigma: Note that Sztildep=Sztilde+sigma
		      siteqnums[2]=Sztilde;
		      int sitestate=SingleSite.GetBlockFromQNumbers(siteqnums);
		      for (double sigma=-0.5;sigma<=0.5;sigma+=1.0)
			{
			  siteqnumsp[2]=Sztilde+sigma;
			  // Test Find site state
			  int sitestatep=SingleSite.GetBlockFromQNumbers(siteqnumsp);

			  //cout << "Sz~ = " << Sztilde 
			  //   << " sigma = " << sigma << endl;
	
			  auxCG[0]=CGordan(Scf,Si-Sztilde,
					   Stilde,Sztilde,Si,Si);
			  //cout << "CG1( " << Scf << "," << Si-Sztilde << "," 
			  //              <<Stilde<< "," <<Sztilde<< "," 
			  //              <<Si<< "," <<Si<< ") = " 
			  //   << auxCG[0] << endl;
			  auxCG[1]=CGordan(Scfp,Si-Sztilde-sigma,
					   Stildep,Sztilde+sigma,Si,Si);
			  // cout << "CG2( " << Scfp << "," << Si-Sztilde-sigma << "," 
			  //              <<Stildep<< "," <<Sztilde+sigma<< "," 
			  //               <<Si<< "," <<Si<< ") = " 
			  //   << auxCG[1] << endl;
			  //cout << "CG3( " << Scfp << "," << Si-Sztilde-sigma << "," 
			  //     <<0.5<< "," <<sigma<< "," 
			  //     <<Scf<< "," <<Si-Sztilde<< ") = ";
			  auxCG[2]=CGordan(Scfp,Si-Sztilde-sigma,0.5,sigma,Scf,Si-Sztilde);
			  //cout << auxCG[2] << endl;
			  //cout << "typep = "<< typep << "type = " << type << endl;
			  //cout << "sitestp = "<< sitestatep 
			  //     << " sitest = " << sitestate << endl;
 
			  auxCG[3]=fd_table((int)(sigma/0.5),sitestatep,sitestate);
			  //cout << "<sitestp|f|sitest> = " << auxCG[3] << endl;
			  CoefMatEl+=auxCG[0]*auxCG[1]*auxCG[2]*auxCG[3];
			;			
			}
		    }

		  //auxEl*=CoefMatEl;
		  auxEl*=CoefMatEl*FermiSign;
		  //cout << " auxEl = " << auxEl << endl;
		  auxEl*=chi_N;
		  HN.MatEl.push_back(auxEl);
		}
	      // end if ist==jst
	      jbl++;
	    }
	  // end loop in jst
	  ibl++;
	}
      // end loop in ist
      cout << " Block : " << ii 
	   << " size: " << sizebl 
	   <<"  Time for set-up: " << t.elapsed() << endl;
    }
  // end loop in blocks

   cout << "Num Blocks in HN     : " << HN.NumBlocks() << endl;


//   cout << "HN : " << endl;

//   int jj=0;
//   for (int ii=0;ii<HN.MatEl.size();ii++)
//      {
//        cout << "Block = " << HN.MatBlockMap[jj] << " " << HN.MatBlockMap[jj+1];
//        cout << "  i in block = " << HN.vec[ii].ibl 
// 	    << "  j in block = " << HN.vec[ii].jbl;
//        cout << "  i = " << HN.vec[ii].ist << "  j = " << HN.vec[ii].jst;
//        cout << "  val = " << scientific << HN.MatEl[ii] << endl;
//        jj+=2;

//      }

  // Set FilterMap_BegEnd:
  // Takes tooo much time!
  //it will get rid of repetitions 
  //in the block map and set MatBlockBegEnd

//    cout << "HN MatBlockMap : ";
//    for (int ii=0;ii<HN.MatBlockMap.size();ii++)
//      cout << HN.MatBlockMap[ii] << " ";
//    cout << endl;

   //HN.FilterMap_SetBegEnd();

//    cout << "Filtered HN MatBlockMap : ";
//    for (int ii=0;ii<HN.MatBlockMap.size();ii++)
//      cout << HN.MatBlockMap[ii] << " ";
//    cout << endl;
//    cout << "HN MatBlockBegEnd : ";
//    for (int ii=0;ii<HN.MatBlockBegEnd.size();ii++)
//      cout << HN.MatBlockBegEnd[ii] << " ";
//    cout << endl;




  // Some debugging:

   if (display){
     cout << "Basis and HN for Nshell = " << Abasis.Nshell << endl;
     for (int ibl=0;ibl<Abasis.NumBlocks();ibl++){
       Abasis.PrintBlockBasis(ibl);
       HN.PrintMatBlock(ibl,ibl);
     }
   }
   // end if display

  // Update Aeig

  cout << " Updating Aeig. " << endl;

  Aeig.ClearAll();

  // This works!
  Aeig=HN;
  
  Aeig.Nshell=Abasis.Nshell;
  // Diagonalize blocks
  Aeig.dEn.clear();
  Aeig.dEigVec.clear();
  cout << "HN eigenvalues for Nshell = " << Abasis.Nshell << endl;
  for (int ii=0;ii<HN.NumBlocks();ii++)
    {
      HN.DiagBlock(ii, Aeig.dEn, Aeig.dEigVec);
      if (display){Aeig.PrintBlock(ii);}
    }

  //Set min to zero
  Aeig.SetE0zero();

  return(0);
}


// 		  if ( ((type_i==1)&&(type_j==2))||
// 		       ((type_i==2)&&(type_j==1))||
// 		       ((type_i==1)&&(type_j==3))||
// 		       ((type_i==3)&&(type_j==1))||
// 		       ((type_i==2)&&(type_j==4))||
// 		       ((type_i==4)&&(type_j==2))||
// 		       ((type_i==3)&&(type_j==4))||
// 		       ((type_i==4)&&(type_j==3)) )
// 		    {
// 		      auxEl=OneChQS_prefactorHN(type_i,type_j,Si)*
// 		       ( Qm1fNQ.GetMatEl(Abasis.StCameFrom[ist],
// 					Abasis.StCameFrom[jst])+
// 			 Qm1fNQ.GetMatEl(Abasis.StCameFrom[jst],
// 					 Abasis.StCameFrom[ist]) );
// 		      auxEl*=chi_N;
// 		      HN.MatEl.push_back(auxEl);
// 		    }
// 		  else  HN.MatEl.push_back(0.0);
