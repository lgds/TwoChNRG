
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include <boost/timer.hpp>

using namespace std;

//CNRGbasisarray CutStates(CNRGarray* Ain, int Ncutoff){

void CutStates(CNRGarray &Ain, CNRGbasisarray &Acut, int Ncutoff){


  //
  // Function that reduces the number of states in Ain
  //

  // No longer like this
  //CNRGbasisarray Acut(*Ain);
  Acut.ClearAll();
  Acut.SyncNRGarray(Ain);
  Acut.SetSCFfromdEn();
  Acut.BlockBegEndBC=Ain.BlockBegEnd;
  Acut.SetBlockBegEndEigVec();

  if (Ncutoff>=Ain.Nstates()) return; // Does nothing


  // Get Ecutoff

  double Ecutoff=Ain.Ecut(Ncutoff);  // Equivalent to (*Ain).Ecut(Ncutoff)

  cout << "Ecutoff = " << Ecutoff << endl;

  Acut.SetKept(Ncutoff);
  //cout << "No blocks before cut: " << Acut.NumBlocks() << endl;
  //
  // Identify which blocks are eliminated
  //

  cout << "Removing blocks..." << endl;

  for (int ibl=0;ibl<Acut.NumBlocks();ibl++){
    double Eblmin=Acut.GetBlockEmin(ibl);
    if (Eblmin>Ecutoff){cout << "ibl = " << ibl << " : Eblmin>Ecut" <<endl;}
    if (Eblmin>Ecutoff){ // not enough. Make sure you're not slicing a block
    //if ( dGTPrec(Eblmin,Ecutoff,0.01/fabs(Ecutoff)) ){
      cout << " Removing block ibl= " << ibl 
	   << " Eblmin = " << Eblmin 
	   << endl;
      Acut.RemoveBlock(ibl);
      ibl--; // Needs this!!
    }
  }

  cout << "...done removing blocks." << endl;

  //
  // Remove all remaining states
  //

//   Acut.PrintEn();

  // Check time
  boost::timer t;

  t.restart();

  cout << "Removing states..." << endl;

  for (int ibl=0;ibl<Acut.NumBlocks();ibl++){ 
    double Eblmax=Acut.GetBlockEmax(ibl);
    double Eblmin=Acut.GetBlockEmin(ibl);

    if (Eblmax>Ecutoff){ // Make sure you're not slicing a block
    //if ( dGTPrec(Eblmax,Ecutoff,0.01/fabs(Ecutoff)) ){

//       cout << " Block " << ibl << " of " << Acut.NumBlocks() << 
// 	" has Emax = " << Eblmax << endl;
//       cout << " and Emin = " << Eblmin << " <= Ecut = " << Ecutoff << endl;
//       cout << " Block begins at " << Acut.GetBlockLimit(ibl,0) << 
// 	" ends at " << Acut.GetBlockLimit(ibl,1) << endl;
//       cout << " First energy: " << Acut.dEn[Acut.GetBlockLimit(ibl,0)] << endl;

      int ipos=Acut.GetBlockEcutPos(ibl,Ecutoff);
      // Assumes states are energy-ordered within block
	Acut.RemoveStatesFromBlock(ibl,ipos,Acut.GetBlockLimit(ibl,1));
    }
    //End if Eblmax > Ecutoff

  }
  // end loop in remaining blocks

  cout << "...done removing states. Total time =" << t.elapsed() << endl;

}


/////////////////////////////////
/////////////////////////////////
/////     OLD version       /////
/////////////////////////////////
/////////////////////////////////


CNRGbasisarray CutStates(CNRGarray* pAin, int Ncutoff){


  //
  // Function that reduces the number of states in Ain
  //

  // No longer like this
  CNRGbasisarray Acut(*pAin);
  Acut.SetSCFfromdEn();
  Acut.BlockBegEndBC=pAin->BlockBegEnd;
  Acut.SetBlockBegEndEigVec();

  if (Ncutoff>=pAin->Nstates()) return(Acut); // Does nothing

  // Get Ecutoff

  double Ecutoff=pAin->Ecut(Ncutoff);  // Equivalent to (*Ain).Ecut(Ncutoff)

  cout << "Ecutoff = " << Ecutoff << endl;
  //cout << "No blocks before cut: " << Acut.NumBlocks() << endl;
  //
  // Identify which blocks are eliminated
  //

  cout << "Removing blocks..." << endl;

  for (int ibl=0;ibl<Acut.NumBlocks();ibl++)
    {
      double Eblmin=Acut.GetBlockEmin(ibl);
      if (Eblmin>Ecutoff)
	{
	 Acut.RemoveBlock(ibl);
	 ibl--; // Needs this!!
	}
    }

  cout << "...done removing blocks." << endl;

  //
  // Remove all remaining states
  //

//   Acut.PrintEn();

  // Check time
  boost::timer t;

  t.restart();

  cout << "Removing states..." << endl;

  for (int ibl=0;ibl<Acut.NumBlocks();ibl++)
    { 
      double Eblmax=Acut.GetBlockEmax(ibl);
      double Eblmin=Acut.GetBlockEmin(ibl);

//       cout << " Block " << ibl << " of " << Acut.NumBlocks() << 
// 	" has Emax = " << Eblmax << endl;
//       cout << " and Emin = " << Eblmin << " <= Ecut = " << Ecutoff << endl;
//       cout << " Block begins at " << Acut.GetBlockLimit(ibl,0) << 
// 	" ends at " << Acut.GetBlockLimit(ibl,1) << endl;
//       cout << " First energy: " << Acut.dEn[Acut.GetBlockLimit(ibl,0)] << endl;

      if (Eblmax>Ecutoff)
	{

	  int ipos=Acut.GetBlockEcutPos(ibl,Ecutoff);
	  // Assumes states are energy-ordered within block
	  Acut.RemoveStatesFromBlock(ibl,ipos,Acut.GetBlockLimit(ibl,1));

//  	  for (int ist=Acut.GetBlockLimit(ibl,0);ist<=Acut.GetBlockLimit(ibl,1);ist++)
//  	    {
//  	      if (Acut.dEn[ist]>Ecutoff)
//  		{
//  		  cout << " Should remove ist= " << ist << " > ipos = " << ipos << " E = " << Acut.dEn[ist] <<endl;
//  		}
// 	      // end Ei>Ecutoff
//  	    }
// 	  // for inside block

	}
      //End if Eblmax > Ecutoff

    }

    cout << "...done removing states. Total time =" << t.elapsed() << endl;


  // Updates StCameFrom

  return(Acut);

}


/////////////////////////////////


