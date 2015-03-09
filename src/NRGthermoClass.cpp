

#include "NRGclasses.hpp"
#include <sys/stat.h>
//#include "NRGfunctions.hpp"
using namespace std;


//////////////
void CNRGthermo::AddValue(vector<double> Params, 
			  CNRGarray* pAeig, 
			  int sqnumber,bool totalS,
			  double Temp){

  double aux=Calc(Params, pAeig,sqnumber,totalS);


  if (CalcChain) aux-=dImpValue;

  dValues.push_back(aux);
  dTempValues.push_back(Temp);

}
// end add value


//////////////

int CNRGthermo::ReadChain(){
  
  // instream
  ifstream InFile;
  double daux,Sus;
  int nlines=0;

  if (CalcChain) return(-1);

  dChainValues.clear();
  
  InFile.open(ChainArqName);
  if (InFile.fail())
    {
      cout << "ReadChain: Can't find chain file " << ChainArqName << endl;
      cout << " Setting " << ArqName << " as " << ChainArqName << endl;
      CalcChain=true;
      strcpy(ArqName,ChainArqName);
      return(-1);
    }
  else
    {
      while (!InFile.eof())
	{
	  InFile >> daux >> daux >> daux >> Sus;
	  dChainValues.push_back(Sus);
	}
      dChainValues.pop_back();
      InFile.close();
    }

  cout << " ReadChain: Read " << dChainValues.size() 
       << " values" << endl;
  return(0);
}
// end ReadChain

//////////////

void CNRGthermo::ReadNChainValue(int Nsites,int Nsites0){
  
  // instream
  ifstream InFile;
  double daux,Sus0;
  int nlines=0;

  if (CalcChain) return;
  if (Nsites0>Nsites) return;

  InFile.open(ChainArqName);
  InFile.clear(); 
  InFile.seekg(0, ios::beg); // rewind
  if (InFile.fail())
    {
      cout << "Can't find " << ChainArqName << endl;
      return;
    }
  else
    {
      while ( (!InFile.eof())&&(nlines<=(Nsites-Nsite0Chain)) ) 
	{
	  InFile >> daux >> daux >> daux >> Sus0;
	  nlines++;
	}
    }
  InFile.close();

  if ((Nsites-Nsites0)<dChainValues.size())
    dChainValues[Nsites-Nsites0]=Sus0;
  else
    if ((Nsites-Nsites0)==dChainValues.size())
      dChainValues.push_back(Sus0);
    else
      {
	cout << " ReadNChain: Nsites-Nsites0 longer than current vector size " << endl;
	return;
      }
}
// end ReadChainValue


//////////////

void CNRGthermo::SaveNValue(int Nsites, int Nsites0){

  // outstream
  ofstream OutFile;

  char outarqname[64];

  double outvals[3];

  if ( (Nsites<Nsites0)||(Nsites-Nsites0)>=dValues.size() ) return;
  
  outvals[0]=dTempValues[Nsites-Nsites0];
  outvals[1]=dValues[Nsites-Nsites0];


  if (CalcChain)
    {
      strcpy(outarqname,ChainArqName);
      outvals[2]=0.0;
    }
  else
    {
      strcpy(outarqname,ArqName);
      outvals[2]=dChainValues[Nsites-Nsites0];
    }

  cout << "Thermo: Nsites = " << Nsites << " - Saving value in " << outarqname << endl;

  if (Nsites==Nsites0) OutFile.open(outarqname);
  else OutFile.open(outarqname,ofstream::app);
  OutFile.precision(20);
  OutFile << scientific 
	  << outvals[0] << " " 
	  << outvals[1]  << " " 
	  << outvals[2] << " " 
	  << outvals[1]-outvals[2] << endl;
  OutFile.close();


}
// end SaveNValue


//////////////

void CNRGthermo::SaveToFile(){

  int ii=0;

  // outstream
  ofstream OutFile;

  char outarqname[64];

  if (CalcChain)
    strcpy(outarqname,ChainArqName);
  else
    strcpy(outarqname,ArqName);


  OutFile.open(outarqname);
  OutFile.precision(20);
  for (int ii=0;ii<dValues.size();ii++)
    {
      OutFile << scientific 
	      << dTempValues[ii] << " " 
	      << dValues[ii]+dChainValues[ii] << " " 
	      << dChainValues[ii] << " " 
	      << dValues[ii] << endl;
    }
  OutFile.close();
      

}
// end SaveToFile


