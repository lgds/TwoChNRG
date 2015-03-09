
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include <algorithm>
#include <math.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
using namespace std;

//////////////

int CNRGarray::NumBlocks(){
  
  return((int)(QNumbers.size()/NQNumbers));

}


//////////////
double CNRGarray::EniBlockj(int i, int j){

  // returns the energy of the i-th state in block j

  double aux=0.0;
  
  if ( (2*j-2)<BlockBegEnd.size() )
    {
      if ( (BlockBegEnd[2*j-2]+(i-1))<dEn.size() )
	{
	  aux=dEn[BlockBegEnd[2*j-2]+(i-1)];
	}
    }
  
  return(aux);

}
//////////////

//////////////
void CNRGarray::PrintQNumbers(){


  int ii;
  cout << "Qnumbers : " << endl;
  for (ii=0;ii<QNumbers.size();ii++)
    {
      cout << fixed << QNumbers[ii] << "  ";
      if ((ii+1) % NQNumbers == 0) cout << endl;
    }
  cout << endl;
  cout << resetiosflags (ios_base::floatfield);


}


//////////////
void CNRGarray::PrintEn(){


  int ii;
  int iblock,pos_in_block;
  for (ii=0;ii<dEn.size();ii++){
    iblock=GetBlockFromSt(ii,pos_in_block);
    cout << "En = " << dEn[ii] << " ---| " ;
    for (int iQn=0;iQn<NQNumbers;iQn++){
      cout << GetQNumber(iblock,iQn) << "  ";
    }
    if (ii<Kept.size()){if(Kept[ii]){cout << "(K)";}}
    cout << ">_N= ";
    cout << Nshell << endl;
  }
}

//////////////
void CNRGarray::PrintBlockQNumbers(int iblock){

  cout << "Block no. " << iblock << " : ";
  for (int ii=0;ii<NQNumbers;ii++) 
    cout << fixed << GetQNumber(iblock,ii) << " ";  
  cout << endl;
  cout << resetiosflags (ios_base::floatfield);

}

//////////////
void CNRGarray::PrintBlockEn(int iblock){

  PrintBlockQNumbers(iblock);

  int iBeg=GetBlockLimit(iblock,0);
  int iEnd=GetBlockLimit(iblock,1);
  cout << " dEn : " << endl;
  if (dEn.size()>0)
    {
      for (int ii=iBeg;ii<=iEnd;ii++) 
	cout << ii << " : " << scientific << dEn[ii] << endl;  
    }
  //end if
  cout << resetiosflags (ios_base::floatfield);
  cout << "Kept States: "<< endl;
  if (Kept.size()>0){
    for (int ii=iBeg;ii<=iEnd;ii++){
      cout << Kept[ii] << " ";
    }
    cout << endl;
  }
  // end if

}

//////////////
void CNRGarray::PrintBlock(int iblock){

  PrintBlockEn(iblock);

  int iBeg=GetBlockLimitEigv(iblock,0);
  int iEnd=GetBlockLimitEigv(iblock,1);

  int NstBl=GetBlockLimit(iblock,1)-GetBlockLimit(iblock,0)+1;
  cout << " dEigVec (rows): " << endl;;
  if (dEigVec.size()>0){
    int i1=1;
    for (int ii=iBeg;ii<=iEnd;ii++){
      cout << scientific << dEigVec[ii] << "  ";
      if (i1%NstBl==0) cout << endl;
//       // complex
//       cout << scientific << cEigVec[ii] << "  ";
//       if (i1%NstBl==0) cout << endl;
      i1++;
    }
  } else if (cEigVec.size()>0){
    int i1=1;
    for (int ii=iBeg;ii<=iEnd;ii++){
      // complex
      cout << scientific << cEigVec[ii] << "  ";
      if (i1%NstBl==0) cout << endl;
      i1++;
    }
  } 
  //end if
  
  cout << resetiosflags (ios_base::floatfield);

}


//////////////
void CNRGarray::PrintAll(){


  int ii;
  //void CNRGarray::PrintQNumbers();

  PrintQNumbers();

  cout << " Any totalS QNs? = " << totalS << endl;
  if (Sqnumbers.size()>0){
    cout << "  totalS QNs positions : ";
    for (ii=0;ii<Sqnumbers.size();ii++)
      cout << Sqnumbers[ii] << " ";
    cout << endl;
  }
  cout << "dEn : " << endl;
  for (ii=0;ii<dEn.size();ii++)
    cout << scientific << dEn[ii] << " ";
  cout << endl;

  cout << "dEigVec (rows): " << endl;
  for (ii=0;ii<dEigVec.size();ii++)
    cout << scientific << dEigVec[ii] << " ";
  cout << endl;
  for (ii=0;ii<cEigVec.size();ii++)
    cout << scientific << cEigVec[ii] << " ";
  cout << endl;


  cout << "BlockBegEnd : " << endl;
  for (ii=0;ii<BlockBegEnd.size();ii+=2)
    cout << BlockBegEnd[ii] << " "<< BlockBegEnd[ii+1] << ", ";
  cout << endl;

  cout << "Size Blocks : " << endl;
  for (ii=0;ii<NumBlocks();ii++)
    cout << GetBlockSize(ii) << " ";
  cout << endl;

  cout << "Child States (size= "<< ChildStates.size() <<" ):"<< endl;
  // Why clear?
  //ChildStates.clear();
  if (ChildStates.size()>0){
    for (ii=0;ii<ChildStates.size();ii++){
      if (ChildStates[ii].size()>0){
	for (int jj=0;jj<ChildStates[ii].size();jj++){
	  cout << ChildStates[ii][jj] << " ";
	}
      cout << endl;
      }
      // end if
    }
    cout << endl;
  }
  // end if


  cout << "Kept States (size= "<< Kept.size() <<" ):"<< endl;
  if (Kept.size()>0){
    for (ii=0;ii<Kept.size();ii++){
      cout << Kept[ii] << " ";
    }
    cout << endl;
  }
  // end if





  cout << resetiosflags (ios_base::floatfield);
}
//////////////
int CNRGarray::Nstates(){

  // Changing this:
  //return(BlockBegEnd[NQNumbers*(NumBlocks()-1)+1]-BlockBegEnd[0]+1);
  return(BlockBegEnd.back()-BlockBegEnd.front()+1);
  
}

//////////////
void CNRGarray::SetE0zero(){

  double Emin=*min_element(dEn.begin(),dEn.end());
  cout <<  "Emin = " << Emin << endl;

  for (int ii=0;ii<dEn.size();ii++) dEn[ii]-=Emin;

 
  

}

//////////////
void CNRGarray::SetDegen(){

  int nst,ii,jj;

  // Needs work
  
//   iDegen.push_back(1);
//   ii=1;
//   while ( ii<dEn.size() )
//     {
//       if (fabs(dEn[ii]-dEn[ii-1])<1.0E-10) nst=iDegen[ii-1]+1;
//       else nst=1;
//       iDegen.push_back(nst);
//       ii++;
//     }

}

///////////////////////
void CNRGarray::FilterQNumbers(){

  int equal=0;
  int ii,i1;

  vector<double>::iterator qnums_iter,qnums_iter2,qnums_ii1,qnums_ii2;

  for (qnums_iter=QNumbers.begin(); qnums_iter<QNumbers.end(); 
       qnums_iter+=NQNumbers)
    {

      //cout << " QN1 = " << (*qnums_iter) << " , " << *(qnums_iter+1)<< endl;
      ii=1;
      qnums_iter2=qnums_iter+ii*NQNumbers;
      while (qnums_iter2<QNumbers.end())
	{
	  //cout << " QN2 = " << (*qnums_iter2) << " , " << *(qnums_iter2+1)<< endl;
	  equal=1;
	  i1=0;
	  for (qnums_ii2=qnums_iter2;
	       qnums_ii2<qnums_iter2+NQNumbers;qnums_ii2++)
	    {
	      qnums_ii1=qnums_iter+i1;
	      if ( (*qnums_ii2==*qnums_ii1)&&(equal==1) ) equal=1;
	      else equal=0;
	      i1++;
	    }
	  if (equal==1)
	    {
	      //cout << "Found equal " << endl;
	      QNumbers.erase(qnums_iter2,qnums_iter2+NQNumbers);
	      ii--;
	    }
	  ii++;
	  qnums_iter2=qnums_iter+ii*NQNumbers;
	}
    }

}

//////////////
void CNRGarray::ClearAll(){

  QNumbers.clear();
  dEn.clear();
  dEigVec.clear();
  iDegen.clear(); 

  BlockBegEnd.clear();
  
  ChildStates.clear();

  Kept.clear();

}

//////////////
double CNRGarray::GetQNumber(int iblock, int whichqn){

  int iq=iblock*NQNumbers+whichqn;
      
  return(QNumbers[iq]);
}

//////////////
int CNRGarray::GetBlockLimit(int iblock, int whichlimit){

  if(iblock<0)
    return(0);
  else
    return(BlockBegEnd[2*iblock+whichlimit]);
  
}

//////////////
int CNRGarray::GetBlockSize(int iblock){

  return(BlockBegEnd[2*iblock+1]-BlockBegEnd[2*iblock]+1);
  
}

//////////////
int CNRGarray::GetBlockSize(int iblock, bool kp){

  int Nstkp=0;

  if (Kept.size()<=BlockBegEnd[2*iblock+1]) return(0);

  for (int ii=BlockBegEnd[2*iblock];ii<=BlockBegEnd[2*iblock+1];ii++){

    if (Kept[ii]==kp) Nstkp++;

  }

  return(Nstkp);
  
}

//////////////
int CNRGarray::GetBlockFromSt(int ist,int &pos_in_block){

  // Find which block ist belongs to:
  // FInds iblock such that 
  // BlockBegEnd[2*iblock]<=ist<=BlockBegEnd[2*iblock+1]

  int iblock=0;

  while (ist>BlockBegEnd[2*iblock+1]) {iblock++;}

  pos_in_block=ist-BlockBegEnd[2*iblock];

  // if state is not there (holes in BlockBegEnd, pos_in_block is negative.


  return(iblock);
}
/////////////
int CNRGarray::GetBlockLimitEigv(int iblock, int whichlimit){

  // Returns position in dEigVec corresponding to block iblock

  int posEig=0;
  int NstBl=0;

  for (int ii=0;ii<iblock;ii++)
    {
      NstBl=GetBlockSize(ii);
      posEig+=NstBl*NstBl;
    }

  NstBl=GetBlockSize(iblock);
  posEig+=whichlimit*(NstBl*NstBl - 1);

  return(posEig);

}

/////////////
double CNRGarray::GetBlockEmin(int iblock){

  double Emin = *min_element(dEn.begin()+GetBlockLimit(iblock,0),
			     dEn.begin()+GetBlockLimit(iblock,1)+1);

  return(Emin);
}

/////////////
double CNRGarray::GetBlockEmax(int iblock){

  double Emax = *max_element(
			     dEn.begin()+GetBlockLimit(iblock,0),
			     dEn.begin()+GetBlockLimit(iblock,1)+1);
  
  return(Emax);
}

/////////////
/////////////

// bool GTEcut (double E, double Ecut) {
//   return (E>Ecut);
// }


int CNRGarray::GetBlockEcutPos(int iblock,double Ecut){


  vector<double>::iterator FindEnGTEcut=dEn.begin()+GetBlockLimit(iblock,0);

   while ( (*FindEnGTEcut<=Ecut)
 	  &&(FindEnGTEcut<=dEn.begin()+GetBlockLimit(iblock,1)) )
//  to allow close-by states to be retained (did not work)
//   while ( ( dLTPrec(*FindEnGTEcut,Ecut,0.01/fabs(Ecut)) )
// 	  &&(FindEnGTEcut<=dEn.begin()+GetBlockLimit(iblock,1)) )
    {FindEnGTEcut++;}

//   =find_if(dEn.begin()+GetBlockLimit(iblock,0),
// 	      dEn.begin()+GetBlockLimit(iblock,1),
// 	      GTEcut); // doesnt work

  
  return((int)(FindEnGTEcut-dEn.begin()));
}


/////////////
void CNRGarray::SetKept(int Ncutoff){


  double dEcut=Ecut(Ncutoff);

  Kept.clear();
  vector<double>::iterator it;

  cout << "SetKept: Ncutoff = " << Ncutoff 
       << " Ecut = " << dEcut << endl;

  int Nkept=0;
  int Ndisc=0;

  for (it=dEn.begin();it<dEn.end();it++){
    if ((*it)<=dEcut){
      Kept.push_back(true); Nkept++;
    }else{Kept.push_back(false);Ndisc++;}

  }


  cout << "SetKept: Nkept = " << Nkept << " Ndisc = " << Ndisc << endl;


}
////////////
bool CNRGarray::CheckKept(int ist, bool kp){
  //
  // Checks if Keep[ist] matches kp


  if (Kept.size()<=ist){return(false);}
  else{
    if (Kept[ist]==kp) return(true);
    else return(false);
  }
}

bool CNRGarray::CheckKept(int iblock, int istbl, bool kp){

  int ist=GetBlockLimit(iblock,0)+istbl;

  return(CheckKept(ist,kp));
}


/////////////
int CNRGarray::GetBlockFromQNumbers(double *qnums){

  // uses predicate "dEqual" defined NRGfunctions.hpp

  vector <double>::iterator it=QNumbers.begin()-1;
  bool ok=false;

  // watch out for false alarms
  while (!ok)
   {
  it = search (it+1, QNumbers.end(), qnums, qnums+NQNumbers,dEqual);
  
    ok=(((int)(it-QNumbers.begin()))%NQNumbers==0)||((int)(it-QNumbers.end())==0)?true:false;
   }

  //return (int(it-QNumbers.begin())/NQNumbers);
  // Adding safeguards
  int ibl=(int)(it-QNumbers.begin())/NQNumbers;
  if (ibl>NumBlocks()-1){
    //cout << "GetBlockFromQNumbers: ibl larger than NumBlocks()-1. Returning -1" << endl;
    ibl=-1;
  }

  return (ibl);

  // Adding safeguards


}


///////////// April 08
double CNRGarray::GetQNumberFromSt(int ist, int whichqn){

  int pos=0;
  int ibl=GetBlockFromSt(ist,pos);

  return(GetQNumber(ibl,whichqn));

}
/////////////
boost::numeric::ublas::matrix<double> CNRGarray::EigVec2BLAS(int iblock){

  // Converts the section in dEigVec corresponding to block iblock to uBLAS
  // Format: each eigenvector in a ROW!

  
  int Nst=GetBlockSize(iblock);
  int pos=GetBlockLimitEigv(iblock,0);

  //cout << "NstBlock = " << Nst << "  Pos in dEigVec = " << pos << endl;

  boost::numeric::ublas::matrix<double> Maux(Nst,Nst);

//   copy(dEigVec.begin()+pos,dEigVec.begin()+pos+Nst*Nst,
// 	   Maux.begin2());

  // Let me try this:

  vector<double>::iterator it;

  it=dEigVec.begin()+pos;
  for (int ii=0;ii<Nst;ii++)
    {
      for (int jj=0;jj<Nst;jj++)
	{
	  Maux(ii,jj)=*it;
	  it++;
	}
    }

  //return(boost::numeric::ublas::trans(Maux)); // In Columns
  return(Maux);          // In Rows
  

}
/////////////
boost::numeric::ublas::matrix<complex<double> > CNRGarray::cEigVec2BLAS(int iblock){

  // Converts the section in cEigVec corresponding to block iblock to uBLAS
  // Format: each eigenvector in a ROW! COMPLEX MATRIX

  
  int Nst=GetBlockSize(iblock);
  int pos=GetBlockLimitEigv(iblock,0);

  //cout << "NstBlock = " << Nst << "  Pos in dEigVec = " << pos << endl;

  boost::numeric::ublas::matrix<complex<double> > Maux(Nst,Nst);

//   copy(dEigVec.begin()+pos,dEigVec.begin()+pos+Nst*Nst,
// 	   Maux.begin2());

  // Let me try this:

  vector<complex<double> >::iterator it;

  it=cEigVec.begin()+pos;
  for (int ii=0;ii<Nst;ii++)
    {
      for (int jj=0;jj<Nst;jj++)
	{
	  Maux(ii,jj)=*it;
	  it++;
	}
    }

  //return(boost::numeric::ublas::trans(Maux)); // In Columns
  return(Maux);          // In Rows
 
}
///////


///////////// April 09
double CNRGarray::GetEigVecComponent(int ist,int istbasis){

  // For "uncut" vectors only!
  // A CNRGbasisarray version for this one and for
  // GetBlockLimitEigv should do the job
  // in cut objects.
  // 

  if (dEigVec.size()==0){
    cout << "GetEigVecComponent error: No dEigVec. Return 0." << endl;
    return(0.0);
  }

  int posBl;
  int posBlbasis;
  int ibl=GetBlockFromSt(ist,posBl);
  int ibl_basis=GetBlockFromSt(istbasis,posBlbasis);
  int bl_size=GetBlockSize(ibl);
  int ieigv0=GetBlockLimitEigv(ibl,0);

  if (ibl==ibl_basis){
    vector<double>::iterator it;
    it=dEigVec.begin()+ieigv0+posBl*bl_size+posBlbasis;

    return(*it);
  }
  else{return(0.0);}

}


/////////////
double CNRGarray::Ecut(int Ncut){

  vector<double> Eaux=dEn;
  double aux,aux2;

  sort(Eaux.begin(),Eaux.end());

  
  //  if (Ncut>Eaux.size()) aux=Eaux.back();
  //  else aux=Eaux[Ncut];

// Let's try this
  if ( Ncut>=(Eaux.size()-1) ) aux=Eaux.back();
  else{
    aux=Eaux[Ncut-1]; // state Ncut
    aux2=Eaux[Ncut];  // state Ncut+1
    cout << " CNRGarray::Ecut : Original Ecut = "<< aux 
	 << " Ecut+1 = " << aux2 << endl;
    //check for near-degenerate levels
    int ii=1;
    while ( (aux2>0.0)&&(fabs((aux2-aux)/aux)<0.001)
	    &&(ii<Eaux.size()) ){aux2=Eaux[Ncut+ii];ii++;}
    if (ii>1) aux=Eaux[Ncut+ii-2];
    cout << " CNRGarray::Ecut : Ecut = "<< aux << " Nkept = " << Ncut+ii-1 << endl;
  }

  return(aux);
}



///////////// May 08
int CNRGarray::GetiGS(){

  // returns *position* of the minimum element in dEn 

  vector<double>::iterator it=min_element(dEn.begin(),dEn.end());

  int iGS=0;
  while (it>dEn.begin()){it--;iGS++;}


  return(iGS);

}

void CNRGarray::SaveQSParameters(){

  cout << "Saving in old format " << endl;

  // outstream
  ofstream OutFile;
  char arqname[32];
  char CNsites[8];
  
  sprintf(CNsites,"%d",Nshell);

  //
  // Saving QS parameters
  //

  strcpy(arqname,"DataQSN");
  strcat(arqname,CNsites);
  strcat(arqname,".dat");

  cout << "File is " << arqname << endl;

  OutFile.open(arqname);

  int iGS=GetiGS();
  int pos=0;
  int blockGS=GetBlockFromSt(iGS,pos);

  cout << " Nsites " << Nshell+1 << " Nstates " << Nstates() 
       << " iGS = " << iGS << " BlockGS = " << blockGS
       << " NumQSblocks = " << NumBlocks()
       << endl;

  OutFile << Nshell+1 << " " << Nstates() 
       << " " << iGS << " " << blockGS
       << " " << NumBlocks()
       << endl;

  for (int iblock=0;iblock<NumBlocks();iblock++)
    {
      OutFile << fixed << GetQNumber(iblock,0) << "  " << GetQNumber(iblock,1) << endl; 
    }
  OutFile << resetiosflags (ios_base::floatfield);
  for (int iblock=0;iblock<NumBlocks();iblock++)
    {
      OutFile << GetBlockLimit(iblock,0) << "  " << GetBlockLimit(iblock,1) << endl; 
    }

  OutFile.close();

  //
  // Saving Energies
  //

  strcpy(arqname,"DataEQSN");
  strcat(arqname,CNsites);
  strcat(arqname,".bin");

  // Using C format...
  FILE *pFile;
  
  pFile=fopen(arqname,"wb");
  fwrite (&dEn[0], sizeof(double), dEn.size() ,pFile); 
  fclose(pFile);


  

}
//

void CNRGarray::SaveBin(char arqname[]){

  vector<double>::iterator dit;
  vector<int>::iterator iit;
  vector<bool>::iterator boolit;
  vector<complex<double> >::iterator cit;


  ofstream OutFile;
  double daux;
  int iaux;
  bool auxbool;

  complex<double> caux;

  OutFile.open(arqname, ios::out | ios::binary);

  if (!OutFile){cout << "SaveBin: Cannot save data in " << arqname << endl; return;}


  // Save Nshell
  OutFile.write((char*)&Nshell, sizeof(int));
  // Save NQNumbers
  OutFile.write((char*)&NQNumbers, sizeof(int));
  // Save QNumbers
  iaux=QNumbers.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(dit=QNumbers.begin();dit<QNumbers.end();dit++){
    daux=*dit;
    OutFile.write(reinterpret_cast<char*>(&daux), sizeof(double));
  }
  // Save dEn
  iaux=dEn.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(dit=dEn.begin();dit<dEn.end();dit++){ 
    // (&(*dit)) is weird but it is right...
    OutFile.write(reinterpret_cast<char*>(&(*dit)), sizeof(double));
  }
  // Save dEigVec
  iaux=dEigVec.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(dit=dEigVec.begin();dit<dEigVec.end();dit++){ 
    // (&(*dit)) is weird but it is right...
    OutFile.write(reinterpret_cast<char*>(&(*dit)), sizeof(double));
  }
  // Save cEigVec
  iaux=cEigVec.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(cit=cEigVec.begin();cit<cEigVec.end();cit++){ 
    // (&(*cit)) is weird but it is right...
    OutFile.write(reinterpret_cast<char*>(&(*cit)), sizeof(complex<double> ));
  }
  // Save BlockBegEnd
  iaux=BlockBegEnd.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=BlockBegEnd.begin();iit<BlockBegEnd.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }
  // Save iDegen
  iaux=iDegen.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=iDegen.begin();iit<iDegen.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }
  // Save Kept
  iaux=Kept.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(boolit=Kept.begin();boolit<Kept.end();boolit++){ 
    //OutFile.write(reinterpret_cast<char*>(&(*boolit)), sizeof(bool));
    auxbool=*boolit;
    OutFile.write((char*)&auxbool, sizeof(bool));
  }
  // Save totalS,Sqnumbers
  OutFile.write((char*)&totalS, sizeof(bool));
  iaux=Sqnumbers.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=Sqnumbers.begin();iit<Sqnumbers.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }


  OutFile.close();
  
  if (!OutFile.good()){cout << "SaveBin: error saving" << endl;}


}
///

void CNRGarray::ReadBin(char arqname[]){

  // outstream
  ifstream InFile;

  InFile.open(arqname, ios::in | ios::binary);

  ReadBinInStream(InFile);

  InFile.close();

  if (!InFile.good()){cout << "ReadBin: error reading" << endl;}


}
///

void CNRGarray::ReadBinInStream(ifstream &InFile){

  int isize,iaux, iaux2;
  double daux;
  bool boolaux;
  complex<double> caux;


  // Read Nshell
  InFile.read((char*)&Nshell, sizeof(int));
  // Read NQNumbers
  InFile.read((char*)&NQNumbers, sizeof(int));
  // Read QNumbers
  QNumbers.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&daux), sizeof(double));
    QNumbers.push_back(daux);
  }
  // Read dEn
  dEn.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&daux), sizeof(double));
    dEn.push_back(daux);
  }
  // Read dEigVec
  dEigVec.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&daux), sizeof(double));
    dEigVec.push_back(daux);
  }
  // Read cEigVec
  cEigVec.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&caux), sizeof(complex<double> ));
    cEigVec.push_back(caux);
  }
  // Read BlockBegEnd
  BlockBegEnd.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    BlockBegEnd.push_back(iaux2);
  }
  // Read iDegen
  iDegen.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    iDegen.push_back(iaux2);
  }
  // Read Kept
  Kept.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    //InFile.read(reinterpret_cast<char*>(&boolaux), sizeof(bool));
    InFile.read((char*)(&boolaux), sizeof(bool));
    Kept.push_back(boolaux);
  }
  // Read totalS,Sqnumbers
  InFile.read((char*)&totalS, sizeof(bool));
  Sqnumbers.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    Sqnumbers.push_back(iaux2);
  }


}

//////////////
double CNRGarray::PartitionFuncTeq0(){

  // Works only for Q,Sz basis (not anymore)
  // 04/2010: Including 2S+1 degeneracy factor for SU(2) symmetries
  //

  double sum=0.0;

//   vector<double>::iterator dit;

//   for (dit=dEn.begin();dit<dEn.end();dit++){
//     if (dEqualPrec(abs((*dit)),0.0,1e-10)) sum+=1.0;
//   }

  for (int ibl=0;ibl<NumBlocks();ibl++){
    double Si=0.0;
    if (totalS){
      Si=GetQNumber(ibl,Sqnumbers[0]); // only a single SU(2) for now
    }
    
    for (int ist=GetBlockLimit(ibl,0);ist<=GetBlockLimit(ibl,1);ist++){
      if ( dEqualPrec(abs(dEn[ist]),0.0,1e-10) ) sum+=2.0*Si+1.0;
    }
    // end loop in states
  }
  // end loop in blocks

//
 return(sum);


}



//////////////
double CNRGarray::PartitionFunc(double betabar){

  // Works only for Q,Sz basis (not anymore)
  // 04/2010: Including 2S+1 degeneracy factor for SU(2) symmetries
  //


  double sum=0.0;

//   vector<double>::iterator dit;

//   for (dit=dEn.begin();dit<dEn.end();dit++){
//     sum+=exp(-betabar*(*dit));
//   }


  for (int ibl=0;ibl<NumBlocks();ibl++){
    double Si=0.0;
    if (totalS){
      Si=GetQNumber(ibl,Sqnumbers[0]); // only a single SU(2) for now
    }
    
    for (int ist=GetBlockLimit(ibl,0);ist<=GetBlockLimit(ibl,1);ist++){
      sum+=(2.0*Si+1.0)*exp(-betabar*(dEn[ist]));
    }
    // end loop in states
  }
  // end loop in blocks

  return(sum);
}


//////////////
void CNRGarray::SetEigVecToOne(){

  // Sets eigenvectors to (1 0 0, 0 1 0, 0 0 1... on each block in the basis)
  // Creates an "uncut" dEigVec with Nbl x Nbl block components 
  dEigVec.clear();

  for(int ibl=0;ibl<NumBlocks();ibl++){
    // within each block
    int ist0=GetBlockLimit(ibl,0);
    int istF=GetBlockLimit(ibl,1);
    for (int ist1=ist0;ist1<=istF;ist1++){
      for (int ist2=ist0;ist2<=istF;ist2++){
	if (ist1==ist2)
	  dEigVec.push_back(1.0); // diagonal
	else
	  dEigVec.push_back(0.0); // off-diagonal
      }
    }
    // end loops in block states
  }
  // end block loops

}

///////////////


boost::numeric::ublas::matrix<double> CNRGarray::ExpEi2BLAS(int iblock, double betabar){


  int Nst1=GetBlockSize(iblock);

  int ist=GetBlockLimit(iblock,0);

  // NO NEED! 
//   double Si=0.0;
//   if (totalS){
//     Si=GetQNumber(iblock,Sqnumbers[0]); 
//     // only a single SU(2) for now
//   }

  boost::numeric::ublas::matrix<double> dMatAux(Nst1,Nst1);

  for (int ii=0;ii<Nst1;ii++){
    for (int jj=0;jj<Nst1;jj++){
      if (ii==jj)
	// WRONG
// 	dMatAux(ii,jj)=(2.0*Si+1.0)*exp(-betabar*dEn[ist]);
	dMatAux(ii,jj)=exp(-betabar*dEn[ist]);
      else dMatAux(ii,jj)=0.0;
    }
    ist++;
  }
  // end loop

  return(dMatAux);

}
//

///////////////

bool CNRGarray::CheckComplex(){

  bool result=false;

  if ( (cEigVec.size()>0)&&(dEigVec.size()==0) ) result=true;

  return(result);

}


///////////////

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//                                              //
//                                              //
//     Class CNRGbasisarray functions           //
//                                              //
//                                              //
//////////////////////////////////////////////////
//////////////////////////////////////////////////


//////////////
//void CNRGbasisarray::ClearVecsBasis(){
void CNRGbasisarray::ClearAll(){

  CNRGarray::ClearAll();

  iType.clear();
  StCameFrom.clear();
  BlockBegEndBC.clear();
  BlockBegEndEigVec.clear();
  
  StCameFromQNumbers.clear(); 

}
//////////////
void CNRGbasisarray::PrintAll(){

  CNRGarray::PrintAll();

  cout << "Type : " << endl;
  for (int ii=0;ii<iType.size();ii++)
    cout << iType[ii] << " ";
  cout << endl;

  cout << "StCameFrom : " << endl;
  for (int ii=0;ii<StCameFrom.size();ii++)
    cout << StCameFrom[ii] << " ";
  cout << endl;

  cout << "iDegen : " << endl;
  for (int ii=0;ii<iDegen.size();ii++)
    cout << iDegen[ii] << " ";
  cout << endl;
  
  cout << "BlockBegEndBC : " << endl;
  for (int ii=0;ii<BlockBegEndBC.size();ii+=2)
    cout << BlockBegEndBC[ii] << " "<< BlockBegEndBC[ii+1] << ", ";
  cout << endl;

  cout << "BlockBegEndEigVec : " << endl;
  for (int ii=0;ii<BlockBegEndEigVec.size();ii++)
    cout << BlockBegEndEigVec[ii] << " ";
  cout << endl;

  cout << "StCameFromQNumbers : " << endl;
  for (int ii=0;ii<StCameFromQNumbers.size();ii++)
    {
      if ( (ii%NQNumbers_stcf==0) ) cout << ii/NQNumbers_stcf << ": ";
      //cout << ii << " -  " << StCameFromQNumbers[ii] << " ";
      cout << StCameFromQNumbers[ii] << " ";
      if ( ((ii+1)%NQNumbers_stcf==0) ) cout << endl;
    }

}
///////////////////
//////////////
void CNRGbasisarray::PrintBlockBasis(int iblock){

  CNRGarray::PrintBlockQNumbers(iblock);

  int iBeg=GetBlockLimit(iblock,0);
  int iEnd=GetBlockLimit(iblock,1);
  cout << " Beg: " << iBeg << " End: " << iEnd << endl;
  for (int ii=iBeg;ii<=iEnd;ii++){
    cout << ii << " ";
    if (dEn.size()>0)
      cout << scientific << "Eold = " << dEn[ii] << "  "; 
    if (iType.size()>0)
      cout << fixed << "type = " << iType[ii] << "  ";
    if (StCameFrom.size()>0)  
      cout << fixed << "StCFrom = " << StCameFrom[ii]  << "  ";
    if (iDegen.size()>ii)
      cout << fixed << "iDeg = " << iDegen[ii] << "  ";
    cout << resetiosflags (ios_base::floatfield);
    if (StCameFromQNumbers.size()>0)
      {
	cout << "SCF_QNbrs: ";
	for (int jj=ii*NQNumbers_stcf;jj<(ii+1)*NQNumbers_stcf;jj++)
	  cout << StCameFromQNumbers[jj] << " ";
      }
    cout << endl;

  }
  // end loop in ii

}
/////////////
void  CNRGbasisarray::SetSCFfromdEn(){

  // Sets StCameFrom

  StCameFrom.clear();
  for (int ii=0;ii<dEn.size();ii++) StCameFrom.push_back(ii);



}

/////////////

//////////////
int CNRGbasisarray::GetBlockSizeBC(int iblock){

  if(BlockBegEndBC.size()==0){return(GetBlockSize(iblock));}

  return(BlockBegEndBC[2*iblock+1]-BlockBegEndBC[2*iblock]+1);
  
}

//////////////
int CNRGbasisarray::GetBlockFromStBC(int ist,int &pos_in_block){

  // Find which block ist belongs in the BC structure:
  // FInds iblock such that 
  // BlockBegEnd[2*iblock]<=ist<=BlockBegEnd[2*iblock+1]

  
  if (BlockBegEndBC.size()==0){
    int auxst=CNRGarray::GetBlockFromSt(ist,pos_in_block);
    return(auxst);
  }
  else{

    int iblock=0;
    
    while (ist>BlockBegEndBC[2*iblock+1]) {iblock++;}
    
    pos_in_block=ist-BlockBegEndBC[2*iblock];

    return(iblock);
  }
  // end if else
}



//////////////
void CNRGbasisarray::SetBlockBegEndEigVec(){

  BlockBegEndEigVec.clear();

  int Nbls=(int)BlockBegEndBC.size()/2;
  if (Nbls==0) return;

  int aux=0;
  for (int ibl=0;ibl<Nbls;ibl++)
    {
      BlockBegEndEigVec.push_back(aux);
      aux+=GetBlockSizeBC(ibl)*GetBlockSizeBC(ibl)-1;
      BlockBegEndEigVec.push_back(aux);
      aux++;
    }

}




/////////////
void CNRGbasisarray::RemoveBlock(int iblock){

  // Checks StCameFrom

  if ( (dEn.size()==0)||( (dEigVec.size()==0)&&(cEigVec.size()==0) )
       ||(StCameFrom.size()==0) ){
      cout << "Cant remove block: " << iblock << endl;
      return;
    }
  //CNRGbasisarray::PrintAll();

  // Erase elements in dEn, dEigVec, StCameFrom, Kept

  dEn.erase(dEn.begin()+GetBlockLimit(iblock,0),
	    dEn.begin()+GetBlockLimit(iblock,1)+1);

  StCameFrom.erase(StCameFrom.begin()+GetBlockLimit(iblock,0),
		   StCameFrom.begin()+GetBlockLimit(iblock,1)+1);


  Kept.erase(Kept.begin()+GetBlockLimit(iblock,0),
	    Kept.begin()+GetBlockLimit(iblock,1)+1);


  // Watch out
  if (dEigVec.size()>0)
  dEigVec.erase(dEigVec.begin()+BlockBegEndEigVec[2*iblock],
		dEigVec.begin()+BlockBegEndEigVec[2*iblock+1]+1);

  if (cEigVec.size()>0)
  cEigVec.erase(cEigVec.begin()+BlockBegEndEigVec[2*iblock],
		cEigVec.begin()+BlockBegEndEigVec[2*iblock+1]+1);


  // Erase elements in QNumbers

  QNumbers.erase(QNumbers.begin()+NQNumbers*iblock,
		 QNumbers.begin()+NQNumbers*(iblock+1));


  // Re-arrange BlockBegEnd
  // Re-arrange BlockBegEndEigVec

  for (int ibl1=iblock;ibl1<NumBlocks();ibl1++)
    {
      int Size_next=GetBlockSize(ibl1+1);
      BlockBegEnd[2*ibl1+1]=BlockBegEnd[2*ibl1]+Size_next-1;
      BlockBegEnd[2*ibl1+2]=BlockBegEnd[2*ibl1]+Size_next;


      int Size_next2=Size_next*Size_next;
      BlockBegEndEigVec[2*ibl1+1]=BlockBegEndEigVec[2*ibl1]+Size_next2-1;
      BlockBegEndEigVec[2*ibl1+2]=BlockBegEndEigVec[2*ibl1]+Size_next2;

    }
  BlockBegEnd.pop_back();
  BlockBegEnd.pop_back();

  BlockBegEndEigVec.pop_back();
  BlockBegEndEigVec.pop_back();

  // Erase elements in BlockBegEndBC

  BlockBegEndBC.erase(BlockBegEndBC.begin()+2*iblock,
		      BlockBegEndBC.begin()+2*(iblock+1));

}

/////////////


/////////////
void CNRGbasisarray::RemoveState(int ist){

  // Now it does cut egenvectors as well!!

  // Checks vectors

  if ( (dEn.size()==0)||(StCameFrom.size()==0) )
    return;



  int pos_in_bl;
  int iblock=CNRGarray::GetBlockFromSt(ist,pos_in_bl);
  int bl_size=CNRGarray::GetBlockSize(iblock);
  int bl_sizeBC=GetBlockSizeBC(iblock);


  if (BlockBegEnd[2*iblock]==
      BlockBegEnd[2*iblock+1]) CNRGbasisarray::RemoveBlock(iblock);
  else{

    dEn.erase(dEn.begin()+ist);

    StCameFrom.erase(StCameFrom.begin()+ist);

    Kept.erase(Kept.begin()+ist);


    // This doesn't work: 
    //       vector<double>::iterator it=dEigVec.begin()+GetBlockLimitEigv(iblock,0)+bl_size*pos_in_bl;

    //       dEigVec.erase(it,it+bl_size);
    //
    //
    // bl_size should be the original block size but
    // no information on the original block length is available anymore.
    //
    // Needs a "BlockBegEndBC" or something.
    //
    // But we won't be needing the eigenvectors to generate the basis anyway.
    //

    // Now this should work:

    vector<double>::iterator it=dEigVec.begin()+
      BlockBegEndEigVec[2*iblock]+bl_sizeBC*pos_in_bl;

    vector<complex<double> >::iterator cit=cEigVec.begin()+
      BlockBegEndEigVec[2*iblock]+bl_sizeBC*pos_in_bl;


//     if (it+bl_sizeBC>dEigVec.end()){
    if ( (it+bl_sizeBC>dEigVec.end())&&(cit+bl_sizeBC>cEigVec.end()) ) {
      cout << " Ops, Problems in removing eigvec! ";
      cout << "Removing section of EigVec. bl_sizeBC = " << bl_sizeBC;
      cout << " BegEigVec at : " << BlockBegEndEigVec[2*iblock];
      cout << " EndEigVec at : " << BlockBegEndEigVec[2*iblock+1];
      cout << " bl_sizeBC*pos_in_bl : " << bl_sizeBC*pos_in_bl;
      cout << endl;
      exit(0);
    }


    if (it+bl_sizeBC<=dEigVec.end())
      dEigVec.erase(it,it+bl_sizeBC);
    if (cit+bl_sizeBC<=cEigVec.end())
      cEigVec.erase(cit,cit+bl_sizeBC);


    //Update BlockBegEnd
    //Update BlockBegEndEigVec
    BlockBegEnd[2*iblock+1]-=1;
    BlockBegEndEigVec[2*iblock+1]-=bl_sizeBC;

    for (int ibl1=iblock+1;ibl1<NumBlocks();ibl1++){
      BlockBegEnd[2*ibl1]-=1;
      BlockBegEnd[2*ibl1+1]-=1;

      BlockBegEndEigVec[2*ibl1]-=bl_sizeBC;
      BlockBegEndEigVec[2*ibl1+1]-=bl_sizeBC;
    }

  }

}


/////////////
void CNRGbasisarray::RemoveStatesFromBlock(int iblock, int ist1,int ist2){

  // Now it does cut egenvectors as well!!
  // Similar as RemoveState but does that in blocks

  // Checks vectors

  if ( (dEn.size()==0)||(StCameFrom.size()==0) )
    return;

  int NstGone=ist2-ist1+1;

  // Make sure that ist1 and ist2 belong to iblock
  int pos_in_bl[2];
  int ibl_check1=CNRGarray::GetBlockFromSt(ist1,pos_in_bl[0]);
  int ibl_check2=CNRGarray::GetBlockFromSt(ist2,pos_in_bl[1]);
  if ( (ibl_check1!=iblock)||(ibl_check2!=iblock) ){
    cout << " Error in RemoveStatesFromBlock: ist1, ist2  are not within block " << endl;
    return;
  }

  int bl_size=CNRGarray::GetBlockSize(iblock);
  int bl_sizeBC=GetBlockSizeBC(iblock);


  if (BlockBegEnd[2*iblock]==
      BlockBegEnd[2*iblock+1]) CNRGbasisarray::RemoveBlock(iblock);
  else{

    dEn.erase(dEn.begin()+ist1,dEn.begin()+ist2+1); // Check if works

    StCameFrom.erase(StCameFrom.begin()+ist1,StCameFrom.begin()+ist2+1);

    Kept.erase(Kept.begin()+ist1,Kept.begin()+ist2+1);
    //
    // Now the tricky part: the eigenvector
    //

    //
    // Iterator set at begining of first vector 
    //
    vector<double>::iterator it=dEigVec.begin()+
      BlockBegEndEigVec[2*iblock]+bl_sizeBC*pos_in_bl[0];
    vector<complex<double> >::iterator cit=cEigVec.begin()+
      BlockBegEndEigVec[2*iblock]+bl_sizeBC*pos_in_bl[0];
   
//     if (it+bl_sizeBC>dEigVec.end()){
    if ( (it+bl_sizeBC>dEigVec.end())&&(cit+bl_sizeBC>cEigVec.end()) ){
      cout << " Ops, Problems in removing eigvec! ";
      cout << "Removing section of EigVec. bl_sizeBC = " << bl_sizeBC;
      cout << " BegEigVec at : " << BlockBegEndEigVec[2*iblock];
      cout << " EndEigVec at : " << BlockBegEndEigVec[2*iblock+1];
      cout << " bl_sizeBC*pos_in_bl[0] : " << bl_sizeBC*pos_in_bl[0];
      cout << endl;
      return;
    }

    if (it+bl_sizeBC<=dEigVec.end())
      dEigVec.erase(it,it+bl_sizeBC*NstGone); // Should work

    if (cit+bl_sizeBC<=cEigVec.end())
      cEigVec.erase(cit,cit+bl_sizeBC*NstGone); // Should work


    // Re-struct
    //Update BlockBegEnd
    //Update BlockBegEndEigVec
    BlockBegEnd[2*iblock+1]-=NstGone;
    BlockBegEndEigVec[2*iblock+1]-=bl_sizeBC*NstGone;
   
    for (int ibl1=iblock+1;ibl1<NumBlocks();ibl1++){
      BlockBegEnd[2*ibl1]-=NstGone;
      BlockBegEnd[2*ibl1+1]-=NstGone;
       
      BlockBegEndEigVec[2*ibl1]-=bl_sizeBC*NstGone;
      BlockBegEndEigVec[2*ibl1+1]-=bl_sizeBC*NstGone;
    }

  }
  // If not RemoveBlock

}
///////////// April 09
int CNRGbasisarray::GetBlockLimitEigv(int iblock, int whichlimit){

  // Returns position in dEigVec corresponding to block iblock
  // USES BlockBegEndEigVec if available!!


  if (BlockBegEndEigVec.size()==0){
    return(CNRGarray::GetBlockLimitEigv(iblock,whichlimit));
  }
  else{
//    cout << "Overloaded GetBlockLimitEigv (CNRGbasisarray)" << endl;
    return(BlockBegEndEigVec[2*iblock+whichlimit]);
  }

}
///////////////
///////////// April 09
double CNRGbasisarray::GetEigVecComponent(int ist,int istbasis){

  // Overloads in the case of "cut" vectors
  // For "uncut" vectors only!
  // A CNRGbasisarray version for this one and for
  // GetBlockLimitEigv should do the job
  // in cut objects.
  // 

  if  (dEigVec.size()==0){
    cout << "GetEigVecComponent: No dEigVec. Return 0." << endl;
    return(0.0);
  }


  int posBl;
  int posBlbasis;
  int ibl=GetBlockFromSt(ist,posBl);
  int iblbasis=GetBlockFromStBC(istbasis,posBlbasis);
  int bl_sizeBC=GetBlockSizeBC(ibl);

  //  cout << "GetEigVecComponent: Using overloaded version" << endl;

  if (ibl==iblbasis){
    int ieigv0=GetBlockLimitEigv(ibl,0);
    vector<double>::iterator it;
    it=dEigVec.begin()+ieigv0+posBl*bl_sizeBC+posBlbasis;
    return(*it);
  }
  else{return(0.0);}

}



///////////////
///////////////
boost::numeric::ublas::matrix<double> CNRGbasisarray::EigVecCut2BLAS(int iblock){

  // Converts the section in dEigVec corresponding to block iblock to uBLAS
  // Format: each eigenvector in a ROW!

  
  int Nst=GetBlockSize(iblock);
  int pos=BlockBegEndEigVec[2*iblock];

  int NstBC=GetBlockSizeBC(iblock);

  //cout << "NstBlock = " << Nst << "  Pos in dEigVec = " << pos << endl;

  boost::numeric::ublas::matrix<double> Maux(Nst,NstBC);


  // Let me try this:

  vector<double>::iterator it;

  it=dEigVec.begin()+pos;
  for (int ii=0;ii<Nst;ii++){
    for (int jj=0;jj<NstBC;jj++){
      Maux(ii,jj)=*it;
      it++;
    }
  }
  //return(boost::numeric::ublas::trans(Maux)); // In Columns
  return(Maux);          // In Rows
  
}
/// Now the complex one
boost::numeric::ublas::matrix<complex<double> > CNRGbasisarray::cEigVecCut2BLAS(int iblock){

  // Converts the section in cEigVec corresponding to block iblock to uBLAS
  // Format: each eigenvector in a ROW!

  
  int Nst=GetBlockSize(iblock);
  int pos=BlockBegEndEigVec[2*iblock];

  int NstBC=GetBlockSizeBC(iblock);

  //cout << "NstBlock = " << Nst << "  Pos in dEigVec = " << pos << endl;

  boost::numeric::ublas::matrix<complex<double> > Maux(Nst,NstBC);

  vector<complex<double> >::iterator it;

  it=cEigVec.begin()+pos;
  for (int ii=0;ii<Nst;ii++){
    for (int jj=0;jj<NstBC;jj++){
      Maux(ii,jj)=*it;
      it++;
    }
  }
  //return(boost::numeric::ublas::trans(Maux)); // In Columns
  return(Maux);          // In Rows (complex)
  
}


///////////////
boost::numeric::ublas::matrix<double> CNRGbasisarray::EigVecCut2BLAS(int iblock, 
								     bool kp){

  // Converts the section in dEigVec corresponding to block iblock to uBLAS
  // Format: each eigenvector in a ROW!

  // The assumption: NO CUTTING WAS DONE in the basisarray object

  
  int Nst=GetBlockSize(iblock);
  int pos=BlockBegEndEigVec[2*iblock];

  int NstBC=GetBlockSizeBC(iblock); // Should be equal to Nst

  int Nstkp=CNRGarray::GetBlockSize(iblock,kp);

  //cout << "NstBlock = " << Nst << "  Pos in dEigVec = " << pos << endl;

  boost::numeric::ublas::matrix<double> Maux(Nstkp,Nstkp);
  // if Nst != NstBC, then we will hake a problem here: 
  // this matrix will not be square

  vector<double>::iterator it;

  it=dEigVec.begin()+pos;
  int istkp=0;
  for (int ii=0;ii<Nst;ii++){
    int jstkp=0;
    for (int jj=0;jj<NstBC;jj++){
      if (CheckKept(iblock,ii,kp)){
	Maux(istkp,jstkp)=*it;
	jstkp++;
      }
      // if state kept==kp
      it++;
    }
    // end loop in components
    if (CheckKept(iblock,ii,kp)){istkp++;}
  }
  // end loop in states


  //return(boost::numeric::ublas::trans(Maux)); // In Columns
  return(Maux);          // In Rows
  

}



///////////////
///////////////
///////////////
///////////////
///////////// April 08
double CNRGbasisarray::GetStCameFromQNumberFromSt(int ist, int whichqn){

  int stcf=StCameFrom[ist];

  return(StCameFromQNumbers[NQNumbers_stcf*ist+whichqn]);

}
/////////////


///////////////
///////////// June 08
/////////////
void CNRGbasisarray::CopyBlock(int iblock,int NoCopies){

  // Checks StCameFrom

  if ( (NoCopies<1)||(NoCopies>100) )
    {
      cout << " CopyBlock : NoCopies not valid " << endl;
      return;
    }


  // Set iDegen if not set yet
  if (iDegen.size()==0)
    for (int ii=0;ii<Nstates();ii++) iDegen.push_back(0);

  // Set iType if not set yet
  if (iType.size()==0)
    for (int ii=0;ii<Nstates();ii++) iType.push_back(0);

  

  int BlSize=GetBlockSize(iblock);

  vector<double> dAux;
  vector<int> iAux,iAux2;

//   cout << " Bl Size = " << BlSize 
//        << " Limit 0: " << GetBlockLimit(iblock,0)
//        << " Limit 1: " << GetBlockLimit(iblock,1)
//        << endl;


  dAux.assign(dEn.begin()+GetBlockLimit(iblock,0),
 	      dEn.begin()+GetBlockLimit(iblock,1)+1);

  iAux.assign(iType.begin()+GetBlockLimit(iblock,0),
 	      iType.begin()+GetBlockLimit(iblock,1)+1);

  iAux2.assign(StCameFrom.begin()+GetBlockLimit(iblock,0),
 	      StCameFrom.begin()+GetBlockLimit(iblock,1)+1);



  for (int ic=1;ic<=NoCopies;ic++)
    {
      // Copy to dEn
      dEn.insert(dEn.begin()+GetBlockLimit(iblock,1)+1,dAux.begin(),dAux.end());
      // Copy to iType
      iType.insert(iType.begin()+GetBlockLimit(iblock,1)+1,iAux.begin(),iAux.end());
      // Copy to StCameFrom
      StCameFrom.insert(StCameFrom.begin()+GetBlockLimit(iblock,1)+1,iAux2.begin(),iAux2.end());

      // Update iDegen

      iDegen.insert(iDegen.begin()+GetBlockLimit(iblock,1)+1,BlSize,ic);

      // Re-arrange BlockBegEnd, Get ready for next copy
      BlockBegEnd[2*iblock+1]+=BlSize;
      for (int ibl1=iblock+1;ibl1<NumBlocks();ibl1++)
	{
	  BlockBegEnd[2*ibl1]+=BlSize;
	  BlockBegEnd[2*ibl1+1]+=BlSize;
	}
      // end loop in blocks

    }
  iAux.clear();
  dAux.clear();

}

///////////////
/////////////
///////////// September 08
/////////////
void CNRGbasisarray::SetLastQNumber(vector<double> InputVec){

  // Based on the iDegen values, add QNumber and rearrange
  // dEn, QNumbers and BlockBegEnd. Good for phonons.
  // It will re-order the blocks according with iDegen values
  // (it may come down to a single state per block)


  if (InputVec.size()!=Nstates()){
    cout << "SetLastQNumber : Input vec size does not match" << endl; 
    return;
  }

  // assume it is the latest qn

  vector<double> auxdEn;
  vector<double> auxQNumbers;
  vector<int> auxBlockBegEnd;
  vector<int> auxiDegen;

  vector<int> auxiType;
  vector<int> auxStCameFrom;

  vector<double>::iterator it0;
  vector<double>::iterator it1;

  vector<int>::iterator it_ideg;



  double *auxdQN=new double [NQNumbers];


  int NoBlocksOri=NumBlocks(); // will increase QNumbers so need this one
  int Nadded=0;

  // Add new blocks to QNumbers

  auxQNumbers=QNumbers;

  auxdQN[0]=0.0;
  for (int ibl=0;ibl<NoBlocksOri;ibl++){
    int ist0=GetBlockLimit(ibl,0);
    int ist1=GetBlockLimit(ibl,1);
      
    for (int ist=ist0;ist<=ist1;ist++){
      it0=QNumbers.begin()+(ibl+Nadded)*NQNumbers;
      it1=QNumbers.begin()+(ibl+Nadded+1)*NQNumbers;
      // 	  auxdQN[0]=(double)iDegen[ist];
      auxdQN[0]=InputVec[ist];
      QNumbers[(Nadded*NQNumbers)+((ibl+1)*NQNumbers)-1]=auxdQN[0];
      if (ist<ist1)
	{QNumbers.insert(it0,it0,it1);Nadded++;}
    }
    // end loop in states
  }
  // end loop in blocks

  FilterQNumbers();

  // Update Block structure

  auxdEn.clear();
  auxiDegen.clear();
  auxBlockBegEnd.clear();
  auxiType.clear();
  auxStCameFrom.clear();

  Nadded=0;

  // Loop in new block structure
  for (int ibl=0;ibl<NumBlocks();ibl++){
    for (int iqn=0;iqn<NQNumbers;iqn++){auxdQN[iqn]=GetQNumber(ibl,iqn);}

    // Loop in the old block structure      
    auxBlockBegEnd.push_back(Nadded);
    for (int ibl1=0;ibl1<NoBlocksOri;ibl1++){
      bool ok_match=true;
	
      int ist0=GetBlockLimit(ibl1,0);
      int ist1=GetBlockLimit(ibl1,1);
	
      for (int iqn=0;iqn<NQNumbers-1;iqn++)
	if(dNEqual(auxdQN[iqn],auxQNumbers[ibl1*NQNumbers+iqn]))
	  ok_match=false;
      if (ok_match){
	//cout << " Blocks Match: " << " ibl1= " << ibl1 << " ibl= " << ibl  << endl;
	for (int ist=ist0;ist<=ist1;ist++){
	  // 		  if (dEqual(auxdQN[NQNumbers-1],(double)iDegen[ist]))
	  if (dEqual(auxdQN[NQNumbers-1],InputVec[ist])){

	    if (ist<dEn.size()) auxdEn.push_back(dEn[ist]);
	    if (ist<iDegen.size()) auxiDegen.push_back(iDegen[ist]);
	    if (ist<iType.size()) auxiType.push_back(iType[ist]);
	    if (ist<StCameFrom.size()) auxStCameFrom.push_back(StCameFrom[ist]);
	    
	    Nadded++;
	    //cout << " Adding ist = " << ist << " to block " << ibl << endl;
	  }
	}
	// end loop in ist in matchin blocks
      }
      // if match
    }
    // loop in old blocks

    auxBlockBegEnd.push_back(Nadded-1);
  }
  // loop in new blocks

  // Update other vectors
  dEn.clear();
  BlockBegEnd.clear();
  iDegen.clear();
  iType.clear();
  StCameFrom.clear();

  dEn=auxdEn;
  BlockBegEnd=auxBlockBegEnd;
  iDegen=auxiDegen;

  iType=auxiType;
  StCameFrom=auxStCameFrom;

  delete[] auxdQN;

}

///////////


//
/////////////
///////////// February 09
/////////////


void CNRGbasisarray::SaveBin(char arqname[]){

  CNRGarray::SaveBin(arqname);

  vector<double>::iterator dit;
  vector<int>::iterator iit;

  ofstream OutFile(arqname, ios::out | ios::binary | ios::app);
  double daux;
  int iaux;

  if (!OutFile){cout << "SaveBin: Cannot save data in " << arqname << endl; return;}

  // Save iType
  iaux=iType.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=iType.begin();iit<iType.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }

  // Save StCameFrom
  iaux=StCameFrom.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=StCameFrom.begin();iit<StCameFrom.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }

  // Save BlockBegEndBC
  iaux=BlockBegEndBC.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=BlockBegEndBC.begin();iit<BlockBegEndBC.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }

  // Save BlockBegEndEigVec
  iaux=BlockBegEndEigVec.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(iit=BlockBegEndEigVec.begin();iit<BlockBegEndEigVec.end();iit++){ 
    OutFile.write(reinterpret_cast<char*>(&(*iit)), sizeof(int));
  }
  // Save NQNumbers_stcf
  OutFile.write((char*)&NQNumbers_stcf, sizeof(int));
  // Save QNumbers_stcf
  iaux=StCameFromQNumbers.size();
  OutFile.write((char*)&iaux, sizeof(int));
  for(dit=StCameFromQNumbers.begin();dit<StCameFromQNumbers.end();dit++){
    daux=*dit;
    OutFile.write(reinterpret_cast<char*>(&daux), sizeof(double));
  }

  OutFile.close();

  if (!OutFile.good()){cout << "SaveBin: error saving" << endl;}


}
///

void CNRGbasisarray::ReadBin(char arqname[]){

  // ifstream
  ifstream InFile;

  InFile.open(arqname, ios::in | ios::binary);
  //InFile.seekg(0,ios::beg);

  if (!InFile){cout << "ReadBin: Cannot read data in " << arqname << endl; return;}

  CNRGbasisarray::ReadBinInStream(InFile);

  InFile.close();

  if (!InFile.good()){cout << "ReadBin: overall error reading" << endl;}


}

/////
void CNRGbasisarray::ReadBinInStream(ifstream &InFile){


  int isize,iaux, iaux2;
  double daux;


  CNRGarray::ReadBinInStream(InFile);
  
   // Read iType
  iType.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    iType.push_back(iaux2);
  }
  // Read StCameFrom
  StCameFrom.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    StCameFrom.push_back(iaux2);
  }
  // Read BlockBegEndBC BlockBegEigVec
  BlockBegEndBC.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    BlockBegEndBC.push_back(iaux2);
  }
  // Read BlockBegEndEigVec
  BlockBegEndEigVec.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&iaux2), sizeof(int));
    BlockBegEndEigVec.push_back(iaux2);
  }
  // Read NQNumbers_stcf
  InFile.read((char*)&NQNumbers_stcf, sizeof(int));
  // Read QNumbers_stcf
  StCameFromQNumbers.clear();
  InFile.read((char*)&isize, sizeof(int));
  for(iaux=0;iaux<isize;iaux++){
    InFile.read(reinterpret_cast<char*>(&daux), sizeof(double));
    StCameFromQNumbers.push_back(daux);
  }


}
/////////////



///////////////
