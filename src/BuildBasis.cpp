
#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include <iostream>
#include <vector>

//void BuildBasis_RedSymmetry(vector<int> CommonQNs,
void BuildBasis(vector<int> CommonQNs, vector<int> totSpos,
		CNRGbasisarray* pAeigCut,  CNRGbasisarray* pAbasis, 
			    CNRGbasisarray* pSingleSite, int BefCut){

  //
  // Retrieving data from CommonQNs
  //
  // structure:
  //  CommonQNs[0] - Number of common Qnumbers in both basis
  //  CommonQNs[1..NCommon] - position of common QNumbers in the "old" basis
  //  CommonQNs[NCommon+1..2*NCommon] - position in the single site basis
  // 
  //  totSpos[ith] - position of the ith total S qnumber in the "old" basis 
  //             (should be -1 if none)
  //
  // Addition (Sep08): Parity QN
  //
  //  CommonQNs[2*NCommon+1],
  //  CommonQNs[2*NCommon+2], etc. - Positions of "parity" QNs in the "old" basis
  //  
  //
  //
  // Examples: 
  // A)  |Q Sz m1> -> |Q Sz> = |Qold Szold m1>x|Q~ Sz~>
  // then
  //
  // CommonQNs = 2  0 1  0 1
  // totSpos (empty)
  //
  // B)  |Q1 Q2 S m1> -> |Q1 Q2 S> = |Q1old Q2old Sold m1>x|Q1~ Q2~ S~ Sz~>
  // then
  //
  // CommonQNs = 3  0 1 2  0 1 2
  // totSpos = 2
  //
  // C)  |m1 Q S m2> -> |Q S> = |m1 Qold Sold m2>x|Q~ S~ Sz~>
  // then
  //
  // CommonQNs = 2  1 2  0 1
  // totSpos = 2
  //
  // D)  |Q S m1> -> |Q S m1> = |Qold Sold m1>x|Q~ S~ Sz~>
  // then
  //
  // CommonQNs = 3  0 1 2  0 1 -1
  // totSpos = 1
  //(will just repeat the "-1" qnumber)
  //
  // E)  |I1 I2 S m1> -> |I1 I2 S> = |I1old I2old Sold m1>x|I1~ I1z~I2~ I2z~ S~ Sz~>
  // then
  //
  // CommonQNs = 3  0 1 2 0 2 4
  // totSpos = 0 1 2
  //(all 3 qnumbers are "total S-like")
  //
  // F)  |Sz> -> |Q1 Q2 Sz> = |Szold>x|Q1~ Q2~ Sz~>
  // then
  //
  // CommonQNs = 1  0 2
  // totSpos (empty)
  //
  //
  //
  // G) |Q S P m> -> |Q S P> =  |Qold Sold Pold m1>x|Q~ S~ Sz~ P~>
  //
  // CommonQNs= 3  0 1 2  0 1 3  2
  // totSpos = 1
  //
  //
  // H) |Pup Pdn m> -> |Pup Pdn> =  |Pupold Pdnold m1>x|Pup~ Pdn~>
  //
  // CommonQNs= 2  0 1  0 1  0 1
  // totSpos (empty)
  // Notice that there are TWO parity QNS at positions 2*Ncommon+1,2 = 5 and 6 

  if (CommonQNs.size()==0) return;

  int NoCommon=CommonQNs[0];

  if (NoCommon==0){
    cout << "Cant build basis: No common Qnumbers" <<endl;
    exit(0);
  }


  if (CommonQNs.size()<2*NoCommon+1) return;


  //
  // Allocate old QNs
  //
  //
  // Allocate site QNs, new QNs
  //


  double* QNold = new double [NoCommon];
  double* QNsite=new double [NoCommon];
  double* QNnew=new double [NoCommon];

  int* PosInOld = new int [NoCommon];
  int* PosInSingleSite = new int [NoCommon];
  int* TotalS=new int [NoCommon];

  int* MatchOld=new int [NoCommon];

  int* ParityQN=new int [NoCommon];

  int NoQNumsChain=pSingleSite->NQNumbers;

  // Reduction in QNs?

  bool RedSym=false;
  if (pAeigCut->NQNumbers > NoCommon) RedSym=true;


  // Any parity ops? (Sep 2008)
  bool UseParity=false;
  int NoParityQNs=0;
  if (CommonQNs.size()>2*NoCommon+1){
    UseParity=true;
    NoParityQNs=CommonQNs.size()-(2*NoCommon+1);
    cout << " Use parity: " << UseParity << "  " << NoParityQNs << endl;
  }


  // Set the matching QNumbers
  // MatchOld: whether a common QN in AeigOld matches one in SingleSite
  int itotS=0;
  int iParity=0;
  for (int ii=0;ii<NoCommon;ii++){
    PosInOld[ii]=CommonQNs[ii+1];
    PosInSingleSite[ii]=CommonQNs[ii+1+NoCommon];

    // TotalS
    TotalS[ii]=0;
    if ( (itotS<totSpos.size())&&(PosInOld[ii]==totSpos[itotS]) )
      {TotalS[ii]=1;itotS++;}

    // Parity
    ParityQN[ii]=0;
    if ((UseParity)&&(iParity<NoParityQNs)){
      if (PosInOld[ii]==CommonQNs[iParity+2*NoCommon+1])
	{ParityQN[ii]=1;iParity++;}
    }

    // Match
    if ((PosInSingleSite[ii]>=0)&&(PosInSingleSite[ii]<NoQNumsChain) )
      MatchOld[ii]=1;
    else
      MatchOld[ii]=0;      
  }
  //end set MatchOld 

  cout << " NoCommon = " << NoCommon << endl;

  cout << " ParityQN : ";
  for (int ii=0;ii<NoCommon;ii++) cout << ParityQN[ii] << " ";
  cout << endl;

  cout << " PosInOld : ";
  for (int ii=0;ii<NoCommon;ii++) cout << PosInOld[ii] << " ";
  cout << endl;

  cout << " PosInSingleSite : ";
  for (int ii=0;ii<NoCommon;ii++) cout << PosInSingleSite[ii] << " ";
  cout << endl;


  cout << " MatchOld : ";
  for (int ii=0;ii<NoCommon;ii++) cout << MatchOld[ii] << " ";
  cout << endl;


  // Reducing symmetry!
  //int posTotSchain=-1;
  //int NoQNnew=NoCommon;
  // To do: check if total S is a quantum number 
  //        check if total S is in Old and single site
  //   If that is the case,
  //   NoQNnew=NoQNumsChain-1 
  //   because SingleSite will have both S and Sz.
  //  but the new QNs will have only S

  pAbasis->Nshell=pAeigCut->Nshell+1;
  pAbasis->NQNumbers=NoCommon;

  if (RedSym)
    pAbasis->NQNumbers_stcf=pAeigCut->NQNumbers;

  // new (April 2010): totalS and Sqnumbers
  pAbasis->Sqnumbers.clear();
  vector<int>::iterator iit;
  if ( totSpos.size()>0 ){
    pAbasis->totalS=true;
    int pos=0;
    for (iit=totSpos.begin();iit<totSpos.end();iit++){
      //if (*iit!=0) // WHY??? THIS CREATES A PROBLEM FOR OneChS!
      pAbasis->Sqnumbers.push_back(*iit);
      pos++;
    }
  }
  else pAbasis->totalS=false;

  pAbasis->ClearAll();

  //
  // Generate basis
  //

  for (int ibl=0;ibl<pAeigCut->NumBlocks();ibl++){
    //
    // Get old QNumbers
    //
    for (int ii=0;ii<NoCommon;ii++)
      QNold[ii]=pAeigCut->GetQNumber(ibl,PosInOld[ii]);

    //
    // Get site qnumbers
    //

    for (int ibls=0; ibls<pSingleSite->NumBlocks(); ibls++){

      for (int ii=0;ii<NoCommon;ii++){
	if (MatchOld[ii]==1)
	  QNsite[ii]=pSingleSite->GetQNumber(ibls,PosInSingleSite[ii]);
	else{
	  if (ParityQN[ii]==1) QNsite[ii]=1.0; else QNsite[ii]=0.0; 
	}
	// Changed to accomodate parity
	// 	      QNnew[ii]=QNsite[ii]+QNold[ii];
	// 	      pAbasis->QNumbers.push_back(QNsite[ii]+QNold[ii]);
	if (ParityQN[ii]==1)
	  QNnew[ii]=QNsite[ii]*QNold[ii];
	else
	  QNnew[ii]=QNsite[ii]+QNold[ii];
	pAbasis->QNumbers.push_back(QNnew[ii]);
      }
      // end loop in NoCommon
      // If there are SU(2) symmetries, add extra blocks
      // right now works for only one SU(2) symmetry!
      if (totSpos.size()>0){
	for (int ii=0;ii<NoCommon;ii++){
	  if (TotalS[ii]==1){
	    double Snew=QNold[ii]+QNsite[ii]-1.0;
	    while (Snew>=fabs(QNold[ii]-QNsite[ii])){
	      QNnew[ii]=Snew;
	      for (int ii1=0;ii1<NoCommon;ii1++)		       
		pAbasis->QNumbers.push_back(QNnew[ii1]);
	      Snew-=1.0;
	    }
	    // end loop in Snew
	  }
	  // end if TotalS==1
	}
	// end 1st loop in iiNcommon
      }
      // end if SU(2) symmtries are present
    }
    // end loop in pSingleSite
  }
  //End loop in AeigCutblocks

  //pAbasis->PrintQNumbers();

  // Filter blocks
  pAbasis->FilterQNumbers();
  
  //pAbasis->PrintQNumbers();

  // Looks ok. Going on.

  // Set Type, StCameFrom, dEn, BlockBegEnd

  int istnew=0;
  for (int ibl=0;ibl<pAbasis->NumBlocks();ibl++){
    // set Qnew
    for (int inew=0;inew<pAbasis->NQNumbers;inew++)
      QNnew[inew]=pAbasis->GetQNumber(ibl,inew);

    // Mark block begin
    pAbasis->BlockBegEnd.push_back(istnew);

    // Loop in old block
    for (int ibl_old=0;ibl_old<pAeigCut->NumBlocks();ibl_old++){
      //
      // Get old QNumbers
      //
      for (int ii=0;ii<NoCommon;ii++)
	QNold[ii]=pAeigCut->GetQNumber(ibl_old,PosInOld[ii]);
      //
      // Get SingleSite Numbers and check if blocks match
      //
      for (int ibls=0; ibls<pSingleSite->NumBlocks(); ibls++){
	int ok_match=1;
	for (int icom=0;icom<NoCommon;icom++)
	  {
	    if (MatchOld[icom]==1)
	      QNsite[icom]=pSingleSite->GetQNumber(ibls,PosInSingleSite[icom]);
	    else
	      QNsite[icom]=0.0;

	    // 		  cout << "icom = " << icom 
	    // 		       << " Qold = " << QNold[icom]
	    // 		       << " Qsite = " << QNsite[icom]
	    // 		       << " Qnew = " << QNnew[icom]
	    // 		       << " TotalS = " << TotalS[icom]
	    // 		       << " Parity = " <<  ParityQN[icom]
	    // 		       << endl;

	    if (TotalS[icom]==0){
	      // 		      if (!dEqual(QNnew[icom],QNsite[icom]+QNold[icom]))
	      // 			ok_match=0;
	      if (ParityQN[icom]==1){
		if (!dEqual(QNnew[icom],QNsite[icom]*QNold[icom]))
		  ok_match=0;
	      }
	      else{
		if (!dEqual(QNnew[icom],QNsite[icom]+QNold[icom]))
		  ok_match=0;
	      }
	    }
	    else{
	      double Sold=QNold[icom];
	      double Ssite=pSingleSite->GetQNumber(ibls,PosInSingleSite[icom]);
	      // Take Sz
	      double Szsite=pSingleSite->GetQNumber(ibls,PosInSingleSite[icom]+1);
	      double Snew=QNnew[icom];

	      if ( (!dEqual(Snew,Szsite+Sold))||
		   (!TriangIneq(Sold,Ssite,Snew)))
		ok_match=0;
	      // 		    cout << " Sold = " << Sold
	      // 			 << " Ssite = " << Ssite
	      // 			 << " Szsite = " << Szsite
	      // 			 << " Snew = " << Snew
	      // 			 << " ok_match = " << ok_match
	      // 			 << endl;

	    }
	    // end if totalS==0

	  }
	// end loop in QNums



	if (ok_match==1){
	  //  		  cout << " Got a match! " << endl;
	  //  		  cout << " Old QNs : ";
	  //  		  for (int iqn = 0; iqn<pAeigCut->NQNumbers; iqn++)
	  //  		    cout << pAeigCut->GetQNumber(ibl_old,iqn) << "  ";
	  //  		  cout << endl;
	  //  		  cout << "  -- QNold: ";
	  //  		  for (int ii=0;ii<NoCommon;ii++) cout << QNold[ii] << " ";
	  //  		  cout << endl;
	  //  		  cout << " Site QNs : ";
	  //  		  for (int iqn = 0; iqn<pSingleSite->NQNumbers; iqn++)
	  //  		    cout << pSingleSite->GetQNumber(ibls,iqn) << "  ";
	  //  		  cout << endl;
	  // 		  cout << "  -- QNsite: ";
	  // 		  for (int inew=0;inew<NoQNnew;inew++) cout << QNsite[PosInSingleSite[inew]] << " ";
	  // 		  cout << endl;
	  //  		  cout << " New QNs : ";
	  //  		  for (int iqn = 0; iqn<pAbasis->NQNumbers; iqn++)
	  //  		    cout << pAbasis->GetQNumber(ibl,iqn) << "  ";
	  //  		  cout << endl;
	  // 		  cout << "  -- QNnew: ";
	  // 		  for (int inew=0;inew<NoQNnew;inew++) cout << QNnew[inew] << " ";
	  // 		  cout << endl;



	  for (int istbl=pAeigCut->BlockBegEnd[2*ibl_old];
	       istbl<=pAeigCut->BlockBegEnd[2*ibl_old+1];istbl++){
	    for (int ists=pSingleSite->GetBlockLimit(ibls,0);
		 ists<=pSingleSite->GetBlockLimit(ibls,1);ists++){
	      pAbasis->iType.push_back(pSingleSite->iType[ists]);
	      pAbasis->dEn.push_back(pAeigCut->dEn[istbl]);
	      // Has to be in the uncut basis if BefCut=1!
	      if (BefCut==1)
		pAbasis->StCameFrom.push_back(pAeigCut->StCameFrom[istbl]);
	      else
		pAbasis->StCameFrom.push_back(istbl);
	      // New: if RedSym, keep old QNumbers

	      if (RedSym){
		for (int whichqn=0;
		     whichqn<pAeigCut->NQNumbers;
		     whichqn++)
		  pAbasis->StCameFromQNumbers.push_back(pAeigCut->GetQNumberFromSt(istbl,whichqn));
	      }
	      istnew++;
	    }
	    // end loop in SingleSite block states
	  }
	  // end in loop in old states
	}
	// end if match is good
      }
      // end loop in pSingleSite
    }
    //Loop in AeigCut blocks
    pAbasis->BlockBegEnd.push_back(istnew-1);     

  }
  // Loop in new block structure


  cout << "No blocks = " << pAbasis->NumBlocks() << endl;

  cout << "No states = " << pAbasis->Nstates() << endl;


  delete[] ParityQN;
  delete[] QNold;
  delete[] PosInOld;
  delete[] PosInSingleSite;
  delete[] QNsite;
  delete[] MatchOld;
  delete[] TotalS;
  delete[] QNnew;


}
