
  // Parameters for command-line passing (GetOpt)
  // To be included in the main code 

int c;
opterr = 0;

//   char ModelOption[32]="Anderson";
//   int ModelNo=0;
//   char ModelSymmetry[32]="OneChQSz";
//   int SymNo=0;
// // new
//   char BandType[32]="SquareWilson";
//   int BandNo=0;

strcpy(ThisCode.Symmetry,"OneChQSz");
strcpy(ThisCode.ModelOption,"Anderson");
strcpy(ThisCode.BandType,"SquareWilson");
ThisCode.ModelNo=0;
ThisCode.SymNo=0;
ThisCode.BandNo=0;

while ((c = getopt (argc, argv, "h ? m: b: D: z:")) != -1)
  switch (c){
  case 'm':
    // 	strcpy(ModelOption,optarg);
    strcpy(ThisCode.ModelOption,optarg);
    break;
  case 'b':
    // 	strcpy(BandType,optarg);
    strcpy(ThisCode.BandType,optarg);
    break;
  case 'h':
    cout << "Usage: " << argv[0] << " (-m model -b band_type -D disc_type)" << endl; 
    cout << " Implemented models/Symmetries: " << endl; 
    cout << "   -m Anderson, 1chQS_Anderson (in |Q S> basis) " << endl; 
    cout << "   -m 1chQSz_Anderson " << endl; 
    cout << "   -m 1chQ_Anderson (in |Q> basis) " << endl; 
    cout << "   -m Kondo, 1chQS_Kondo " << endl; 
    cout << "   -m Phonon, 1chQS_AndersonPhonon (Anderson with phonons) " << endl; 
    cout << "   -m CMphonon, 2chQSP_CMphonon (2 ch: Anderson with CM phonons) " << endl;
    cout << "   -m SMM, 1chQ_SMM  " << endl;
    cout << "   -m DQD, 1chQS_DQD  " << endl;
    cout << "   -m 1chQSz_DQD (with Zeeman in both dots) " << endl;
    cout << "   -m 2chQS_Kondo (implenting) " << endl;
    cout << "   -m 1chNupPdn_Majorana " << endl;
    cout << "   -m 1chS_AndersonSC " << endl;
    cout << "   -m 1chSz_AndersonSC " << endl;
    cout << " Models with chains only (no impurity): " << endl;
    cout << "   -m 1chQS_chain  " << endl;
    cout << "   -m 2chQS_chain  " << endl;
    cout << "   -m 2chQSP_chain  " << endl;
    cout << "   -m 1chQ_chain  " << endl;
    cout << " Band types implemented " << endl;
    cout << "   -b SquareWilson (default): metallic phs square band, no z-trick" << endl; 
    cout << "   -b SideDot : Lorentzian effective HybFunc." << endl;
    cout << "   -b Cavity : many-level cavity attached (effective HybFunc)." << endl;
    cout << "   -b Const : same as square band but using Lanczos procedure with z-trick" << endl;
    cout << "   -b PowerLaw : DoS(w)=C*|w|^r for |w|>Gap, zero otherwise. (needs lanc.in)" << endl; 
    cout << "   -b FromFile : Reads Gamma(w) from file HybFunc.dat (format: w, DoS). Can use z-trick" << endl; 
    cout << "   -b SumDeltas : Reads Gamma(w)=pi.sum_i t^2_i.delta(w-e_i) from file HybDeltas.dat (format: e_i, pi*t(e_i)^2)." << endl; 
    cout << "   -z (z value) : For instance, \"-z 0.5\". Default is z=1.0. " << endl;
    cout << "   NOTE: for z !=1, use \"-b Const\" or \"-b FromFile\" (not SquareWilson)  " << endl;
    cout << " Alternative discretization schemes:  " << endl;
    cout << "   -D (integer value): 0: usual (default) 1: CampoOliveira  " << endl;
    cout << "   NOTE: for alternative schemes, use -b Const or -b FromFile " << endl;
    return 1;
   case 'z':
     //cout << " Warning: z-trick is now in the 2nd line of lanc.in!! Use -b Const" << endl;
     //return 1;
     ThisCode.chain.z_twist=atof(optarg);
//      if (dEqual(ThisCode.chain.z_twist,0.0)) ThisCode.chain.z_twist=1.0;
     if (ThisCode.chain.z_twist<0.0) ThisCode.chain.z_twist=1.0;
     cout << " Using z-trick: z = " << ThisCode.chain.z_twist << endl;
     break;
  case 'D':
    ThisCode.chain.DiscScheme=atof(optarg);
    if (ThisCode.chain.DiscScheme>1) ThisCode.chain.DiscScheme=0;
    cout << " Using Discretization Scheme= " << ThisCode.chain.DiscScheme << endl;
    cout << " 0- usual (default) 1- CampoOliveira  " << endl;
    break;
  case '?':
    cout << "Usage: " << argv[0] << " (-m model -b band_type)" << endl; 
    cout << "  try \"-h\" for more details. " << endl; 
    return 1;
  default:
    return 1;
  }


  cout << " Got Model = " << ThisCode.ModelOption;

// At some point, the number of combinations will be too large 
// and we will have to do this:
//
//  -m model (Kondo, Anderson, DQD, etc.)
//  -s Symmetry (OneChQS, OneChQSz, TwoChQS, etc.)
//  -b Band (Flat, etc.)
//

if ( (strcmp(ThisCode.ModelOption,"Anderson") == 0)||
     (strcmp(ThisCode.ModelOption,"1chQS_Anderson") == 0) )
    {
      ThisCode.ModelNo=0;
      strcpy(ThisCode.Symmetry,"OneChQS");ThisCode.SymNo=0;}

 else
    if (strcmp(ThisCode.ModelOption,"1chQSz_Anderson") == 0)
      {ThisCode.ModelNo=0;
	strcpy(ThisCode.Symmetry,"OneChQSz");ThisCode.SymNo=2;}
    else
      if (strcmp(ThisCode.ModelOption,"1chQ_Anderson") == 0)
	{ThisCode.ModelNo=0;
	  strcpy(ThisCode.Symmetry,"OneChQ");ThisCode.SymNo=3;}
      else
	if ( (strcmp(ThisCode.ModelOption,"Kondo") == 0)||
	     (strcmp(ThisCode.ModelOption,"1chQS_Kondo") == 0) )
	  {ThisCode.ModelNo=1;
	    strcpy(ThisCode.Symmetry,"OneChQS");ThisCode.SymNo=0;}
	else
	  if ( (strcmp(ThisCode.ModelOption,"Phonon") == 0)||
	       (strcmp(ThisCode.ModelOption,"1chQS_AndersonPhonon") == 0) )
	    { ThisCode.ModelNo=2;
	      strcpy(ThisCode.Symmetry,"OneChQS");ThisCode.SymNo=0;}
	  else
	    if ( (strcmp(ThisCode.ModelOption,"Chain") == 0)||
		 (strcmp(ThisCode.ModelOption,"1chQS_chain") == 0) )
	      { ThisCode.ModelNo=3;
		strcpy(ThisCode.Symmetry,"OneChQS");ThisCode.SymNo=0;}
	    else
	      if ( (strcmp(ThisCode.ModelOption,"CMphonon") == 0)||
		   (strcmp(ThisCode.ModelOption,"2chQSP_CMphonon") == 0) )
		{cout << " Got 2ch CM Phonon " << endl; 
		  ThisCode.ModelNo=4;
		  strcpy(ThisCode.Symmetry,"TwoChQSP");ThisCode.SymNo=4;}
	      else
		if ( (strcmp(ThisCode.ModelOption,"SMM") == 0)||
		     (strcmp(ThisCode.ModelOption,"1chQ_SMM") == 0) )
		  {cout << " Got SMM 1ch " << endl; 
		    ThisCode.ModelNo=5;
		    strcpy(ThisCode.Symmetry,"OneChQ");ThisCode.SymNo=3;}
		else
		  if (strcmp(ThisCode.ModelOption,"1chQ_chain") == 0)
		    {cout << " Got 1ch Q Chain " << endl;
		      ThisCode.ModelNo=3;
		      strcpy(ThisCode.Symmetry,"OneChQ");ThisCode.SymNo=3;}
		  else
		    if (strcmp(ThisCode.ModelOption,"2chQS_chain") == 0)
		      {cout << " Got 2ch QS Chain " << endl;
			ThisCode.ModelNo=3;
			strcpy(ThisCode.Symmetry,"TwoChQS");ThisCode.SymNo=1;}
		    else
		      if (strcmp(ThisCode.ModelOption,"2chQSP_chain") == 0)
			{cout << " Got 2ch QSP Chain " << endl;
			  ThisCode.ModelNo=3;
			  strcpy(ThisCode.Symmetry,"TwoChQSP");ThisCode.SymNo=4;}
		      else
			if ( (strcmp(ThisCode.ModelOption,"DQD") == 0)||
			     (strcmp(ThisCode.ModelOption,"1chQS_DQD") == 0) )
			  {cout << " Got 1chQS Double Dot " << endl;
			    ThisCode.ModelNo=6;
			    strcpy(ThisCode.Symmetry,"OneChQS");ThisCode.SymNo=0;}
		      else
			if ( strcmp(ThisCode.ModelOption,"2chQS_Kondo") == 0)
			  {cout << " Got 2chQS Kondo " << endl;
			    ThisCode.ModelNo=1;
			    strcpy(ThisCode.Symmetry,"TwoChQS");
			    ThisCode.SymNo=1;}
		      else
			if ( strcmp(ThisCode.ModelOption,"1chQSz_DQD") == 0)
			  {cout << " Got 1chQSz DQD (Zeeman in both dots) " << endl;
			    ThisCode.ModelNo=6;
			    strcpy(ThisCode.Symmetry,"OneChQSz");
			    ThisCode.SymNo=2;}
		      else
			if ( strcmp(ThisCode.ModelOption,"1chNupPdn_Majorana") == 0)
			  {cout << " Got 1chNupPdn Majorana " << endl;
			    ThisCode.ModelNo=7;
			    strcpy(ThisCode.Symmetry,"OneChNupPdn");
			    ThisCode.SymNo=6;}

		      else
			if ( strcmp(ThisCode.ModelOption,"1chS_AndersonSC") == 0)
			  {cout << " Got 1chS_AndersonSC " << endl;
			    ThisCode.ModelNo=0;
			    strcpy(ThisCode.Symmetry,"OneChS");
			    ThisCode.SymNo=7;}
		      else
			if ( strcmp(ThisCode.ModelOption,"1chSz_AndersonSC") == 0)
			  {cout << " Got 1chSz_AndersonSC " << endl;
			    ThisCode.ModelNo=0;
			    strcpy(ThisCode.Symmetry,"OneChSz");
			    ThisCode.SymNo=8;}
		      else

			{cout << " Invalid model. Exiting... " << endl;
			  exit(0);}


//
// Band type (BandNo follows the convention from old lanc.in files)
//  (4 - side dot; 41 - Cavity, 5 - parallel dot w/ pseudogap, etc...)
//

if (strcmp(ThisCode.BandType,"SquareWilson")==0)
  ThisCode.BandNo=0;   
 else
   if (strcmp(ThisCode.BandType,"SideDot") == 0)
     ThisCode.BandNo=4;
   else
   if (strcmp(ThisCode.BandType,"Cavity") == 0)
     ThisCode.BandNo=41;
   else
     if (strcmp(ThisCode.BandType,"Const") == 0)
       ThisCode.BandNo=11;
     else
       if (strcmp(ThisCode.BandType,"PowerLaw") == 0)
	 ThisCode.BandNo=12;
       else
	 if (strcmp(ThisCode.BandType,"FromFile") == 0)
	   ThisCode.BandNo=13;
	 else
	   if (strcmp(ThisCode.BandType,"SumDeltas") == 0)
	     ThisCode.BandNo=14;
	   else{
	     cout << " Invalid Band type. Exiting... " << endl;
	     exit(0);
	 }


cout << " - Model No: " << ThisCode.ModelNo << endl;
cout << " Symmetry : " << ThisCode.Symmetry << " - SymNo : " << ThisCode.SymNo << endl;
cout << " BandType : " << ThisCode.BandType << " - BandNo : " << ThisCode.BandNo << endl;

// SymNo:
//  0 - OneChQS
//  1 - TwoChQS
//  2 - OneChQSz
//  3 - OneChQ
//  4 - TwoChQSP
//  5 - TwoChQSz
//  6 - OneChNupPdn
//  7 - OneChS
//  8 - OneChSz


// ModelNo:
//  0 - Anderson
//  1 - Kondo
//  2 - Anderson w Holstein phonons
//  3 - Chain only
//  4 - CM phonons (2ch only)
//  5 - SMM (1chQ only)
//  6 - DQD
//  7 - Anderson+Local Majorana (1ch NupPdn only)


// BandNo - BandType (BandNo values come from old lanczosNRG.for code):
//  0 - SquareWilson
//  1* - PowerLaw : w^r 
//  2* - (empty)
//  3* - (empty)
//  4 - SideDot : side-dot case
//  41 - Cavity : many-level cavity.
//  5* - ParallelDot : parallel dot case.
//  6* - ?? : Fully connected dot 
//  7* - ?? :(Case 6 centered at peak) 
//  8* - ?? :(Case 6 centered at dip) 
//  9* - ?? :Parallel dot with a mag field
// 10* - ?? :Semi-elliptical
// 11 - Const : Constant DoS/Constant gamma but uses the Lanczos (e.g., z-trick)
// 12 - PowerLaw : with or without a gap
// 13 - FromFile : Read input file HybFunc.dat
// 14 - SumDeltas : Read input file HybDeltas.dat
// * - Not implemented yet. 
