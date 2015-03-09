
  // Parameters for command-line passing (GetOpt)
  // To be included in the main code 

  int c;
  opterr = 0;

  char ModelOption[]="Anderson";
  int ModelNo=0;
  char ModelSymmetry[]="OneChQ";


  while ((c = getopt (argc, argv, "h m: S:")) != -1)
    switch (c)
      {
      case 'm':
	strcpy(ModelOption,optarg);
	break;
      case 'S':
	strcpy(ModelSymmetry,optarg);
	break;
      case 'h':
	cout << "Usage: " << argv[0] << " (-m model -S Symmetry)" << endl; 
	cout << " Implemented models: " << endl; 
	cout << "   -m Anderson (default) " << endl; 
	cout << "   -m Kondo " << endl; 
	cout << "   -m Phonon (Anderson with phonons) " << endl; 
	cout << "   -m CMphonon (2 ch: Anderson with CM phonons) " << endl;
	cout << "   -m SMM (1 ch: SMM model) " << endl;
	cout << " Implemented Symmetries: " << endl; 
	cout << "   -S OneChQ (default) " << endl;
	cout << "   -S OneChQS (to be implemented)  " << endl;  
	cout << "   -S TwoChQS (to be implemented)  " << endl;  
	return 1;
       case '?':
 	cout << "Usage: " << argv[0] << " (-m model)" << endl; 
 	cout << "  try \"-h\" for more details. " << endl; 
 	return 1;
      default:
	return 1;
      }

  cout << " Model = " << ModelOption << endl;


  if (strcmp(ModelOption,"Anderson") == 0)
    {cout << " Got Anderson " << endl;ModelNo=0;}
  else
    if (strcmp(ModelOption,"Kondo") == 0)
      {cout << " Got Kondo " << endl;ModelNo=1;}
    else
      if (strcmp(ModelOption,"Phonon") == 0)
	{cout << " Got Phonon " << endl;ModelNo=2;}
      else
	if (strcmp(ModelOption,"Chain") == 0)
	  {cout << " Got Chain " << endl;ModelNo=3;}
	else
	  if (strcmp(ModelOption,"CMphonon") == 0)
	    {cout << " Got CM Phonon " << endl; ModelNo=4;}
          else
	    if (strcmp(ModelOption,"SMM") == 0)
	     {cout << " Got SMM 1ch " << endl; ModelNo=5;
	     strcpy(ModelSymmetry,"OneChQ");}
            else
	      {cout << " Invalid model. Exiting... " << endl;exit(0);}
	


//   if (strcmp(ModelOption,"Anderson") == 0)
//     {
//       cout << " Got Anderson " << endl;
//       ModelNo=0;
//     }
//   else
//     {
//       if (strcmp(ModelOption,"Kondo") == 0)
// 	{
// 	  cout << " Got Kondo " << endl;
// 	  ModelNo=1;
// 	}
//       else
// 	{
// 	  if (strcmp(ModelOption,"Phonon") == 0)
// 	    {
// 	      cout << " Got Phonon " << endl;
// 	      ModelNo=2;
// 	    }
// 	  else
// 	    {
// 	      cout << " Invalid model. Exiting... " << endl;
// 	      exit(0);
// 	    }
// 	}
//     }      

