

double fd_table(int channel, int sigma, int type_i, int type_j){

  // Table for fd_{channel, sigma} Dagger only!
  // type runs from 0 to 15

  double Mat[16][16]={0.0};


  switch (channel)
    {
    case(1):
      // channel 1
      switch (sigma)
	{
	case(1): // fd_(1 up)
	  Mat[4][0]=1.0;
	  Mat[5][1]=1.0;
	  Mat[6][2]=1.0;
	  Mat[7][3]=1.0;

	  Mat[12][8]=1.0;
	  Mat[13][9]=1.0;
	  Mat[14][10]=1.0;
	  Mat[15][11]=1.0;
	  break;
	case(-1): // fd_(1 dn)
	  Mat[8][0]=1.0;
	  Mat[9][1]=1.0;
	  Mat[10][2]=1.0;
	  Mat[11][3]=1.0;

	  Mat[12][4]=-1.0;
	  Mat[13][5]=-1.0;
	  Mat[14][6]=-1.0;
	  Mat[15][7]=-1.0;
	  break;
	default:
	  return(0.0);
	}
      break;
      // channel 2
    case(2):
      switch (sigma)
	{
	case(1): // fd_(2 up)
	  Mat[1][0]=1.0;
	  Mat[3][2]=1.0;
	  Mat[5][4]=-1.0;
	  Mat[7][6]=-1.0;

	  Mat[9][8]=-1.0;
	  Mat[11][10]=-1.0;
	  Mat[13][12]=1.0;
	  Mat[15][14]=1.0;
	  break;
	case(-1): // fd_(2 dn)
	  Mat[2][0]=1.0;
	  Mat[3][1]=-1.0;
	  Mat[6][4]=-1.0;
	  Mat[7][5]=1.0;

	  Mat[10][8]=-1.0;
	  Mat[11][9]=1.0;
	  Mat[14][12]=1.0;
	  Mat[15][13]=-1.0;
	  break;
	default:
	  return(0.0);
	}
      break;
    default:
      return(0.0);
    }
  // END switch ich

  return(Mat[type_i][type_j]);

}
