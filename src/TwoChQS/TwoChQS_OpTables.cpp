#include <cmath>
#include <iostream>
using namespace std;
//////////////



//////////////
//////////////
//////////////


//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////



double TwoChQS_fd_table(int channel, int sigma, int type_i, int type_j){

  // Table for fd_{channel, sigma} Dagger only!
  // type runs from 0 to 15

  // Table for f^{\dagger}_{sigma}
  // order : |up dn up dn>= fd_(1 up)fd_(1 dn) fd_(2 up)fd_(2 dn)|0>
  //
  // type   States: 
  // 0    -  |0 0> 
  //
  // 1    -  |0 up>
  // 2    -  |up 0> 
  //
  // 3    -  |0 dn>
  // 4    -  |dn 0> 
  //
  // 5    -  1/sqrt(2)(|up dn>-|dn up>)
  // 6    -  |0 up dn>
  // 7    -  |up dn 0> 
  //
  // 8    -  |dn dn>
  //
  // 9    -  1/sqrt(2)(|up dn>+|dn up>)
  //
  // 10    -  |up up>
  //
  // 11    -  |up  up dn>
  // 12    -  |up dn  up>
  //
  // 13    -  |dn  up dn>
  // 14    -  |up dn  dn>
  //
  // 15    -  |up dn  up dn>


  double Mat[16][16]={{0.0}};
  Mat[type_i][type_j]=0.0;

  switch (channel)
    {
    case(1):
      // channel 1
      switch (sigma)
	{
	case(1): // fd_(1 up)
	  Mat[2][0]=1.0;
	  Mat[10][1]=1.0;
	  Mat[9][3]=1.0/sqrt(2.0); // Changes
	  Mat[11][6]=1.0;
	  Mat[5][3]=1.0/sqrt(2.0); // Add

	  Mat[7][4]=1.0;
	  Mat[12][9]=1.0/sqrt(2.0); // Add
	  Mat[12][5]=-1.0/sqrt(2.0); // Changes
	  Mat[14][8]=1.0;
	  Mat[15][13]=1.0;
	  break;
	case(-1): // fd_(1 dn)
	  Mat[4][0]=1.0;
	  Mat[9][1]=1.0/sqrt(2.0); // Add
	  Mat[5][1]=-1.0/sqrt(2.0); // Changes
	  Mat[8][3]=1.0;
	  Mat[13][6]=1.0;

	  Mat[7][2]=-1.0;
	  Mat[12][10]=-1.0;
	  Mat[14][9]=-1.0/sqrt(2.0); // Changes
	  Mat[14][5]=-1.0/sqrt(2.0); // Add
	  Mat[15][11]=-1.0;
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
	  Mat[6][3]=1.0;
	  Mat[10][2]=-1.0;
	  Mat[11][9]=-1.0/sqrt(2.0); // Changes
	  Mat[11][5]=-1.0/sqrt(2.0); // Add

	  Mat[9][4]=-1.0/sqrt(2.0); // Add
	  Mat[5][4]=1.0/sqrt(2.0); // Changes
	  Mat[13][8]=-1.0;
	  Mat[12][7]=1.0;
	  Mat[15][14]=1.0;
	  break;
	case(-1): // fd_(2 dn)
	  Mat[3][0]=1.0;
	  Mat[6][1]=-1.0;
	  Mat[9][2]=-1.0/sqrt(2.0); // Changes
	  Mat[11][10]=1.0;
	  Mat[5][2]=-1.0/sqrt(2.0); // Add

	  Mat[8][4]=-1.0;
	  Mat[13][9]=1.0/sqrt(2.0); // Add
	  Mat[13][5]=-1.0/sqrt(2.0); // Changes
	  Mat[14][7]=1.0;
	  Mat[15][12]=-1.0;
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





//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////




double TwoChQS_SpSm_table(int type_i, int type_j){

  //
  // Non-zero for states connected by S+ + S- within the same channel
  //

  double Mat[16][16];


  Mat[type_i][type_j]=0.0;

  // Diagonal terms are 1, others are zero
  for (int ii=0;ii<16;ii++) 
    Mat[ii][ii]=1.0;

  Mat[1][3]=1.0;
  Mat[3][1]=1.0;

  Mat[2][4]=1.0;
  Mat[4][2]=1.0;

  Mat[11][13]=1.0;
  Mat[13][11]=1.0;

  Mat[12][14]=1.0;
  Mat[14][12]=1.0;

  Mat[8][9]=1.0;
  Mat[9][8]=1.0;

  Mat[9][10]=1.0;
  Mat[10][9]=1.0;

  Mat[8][10]=1.0;
  Mat[10][8]=1.0;


  return(Mat[type_i][type_j]);


}



//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////



double TwoChQS_fdupfdn_table(int channel, int type_i, int type_j){

  //
  // Table for the spin-flip operator f+_{i up} f_{i dn} in channel i
  //

  double Mat[16][16]={0.0};
  Mat[type_i][type_j]=0.0;


  switch (channel)
    {
    case(1):
      // channel 1
      Mat[2][4]=1.0;
      Mat[11][13]=1.0;
      Mat[10][5]=-sqrt(2.0)/2.0;
      Mat[5][8]=sqrt(2.0)/2.0;
      Mat[10][9]=sqrt(2.0)/2.0;
      Mat[9][8]=sqrt(2.0)/2.0;
      break;
    case(2):
      // channel 2
      Mat[1][3]=1.0;
      Mat[12][14]=1.0;
      Mat[10][5]=sqrt(2.0)/2.0;
      Mat[5][8]=-sqrt(2.0)/2.0;
      Mat[10][9]=sqrt(2.0)/2.0;
      Mat[9][8]=sqrt(2.0)/2.0;
      break;
    default:
      return(0.0);
    }

      
  return(Mat[type_i][type_j]);


}


//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////



double TwoChQS_Szch_table(int channel, int type_i){


  double Mat[16]={0.0};
  Mat[type_i]=0.0;


  switch (channel)
    {
    case(1):
      // channel 1
      Mat[2]=1.0;
      Mat[4]=-1.0;
      Mat[8]=-1.0;
      Mat[10]=1.0;
      Mat[11]=1.0;
      Mat[13]=-1.0;
      break;
    case(2):
      // channel 2
      Mat[1]=1.0;
      Mat[3]=-1.0;
      Mat[8]=-1.0;
      Mat[10]=1.0;
      Mat[12]=1.0;
      Mat[14]=-1.0;
      break;
    default:
      return(0.0);
    }
      
  return(0.5*Mat[type_i]);

}

//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////

