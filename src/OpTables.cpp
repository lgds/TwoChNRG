
#include <cmath>
//////////////
//////////////
//////////////


double OneCh_fd_table(int sigma, int type_i, int type_j){

  // Table for f^{\dagger}_{sigma}
  // type   States: 
  // 0    -  |0> 
  // 1    -  |up>=fd_(up)|0>
  // 2    -  |dn>=fd_(dn)|0>
  // 3    -  |up dn>=fd_(up)fd_(dn)|0> (notice order)
  switch (sigma)
    { 
    case (1):
      if ( (type_i==1)&&(type_j==0) ) return(1.0); 
      else
      if ( (type_i==3)&&(type_j==2) ) return(1.0);
      else return(0.0);
      break;
    case (-1):
      if ( (type_i==2)&&(type_j==0) ) return(1.0);
      else
      if ( (type_i==3)&&(type_j==1) ) return(-1.0);
      else return(0.0); 
      break;
    default:
      return(0.0);
    }
}

//////////////
//////////////
//////////////



double TwoCh_fd_table(int channel, int sigma, int type_i, int type_j){

  // Table for fd_{channel, sigma} Dagger only!
  // type runs from 0 to 15

  // Table for f^{\dagger}_{sigma}
  // order : |up dn up dn>= fd_(1 up)fd_(1 dn) fd_(2 up)fd_(2 dn)|0>
  //
  // Oct 11: CHANGED ORDER OF THE STATES!! (I did check and this routine was 
  //  not being called by TwoCh Codes BUT watch out...)
  // Using the order set in TwoChQSz_SetSingleSite

  // type   States: 
  // 0    -  |0 0> 

  // 1    -  |0 up>
  // 2    -  |up 0> 

  // 3    -  |0 dn>
  // 4    -  |dn 0> 

  // 5    -  |up dn>
  // 6    -  |0  updn>
  // 7    -  |updn 0> 
  // 8    -  |dn up>

  // 9    -  |up up>

  // 10    -  |dn dn>

  // 11    -  |up  updn>
  // 12    -  |updn up>

  // 13    -  |dn  updn>
  // 14    -  |updn dn>

  // 15    -  |updn updn>


  double Mat[16][16]={0.0};


  switch (channel)
    {
    case(1):
      // channel 1
      switch (sigma)
	{
	case(1): // fd_(1 up)
	  Mat[2][0]=1.0;
	  Mat[9][1]=1.0;
	  Mat[5][3]=1.0;
	  Mat[7][4]=1.0;

	  Mat[11][6]=1.0;
	  Mat[12][8]=1.0;
	  Mat[14][10]=1.0;
	  Mat[15][13]=1.0;
	  break;
	case(-1): // fd_(1 dn)
	  Mat[4][0]=1.0;
	  Mat[8][1]=1.0;
	  Mat[7][2]=-1.0;
	  Mat[10][3]=1.0;

	  Mat[14][5]=-1.0;
	  Mat[13][6]=1.0;
	  Mat[12][9]=-1.0;
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
	  Mat[9][2]=-1.0;
	  Mat[6][3]=1.0;
	  Mat[8][4]=-1.0;

	  Mat[11][5]=-1.0;
	  Mat[12][7]=1.0;
	  Mat[13][10]=-1.0;
	  Mat[15][14]=1.0;
	  break;
	case(-1): // fd_(2 dn)
	  Mat[3][0]=1.0;
	  Mat[6][1]=-1.0;
	  Mat[5][2]=-1.0;
	  Mat[10][4]=-1.0;

	  Mat[14][7]=1.0;
	  Mat[13][8]=1.0;
	  Mat[11][9]=1.0;
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

////////////////////
////////////////////
////////////////////



double TwoDotQS_fd_reduced_table(int channel, int type_i, int type_j){

  // Table for reduced matrix elements for fd_{channel} Dagger only!
  // type runs from 0 to 9

  // Table for f^{\dagger}:
  // order : |up dn up dn>= fd_(1 up)fd_(1 dn) fd_(2 up)fd_(2 dn)|0>
  //
  // type   States: Q S Sz 
  // 0    -  |-2 0 0 >  - | 0 0>
  // 1    -  |-1 1/2 1/2> - |0 up> 
  // 2    -  |-1 1/2 1/2> - |up 0>

  // not inc (3)   -  |-1 1/2 -1/2> - |0 dn> 
  // not inc (4)    -  |-1 1/2 -1/2> - |dn 0>

  // 3 (5)   -  |0 0 0> - |up dn>-|dn up>/sqrt(2.0)
  // 4 (6)   -  |0 0 0> - |0 updn>
  // 5 (7)    -  |0 0 0> - |updn 0>
  //
  // not inc (8)   -  |0 1 -1> - |dn dn>
  // not inc (9)    -  |0 1 0> - |up dn>+|dn up>/sqrt(2.0)
  // 6 (10)   -  |0 1 1> - |up up>

  // 7 (11)   -  |1 1/2 1/2> - |up updn> 
  // 8 (12)   -  |1 1/2 1/2> - |updn up>

  // not inc (13)    -  |1 1/2 -1/2> - |dn updn> 
  // not inc (14)   -  |1 1/2 -1/2> - |updn dn>

  // 9 (15)    -  |2 0 0> - |updn updn>

  double Mat[10][10]={0.0};
  // Clebsch Gordan coefs:

  // up 
  // <0 0 1/2 1/2|1/2 1/2> = 1.0  
  // <1/2 -1/2 1/2 1/2|0 0> = -1.0/sqrt(2.0)
  // <1 -1 1/2 1/2|1/2 -1/2> = -sqrt(2.0/3.0)
  // <1 0 1/2 1/2|1/2 1/2> = -1.0/sqrt(3.0)
  // <1/2 1/2 1/2 1/2 |1 1> = 1.0
  // <1/2 -1/2 1/2 1/2 |1 0> = 1.0/sqrt(2.0)


  // dn 
  // <0 0 1/2 -1/2|1/2 -1/2> = 1.0
  // <1/2 1/2 1/2 -1/2|0 0> = 1.0/sqrt(2.0)
  // <1 1 1/2 -1/2|1/2 1/2> = sqrt(2.0/3.0)
  // <1 0 1/2 -1/2|1/2 -1/2> = 1.0/sqrt(3.0)
  // <1/2 -1/2 1/2 -1/2 |1 -1> = 1.0
  // <1/2 1/2 1/2 -1/2 |1 0> = 1.0/sqrt(2.0)


  switch (channel){
  case(1):
    // channel/dot 1
    // fd_(1 up)   
    Mat[2][0]=1.0;
    Mat[6][1]=1.0; 
    Mat[8][3]=-1.0/sqrt(2.0);
    Mat[7][4]=1.0;
    // fd_(1 dn)
    Mat[3][1]=-1.0; 
    Mat[5][2]=-sqrt(2.0);
    Mat[8][6]=-sqrt(3.0/2.0);
    Mat[9][7]=-sqrt(2.0);
//     // fd_(1 up)   
//     Mat[2][0]=1.0;
//     Mat[10][1]=1.0; 
//     Mat[5][3]=-1.0;
//     Mat[9][3]=-1.0;
//     Mat[7][4]=-sqrt(2.0);
    
//     Mat[12][5]=-1.0/sqrt(2.0);
//     Mat[11][6]=1.0;
//     Mat[14][8]=-sqrt(3.0/2.0);
//     Mat[12][9]=-sqrt(3.0/2.0);
//     Mat[15][13]=-sqrt(2.0);
//     // fd_(1 dn)
//     Mat[4][0]=1.0;
//     Mat[5][1]=-1.0; 
//     Mat[9][1]=1.0;
//     Mat[7][2]=-sqrt(2.0);
//     Mat[8][3]=1.0;
    
//     Mat[14][5]=-1.0/sqrt(2.0);
//     Mat[13][6]=1.0;
//     Mat[14][9]=-sqrt(3.0/2.0);
//     Mat[12][10]=-sqrt(3.0/2.0);
//     Mat[15][11]=-sqrt(2.0);
    break;
    // channel/dot 2
  case(2):
    // fd_(2 up)
    Mat[1][0]=1.0;
    Mat[6][2]=-1.0; 
    Mat[7][3]=-1.0/sqrt(2.0);
    Mat[8][5]=1.0;
    // fd_(2 dn)
    Mat[4][1]=-sqrt(2.0);
    Mat[3][2]=-1.0; 
    Mat[7][6]=sqrt(3.0/2.0);
    Mat[9][8]=-sqrt(2.0);

//     // fd_(2 up)
//     Mat[1][0]=1.0;
//     Mat[10][2]=-1.0; 
//     Mat[6][3]=-sqrt(2.0);
//     Mat[5][4]=-1.0;
//     Mat[9][4]=1.0;
    
//     Mat[11][5]=-1.0/sqrt(2.0);
//     Mat[12][7]=1.0;
//     Mat[13][8]=sqrt(3.0/2.0);
//     Mat[11][9]=sqrt(3.0/2.0);
//     Mat[15][14]=-sqrt(2.0);
//     // fd_(2 dn)
//     Mat[3][0]=1.0;
//     Mat[6][1]=-sqrt(2.0);
//     Mat[5][2]=-1.0; 
//     Mat[9][2]=-1.0;
//     Mat[8][4]=-1.0;
    
//     Mat[13][5]=-1.0/sqrt(2.0);
//     Mat[14][7]=1.0;
//     Mat[13][9]=sqrt(3.0/2.0);
//     Mat[11][10]=sqrt(3.0/2.0);
//     Mat[15][12]=-sqrt(2.0);
    break;
  default:
    return(0.0);
  }
  // END switch ich

return(Mat[type_i][type_j]);

}

////////////////////

double TwoDotQSSz_fd_table(int channel, int sigma, int type_i, int type_j){

  // Table for reduced matrix elements for fd_{channel} Dagger only!
  // type runs from 0 to 9

  // Table for f^{\dagger}:
  // order : |up dn up dn>= fd_(1 up)fd_(1 dn) fd_(2 up)fd_(2 dn)|0>
  //
  // type   States: Q S Sz 
  // 0    -  |-2 0 0 >  - | 0 0>
  // 1    -  |-1 1/2 1/2> - |0 up> 
  // 2    -  |-1 1/2 1/2> - |up 0>

  // 3    -  |-1 1/2 -1/2> - |0 dn> 
  // 4    -  |-1 1/2 -1/2> - |dn 0>

  // 5    -  |0 0 0> - |up dn>-|dn up>/sqrt(2.0)
  // 6    -  |0 0 0> - |0 updn>
  // 7    -  |0 0 0> - |updn 0>
  //
  // 8    -  |0 1 -1> - |dn dn>
  // 9    -  |0 1 0> - |up dn>+|dn up>/sqrt(2.0)
  // 10   -  |0 1 1> - |up up>

  // 11   -  |1 1/2 1/2> - |up updn> 
  // 12   -  |1 1/2 1/2> - |updn up>

  // 13   -  |1 1/2 -1/2> - |dn updn> 
  // 14   -  |1 1/2 -1/2> - |updn dn>

  // 15   -  |2 0 0> - |updn updn>

  double Mat[16][16]={0.0};

  switch (channel){
  case(1):
    switch (sigma){
    case(1): 
      // fd_(1 up)   
      Mat[2][0]=1.0;
      Mat[10][1]=1.0; 
      Mat[5][3]=1.0/sqrt(2.0);
      Mat[9][3]=1.0/sqrt(2.0);
      Mat[7][4]=1.0;
    
      Mat[12][5]=-1.0/sqrt(2.0);
      Mat[11][6]=1.0;
      Mat[14][8]=1.0;
      Mat[12][9]=1.0/sqrt(2.0);
      Mat[15][13]=1.0;
      break;
//     // fd_(1 dn)
    case(-1): 
      Mat[4][0]=1.0;
      Mat[5][1]=-1.0/sqrt(2.0); 
      Mat[9][1]=1.0/sqrt(2.0);
      Mat[7][2]=-1.0;
      Mat[8][3]=1.0;
    
      Mat[14][5]=-1.0/sqrt(2.0);
      Mat[13][6]=1.0;
      Mat[14][9]=-1.0/sqrt(2.0);
      Mat[12][10]=-1.0;
      Mat[15][11]=-1.0;
      break;
    default:
      return(0.0);
    }
    break;
    // channel/dot 2
  case(2):
    switch (sigma){
    case(1): 
      // fd_(2 up)
      Mat[1][0]=1.0;
      Mat[10][2]=-1.0; 
      Mat[6][3]=1.0;
      Mat[5][4]=1.0/sqrt(2.0);
      Mat[9][4]=-1.0/sqrt(2.0);
    
      Mat[11][5]=-1.0/sqrt(2.0);
      Mat[12][7]=1.0;
      Mat[13][8]=-1.0;
      Mat[11][9]=-1.0/sqrt(2.0);
      Mat[15][14]=1.0;
      break;
    case(-1): 
      // fd_(2 dn)
      Mat[3][0]=1.0;
      Mat[6][1]=-1.0;
      Mat[5][2]=-1.0/sqrt(2.0); 
      Mat[9][2]=-1.0/sqrt(2.0);
      Mat[8][4]=-1.0;
    
      Mat[13][5]=-1.0/sqrt(2.0);
      Mat[14][7]=1.0;
      Mat[13][9]=1.0/sqrt(2.0);
      Mat[11][10]=1.0;
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

////////////////////
