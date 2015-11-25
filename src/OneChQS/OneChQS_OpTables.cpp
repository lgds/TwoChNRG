
double fd_table(int sigma, int type_i, int type_j){

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
