
int r2ij(int Nel, int r, int &i, int &j){

  i=(int)(r/Nel);
  j=r-Nel*i;

  return(0);

}
