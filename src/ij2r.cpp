
int ij2r(int Nel, int i, int j){

  // Assumes:
  // - Square Matrix size Nel
  // - If input is 
  //     i=1,Nel and j=1,Nel THEN returns r=1 through Nel*Nel
  //     i=0,Nel-1 and j=0,Nel-1 THEN returns r=0 through Nel*Nel-1
  //

  int r=Nel*i+j;

  return(r);

}

////////////////

int ij2rNSq(int Neli, int Nelj, int i, int j){


  int r=Nelj*i+j;

  //  int r=(Neli <? Nelj)*i+j;


  return(r);

}

////////////////
// Returns position of element aij of a Nel x Nel matrix on
// a vector arranged as
// vec: a00 a01 a10 a02 a20.. a11 a12 a21 ...
// pos:  0   1   2   3   4 ..  
//
int Pmn(int Nel, int i, int j)
{

  int Pmm=0;
  int add=0;

  if (i>=j) //n=i,m=j
    {
      Pmm=(2*Nel-j)*j;
      add=2*(i-j);
    }
  else // n=j,m=i
    {
      Pmm=(2*Nel-i)*i;
      add=2*(j-i)-1;
    }

  return(Pmm+add);
}

////////////////

////////////////
// Returns position of element aij of a the upper triangular 
// part of a Nel x Nel matrix 

int ij2rUpTr(int Nel, int i, int j)
{
  if ((i>j)||(i<0)||(j<0)) return(0);
  else
    return(i*Nel-i*(i-1)/2+(j-i));
}
