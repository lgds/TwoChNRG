
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
/////////                     //////////////////
///////// hpsort SUBROUTINES  //////////////////
/////////                     //////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/* (C) Copr. 1986-92 Numerical Recipes Software p,{2. */
/* Modified by Luis da Silva to account for vectors on the form ra[0..n] */
/* Modified by Luis da Silva to account for vectors with different types */


/////////               //////////////////
///////// hpsort FLOAT  //////////////////
/////////               //////////////////

#include <stdlib.h>
#include <math.h>


void hpsort_float(unsigned long n, float ra[])
{
    unsigned long i,ir,j,l;
    float rra;

    if (n < 1) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l-1];
        } else {
            rra=ra[ir-1];
            ra[ir-1]=ra[0];
            if (--ir == 1) {
                ra[0]=rra;
                break;
            }
        }
        i=l;
        /*j=l+l;*/
        j=l;
        while (j <= ir) {
            if (j < ir && fabs(ra[j-1]) < fabs(ra[j])) j++;
            if (fabs(rra) < fabs(ra[j-1])) {
                ra[i-1]=ra[j-1];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i-1]=rra;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software p,{2. */
/* Modified by Luis da Silva to account for vectors on the form ra[0..n] */


////////////////////////////////////////////////
////////////////////////////////////////////////

/////////               //////////////////
///////// hpsort INT    //////////////////
/////////               //////////////////

void hpsort_int(unsigned long n, int ra[],int ra2[])
{


    unsigned long i,ir,j,l;
    int rra;
    int rra2;

    if (n < 1) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l-1];
            rra2=ra2[l-1];
        } else {
            rra=ra[ir-1];
            ra[ir-1]=ra[0];

            rra2=ra2[ir-1];
            ra2[ir-1]=ra2[0];

            if (--ir == 1) {
                ra[0]=rra;
                ra2[0]=rra2;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && abs(ra[j-1]) < abs(ra[j]) ) j++;
            if (abs(rra) < abs(ra[j-1])) {
                ra[i-1]=ra[j-1];
                ra2[i-1]=ra2[j-1];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i-1]=rra;
        ra2[i-1]=rra2;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software p,{2. */
/* Modified by Luis da Silva to account for vectors on the form ra[0..n] */

/////////                  //////////////////
/////////  hpsort DOUBLE   //////////////////
/////////                  //////////////////


void hpsort_double(unsigned long n, double ra[])
{
    unsigned long i,ir,j,l;
    double rra;

    if (n < 1) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l-1];
        } else {
            rra=ra[ir-1];

            ra[ir-1]=ra[0];
            if (--ir == 1) {
                ra[0]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        /*j=l; Jan 07. l+l seems right. Dont rememeber why I had j=l*/
        while (j <= ir) {
            if (j < ir && ra[j-1] < ra[j]) j++;
            if (rra < ra[j-1]) {
                ra[i-1]=ra[j-1];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i-1]=rra;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software p,{2. */
/* Modified by Luis da Silva to account for vectors on the form ra[0..n] */

////////////////////////////////////////////////
////////////////////////////////////////////////


/////////                         //////////////////
/////////  hpsort DOUBLE-DOUBLE   //////////////////
/////////                         //////////////////


void hpsort_doubledouble(unsigned long n, double ra[], double ra2[])
{
    unsigned long i,ir,j,l;
    double rra;
    double rra2;

    if (n < 1) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l-1];
            rra2=ra2[l-1];
        } else {
            rra=ra[ir-1];
            rra2=ra2[ir-1];

            ra[ir-1]=ra[0];
            ra2[ir-1]=ra2[0];
            if (--ir == 1) {
                ra[0]=rra;
                ra2[0]=rra2;
                break;
            }
        }
        i=l;
        j=l+l;
        /*j=l; has to be l+l otherwise it does not work*/
        while (j <= ir) {
            if (j < ir && ra[j-1] < ra[j]) j++;
            if (rra < ra[j-1]) {
                ra[i-1]=ra[j-1];
                ra2[i-1]=ra2[j-1];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i-1]=rra;
        ra2[i-1]=rra2;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software p,{2. */
/* Modified by Luis da Silva to account for vectors on the form ra[0..n-1] */

////////////////////////////////////////////////
////////////////////////////////////////////////


/////////                                     //////////////////
/////////  hpsort DOUBLE-DOUBLE matrix[j]     //////////////////
/////////                                     //////////////////


void hpsort_double_doublematrix(unsigned long n, unsigned long nmat, double ra[], double **ra2)
/* Sorts a double vector ra of size n. */
/* In the process, re-shuffles a matrix ra2 (n lines x nmat columns) in order to keep the correspondence ra[i]<->ra2[i][j] c */
{
    unsigned long i,ir,j,l;
    unsigned long iaux;

    double rra;
    double rra2[nmat];

    if (n < 1) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l-1];
            for(iaux=0;iaux<nmat;iaux++){rra2[iaux]=ra2[l-1][iaux];}
/////         rra2=ra2[l-1];
        } else {
            rra=ra[ir-1];
            for(iaux=0;iaux<nmat;iaux++){rra2[iaux]=ra2[ir-1][iaux];}
/////            rra2=ra2[ir-1];

            ra[ir-1]=ra[0];
            for(iaux=0;iaux<nmat;iaux++){ra2[ir-1][iaux]=ra2[0][iaux];}
/////            ra2[ir-1]=ra2[0];
            if (--ir == 1) {
                ra[0]=rra;
                for(iaux=0;iaux<nmat;iaux++){ra2[0][iaux]=rra2[iaux];}
/////                ra2[0]=rra2;
                break;
            }
        }
        i=l;
        /*j=l+l;*/
        j=l;
        while (j <= ir) {
            if (j < ir && ra[j-1] < ra[j]) j++;
            if (rra < ra[j-1]) {
                ra[i-1]=ra[j-1];
                for(iaux=0;iaux<nmat;iaux++){ra2[i-1][iaux]=ra2[j-1][iaux];}
/////                ra2[i-1]=ra2[j-1];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i-1]=rra;
        for(iaux=0;iaux<nmat;iaux++){ra2[i-1][iaux]=rra2[iaux];}
/////        ra2[i-1]=rra2;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software p,{2. */
/* Modified by Luis da Silva to account for vectors on the form ra[0..n] */

////////////////////////////////////////////////
////////////////////////////////////////////////
/////////                         //////////////////
/////////  hpsort DOUBLE-DOUBLE   //////////////////
/////////                         //////////////////


void hpsort_tripledouble(unsigned long n, double ra[], double ra2[], double ra3[])
{
    unsigned long i,ir,j,l;
    double rra;
    double rra2;
    double rra3;

    if (n < 1) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l-1];
            rra2=ra2[l-1];
            rra3=ra3[l-1];
        } else {
            rra=ra[ir-1];
            rra2=ra2[ir-1];
            rra3=ra3[ir-1];

            ra[ir-1]=ra[0];
            ra2[ir-1]=ra2[0];
            ra3[ir-1]=ra3[0];
            if (--ir == 1) {
                ra[0]=rra;
                ra2[0]=rra2;
                ra3[0]=rra3;
                break;
            }
        }
        i=l;
        j=l+l;
        /*j=l; has to be l+l otherwise it does not work*/
        while (j <= ir) {
            if (j < ir && ra[j-1] < ra[j]) j++;
            if (rra < ra[j-1]) {
                ra[i-1]=ra[j-1];
                ra2[i-1]=ra2[j-1];
                ra3[i-1]=ra3[j-1];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i-1]=rra;
        ra2[i-1]=rra2;
        ra3[i-1]=rra3;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software p,{2. */
/* Modified by Luis da Silva to account for vectors on the form ra[0..n-1] */

////////////////////////////////////////////////
////////////////////////////////////////////////
