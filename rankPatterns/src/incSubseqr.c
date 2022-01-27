//
//  incSubseqr.c converts all arguments in incSubseq to pointers
//  
//
//  Created by Rachel Wang on 2/24/14.
//
//
void incSubseq(double *, int, int, double *);

void incSubseq_r(double *x, int *k, int *n, double *pt){
    incSubseq(x, *k, *n, pt);
}
