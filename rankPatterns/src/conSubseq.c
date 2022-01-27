#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void rank(int*, double*,int);

void conSubseq(double* x, double* y, int k, int n, double* pt) {
    int* rx;
    int* ry;
    int i, j, temp;
    double xsub[k], ysub[k];
    
    *pt = 0;
    for (i=0; i<(n-k+1); i++){
        for (j=i; j<(i+k); j++){
            xsub[j-i] = x[j];
            ysub[j-i] = y[j];
        }
        rx = malloc(k * sizeof(int));
        rank(rx,xsub,k);
        ry = malloc(k * sizeof(int));
        rank(ry,ysub,k);
        temp = 0;
        for (j=0; j<k; j++){
            temp = temp + (rx[j]==ry[j]);
        }
        if (temp == k){
            *pt = *pt + 1;
        }
        free(rx);
        free(ry);
    }
}

void rank(int* r, double* xsub, int k){

    double x;
    int i, j;
    
    for (i=0; i<k; i++){
        x = xsub[i];
        r[i] = 1;
        for (j=0; j<k; j++){
            r[i] += (xsub[j] < x);
        }
    }
}


