#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void incSubseq(double* x, int k, int n, double* pt) {
    double dp[n*k];
    int i, j, p;
    
    for (i=0; i<n; i++){
        dp[i] = 1;
    }
    for (i=n; i<n*k; i++){
        dp[i] = 0;
    }
    
    for (i=1; i<n; i++){
        for (j=0; j<i; j++){
            if (x[i] > x[j]){
                for (p=1; p<k; p++){
                    dp[p*n+i] = dp[(p-1)*n+j] + dp[p*n+i];
                }
            }
        }
    }

    *pt = 0;
    for (i=(k-1)*n; i<k*n; i++){
        *pt = *pt + dp[i];
    }

}


