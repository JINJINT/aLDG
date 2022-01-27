#include <iostream>
#include <map>
#include <vector>

using namespace std;
void incSubseq(double *x, int k, int n, double *cnt);

extern "C" {
    void incSubseqr(double *x, int *k, int *n, double *cnt){
        incSubseq(x, *k, *n, cnt);
    }
}
