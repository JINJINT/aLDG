#include <iostream>
#include <map>
#include <vector>
using namespace std;
void conSubseq(double *x, double *y, int k, int n, double *cnt);

extern "C" {
    void conSubseqr(double *x, double *y, int *k, int *n, double *cnt){
        conSubseq(x, y, *k, *n, cnt);
    }
}
