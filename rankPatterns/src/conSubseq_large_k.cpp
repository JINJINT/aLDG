#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

void conSubseq(double *x, double *y, int k, int n, double *cnt){
    
    map<double,int> rankMap_x, rankMap_y;
    vector<double> xsub, ysub, temp_x, temp_y;
    vector<double> xv (x,x+n);
    vector<double> yv (y,y+n);
    int temp;
    
    *cnt = 0;
    for (int i=0; i<(n-k+1); i++){
        xsub.insert (xsub.begin(), xv.begin()+i, xv.begin() + i+k);
        ysub.insert (ysub.begin(), yv.begin()+i, yv.begin() + i+k);
        temp_x = xsub;
        temp_y = ysub;
        sort(xsub.begin(),xsub.end());
        sort(ysub.begin(),ysub.end());
        for(int j=0;j<k;j++) {
            rankMap_x[xsub[j]]=j+1;
            rankMap_y[ysub[j]]=j+1;
        }
        xsub = temp_x;
        ysub = temp_y;
        temp = 0;
        for (int l=0;l<k;l++){
            temp = temp + (rankMap_x[xsub[l]] == rankMap_y[ysub[l]]);
        }
        if (temp == k){
            *cnt = *cnt+1;
        }
        xsub.clear();
        ysub.clear();
    }
}
 
