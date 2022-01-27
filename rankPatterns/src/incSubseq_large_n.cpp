#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

const int n_max = 100000;
const int k_max = 200;
typedef double T;

T kBIT[k_max][n_max];
T sub[k_max][n_max];

void update(int idx, int k, int n, T val);
T getSum(int idx, int k);

void incSubseq(double *x, int k, int n, double *cnt){
    map<double,int> rankMap;
        
    memset(kBIT,0,sizeof(kBIT));
    memset(sub,0,sizeof(sub));
    vector<double> xv (x,x+n);
    vector<double> temp = xv;
    sort(xv.begin(),xv.end());
    
    for(int i=0; i<n; i++){
        rankMap[xv[i]]=i+1;
    }
    xv = temp;
    
    *cnt = 0;
    for(int i=0; i<n; i++){
        int index = rankMap[xv[i]];
        sub[1][i]=1;
        update(index,1,n,sub[1][i]);
        for(int j=2;j<=k;j++){
            sub[j][i]=getSum(index-1,j-1);
            update(index,j,n,sub[j][i]);
        }
        *cnt = *cnt + sub[k][i];
    }
}

void update(int idx, int k, int n, T val){
    while(idx<n){
        kBIT[k][idx]+=val;
        idx+=(idx& (-idx));
    }
}

T getSum(int idx, int k){
    T ret = 0;
    while(idx>0) {
        
        ret+=kBIT[k][idx];
        
        idx-=(idx& (-idx));
    }
    return ret;
}

