//
// Created by Siyi Yang on 1/5/24.
//

#ifndef QLDPC_CODE_CONSTRUCTION_H
#define QLDPC_CODE_CONSTRUCTION_H


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"
using namespace std;

vector<vector<vector<int> > > GB(int l,int z,vector<int> a,vector<int> b);
// B1[[882,24,18-24]]: l=63,k=7,a=[27,54,0],b=[0,1,6]
// B2[[882,48,16]]: l=63,k=7,a=[27,0,27,18],b=[0,1,6]
vector<vector<vector<int> > > GB(int l,int z,vector<int> a,vector<int> b){
    int n=l*z;
    int sl,s1,s2,s3,s4;
    int i,j,k;

    vector<vector<int> > Hx(n);
    vector<vector<int> > Hz(n);
    for(i=0;i<n;i++){
        vector<int> t1(2*n);
        Hx[i]=t1;
        vector<int> t2(2*n);
        Hz[i]=t2;
    }

    for(i=0;i<a.size();i++){
        sl=a[i];
        for(j=0;j<z;j++){
            s1=((j+i)%z)*l;
            s2=j*l;
            s3=s2;
            s4=n+s1;
            for(k=0;k<l;k++){
                Hx[s1+((k+sl)%l)][s2+k]=1;
                Hz[s3+k][s4+((k+sl)%l)]=1;
            }
        }
    }

    for(i=0;i<b.size();i++){
        sl=b[i];
        for(j=0;j<z;j++){
            s1=s3=s4=j*l;
            s2=n+s1;
            for(k=0;k<l;k++){
                Hx[s1+((k+sl)%l)][s2+k]=1;
                Hz[s3+k][s4+((k+sl)%l)]=1;
            }
        }
    }

    vector<vector<vector<int> > > H(2);
    H[0]=Hx;
    H[1]=Hz;

    return H;
}

#endif //QLDPC_CODE_CONSTRUCTION_H
