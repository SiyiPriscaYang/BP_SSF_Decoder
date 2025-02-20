#ifndef QLDPC_NODE_BASIC_H
#define QLDPC_NODE_BASIC_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"

using namespace std;

class Node{
protected:
    vector<int> neighbors;// store neighboring edge indices
    vector<double> message;// store messages transmitted through the neighboring edges
public:
    Node(){vector<int> neigh(0);vector<double> mes(0);neighbors=neigh;message=mes;}
    Node(vector<int> neigh, vector<double> mes){neighbors=neigh;message=mes;}
    vector<double> get_message(){return this->message;}
    vector<int> get_neighbors(){return this->neighbors;}
    int operator[](int i){return neighbors[i];}
    void set_message(int ind,double m){message[ind]=m;}
    void print(){
        for(int i=0;i<neighbors.size();i++)
            cout<<neighbors[i]<<",";
        cout<<endl;
    }
};

class CN: public Node{
public:
    CN(){vector<int> neigh(0);vector<double> mes(0);neighbors=neigh;message=mes;}
    CN(vector<int> neigh, vector<double> mes){neighbors=neigh;message=mes;}
    vector<bool> l2d();
    void d2l(vector<bool> sign);
    void update_c2v(int syndrome);// update c2v messages
};


// checked
vector<bool> CN::l2d(){
    vector<bool> sign(message.size());
    double llr;
    for(int i=0;i<message.size();i++){
        llr=message[i];
        message[i]=(llr>0)? phi(llr):phi(-llr);
        sign[i]=(llr<0);
    }
    return sign;
}

// checked
void CN::d2l(vector<bool> sign){
    for(int i=0;i<message.size();i++)
        message[i]=sign[i]? -phi(message[i]):phi(message[i]);
}

// checked
void CN::update_c2v(int syndrome) {
    vector<bool> sign=l2d();
    double sum=0;
    bool sign_sum= (syndrome>0);
    for(int i=0;i<message.size();i++){
        sum+=message[i];
        sign_sum=(sign[i])? (!sign_sum):sign_sum;
    }
    for(int i=0;i<message.size();i++){
        message[i]=sum-message[i];
        sign[i]=(sign[i])? (!sign_sum):sign_sum;
    }
    d2l(sign);
}



#endif //QLDPC_NODE_BASIC_H
