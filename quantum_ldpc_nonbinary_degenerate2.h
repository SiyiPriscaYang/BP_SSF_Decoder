//
// Created by Siyi Yang on 7/20/24.
//

#ifndef QLDPC_DEGENERATE_QUANTUM_LDPC_NONBINARY_DEGENERATE2_H
#define QLDPC_DEGENERATE_QUANTUM_LDPC_NONBINARY_DEGENERATE2_H


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"
#include "node_basic.h"

class qVN: public Node{
protected:
    vector<int> stabilizer;
    vector<vector<double> >  qmessage;
    vector<int> x_neighbors;
    vector<int> z_neighbors;
public:
    qVN(){vector<int> neigh(0);vector<double> mes(0);vector<int> stab(0);vector<vector<double> >  qmes(0);vector<int> x_neigh(0);vector<int> z_neigh(0);neighbors=neigh;message=mes;stabilizer=stab;qmessage=qmes;x_neighbors=x_neigh;z_neighbors=z_neigh;}
    qVN(vector<int> neigh, vector<double> mes,vector<int> stab,vector<vector<double> > qmes,vector<int> x_neigh,vector<int> z_neigh){neighbors=neigh;message=mes;stabilizer=stab;qmessage=qmes;x_neighbors=x_neigh;z_neighbors=z_neigh;}
    void set_qmessage(int ind,double x,double y,double z){qmessage[ind][0]=x;qmessage[ind][1]=y;qmessage[ind][2]=z;}
    vector<int> get_x_neighbors(){return this->x_neighbors;}
    vector<int> get_z_neighbors(){return this->z_neighbors;}
    double update_v2c_x(double lh);
    double update_v2c_z(double lh);
    vector<double> update_v2c(vector<double> lh);// update v2c messages
    void print2(){
        for(int i=0;i<qmessage.size();i++) {
            cout<<i<<':';
            for (int j = 0; j < 3; j++)
                cout << qmessage[i][j] << ',';
            cout<<endl;
        }
    }
};

class qCode{
private:
    vector<qVN> vn;// store VNs
    vector<CN> cn;// store CNs
    vector<vector<int> > vc;//store VN indices adjacent to CNs
    vector<vector<int> > cv;//store CN indices adjacent to VNs
    vector<vector<int> > cv_x;//store X stabilizer indices adjacent to VNs
    vector<vector<int> > cv_z;//store Z stabilizer indices adjacent to VNs
    vector<vector<int> >  sc;//store VN stabilizer values adjacent to CNs
    vector<double> lv2c;// store v2c messages of all edges
    vector<double> lc2v;// store c2v messages of all edges
    vector<vector<int> > bx;
    vector<vector<int> > bz;
    vector<int> idx_cx;
    vector<int> idx_cz;
public:
    qCode(vector<vector<int> >  H);
    void setxz(vector<vector<int> > H);
    vector<int> get_syndrome(vector<int> x);
    bool is_satisfy(vector<int> x,vector<int> syndrome);
    bool satisfy_bin(vector<int> x,vector<int> syndrome,bool x_satisfy,bool z_satisfy);
    bool x_satisfy(vector<int> x,vector<int> syndrome);
    bool z_satisfy(vector<int> x,vector<int> syndrome);
    vector<vector<int> > unsatisfied_CN(vector<int> x,vector<int> syndrome);
    vector<int> unsatisfied_CN_bin(vector<int> x,vector<int> syndrome,bool x_satisfy,bool z_satisfy);
    vector<int> BF(vector<int> uscn,bool is_x, vector<int> vn_flag);
    vector<int> SI(vector<int> uscn, bool is_x);
    vector<int> SP(vector<int> uscn, bool is_x, vector<int> vn_flag);
    void update_c2v(vector<int> syndrome,bool satisfy_x,bool satisfy_z);
    vector<vector<double> > update_v2c(vector<vector<double> >  lch);
    vector<double> update_v2c_x(vector<double> lch);
    vector<double> update_v2c_z(vector<double> lch);
    vector<vector<int> > decode(vector<int> x, vector<vector<double> >  lch, double p, int num_iter, double th);
    bool degenerate(vector<int> x);
    void print(){
        for(int i=0;i<cn.size();i++) {
            cout<<"cn-"<<i<<':';
            cn[i].print();
        }
        for(int i=0;i<vn.size();i++) {
            cout<<"vn-"<<i<<':';
            vn[i].print();
        }
    }
};

// checked
double qVN::update_v2c_x(double lh) {
    double sum=lh;
    for(int i=0;i<x_neighbors.size();i++)
        sum+=message[i];
    for(int i=0;i<x_neighbors.size();i++)
        message[i]=sum-message[i];
    return sum;
}

double qVN::update_v2c_z(double lh) {
    double sum=lh;
    for(int i=0;i<z_neighbors.size();i++)
        sum+=message[i];
    for(int i=0;i<z_neighbors.size();i++)
        message[i]=sum-message[i];
    return sum;
}

// checked
vector<double> qVN::update_v2c(vector<double> lh) {
    vector<double> sum(lh);
    int i,j;
    double a;
    for(i=0;i<message.size();i++) {
        a=message[i];
        if(stabilizer[i]==1)
            set_qmessage(i,0,-a,-a);// a=log(p0/p1)=log(pi+px/py+pz)=log(pi/py,pz)=log(px/py,pz)
        if(stabilizer[i]==2)
            set_qmessage(i,-a,0,-a);// a=log(p0/p1)=log(pi+py/px+pz)=log(pi/px,pz)=log(py/px,pz)
        if(stabilizer[i]==3)
            set_qmessage(i,-a,-a,0);
        for (j = 0; j < 3; j++)
            sum[j] += qmessage[i][j];
        //cout<<"stabilizer="<<stabilizer[i]<<','<<a<<':';
        //print_double(sum);
    }

    for(i=0;i<message.size();i++) {
        for (j = 0; j < 3; j++)
            qmessage[i][j] = sum[j] - qmessage[i][j];
        if(stabilizer[i]==1)
            a=q2b(0,qmessage[i][0])-q2b(qmessage[i][1],qmessage[i][2]);
        if(stabilizer[i]==2)
            a=q2b(0,qmessage[i][1])-q2b(qmessage[i][2],qmessage[i][0]);
        if(stabilizer[i]==3)
            a=q2b(0,qmessage[i][2])-q2b(qmessage[i][0],qmessage[i][1]);
        set_message(i,a);
    }
    return sum;
}


void qCode::setxz(vector<vector<int> > H) {
    int m=H.size();
    int n=H[0].size();
    int t;
    vector<vector<int> > Hx(m, vector<int>(n, 0));
    vector<vector<int> > Hz(m, vector<int>(n, 0));
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            t=H[i][j];
            if(t){
                Hx[i][j]=((t==1)||(t==2))? 1:0;
                Hz[i][j]=((t==2)||(t==3))? 1:0;
            }
        }
    }

    bx= generator(Hx);
    bz= generator(Hz);
}

bool qCode::degenerate(vector<int> x) {
    int len=x.size();
    vector<int> vx(len);
    vector<int> vz(len);
    int sys=bx.size();
    vector<int> sysx(0);
    vector<int> sysz(0);
    int t;

    for(int i=0;i<len;i++){
        t=x[i];
        vx[i]=((t==1)||(t==2))? 1:0;
        vz[i]=((t==2)||(t==3))? 1:0;
    }

    vector<int> idx= is_in_row_space(bx,vx);
    if(idx.size())
        if(idx[0]==-1)
            return false;
    vector<int> idz= is_in_row_space(bz,vz);
    if(idz.size())
        if(idz[0]==-1)
            return false;
/*
    vector<int> idx1, idz1;
    set_difference(idx.begin(), idx.end(), idz.begin(), idz.end(), back_inserter(idx1));
    set_difference(idz.begin(), idz.end(), idx.begin(), idx.end(), back_inserter(idz1));

    for(int i: idx1){
        for(int j=0;j<len;j++){
            if(bx[i][j]&&(!bz[i][j]))
                return false;
        }
    }
    for(int i: idz1){
        for(int j=0;j<len;j++){
            if(bz[i][j]&&(!bx[i][j]))
                return false;
        }
    }
*/

    return true;
}


qCode::qCode(vector<vector<int> >  H){
    int num_cn=H.size();
    int num_vn=H[0].size();
    int num_e=0;
    int e=0;
    vector<vector<int> >  Ht=H;
    vector<vector<int> >  vct(num_cn);
    vector<vector<int> >  cvt(num_vn);
    vector<vector<int> >  cvt_x(num_vn);
    vector<vector<int> >  cvt_z(num_vn);
    vector<vector<int> >  sct(num_cn);
    vector<qVN> vt(num_vn);
    vector<CN> ct(num_cn);
    vector<int> idx_cx_t(0);
    vector<int> idx_cz_t(0);
    bool xz_allocate;
    int t;
    for(int c=0;c<num_cn;c++){
        vector<int> neigh(0);// neigh_cn in e
        vector<int> neigh_c(0);// neigh_cn in v
        vector<int> neigh_sc(0);// neigh_cn in stabilizer
        xz_allocate=false;
        for(int v=0;v<num_vn;v++) {
            if (H[c][v]) {
                if(!xz_allocate){
                    if(H[c][v]==1)
                        idx_cx_t.push_back(c);
                    else
                        idx_cz_t.push_back(c);
                    xz_allocate= true;
                }
                neigh.push_back(num_e);// e
                neigh_c.push_back(v);
                neigh_sc.push_back(H[c][v]);
                Ht[c][v] = num_e;// Ht records edge number
                num_e++;
                //cout<<H[c][v]<<',';
            }
        }
        //cout<<endl;
        vct[c]=neigh_c;// in v
        sct[c]=neigh_sc;// in stabilizer value
        vector<double> mes(neigh.size());
        CN node(neigh,mes);// check
        ct[c]=node;
    }

    for(int v=0;v<num_vn;v++){
        vector<int> neigh(0);//neigh in e
        vector<int> x_neigh(0);
        vector<int> z_neigh(0);
        vector<int> neigh_v(0);// neigh_vn in v
        vector<int> neigh_v_x(0);// neigh_vn in v
        vector<int> neigh_v_z(0);// neigh_vn in v
        vector<int> stab(0);
        for(int c=0;c<num_cn;c++) {
            t=H[c][v];
            if (t) {
                e = Ht[c][v];
                neigh.push_back(e);
                if(t==1) {
                    x_neigh.push_back(e);
                    neigh_v_x.push_back(c);
                }
                else {
                    z_neigh.push_back(e);
                    neigh_v_z.push_back(c);
                }
                neigh_v.push_back(c);
                stab.push_back(t);
            }
        }
        cvt[v]=neigh_v;
        cvt_x[v]=neigh_v_x;
        cvt_z[v]=neigh_v_z;
        vector<double> mes(neigh.size());
        vector<vector<double> >  qmes(neigh.size());
        vector<double> qmest(3);
        for(int i=0;i<neigh.size();i++)
            qmes[i]=qmest;
        qVN node(neigh,mes,stab,qmes,x_neigh,z_neigh);
        vt[v]=node;
    }

    vector<double> c2v(num_e);
    vector<double> v2c(num_e);

    vn=vt;
    cn=ct;
    lv2c=v2c;
    lc2v=c2v;
    vc=vct;
    cv=cvt;
    cv_x=cvt_x;
    cv_z=cvt_z;
    sc=sct;
    idx_cx=idx_cx_t;
    idx_cz=idx_cz_t;
    setxz(H);
}


// checked
vector<int> qCode::get_syndrome(vector<int> x) {
    vector<int> syndrome(vc.size());
    vector<int> lv,ls;
    int parity,ind,t,s;
    for(int c=0;c<vc.size();c++){
        parity=0;
        lv=vc[c];
        ls=sc[c];
        for(int v=0;v<lv.size();v++) {
            ind=lv[v];
            s=ls[v];// edge pauli
            t=x[ind];// vn pauli
            if(t&&(t!=s))
                parity = parity ^ 1;
            //if(t)
            //    cout<<v<<':'<<s<<','<<x[ind]<<';'<<"parity="<<parity<<',';
        }
        syndrome[c]=parity;
        //if(syndrome[c])
        //    cout<<"syndrome["<<c<<"]="<<parity<<';';
    }
    //cout<<endl;
    return syndrome;
}

// checked
bool qCode::is_satisfy(vector<int> x,vector<int> syndrome) {
    vector<int> lv,ls;
    int parity,t,s,c,v;
    for(c=0;c<vc.size();c++){
        parity=0;
        lv=vc[c];
        ls=sc[c];
        for(v=0;v<lv.size();v++) {
            s=ls[v];
            t=x[lv[v]];
            if(t&&(t!=s))
                parity = parity ^ 1;
        }
        if(parity!=syndrome[c])
            return false;
    }
    return true;
}


bool qCode::satisfy_bin(vector<int> x,vector<int> syndrome,bool x_satisfy,bool z_satisfy) {
    vector<int> lv;
    int parity,t,c,k,v;
    vector<int> t_cn;
    if(x_satisfy)//update Z
        t_cn=idx_cz;
    else//update X
        t_cn=idx_cx;

    for(k=0;k<t_cn.size();k++){
        c=t_cn[k];
        parity=0;
        lv=vc[c];
        for(v=0;v<lv.size();v++) {
            t=x[lv[v]];
            if (t)
                parity = parity ^ 1;
        }
        if(parity!=syndrome[c])
            return false;
    }
    return true;
}

bool qCode::x_satisfy(vector<int> x,vector<int> syndrome) {
    vector<int> lv;
    int parity,t,c,k,v;
    for(k=0;k<idx_cx.size();k++){
        c=idx_cx[k];
        parity=0;
        lv=vc[c];
        for(v=0;v<lv.size();v++) {
            t=x[lv[v]];
            if(t&&(t!=1))// t=2,3 (Y,Z)
                parity = parity ^ 1;
        }
        if(parity!=syndrome[c])
            return false;
    }
    return true;
}

bool qCode::z_satisfy(vector<int> x,vector<int> syndrome) {
    vector<int> lv;
    int parity,t,c,k,v;
    for(k=0;k<idx_cz.size();k++){
        c=idx_cz[k];
        parity=0;
        lv=vc[c];
        for(v=0;v<lv.size();v++) {
            t=x[lv[v]];
            if(t&&(t!=3))// t=1,2 (X,Y)
                parity = parity ^ 1;
        }
        if(parity!=syndrome[c])
            return false;
    }
    return true;
}

//only works for CSS codes
vector<vector<int> > qCode::unsatisfied_CN(vector<int> x,vector<int> syndrome) {
    vector<vector<int> > uscn(2);
    vector<int> uscn_x(0);
    vector<int> uscn_z(0);
    vector<int> lv,ls;
    int parity,t,s,c,v;
    for(c=0;c<vc.size();c++){
        parity=0;
        lv=vc[c];
        ls=sc[c];
        for(v=0;v<lv.size();v++) {
            s=ls[v];
            t=x[lv[v]];
            if(t&&(t!=s))
                parity = parity ^ 1;
        }
        if(parity!=syndrome[c]){
            if(s==1)//assume CSS code, then all elements in ls are identical
                uscn_x.push_back(c);
            else
                uscn_z.push_back(c);
        }
    }

    uscn[0]=uscn_x;
    uscn[1]=uscn_z;
    return uscn;
}

// one type of stabilizers all satisfied, x=xt_bin
vector<int> qCode::unsatisfied_CN_bin(vector<int> x,vector<int> syndrome,bool x_satisfy,bool z_satisfy){
    int num_vn=x.size();
    vector<int> uscn(0);
    vector<int> t_cn;
    vector<int> lv;
    int parity,k,c,v;

    if(x_satisfy)// X stabilizers all satisfied, find unsatisfied Z stabilizers
        t_cn=idx_cz;
    else
        t_cn=idx_cx;

    for(k=0;k<t_cn.size();k++){
        c=t_cn[k];
        parity=0;
        lv=vc[c];
        for(v=0;v<lv.size();v++){
            if(x[lv[v]])
                parity=parity^1;
        }
        if(parity!=syndrome[c])
            uscn.push_back(c);
    }

    return uscn;
}

// is_x=true (satisfy_z): uscn= unsatisfied X stabilizers
// is_x=false (satisfy_z): uscn= unsatisfied Z stabilizers
vector<int> qCode::BF(vector<int> uscn, bool is_x, vector<int> vn_flag){
    vector<int> flip_vn(0);
    vector<int> lv;
    int k,c,t,v;

    // X(Z) neighbor of VNs
    vector<vector<int> > t_cv;
    vector<int> t_cv_p;
    if(is_x) {// update X stabilizer
        t_cv = cv_x;// X stabilizers adjacent to VNs
        t_cv_p = idx_cz;// Z stabilizers
    }
    else {
        t_cv = cv_z;
        t_cv_p = idx_cx;
    }

    // flip_flag records the number of unsatisfied X(Z) CNs adjacent to each VN
    vector<int> flip_flag(t_cv.size());
    for(k=0;k<uscn.size();k++){
        c=uscn[k];
        lv=vc[c];//VN neighbors of unsatisfied CN
        for(t=0;t<lv.size();t++){
            v=lv[t];
            flip_flag[v]++;
        }
    }

    // find the VNs that there are more adjacent uscns than scns
    double maxd=0;
    int maxt=0;
    double temp=0;
    int v_opt=0;
    for(v=0;v<t_cv.size();v++){
        k = t_cv[v].size();
        t = flip_flag[v];
        temp = ((double) t / (double) k);
        if (temp > 0.5)
            if (t > maxt) {
                v_opt = v;
                maxd = temp;
                maxt = t;
            }
    }
    if(maxt) {
        flip_vn.push_back(v_opt);
        return flip_vn;
    }

    // find the Z(X) stabilizer that is mostly likely to be flipped
    int max=0;
    int deg=0;
    for(k=0;k<t_cv_p.size();k++){// loop through Z(X) stabilizers
        c=t_cv_p[k];
        lv=vc[c];// the VN neighbors of the Z(X) stabilizer
        deg=0;// the number of VNs that are connected to unsatisfied X(Z) CNs adjacent to Z(X) stabilizers
        for(t=0;t<lv.size();t++){
            v=lv[t];// loop through the VN neighbors of the Z(X) stabilizer
            if(flip_flag[v])
                deg++;
        }
        if((2*deg)>lv.size()) {
            if (deg > max) {
                max = deg;
                flip_vn = lv;
            }
        }
    }

    return flip_vn;
}

vector<int> qCode::SI(vector<int> uscn, bool is_x){
    vector<int> flip_vn(0);
    vector<int> lv;
    int k,c,t,v;

    // X(Z) neighbor of VNs
    vector<vector<int> > t_cv;
    vector<int> t_cv_p;
    if(is_x) {// update X stabilizer
        t_cv = cv_x;// X stabilizers adjacent to VNs
        t_cv_p = idx_cz;// Z stabilizers
    }
    else {
        t_cv = cv_z;
        t_cv_p = idx_cx;
    }

    // flip_flag records the number of unsatisfied X(Z) CNs adjacent to each VN
    vector<int> flip_flag(t_cv.size());
    for(k=0;k<uscn.size();k++){
        c=uscn[k];
        lv=vc[c];//VN neighbors of unsatisfied CN
        for(t=0;t<lv.size();t++){
            v=lv[t];
            flip_flag[v]++;
        }
    }

    // find the Z(X) stabilizer that is mostly likely to be flipped
    int max=0;
    int deg=0;
    for(k=0;k<t_cv_p.size();k++){// loop through Z(X) stabilizers
        c=t_cv_p[k];
        lv=vc[c];// the VN neighbors of the Z(X) stabilizer
        deg=0;// the number of VNs that are connected to unsatisfied X(Z) CNs adjacent to Z(X) stabilizers
        for(t=0;t<lv.size();t++){
            v=lv[t];// loop through the VN neighbors of the Z(X) stabilizer
            if(flip_flag[v])
                deg++;
        }
        if(deg>max)
            if((2*deg)>lv.size()){
                max=deg;
                flip_vn=lv;
            }
    }

    return flip_vn;
}

vector<int> qCode::SP(vector<int> uscn, bool is_x, vector<int> vn_flag){
    vector<int> flip_vn(0);
    int num_vn=vn.size();
    int num_cn=cn.size();
    vector<int> dis_vn(num_vn);
    vector<int> dis_cn(num_cn);
    vector<int> nb_vn(num_vn);
    vector<int> nb_cn(num_cn);
    int i,c,v,t;
    vector<vector<int> > sp_vn(0);
    vector<vector<int> > sp_cn(0);
    vector<int> lv(0);
    vector<int> lc(0);
    vector<vector<int> > cv_t(0);
    if(is_x)
        cv_t=cv_x;
    else
        cv_t=cv_z;

    vector<int> new_vn(0);
    vector<int> new_cn(0);
    vector<int> new_vn_nb(0);
    vector<int> new_cn_nb(0);

    vector<int> zerovector(0);

    for(i=0;i<uscn.size();i++)
        dis_cn[uscn[i]]=1;
    new_cn=uscn;
    sp_cn.push_back(new_cn);

    int s_dis=1;
    bool find_sp=false;
    bool find_vn=false;
    bool find_cn=false;
    int mid_vn=0;
    int mid_cn=0;
    int nb_vn1=0;
    int nb_vn2=0;
    int nb_cn1=0;
    int nb_cn2=0;
    vector<int> t_cn(0);
    vector<int> t_vn(0);

    //cout<<"USCN: ";
    //print_int(uscn);
    //cn-vn-cn-vn-cn-vn-cn
    //1-1-2-2-2-1-1
    //cn-vn-cn-vn-cn-vn-cn-vn-cn
    //1-1-2-2-3-2-2-1-1
    do{
        t_cn=sp_cn[s_dis-1];
        new_vn=zerovector;
        new_cn=zerovector;
        new_vn_nb=zerovector;
        new_cn_nb=zerovector;
        for(i=0;i<t_cn.size();i++){
            c=t_cn[i];
            lv=vc[c];
            //cout<<"neighbors of CN "<<c<<": ";
            //print_int(lv);
            for(t=0;t<lv.size();t++){
                v=lv[t];
                if(vn_flag[v]<2) {
                    if (!dis_vn[v]) {
                        dis_vn[v] = s_dis;
                        nb_vn[v] = c;
                        new_vn.push_back(v);
                    } else if (dis_vn[v] == s_dis) {
                        find_sp = true;
                        find_vn = true;
                        mid_vn = v;
                        nb_vn1 = c;
                        nb_vn2 = nb_vn[v];
                        //cout << "find sp: " << mid_vn << ", " << nb_vn1 << ", " << nb_vn2 << endl;
                        break;
                    }
                }
            }
            if(find_sp)
                break;
        }
        if(find_sp)
            break;
        sp_vn.push_back(new_vn);
        //cout<<"new_vn: ";
        //print_int(new_vn);

        t_vn=sp_vn[s_dis-1];
        for(i=0;i<t_vn.size();i++){
            v=t_vn[i];
            lc=cv_t[v];
            //cout<<"neighbors of VN "<<v<<": ";
            //print_int(lc);
            for(t=0;t<lc.size();t++){
                c=lc[t];
                if(!dis_cn[c]) {
                    dis_cn[c] = s_dis+1;
                    nb_cn[c] = v;
                    new_cn.push_back(c);
                }
                else if(dis_cn[c]==(s_dis+1)){
                    find_sp=true;
                    find_cn=true;
                    mid_cn=c;
                    nb_cn1=v;
                    nb_cn2=nb_cn[c];
                    //cout<<"find sp: "<<mid_cn<<", "<<nb_cn1<<", "<<nb_cn2<<endl;
                    break;
                }
            }
            if(find_sp)
                break;
        }
        if(find_sp)
            break;
        sp_cn.push_back(new_cn);
        //cout<<"new_cn: ";
        //print_int(new_cn);
        s_dis++;
    }
    while(!find_sp);

    if(find_vn){
        flip_vn.push_back(mid_vn);
        s_dis--;
        while(s_dis>0){
            nb_cn1=nb_cn[nb_vn1];
            nb_cn2=nb_cn[nb_vn2];
            flip_vn.push_back(nb_cn1);
            flip_vn.push_back(nb_cn2);
            nb_vn1=nb_vn[nb_cn1];
            nb_vn2=nb_vn[nb_cn2];
            s_dis--;
        }
    }
    else if(find_cn){
        while(s_dis>0){
            flip_vn.push_back(nb_cn1);
            flip_vn.push_back(nb_cn2);
            nb_vn1=nb_vn[nb_cn1];
            nb_vn2=nb_vn[nb_cn2];
            nb_cn1=nb_cn[nb_vn1];
            nb_cn2=nb_cn[nb_vn2];
            s_dis--;
        }
    }

    return flip_vn;
}


// checked
vector<vector<double> >  qCode::update_v2c(vector<vector<double> >  lch) {
    vector<double> tv2c;
    vector<vector<double> >  lh(vn.size());
    vector<int> lc;
    for(int n=0;n<vn.size();n++){
        lc=vn[n].get_neighbors();// vn neighbor in e
        for(int c=0;c<lc.size();c++)
            vn[n].set_message(c,lc2v[lc[c]]);
        lh[n]=vn[n].update_v2c(lch[n]);// already changed mes
        tv2c=vn[n].get_message();// should change to get qmes
        for(int c=0;c<lc.size();c++)
            lv2c[lc[c]]=tv2c[c];
    }
    return lh;
}

// checked
vector<double> qCode::update_v2c_x(vector<double> lch) {
    vector<double> tv2c;
    vector<double> lh(vn.size());
    vector<int> lc;
    for(int n=0;n<vn.size();n++){
        lc=vn[n].get_x_neighbors();// vn neighbor in e
        for(int c=0;c<lc.size();c++)
            vn[n].set_message(c,lc2v[lc[c]]);
        lh[n]=vn[n].update_v2c_x(lch[n]);// already changed mes
        tv2c=vn[n].get_message();// should change to get qmes
        for(int c=0;c<lc.size();c++)
            lv2c[lc[c]]=tv2c[c];
    }
    return lh;
}

// checked
vector<double> qCode::update_v2c_z(vector<double> lch) {
    vector<double> tv2c;
    vector<double> lh(vn.size());
    vector<int> lc;
    for(int n=0;n<vn.size();n++){
        lc=vn[n].get_z_neighbors();// vn neighbor in e
        for(int c=0;c<lc.size();c++)
            vn[n].set_message(c,lc2v[lc[c]]);
        lh[n]=vn[n].update_v2c_z(lch[n]);// already changed mes
        tv2c=vn[n].get_message();// should change to get qmes
        for(int c=0;c<lc.size();c++)
            lv2c[lc[c]]=tv2c[c];
    }
    return lh;
}

// checked
void qCode::update_c2v(vector<int> syndrome,bool satisfy_x,bool satisfy_z) {
    vector<double> tc2v;
    vector<int> lv;
    int n,v;
    if(satisfy_x||satisfy_z){
        if(satisfy_x){//only update z stabilizers
            for(int k=0;k<idx_cz.size();k++){
                n=idx_cz[k];
                lv = cn[n].get_neighbors();
                for (v = 0; v < lv.size(); v++)
                    cn[n].set_message(v, lv2c[lv[v]]);
                cn[n].update_c2v(syndrome[n]);
                tc2v = cn[n].get_message();
                for (v = 0; v < lv.size(); v++)
                    lc2v[lv[v]] = tc2v[v];
            }
        }
        else{//satisfy_z, only update x stabilizers
            for(int k=0;k<idx_cx.size();k++){
                n=idx_cx[k];
                lv = cn[n].get_neighbors();
                for (v = 0; v < lv.size(); v++)
                    cn[n].set_message(v, lv2c[lv[v]]);
                cn[n].update_c2v(syndrome[n]);
                tc2v = cn[n].get_message();
                for (v = 0; v < lv.size(); v++)
                    lc2v[lv[v]] = tc2v[v];
            }
        }
    }
    else{//update full nonbinary decoder
        for(n=0;n<cn.size();n++) {
            lv = cn[n].get_neighbors();
            for (v = 0; v < lv.size(); v++)
                cn[n].set_message(v, lv2c[lv[v]]);
            cn[n].update_c2v(syndrome[n]);
            tc2v = cn[n].get_message();
            for (v = 0; v < lv.size(); v++)
                lc2v[lv[v]] = tc2v[v];
        }
    }
}


// checked, err{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}}
vector<vector<int> > qCode::decode(vector<int> x, vector<vector<double> >  lch, double p, int num_iter, double th) {
    int num_vn=x.size();
    int num_cn=vc.size();
    vector<vector<int> > err{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
    int l,v,c,t,idx,xzsu,s,sum;
    int xs,xu,zs,zu;

    vector<int> syndrome=get_syndrome(x);
    vector<int> diff(num_vn);

    bool satisfy=false;//all stabilizers satisfied
    bool satisfy_x=false;//x stabilizer all satisfied
    bool satisfy_z=false;//z stabilizer all satisfied
    bool reset_lch=false;//channel llr for the other binary decoder set if X/Z stabilizers are all satisfied
    bool identical;//estimation identical to the vector sent
    bool is_converge=false;//non-binary iteration converges
    bool is_converge_bin=false;//binary iteration converges
    bool suc_decode;
    int update_error=0;//if X(Z) stabilizers are all satisfied, X(Z) error still need to be updated, update_error=1(3)

    vector<int> xt(x.size());//estimation of error vector at each iteration (NB decoder) or Z(X) error vector if X(Z) stabilizers are all satisfied
    vector<int> xt_bin(x.size());//estimation of binary X(Z) error vector in Z(X) stabilizer binary decoder if X(Z) stabilizers are all satisfied
    vector<int> xtt(x.size());//estimation of nonbinary error vector
    vector<int> xd(x.size());
    vector<int> uscn(0);
    vector<int> flip(0);

    vector<double> lch_bin(num_vn);//record the channel llr in binary decoder
    double llr0_bin;
    if(p==0)
        llr0_bin=27.5;
    else
        llr0_bin=log(3*(1-p)/p);

    vector<vector<double> >  llr(num_vn);//record the VN llr in the previous iteration, binary decoder
    vector<vector<double> >  llr_new(num_vn);//record the VN llr in the current iteration, binary decoder
    double llr_temp;
    vector<double>  llr_bin(num_vn);//record the VN llr vector in the previous iteration, non-binary decoder
    vector<double>  llr_bin_new(num_vn);//record the VN llr vector in the current iteration, non-binary decoder
    vector<double> nb_llr_temp(3);
    for(t=0;t<num_vn;t++)
        llr[t]=nb_llr_temp;

    for(l=0;l<num_iter;l++) {
        is_converge = true;
        llr_new = update_v2c(lch);
        update_c2v(syndrome, false,false);
        //update_c2v(syndrome, satisfy_x, satisfy_z);
        for (v = 0; v < num_vn; v++) {
            nb_llr_temp = llr_new[v];
            xt[v] = estimate(nb_llr_temp);
            if (is_converge) {
                for (t = 0; t < 3; t++) {
                    if (abs(nb_llr_temp[t] - llr[v][t]) > th)
                        is_converge = false;
                    if (!is_converge)
                        break;
                }
            }
            llr[v] = nb_llr_temp;
        }
        satisfy_x = x_satisfy(xt, syndrome);
        satisfy_z = z_satisfy(xt, syndrome);
        satisfy=(satisfy_x&&satisfy_z);
        if(satisfy||is_converge)
            break;
    }

    vector<int> sp_flag(num_vn);

    //cout<<"finish 1st round\n";
    vector<int> syndrome_t(0);
    if(!satisfy){
        if(satisfy_x||satisfy_z){
            //obtain new syndrome (extract xt)
            update_error=(satisfy_x)? 1:3;
            syndrome_t = get_syndrome(xt);
            for (c = 0; c < syndrome.size(); c++)
                syndrome[c] ^= syndrome_t[c];

            bool last_flip=true;
            int count=0;
            do {
                //cout<<"Entering BF "<< count<<endl;
                uscn = unsatisfied_CN_bin(xd, syndrome, satisfy_x, satisfy_z);
                if(uscn.empty()) {
                    satisfy=true;
                    break;
                }
                flip = BF(uscn, satisfy_z, sp_flag);
                if (flip.size()) {
                    //print_int(flip);
                    // obtain real xt_bin
                    for (t = 0; t < flip.size(); t++) {
                        v = flip[t];
                        xd[v] ^= 1;
                        sp_flag[v]++;
                    }
                }
                else {
                    //cout<<"Entering SP: "<<count<<endl;
                    flip = SP(uscn, satisfy_z, sp_flag);
                    if (flip.size()) {
                        //print_int(flip);
                        // obtain real xt_bin
                        for (t = 0; t < flip.size(); t++) {
                            v = flip[t];
                            xd[v] ^= 1;
                            sp_flag[v]++;
                        }
                    } else
                        last_flip = false;
                }
                count++;
            }
            while(last_flip&&(count<20));

            //cout<<"finish discrete part\n";


            uscn = unsatisfied_CN_bin(xd, syndrome, satisfy_x, satisfy_z);
            if(uscn.empty())
                satisfy=true;
            else
                satisfy=false;
            xt_bin=xd;
            /*
            if(uscn.empty()){
                satisfy=true;
                xt_bin=xd;
            }
            else{
                //cout<<"entering 2nd round\n";
                satisfy=false;

                xtt=xt;
                syndrome_t = get_syndrome(xd);
                for (c = 0; c < syndrome.size(); c++)
                    syndrome[c] ^= syndrome_t[c];
                for(v=0;v<num_vn;v++)
                    if(xd[v]) {
                        xd[v] = 0;
                        xt[v] = err[xtt[v]][update_error];
                    }


                if (satisfy_x) {//change lch, fix Z errors
                    //update lch and extract Z error vector in xt
                    for (v = 0; v < num_vn; v++) {
                        t = xt[v];
                        if ((t == 0) || (t == 1)) {//I or X
                            //xt[v] = 0;// extract Z error as I
                            lch_bin[v] = llr0_bin;
                            //lch_bin[v] = -llr_new[v][0];//llr(I)-llr(X)
                        } else {// Y or Z
                            //xt[v] = 3;// extract Z error as Z
                            lch_bin[v] = 0.0;
                            //lch_bin[v] = llr_new[v][2] - llr_new[v][1];//llr{Z)-llr(Y)
                        }
                        //xt_bin[v] = 0;
                    }
                }
                if (satisfy_z) {
                    for (v = 0; v < num_vn; v++) {
                        t = xt[v];
                        if ((t == 0) || (t == 3)) {//I or Z
                            //xt[v] = 0;
                            lch_bin[v] = llr0_bin;
                            //lch_bin[v] = -llr_new[v][2];
                        } else {// X or Y
                            //xt[v] = 1;
                            lch_bin[v] = 0.0;
                            //lch_bin[v] = llr_new[v][0] - llr_new[v][1];
                        }
                        //xt_bin[v] = 0;
                    }
                }
                //reset lc2v and lv2c to all zero vector
                for (t = 0; t < lc2v.size(); t++)
                    lc2v[t] = lv2c[t] = 0.0;
                is_converge_bin = false;//restart

                for (l = 0; l < num_iter; l++) {
                    if (satisfy_x)// only update Z stabilizers if X stabilizers are all satisfied
                        llr_bin_new = update_v2c_z(lch_bin);
                    else// only update X stabilizers if Z stabilizers are all satisfied
                        llr_bin_new = update_v2c_x(lch_bin);
                    update_c2v(syndrome, satisfy_x, satisfy_z);
                    is_converge_bin = true;
                    for (v = 0; v < num_vn; v++) {
                        llr_temp = llr_bin_new[v];
                        xt_bin[v] = (llr_temp >= 0) ? 0 : 1;
                        if (abs(llr_temp - llr_bin[v]) > th)
                            is_converge_bin = false;
                        llr_bin[v] = llr_temp;
                    }
                    if (satisfy_x)
                        satisfy_z = satisfy_bin(xt_bin, syndrome, satisfy_x, satisfy_z);
                    else
                        satisfy_x = satisfy_bin(xt_bin, syndrome, satisfy_x, satisfy_z);

                    satisfy = (satisfy_x && satisfy_z);
                    if (satisfy || is_converge_bin)
                        break;
                }
                //cout<<"finish 2nd round\n";
            }
             */
        }
    }

    identical=true;
    xtt=xt;
    if(update_error) {
        // xt+xt_bin*update_error
        for(v=0;v<vn.size();v++){
            if(xt_bin[v])
                xt[v] = err[xtt[v]][update_error];
        }
    }

    for(v=0;v<vn.size();v++){
        if(xt[v]!=x[v]) {
            diff[v] = err[xt[v]][x[v]];
            if(identical)
                identical=false;
        }
        else
            diff[v]=0;
    }

    // find error type
    if(satisfy) {
        if (identical) {
            suc_decode=true;
            //cout << "Decoder succeed\n";
            diff.push_back(0);
        }
        else{
            suc_decode = degenerate(diff);
            if (suc_decode) {
                //cout << "Decoder succeed\n";
                diff.push_back(0);
            }
            else {
                //cout << "Decoder failed: convergent error\n";
                diff.push_back(1);// convergent error
            }
        }
    }
    else{
        suc_decode=false;
        //cout<<"Decoder failed: non-convergent error\n";
        diff.push_back(2);// non-convergent error
    }

    vector<vector<int> > io_pair(2);
    io_pair[0]=diff;
    io_pair[1]=x;

    return io_pair;
}

#endif //QLDPC_DEGENERATE_QUANTUM_LDPC_NONBINARY_DEGENERATE2_H
