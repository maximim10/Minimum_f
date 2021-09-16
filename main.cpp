#include <iostream>
#include <vector>
#include <cmath>
const int N=2;
const double E=0.001;
using namespace std;
double func_at_point(vector <double> *point){
    double c=((*point)[0]+1000)*((*point)[0]+1000)+(*point)[1]*(*point)[1];
    return c;
};
void grad_at_point(vector <double> *point, vector <double> *grad){
    vector <double> point_temp1=*point;
    vector <double> point_temp2=*point;
    for (int i=0;i<N;i++){
        point_temp1[i]+=E;
        point_temp2[i]-=E;
        (*grad)[i]=(func_at_point(&point_temp1)-func_at_point(&point_temp2))/(2*E);
        point_temp1=*point;
        point_temp2=*point;
    }
};
void grad2_at_point(vector <double> *point, vector < vector <double> > *grad2){
    vector <double> point_temp11=*point;
    vector <double> point_temp12=*point;
    vector <double> point_temp21=*point;
    vector <double> point_temp22=*point;
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            point_temp11[i]+=E; point_temp11[j]+=E;
            point_temp12[i]+=E; point_temp12[j]-=E;
            point_temp21[i]-=E; point_temp21[j]+=E;
            point_temp22[i]-=E; point_temp22[j]-=E;
            (*grad2)[i][j]=(func_at_point(&point_temp11)-func_at_point(&point_temp21)
                        -func_at_point(&point_temp12)-func_at_point(&point_temp22))/(4*E*E);
            point_temp11=point_temp12=point_temp21=point_temp22=*point;
        }
    }
};
double det(vector < vector <double> > *matr){
    vector < int > perest (N), perest_temp (N);
    int z,t,d=0,d_temp;
    for (int i=0;i<N;i++) { perest[i]=i; }
    while (perest[0]<N){
        perest_temp=perest; z=1; d_temp=1;
        for (int i=0;i<N;i++){
            for (int j=0;j<N-1;j++){
                if (perest_temp[j]==perest_temp[j+1]){
                    z=0;
                    break;
                } else {
                    if (perest_temp[j]>perest_temp[j+1]){
                        t=perest_temp[j];
                        perest_temp[j]=perest_temp[j+1];
                        perest_temp[j+1]=t;
                        if (z==1) { z=-1; } else { z=1; }
                    }
                }
            }
            if (z==0) { break; }
        }
        d_temp*=z;
        for (int i=0;i<N;i++){
            d_temp*=(*matr)[i][perest[i]];
        }
        d+=d_temp;
        perest[N-1]+=1;
        for (int i=N-1;i>0;i--){
            if (perest[i]==N) { perest[i-1]+=1; }
        }
    }
    return d;
}
double Alg_dop(int a, int b, vector < vector <double> > *matr){                                      // Если a, b равны N - имеем определитель
    int M=N-1;
    if (a==N||b==N){ M++; }
    vector < int > perest (M), perest_temp (M);
    int z,t,d=0,d_temp,k,l;
    for (int i=0;i<M;i++) { perest[i]=i; }
    while (perest[0]<M){
        perest_temp=perest; z=1; d_temp=1; k=0; l=0;
        for (int i=0;i<M;i++){
            for (int j=0;j<M-1;j++){
                if (perest_temp[j]==perest_temp[j+1]){
                    z=0;
                    break;
                } else {
                    if (perest_temp[j]>perest_temp[j+1]){
                        t=perest_temp[j];
                        perest_temp[j]=perest_temp[j+1];
                        perest_temp[j+1]=t;
                        if (z==1) { z=-1; } else { z=1; }
                    }
                }
            }
            if (z==0) { break; }
        }
        d_temp*=z;
        k=0;l=0;
        for (int i=0;i<M;i++){
            if (a==i) {k=1;}
            if (perest[i]==b) {l=1;}
            cout<<i+k<<" "<<perest[i]+l<<endl;
            d_temp*=(*matr)[i+k][perest[i]+l];
        }cout<<endl<<endl;
        d+=d_temp;
        perest[M-1]+=1;
        for (int i=M-1;i>0;i--){
            if (perest[i]==M) { perest[i]=0; perest[i-1]+=1; }
        }
    }
    d*=pow(-1,a+b);
    return d;
}
void rev(vector < vector <double> > *matr, vector < vector <double> > *matr_rev){
    double d=0;
    d=Alg_dop(N,N,matr);
    cout<<111;
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){

            (*matr_rev)[i][j]=Alg_dop(i,j,matr)/d;
        }
    }
}
double find_Min (double a, vector <double> *point, vector <double> *grad){

    double Min=func_at_point(point), Min_temp=Min, Min_comp=Min;
    vector <double> point_temp1=(*point), point_temp2=(*point), grad_temp (N);
    grad_at_point(point, grad);
    for (int i=0;i<N;i++) {
        (point_temp1)[i]-=a*(*grad)[i];
    }
    for (int j=0;j*E<2*a;j++){
        for (int i=0;i<N;i++) { (point_temp1)[i]+=E*(*grad)[i]; }
        Min_temp=func_at_point(&point_temp1);
        if (Min_temp<Min){
            Min=Min_temp;
            point_temp2=point_temp1;
        }
    }
    Min_comp=Min_temp;
    for (int i=0;i<N;i++) {
        (point_temp1)[i]-=3*a*(*grad)[i];
    }
    for (int j=0;j*E<4*a;j++){
        for (int i=0;i<N;i++) { (point_temp1)[i]+=E*(*grad)[i]; }
        Min_temp=func_at_point(&point_temp1);
        if (Min_temp<Min_comp){
            Min_comp=Min_temp;
            point_temp2=point_temp1;
        }
    }
    if (Min_comp<Min){
        return find_Min((a+1)*2, point, grad);
    } else {
        *point=point_temp2;
        return Min;
    }
}
double Find_Min_2(vector <double> *point, vector <double> *grad, vector < vector <double> > *grad2){
    grad2_at_point(point, grad2);
    grad_at_point(point, grad);
    vector < vector <double> > grad2_rev=*grad2;
    double a=0;
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){

        }
    }
    vector <double> gr2_m_gr (N);
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            a+=(*grad2)[i][j];
        }
    }
}
int main()
{
    vector <double> point_start (N);
    for (int i=0;i<N;i++) {point_start[i]=3;}
    vector <double> point=point_start;
    vector <double> grad (N);
    vector < vector <double> > grad2 (N);
    vector < vector <double> > matr (N);
    vector < vector <double> > matr_rev (N);
    for (int i=0;i<N;i++) { for (int j=0;j<N;j++) {matr[i].push_back(0);} }
    for (int i=0;i<N;i++) { for (int j=0;j<N;j++) {matr_rev[i].push_back(0);} }

    matr[0][0]=2; matr[1][1]=2; matr[0][1]=2; matr[1][0]=0;
    rev(&matr, &matr_rev);
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            cout<<matr_rev[i][j]<<" ";
        }
        cout<<endl;
    }
    for (int i=0;i<N;i++) { for (int j=0;j<N;j++) {grad2[i].push_back(0);} }
    double a=find_Min(E, &point, &grad),n=1;
    double b=a;
    cout<<a<<" "<<b<<endl;
    for(;;){
        b=find_Min(E, &point, &grad);
        cout<<b<<" "<<n<<endl;
        if (b-a<E) { break; }
        a=b;
    }

    for (int i=0;i<N;i++) { cout<<point[i]<<" "; } cout<<endl;
    cout << "Hello world!" << endl;
    return 0;
}
