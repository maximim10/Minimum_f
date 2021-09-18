#include <iostream>
#include <vector>
#include <cmath>
const int N=2;
const double E=0.000001;
using namespace std;
double func_at_point(vector <double> *point){
    double c=pow(1-(*point)[0],2)+10*pow((*point)[1]-pow((*point)[0],2),2)+2;
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
                        -func_at_point(&point_temp12)+func_at_point(&point_temp22))/(4*E*E);
            point_temp11=point_temp12=point_temp21=point_temp22=*point;
        }
    }
};
double Alg_dop(int a, int b, vector < vector <double> > *matr){                                      // If a, b equals N - this is a determinant
    int M=N-1;
    if (a==N||b==N){ M++; }
    vector < int > perest (M), perest_temp (M);
    int z,t,k,l;
    double d=0,d_temp=0;
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
                        z*=-1;
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
            d_temp*=(*matr)[i+k][perest[i]+l];
        }
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
    if (d==0){
        cout<<"E pobolshe nado, grad2 takoe ne perevarivaet"<<endl;
    }
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            (*matr_rev)[i][j]=Alg_dop(i,j,matr)/d;
        }
    }
}
double Find_Min (vector <double> *point, vector <double> *grad){
    double f=func_at_point(point);
    double Min=func_at_point(point), Min_temp=Min, Min_comp=Min;
    vector <double> point_temp1=(*point), point_temp2=(*point), point_temp, grad_temp (N);
    grad_at_point(point, grad);
    double a=1,b=1,m=2,l=0.75,f_prev=f;
    for (int b=0,d=1;b<1;){
        point_temp=(*point);
        for (int i=0;i<N;i++){
            point_temp[i]=(*point)[i]-a*(*grad)[i];
        }
        if (func_at_point(&point_temp)<f){
            a*=m;
            continue;
        }
        for (;;d=1){
            a*=l;

            point_temp=(*point);
            for (int i=0;i<N;i++){
                point_temp[i]-=a*(*grad)[i];
            }
            f_prev=func_at_point(&point_temp);
            if (f_prev<f){
                (*point)=point_temp;
                b=1;
                break;
            }
            for (int i=0;i<N;i++){
                if ((f_prev-f!=0)||(a*(*grad)[i]>pow(E/4,0.5))) { d=0; }
            }
            if (d==1){
                (*point)=point_temp;
                b=1;
                break;
            }
        }
    }
    return func_at_point(point);
}
double Find_Min_2(vector <double> *point, vector <double> *grad, vector < vector <double> > *grad2){

    grad2_at_point(point, grad2);
    grad_at_point(point, grad);
    double f=func_at_point(point);
    vector < vector <double> > grad2_rev=(*grad2);
    vector < double > point_temp1=(*point), point_temp2=(*point);
    rev(grad2, &grad2_rev);
    double d=0,e=1,a=-1,b=0,c=0;
    for (int i=0;i<N;i++){ for (int j=0;j<N;j++){ cout<<(*grad2)[i][j]<<" ";}cout<<endl;}
    vector < double > move_vec (N);
    for (int i=0;i<N;i++){
        move_vec[i]=0;
        for (int j=0;j<N;j++){
            move_vec[i]-=grad2_rev[i][j]*(*grad)[j];
        }
    }
    grad_at_point(point,grad);
    e=0; for (int i=0;i<N;i++) {e+=(*grad)[i]*move_vec[i];}
    if (e<0) {
        a=+1;
    }
    for (;;a*=2){cout<<c<<"  "<<e<<endl;
        for (int i=0;i<N;i++){
                point_temp1[i]=(*point)[i]+a*(move_vec)[i];
        }
        grad_at_point(&point_temp1,grad);
        c=0; for (int i=0;i<N;i++) {c+=(*grad)[i]*move_vec[i];}
        if (c*e<0){
            break;
        }
    }
    d=(func_at_point(&point_temp1)-func_at_point(point)+e*0-c*a)/(e-c);
    for (int i=0;i<N;i++){
            (*point)[i]=(*point)[i]+d*(move_vec)[i];
    }
    return func_at_point(point);
}
int main()
{
    vector <double> point_start (N);
    for (int i=0;i<N;i++) { point_start[i]=-3; }
    vector <double> point=point_start;
    vector <double> grad (N);
    vector <double> point_temp=point;
    vector < vector <double> > grad2 (N);
    for (int i=0;i<N;i++) { for (int j=0;j<N;j++) {grad2[i].push_back(0);} }
    double a=func_at_point(&point),n=0;
    double b=a,c,d;
    for(;;){
        point_temp=point;
        b=Find_Min(&point, &grad);
        cout<<b<<" at step "<<++n<<" at point ("; for (int i=0;i<N;i++) { cout<<point[i]<<", "; } cout<<")"<<endl;
        c=0; for (int i=0;i<N;i++) {c+=pow(point[i]-point_temp[i],2);} c=pow(c,0.5);
        grad_at_point(&point,&grad);
        d=0; for (int i=0;i<N;i++) {d+=pow(grad[i],2);} d=pow(d,0.5);
        if ((b-a<pow(E,0.5))&&(c<pow(E,0.5))&&(d<pow(E,0.5))) { break; }
        a=b;
    }
    cout<<"End of 1 path"<<endl;
    n=0;
    for(;;){
        point_temp=point;
        b=Find_Min_2(&point, &grad, &grad2);
        cout<<b<<" at step "<<++n<<" at point ("; for (int i=0;i<N;i++) { cout<<point[i]<<", "; } cout<<")"<<endl;
        c=0; for (int i=0;i<N;i++) {c+=pow(point[i]-point_temp[i],2);} c=pow(c,0.5);
        grad_at_point(&point,&grad);
        d=0; for (int i=0;i<N;i++) {d+=pow(grad[i],2);} d=pow(d,0.5);
        if ((b-a<E)&&(c<E)&&(d<E)) { break; }
        a=b;
    }
    cout<<endl;
    cout<<"f = "<<func_at_point(&point)<<" at point (";
    for (int i=0;i<N;i++) { cout<<point[i]<<", "; }
    cout<<")"<<endl;
    return 0;
}
