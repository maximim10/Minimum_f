#include <iostream>
#include <vector>
#include <cmath>
const int N=2;
const double E=0.001;
using namespace std;
double func_at_point(vector <double> *point){
    double c=pow((*point)[0],2)+pow((*point)[1],6);
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
            //cout<<i+k<<" "<<perest[i]+l<<endl;
            d_temp*=(*matr)[i+k][perest[i]+l];
        }//cout<<endl<<endl;
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
    //cout<<111;
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){

            (*matr_rev)[i][j]=Alg_dop(i,j,matr)/d;
        }
    }
}
double Find_Min (vector <double> *point, vector <double> *grad){
    int f=func_at_point(point);
    double Min=func_at_point(point), Min_temp=Min, Min_comp=Min;
    vector <double> point_temp1=(*point), point_temp2=(*point), point_temp, grad_temp (N);
    grad_at_point(point, grad);
    double a=1,b=1,m=2,l=0.75;
    for (int b=0;b<1;){
        point_temp=(*point);
        for (int i=0;i<N;i++){
            point_temp[i]=(*point)[i]-a*(*grad)[i];
        }
        if (func_at_point(&point_temp)<f){
            a*=m;
            continue;
        }
        for (;;){cout<<a;
            a*=l;
            point_temp=(*point);
            for (int i=0;i<N;i++){
                point_temp[i]=(*point)[i]-a*(*grad)[i];
            }
            if (func_at_point(&point_temp)<f){
                (*point)=point_temp;
                b=1;
                break;
            }
            if (a<4e-306){
                b=1;
                break;
            }
        }
    }

}
double Find_Min_2(vector <double> *point, vector <double> *grad, vector < vector <double> > *grad2){
    grad2_at_point(point, grad2);
    grad_at_point(point, grad);
    int f=func_at_point(point);
    vector < vector <double> > grad2_rev=*grad2;
    vector < double > point_temp=(*point);
    double a=1,b=1,m=2,l=0.75;
    rev(grad2, &grad2_rev);
    vector < double > move_vec (N);
    //cout<<(*point)[1]<<endl;
    //cout<<(*grad2)[0][0]<<" "<<(*grad2)[0][1]<<endl; cout<<(*grad2)[1][0]<<" "<<(*grad2)[1][1]<<endl;
    //cout<<(grad2_rev)[0][0]<<" "<<(grad2_rev)[0][1]<<endl; cout<<(grad2_rev)[1][0]<<" "<<(grad2_rev)[1][1]<<endl;
    for (int i=0;i<N;i++){
        move_vec[i]=0;
        for (int j=0;j<N;j++){
            move_vec[i]-=grad2_rev[i][j]*(*grad)[j];
            //cout<<move_vec[i]<<endl;
        }
    }
    //cout<<move_vec[0]<<" "<<move_vec[1]<<endl;
    for (int b=0;b<1;){
        point_temp=(*point);
        for (int i=0;i<N;i++){
            point_temp[i]=(*point)[i]+a*move_vec[i];
        }
        if (func_at_point(&point_temp)<f){
            a*=m;
            continue;
        }
        for (;;){
            a*=l;
            point_temp=(*point);
            for (int i=0;i<N;i++){
                point_temp[i]=(*point)[i]+a*move_vec[i];
            }
            if (func_at_point(&point_temp)<f){
                (*point)=point_temp;
                b=1;
                break;
            }
            if (a<4e-306){
                b=1;
                break;
            }
        }
    }
    return func_at_point(point);

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
    for (int i=0;i<N;i++) { point_start[i]=3; }
    vector <double> point=point_start;
    vector <double> grad (N);
    vector <double> point_temp=point;
    vector < vector <double> > grad2 (N);
    //vector < vector <double> > matr (N);
    //vector < vector <double> > matr_rev (N);
    //for (int i=0;i<N;i++) { for (int j=0;j<N;j++) {matr[i].push_back(0);} }
    //for (int i=0;i<N;i++) { for (int j=0;j<N;j++) {matr_rev[i].push_back(0);} }

    //matr[0][0]=2; matr[1][1]=2; matr[0][1]=2; matr[1][0]=0;
    //rev(&matr, &matr_rev);
    //for (int i=0;i<N;i++){
        //for (int j=0;j<N;j++){
    //        cout<<matr_rev[i][j]<<" ";
    //    }
    //    cout<<endl;
    //}
    for (int i=0;i<N;i++) { for (int j=0;j<N;j++) {grad2[i].push_back(0);} }
    double a=func_at_point(&point),n=0;
    double b=a,c,d;

    cout<<a<<endl;
    for(;;){
        point_temp=point;
        b=Find_Min(&point, &grad);
        cout<<b<<" at step "<<++n<<" at point ("; for (int i=0;i<N;i++) { cout<<point[i]<<", "; } cout<<")"<<endl;
        c=0; for (int i=0;i<N;i++) {c+=pow(point[i]-point_temp[i],2);} c=pow(c,0.5);
        grad_at_point(&point,&grad);
        cout<<grad[0]<<" "<<grad[1]<<endl;
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
