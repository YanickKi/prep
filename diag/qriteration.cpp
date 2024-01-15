#include <iostream>
#include </home/yanick/.include/Eigen/Dense>    
#include <cmath>
#include <vector>
using  namespace std;


tuple<Eigen::MatrixXd, Eigen::MatrixXd> QRdecomp(Eigen::MatrixXd A, unsigned int N){

    // This is a function to determine the QR decomposition Q,R of a tridiagonal matrix A 
    double r,c,s;
    Eigen::MatrixXd G = Eigen::MatrixXd::Identity(N,N), Q = Eigen::MatrixXd::Identity(N,N);

    for (int i = 0; i<N-1; i++){
        r = sqrt(pow(A(i+1,i),2)+pow(A(i,i),2));
        c = A(i,i)/r;
        s = -A(i+1,i)/r;    
        G.block(i,i,2,2) << c,-s,
                            s,c;
        
        Q = Q*G.transpose();
        A = G*A;
        cout << "Q " << i << ":" << endl << Q << endl; 
        cout << "G: " << i << ": " << endl << G << endl;
        cout << "QQT: " << i << endl << Q.transpose()*Q << endl;    
        G.block(i,i,2,2) = Eigen::MatrixXd::Identity(2,2);
    }
    
    //Eigen::MatrixXd R = Q.transpose()*A;
    cout << "A:" << endl << A << endl;
    cout << "Q:" << endl << Q << endl;
    //cout << "R:" << endl << R << endl; 

    return {Q,A};
}


int sign(double x){
    if(x<0){
        return -1;
    }
    else{
        return +1;
    }
}


int main(){
    double l = 1, k = 5;
    unsigned int N = 3; //dimension of matrix
    Eigen::MatrixXd A (N,N);
    //Eigen::MatrixXd A = Eigen::MatrixXd::Constant(N,N,k+l);


    //for(int i = 0; i < N; i++){
    //    A(i,i) += k; 
    //}


    A <<    6,5,0,
            5,1,4,
            0,4,3;

    double kn;

    cout << "Householder transformation" << endl;

    //BEGIN Part a) Householder transformation

    for(int i = 0; i < N-2; i++){


        //big loop for going from dimension 1 to N-2

        Eigen::VectorXd v(N-i-1), u(N-i-1), unitVec = Eigen::ArrayXd::Zero(N-i-1);
        Eigen::MatrixXd P = Eigen::MatrixXd::Zero(N,N); 

        unitVec(0) = 1;

        // now get the value kn for the current dimension and the vector u
        v = A.col(i).tail(N-i-1);
        kn = - sign(A(i+1,i)) * v.norm();
        u = (v - kn * unitVec);
        u.normalize();

        //top left block from the transformation matrix

        P.block(0,0,i+1, i+1) = Eigen::MatrixXd::Identity(i+1,i+1);

        //bottom right of the transformation matrix
        P.block(i+1,i+1, N-i-1, N-i-1) = Eigen::MatrixXd::Identity(N-i-1,N-i-1) - 2 * u*u.transpose();
        //cout << "P" << i << endl << P << endl;
        A = P.transpose() * A * P;
    }

    cout << A << endl << endl;
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(abs(A(i,j)) < 1e-10){
                A(i,j) = 0; 
            }
        }
    }
    //cout << "Nach 0 setzen: " << endl << A << endl;
    //END Part a) Householder transformation

    //BEGIN Part b) QR iteration

    cout << "QR iteration" << endl;

    A <<    6,5,0,
            5,1,4,
            0,4,3;

    for(int i = 0; i<10000; i++){
        Eigen::MatrixXd Q(N,N),R(N,N); 
        tie(Q, R) = QRdecomp(A,N);
        A = R*Q;
    }


    cout << "Diagonal Matrix: " << endl << A;
    //END Part b) QR iteration

    return 0;
}