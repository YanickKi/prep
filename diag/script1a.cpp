#include <iostream>
#include </home/yanick/.include/Eigen/Dense>    
#include<cmath>
#include <vector>
using  namespace std;


void CompEigenvalueAndEigenvector(Eigen::Matrix4d M, Eigen::Vector4d & Eigenvector, double & Eigenvalue){
    double error = 1;
    double TempEigenvalue = Eigenvector.transpose() * M * Eigenvector; //initial eigenvalue

    while(error > 1e-5){
        Eigenvector = M*Eigenvector;
        Eigenvector.normalize();
        Eigenvalue = Eigenvector.transpose() * M * Eigenvector;
        error = abs(Eigenvalue - TempEigenvalue);
        TempEigenvalue = Eigenvalue;
    }
}

Eigen::Matrix4d NewMatrix(Eigen::Matrix4d M, Eigen::Vector4d Eigenvector, double Eigenvalue){

    return M - Eigenvalue * Eigenvector * Eigenvector.transpose();
}

int main(){

Eigen::Matrix4d A;

A << 1, 2, 3, 4,
     2, 5, 2, 1,
     3, 2, 5, 3,
     4, 1, 3, 4;

Eigen::VectorXcd eivals = A.eigenvalues();
cout << "EVs in a): " << endl << eivals << endl;

// part b

double Eigenvalue; 
Eigen::Vector4d Eigenvector (1,1,1,1);

cout << "Eivs in b): ";

for(int i = 0; i < 4; i++){
    CompEigenvalueAndEigenvector(A, Eigenvector, Eigenvalue);
    A = NewMatrix(A, Eigenvector, Eigenvalue);
    cout << Eigenvalue << " ";
}



//ex. 2

cout << endl << "Exercise 2" << endl;

unsigned int L = 10; // distance Left and right from origin \xi \in [-L,L]
double Dxi = 0.1, lambda = 0;
int N = 2*L/Dxi;
cout << "N ist: " << N << endl;

Eigen::MatrixXd deltam1 = Eigen::ArrayXXd::Zero(N+1,N+1); //\delta_{n,m-1}
Eigen::MatrixXd delta   = Eigen::ArrayXXd::Zero(N+1,N+1); //\delta_{n,m}
Eigen::MatrixXd deltap1(N+1,N+1);
Eigen::MatrixXd nMatrix = Eigen::ArrayXXd::Zero(N+1,N+1);
Eigen::MatrixXd H(N+1,N+1);
Eigen::VectorXd eigenvalues(N+1);

//diagonal matrices

for(int n = 0; n < N+1; n++){
    delta(n,n) = 1;
    nMatrix(n,n) = L/Dxi*(2.*n/N-1); // integer indices but n \in [-L/d\xi, ..., L/d/xi]
}

//off diagonal matrices

for(int n = 0; n < N; n++){
    deltam1(n,n+1) = 1;
}

deltap1 = deltam1.transpose();

//now construct the actual Hamilton matrix

H = -1/pow(Dxi,2) * (deltam1 + deltap1 - 2*delta)+ (pow(Dxi,2)*nMatrix*nMatrix + lambda * pow(Dxi,4)*nMatrix*nMatrix*nMatrix*nMatrix);

eigenvalues = H.eigenvalues().real();

sort(begin(eigenvalues), end(eigenvalues));

cout << nMatrix * nMatrix << endl;

cout << "10 Smallest eigenvalues are: ";

for(int i = 0; i < 10; i++){
    cout << eigenvalues(i) << " "; 
}

cout << endl;

return 0;
}