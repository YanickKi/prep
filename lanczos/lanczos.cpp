#include <iostream>
#include </home/yanick/.include/Eigen/Dense>    
#include <cmath>
#include <vector>
using  namespace std;
using namespace Eigen;


Eigen::VectorXd linspace(double start, double end, int steps){
    Eigen::VectorXd vec(steps);
    double delta = (end-start)/(steps-1);
    for(int i = 0; i < steps; i++){
        vec(i) = start+ i*delta;
    }
    return vec;
}


MatrixXd lanczos(MatrixXd A, double dim, double iter){    

     //iter will be dimension of subspace and dim the dimension of the full problem
    
    MatrixXd A (dim,dim);

    cout << "A: " << endl << A << endl << endl;

    MatrixXd Q(dim,iter); //Matrix with Lanczos vectors as rows
    Q.col(0) = VectorXd::Random(dim);
    Q.col(0).normalize();
    VectorXd v(dim);
    double alpha, beta;

    //do the actual lanczos iteration

    for(int i = 0; i < iter-1; i++){
        v = A*Q.col(i);
        alpha = Q.col(i).transpose()*v;
        if(i == 0){
            v = v - alpha *Q.col(i);
        }
        else{
            v = v - beta*Q.col(i-1) - alpha *Q.col(i);
        }
        beta = v.norm();
        Q.col(i+1) = v/beta;
        
    }
    cout << "Q: " << endl << Q << endl << endl;

    MatrixXd T = Q.transpose() * A * Q;

    cout << "T: " << endl << T << endl << endl;
    
    cout << "Eigenvalues of A: " << endl << A.eigenvalues().real() << endl << endl;

    cout << "First " << iter << " approximate eigenvalues of A: " << endl << T.eigenvalues().real() << endl << endl;

    cout << "Orthogonality check:" << endl;

    for(int i = 0; i < Q.cols()-1; i++){
        cout << i << ", " << i+1 << ": " << Q.col(i).transpose()*Q.col(i+1) << endl;
    }
}
