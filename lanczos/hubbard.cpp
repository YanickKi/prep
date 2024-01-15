#include <iostream>
#include </home/yanick/eigen/Eigen/Dense>
#include </home/yanick/eigen/unsupported/Eigen/CXX11/Tensor>    
#include <cmath>
#include <vector>
using  namespace std;
using namespace Eigen;

unsigned int N = 2; //number of sites
int numStates = 6;


int overlap(VectorXd bra, VectorXd ket){

    int overlap = 0;
    for(int i = 2; i < 2*N+2; i++){
        overlap += bra(i) * ket(i);
    }
    
    if(overlap == 2){
        return 1; 
    }
    else{
        return 0;
    }

}

int numoccbefore(VectorXd state, int index){
    int sum = 0;
    for(int i = 2; i < index; i++){
        sum += state(i);
    }
    return sum;
}

tuple <VectorXd, double> hopping(VectorXd state, int annhilate, int create){
    int sign = 1;
    int 
    if(state(annhilate) == 0 || state(create) == 1){
        return {VectorXd::Zero(2*N+4),0}; 
    }

    sign *= pow((-1),numoccbefore(state, annhilate));
    state(annhilate) = 0;

    if(annhilate < 2){
        state(annhilate+2*N) = state(annhilate);
    }
    if(annhilate > 2*N+1){
        state(annhilate-2*N) = state(annhilate);
    }


    sign *= pow((-1),numoccbefore(state, create));
    state(create) = 1;

    if(create < 2){
        state(create+2*N) = state(create);
    }
    if(annhilate > 2*N+1){
        state(create-2*N) = state(create);
    }
    //state(2) = state(2*N+2);
    //state(3) = state(2*N+3);
    //state(2*N+1) = state(0);
    //state(2*N+2) = state(1);

    return {state, sign};
}

double Matrixelement(VectorXd initial, VectorXd final){
    double U=10, t = 1; 
    double M = 0; 
    VectorXd result = VectorXd::Zero(2*N+4);
    int sign;
    // this function acts the Hamiltonian on a state and returns the outcoming state and the value 
    // H |i>  = E_ji |j> 

    for(int i = 2; i < 2*N+2; i += 2){
        M += U*initial(i)*initial(i+1)*overlap(final, initial);
    }
    for(int i = 2; i < 2*N+2; i += 2){
        
        tie(result, sign) = hopping(initial,i,i+2);
        M += -t*sign*overlap(final, result);

        tie(result, sign) = hopping(final,i-1,i+1);
        M += -t*sign*overlap(final, result);

        tie(result, sign) = hopping(initial,i,i-2);
        M += -t*sign*overlap(final, result);
        
        tie(result, sign) = hopping(initial,i+1,i-1);
        M += -t*sign*overlap(final, result);
    }

    return M;
}

MatrixXd init(){

    MatrixXd basis(2*N+4, numStates); // 6 states, N sites, 2 for spin up and down per site + 4 for PBC
    int count = 0;
    for(int i = 2; i < 2*N+1; i++){
        for(int j = i+1; j < 2*N+2; j++){
            basis(i, count) = 1;
            basis(j, count) = 1;
            count += 1;
        }
    }
    basis.block(0,0,2,numStates) = basis.block(2*N+1,0,2,numStates);
    basis.block(2*N+2,0,2,numStates) = basis.block(2,0,2,numStates);

    return basis;
}

int main(){

    MatrixXd basis = init();
    MatrixXd H(numStates, numStates);
    cout << "Basis vectors: " << endl << basis << endl;
    for(int i = 0; i<numStates; i++){
        for(int j = 0; j < numStates; j++){
            H(i,j) = Matrixelement(basis.col(i), basis.col(j));
        }
    }
    cout << "H: " << endl << H << endl;

    VectorXd test(2*N+4);

    test << 0,1,1,0,0,1,1,0; 
    VectorXd result(2*N+4);
    int sign;
    cout << "test2: " << test(2) << " test5: " << test(4) << endl;
    tie(result, sign) = hopping(test, 2, 4);
    cout << "TEST RESULT: " << sign*result << endl;
    return 0;
}