#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using  namespace std;


double func(double x){
    return pow(x,2)-2;
}


double derivative(double ( *funptr )(double), double x, double h){
    return (funptr(x+h) - funptr(x-h))/(2*h);
}

double secondDerivative(double ( *funptr )(double), double x, double h){
    return (funptr(x+2*h) + funptr(x-2*h) - 2*funptr(x))/(4*pow(h,2));
}

void NewIntervall(double u, double &inner, double &outer, double &rest, double (*funptr) (double)){
    if(funptr(u) < funptr(inner)){
        rest = inner; 
        inner = u;
        outer = outer;
    }
    else{
        rest = rest; 
        inner = inner;
        outer = u;
    }
}

void intervallHalbierung(double ( *funptr )(double), double x, double y, double z, double error){
    ofstream myfile;
    myfile.open ("intervallHalb.txt");
    myfile << x << "\t" << y << "\t" << z << endl; 
    
    while(z-x >= error){
        if(y-x < z-y){
            NewIntervall((y+z)/2, y, z, x, *funptr);
        }
        else {
            NewIntervall((x+y)/2, y, x, z, *funptr);
        }
    myfile << x << "\t" << y << "\t" << z << endl; 
    }
    myfile.close();
    
}

void newton(double ( *funptr ) (double), double x, double error){
    double xtemp = x + error + 1;
    double h = 1e-6;

    ofstream myfile;
    myfile.open ("newton.txt");
    myfile << x;

    while(abs(x-xtemp) > error){
        xtemp = x;
        x = x + - (derivative(funptr, x , h))/(secondDerivative(funptr, x , h));
        myfile << x << endl; 
    }
    
    myfile.close();
}

int main(){
    double min;
    double x0 = -0.5, y0 = -0.1, z0 = 2, error = 1e-9;

    intervallHalbierung(func, x0, y0, z0, error);
    x0 = 1; // newton only needs one inital value
    newton(func, x0, error);
    //cout << "Minimum is at" << min << " with value f(min) = " << func(min) << endl;
    return 0;
}