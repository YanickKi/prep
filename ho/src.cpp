#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;


/* df/dx = (f(x+h)-f(x-h))/2h = f'
    f'' = (f'(x+h)-f'(x-h))/2h 
        = (f(x+2h)-f(x)-f(x) + f(x-2h))4h**2
        = (f(x+2h)-2f(x)+f(x-2h))/4h**2

Initial condition :
at i = 0 there is f[0] = 1;
the function does not change there
=> f[1] = 1 

DEQ of harm osc: d^2f/dt^2 + w^2 f(x) = 0

(f(x+2h)-2f(x)+f(x-2h))/4h**2 + w^2 f(x) = 0


update equation: f[i+1] = 2f[i] - f(i-1) - 4h**2 w^2 f[i]

solve with f(x=0)=0 => starts with cosine 

*/

int main(){

  double h = 0.001, w = 5, STEPS = 1000;

  vector<double> f{1, 1};

  for(int i = 1; i < STEPS  ; i++){
    f.push_back(2* f[i] - f[i-1] - 4*h*h*w*w *f[i]);
  }


  ofstream myfile ("data.txt");

  if (myfile.is_open())
  {
    for(int i = 0; i < f.size(); i++){
      myfile << h * i;
      myfile << "\t" << f[i] << endl;//HIER CHECKEN OB ES DIE RICHTIGEN ELEMENTE SIND
      }
      myfile.close();
  }
  else cout << "Unable to open file";

return 0;
}
