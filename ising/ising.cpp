#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include<cstdlib>
using namespace std;



/* DE = E - E'  E' = E - DE*/

/* this is short program to get the ground state configuration of the ising model in 1D*/

double eV = 1.602176634e-19;

double  J = 1*eV; // coupling constant
unsigned int N = 1000; // number of spins 
double kb = 1.380658e-23; 

double Energy(int *spins){
  double energy = 0;
  for(int i = 0; i < N-1; i++){
  energy -= spins[i] * spins[i+1];
  }
  energy *= J;
  return energy;
}

double EnergyDifference(int spins[], unsigned int index){
  double DE = 0;
  if(index != N-1 and index != 0){
    DE += -2*J*(spins[index-1] * spins[index] + spins[index] * spins[index+1]);
  }
  else if(index == N-1){
    DE += -2*J*(spins[index-1] * spins[index]);
  }
  else if(index == 0){
    DE += -2*J*(spins[index] * spins[index+1]);
  }
  return DE;
}

void initializeSpins(int spins[]){

	// Loop to assign sing of spin randomly
	for(int i=0; i<N; i++){
		
    srand((unsigned) time(NULL));

	  // Retrieve a random number
	  int random = rand();

    spins[i] = pow(-1,random);
  }
}


double getMagnetization(int spins[]){
  double m = 0;

  // magnetization for ferromagnetic coupling
  for(int i = 0; i<N; i++){
      m += spins[i];
    }

  m /= N;
  return m;
}


void getGroundstateconfig(double T, int spins[]){

  // initialize spins
  initializeSpins(spins);

  // get the Energy of initial configuration
  double E = Energy(spins), DE; 

  // Number of times a spins is randomly flipped
  unsigned int iterations = 10000;
  
  // Providing a seed value
	srand((unsigned) time(NULL));

	// Loop to get the index of spin which is flipped
	for(int i=1; i<=iterations; i++){
		
	  // Retrieve a random number between 0 and number of spins
	  int random = (rand() % N);

    //generate random threshold for thermal flucutations

    double Z = (double)rand() / ((double)RAND_MAX + 1);
    
    // calculate Energy difference between current and possible new configuration
    DE = EnergyDifference(spins, random);

  // check if new configuration is favourable
    if(exp((DE)/(kb * T)) > Z){
      spins[random] *= -1;
      E -= DE;
    }
  }

  // print the ground state config
  for(int i = 0; i<N;i++){
    cout << spins[i] << " ";
  }
  cout << endl << endl;
}

int main(){

  unsigned int numTemp = 100;
  double temperatures[numTemp], magnetization[numTemp];

  int spins[N];

  /* loop through different temperatures*/

  for(int i = 0; i < numTemp; i++){
    temperatures[i] = exp(0.35*i-10); // make array with temperatures
    getGroundstateconfig(temperatures[i], spins); //get the ground state config at fixed temperature
    magnetization[i] = getMagnetization(spins); //calculate the magnetization at fied temperature
  }

  ofstream myfile ("magnetization.txt");

  if (myfile.is_open())
  {
    for(int i = 0; i < numTemp; i++){
      myfile << temperatures[i] << "\t" << magnetization[i] << endl;
      }
      myfile.close();
  }
  else cout << "Unable to open file";

return 0;
}
