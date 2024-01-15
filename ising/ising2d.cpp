#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include<cstdlib>
using namespace std;



/* DE = E' - E  <=> E' = E + DE*/

/* this is short program to get the ground state configuration of the ising model in 2D*/

double eV = 1.602176634e-19, kb = 1.380658e-23;

double  J = 1*eV, h = 0; // coupling constant
unsigned const int N = 100; // number of spins and temperatures 
double Tc = 2*abs(J)/(kb * log(1+sqrt(2)));
int numTemp = 200;
unsigned int eqSteps = 500, mcSteps = 500;
	
void writeToFile(double array1[], double array2[], unsigned int size ,string filename){
  
  ofstream file (filename);

  if (file.is_open())
  {
      for(int i = 0; i < size; i++){
        file << array1[i] << " " << array2[i] << endl;
      }
      file.close();
  }
  else cout << "Unable to open file";
}

double Energy(int spins[][N])
{
  double energy = 0;
  // assuming periodic boundary conditions
  for(int i = 0; i < N;i++){
    for(int j = 0; j < N; j++){
      energy -= J*(spins[(j+1)%N][i]*spins[j][i]+spins[j][(i+1)%N]+spins[j][i]);
      energy -= h*(spins[j][i]);
    }

  }
  return energy;
}

double getMagnetization(int spins[][N])
{
  double m = 0;

  for (int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++){
    m += spins[i][j];
  }
  }
  return m;
}

double EnergyDifference(int spins[][N], unsigned int x, unsigned int y){
  double DE = 0;
  unsigned int NNx[2], NNy[2]; //nearest neighbours in x and y direction

  if(x == 0){
    NNx[0] = N-1; 
    NNx[1] = 1;
  }
  else if(x == N-1){
    NNx[0] = N-2;
    NNx[1] = 0;
  }
  else{
    NNx[0] = x-1;
    NNx[1] = x+1;
  }

  if(y == 0){
    NNy[0] = N-1;
    NNy[1] = 1;
  }
  else if(y == N-1){
    NNy[0] = N-2;
    NNy[1] = 0;
  }
  else{
    NNy[0] = y-1;
    NNy[1] = y+1;
  }

  for(int i = 0; i < 2; i++){
  DE += 2*J*(spins[x][y] * spins[NNx[i]][y]);
  DE += 2*J*(spins[x][y] * spins[x][NNy[i]]);
  }
  DE += 2*h*spins[x][y];

  return DE;
}

void initializeSpins(int spins[][N])
{

  // Loop to assign sign of spin randomly
  for (int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++){
    
    srand((unsigned)time(NULL));

    // Retrieve a random number
    //int random = rand();

    spins[i][j] = 1;
   
    }
  }
}

void mcmove(double T, int spins[][N]){

  // get the enegie of the inital config at a given Temperature T
  double DE; 

  // Number of times a spin  is randomly flipped
  unsigned int iterations = N*N;
  
  // Providing a seed value
	srand((unsigned) time(NULL));

	// Loop to get the index of spin which is flipped
	for(int i=1; i<=iterations; i++){
		
	  // Retrieve a random number between 0 and number of spins
	  unsigned int randomx = (rand() % N);
    unsigned int randomy = (rand() % N);

    //generate random threshold for thermal flucutations
    
    double Z = (double)rand() / ((double)RAND_MAX + 1);

    // calculate Energy difference between current and possible new configuration
    DE = EnergyDifference(spins, randomx, randomy);

  // check if new configuration is favourable
    if(DE <= 0){
      spins[randomx][randomy] *= -1;
    }

    else if(exp((-DE)/(kb * T)) > Z){
  
      spins[randomx][randomy] *= -1;

      }  
  }


}

int main(){

  double temperatures[numTemp], magnetization[numTemp], energies[numTemp], specificHeat[numTemp], susceptibility[numTemp];
  double Ene, Mag, n1 = 1.0/(mcSteps * N * N), n2 = 1.0/(mcSteps * mcSteps * N * N)  ;
  int spins[N][N];
  initializeSpins(spins);
  int tempSpins[N][N];

  // assign temperature array with temperatures

  for(int i = 0; i < numTemp;i++){
    temperatures[i] = 5 * Tc* i/(numTemp -1);
  }

  ofstream gsfile ("groundstates2d.txt");

  if (gsfile.is_open()){

  for(int i = 0; i < numTemp; i++){
    double E1 =0, E2 = 0, M1 = 0, M2 = 0;

    for(int j = 0; j < eqSteps; j++){
      mcmove(temperatures[i], spins); //equilibrate
    }

    copy(&spins[0][0], &spins[0][0]+N*N,&tempSpins[0][0]);    
    // write configuration in equilibrium in file
      for(int k = 0; k<N;k++){
        for(int l=0; l<N; l++){
          gsfile << spins[k][l] << " ";
        }
        gsfile << endl; 
      }
      gsfile << endl; 
    
  // now actually measure in equilibrium
    for(int j = 0; j < mcSteps; j++){
      mcmove(temperatures[i], tempSpins);
      Ene = Energy(tempSpins);
      Mag = getMagnetization(tempSpins);

      E1 += Ene;
      E2 += Ene*Ene;
      M1 += Mag;
      M2 += Mag*Mag;
    }
    
  energies[i] = E1*n1;
  magnetization[i] = M1*n1;
  specificHeat[i] = (n1*E2 - n2*E1*E1) / (kb * temperatures[i] * temperatures[i]);
  susceptibility[i] = 1/(kb*temperatures[i]) * (n1*M2-n2*M1*M1);
  }

    gsfile.close();
  }
  else cout << "Unable to open file";

  writeToFile(temperatures, magnetization,numTemp, "magnetization2d.txt");
  writeToFile(temperatures, energies,     numTemp,"energies2d.txt");
  writeToFile(temperatures, specificHeat,     numTemp,"specificHeat2d.txt");
  writeToFile(temperatures, susceptibility,     numTemp,"susceptibility2d.txt");

  ofstream InputData ("InputData.txt");

  if (InputData.is_open())
  {
      InputData << N  << endl << numTemp << endl << Tc;
      InputData.close();
  }
  else cout << "Unable to open file";

return 0;


}