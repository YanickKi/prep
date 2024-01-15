#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include<cstdlib>
using namespace std;



/* DE = E' - E  <=> E' = E + DE*/

/* this is short program to get the ground state configuration of the ising model in 1D*/

double eV = 1.602176634e-19;

double  J = 1e-3*eV, h = J/100; // coupling constant
unsigned const int N = 100, numTemp = 50; // number of spins and temperatures 
double kb = 1.380658e-23;
double Tc = 2*J/(kb * log(1+sqrt(2)));

	
double Energy(int spins[][N])
{
  double energy = 0;

  //iterate through the c enter
  for (int i = 0; i < N-1; i++){
    for(int j = 0; j < N-1; j++){
      energy -= J*(spins[i][j] * spins[i+1][j]);
      energy -= J*(spins[i][j] * spins[i][j+1]);
      energy -= h* spins[i][j];
    }
    // iterate through the borders
    energy -= J*(spins[i][N] * spins[i+1][N]);
    energy -= J*(spins[N][i] * spins[N][i+1]);
    energy -= h*(spins[i][N] + spins[N][i]);
  }
  energy -= h*spins[N][N];

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
  m /= N*N;
  return m;
}

double EnergyDifference(int spins[][N], unsigned int x, unsigned int y){
  double DE = 0;
  if(x==0){
    if(y == 0){
      DE += -2*J*(spins[x+1][y]*spins[x][y] + spins[x][y]*spins[x][y+1]); // above and right
    }
    else if(y == N-1){
      DE += -2*J*(spins[x+1][y]*spins[x][y] + spins[x][y-1]*spins[x][y]); // above and left
    }
    else{
      DE += -2*J*(spins[x+1][y]*spins[x][y] + spins[x][y-1]*spins[x][y] + spins[x][y]*spins[x][y+1]); // above left and right 
    }
  }
  else if(x==N-1){
    if(y == 0){
      DE += -2*J*(spins[x-1][y]*spins[x][y] + spins[x][y]*spins[x][y+1]); // under and right
    }
    else if(y == N-1){
      DE += -2*J*(spins[x-1][y]*spins[x][y] + spins[x][y-1]*spins[x][y]); // under and left
    }
    else{
      DE += -2*J*(spins[x-1][y]*spins[x][y] + spins[x][y-1]*spins[x][y] + spins[x][y]*spins[x][y+1]); // under left and right 
    }
  }
  else if(y==0){
    DE += -2*J*(spins[x-1][y]*spins[x][y] + spins[x+1][y]*spins[x][y] + spins[x][y+1]*spins[x][y]); // left right and up
  }
  else if(y==N-1){
    DE += -2*J*(spins[x-1][y]*spins[x][y] + spins[x+1][y]*spins[x][y] + spins[x][y-1]*spins[x][y]); // left right and down
  }

  else {
    DE += -2*J*(spins[x-1][y] * spins[x][y] + spins[x][y] * spins[x+1][y]); // interaction left and right
    DE += -2*J*(spins[x][y-1] * spins[x][y] + spins[x][y] * spins[x][y+1]); // interaction above and under
  }
  DE += -2*h*spins[x][y];
  return -DE;
}

void initializeSpins(int spins[][N])
{

  // Loop to assign sign of spin randomly
  for (int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++){
    
    srand((unsigned)time(NULL));

    // Retrieve a random number
    int random = rand();

    spins[i][j] = pow(-1, random);
   
    }
  }
}


void getGroundstateconfig(double T, int spins[][N], ofstream& filename){

  // initialize spins
  initializeSpins(spins);

  // get the Energy of initial configuration
  double E = Energy(spins), DE; 

  // Number of times a spins is randomly flipped
  unsigned int iterations = 500;
  
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
      E += DE;
    }

    else {
      if(exp((-DE)/(kb * T)) > Z){
  
      spins[randomx][randomy] *= -1;
  
      E += DE;

      }
    }
  
  }

  for(int i = 0; i<N;i++){
    for(int j=0; j<N; j++){
      filename << spins[i][j] << " ";
    }
    filename << endl; 
  }
    filename << endl; 
}

int main(){

  double temperatures[numTemp], magnetization[numTemp];

  int spins[N][N];

  /* loop through different temperatures*/

  ofstream gsfile ("groundstates2d.txt");

  if (gsfile.is_open()){
    
  // make array with temperatures

  temperatures[0]=1e-6;
  
  for(int i = 1; i < numTemp;i++){
    temperatures[i] = 1000/kb * i/(numTemp-1);
  }

  for(int i = 0; i < numTemp; i++){
    getGroundstateconfig(temperatures[i], spins, gsfile); //get the ground state config at fixed temperature
    magnetization[i] = getMagnetization(spins); //calculate the magnetization at fixed temperature
  }

    gsfile.close();
  }
  else cout << "Unable to open file";

  ofstream FileMagn ("magnetization2d.txt");

  if (FileMagn.is_open())
  {
    for(int i = numTemp-1; i >= 0; i--){
      FileMagn << temperatures[i] << "\t" << magnetization[i] << endl;
      }
      FileMagn.close();
  }
  else cout << "Unable to open file";

  ofstream InputData ("InputData.txt");

  if (InputData.is_open())
  {
      InputData << N << "\t" << "Number of spins" << endl
                << numTemp << "\t" << "Number of temperatures";
      InputData.close();
  }
  else cout << "Unable to open file";

return 0;


}