#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <chrono>

using namespace std;
using Clock = std::chrono::steady_clock;

const int Nx = 4; //Spins in direction x
const int Ny = 4;
const int N2 = 16; // all Spins
int MCSteps = 2000; // number of Metropolis Steps
const int k_B = 1; // Boltzmann constant
//vector<int> h; //magnetic field	
double h;
const int J = 0.2; // ferromagnetic coupling
const int mu = 1; //permeability
double beta; // invers temperatur
double T; // temperatur
int s[Nx][Ny]; //spin
int li[N2] , lj[N2] , ri[N2] , rj[N2]; //spin neightbours
int sum_s;
int i; //spin position in direction x
int j; // spin position in direction y
double dE; // energy difference (E_new - E_old)
double E; // energy
double Boltz_prob; // transmission probability
double Random_prob; // random probability betwenn 0 and 1
double mag; // magnetisation
double med_mag; // mean magnetisation
double med_E; // mean energy
int count1; //counter
int count2; //counter


void Initialization () //inizialisation: random start +1 oder -1
{
	random_device seed;
	mt19937 engine(seed());
	uniform_real_distribution<double> values(0.0, 1.0);	// oh, this is double

	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
			
			if(values(engine) < 0.5){
				s[i][j] = 1;
			}
			else{
			    	s[i][j] = -1;
			}
		}
	}	
}

double Magnetisation() 
{
	int Magnetisation = 0;
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			
			Magnetisation += s[i][j];
		}
	}
	return Magnetisation;
}

double Energy() 
{
	int sumN = 0;
	//int Energy = 0;
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			
		mag += s[i][j];

		int ri = i == Nx-1 ? 0 : i+1;	// periodic condition
		int rj = j == Ny-1 ? 0 : j+1;
		//int li = i == 0 ? Nx-1 : i-1;
		//int lj = j == 0 ? Ny-1 : j-1;
		
		//int neighbor = s[ri][j] + s[li][j] + s[i][rj] + s[i][lj]; 	
		//int neighbor = s[ri][j] + s[i][rj]; 	
		//Energy -=  s[i][j] * (J* neighbor + h);
		sumN += s[i][j] * (s[ri][j] + s[i][rj]);
		 
		}
	}
	//return Energy;	
	return -(J*sumN + h*mu*mag);	
}

void OneMetropolisStep()
{
	random_device seed;
	mt19937 engine(seed());	// already in initialization	/ seed the engine
	uniform_int_distribution<int> value1(0, Nx);	// btwn 0, Nx

	int i = value1(engine);
	
	uniform_int_distribution<int> value2(0, Ny);
	int j = value2(engine);

	// nearest neightbours
	int li;		
	if(i==0)
		li = Nx-1;	// ---(i-1)(i=0)(i+1)--- then at Nx-1! Perioic
	else li = i-1;		// else, ---(i-1)(i)(i+1)---- just

	int ri;
	if(i == Nx-1)		// ---(i-1)(i=Nx-1)(i+1)--- 
		ri = 0;		// else just
	else ri = i+1;

	int lj;			// same as i
	if(j==0)
		lj= Ny-1;
	else lj = j-1;

	int rj;
	if(j == Ny-1)
		rj = 0;
	else rj = j+1;

	// Perioduc function before mag & E, so it looks fine...maybe?

	mag = Magnetisation();
	E = Energy();

	beta = 1/(k_B);
	//sum_s = (s[ri][j] + s[li][j] + s[i][rj] + s[i][lj]);
	//dE =  2 *E;
	dE =  - 2* s[i][j]*(J *sum_s + mu *h);

// E =  -(J*(s[i][j] * (s[ri][j] + s[i][rj])) + B* s[i][j]*1)
	if(dE <= 0){		// flip spin!
		s[i][j] = -s[i][j];
		mag += 2*s[i][j];
		//mag += s[i][j];
		//dE =  - s[i][j]*(J *sum_s + mu *B);	//here
		E += dE;
	}
	else{ 			// make a random number
		random_device seed;
		mt19937 engine(seed());
		uniform_real_distribution<double> values(0.0, 1.0);

		Boltz_prob = exp(-beta* dE);
		Random_prob = values(engine);
	
	if (Random_prob < Boltz_prob){
		s[i][j] = -s[i][j];
		mag += 2*s[i][j];
		//mag += s[i][j];
		//dE =  - s[i][j]*(J *sum_s + mu *B);	//here
		E += dE;
	}
	else{
		s[i][j] = s[i][j];
		mag = mag;
		E = E;
	}
	}
}
void oneMetropolisStepPerSpin()
{
	for (int i = 0; i < N2; i++){
		OneMetropolisStep();
	}
	mag = mag/ N2;
	E = E/ N2;
}

// main function!!
int main(int argc, const char * argv[])
{
	
	if (argc!=2){
		cout << "Need one argument: the output file name" << endl;
		return 1;
	}
	chrono::duration<double, milli> randomise_time, mag_time, E_time; 
		

	ofstream ofs(argv[1]);

	for(double h_i = -2.0; h_i <= 2.1; h_i+= 0.1){
		/*		
		if(h_i < 0.01 and h_i > -0.01){		
			cout << "0 point" << endl;
			h_i = 0.0;		
		}
		*/
		h = h_i;

		med_mag =0;	// reset
		med_E =0;

		for(count1 = 0; count1 < MCSteps; count1++){
		
			auto const t0 = Clock::now();
			Initialization();
			auto const t1 = Clock::now();

			oneMetropolisStepPerSpin();
			med_mag += mag;
			auto const t2 = Clock::now();
			med_E += E;
			auto const t3 = Clock::now();
		
			randomise_time += t1 -t0;
			mag_time += t2 - t1;
			E_time += t3 -t2;	
		}
	
		med_mag = med_mag / MCSteps;
		med_E = med_E / MCSteps;


		med_mag = - med_mag;
		//med_E = - med_E;


		//cout << "mc step = " << MCSteps << endl;
		ofs << h << "   " << med_mag << "   " << med_E << endl;
		// analytical solution:
		//cout << "Randomisation time = " << randomise_time.count() << " ms" << endl;
		//cout << "Magnetisation time = " << mag_time.count() << " ms" << endl;
		//cout << "Energy time = " << E_time.count() << " ms" << endl;
		cout << "The external magenetic field h = "<< h <<'\n';
		cout <<"<m> = "<< med_mag <<'\n';
		cout <<"<E> = "<<med_E<<'\n';
		cout << "-------------------------------------------" << '\n';
	

	}

	return 0;

}







