#include <iostream>
#include <vector>
#include <random>
#include <cmath>

const double DELTA = 1.0;
const int N =  350;
const int Npart =  60; //N particles = 60*60
const int lmin = (N/2)-(Npart/2);
const int lmax = (N/2)+(Npart/2);
const int NSTEPS = 10000;
const int Np = N/10;

typedef std::vector<double> vector;

void initial_conditions(vector & data, int N);
void cream_in_coffee(vector & data, int N);
void evolve(vector & data, vector & prob, int Np, int N, int nsteps);
void entropy(vector & data, int N);
double grid(vector & data, vector & prob, int Np, int N);
void print_screen(const vector & data, int N);
void start_gnuplot(void);
void print_gnuplot(const vector & data, int N);


int main(int argc, char **argv)
{
    // declare data structures
        vector potential(N*N); // [ii, jj] -> ii* + jj
		vector prob(Np*Np);
		//int Np = std::atoi(argv[1]);
    // set initial and boundary conditions
    initial_conditions(potential, N);
		initial_conditions(prob, Np);
    cream_in_coffee(potential, N);

    // evolve and print
    evolve(potential,prob,Np, N, NSTEPS);

    return 0;
}

void initial_conditions(vector & data, int N)
{
    for(int ix = 0; ix < N; ++ix) {
        for(int iy = 0; iy < N; ++iy) {
            data[ix*N + iy] = 0.0;
        }
    }
}
void cream_in_coffee(vector & data, int N)
{
    for(int ix = lmin; ix < lmax; ++ix) {
        for(int iy = lmin; iy < lmax; ++iy) {
        data[ix*N + iy] = 1.0;
        }
    }  
}

void evolve(vector & data, vector & prob, int Np, int N, int nsteps)
{
    start_gnuplot();
    print_gnuplot(data, N);
    //print_screen(data, N);
    //grid(data,prob,Np, N);
    initial_conditions(prob, Np);
    for(int istep = 0; istep < nsteps; istep += 1) {
        entropy(data, N);
        //print_screen(data, N);
		//std::cout << istep << "\t" << grid(data,prob,Np, N) << std::endl;
        //initial_conditions(prob, Np);
        print_gnuplot(data, N);
    }
}

void entropy(vector & data, int N)
{
    std::mt19937 gen(5);
    std::uniform_real_distribution<double> dis(0, 4.0);
    for(int ix = 0; ix < N; ++ix) {
        for(int iy = 0; iy < N; ++iy) {
            if (data[ix*N + iy] != 0.0){
                double a = dis(gen);
                //Move up
                if (a < 1.0){
                    if(ix == 0){
                        continue;
                            //data[ix*N + iy] = 0.0;
                            //data[(ix+1)*N + iy] = 1.0;
                        }
                    else if (data[(ix-1)*N + iy] == 0.0){
                            data[ix*N + iy] = 0.0;
                            data[(ix-1)*N + iy] = 1.0;
                    }
                    else {
                        continue;
                    }
                }
                //Move right 
                if (1.0 <= a && a < 2.0){
                    if(iy == N-1){
                        continue;
                            //data[ix*N + iy] = 0.0;
                            //data[ix*N + iy-1] = 1.0;
                        }
                    else if (data[ix*N + iy+1] == 0.0){
                            data[ix*N + iy] = 0.0;
                            data[ix*N + iy+1] = 1.0;
                    } 
                    else {
                        continue;
                    } 
                }
                //Move left
                if (2.0 <= a && a < 3.0){
                    if(iy == 0){
                        continue;
                            //data[ix*N + iy] = 0.0;
                            //data[ix*N + iy+1] = 1.0;
                        }
                    else if (data[ix*N + iy-1] == 0.0){
                            data[ix*N + iy] = 0.0;
                            data[ix*N + iy-1] = 1.0;
                    } 
                    else {
                        continue;
                    }
                }
                //Move down
                if (3.0 < a){
                    if(ix == N-1){
                        continue;
                            //data[ix*N + iy] = 0.0;
                            //data[(ix-1)*N + iy] = 1.0;
                        }
                    else if (data[(ix+1)*N + iy] == 0.0){
                        data[ix*N + iy] = 0.0;
                        data[(ix+1)*N + iy] = 1.0;
                    } 
                    else {
                        continue;
                    }
                }
            }
            else {
                continue;
            }
        }
    }

}

double grid(vector & data, vector & prob, int Np, int N)
{
	for(int ii = 0; ii < Np; ++ii) {
		for(int jj = 0; jj < Np; ++jj) {
		    for(int ix = (N/Np)*ii; ix < (N/Np)*(ii+1); ++ix) {
                for(int iy = (N/Np)*jj; iy < (N/Np)*(jj+1); ++iy) {
                    prob[ii*Np+jj] += data[ix*N + iy]; 
                }
			}
      prob[ii*Np+jj] = prob[ii*Np+jj]/(Npart*Npart);
         }    
	} 
    double total = 0;
	for(int kk=0; kk< Np*Np; ++kk){
		if (prob[kk] != 0){
		total +=  prob[kk] * - std::log(prob[kk]) ; 
		}
		else {
			continue;
			}
	}
    return total;

}


void print_screen(const vector & data, int N)
{
    for(int ix = 0; ix < N; ++ix) {
        for(int iy = 0; iy < N; ++iy) {
            std::cout << data[ix*N + iy] << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void start_gnuplot(void)
{
    std::cout << "set term gif animate delay 10\n";
    std::cout << "set output 'anim_03.gif'\n";
		std::cout << "set xrange [0.0:"<<N<<".0]\n";

		std::cout << "set yrange [0.0:"<<N<<".0]\n";
		
}

void print_gnuplot(const vector & data, int N)
{
    std::cout << "plot '-' w p ls 2 \n";
    for(int ix = 0; ix < N; ++ix) {
        double x =  ix*DELTA;
        for(int iy = 0; iy < N; ++iy) {
					if (data[ix*N + iy] != 0.0){
            double y = iy*DELTA;
            std::cout << x << "\t" << y << "\n";
        }
				else{
        continue;
        }
				}
        std::cout << "\n";
        }
    std::cout << "e\n";
}
