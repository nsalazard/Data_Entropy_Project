#include <iostream>
#include <vector>
#include <random>
#include <cmath>
typedef std::vector<double> vector;
void initial_conditions(vector & data, int N);
void cream_in_coffee(vector & data, int N,int lmin, int lmax,double DELTA);
void evolve(vector & data, vector & prob, int Np, int N, int nsteps,double Xmin,double Ymin,double DELTA,int Npart,int u,int hole);
void entropy(vector & data, int N, int b, double DELTA,int hole,int &Npart);
double grid(vector & data, vector & prob, int Np, int N,int Npart);
void print_screen(const vector & data, int N);
void start_gnuplot(double Xmin,double Ymin);
void print_gnuplot(const vector & data, int N, double DELTA, double Xmin,double Ymin);
int main(int argc, char **argv)
{
	int N = std::atoi(argv[1]); // Matrix  (500)
	int Npart = std::atoi(argv[2]); // Particles = Npart*Npart (20)
	int NSTEPS = std::atoi(argv[3]); //Steps
	int u = std::atoi(argv[4]); //Steps
	double DELTA = 1.0;
	double Xmin = -N/2;
	double Ymin = -N/2;
	int lmin = -Xmin-(Npart/2);
	int lmax = -Xmin+(Npart/2);
	int Np = N/10;
	int hole = 10;
    // declare data structures
  vector potential(N*N); 
	vector prob(Np*Np);
    // set initial and boundary conditions
  initial_conditions(potential, N);
	initial_conditions(prob, Np);
  cream_in_coffee(potential, N, lmin,lmax, DELTA);
    // evolve and print
  evolve(potential,prob,Np, N, NSTEPS, Xmin, Ymin, DELTA,Npart,u, hole);
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
void cream_in_coffee(vector & data, int N,int lmin, int lmax,double DELTA)
{
    for(int ix = lmin; ix < lmax; ++ix) {
        for(int iy = lmin; iy < lmax; ++iy) {
        data[ix*N + iy] = DELTA;
        }
    }  
}
void evolve(vector & data, vector & prob, int Np, int N, int nsteps,double Xmin,double Ymin,double DELTA,int Npart,int u,int hole)
{
	if (u == 0){  //Create a Gif using Gnuplot
    start_gnuplot(Xmin,Ymin);
    print_gnuplot(data, N, DELTA,Xmin,Ymin);
    for(int istep = 0; istep < nsteps; istep += 1) {
        entropy(data, N,istep, DELTA, hole,Npart);
	if (istep%100 == 0) {
		if (Npart == 0) {
			exit (EXIT_FAILURE);
		}
        print_gnuplot(data, N, DELTA,Xmin,Ymin);
			}
    }
	}
	if (u == 1){  //Print the entropy
    //print_screen(data, N);
    std::cout << 0 << "\t" << grid(data,prob,Np, N, Npart) << std::endl;
    initial_conditions(prob, Np);
    for(int istep = 1; istep <= nsteps; istep += 1) {
        entropy(data, N,istep, DELTA, hole,Npart);
        //print_screen(data, N);
        if (istep%100 == 0) {
		if (Npart == 0) {
			exit (EXIT_FAILURE);
		}
		  	std::cout << istep << "\t" << grid(data,prob,Np, N, Npart) << std::endl;
		  	initial_conditions(prob, Np);
    		}
			}
		}
	else{
		
	std::cout << 0 << "\t" << Npart << std::endl;
    	for(int istep = 1; istep <= nsteps; istep += 1) {
        entropy(data, N,istep, DELTA, hole,Npart);
		if (Npart == 0) {
			exit (EXIT_FAILURE);
		}
        //print_screen(data, N);
        	if (istep%100 == 0) {
		  	std::cout << istep << "\t" <<  Npart << std::endl;
    		}
	}	
	}
}
void entropy(vector & data, int N,int b, double DELTA,int hole,int &Npart){
    std::mt19937 gen(b);
    std::uniform_real_distribution<double> dis(0, 4.0);
    for(int ix = 0; ix < N; ++ix) {
        for(int iy = 0; iy < N; ++iy) {
            if (data[ix*N + iy] != 0.0){
                double a = dis(gen);
                //Move up
                if (a < 1.0){
                    if(ix == 0){
                        continue;
                        }
                    else if (data[(ix-1)*N + iy] == 0.0){
                        data[ix*N + iy] = 0.0;
                        data[(ix-1)*N + iy] = DELTA;
                    } 
                    else {
                        continue;
                    }
                }
                //Move right 
                if (1.0 <= a && a < 2.0){
                    if(iy == N-1){
                        continue;
                        }
                    else if (data[ix*N + iy+1] == 0.0){
                            data[ix*N + iy] = 0.0;
                            data[ix*N + iy+1] = DELTA;
                    } 
                    else {
                        continue;
                    } 
                }
                //Move left
                if (2.0 <= a && a < 3.0){
                    if(iy == 0){
                        continue;
                        }
                    else if (data[ix*N + iy-1] == 0.0){
                            data[ix*N + iy] = 0.0;
                            data[ix*N + iy-1] = DELTA;
                    } 
                    else {
                        continue;
                    }
                }
                //Move down
                if (3.0 < a){
                    if(ix == N-1){
			if (N/2 - hole/2 <= iy <= N/2 + hole/2){
				data[ix*N + iy] = 0.0;
				Npart -= 1;
				}
			else{
                        continue;
			}
                        }
                    else if (data[(ix+1)*N + iy] == 0.0){
                        data[ix*N + iy] = 0.0;
                        data[(ix+1)*N + iy] = DELTA;
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
double grid(vector & data, vector & prob, int Np, int N,int Npart){
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
void print_screen(const vector & data, int N){
    for(int ix = 0; ix < N; ++ix) {
        for(int iy = 0; iy < N; ++iy) {
            std::cout << data[ix*N + iy] << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}
void start_gnuplot(double Xmin,double Ymin)
{
    std::cout << "set term gif animate delay 10\n";
    std::cout << "set output 'Entropy_03.gif'\n";
		std::cout << "set xrange ["<<Xmin<<":"<<-Xmin<<".0]\n";
		std::cout << "set yrange ["<<Ymin<<":"<<-Ymin<<".0]\n";
}
void print_gnuplot(const vector & data, int N, double DELTA, double Xmin,double Ymin)
{
    std::cout << "plot '-' w p ls 3 \n";
    for(int ix = 0; ix < N; ++ix) {
        double x =  Xmin + ix*DELTA;
        for(int iy = 0; iy < N; ++iy) {
	 if (data[ix*N + iy] != 0.0){
            double y = Ymin + iy*DELTA;
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
