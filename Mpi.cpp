#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <mpi.h>

typedef std::vector<double> vector;

void initial_conditions(vector & data, int nx, int ny, int pid, int np);
void cream_in_coffee(vector & data, int nx,int ny, int Npart,double DELTA, int pid, int np);
void evolve(vector & data, vector & prob , int nx, int ny, int nsteps,double Xmin,double Ymin,double DELTA,int Np,int u, int pid, int np);
void entropy(vector & data, int nx, int ny,int b, double DELTA,int pid, int np);
double grid(vector & data, vector & prob, int Ng , int nx, int ny,int Np);
void print_screen(const vector & data, int nx, int ny);
void print_screen(const vector & data, int nx, int ny, int pid, int np);
void start_gnuplot(double Xmin,double Ymin,int N,int pid, int np);
void print_gnuplot_slice(const double * data, int nx, int ny,double Xmin,double Ymin,double DELTA, int pid, int np);
void print_gnuplot(const vector & data, int nx, int ny,double Xmin,double Ymin, double DELTA,int istep, int pid, int np);
void communication(vector & data, int nx, int ny, int pid, int np);
double spread(vector & data, int N, double DELTA, double Xmin,double Ymin, int istep, int Np);

int main(int argc, char **argv)
{
	int N = std::atoi(argv[1]); // Matrix  (500)
	int Npart = std::atoi(argv[2]); // Particles = Npart*Npart (20)
	int NSTEPS = std::atoi(argv[3]); //Steps
	int u = std::atoi(argv[4]); //Steps
    // MPI //

    int pid = 0, nproc = 0;
    // init mpi environment
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int NY = N;
    int NX = N/nproc;
    int NXlocal = NX + 2;
    double DELTA = 1.0;
	double Xmin = -NX/2;
	double Ymin = -NY/2;
	int lmin = -Xmin-(Npart/2);
	int lmax = -Xmin+(Npart/2);
	int Np = Npart * Npart;



    // declare data structures
  vector matrix(NXlocal*NY);
	vector prob(NX);

    // set initial and boundary conditions
    initial_conditions(matrix, NXlocal, NY, pid, nproc);
    cream_in_coffee(matrix,NXlocal, NY, Npart, DELTA, pid, nproc);
    communication(matrix, NXlocal, NY, pid, nproc);

    // evolve and print
    evolve(matrix, prob, NXlocal,NY, NSTEPS,Xmin,Ymin, DELTA,Np,u,pid, nproc);
    //print_screen(matrix, NXlocal, NY, pid, nproc);

    // close mpi environment
    MPI_Finalize();




  return 0;
}

void initial_conditions(vector & data, int nx, int ny, int pid, int np)
{
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            data[ix*ny + iy] = 0.0;
        }
    }
}

void cream_in_coffee(vector & data, int nx,int ny, int Npart,double DELTA, int pid, int np)
{
    int xmin = (nx/2)-(Npart/2);
    int xmax = (nx/2)+(Npart/2);
    int ymin = (ny/2)-(Npart/2);
    int ymax = (ny/2)+(Npart/2);
    if (pid == np/2) {
    for(int ix = xmin; ix < xmax; ++ix) {
        for(int iy = ymin; iy < ymax; ++iy) {
        data[ix*ny + iy] = DELTA;
        }
    }  
    }
}

void evolve(vector & data, vector & prob , int nx, int ny, int nsteps,double Xmin,double Ymin,double DELTA,int Np,int u, int pid, int np)
{
	double s = 5;
	if (u == 0){  //Create a Gif using Gnuplot
    start_gnuplot(Xmin,Ymin,ny, pid,np);
    print_gnuplot( data, nx, ny, Xmin,Ymin,DELTA,0, pid, np);
    for(int istep = 1; istep <= nsteps; istep += 1) {

        entropy(data,nx, ny, istep, DELTA, pid,np);
	communication(data, nx, ny, pid, np);

	if (istep%100 == 0) {
        print_gnuplot( data, nx, ny, Xmin,Ymin,DELTA,istep, pid, np);
			}
    }
	}
	
    /*
	if (u == 1){  //Print the entropy
    //print_screen(data, N);
    std::cout << 0 << "\t" << grid(data,prob,Ng , N, Np) << std::endl;
    initial_conditions(prob, Ng );
    double a = 0;
    for(int istep = 1; istep <= nsteps; istep += 1) {
        entropy(data,nx, ny, istep, DELTA, pid,np);
        //print_screen(data, N);
        if (istep%100 == 0) {
		  	std::cout << istep << "\t" << grid(data,prob,Ng , N, Np) << std::endl;
		  	initial_conditions(prob, Ng );
    		}
			}
		}
        
    else if (u == 2 || u ==3){  //Print r  // 
    //print_screen(data, N);
		double r = 0.0;
    double a = 0;
    std::cout << 0 << "\t" << spread( data,N, DELTA, Xmin, Ymin,0,Np) << std::endl;
    for(int istep = 1; istep <= nsteps; istep += 1) {
        entropy(data,nx, ny, istep, DELTA, pid,np);
				r = spread( data,N, DELTA, Xmin, Ymin,istep,Np);
        //print_screen(data, N);
				 if (istep%100 == 0) {
		  	std::cout << istep << "\t" << r << std::endl;
				 }
			}
		}
*/
}

void entropy(vector & data, int nx, int ny,int b, double DELTA,int pid, int np){
	std::mt19937 gen(b);
    std::uniform_real_distribution<double> dis(0, 4.0);
    for(int ix = 1; ix <= nx-2; ++ix) {
        for(int iy = 1; iy < ny; ++iy) {
            if (data[ix*ny + iy] != 0.0){
                if (0 == pid and 1 == ix) continue;
                if (np-1 == pid and nx-2 == ix) continue;
                double a = dis(gen);
                //Move up
                if (a < 1.0){
                    if(ix == 1){
                        continue;
                        }
                    else if (data[(ix-1)*ny + iy] == 0.0){
                        data[ix*ny + iy] = 0.0;
                        data[(ix-1)*ny + iy] = DELTA;
                    } 
                    else {
                        continue;
                    }
                }
                //Move right 
                if (1.0 <= a && a < 2.0){
                    if(iy == ny-1){
                        continue;
                        }
                    else if (data[ix*ny + iy+1] == 0.0){
                            data[ix*ny + iy] = 0.0;
                            data[ix*ny + iy+1] = DELTA;
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
                    else if (data[ix*ny + iy-1] == 0.0){
                            data[ix*ny + iy] = 0.0;
                            data[ix*ny + iy-1] = DELTA;
                    } 
                    else {
                        continue;
                    }
                }
                //Move down
                if (3.0 < a){
                    if(ix == nx-2){
                        continue;
                        }
                    else if (data[(ix+1)*ny + iy] == 0.0){
                        data[ix*ny + iy] = 0.0;
                        data[(ix+1)*ny + iy] = DELTA;
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

double grid(vector & data, vector & prob, int Ng , int nx, int ny,int Np){
    int kk = 0;
	for(int ix = 1; ix <= nx-2; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            prob[kk] += data[ix*ny + iy];  
         } 
         prob[kk] = prob[kk]/(Np);
         kk += 1;   
	} 
    double total = 0;
	for(int kk=0; kk< nx-2; ++kk){
		if (prob[kk] != 0){
		total +=  prob[kk] * - std::log(prob[kk]) ; 
		}
		else {
			continue;
			}
	}
    return total;
}

void print_screen(const vector & data, int nx, int ny)
{
    for(int ix = 0; ix < nx; ++ix) {
        if (0 == ix or nx-1 == ix) {
            std::cout << "G: " ;
        } else {
            std::cout << "   " ;
        }
        for(int iy = 0; iy < ny; ++iy) {
            std::cout << data[ix*ny + iy] << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void print_screen(const vector & data, int nx, int ny, int pid, int np)
{
    int tag = 0;
    if (0 == pid) {
        print_screen(data, nx, ny);
        std::vector<double> buffer(nx*ny);
        for (int src = 1; src < np; ++src) {
            MPI_Recv(&buffer[0], nx*ny, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            print_screen(buffer, nx, ny);
        }
    } else {
        int dest = 0;
        MPI_Send(&data[0], nx*ny, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    }
}

void start_gnuplot(double Xmin,double Ymin,int N,int pid, int np)
{
    if (0 == pid) {

        std::cout << "set term gif animate delay 10\n";
        std::cout << "set output 'Anim_MPI_"<<N<<".gif'\n";
		std::cout << "set xrange ["<<Xmin<<":"<<-Xmin<<".0]\n";
		std::cout << "set yrange ["<<Ymin<<":"<<-Ymin<<".0]\n";

    }
}

void print_gnuplot_slice(const double * data, int nx, int ny,double Xmin,double Ymin,double DELTA, int pid, int np)
{
    for(int ix = 0; ix < nx; ++ix) {
        double x = pid*DELTA*nx + Xmin + ix*DELTA;
        for(int iy = 0; iy < ny; ++iy) {
            double y = Ymin + iy*DELTA;
            std::cout << x << "  " << y << "  " << data[ix*ny + iy] << "\n";
        }
        std::cout << "\n";
    }
}

void print_gnuplot(const vector & data, int nx, int ny,double Xmin,double Ymin, double DELTA,int istep, int pid, int np)
{
    int tag = 0;
    if (0 == pid) {
        std::cout << "plot '-' w p ls 3  t '"<< istep <<"' \n";
	print_gnuplot_slice(&data[ny], nx-2, ny, Xmin,Ymin,DELTA, pid, np);
        std::vector<double> buffer(nx*ny);
        for (int src = 1; src < np; ++src) {
            MPI_Recv(&buffer[0], nx*ny, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            print_gnuplot_slice(&buffer[ny], nx-2, ny, Xmin,Ymin,DELTA, src, np);
        }
        std::cout << "e\n";
    } else {
        int dest = 0;
        MPI_Send(&data[0], nx*ny, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    }
}

void communication(vector & data, int nx, int ny, int pid, int np)
{
    int tag = 0;
    if (pid != 0) {
        MPI_Send(&data[1*ny], ny*1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD);
        MPI_Recv(&data[(0)*ny], ny*1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (np-1 != pid) {
        MPI_Recv(&data[(nx-1)*ny], ny*1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&data[(nx-2)*ny], ny*1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD);
    }


}

double spread(vector & data, int N, double DELTA, double Xmin,double Ymin, int istep, int Np)
{
	double r = 0.0;
    for(int ix = 0; ix < N; ++ix) {
        double x =  Xmin + ix*DELTA;
        for(int iy = 0; iy < N; ++iy) {
					if (data[ix*N + iy] != 0.0){
            double y = Ymin + iy*DELTA;
						r += (x*x) + (y*y);
        }
				else{
        continue;
        }
				}
		}
		return r/Np; // Añadir el cuadrado cuando haya terminado std::sqrt(r/Np)
}
