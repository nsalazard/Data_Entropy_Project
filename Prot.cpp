#include <cmath>
#include <iostream>
#include <random>
#include <vector>

typedef std::vector<double> vector;

void initial_conditions(vector &data, int N);
void cream_in_coffee(vector &data, int N, int lmin, int lmax, double DELTA);
void evolve(vector &data, vector &prob, int Ng, int N, int nsteps, double Xmin,
            double Ymin, double DELTA, int Np, int u);
void entropy(vector &data, int N, int b, double DELTA);
double grid(vector &data, vector &prob, int Ng, int N, int Np);
void print_screen(const vector &data, int N);
void start_gnuplot(double Xmin, double Ymin, int N);
void print_gnuplot(const vector &data, int N, double DELTA, double Xmin,
                   double Ymin, int istep);
double spread(vector &data, int N, double DELTA, double Xmin, double Ymin,
              int istep, int Np);

int main(int argc, char **argv) {
  int N = std::atoi(argv[1]);      // Matrix  (500)
  int Npart = std::atoi(argv[2]);  // Particles = Npart*Npart (20)
  int NSTEPS = std::atoi(argv[3]); // Steps
  int u = std::atoi(argv[4]);      // Steps
  double DELTA = 1.0;
  double Xmin = -N / 2;
  double Ymin = -N / 2;
  int lmin = -Xmin - (Npart / 2);
  int lmax = -Xmin + (Npart / 2);
  int Ng = N / 10;
  int Np = Npart * Npart;
  // declare data structures
  vector matrix(N * N);
  vector prob(Ng * Ng);
  // set initial and boundary conditions
  initial_conditions(matrix, N);
  if (u == 3) {
    matrix[N / 2 * N + N / 2] = DELTA;
  }
  if (u != 3) {
    initial_conditions(prob, Ng);
    cream_in_coffee(matrix, N, lmin, lmax, DELTA);
  }
  // evolve and print
  evolve(matrix, prob, Ng, N, NSTEPS, Xmin, Ymin, DELTA, Np, u);
  return 0;
}
void initial_conditions(vector &data, int N) {
  for (int ix = 0; ix < N; ++ix) {
    for (int iy = 0; iy < N; ++iy) {
      data[ix * N + iy] = 0.0;
    }
  }
}
void cream_in_coffee(vector &data, int N, int lmin, int lmax, double DELTA) {
  for (int ix = lmin; ix < lmax; ++ix) {
    for (int iy = lmin; iy < lmax; ++iy) {
      data[ix * N + iy] = DELTA;
    }
  }
}
void evolve(vector &data, vector &prob, int Ng, int N, int nsteps, double Xmin,
            double Ymin, double DELTA, int Np, int u) {
  double s = 5;
  if (u == 0) { // Create a Gif using Gnuplot
    start_gnuplot(Xmin, Ymin, N);
    print_gnuplot(data, N, DELTA, Xmin, Ymin, 0);
    for (int istep = 1; istep <= nsteps; istep += 1) {
      entropy(data, N, istep, DELTA);
      if (istep >= 16000) {
      if (istep % 1000 == 0) {
        print_gnuplot(data, N, DELTA, Xmin, Ymin, istep);
      }
      }
    }
  }
  if (u == 1) { // Print the entropy
    // print_screen(data, N);
    std::cout << 0 << "\t" << grid(data, prob, Ng, N, Np) << std::endl;
    initial_conditions(prob, Ng);
    double a = 0;
    for (int istep = 1; istep <= nsteps; istep += 1) {
      entropy(data, N, istep, DELTA);
      // print_screen(data, N);
      if (istep % 100 == 0) {
        std::cout << istep << "\t" << grid(data, prob, Ng, N, Np) << std::endl;
        initial_conditions(prob, Ng);
      }
    }
  } else if (u == 2 || u == 3) { // Print r  //
    // print_screen(data, N);
    double r = 0.0;
    double a = 0;
    std::cout << 0 << "\t" << spread(data, N, DELTA, Xmin, Ymin, 0, Np)
              << std::endl;
    for (int istep = 1; istep <= nsteps; istep += 1) {
      entropy(data, N, istep, DELTA);
      r = spread(data, N, DELTA, Xmin, Ymin, istep, Np);
      // print_screen(data, N);
      if (istep % 100 == 0) {
        std::cout << istep << "\t" << r << std::endl;
      }
    }
  }
}

void entropy(vector &data, int N, int b, double DELTA) {
  std::mt19937 gen(b);
  std::uniform_real_distribution<double> dis(0, 4.0);
  for (int ix = 0; ix < N; ++ix) {
    for (int iy = 0; iy < N; ++iy) {
      if (data[ix * N + iy] != 0.0) {
        double a = dis(gen);
        // Move up
        if (a < 1.0) {
          if (ix == 0) {
            continue;
            // data[ix*N + iy] = 0.0;
            // data[(ix+1)*N + iy] = 1.0;
          } else if (data[(ix - 1) * N + iy] == 0.0) {
            data[ix * N + iy] = 0.0;
            data[(ix - 1) * N + iy] = DELTA;
          } else {
            continue;
          }
        }
        // Move right
        if (1.0 <= a && a < 2.0) {
          if (iy == N - 1) {
            continue;
            // data[ix*N + iy] = 0.0;
            // data[ix*N + iy-1] = 1.0;
          } else if (data[ix * N + iy + 1] == 0.0) {
            data[ix * N + iy] = 0.0;
            data[ix * N + iy + 1] = DELTA;
          } else {
            continue;
          }
        }
        // Move left
        if (2.0 <= a && a < 3.0) {
          if (iy == 0) {
            continue;
            // data[ix*N + iy] = 0.0;
            // data[ix*N + iy+1] = 1.0;
          } else if (data[ix * N + iy - 1] == 0.0) {
            data[ix * N + iy] = 0.0;
            data[ix * N + iy - 1] = DELTA;
          } else {
            continue;
          }
        }
        // Move down
        if (3.0 < a) {
          if (ix == N - 1) {
            continue;
            // data[ix*N + iy] = 0.0;
            // data[(ix-1)*N + iy] = 1.0;
          } else if (data[(ix + 1) * N + iy] == 0.0) {
            data[ix * N + iy] = 0.0;
            data[(ix + 1) * N + iy] = DELTA;
          } else {
            continue;
          }
        }
      } else {
        continue;
      }
    }
  }
}

double grid(vector &data, vector &prob, int Ng, int N, int Np) {
  for (int ii = 0; ii < Ng; ++ii) {
    for (int jj = 0; jj < Ng; ++jj) {
      for (int ix = (N / Ng) * ii; ix < (N / Ng) * (ii + 1); ++ix) {
        for (int iy = (N / Ng) * jj; iy < (N / Ng) * (jj + 1); ++iy) {
          prob[ii * Ng + jj] += data[ix * N + iy];
        }
      }
      prob[ii * Ng + jj] = prob[ii * Ng + jj] / (Np);
    }
  }
  double total = 0;
  for (int kk = 0; kk < Ng * Ng; ++kk) {
    if (prob[kk] != 0) {
      total += prob[kk] * -std::log(prob[kk]);
    } else {
      continue;
    }
  }
  return total;
}
void print_screen(const vector &data, int N) {
  for (int ix = 0; ix < N; ++ix) {
    for (int iy = 0; iy < N; ++iy) {
      std::cout << data[ix * N + iy] << "  ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}
void start_gnuplot(double Xmin, double Ymin, int N) {
  std::cout << "set term gif animate delay 10\n";
  std::cout << "set output 'Entropy_" << N << ".gif'\n";
  std::cout << "set xrange [" << Xmin << ":" << -Xmin << ".0]\n";
  std::cout << "set yrange [" << Ymin << ":" << -Ymin << ".0]\n";
}
void print_gnuplot(const vector &data, int N, double DELTA, double Xmin,
                   double Ymin, int istep) {
  std::cout << "plot '-' w p ls 3  t '" << istep << "' \n";
  for (int ix = 0; ix < N; ++ix) {
    double x = Xmin + ix * DELTA;
    for (int iy = 0; iy < N; ++iy) {
      if (data[ix * N + iy] != 0.0) {
        double y = Ymin + iy * DELTA;
        std::cout << x << "\t" << y << "\n";
      } else {
        continue;
      }
    }
    std::cout << "\n";
  }
  std::cout << "e\n";
}

double spread(vector &data, int N, double DELTA, double Xmin, double Ymin,
              int istep, int Np) {
  double r = 0.0;
  for (int ix = 0; ix < N; ++ix) {
    double x = Xmin + ix * DELTA;
    for (int iy = 0; iy < N; ++iy) {
      if (data[ix * N + iy] != 0.0) {
        double y = Ymin + iy * DELTA;
        r += (x * x) + (y * y);
      } else {
        continue;
      }
    }
  }
  return std::sqrt(r)/Np;
}
