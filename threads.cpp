#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <cmath>

// receives a vector, a set of indexes, and stores the vector sum on result
void func(const std::vector<double> & data, int iimin, int iimax, double & result);


int main(int argc, char **argv) {
  const int N = std::atoi(argv[1]);
  const int NTHREADS = std::atoi(argv[2]);
  
  std::vector<double> values(N);
  std::iota(values.begin(), values.end(), 1); // fill the vector with values 1, 2, 3 , ..., N-1 

  /*
  // only one thread
  double val1 = 0.0;
  std::thread th1(&func, values, 0, N, std::ref(val1));
  th1.join();
  std::cout << val1 << std::endl;
  */

  // several threads
  std::vector<std::thread> mythreads(NTHREADS);
  std::vector<double> vals(NTHREADS, 0.0); // store here partial sums
  for (int ii = 0; ii < NTHREADS; ++ii) {
    const int localsize = N/NTHREADS; // how much size for each thread, NTHREADS must be divisor of N
    const int iimin = ii * localsize; // FILL HERE minimum global index for this thread
    const int iimax = (ii+1) * localsize; // FILL HERE maximim global index for this thread
    mythreads[ii] = std::thread(&func, values, iimin, iimax, std::ref(vals[ii]));
  }
  for (int ii = 0; ii < NTHREADS; ++ii) {
    mythreads[ii].join();
  }
  std::cout << std::accumulate(vals.begin(), vals.end(), 0.0) << "\n"; // print the accumulated sum
  
  return 0;
}

void func(const std::vector<double> & data, int iimin, int iimax, double & result)
{
  result = 0;
  for (int ii = iimin; ii < iimax; ++ii) {
    result += std::sqrt(std::sqrt(data[ii]));
  }
}
