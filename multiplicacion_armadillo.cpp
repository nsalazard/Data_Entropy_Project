# include <iostream>
# include <cstdio>
# include <cstdlib>
#include <armadillo>
# include <vector>
# include "papi.h"
int code_to_be_measured(const arma::mat & A, const arma::mat & B, arma::mat & C);

typedef std::vector<int> VEC;

int main(int argc, char **argv)
{
  float real_time, proc_time,mflops;
    long long flpops;
    float ireal_time, iproc_time, imflops;
    long long iflpops;
    int retval;
// PERFOMANCE MEASURE
    
   // VEC data = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};
      VEC data = {2, 4, 8,  16, 32, 64, 128, 256};
    
      std::cout<< "N" << "\t" <<  "Time CPU " << "\t" <<  "Total Flops" << "\t" << "MFLOPS" << "\t"<< "C(1,1)" << "\n";
      for (auto &N : data) {
      real_time = 0.0;
      proc_time = 0.0;
      mflops = 0.0;
      flpops = 0;
      ireal_time = 0.0;
      iproc_time = 0.0;
      imflops = 0.0;
      iflpops = 0;
      retval = 0;
    
arma::mat A(N, N, arma::fill::randu);
arma::mat B(N, N, arma::fill::randu);
arma::mat C(N, N);

// start PAPI counters
    if((retval=PAPI_flops_rate(PAPI_FP_OPS,&ireal_time,&iproc_time,&iflpops,&imflops)) < PAPI_OK)
    {
        printf("Could not initialise PAPI_flops \n");
        printf("Your platform may not support floating point operation event.\n");
        printf("retval: %d\n", retval);
        exit(1);
    }
    code_to_be_measured(A, B, C);
    if((retval=PAPI_flops_rate(PAPI_FP_OPS,&real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
    {
        printf("retval: %d\n", retval);
        exit(1);
    }
    
    std::cout <<  N  << "\t" <<  proc_time << "\t" <<  flpops << "\t" <<  mflops << "\t" << C(1,1) << "\n";
    // Do something here, like computing the average of the resulting matrix, to avoid the optimizer deleting the code
      }
           return 0;
}

int code_to_be_measured(const arma::mat & A, const arma::mat & B, arma::mat & C)
{
    C = A*B;
    return 0;
}
