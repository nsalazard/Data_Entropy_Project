# include <iostream>
# include <cstdio>
# include <cstdlib>
# include <vector>
# include <eigen3/Eigen/Dense>
# include "papi.h"
int code_to_be_measured(const Eigen::MatrixXd & M, Eigen::MatrixXd & N, Eigen::MatrixXd & R);

typedef std::vector<int> VEC;

int main(int argc, char **argv)
{
// PAPI vars
    float real_time, proc_time,mflops;
    long long flpops;
    float ireal_time, iproc_time, imflops;
    long long iflpops;
    int retval;
// PERFOMANCE MEASURE
    
   // VEC data = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};
      
      VEC data = {2, 4, 8,  16, 32, 64, 128, 256};

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
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(N, N);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(N, N);
    Eigen::MatrixXd C(N,N);
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
   
// Do something here, like computing the average of the resulting matrix, to avoid the optimizer deleting the code
     std::cout<< "N" << "\t\t" <<  "Time CPU " << "\t\t" <<  "Total Flops" << "\t\t" << "MFLOPS" << "\t\t"<< "C[3]" << "\n";
     std::cout <<  N  << "\t\t" <<  proc_time << "\t\t" <<  flpops << "\t\t" <<  mflops << "\t\t" << C.sum() << "\n"; 
      }
    return 0;
}

int code_to_be_measured(const Eigen::MatrixXd & M, Eigen::MatrixXd & N, Eigen::MatrixXd & R)
{
    R = M*N;
    return 0;
}
