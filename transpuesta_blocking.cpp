#include "papi.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>

typedef std::vector<double> VEC;
typedef std::vector<int> IVEC;
void fill(VEC &D);
int code_to_be_measured(const int N, const int blockSize, const VEC &M, VEC &MT);
void papi(float real_time, float proc_time, float mflops, long long flpops,float ireal_time,float iproc_time,float imflops, 
long long iflpops, int retval, IVEC & data,const int N, const int blockSize, const VEC &M, VEC &MT, int x);


int main(int argc, char **argv) {
  const int P = std::atoi(argv[1]);
  const int Q = std::atoi(argv[2]);

  // PAPI vars
  float real_time, proc_time, mflops;
  long long flpops;
  float ireal_time, iproc_time, imflops;
  long long iflpops;
  int retval;
  // PERFOMANCE MEASURE
  if (P == 0) { // Blocksize
 
    IVEC data = {1, 2, 4, 8, 16, 32, 64, 128, 256,  512, 1024, 2048, 4096};
    std::cout<< "Blocksize" << "\t" <<  "Time CPU " << "\t" <<  "Total Flops" << "\t" << "MFLOPS" << "\t" << "MT[3]" << "\n";
    if (Q == 0) {
      int N = 2048;
      VEC M(N*N);
      VEC MT(N*N);
      // initialize matrices
      fill(M);
      
      for (auto &blocksize : data) {
	      if (blocksize > N) return 0;
      real_time = 0.0;
      proc_time = 0.0;
      mflops = 0.0;
      flpops = 0;
      ireal_time = 0.0;
      iproc_time = 0.0;
      imflops = 0.0;
      iflpops = 0;
      retval = 0;
      papi(real_time,proc_time, mflops, flpops,ireal_time,iproc_time, imflops, iflpops, retval, data,N, blocksize, M,MT, blocksize);

				}
      }
      else if (Q == 1) {
      int N = 4096;
      VEC M(N*N);
      VEC MT(N*N);
      // initialize matrices
      fill(M);
        for (auto &blocksize : data) {
          real_time = 0.0;
      proc_time = 0.0;
      mflops = 0.0;
      flpops = 0;
      ireal_time = 0.0;
      iproc_time = 0.0;
      imflops = 0.0;
      iflpops = 0;
      retval = 0;
        papi(real_time,proc_time, mflops, flpops,ireal_time,iproc_time, imflops, iflpops, retval, data, N, blocksize, M,MT, blocksize);
      }
    }
  }

    else if (P == 1) { // for (auto & N : data)
    
      std::cout<< "N" << "\t" <<  "Time CPU " << "\t" <<  "Total Flops" << "\t" << "MFLOPS" << "\t" << "MT[3]" << "\n";
      IVEC data = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,2048, 4096, 8192, 16384};
     

      for (auto &N : data) {
        VEC M(N*N);
        VEC MT(N*N);
      // initialize matrices
        fill(M);
        real_time = 0.0;
        proc_time = 0.0;
        mflops = 0.0;
        flpops = 0;
        ireal_time = 0.0;
        iproc_time = 0.0;
        imflops = 0.0;
        iflpops = 0;
        retval = 0;
        int blocksize = 64;
	      if (64 > N){
	      blocksize = N;
	      }
      papi(real_time,proc_time, mflops, flpops,ireal_time,iproc_time, imflops, iflpops, retval, data, N, blocksize, M,MT, N);
      }


    }
    else {
      std::cout << "Not a valid argument" << std::endl;
      return 0;
    }


    // start PAPI counters

    return 0;
  }

  // Functions
  int code_to_be_measured(const int N, const int blockSize, const VEC &M, VEC &MT) {
    for (int ii = 0; ii<N; ii+=blockSize) {
    for (int jj=0; jj<N; jj+=blockSize) {
      for(int bi = ii; bi < ii +blockSize ; bi++){
        for(int bj = jj; bj < jj + blockSize; bj++){
      MT [bi*N+bj] =  M [bj*N+bi]*12.50/6.25*0,5;
        }
      }
    }
  }
  return 0;
  }

  void fill(VEC &D){
	std::mt19937 gen(0);
  std::normal_distribution<double> dis(1.0, 1.5);
  	for (auto & x : D) {
   	 x = dis(gen);
  }
    }

void papi(float real_time, float proc_time, float mflops, long long flpops,float ireal_time,float iproc_time,float imflops, 
long long iflpops, int retval, IVEC & data,const int N, const int blockSize, const VEC &M, VEC &MT, int x){



      if ((retval = PAPI_flops_rate(PAPI_FP_OPS, &ireal_time, &iproc_time,
                                    &iflpops, &imflops)) < PAPI_OK) {
        printf("Could not initialise PAPI_flops \n");
        printf(
            "Your platform may not support floating point operation event.\n");
        printf("retval: %d\n", retval);
        exit(1);
      }

      code_to_be_measured(N,blockSize, M,MT);

      if ((retval = PAPI_flops_rate(PAPI_FP_OPS, &real_time, &proc_time,
                                    &flpops, &mflops)) < PAPI_OK) {
        printf("retval: %d\n", retval);
        exit(1);
      }
      std::cout<< x  << "\t" <<  proc_time << "\t" <<  flpops << "\t" <<  mflops << "\t" << MT[3] << "\n";
      // Do something here, like computing the average of the resulting matrix,
      // to avoid the optimizer deleting the code
      MT [3]= {0.0};

    }
    // End for
