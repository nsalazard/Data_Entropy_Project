#include "papi.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>

typedef std::vector<int> VEC;
typedef std::vector<int> IVEC;
void fill(VEC &D);
int code_to_be_measured(const int N, const int blockSize, const VEC &A,const VEC &B, VEC &C);
void papi(float real_time, float proc_time, float mflops, long long flpops,float ireal_time,float iproc_time,float imflops, 
long long iflpops, int retval, IVEC & data,const int N, const int blockSize, const VEC &A,const VEC &B, VEC &C);


int main(int argc, char **argv) {
  const int P = std::atoi(argv[1]);
  const int Q = std::atoi(argv[2]);
  const int N = std::atoi(argv[3]);

  // PAPI vars
  float real_time, proc_time, mflops;
  long long flpops;
  float ireal_time, iproc_time, imflops;
  long long iflpops;
  int retval;
  // PERFOMANCE MEASURE
  if (P == 0) { // Blocksize
    // int data[13] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    IVEC data = {1, 2, 4, 8, 16, 32, 64, 128, 256};
    //std::cout<< "Blocksize" << "\t" <<  "Time CPU " << "\t" <<  "Total Flops" << "\t" << "MFLOPS" << "\t" << "C[3]" << "\n";
    if (Q == 0) {
      VEC A(N*N);
      VEC B(N*N);
      VEC C(N*N);
      // initialize matrices
      fill(A);
      fill(B);
      
      for (auto &blocksize : data) {
	      if (blocksize <= N){
        real_time = 0.0;
      proc_time = 0.0;
      mflops = 0.0;
      flpops = 0;
      ireal_time = 0.0;
      iproc_time = 0.0;
      imflops = 0.0;
      iflpops = 0;
      retval = 0;
      papi(real_time,proc_time, mflops, flpops,ireal_time,iproc_time, imflops, iflpops, retval, data,N, blocksize, A,B,C);

	}
	      else {return 0;}
      }
      }
      else if (Q == 1) {
        VEC A(N*N);
        VEC B(N*N);
        VEC C(N*N);
      // initialize matrices
        fill(A);
        fill(B);
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
        papi(real_time,proc_time, mflops, flpops,ireal_time,iproc_time, imflops, iflpops, retval, data, N, blocksize, A,B,C);
      }
    }
  }

    else if (P == 1) { // for (auto & N : data)
      // int data[14] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
      // 8192, 16384};
      std::cout<< "N" << "\t" <<  "Time CPU " << "\t" <<  "Total Flops" << "\t" << "MFLOPS" << "\t" << "C[3]" << "\n";
      IVEC data = {2, 4, 8};
     // int cache = 32;

      for (auto &N : data) {
        VEC A(N*N);
        VEC B(N*N);
        VEC C(N*N);
        fill(A);
        fill(B);
        real_time = 0.0;
        proc_time = 0.0;
        mflops = 0.0;
        flpops = 0;
        ireal_time = 0.0;
        iproc_time = 0.0;
        imflops = 0.0;
        iflpops = 0;
        retval = 0;
        int blocksize = 8;
      papi(real_time,proc_time, mflops, flpops,ireal_time,iproc_time, imflops, iflpops, retval, data, N, blocksize, A,B,C);
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
  int code_to_be_measured(const int N, const int blockSize, const VEC &A,const VEC &B, VEC &C) {
    for (int bii = 0; bii < N; bii += blockSize) {
      for (int bjj = 0; bjj < N; bjj += blockSize) {
        for (int bkk = 0; bkk < N; bkk += blockSize) {
          for (int ii = 0; ii < blockSize; ii++) {
            for (int jj = 0; jj < blockSize; jj++) {
              for (int kk = 0; kk < blockSize; kk++) {
                C[N * (bii + ii) + (bjj + jj)] +=
                    A[N * (bii + ii) + (bkk + kk)] *
                    B[N * (bkk + kk) + (bjj + jj)];
              }
            }
          }
        }
      }
    }

    return 0;
  }

  void fill(VEC &D){
	//std::mt19937 gen(0);
  //std::normal_distribution<int> dis(1,2);
  	for (auto & x : D) {
   	 x = 1;//dis(gen);
  }
    }


void papi(float real_time, float proc_time, float mflops, long long flpops,float ireal_time,float iproc_time,float imflops,
 long long iflpops, int retval, IVEC & data, const int N, const int blocksize, const VEC &A,const VEC &B, VEC &C){



      if ((retval = PAPI_flops_rate(PAPI_FP_OPS, &ireal_time, &iproc_time,
                                    &iflpops, &imflops)) < PAPI_OK) {
        printf("Could not initialise PAPI_flops \n");
        printf(
            "Your platform may not support floating point operation event.\n");
        printf("retval: %d\n", retval);
        exit(1);
      }

      code_to_be_measured(N,blocksize, A,B,C);

      if ((retval = PAPI_flops_rate(PAPI_FP_OPS, &real_time, &proc_time,
                                    &flpops, &mflops)) < PAPI_OK) {
        printf("retval: %d\n", retval);
        exit(1);
      }
      std::cout<< "\t" <<  blocksize  << "\t" <<  proc_time << "\t" <<  flpops << "\t" <<  mflops << "\t" << C[3] << "\n";
      // Do something here, like computing the average of the resulting matrix,
      // to avoid the optimizer deleting the code
      C [3]= {0};

    }
    // End for
