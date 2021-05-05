#include "papi.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

typedef std::vector<double> VEC;

int code_to_be_measured(const int N, const int blockSize, const VEC &A,const VEC &B, VEC &C);
void fill(VEC & data);
void papi(float real_time, float proc_time, float mflops, long long flpops,float ireal_time,float iproc_time,float imflops, long long iflpops, int retval, VEC & data,const int N, const int blockSize, const VEC &A,const VEC &B, VEC &C);


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
    // int data[13] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    VEC data(4) = {1, 2, 4, 8};
    if (Q == 0) {
      // int N = 2048;
      int N = 2048;
      VEC A(N*N);
      VEC B(N*N);
      VEC C(N*N);
      // initialize matrices
      fill(A);
      fill(B);
      for (auto &Var : data) {
      papi(real_time,proc_time, mflops, flpops,ireal_time,iproc_time, imflops, iflpops, retval, data,N, Var, A,B,C);

				}
      }
      if (Q == 1) {
        int N = 4096;
      VEC A(N*N);
      VEC B(N*N);
      VEC C(N*N);
      // initialize matrices
      fill(A);
      fill(B);
      for (auto &Var : data) {
      papi(real_time,proc_time, mflops, flpops,ireal_time,iproc_time, imflops, iflpops, retval, data, N, Var, A,B,C);
      }
    }
		/*
    else if (P == 1) { // for (auto & N : data)
      // int data[14] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
      // 8192, 16384};
      int data[2] = {4, 8};
      int blocksize = 2; //Ejemplo para pruebas
    }
    else {
      std::cout << "Not a valid argument" << std::endl;
      return 0;
    }
		*/

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

  void fill(VEC data){
    for (long unsigned int ii = 0; ii < data.size(); ++ii) {
      data[ii] = ii;
    }
  }


void papi(float real_time, float proc_time, float mflops, long long flpops,float ireal_time,float iproc_time,float imflops, long long iflpops, int retval, VEC & data,const int N, const int blockSize, const VEC &A,const VEC &B, VEC &C){


      real_time = 0.0;
      proc_time = 0.0;
      mflops = 0.0;
      flpops = 0;
      ireal_time = 0.0;
      iproc_time = 0.0;
      imflops = 0.0;
      iflpops = 0;
      retval = 0;


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
      printf("Real_time: %f Proc_time: %f Total flpops: %lld MFLOPS: %f\n",
             real_time, proc_time, flpops, mflops);
      // Do something here, like computing the average of the resulting matrix,
      // to avoid the optimizer deleting the code
      printf("%.15e\n", C[0]);
    }
    // End for

