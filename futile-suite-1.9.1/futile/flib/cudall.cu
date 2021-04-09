#include <iostream>
#include <stdio.h>
#include "cufft.h"
#include "cuda.h"
#include "cublas_v2.h"
#include "cuda_runtime_api.h"
#include "config.h"
 
#define DOUBLE

#ifdef DOUBLE
#define Complex  cufftDoubleComplex
#define Real double
#define Transform CUFFT_Z2Z
#define TransformExec cufftExecZ2Z
#else
#define Complex  cufftComplex
#define Real float
#define Transform CUFFT_C2C
#define TransformExec cufftExecC2C
#endif

#define TILE_DIM  8


static const char *_cublasGetErrorString(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
        case CUBLAS_STATUS_NOT_SUPPORTED:
            return "CUBLAS_STATUS_NOT_SUPPORTED";
#if CUDA_VERSION >= 6500
        case CUBLAS_STATUS_LICENSE_ERROR:
            return "CUBLAS_STATUS_LICENSE_ERROR";
#endif
    }
    return "<unknown>";
}

static const char *_cufftGetErrorString(cufftResult error)
{
    switch (error)
    {
        case CUFFT_SUCCESS:
            return "CUFFT_SUCCESS";
        case CUFFT_INVALID_PLAN:
            return "CUFFT_INVALID_PLAN";
        case CUFFT_ALLOC_FAILED:
            return "CUFFT_ALLOC_FAILED";
        case CUFFT_INVALID_TYPE:
            return "CUFFT_INVALID_TYPE";
        case CUFFT_INVALID_VALUE:
            return "CUFFT_INVALID_VALUE";
        case CUFFT_INTERNAL_ERROR:
            return "CUFFT_INTERNAL_ERROR";
        case CUFFT_EXEC_FAILED:
            return "CUFFT_EXEC_FAILED";
        case CUFFT_SETUP_FAILED:
            return "CUFFT_SETUP_FAILED";
        case CUFFT_INVALID_SIZE:
            return "CUFFT_INVALID_SIZE";
        case CUFFT_UNALIGNED_DATA:
            return "CUFFT_UNALIGNED_DATA";
    }
    return "<unknown>";
}


cudaStream_t stream1=NULL;
cublasHandle_t handle1=NULL;
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


#define cufftErrchk(ans) { __cufftAssert((ans), __FILE__, __LINE__); }

inline void __cufftAssert(cufftResult code, const char *file, const int line, bool abort=true)
{
   if(code != CUFFT_SUCCESS) 
   {
      fprintf(stderr, "cufftAssert : %s %s %d.\n",
      _cufftGetErrorString(code), file, line);
      if (abort) exit(-1);
   }
}

#define cublasErrchk(ans) { __cublasAssert((ans), __FILE__, __LINE__); }

inline void __cublasAssert(cublasStatus_t code, const char *file, const int line, bool abort=true)
{
   if(code !=CUBLAS_STATUS_SUCCESS) 
   {
      fprintf(stderr, "cublasAssert : %s %s %d.\n",
      _cublasGetErrorString(code), file, line);
      if (abort) exit(-1);
   }
}

// synchronize blocks
extern "C" void FC_FUNC(synchronize, SYNCHRONIZE)() {
   cudaStreamSynchronize(stream1);
}


// allocate device memory
extern "C" void FC_FUNC(cudamalloc, CUDAMALLOC) (int *size, Real **d_data,int *ierr) {

  *ierr = cudaMalloc((void**)d_data, sizeof(Real)*(*size));
  //errors should be treated in the fortran part
}

// allocate device memory
extern "C" void FC_FUNC(cudamemset, CUDAMEMSET) (Real **d_data, int* value, int* size,int *ierr) {

  *ierr = cudaMemsetAsync((void*)*d_data, *value, sizeof(Real)*(*size),stream1);
}

extern "C" void FC_FUNC(cudafree, CUDAFREE) (Real **d_data) {
  cudaFree(*d_data);
}

// set device memory
extern "C" void FC_FUNC_(reset_gpu_data, RESET_GPU_DATA)(int *size, Real* h_data, Real **d_data){
  cudaMemcpyAsync(*d_data, h_data, sizeof(Real)*(*size),
         cudaMemcpyHostToDevice,stream1);
  gpuErrchk( cudaPeekAtLastError() );
}

// copy data on the card
extern "C" void FC_FUNC_(copy_gpu_data, COPY_GPU_DATA)(int *size, Real** dest_data, Real **send_data){
  cudaMemcpyAsync(*dest_data, *send_data, sizeof(Real)*(*size),
         cudaMemcpyDeviceToDevice,stream1);
  gpuErrchk( cudaPeekAtLastError() );
}


// read device memory
extern "C" void FC_FUNC_(get_gpu_data, GET_GPU_DATA)(int *size, Real *h_data, Real **d_data) {
  cudaMemcpyAsync(h_data, *d_data, sizeof(Real)*(*size),
         cudaMemcpyDeviceToHost,stream1);
  gpuErrchk( cudaPeekAtLastError() );
}

extern "C" void FC_FUNC_(cuda_get_mem_info, CUDA_GET_MEM_INFO)(size_t* freeSize, size_t* totalSize){
 gpuErrchk(cudaMemGetInfo(freeSize,totalSize));
}

// set device memory
extern "C" void FC_FUNC_(poisson_cublas_daxpy, POISSON_CUBLAS_DAXPY)(int *size, const double* alpha, Real** d_x,int* facx, Real ** d_y, int* facy,int* offset_y){

  cublasSetStream(handle1, stream1);
  cublasErrchk(cublasDaxpy(handle1,*size,alpha,*d_x,*facx,*d_y+*offset_y,*facy));
//  gpuErrchk( cudaPeekAtLastError() );
}


extern "C" void FC_FUNC_(cudagetdevicecount, CUDAGETDEVICECOUNT)(int* num_devices){
  gpuErrchk(cudaGetDeviceCount(num_devices));
}

extern "C" void FC_FUNC_(cudasetdevice, CUDASETDEVICE)(int* device){
  gpuErrchk(cudaSetDevice(*device));
}

extern "C" void FC_FUNC_(cudadevicereset, CUDADEVICERESET)(int* device){
  gpuErrchk(cudaDeviceReset());
}

