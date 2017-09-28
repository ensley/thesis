#include <math.h>
#include <R.h>

#define TPB 1024

__device__ double mcCalc(double x, double *d_samps, int S)
{
	double total = 0.0f;
	for (int i = 0; i < S; i++)
	{
		total += cos(x * d_samps[i]);
	}
	return total / S;
}

__global__ void mcKernel(double *d_vec, double *d_samps, double *d_mat, int N, int S)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= N) return;
	const double x = d_vec[i];
	d_mat[i] = mcCalc(x, d_samps, S);
}


// Helper function for using CUDA to add vectors in parallel.
extern "C" void mcCuda(double *vec, double *samps, double *mat, int N, int S)
{

    double *d_vec = 0;
    double *d_samps = 0;
    double *d_mat = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        error("cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output).
    cudaStatus = cudaMalloc((void**)&d_vec, N * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        error("cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&d_samps, S * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        error("cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&d_mat, N * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        error("cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(d_vec, vec, N * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        error("cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(d_samps, samps, S * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        error("cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with TPB threads per block.
    mcKernel<<<(N+TPB-1)/TPB, TPB>>>(d_vec, d_samps, d_mat, N, S);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        error("mcKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        error("cudaDeviceSynchronize returned error code %d after launching mcKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(mat, d_mat, N * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        error("cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(d_vec);
    cudaFree(d_samps);
    cudaFree(d_mat);

    cudaThreadExit();

}
