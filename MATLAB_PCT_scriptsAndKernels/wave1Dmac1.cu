__global__ void wave1Dmac1(double * f_tmp1, double * f_in,
				double u, double dt, double dx,
				int N){
  int tid = threadIdx.x+blockIdx.x*blockDim.x;
  if(tid<N){
    int x_p = tid+1;
    if(x_p == N) x_p = 0;
    
    double f_tmp = f_in[tid];
    f_tmp1[tid]= f_tmp - u*(dt/dx)*(f_in[x_p] - f_tmp);

  }
}
