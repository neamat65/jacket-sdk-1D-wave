__global__ void wave1Dmac2(double * f_next, double * f_tmp1,
				double * f_in, double u, double dt,
				double dx, int N){
  int tid = threadIdx.x+blockIdx.x*blockDim.x;
  if(tid<N){
    int x_m = tid-1;
    if(x_m <0) x_m = N-1;

    double ft1_tmp = f_tmp1[tid];
    f_next[tid]=0.5*(f_in[tid]+ft1_tmp - u*(dt/dx)*(ft1_tmp-f_tmp1[x_m]));

  }
}
