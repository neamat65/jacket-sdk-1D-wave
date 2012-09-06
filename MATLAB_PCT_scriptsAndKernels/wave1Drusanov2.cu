__global__ void wave1Drusanov2(double * f_tmp,double * f_nm, 
				double * f_in, double nu, int N){
  int tid = threadIdx.x+blockIdx.x*blockDim.x;
  if(tid<N){
    int x_m = tid-1;
    if(x_m<0) x_m = (N-1);
    f_tmp[tid]=f_in[tid]-(2.*nu/3.)*(f_nm[tid]-f_nm[x_m]);

  }
}
