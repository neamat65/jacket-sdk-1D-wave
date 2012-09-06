__global__ void wave1Dlax(double * f_next, double * f, double u, 
			  double dt, double dx, int N){
  int tid = threadIdx.x+blockIdx.x*blockDim.x;
  if(tid<N){
    int x_p = tid+1;
    if(x_p ==N)
      x_p = 0;
    int x_m = tid-1;
    if(x_m<0)
      x_m = N-1;

    
    double f_p = f[x_p];
    double f_m = f[x_m];

    f_next[tid]=0.5*(f_p + f_m) - (u*dt/(2.*dx))*(f_p - f_m);

  }
}
