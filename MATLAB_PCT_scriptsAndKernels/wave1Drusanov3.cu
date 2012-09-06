__global__ void wave1Drusanov3(double * f_next,double * f_tmp, 
				double * f_in, double nu,
				double omega, int N){
  int tid=threadIdx.x+blockIdx.x*blockDim.x;
  if(tid<N){
    int x_2m=tid-2;
    if(x_2m<0) x_2m+=N;
    int x_m = tid-1;
    if(x_m<0) x_m+=N;

    int x_p = tid+1;
    if(x_p>(N-1)) x_p-=N;

    int x_2p = tid+2;
    if(x_2p>(N-1)) x_2p-=N;

    double f_2m = f_in[x_2m];
    double f_m = f_in[x_m];
    double f = f_in[tid];
    double f_p = f_in[x_p];
    double f_2p = f_in[x_2p];

    f_next[tid]=f-(nu/24.)*(-2.*f_2p+7.*f_p - 7.*f_m+2.*f_2m)
      -(3.*nu/8.)*(f_tmp[x_p]-f_tmp[x_m])
      -(omega/24.)*(f_2p - 4.*f_p + 6.*f - 4.*f_m + f_2m);


  }
}
