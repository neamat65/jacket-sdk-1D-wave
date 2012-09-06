__global__ void wave1Drusanov1(double * f_nm, double * f_in, 
				double nu,int N){

  int tid=threadIdx.x+blockIdx.x*blockDim.x;
  if(tid<N){
    int x_p = tid+1;
    if(x_p==N) x_p=0;

    double fp = f_in[x_p];
    double f = f_in[tid];
    f_nm[tid]=0.5*(fp+f)-(nu/3.)*(fp-f);

  }
}
