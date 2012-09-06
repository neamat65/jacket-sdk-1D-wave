% moving_wave_PCT_rusanov.m

clear
%clc
close('all')

N = 1000000;
u = 1;

plot_freq = 200;
plot_switch = 0;
x_left = -10;
x_right = 10;
x_space = linspace(x_left,x_right,N);

dx = x_space(2)-x_space(1);
dt = 0.6*(dx)/u;
nu = u*dt/dx;
omega = (4*nu*nu+1)*(4-nu*nu)/5;
Num_ts = min(5000,ceil(15/dt));
% set initial condition
%f = zeros(N,1);
x_space = linspace(x_left,x_right,N);

f_l = 1;
f = f_l*exp(-(x_space.*x_space));
f((x_space < -5) & (x_space > -7)) = 1;


f_tmp = zeros(N,1);
f_tmp1 = zeros(N,1);
f_nm = zeros(N,1);

% plot initial condition
plot(x_space,f,'-b');
axis([x_left x_right 0 1.1*f_l]);
title('\bf{Initial Condition}');
drawnow

tic;


f_nm=gpuArray(f_nm);
f_tmp = gpuArray(f_tmp);
f_tmp1= gpuArray(f_tmp1);
f = gpuArray(f);

TPB=128;
wave1Drusanov1 = parallel.gpu.CUDAKernel('wave1Drusanov1.ptx','wave1Drusanov1.cu');
wave1Drusanov2 = parallel.gpu.CUDAKernel('wave1Drusanov2.ptx','wave1Drusanov2.cu');
wave1Drusanov3 = parallel.gpu.CUDAKernel('wave1Drusanov3.ptx','wave1Drusanov3.cu');

wave1Drusanov1.ThreadBlockSize=[TPB 1 1];
wave1Drusanov2.ThreadBlockSize=[TPB 1 1];
wave1Drusanov3.ThreadBlockSize=[TPB 1 1];

wave1Drusanov1.GridSize=[ceil(N/TPB) 1];
wave1Drusanov2.GridSize=[ceil(N/TPB) 1];
wave1Drusanov3.GridSize=[ceil(N/TPB) 1];



for ts = 1:Num_ts
    
    if(mod(ts,500)==0)
        fprintf('Executing time step number %d.\n',ts);
    end
    
    %f_nm = (0.5).*(f(x_p) + f) - (nu/3).*(f(x_p) - f);
    f_nm = feval(wave1Drusanov1,f_nm,f,nu,int32(N));
    %f_tmp = f - (2*nu/3).*(f_nm - f_nm(x_m));
    f_tmp=feval(wave1Drusanov2,f_tmp,f_nm,f,nu,int32(N));
%     f = f - (nu/24).*(-2*f(x_2p) + 7*f(x_p) - 7*f(x_m) + 2*f(x_2m)) ...
%         -(3*nu/8).*(f_tmp(x_p) - f_tmp(x_m)) ...
%         -(omega/24).*(f(x_2p) - 4*f(x_p) + 6*f - 4*f(x_m) + f(x_2m));
    f_tmp1=feval(wave1Drusanov3,f_tmp1,f_tmp,f,nu,omega,int32(N));
    
    f=f_tmp1;
    
    if(plot_switch==1)
        if(mod(ts,plot_freq)==0)
            plot(x_space,f,'-b')
            axis([x_left x_right 0 1.1*f_l]);
            drawnow
        end
    end
    
end
f = double(f);
ex_time = toc;

plot(x_space,f,'-b');
axis([x_left x_right 0 1.1*f_l]);
title('\bf{Final Condition}');
drawnow

fprintf('Execution time = %g.\n Average time per DOF*update = %g. \n',ex_time, ex_time/(N*Num_ts));