% moving_wave_jacketSDK_lax.m

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
nu = u*dx/dt;

Num_ts = min(1000,ceil(15/dt));

% set initial condition
%f = zeros(N,1);
x_space = linspace(x_left,x_right,N);

f_l = 1;
f = f_l*exp(-(x_space.*x_space));
f((x_space < -5) & (x_space > -7)) = 1;

f_tmp = zeros(N,1);

% plot initial condition
plot(x_space,f,'-b');
axis([x_left x_right 0 1.1*f_l]);
title('\bf{Initial Condition}');
drawnow

tic;

f = gdouble(f);
f_tmp = gdouble(f);


addpath('/home/srblair/Dropbox/matlab/jacketSDK/wave1D_lax');

for ts = 1:Num_ts
    
    if(mod(ts,100)==0)
       fprintf('Executing time step number %d.\n',ts);
    end
    
    %f_tmp = 0.5.*(f(x_p)+f(x_m))-(u*dt/(2*dx)).*(f(x_p)-f(x_m));
    f = wave1D_lax(f,u,dt,dx);
    
       
    if(plot_switch==1)
        if(mod(ts,plot_freq)==0)
            plot(x_space,f,'-b')
            axis([x_left x_right 0 1.1*f_l]);
            drawnow
        end
    end
    
end

ex_time = toc;

f = double(f);

plot(x_space,f,'-b');
axis([x_left x_right 0 1.1*f_l]);
title('\bf{Final Condition}');
drawnow

fprintf('Execution time = %g.\n Average time per DOF*update = %g. \n',ex_time, ex_time/(N*Num_ts));
fprintf('DOF updates per second = %g.\n',(N*Num_ts)/ex_time);