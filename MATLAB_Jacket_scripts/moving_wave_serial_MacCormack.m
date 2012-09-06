% moving_wave_serial_MacCormack.m


clear
%clc
close('all')

N = 50000;
plot_freq = 200;
plot_switch=0;
u = 1;
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
f_tmp1 = zeros(N,1);

% plot initial condition
plot(x_space,f,'-b');
axis([x_left x_right 0 1.1*f_l]);
title('\bf{Initial Condition}');
drawnow

% time step
tic
for ts = 1:Num_ts
    
    if(mod(ts,100)==0)
       fprintf('Executing time step number %d.\n',ts);
    end
    % slower loop-driven update
    for x = 1:N
       x_p = x+1;
       % enforce periodicity at boundary
       if (x_p >N)
           x_p = 1;
       end
       
       x_m = x-1;
       % enforce periodicity at boundary
       if(x_m < 1)
           x_m = N;
       end
       
       f_tmp1(x)=f(x) - (u*dt/dx)*(f(x_p)-f(x));
            
       
    end
    
    for x = 1:N
        x_p = x+1;
       % enforce periodicity at boundary
       if (x_p >N)
           x_p = 1;
       end
       
       x_m = x-1;
       % enforce periodicity at boundary
       if(x_m < 1)
           x_m = N;
       end
       
       f_tmp(x)=0.5*(f(x)+f_tmp1(x)-(u*dt/dx)*(f_tmp1(x)-f_tmp1(x_m))); 
        
    end
       
    
    f = f_tmp;
   
    if(plot_switch==1)
        if(mod(ts,plot_freq)==0)
            plot(x_space,f,'-b')
            axis([x_left x_right 0 1.1*f_l]);
            drawnow
        end
    end
    
end
ex_time = toc;

plot(x_space,f,'-b');
axis([x_left x_right 0 1.1*f_l]);
title('\bf{Final Condition}');
drawnow

fprintf('Execution time = %g.\n Average time per DOF*update = %g. \n',ex_time, ex_time/(N*Num_ts));