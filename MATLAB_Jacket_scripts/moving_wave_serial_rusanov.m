% moving_wave_serial_rusanov.m


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
nu = u*dt/dx;

Num_ts = min(1000,ceil(15/dt));


% set initial condition
%f = zeros(N,1);
x_space = linspace(x_left,x_right,N);

f_l = 1;
f = f_l*exp(-(x_space.*x_space));
f((x_space < -5) & (x_space > -7)) = 1;

f = f';
f_tmp = zeros(N,1);
f_tmp1 = zeros(N,1);
f_nm = zeros(N,1);

% plot initial condition
plot(x_space,f,'-b');
axis([x_left x_right 0 1.1*f_l]);
title('\bf{Initial Condition}');
drawnow


omega = (4*nu*nu+1)*(4-nu*nu)/5;


% time step
tic
for ts = 1:Num_ts
    
    if(mod(ts,500)==0)
       fprintf('Executing time step number %d.\n',ts);
    end
    % slower loop-driven update
    
    for x = 1:N
       x_p = x+1;
       
       % enforce periodicity at boundary
       if (x_p >N)
           x_p = 1;
       end
       
      f_nm(x) = 0.5*(f(x_p)+f(x))-(nu/3)*(f(x_p)-f(x));
         
      
    end
    
   
    
    for x = 1:N
       x_m = x-1;
       if(x_m<1)
           x_m = N;
       end
       f_tmp1(x)=f(x)-(2*nu/3)*(f_nm(x)-f_nm(x_m));
        
        
    end
    
   
    
    for x = 1:N
        x_p = x+1;
       % enforce periodicity at boundary
       if (x_p >N)
           x_p = 1;
       end
       
       x_2p = x+2;
       if(x_2p>N)
          x_2p = x_2p-N; 
       end
       
             
       x_m = x-1;
       % enforce periodicity at boundary
       if(x_m < 1)
           x_m = N;
       end
       
       x_2m = x-2;
       if(x_2m<1)
           x_2m=x_2m+N;
       end
       
       f_tmp(x)=f(x)-(nu/24)*(-2*f(x_2p)+7*f(x_p)-7*f(x_m)+2*f(x_2m))...
           -(3*nu/8)*(f_tmp1(x_p)-f_tmp1(x_m))...
           -(omega/24)*(f(x_2p)-4*f(x_p)+6*f(x)-4*f(x_m)+f(x_2m));
        
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