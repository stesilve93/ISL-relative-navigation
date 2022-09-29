function [x_sol,t_sol,fevals,cpu_time] = RK4(f,x0,h,tspan)
% This function implements Runge-Kutta 4-th order method (RK4)
% 
% INPUT
% f         = function handle of the ODE's RHS
% x0        = initial condition
% h         = method's stepsize
% tspan     = time interval of integration
% 
% OUTPUT
% x_sol     = solutions history vector in tspan
% t_sol     = time vector
% fevals    = number of RHS evaluations
% cpu_time  = computational time

    tic

    % Avoid problems if the interval is not a multiple of h
    nstep = ceil((max(tspan)-min(tspan))/h);
    t_sol = min(tspan):h:max(tspan);

    % Make sure to catch all the dynamics, until the end of a generic tspan
    if t_sol(end) < max(tspan)
        t_sol = [t_sol, t_sol(end)+h];
    end
    
    n_st = length(x0);
    
    x_c = [x0,zeros(n_st,nstep)];
    k = 0;

    for i = 1:nstep

        k1 = f(x_c(:,i),        t_sol(i));             % 1st fcn eval
        k2 = f(x_c(:,i)+h/2*k1, t_sol(i)+h/2);         % 2nd fcn eval
        k3 = f(x_c(:,i)+h/2*k2, t_sol(i)+h/2);         % 3rd fcn eval
        k4 = f(x_c(:,i)+h*k3,   t_sol(i)+h);           % 4th fcn eval

        x_c(:,i+1) = x_c(:,i)+h/6*(k1+2*k2+2*k3+k4);     % Corrector

        k = k+4;     % 4 fcn evals per step

    end

    x_sol = x_c;
    fevals = k;
    cpu_time = toc;
end
