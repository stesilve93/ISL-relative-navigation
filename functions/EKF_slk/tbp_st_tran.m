function x_m = tbp_st_tran(x_prev)

% State time update (prediction)
dt = 1;
t_span = [0, dt];
% t_span = [0, dt/2, dt];

x_int = [x_prev(1:6);
    x_prev(10:15)];

a_int = [x_prev(7:9)';
    x_prev(16:18)'];

step = 0.1;
odefun = @(x,t) tbp(t,x,a_int);
[x_sol,~,~,~] = RK4(odefun,x_int,step,t_span);

x_m = x_sol(:,end)';

% Empirical acceleration time update
tau = 600;
m_gm = exp(-dt/tau);
a_m = m_gm*a_int;

% State reordering
x_m = [x_m(1:6), a_m(1,:), x_m(7:12), a_m(2,:)]';

end