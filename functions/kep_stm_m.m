function stm = kep_stm_m(state, dt, mi)

n = length(state);

xA = state(1);
yA = state(2);
zA = state(3);

% dxA = state(4);
% dyA = state(5);
% dzA = state(6);
% 
% axA = state(7);
% ayA = state(8);
% azA = state(9);
% 
% xB = state(10);
% yB = state(11);
% zB = state(12);

% dxB = state(13);
% dyB = state(14);
% dzB = state(15);
% 
% axB = state(16);
% ayB = state(17);
% azB = state(18);

% Auxiliary variables
rAq = (xA^2 + yA^2 + zA^2)^(5/2);
% rBq = (xB^2 + yB^2 + zB^2)^(5/2);

J1 = [mi*(2*xA^2 - yA^2 - zA^2)/rAq,  3*mi*xA*yA/rAq,                 3*mi*xA*zA/rAq;
      3*mi*xA*yA/rAq,                 mi*(2*yA^2 - xA^2 - zA^2)/rAq,  3*mi*yA*zA/rAq;
      3*mi*xA*zA/rAq,                 3*mi*zA*yA/rAq,                 mi*(2*zA^2 - xA^2 - yA^2)/rAq];
  
% J2 = [mi*(2*xB^2 - yB^2 - zB^2)/rBq,  3*mi*xB*yB/rBq,                 3*mi*xB*zB/rBq;
%       3*mi*xB*yB/rBq,                 mi*(2*yB^2 - xB^2 - zB^2)/rBq,  3*mi*yB*zB/rBq;
%       3*mi*xB*zB/rBq,                 3*mi*zB*yB/rBq,                 mi*(2*zB^2 - xB^2 - yB^2)/rBq];

c1 = [zeros(3,3);
      J1];
  
c2 = [eye(3);
      zeros(3,3)];
  
% c3 = [zeros(12,3);
%       J2;
%       zeros(3,3)];
%   
% c4 = [zeros(9,6);
%       eye(6);
%       zeros(3,6)];

% Jacobian of the dynamics wrt state vector
% J = [c1, c2, c3, c4];
J = [c1, c2];


stm = eye(n) + J*dt;

% adding Gauss-Markov process terms
% m = eye(3)*exp(-dt/tau);
% 
% stm(7:9,7:9) = m;
% stm(16:18,16:18) = m;
    

end
