function [x,P] = irodr(data, simout)

% [x,P] = irodr(data, simout)
% 
% This function performs Initial Relative Orbit Determination (IROD) by
% processing range only measurements.
% Ref: J. Christian, "Relative Navigation Using Only Intersatellite Range
% Measurements"
% 
% INPUT
% data
% simout
%
% OUTPUT
% x
% P

n_orb = data.sim.n_orb;

if n_orb < 1
    warning('The simulation has been run on less than one orbit.\nResults may be inaccurate.')
end

om = 2*pi/data.satA.orbit.T;
t = simout.tout(simout.tout <= data.satA.orbit.T);
m = length(t);

Z = zeros(6,6,m);

for i = 1:m 
    [~, Z(:,:,i)] = cw_stm(om, t(i)); 
end

% rA = simout.r_A;
% rB = simout.r_B;
% rAB = rB - rA;
% r_rel = vecnorm(rAB, 2, 2);
% r_rel = r_rel(1:m);

r_rel = simout.r_meas(1:m);
y = r_rel.^2;

% Covariance estimate
G = diag(r_rel);
R = data.sens.isr.sig^2 * eye(m);
Ry = 4*G*R*G';
Rinv = inv(Ry);

% Init. state
x0 = [1 1 1 1 1 1]';
k = 0;
maxiter = 100;
A = zeros(m,6);
upd = 1;
tol = 0.001;

while k < maxiter && upd > tol
    for i = 1:m
        A(i,:) = x0'*Z(:,:,i);
    end
    
    z = (A'*Rinv*A)\A'*Rinv*y;
    
    x = 0.5*(x0+z);
    
    upd = abs(mean(x-x0));
    x0 = x;
    k = k + 1;
end

P = 1/4*inv(A'*Rinv*A);

% Mirror solutions computation
% THE FOLLOWING SOLUTION IS UNIQUE ONLY FOR IN-PLANE NON-PERIODIC RELATIVE ORBITS

x = [x, -x]; 

end