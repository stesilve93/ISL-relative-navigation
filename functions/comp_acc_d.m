function [a_d_A, a_d_B] = comp_acc_d(x)

mi = 398600;

rrA = x(:,1:3);
rA = vecnorm(rrA,2,2);
n = length(rA);
a_mainA = zeros(n,3);
taccA = x(:,7:9);

rrB = x(:,10:12);
rB = vecnorm(rrB,2,2);
n = length(rB);
a_mainB = zeros(n,3);
taccB = x(:,16:18);

for i = 1:n
   
a_mainA(i,:) = -mi/rA(i)^3*rrA(i,:);
a_mainB(i,:) = -mi/rB(i)^3*rrB(i,:);

end

a_d_A = taccA - a_mainA;
a_d_B = taccB - a_mainB;

end