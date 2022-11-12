function [xrel, Prelpos_new_rot] = relstate(x, P)
mi = 398600;

rA = x(1:3);
vA = x(4:6);
aA = x(7:9) - mi/norm(rA)^3*rA;

rB = x(10:12);
vB = x(13:15);
aB = x(16:18) - mi/norm(rB)^3*rB;

rot = ECI2LVLH(rA,vA);
h = cross(rA,vA);
r2 = norm(rA)^2;
om = h/r2;
om_dot = -2*dot(rA,vA)/r2*om;

r_rel = rB-rA;
cr_o_r = cross(om,r_rel);
v_rel = vB-vA - cr_o_r;
a_rel = aB-aA - cross(om_dot,r_rel) - cross(om,cr_o_r) - 2*cross(om,v_rel);

r_rel = rot*r_rel';
v_rel = rot*v_rel';
a_rel = rot*a_rel';

xrel = [r_rel', v_rel', a_rel'];

% Covariance of relative state
P_A  = P(  1:3,   1:3);
P_B  = P(10:12, 10:12);
P_AB = P(  1:3, 10:12);
P_BA = P(10:12,   1:3);

P_rel = [P_A, P_AB;
         P_BA, P_B];

reltrn = [-eye(3), eye(3)];

Prelpos_new = reltrn*P_rel*reltrn';
Prelpos_rot = rot*Prelpos_new*rot';

Prelpos_new_rot = sqrt(diag(Prelpos_rot));
end