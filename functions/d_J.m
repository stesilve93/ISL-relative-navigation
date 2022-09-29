function Jd = d_J(x)

Hd = [-eye(3), zeros(3,6), eye(3), zeros(3,6)];
dr = x(10:12)-x(1:3);

Jd = 2*dr'*Hd;

end