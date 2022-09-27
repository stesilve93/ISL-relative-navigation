function rot = ECI2LVLH(r,v)
h = cross(r,v);
rh = normalize(r,'norm');
hh = normalize(h,'norm');



x_hat = rh;
z_hat = hh;
y_hat = cross(z_hat,x_hat);

rot = [x_hat;
       y_hat;
       z_hat];
end