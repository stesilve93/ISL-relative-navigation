function y = gpsA_meas(x)

HA = [eye(6), zeros(6,12)];

y = HA*x;

end