function y = gpsB_meas(x)

HB = [zeros(6,9), eye(6), zeros(6,3)];

y = HB*x;

end