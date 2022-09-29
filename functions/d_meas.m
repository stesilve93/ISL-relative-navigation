function y = d_meas(x)

dr = x(10:12)-x(1:3);

y = dr'*dr;

end