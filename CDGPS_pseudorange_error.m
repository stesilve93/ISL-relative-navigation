clear, close all

r1 = 6771;
r2 = 26571;

th1 = 1.2693166;
th2 = -90:.1:90;


AB = r1*sqrt(2*(1-cosd(th1)));

GA = sqrt(r1^2 + r2^2 - 2*r1*r2*cosd(th2-th1));
GB = sqrt(r1^2 + r2^2 - 2*r1*r2*cosd(th2));

err = AB - abs(GB - GA);

figure
plot(th2, err)