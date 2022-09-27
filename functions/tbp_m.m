function dy = tbp_m(t, x, acc)

    mi = 398600;

    posA = x(1:3);
    velA = x(4:6);
    
    rA = norm(posA);
    
% Set the derivatives of the states
    dy = [velA;
          -mi/rA^3*posA + acc(1,:)'];

end