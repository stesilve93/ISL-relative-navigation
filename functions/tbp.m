function dy = tbp(~, x, acc)

    mi = 398600;

    posA = x(1:3);
    posB = x(7:9);
    
    velA = x(4:6);
    velB = x(10:12);
    
    rA = norm(posA);
    rB = norm(posB);
    
% Set the derivatives of the states
    dy = [velA;
          -mi/rA^3*posA + acc(1,:)';
          velB;
          -mi/rB^3*posB + acc(2,:)'];

end