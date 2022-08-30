function [stm, Z] = cw_stm(om, dt)

    phi_rri = [4-3*cos(om*dt),         0;
               6*(sin(om*dt)-om*dt), 1];
           
    phi_rvi = 1/om * [sin(om*dt),         2*(1-cos(om*dt));
                      2*(cos(om*dt)-1),   4*sin(om*dt)-3*om*dt];
                  
    phi_vri = [3*om*sin(om*dt),      0;
               6*om*(cos(om*dt)-1),  0];
           
    phi_vvi = [cos(om*dt),    2*sin(om*dt);
               -2*sin(om*dt), 4*cos(om*dt)-3];
           
    phi_zzi = [cos(om*dt),      sin(om*dt)/om;
               -om*sin(om*dt),  cos(om*dt)];
           
    PHI_rri = [phi_rri,     zeros(2,1);
               zeros(1,2),  phi_zzi(1,1)];
    
    PHI_rvi = [phi_rvi,     zeros(2,1);
               zeros(1,2),  phi_zzi(1,2)];    
           
    PHI_vri = [phi_vri,     zeros(2,1);
               zeros(1,2),  phi_zzi(2,1)];
           
    PHI_vvi = [phi_vvi,     zeros(2,1);
               zeros(1,2),  phi_zzi(2,2)];
           
    stm     = [PHI_rri, PHI_rvi;
               PHI_vri, PHI_vvi];
           
    Z       = [PHI_rri'; PHI_rvi'] * [PHI_rri, PHI_rvi];

end