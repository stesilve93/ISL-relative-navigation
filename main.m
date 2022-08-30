% ISL for relative navigation
% MSc thesis
% Author: Francesco De Cecio
clear, close all
addpath('./functions')

%% Simulation setup
% Load scenario
data_setup

%% Simulation
simout = sim('ff_model');

%% Output processing
t = simout.tout;
t_orb = t/data.satA.orbit.T;

% Absolute position (ECI frame)
rA = simout.r_A;
rB = simout.r_B;

% Relative position (LVLH frame)
r_rel = simout.relpos;
BL = vecnorm(r_rel, 2, 2);

%% CW (STM)
n = length(t);
om = data.satA.orbit.omv;
x_cw = zeros(n,6);
x_cw(1,:) = [simout.relpos(1,:), simout.relvel(1,:)];

for i = 1:n
    [stm, ~] = cw_stm(om, t(i));
    x_cw(i,:) = stm * x_cw(1,:)';
end

% diff_CW = x_cw(:,1:3) - r_rel_CW;

%% State estimation
% Initial relative orbit determination - range only
[irod_state, irod_cov] = irodr(data,simout);
irod_err = [simout.relpos(1,:), simout.relvel(1,:)]' - irod_state;

% Range only EKF
x0 = irod_state(:,2);
P0 = irod_cov;
t  = simout.tout;
y_m  = simout.r_meas;

% Noise covariance matrices
% Q = zeros(6);   % Process noise (6 states)
G = [1 10 1e-4 1e-4 1e-4 1e-4];
Q = G'*G*1e-5;   % Process noise (6 states)
R = 1e5;          % Measurement noise (scalar)

[x,y_e,P] = EKFr(t,y_m,x0,P0,Q,R,data);

%% Data print

fprintf('\nINITIAL RELATIVE ORBIT DETERMINATION\n\n')
fprintf('Estimate of the initial state and its covariance:\n')
for i = 1:6
    fprintf('%+.2f\t+-\t%.2f\t(Error:\t%+.2f)\n',irod_state(i,1),irod_cov(i,i),irod_err(i,1))
end
fprintf('\nOr its mirror solution:\n')
for i = 1:6
    fprintf('%+.2f\t+-\t%.2f\t(Error:\t%+.2f)\n',irod_state(i,2),irod_cov(i,i),irod_err(i,2))
end

%% Plots
close all

figure(1)
    plot3(rA(:,1), rA(:,2), rA(:,3))
    hold on
    scatter3(rA(1,1), rA(1,2), rA(1,3), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
    sun_pnt = simout.r_sun(1,:)*2500;
    quiver3(rA(1,1), rA(1,2), rA(1,3), sun_pnt(1), sun_pnt(2), sun_pnt(3))
    text(rA(1,1)+sun_pnt(1), rA(1,2)+sun_pnt(2), rA(1,3)+sun_pnt(3), 'To Sun')
    % earth3D
    axset(7000)
    plot3(rB(:,1), rB(:,2), rB(:,3))
    scatter3(rB(1,1), rB(1,2), rB(1,3), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
    axis equal
    grid on

    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')

    
figure(2)
    plot3(r_rel(:,1),r_rel(:,2),r_rel(:,3),'DisplayName','Relative trajectory')
    hold on
    scatter3(0, 0, 0, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r','DisplayName','Target')
    scatter3(r_rel(1,1), r_rel(1,2), r_rel(1,3), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b','DisplayName','Chaser (t_0)')
    grid on

    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')
    
    legend('Location','northeast')
    
    
figure(3)
subplot(2,3,1)
    plot(t_orb,r_rel(:,1))
    hold on
    plot(t_orb,x(:,1))
    grid on
    xlabel('Time (orbits)')
    ylabel('x [km]')
subplot(2,3,4)
    plot(t_orb,(r_rel(:,1)-x(:,1))*1e3)
    grid on
    xlabel('Time (orbits)')
    ylabel('Error (m)')
subplot(2,3,2)
    plot(t_orb,r_rel(:,2))
    hold on
    plot(t_orb,x(:,2))
    grid on
    xlabel('Time (orbits)')
    ylabel('y [km]')
subplot(2,3,5)
    plot(t_orb,(r_rel(:,2)-x(:,2))*1e3)
    grid on
    xlabel('Time (orbits)')
    ylabel('Error (m)')
subplot(2,3,3)
    plot(t_orb,r_rel(:,3))
    hold on
    plot(t_orb,x(:,3))
    grid on
    xlabel('Time (orbits)')
    ylabel('z [km]')
subplot(2,3,6)
    plot(t_orb,(r_rel(:,3)-x(:,3)))
    grid on
    xlabel('Time (orbits)')
    ylabel('Error (km)')
    
% figure(6)
% subplot(2,1,1)
%     plot(t_orb,r_rel(:,1))
%     hold on
%     plot(t_orb,x(:,1))
%     grid on
%     xlabel('Time (orbits)')
%     ylabel('x [km]')
% subplot(2,1,2)
%     plot(t_orb,(r_rel(:,1)-x(:,1))*1e3)
%     grid on
%     xlabel('Time (orbits)')
%     ylabel('Error (m)')
%     
% figure(7)
% subplot(2,1,1)
%     plot(t_orb,r_rel(:,2))
%     hold on
%     plot(t_orb,x(:,2))
%     grid on
%     xlabel('Time (orbits)')
%     ylabel('y [km]')
% subplot(2,1,2)
%     plot(t_orb,(r_rel(:,2)-x(:,2))*1e3)
%     grid on
%     xlabel('Time (orbits)')
%     ylabel('Error (m)')
%     
% figure(8)
% subplot(2,1,1)
%     plot(t_orb,r_rel(:,3))
%     hold on
%     plot(t_orb,x(:,3))
%     grid on
%     xlabel('Time (orbits)')
%     ylabel('z [km]')
% subplot(2,1,2)
%     plot(t_orb,(r_rel(:,3)-x(:,3)))
%     grid on
%     xlabel('Time (orbits)')
%     ylabel('Error (km)')
    
    
figure(4)
    hold on
    plot(t_orb, y_m)
    plot(t_orb, y_e)
    plot(t_orb, BL)
    
    xlabel('Time (orbits)')
    ylabel('Baseline (km)')
    
    legend('Measured','EKF','True')
    
    grid on
    
