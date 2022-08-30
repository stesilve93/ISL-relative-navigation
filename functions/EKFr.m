function [x,y_e,P] = EKFr(t,y_m,x0,P0,Q,R,data)

% [x,y_e,P] = EKFr(t,y_m,x0,P0,Q,R,data)
%
% This function implements an EKF to perform relative navigation based on
% range measurements only.
% Ref: J. Christian, "Relative Navigation Using Only Intersatellite Range
% Measurements"
% 
% INPUT
% t     time vector
% y_m   noisy measurements vector
% x0    initial state estimate
% P0    initial state covariance estimate
% Q     process noise covariance matrix
% R     measurement noise covariance matrix
% data  problem's data
% 
% OUTPUT
% x     estimated state
% y_e   estimated measurement (filtered)
% P     estimated state covariance

% Settings
n = length(t);
om = data.satA.orbit.omv;   % [rad/s]

% Filter initialization
x_p = x0;
P_p = P0;

% Output preallocation
x = zeros(n,6);
y_e = zeros(n,1);
P = zeros(6,6,n);
x(1,:) = x0;
y_e(1)   = x0(1:3)'*x0(1:3);
P(:,:,1) = P0;

for i = 2:n
    % Linearized measurement model
    H = 2*x_p(1:3)'*[eye(3) zeros(3)];
    
    % State and covariance propagation
    dt = t(i)-t(i-1);
    [stm,~] = cw_stm(om,dt);
    x_m = stm*x_p;
    P_m = stm*P_p*stm'+Q;
    
    % Gain update
    K = P_m*H'/(H*P_m*H'+R);
    
    % Measurement estimation
    y_est = x_m(1:3)'*x_m(1:3);
    
    % Correction
    x_p = x_m + K*(y_m(i)^2 - y_est);
    W = eye(6)-K*H;
    P_p = W*P_m*W'+K*R*K';
    
    % Output storage
    x(i,:) = x_p;
    y_e(i)   = x_p(1:3)'*x_p(1:3);
    P(:,:,i) = P_p;
end
    
    % Measurements are squared in the filter
    y_e = sqrt(y_e);

end