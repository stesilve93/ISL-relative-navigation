function [x0,P0] = init_state(data)


%% Data loading
load('sim_data.mat','simout');

%% Pre-processing
% time signatures of measurements
tA = simout.GPS_meas.GPS_rec_A.tA.Data;
tB = simout.GPS_meas.GPS_rec_B.tB.Data;

% number of measurements
nA = length(tA);
nB = length(tB);

% position measurements
posA = simout.GPS_meas.GPS_rec_A.posA.Data;
posB = simout.GPS_meas.GPS_rec_B.posB.Data;
posA = reshape(posA,3,nA)';
posB = reshape(posB,3,nB)';

% velocity measurements
velA = simout.GPS_meas.GPS_rec_A.velA.Data;
velB = simout.GPS_meas.GPS_rec_B.velB.Data;
velA = reshape(velA,3,nA)';
velB = reshape(velB,3,nB)';

% GPS measurement vector
measA = [posA, velA];
measB = [posB, velB];


%% Initialization
% filter time vector
stF = 1;

% filter settings
% Measurement noise
sigp = eye(3)*data.sens.gps.sig_p;
sigv = eye(3)*data.sens.gps.sig_v;
R = [sigp,     zeros(3);
     zeros(3), sigv];
R = R*stF;
R = R*1e2;

% Process noise
tau = 600;
m_gm = exp(-stF/tau);
siga = 1e-5;
qm = siga^2*(1-m_gm^2)*eye(3);

% initial state and covariance
x0 = [measA(1,:), [0, 0, 0], ...
      measB(1,:), [0, 0, 0]];
  
x0 = x0';

P0 = [R, zeros(6,12);
      zeros(3,6), qm, zeros(3,9);
      zeros(6,9), R, zeros(6,3);
      zeros(3,15), qm];

end
