% Scenario parameters setup  
% Spacecraft Formation Flying in LEO

%% Environmental parameters
% Earth
data.env.r_e     = 6371;                          % [km]
data.env.mi      = 398600.44;                     % [km^3/s^2]
data.env.w_e     = 7.2921150e-5;                  % [rad/s] Earth's angular speed

% Ionosphere
[data.env.iono.latID, data.env.iono.lonID, data.env.iono.hID, data.env.iono.IRItab] = iono4data(0);

% Sun
data.env.r_sun   = 149.5978707e6;                 % [km] 1 AU
data.env.T_sun   = 31556952;                      % [s] Earth's orbital period
data.env.n_sun   = 2*pi/data.env.T_sun;           % [rad/s] Earth's mean motion
data.env.eps     = deg2rad(23.45);                % [rad] ecliptic plane tilt
data.env.c       = 299792458;
data.env.srp     = 1358/data.env.c;               % [Pa] solar radiation pressure (1AU)
data.env.th0_sun = deg2rad(-79.75);               % [deg] Sun initial true anomaly
data.env.mi_sun  = 132.712e9;                     % [km^3/s^2]

% Moon
data.env.mi_moon = 4903;                          % [km^3/s^2]

%% Formation Flying
data.ff.b        = 150;                           % [km] Nominal baseline

%% Spacecraft A - Target
% Orbital parameters
data.satA.orbit.h       = 400;                    % [km]
data.satA.orbit.a       = data.env.r_e...
    + data.satA.orbit.h;                          % [km]
data.satA.orbit.e       = 0;                      % [-]
data.satA.orbit.i       = 97.0494;                % [deg]
data.satA.orbit.RAAN    = 180;                    % [deg]
% data.satA.orbit.RAAN = 337.5;                   % [deg]
data.satA.orbit.om      = 0;                      % [deg]
data.satA.orbit.th      = 0;                      % [deg]
data.satA.orbit.T       = 2*pi*sqrt(...
    data.satA.orbit.a^3/data.env.mi);             % [s]
data.satA.orbit.omv     = 2*pi/data.satA.orbit.T;

% Spacecraft properties
data.satA.area          = 0.1^2*6*2 + 0.1^2*sqrt(2)*6;
data.satA.mass          = 24;
data.satA.CD            = 2.2;
data.satA.BC            = data.satA.CD * data.satA.area / data.satA.mass;
data.satA.sig_mp        = 0.005; 
data.satA.lnd           = 1e-8;                 

%% Spacecraft B - Chaser
% Orbital parameters
data.satB.orbit.h       = 400;                    % [km]
data.satB.orbit.a       = data.env.r_e...
    + data.satB.orbit.h;                          % [km]
data.satB.orbit.e       = 0;                      % [-]
data.satB.orbit.i       = 97.0494;                % [deg]
data.satB.orbit.RAAN    = 180;                    % [deg]
data.satB.orbit.om      = 0;                      % [deg]

%%% WIP: initial contitions setup based on specified baseline %%%
% formula for same, circular orbit of the two s/c
th = 2 * asind(data.ff.b/2/data.satA.orbit.a);
%%% END %%%

data.satB.orbit.th      = th;                     % [deg]
clear th

data.satB.orbit.T       = 2*pi*sqrt(...
    data.satB.orbit.a^3/data.env.mi);             % [s]

% Spacecraft properties
data.satB.area          = 0.1^2*6*2 + 0.1^2*sqrt(2)*6;
data.satB.mass          = 24;
data.satB.CD            = 2.2;
data.satB.BC            = data.satB.CD * data.satB.area / data.satB.mass;
data.satB.sig_mp        = 0.005; 
data.satB.lnd           = 1e-8;      

%% Sensors
% ISL 
% - Distance measurement
data.sens.isr.sig   = 0.3;                              % [km] Sensor standard deviation
data.sens.isr.st    = 10;                             % [s] Sample time
 
% % - Receiver clock model: ovenized crystal oscillator
% data.sens.isr.sf = 8e-20/2;
% data.sens.isr.sg = 2*pi^2*4e-23;

% % - Receiver clock model: Rubidium
% data.sens.isr.sf = 2e-20/2;
% data.sens.isr.sg = 2*pi^2*4e-29;

% - Receiver clock model: CSAC
data.sens.isr.sf = 7.2e-21/2;
data.sens.isr.sg = 2*pi^2*2.7e-27;

% - Link design
data.sens.isr.f    = 420e6;                             % [Hz]
data.sens.isr.lmb  = data.env.c/data.sens.isr.f;        % [m]
data.satA.sig_mp   = 0.005/data.sens.isr.lmb; % Spurious [cycles]
data.sens.isr.SNR  = 10^(25/10); % [dB->Hz]
data.sens.isr.bPLL = 15; % [Hz]
data.sens.isr.bDLL = 1; % [Hz]
data.sens.isr.sig_thp = 1/(2*pi)*sqrt(data.sens.isr.bPLL/(2*data.sens.isr.SNR)); % [cycles]
data.sens.isr.tc = 1/1023e6;
data.sens.isr.sig_thc = data.env.c * data.sens.isr.tc*...
    sqrt(data.sens.isr.bDLL * data.sens.isr.tc/2/(2*data.sens.isr.tc * data.sens.isr.SNR) *...
    (1+1/2/data.sens.isr.tc/data.sens.isr.SNR)); % [m]
% - Antenna phase center deviation
data.sens.isr.sig_pc = 0.01/data.sens.isr.lmb;

% - Phase bias
data.satA.ph0 = 130.58;
data.satB.ph0 = 25.42;
% data.satA.ph0 = 0;
% data.satB.ph0 = 0;

% GPS receiver
% data.sens.gps.sig_p   = 10e-3;
data.sens.gps.sig_p   = 1;
data.sens.gps.sig_v   = 0.5e-3;
data.sens.gps.sig_t   = 20e-9;  
data.sens.gps.stA   = 60;
data.sens.gps.stB   = 62;

% Delay (Time of Flight satB->satA)
data.sens.del = data.ff.b*1e3/data.env.c;
data.sens.del = data.sens.del * 2;

%% Filter
data.filt.st = 60;

%% Simulation
data.sim.n_orb  = 1;
data.sim.time   = data.sim.n_orb * data.satA.orbit.T;
data.sim.date   = [2020, 1, 1, 12, 0, 0];
data.sim.jdate  = juliandate(2020, 1, 1, 12, 0, 0);

%% SRP

% Spacecraft properties - 6U CubeSat
data.satAB.dim_sc = [360, 240, 230]*1e-3;     % [m] s/c dimensions (main body)
data.satAB.dim_sa = [360, 600]*1e-3;          % [m] solar arrays dimensions

% Spacecraft surfaces in [m^2]
data.satAB.A_sc   = [data.satAB.dim_sc(2)*data.satAB.dim_sc(3);
                     data.satAB.dim_sc(2)*data.satAB.dim_sc(1);
                     data.satAB.dim_sc(1)*data.satAB.dim_sc(3)];
data.satAB.A_sa   =  data.satAB.dim_sa(1)*data.satAB.dim_sa(2)*ones(4,1);

data.satAB.A_i    = [data.satAB.A_sc(1:3);
                     data.satAB.A_sc(1:3);
                     data.satAB.A_sa]';          % [m^2] surfaces vector

% Normal vector to each surface (Body F.O.R.)
data.satAB.n_sa   = [0; 0; 1];
data.satAB.n_i    = [eye(3), -eye(3), data.satAB.n_sa, -data.satAB.n_sa, data.satAB.n_sa, -data.satAB.n_sa];

% Center of pressure coordinates for each surface in [m] Body F.O.R.
data.satAB.r_i = zeros(3,10);
% Main body
data.satAB.r_i(:,1)   = [data.satAB.dim_sc(1)/2; 0; 0];
data.satAB.r_i(:,2)   = [0; data.satAB.dim_sc(3)/2; 0];
data.satAB.r_i(:,3)   = [0; 0; data.satAB.dim_sc(2)/2];
data.satAB.r_i(:,4:6) = -data.satAB.r_i(:,1:3);
% Solar arrays
data.satAB.r_i(:,7)   = [0; data.satAB.dim_sc(3)/2 + data.satAB.dim_sa(2)/2; data.satAB.dim_sc(2)/2];
data.satAB.r_i(:,8)   = data.satAB.r_i(:,7);
data.satAB.r_i(:,9)   = data.satAB.r_i(:,7).*[1; -1; 1];
data.satAB.r_i(:,10)  = data.satAB.r_i(:,9);

% Optical properties for each surface
data.satAB.rhoS  = [ones(1,6)*0.5, ones(1,4)*0.8];	 % [-] reflectivity
data.satAB.rhoD  = ones(1,10)*0.1;                   % [-] diffusivity

%% EKF
data.ekf.st = 1;

% Process noise
tau = 600;
m_gm = exp(-data.ekf.st/tau);
data.ekf.siga = 1e-5;
qm = data.ekf.siga^2*(1-m_gm^2)*eye(3);

data.ekf.Q = [zeros(6,18);
              zeros(3,6), qm, zeros(3,9);
              zeros(6,18);
              zeros(3,15), qm];
clear('tau','m_gm','qm')

% Measurement noise
data.ekf.sigp = eye(3)*data.sens.gps.sig_p^2;
data.ekf.sigv = eye(3)*data.sens.gps.sig_v^2;

data.ekf.R = [data.ekf.sigp,     zeros(3);
              zeros(3), data.ekf.sigv];
data.ekf.R = data.ekf.R*data.ekf.st;
data.ekf.R = data.ekf.R * 1e-1;

data.ekf.Rd = data.sens.isr.sig^2;
data.ekf.Rd = data.ekf.Rd * 1e4;

[data.ekf.x0, data.ekf.P0] = init_state(data);