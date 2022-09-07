% Scenario parameters setup  
% Spacecraft Formation Flying in LEO

%% Environmental parameters
% Earth
data.env.r_e     = 6371;                          % [km]
data.env.mi      = 398600.44;                     % [km^3/s^2]
data.env.w_e     = 7.2921150e-5;                  % [rad/s] Earth's angular speed

% Sun
data.env.r_sun   = 149.5978707e6;                 % [km] 1 AU
data.env.T_sun   = 31556952;                      % [s] Earth's orbital period
data.env.n_sun   = 2*pi/data.env.T_sun;           % [rad/s] Earth's mean motion
data.env.eps     = deg2rad(23.45);                % [rad] ecliptic plane tilt
data.env.srp     = 1358/299792458;                % [Pa] solar radiation pressure (1AU)
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
data.satA.mass          = 11.4;
data.satA.CD            = 2.2;
data.satA.BC            = data.satA.CD * data.satA.area / data.satA.mass;

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
data.satB.mass          = 11.4;
data.satB.CD            = 2.2;
data.satB.BC            = data.satB.CD * data.satB.area / data.satB.mass;

%% Sensors
% Distance measurement
data.sens.sig   = 1;                              % [km] Sensor standard deviation
data.sens.st    = 10;                             % [s] Sample time

%% Simulation
data.sim.n_orb  = 2;
data.sim.time   = data.sim.n_orb * data.satA.orbit.T;
data.sim.date   = [2020, 1, 1, 12, 0, 0];
data.sim.jdate  = juliandate(2020, 1, 1, 12, 0, 0);

%% WIP

% Spacecraft properties - 6U CubeSat
m_sc    = 11.4;                       % [kg] s/c mass (max)
m_sa    = 0.3;                        % [kg] deployable solar arrays mass
dim_sc  = [366, 100, 226.3]*1e-3;     % [m] s/c dimensions (main body)
dim_sa  = [366, 225]*1e-3;            % [m] solar arrays dimensions
mag_res = [1; 1; 1]*1e-3;             % [A*m^2] residual magnetic moment

% Spacecraft surfaces in [m^2]
A_sc   = [dim_sc(2)*dim_sc(3); dim_sc(2)*dim_sc(1); dim_sc(1)*dim_sc(3)];
A_sa   = dim_sa(1)*dim_sa(2);

A_i    = [A_sc(1); A_sc(2); A_sc(3); A_sc(1); A_sc(2); A_sc(3); ...
          A_sa; A_sa; A_sa; A_sa]';   % [m^2] surfaces vector

% Normal vector to each surface (Body F.O.R.)
n_sa   = [0; 0; 1];
n_i    = [eye(3), -eye(3), n_sa, -n_sa, n_sa, -n_sa];

% Center of pressure coordinates for each surface in [m] Body F.O.R.
% Main body
CG_1  = [dim_sc(1)/2; 0; 0];
CG_2  = [0; dim_sc(3)/2; 0];
CG_3  = [0; 0; dim_sc(2)/2];
CG_4  = -CG_1;
CG_5  = -CG_2;
CG_6  = -CG_3;
% Solar arrays
CG_7  = [0; dim_sc(3)/2 + dim_sa(2)/2; dim_sc(2)/2];
CG_8  = CG_7;
CG_9  = CG_7.*[1; -1; 1];
CG_10 = CG_9;

r_i   = [CG_1, CG_2, CG_3, CG_4, CG_5, CG_6, CG_7, CG_8, CG_9, CG_10];

% Optical properties for each surface
rhoS  = [ones(1,6)*0.5, ones(1,4)*0.8];	 % [-] reflectivity
rhoD  = ones(1,10)*0.1;                  % [-] diffusivity

