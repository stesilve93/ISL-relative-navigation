clear

%% Data loading
addpath('functions')
load('sim_data.mat')
data_setup

%% Mode selector
% Logical switch
% m_upd     Use measurements for state update
% m_sync    Propagate measurements up to filter run time
% m_d       Include ISR measure
m_upd  = 1;
m_sync = 1;
m_d    = 1;

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

% True state extraction 
ttime = simout.tout;
ntrue = length(ttime);

tposA = simout.r_A;
tvelA = simout.v_A;

tposB = simout.r_B;
tvelB = simout.v_B;

taccA  = simout.a_A;
taccB  = simout.a_B;

% Absolute state assembly
tstate = [tposA, tvelA, taccA, tposB, tvelB, taccB];
% Relative state assembly
trelstate = [simout.relpos, simout.relvel, simout.relacc];

% Distance measurements
td = ttime;
measd = simout.r_meas.^2;

%% Initialization
% filter time vector
stF = 1;
tF = 0:stF:simout.tout(end);
tF = tF';

nF = length(tF);

% filter settings
% Measurement noise
sigp = eye(3)*data.sens.gps.sig_p;
sigv = eye(3)*data.sens.gps.sig_v;
R = [sigp,     zeros(3);
     zeros(3), sigv];
R = R*stF;
R = R*1e2;

Rd = data.sens.isr.sig;
Rd = Rd*1e6;

% Process noise
tau = 600;
m_gm = exp(-stF/tau);
siga = 1e-5;
qm = siga^2*(1-m_gm^2)*eye(3);
Q = [zeros(6,18);
     zeros(3,6), qm, zeros(3,9);
     zeros(6,18);
     zeros(3,15), qm];
 
HA = [eye(6), zeros(6,12)];
     
HB = [zeros(6,9), eye(6), zeros(6,3)];

Hd = [-eye(3), zeros(3,6), eye(3), zeros(3,6)];

% initial state and covariance
% x0 = [simout.r_A(1,:), simout.v_A(1,:), [0, 0, 0], ...
%     simout.r_B(1,:), simout.v_B(1,:), [0, 0, 0]];

x0 = [measA(1,:), [0, 0, 0], ...
      measB(1,:), [0, 0, 0]];

% P0 = [R, zeros(6,12);
%        zeros(3,18);
%        zeros(6,9), R, zeros(6,3);
%        zeros(3,18)];

P0 = [R, zeros(6,12);
      zeros(3,6), qm, zeros(3,9);
      zeros(6,9), R, zeros(6,3);
      zeros(3,15), qm];
   
% filter initialization
x_prev = x0;
P_prev = P0;

% output preallocation
x = zeros(nF,18);
P = zeros(18,18,nF);
x(1,:) = x0;
P(:,:,1) = P0;
opt = odeset('RelTol',1e-7,'AbsTol',1e-8);

% Live-time measurement
toll = 1;
[log_measA, ind_measA] = ismembertol(tF, tA, toll, 'DataScale', 1);
% [log_measB, ind_measB] = ismembertol(tF, tB, toll, 'DataScale', 1);
ind_measA_prev = 0;
ind_measB_prev = 0;

toll_d = 1e-3;
[log_measd, ind_measd] = ismembertol(tF, td, toll_d, 'DataScale', 1);
ind_measd_prev = 0;

dt_tol = 1e-4;

countA = 0;
countB = 0;
countd = 0;

%% Delay
c = 3e8;
BL = 150e3;
del = BL/c;
del = 2*del;

% Filter time at signal from B reception
tB_rec = tB + del;
% toll = 2e-3;
[log_measB, ind_measB] = ismembertol(tF, tB_rec, toll, 'DataScale', 1);

%% Extended kalman filter
% RK4_err = zeros(1,nF-1);

for i = 2:nF
    
    % State time update (prediction)
    dt = tF(i)-tF(i-1);
    
    t_span = [tF(i-1), tF(i)];
    
    x_int = [x_prev(1:3)';
             x_prev(4:6)';
             x_prev(10:12)';
             x_prev(13:15)'];
         
    a_int = [x_prev(7:9);
             x_prev(16:18)];
         
    [~, x_m] = ode45(@tbp, t_span, x_int, opt, a_int);
    x_m = x_m(end,:);
    
%     
%     %%%%%%%%%%%%%%%%% RK4 TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     step = 0.1;
%     odefun = @(x,t) tbp(t,x,a_int);
%     [x_sol,t_sol,fevals,cpu_time] = RK4(odefun,x_int,step,t_span);
%     RK4_sing = x_sol(:,end) - x_m';
%     RK4_err(i-1) = sum(abs(RK4_sing));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    
    % Empirical acceleration time update
    a_m = m_gm*a_int;
    
    % State reordering
    x_m = [x_m(1:6), a_m(1,:), x_m(7:12), a_m(2,:)];
    
    % Covariance time update
    stm = kep_stm(x_prev, dt, data.env.mi, tau);
    P_m = stm*P_prev*stm' + Q;
    
    % test stm performace
    absErr = x_m' - stm*x_prev';
    relErr = absErr./x_m'*100;
    
    % Update estimates in case of no measurements
    x_p = x_m;
    P_p = P_m;
    
    if log_measA(i) && ind_measA(i) ~= ind_measA_prev && tF(i)>=tA(ind_measA(i)) && m_upd
        % Measurement synchronization
        y_act = measA(ind_measA(i),:)';
        R_act = R;
        dt_m = tF(i)-tA(ind_measA(i));
        
        if dt_m > dt_tol && m_sync
            t_span = [tA(ind_measA(i)), tF(i)];
            
            y_int = y_act;
            
            a_int = x_m(7:9);
            
            [~, y_act] = ode45(@tbp_m, t_span, y_int, opt, a_int);
            y_act = y_act(end,:)';
            
            % Measurement noise covariance time update
            stm_m = kep_stm_m(y_int, dt_m, data.env.mi);
            R_act = stm_m*R*stm_m';
        end
        
        % Gain update
        K = P_m*HA'/(HA*P_m*HA' + R_act);
        
        % Correction
        y_m = y_act - HA*x_m';
        x_p = x_m' + K*y_m;
        x_p = x_p';
        W = eye(18)-K*HA;
        P_p = W*P_m*W'+K*R_act*K';
        
        % Measurement marked as used
        ind_measA_prev = ind_measA(i);
        
        countA = countA +1;
        
%         tF(i)
%         tA(ind_measA(i))
    end
    
    if log_measB(i) && ind_measB(i) ~= ind_measB_prev && tF(i)>=tB_rec(ind_measB(i)) && m_upd
        
        % Measurement synchronization
        y_act = measB(ind_measB(i),:)';
        R_act = R;
        dt_m = tF(i)-tB(ind_measB(i));
        
        if dt_m > dt_tol && m_sync
            t_span = [tB(ind_measB(i)), tF(i)];
            
            y_int = y_act;
            
            a_int = x_m(16:18);
            
            [~, y_act] = ode45(@tbp_m, t_span, y_int, opt, a_int);
            y_act = y_act(end,:)';
            
            % Measurement noise covariance time update
            stm_m = kep_stm_m(y_int, dt_m, data.env.mi);
            R_act = stm_m*R*stm_m';
        end
        
        % Gain update
        K = P_m*HB'/(HB*P_m*HB' + R_act);
        
        % Correction
        y_m = y_act - HB*x_m';
        x_p = x_m' + K*y_m;
        x_p = x_p';
        W = eye(18)-K*HB;
        P_p = W*P_m*W'+K*R_act*K';
        
        % Measurement marked as used
        ind_measB_prev = ind_measB(i);
        
        countB = countB +1;
        
%         tF(i)
%         tB_rec(ind_measB(i))
    end
    
    if log_measd(i) && ind_measd(i) ~= ind_measd_prev && tF(i)>=td(ind_measd(i)) && m_d
        
        % Measurement linearization
        y_act = measd(ind_measd(i),:)';
        R_act = Rd;
        dr = x_m(10:12)-x_m(1:3);
        
        Hd_act = 2*dr*Hd;
        % Gain update
        K = P_m*Hd_act'/(Hd_act*P_m*Hd_act' + R_act);
        
        % Correction
        y_m = y_act - dr*dr';
        x_p = x_m' + K*y_m;
        x_p = x_p';
        W = eye(18)-K*Hd_act;
        P_p = W*P_m*W'+K*R_act*K';
        
        % Measurement marked as used
        ind_measd_prev = ind_measd(i);
        
        countd = countd +1;
        
        %         tF(i)
        %         td(ind_measd(i))
        
    end
    
    x_prev = x_p;
    P_prev = P_p;
    
    x(i,:) = x_p;
    P(:,:,i) = P_p;
end

%% Post-processing
% Computation of accelerations not included in filter state propagation model
% (main term rm)
[a_d_A, a_d_B] = comp_acc_d(tstate);
% Modified true state for comparison reasons
mtstate = [tposA, tvelA, a_d_A, tposB, tvelB, a_d_B];

% Estimated relative state computation
xrel = zeros(nF,9);
for i = 1:nF
    xrel(i,:) = relstate(x(i,:));
end

%% Plots
list = {'Yes','No'};
[indx,~] = listdlg('PromptString','Close previous plots?','ListString',list,'SelectionMode','single');

switch indx
    case 1
        close all
end

list = {'Absolute','Relative','Absolute - Error','Relative - Error',};
[indx,~] = listdlg('PromptString','Select plot type:','ListString',list,'SelectionMode','single');

list = {'Yes','No'};
[cov_plot,~] = listdlg('PromptString','Include covariance bounds?','ListString',list,'SelectionMode','single');

ylabel_list = {'r_x^A [km]','r_y^A [km]','r_z^A [km]','v_x^A [km/s]','v_y^A [km/s]','v_z^A [km/s]','aD_x^A [km/s^2]','aD_y^A [km/s^2]','aD_z^A [km/s^2]',...
               'r_x^B [km]','r_y^B [km]','r_z^B [km]','v_x^B [km/s]','v_y^B [km/s]','v_z^B [km/s]','aD_x^B [km/s^2]','aD_y^B [km/s^2]','aD_z^B [km/s^2]'};

rel_title_list = {'Relative Position','Relative Velocity','Relative Acceleration'};

switch indx
    
    case 1
        n_pl = 1;
        for i = 1:6
            nfigure(i,3,2)
            for k = 1:3
                subplot(3,1,k)
                hold on
                grid on
                plot(ttime, mtstate(:,n_pl))
                plot(tF, x(:,n_pl))
                xlim([-inf,inf])
                xlabel('Time [s]')
                ylabel(ylabel_list(n_pl))
                n_pl = n_pl + 1;
            end
        end
        
    case 2
        n_pl = 1;
        for i = 1:3
            nfigure(i,3,1)
            for k = 1:3
                subplot(3,1,k)
                hold on
                grid on
                plot(ttime, trelstate(:,n_pl))
                plot(tF, xrel(:,n_pl))
                xlim([-inf,inf])
                xlabel('Time [s]')
                n_pl = n_pl + 1;
            end
            sgtitle(rel_title_list(i))
        end
        
    case 3
        x_ds = interp1(tF,x,ttime);
        n_pl = 1;
        for i = 1:6
            nfigure(i,3,2)
            for k = 1:3
                subplot(3,1,k)
                hold on
                grid on
                plot(ttime, x_ds(:,n_pl)-mtstate(:,n_pl))
                
                switch cov_plot
                    case 1
                        sig_p = sqrt(reshape(P(n_pl, n_pl, :), 1, nF));
                        plot(tF,  sig_p.*[-1; 1], 'r--')
                end
                
                xlim([-inf,inf])
                xlabel('Time [s]')
                ylabel(ylabel_list(n_pl))
                n_pl = n_pl + 1;
            end
        end
        
    case 4
        xrel_ds = interp1(tF,xrel,ttime);
        n_pl = 1;
        for i = 1:3
            nfigure(i,3,1)
            for k = 1:3
                subplot(3,1,k)
                hold on
                grid on
                plot(ttime, xrel_ds(:,n_pl)-trelstate(:,n_pl))
                
                if i == 1
                    switch cov_plot
                        case 1
                            covA  = reshape(P(n_pl,   n_pl,   :), 1, nF);
                            covB  = reshape(P(n_pl+9, n_pl+9, :), 1, nF);
                            covAB = reshape(P(n_pl,   n_pl+9, :), 1, nF);
                            sig_p = sqrt(covA+covB-2*covAB);
                            plot(tF,  sig_p.*[-1; 1], 'r--')
                    end
                end
                
                xlim([-inf,inf])
                xlabel('Time [s]')
                legend
                n_pl = n_pl + 1;
            end
            sgtitle(rel_title_list(i))
        end
        
end
