
% Modeling and Simulation of Aerospace Systems (2020/2021)

% Project - Drag free control of GOCE Spacecraft 

% Group 9

% Authors: Enrico Bassissi
%          Alessandro Colombo
%          Maria Alessandra De Luca

%% -------------------- Setup for default options ------------------------- %%

set(0, 'DefaultTextFontSize', 15);
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendFontSize', 15)
set(0, 'DefaultLineLineWidth', 1.4)
format short

colors = [0    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        0.6350    0.0780    0.1840
        0.7       0         0];

%% ----------------------- INIZIALIZATION ------------------------------- %%
clearvars -except colors; close all; clc;
addpath functions

% --------------------- SPACECRAFT'S SYSTEMS DATA --------------------------%
data = GoceData;

% -------------------- ORBIT INITIAL CONDITIONS ---------------------------%
% Initial keplerian elements
alt = 254.9;                                   % altitude [km] 
rp0 = data.earth.a_earth + alt;                % Perigee radius [km]
e0  = 0.0045;                                  % Eccentricity [-]
a0  = rp0/(1 - e0);                            % Semi_major axis [km]
% i0  = deg2rad(90);                             % Inclination [rad]
i0  = deg2rad(45);                             % Inclination [rad]

h0  = sqrt(data.earth.mu_E*a0*(1 - e0^2));     % Angular momentum [km^2/s] 
ra0 = h0^2/data.earth.mu_E/(1 - e0);           % Apogee radius [km]
OM0 = deg2rad(0);                              % Right-Ascension of the ascending node [rad]
om0 = deg2rad(0);                              % Argument of periapsis [rad]
th0 = deg2rad(0);                              % True anomaly [rad]
T = 2*pi*sqrt(a0^3/data.earth.mu_E);           % Orbital Period [s]

% -------------------- S/C SYSTEMS INITIAL CONDITIONS ------------------------%

% Orbital elements
coe0 = [a0 e0 i0 OM0 om0 th0];   

% Accelerometer device
xa0 = 0; va0 = 0; Vout0 = 0; 
acc0 = [xa0, va0, Vout0]; 

% Xenon Control Valve
Vout_int0 = 0; I0 = 0; xv0 = 10*data.valve.A0 + data.valve.d0; vv0 = 0;  
valve0 = [Vout_int0,I0, xv0, vv0];

% Overall initial conditions vector
ICS = [coe0, acc0, valve0];

% Integration condition
cond = struct(); 
%{
Before running the integration, specify: - cond.type -> choose among {'Nominal','Off-Nominal'}
in case of Off-Nominal, specify:         - cond.failure.time  
                                          - cond.failure.mode -> choose among: {'Grid-erosion', 'Xe-Leakage','PI failure'}
%}
%% ------------------------- INTEGRATION (ode15s) ---------------------- %%
 
% Set the integration options
options = odeset('RelTol',1e-10,'AbsTol',[1e-12*ones(6,1); ...
    1e-16*ones(2,1); 1e-10; 1e-6; 1e-10; 1e-6; 1e-7]);

% Integration scheme chosen: ode15s
tic
cond.type = 'Nominal';
cond.failure.mode = [];
[t_pl, z] = ode15s(@simulation, [0, 5*T], ICS , options, data,cond);
time_ode15s = toc;

%% ---------- EXTRACTION OF ADDITIONAL CHARACTERISTIC PARAMETERS ---------%

parout = zeros(length(t_pl),11);
for k = 1:length(t_pl) 
    [~, parout(k,:)] = simulation(t_pl(k), z(k, :) , data,cond);
end

Th = parout(:,1)*(1e+3);  % [mN]
D  = parout(:,2)*(1e+3);  % [mN]
A_valve = parout(:,3);    % [m^2]
lon = parout(:,4);        % [deg]
lat = parout(:,5);        % [deg]
height = parout(:,6);     % [km]
position = parout(:,7);   % [km]
velocity = parout(:,8);   % [km/s]
accel_at = parout(:,9);   % [km/s^2]
accel_an = parout(:,10);  % [km/s^2]
accel_ah = parout(:,11);  % [km/s^2]

%% ------------------------------- PLOTS --------------------------------- %%

%{
Plots performed:
(1) 3D orbit evolution, (2) Ground track ,(3) Thrust vs Drag, (4) Keplerian elements, 
(5) Altitude, (6) Proof Mass Displacement (with zoom out), (7) Proof mass
Velocity (with zoom out),(8) Valve control voltage, (9) Integral of valve control voltage, 
(10) Soleonid current, (11) Valve Spool Displacement, (12) Valve Spool Velocity
(13) Non keplerian state evolution design for latex report
%}

plot_orbit(z, parout,data)
plot_diagnostic(z, parout, t_pl, data, colors)

%% ----------------- INTEGRATION COMPARISON (ode15s vs ode23tb) -------- %%

tic
cond.type = 'Nominal';
cond.failure.mode = [];
[t_pl2, z2] = ode23tb(@simulation, [0, 5*T], ICS , options, data,cond);
time_ode23tb = toc;

Integrators_time = table(time_ode15s, time_ode23tb);
disp('ode15s is the fastest method: it will be used for the following integrations')

%% ------------------------- LINEARIZATION ----------------------------- %%

% Definition of the reference condition for linearization (chosen as the initial one)
ICS2 = [position(1) velocity(1) accel_at(1) accel_an(1) accel_ah(1) h0 data.earth.mu_E data.sc.M  D(1)*(1e-3)];
ICS3 = [cell2mat(struct2cell(data.accelerometer))' data.epsilon];
ICS_valve = struct2cell(data.valve);
ICS4 = [data.valve.A0 data.valve.r0 cell2mat(ICS_valve(6:end))'];
ICS5 = cell2mat(struct2cell(data.xenon))';
ICS6 = cell2mat(struct2cell(data.thruster))';

% Computation of the state matrices A linearizing GPEs and DFACs equations
% (both in symbolic form and already evaluated in the reference condition)
[A_matrix_sym_orb, A_matrix_values_orb,A_matrix_sym_con,A_matrix_values_con] = linearization(ICS, ICS2, ICS3, ICS4, ICS5, ICS6);

% Recovering of eigenvalues 
lambda_orb = eig(A_matrix_values_orb);
lambda_con = eig(A_matrix_values_con);
lambda = [lambda_orb; lambda_con];

% Plot of eigenvalues
figure('Name','Eigenvalues')
subplot(2,1,1)
plot(real(lambda_orb),imag(lambda_orb),'o','Color',colors(1,:))
xlabel('Re $ \{\lambda \} $'); ylabel('Im $ \{\lambda \} $')
title('Orbital mechanics eigenvalues')
xlim([-4*1e-9 0.2*1e-9])

subplot(2,1,2)
plot(real(lambda_con),imag(lambda_con),'o','Color',colors(2,:))
xlabel('Re $ \{\lambda \} $'); ylabel('Im $ \{\lambda \} $')
title('Control system eigenvalues')

% Recover step-size used by ode15s to integrate 
for i = 2:length(t_pl)
    h_step(i-1) = t_pl(i) - t_pl(i-1);
end
figure('Name', 'Stepsize')
plot(h_step,'o-','markersize',3)
xlabel('Step No.'); ylabel('Step size [s]'); grid on;

h_max = max(h_step);
h_min = min(h_step);

% Plot ode15s regions of stability and step-size hisogram 
% (Stability Regions of ODE Formulas, Nick Trefethen, February 2011, Mathworks)

figure('Name', 'BDF and step-size')
subplot(2,1,1)
x_rs = [0 0]; y_rs = [-8 8]; 
t = linspace(0,2*pi,1000);
z_rs = exp(1i*t);

plot(8*y_rs,x_rs,'k'), hold on, plot(x_rs,8*y_rs,'k')
d_rs = 1-1./z_rs; r_rs = 0;
for i = 1:5
  r_rs = r_rs+(d_rs.^i)/i;
  p(i) = plot(r_rs);
end
axis([-15 35 -25 25]), 
axis equal
ylim([-13 13])
hold on
title('BDF orders 1-5')
xlabel('Re $ \{h \lambda \} $'); ylabel('Im $ \{h \lambda \} $')
legend([p(1) p(2) p(3) p(4) p(5)],{'ord.1','ord.2','ord.3','ord.4','ord.5'},'Location','West')

subplot(2,1,2)
histogram(h_step)
xlabel('step-size [s]'); ylabel('BinCount')

%% ----------------------- OFF-NOMINAL CONDITIONS ------------------------- %% 

cond.failure.time = T + rand*2*T;

cond.type = 'Nominal';
cond.failure.mode = [];
[t_1, z_1] = ode15s(@simulation, [0 cond.failure.time],ICS , options, data,cond);
parout_1 = zeros(length(t_1), 11);
for k = 1:length(t_1) 
    [~, parout_1(k,:)] = simulation(t_1(k), z_1(k, :) ,data,cond);
end

%------------------------- OFF-NOMINAL 1 ---------------------------------%

cond.type = 'Off-Nominal';
cond.failure.mode = 'Grid-erosion';
[t_off1, z_off1] = ode15s(@simulation, [cond.failure.time 5*T], z_1(end,:) , options, data,cond);
parout_off1 = zeros(length(t_off1), 11);
for k = 1:length(t_off1) 
    [~, parout_off1(k,:)] = simulation(t_off1(k), z_off1(k, :) , data,cond);
end

%------------------------- OFF-NOMINAL 2 ---------------------------------%

cond.type = 'Off-Nominal';
cond.failure.mode = 'Xe-Leakage';
[t_off2, z_off2] = ode15s(@simulation, [cond.failure.time 5*T], z_1(end,:) , options, data,cond);
parout_off2 = zeros(length(t_off2), 11);
for k = 1:length(t_off2) 
    [~, parout_off2(k,:)] = simulation(t_off2(k), z_off2(k, :), data,cond);
end

%--------------------- OFF-NOMINAL 3 with RECOVERY -----------------------%
% Failure
cond.failure.time2 = cond.failure.time + rand*T;
cond.type = 'Off-Nominal';
cond.failure.mode = 'PI failure';
options_off3 = odeset('RelTol',1e-10,'AbsTol',[1e-12*ones(6,1); ...
    1e-16*ones(2,1); 1e-10; 1e-6; 1e-10; 1e-6; 1e-7],'event',@event_off3);
[t_off31, z_off31,te,xve,ie] = ode15s(@simulation, [cond.failure.time cond.failure.time2],z_1(end,:), options_off3,data,cond);

parout_off31 = zeros(length(t_off31), 11);
for k = 1:length(t_off31) 
    [~, parout_off31(k,:)] = simulation(t_off31(k), z_off31(k, :) , data,cond);
end

if ie == 1 % event occurred
    z_off31_new = z_off31;
    z_off31_new(end,11) = 0;
    z_off31_new(end,13) = 0;
    [t_off32, z_off32] = ode15s(@simulation, [te cond.failure.time2],z_off31_new(end,:), options,data, cond);
else % event didn't occur
    t_off32 = zeros(size(t_off31));
    z_off32 = zeros(size(z_off31));
end
    
parout_off32 = zeros(length(t_off32), 11);
for k = 1:length(t_off32) 
    [~, parout_off32(k,:)] = simulation(t_off32(k), z_off32(k, :) , data,cond);
end

% Recovery
cond.type = 'Nominal';
cond.failure.mode = [];
[t_off33, z_off33] = ode15s(@simulation, [cond.failure.time2 5*T],z_off32(end,:),options,data,cond);

parout_off33 = zeros(length(t_off33), 11);
for k = 1:length(t_off33) 
    [~, parout_off33(k,:)] = simulation(t_off33(k), z_off33(k, :) , data,cond);
end

t_off3 = [t_off31;t_off32(2:end,:);t_off33(2:end,:)];
z_off3 = [z_off31;z_off32(2:end,:); z_off33(2:end,:)];
parout_off3 = [parout_off31; parout_off32(2:end,:);parout_off33(2:end,:)];

% ------------------------------- PLOTS --------------------------------- %

plot_diagnostic_off(z,t_pl,parout,[z_1;z_off1(2:end,:)],[t_1;t_off1(2:end)],[parout_1; parout_off1(2:end,:)],...
              [z_1;z_off2(2:end,:)],[t_1;t_off2(2:end)],[parout_1; parout_off2(2:end,:)],[z_1;z_off3(2:end,:)],...
              [t_1;t_off3(2:end)],[parout_1; parout_off3(2:end,:)], data,colors)

% Plot Off-Nominal condition used in the report
plot_diagnostic_off_report(z,t_pl,parout,[z_1;z_off1(2:end,:)],[t_1;t_off1(2:end)],[parout_1; parout_off1(2:end,:)],...
             [z_1;z_off2(2:end,:)],[t_1;t_off2(2:end)],[parout_1; parout_off2(2:end,:)],[z_1;z_off3(2:end,:)],...
             [t_1;t_off3(2:end)],[parout_1; parout_off3(2:end,:)], data, colors)

%% ------------------- SENSITIVITY ANALYSIS ----------------------------- %%

% Define the parameters used for the sensitivity analysis
x0 = [data.accelerometer.Kpa, data.accelerometer.Kda, data.valve.mfcv, data.valve.Ki, ...
    data.valve.Kpv, data.valve.Kiv, data.valve.kfcv, data.xenon.p2, data.xenon.T2, data.valve.R, data.valve.L];

% Define the range of variations of the parameters
x1 = x0*0.4; % -60% di x0
x2 = x0*0.6; % -40% di x0
x3 = x0*0.8; % -20% di x0
x4 = x0;     % nominal
x5 = x0*1.2; % +20% di x0
x6 = x0*1.4; % +40% di x0
x7 = x0*1.6; % +60% di x0
x = [x1;x2;x3;x4;x5;x6;x7];

[xrows,xcols] = size(x);

best_mod = zeros(1,xcols); % array in which I will change the value of parameters
fji = zeros(xrows, xcols); % array in which I will put the evaluation of th obj function for varying i

% Initialization of the waitbar, to be updated in the for loop
barra1 = waitbar(0);
iter_bar = 0;

% For loop to evaluate the objective function at each individual change of
% parameters
for i =1:xcols % for parameters
    best_mod = x0;
    for j = 1:xrows % for variations of nominal value
        iter_bar = iter_bar + 1; % index for the waitbar
        waitbar(iter_bar/(xrows*xcols),barra1,sprintf('Sensitivity step %d of %d', ...
                iter_bar,(xrows*xcols)))
        best_mod(i) = x(j,i); % set of parameters, one by one changed
        fji(j,i) = ps_obj_sens(best_mod, T, ICS, options, data); %evaluation objective function
    end
end

% Standard deviation of obj function evaluations applied by columns
std_fji = std(fji);
% Mean value of obj function evaluations applied by columns
mean_fji = mean(fji);

% Plot showing results of sensitivity analysis, mean value in red, all
% values in blue, box representing the std in yellow
figure('Name','Std Dist')
for i = 1:xcols
    h_sens1 = plot(i,fji(:,i),'.','markersize',15,'Color',colors(1,:));
    hold on;
    h_sens2 = plot(i,mean_fji(i),'.','markersize',15,'Color',colors(2,:));
    h_sens_line = line([i-1/5,i+1/5],[mean_fji(i)+std_fji(i),mean_fji(i)+std_fji(i)],'Color',colors(3,:));
    line([i-1/5,i+1/5],[mean_fji(i)-std_fji(i),mean_fji(i)-std_fji(i)],'Color',colors(3,:))
    line([i-1/5,i-1/5],[mean_fji(i)-std_fji(i),mean_fji(i)+std_fji(i)],'Color',colors(3,:))
    line([i+1/5,i+1/5],[mean_fji(i)-std_fji(i),mean_fji(i)+std_fji(i)],'Color',colors(3,:))
end
ylabel('Objective Function'); grid on;
legend([h_sens1(1), h_sens2, h_sens_line], 'Evaluations', 'Mean Value', 'std')
set(gca, 'XTick', 1:xcols)
xticklabels({'$K_{pa}$', '$K_{da}$', '$m_{fcv}$', '$K_i$', '$K_{pv}$',...
    '$K_{iv}$','$k_{fcv}$', '$P_2$', '$T_2$','$R_v$', '$L_v$'})

%% -------------------- OPTIMIZATION PATTERNSEARCH ------------------------ %%
data = GoceData;

% Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

% Patternsearch options
opt_ps = optimoptions('patternsearch','Cache','on' ,'Display','iter','PlotFcn',@psplotbestf, ...
                      'UseParallel', true, 'UseCompletePoll', true, 'UseCompleteSearch',true);

% Function to be optimized
fun = @(x) ps_obj(x, T, ICS, options, data);

% Nominal conditions of the parameters to be optimized
x0 = [data.accelerometer.Kpa, data.valve.Ki, data.valve.Kiv, data.valve.R];

% Lower - upper boundaries
lb = x0*0.5e-1;
ub = x0.*0.2e1;
A = [];
b = [];
Aeq = [];
beq = [];

% Call of patternsearch
[x_opt,fval,exitflag,output] = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub, opt_ps);

INTEST = {'ps'; 'baseline'};
KPA = [x_opt(1); data.accelerometer.Kpa];
KI = [x_opt(2); data.valve.Ki];
KIV = [x_opt(3); data.valve.Kiv];                                
RV = [x_opt(4); data.valve.R];

Tbl = table(INTEST, KPA, KI, KIV, RV);
Tbl.Properties.VariableNames = {'-','$K_{pa}','$K_i$','$K_{iv}$','$R_v$'};

%% ------------ ANALYSIS USING OPTIMISED PARAMETERS -------------------- %%
% Substitute optimised data 
data.accelerometer.Kpa = x_opt(1);   
data.valve.Ki = x_opt(2);     
data.valve.Kiv = x_opt(3);
data.valve.R = x_opt(4);

% Integration with optimised data
cond.type = 'Nominal';
cond.failure.mode = [];
options = odeset('RelTol',1e-8,'AbsTol',[1e-8*ones(6,1); ...
    1e-16*ones(2,1); 1e-10; 1e-6; 1e-10; 1e-6; 1e-7]);  

[t_pl_opt, z_opt] = ode15s(@simulation, [0, 5*T], ICS , options, data,cond);

parout_opt = zeros(length(t_pl_opt),11);
for k = 1:length(t_pl_opt) 
    [~, parout_opt(k,:)] = simulation(t_pl_opt(k), z_opt(k, :) , data,cond);
end

% ------------------------------- PLOTS --------------------------------- %
plot_diagnostic(z_opt, parout_opt, t_pl_opt, data, colors);

%% ------------------------- ERROR DFAC MAX ---------------------------- %%

% To understand the effect of the optimization, it's shown the error
% between the thrust and the drag comparing between nominal and optimised case

% To eliminate the transient, the comparison analysis start from about 1000 s
% after the starting moment
diff_opt = abs(parout_opt(415:end,1)*1e3 + parout_opt(415:end,2)*1e3); 
diff_nom = abs(Th(510:end) + D(510:end)); 

figure('Name','ERROR DFAC MAX')
plot(t_pl_opt(415:end)/3600,diff_opt, t_pl(510:end)/3600, diff_nom)
legend('Optimal','Nominal')
grid on; xlabel('$t$ [h]'); ylabel('$|$T-D$|$ [mN]');

% maximum error evaluation
max_diff = [max(diff_opt); max(diff_nom)];
INTEST = {'Opt'; 'Nom'};
Tbl = table(INTEST, max_diff);
Tbl.Properties.VariableNames = {'-','max(err)'}

%% --------------------- VIDEO SECTION ----------------------------------%%
%{ 
Nominal condition video evolution
Video performed : (1) 3D Orbit, (2) Valve displacement, (3) Proof mass dislacement
%}
movie_orbit(t_pl, T, z, data, colors)
evolutions_video(t_pl, z, data)

%{ 
Off-Nominal condition video evolution
Only for PI failure condition
Video performed : (1) Valve displacement, (2) Proof mass dislacement
%}
evolutions_video([t_1; t_off3], [z_1; z_off3], data)

%{ 
Optimised parameters condition video evolution
Video performed : (1) Valve displacement, (2) Proof mass dislacement
%}
evolutions_video(t_pl_opt, z_opt, data)
