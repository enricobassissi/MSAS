function data = GoceData

%{ 
Function that returns all the GOCE parameters.

 INPUT:  -

 OUTPUT: -

 FUNCTIONS REQUIRED: -

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}

% ------------------------- EARTH DATA ----------------------------- %
data.earth.RE = astroConstants(23);            % Earth’s mean radius [km]
data.earth.mu_E = astroConstants(13);          % Earth’s gravitational parameter [km^3/s^2]
data.earth.a_earth = 6378.16;                  % equatorial radius   [km]
data.earth.b_earth = 6356.778;                 % polar radius        [km]

% ------------------------- SPACECRAFT DATA ----------------------------- %

data.initial_date = [2009, 1, 1];
data.sc.Cd = 3;
data.sc.A  = 1.1;
data.sc.M  = 300;
data.sc.B  = 1/300;

% -------------------------- ACCELEROMETER DATA ------------------------- %

data.accelerometer.m     = 0.32;      % [kg]
data.accelerometer.gap   = 5e-4;      % [m]
data.accelerometer.Vbias = 10;        % [V]
data.accelerometer.Cf    = 2e-12;     % [F] op amp capacitance
data.accelerometer.Kpa   = 1e6;       % optimizable, proportional accelemeter control
data.accelerometer.Kda   = 5e4;       % optimizable, derivative gain accel control
data.accelerometer.A     = 1.6e-3;    % [m^2],  accelerometer seismic mass section
data.epsilon             = 8.854e-12; % [C/(V m)] dielectric constant vacuum


% ---------------------------- VALVE DATA ------------------------------- %

data.valve.A0 = 1e-5;                     % [m^2]
data.valve.d0 = sqrt(data.valve.A0*4/pi); % [m]
data.valve.r0 = data.valve.d0/2;          % [m]
data.valve.alpha   = @(zeta) 2*acos(zeta/data.valve.r0 - 1); 
data.valve.A_valve = @(zeta) data.valve.r0^2/2*( zeta - sin(zeta)); %[m^2]
data.valve.mfcv = 2e-1; % [kg], optimizable
data.valve.kfcv = 7e3;  % [N/m]
data.valve.c   = 30;    % [Ns/m]
data.valve.Ki  = 0.2;   % Gain optimizable
data.valve.Kpv = 0.1;   % Gain optimizable
data.valve.L   = 1e-3;  % [H]
data.valve.R   = 0.5;   % [ohm]
data.valve.Kiv = 3;     % optimizable


% ---------------------------- XENON DATA ------------------------------- %

data.xenon.T2 = 240;     % [K]
data.xenon.p2 = 2e5;     % [Pa]
data.xenon.gamma = 1.66; 
data.R = 8314.47;        % [J/(kmol K)]
data.xenon.MW  = 131.29; % [kg/kmol]
data.xenon.RXe = data.R/data.xenon.MW; % [J/(kg K)]
data.xenon.mXe = 2.188e-25; % [kg]
data.xenon.rho = data.xenon.p2/(data.xenon.T2*data.xenon.RXe); % density [kg/m^3]
data.xenon.storageT = 40+273.15; % [K]
data.xenon.storageP =  data.xenon.p2*(data.xenon.T2/data.xenon.storageT)^(data.xenon.gamma/(1-data.xenon.gamma)); % [Pa]

% ------------------------- THRUSTER DATA ------------------------------- %

data.thruster.e = 1.6e-19; % [C]
data.thruster.dV = 2e3;    % [V]

end