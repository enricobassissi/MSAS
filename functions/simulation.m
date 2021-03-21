function [dydt, parout] = simulation(t, y, data, cond)
%{
 Function that integrates the GOCE satellite motion and its components by keeping 
 into account Earth J2 and Atmospheric Drag perturbations (modeled through Gauss 
 variational principle) and the thrust produced by the ion thruster propulsion system.
 

 INPUT:  1. t: Time (can be omitted, as the system is autonomous) (s) 
         2. y: State variables: 
               % Keplerian Elements: [a e i OM om th] 
               % Accelerometer variables: mass position and velocity [xa va] -  Voltage [Vout] 
               % Control valve variables: Integral of the Voltage [Vout_int] - current [I] - position and velocity of the solenoidal actuator [xv  vv]
        
 OUTPUT: 1. dydt: Derivative of the states [a_dot e_dot i_dot OM_dot om_dot th_dot xa_dot va_dot Vout_dot Vout I_dot xv_dot vv_dot]
         2. parout: Additional quantities of interest [Th D A_valve lon lat height r v at an ah]

 FUNCTIONS REQUIRED: astroConstants, sv_from_coe, lonlat 

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}

% STATES
 % Orbital mechanics
 a     = y(1);     %[km]
 e     = y(2);     %[-]
 i     = y(3);     %[rad]
 OM    = y(4);     %[rad]
 om    = y(5);     %[rad]
 th    = y(6);     %[rad]
 % Accelerometer 
 xa    = y(7);     %[m]
 va    = y(8);     %[m/s]
 Vout  = y(9);     %[V]
 % Control valve
 Vout_int = y(10); %[V/s] 
 I     = y(11);    %[A]
 xv    = y(12);    %[m]
 vv    = y(13);    %[m/s]


% ----------------------- ORBITAL MECHANICS PART ------------------------ % 

% Parameters
J2 = 0.00108263;               % J2 effect
Re = astroConstants(23);       % Mean radius of the Earth [km]
phi = om + th;                 % Argument of latitude [rad]
%phi = mod(phi,2*pi);           
mu_E= astroConstants(13);      % Earth planetary constant [km^3/s^2]
h = sqrt(mu_E*a*(1-e^2));      % Angular momentum [km^2/s]

% Position and velocity of the s/c in ECI
[R, V] = sv_from_coe(y,mu_E); 

% Definition of {tnh} reference frame
r = norm(R);
v = norm(V);
H = cross(R,V); 

t_hat = V/v;  
h_hat = H/norm(H);
n_hat = cross(h_hat,t_hat);
Rot   = [ t_hat  n_hat  h_hat ];

% Computation of longitude and latitude
[lon, lat] = lonlat(R,t);

% Radius of the earth changing at each lat position - Earth modeled as an
% ellipsoid
a_earth = 6378.137;  % equatorial radius [km]
b_earth = 6356.7523; % polar radius      [km]
lat_rad = deg2rad(lat);
R_lat = sqrt(((a_earth^2.*cos(lat_rad)).^2+(b_earth^2.*sin(lat_rad)).^2)./...
            ((a_earth.*cos(lat_rad)).^2+(b_earth.*sin(lat_rad)).^2));

% Perturbation Atm Drag components in the rsw frame: 
% Calling of Atmospheric Model function

% Magnetic Field parameters
aph  = [7.3 6.0 9.0 12.0 9.0 10.1 0.5];
f107 = 67.0;
f107a =  67.3;
flags = zeros(1,23);

% s/c altitude
height = r - R_lat; 
 
% [T, RHO] = ATMOSNRLMSISE00( H (m), LAT (deg), LON, YEAR, DOY, SEC, LST, F107A,
% F107, APH, FLAGS, ITYPE, OTYPE, ACTION )
% Atmosphere density at the given altitude
[~, rho_v] = atmosnrlmsise00(height*1e3, lat, lon, data.initial_date(1), data.initial_date(3), t, f107a, f107, aph, flags, 'NoOxygen', 'Warning');
rho = rho_v(6); % [kg/m^3]

% Earth angular velocity
we = (2*pi + 2*pi/365.26)/(24*3600); %[km/s]
WE = [0;0;we];

% Rotation matrix from ECI to ECEF frame
Rot2 = [cos(we) sin(we) 0;
      -sin(we) cos(we) 0;
      0 0 1];
        
% s/c relative velocity 
V_REL = V - cross(WE,R); 
v_rel = norm(V_REL);

% Perturbing acceleration due to atmospheric drag: absolute value and
% components
aD_norm = -1/2*rho*data.sc.B*v_rel^2*(1e+3);     % [km/s^2] 
aD_1 = -1/2*rho*data.sc.B*v_rel*V_REL(1)*(1e+3); % [km/s^2] 
aD_2 = -1/2*rho*data.sc.B*v_rel*V_REL(2)*(1e+3); % [km/s^2] 
aD_3 = -1/2*rho*data.sc.B*v_rel*V_REL(3)*(1e+3); % [km/s^2] 

% Drag force acting on the sc
D = data.sc.M*aD_norm*(1e+3); % [N]

% Position of the s/c in ECEF frame
s_sc = Rot2*R;  
        
aJ2_1_ecef = 3/2*J2*mu_E*Re^2/r^4*(s_sc(1)/r *(5 *s_sc(3)^2/r^2 -1)); %km/s^2
aJ2_2_ecef = 3/2*J2*mu_E*Re^2/r^4*(s_sc(2)/r *(5 *s_sc(3)^2/r^2 -1)); %km/s^2
aJ2_3_ecef = 3/2*J2*mu_E*Re^2/r^4*(s_sc(3)/r *(5 *s_sc(3)^2/r^2 -3)); %km/s^2

% Rotation in ECI frame:
aJ2 = Rot2'*[aJ2_1_ecef; aJ2_2_ecef; aJ2_3_ecef;];
aJ2_1 = aJ2(1);
aJ2_2 = aJ2(2);
aJ2_3 = aJ2(3);

% ------------------- ORBITS TO BE CONTINUED LATER ---------------------- %

% ------------------------ ACCELEROMETER PART --------------------------- %

% Check on the position of the mass
if xa > data.accelerometer.gap
    xa = data.accelerometer.gap;
end
if xa >= data.accelerometer.gap && va > 0
    xa = data.accelerometer.gap;
    va = 0;
end
if xa < -data.accelerometer.gap
    xa = -data.accelerometer.gap;
end
if xa <= -data.accelerometer.gap && va < 0
    xa = -data.accelerometer.gap;
    va = 0;
end

% Voltage drop due to the displacement of the mass
Vx = xa/data.accelerometer.gap*data.accelerometer.Vbias;

% Readout voltage
Vout_dot = -1/data.accelerometer.Cf*2*data.accelerometer.Vbias*va*data.epsilon*data.accelerometer.A*...
    (data.accelerometer.gap^2 + xa^2)/((data.accelerometer.gap^2 - xa^2)^2);

% Voltage obtained from the PD controller
Vc = data.accelerometer.Kpa*Vout + data.accelerometer.Kda*Vout_dot;

% Electric forces
deltaV1 = data.accelerometer.Vbias - Vc - 1/2*Vx;
deltaV2 = data.accelerometer.Vbias + Vc + 1/2*Vx;

Fel1 = 1/2*data.epsilon*data.accelerometer.A/((data.accelerometer.gap - xa)^2)*deltaV1^2;
Fel2 = 1/2*data.epsilon*data.accelerometer.A/((data.accelerometer.gap + xa)^2)*deltaV2^2;

% ----------------- ACCELEROMETER TO BE CONTINUED LATER ----------------- %

% ----------------------------- VALVE PART ------------------------------ %

% Voltage obtained from the PI controller + introduction of possible off-nominal conditions (PI failure)
if strcmp(cond.type,'Nominal')
    VI = data.valve.Kpv*Vout + data.valve.Kiv*Vout_int;
elseif strcmp(cond.type,'Off-Nominal')
      if strcmp(cond.failure.mode,'Grid-erosion') ||  strcmp(cond.failure.mode,'Xe-Leakage')
          VI = data.valve.Kpv*Vout + data.valve.Kiv*Vout_int;
      elseif strcmp(cond.failure.mode,'PI failure') 
          VI = 0;
      end
end

% Solenoidal current state equation
I_dot = VI/data.valve.L - data.valve.R/data.valve.L*I;

% Computation of the orifice area
if xv < 10*data.valve.A0
    xv = 10*data.valve.A0;
end
if xv <= 10*data.valve.A0 && vv < 0
    xv = 10*data.valve.A0;
    vv = 0;
end    
if xv > (10*data.valve.A0 + data.valve.d0)
    xv = 10*data.valve.A0 + data.valve.d0;
end    
if xv >= (10*data.valve.A0 + data.valve.d0) && vv > 0
    xv = 10*data.valve.A0 + data.valve.d0;
    vv = 0;
end

alpha = data.valve.alpha(xv-10*data.valve.A0);
% Control on alpha angle
if sin(alpha/2) >= 0
    alpha = alpha;
else
    alpha = 2*pi - alpha;
end

A_valve = data.valve.A_valve(alpha); %[m^2]

% Valve state equations
xv_dot = vv;
vv_dot = 1/data.valve.mfcv*(data.valve.kfcv*(data.valve.d0 + 10*data.valve.A0 - xv)...
    - data.valve.Ki*I - data.valve.c*vv);

dvalve_dot = [I_dot, xv_dot, vv_dot];

% --------------------------- XENON TANK PART --------------------------- %

% Introduction of possible off-nominal condition (Xe-Leakage)
if strcmp(cond.type,'Off-Nominal')
      if strcmp(cond.failure.mode,'Xe-Leakage')
           P1 =  data.xenon.storageP*exp((cond.failure.time - t)/1e4);
           data.xenon.p2 = P1.*(data.xenon.storageT/data.xenon.T2)^(data.xenon.gamma/(1-data.xenon.gamma));
           data.xenon.rho = data.xenon.p2/(data.xenon.RXe*data.xenon.T2);
      end
end

% Xenon mass flow rate
GAMMA = sqrt(data.xenon.gamma)*(2/(data.xenon.gamma + 1))^((data.xenon.gamma+1)/(2*(data.xenon.gamma-1)));

m_dot = GAMMA*data.xenon.p2/sqrt(data.xenon.RXe*data.xenon.T2)*A_valve;


% --------------------------- THRUSTER PART ----------------------------- %

% Introduction of possible off-nominal condition (Grid erosion)
if strcmp(cond.type,'Off-Nominal')
      if strcmp(cond.failure.mode,'Grid-erosion')
           data.thruster.dV = data.thruster.dV - ((t - cond.failure.time)/70).^(4/3); % generic power law decay
      end
end

% Computation of the thrust
Th = m_dot*sqrt(2*data.thruster.dV*data.thruster.e/data.xenon.mXe);

% --------------------- CONTINUE ACCELEROMETER PART --------------------- %

% External force
Fext = (Th + D)/data.sc.M*data.accelerometer.m; 

% Proof mass state equations
xa_dot = va;
va_dot = 1/data.accelerometer.m*(Fext - Fel1 + Fel2);

dacc_dt = [xa_dot, va_dot, Vout_dot, Vout];

% ------------- ORBITAL PART CONTINUE WITH THRUST ------------------------%

% Propulsion system acceleration components 
aT_1 = Th*V_REL(1)/v_rel/data.sc.M*(1e-3); % [km/s^2] 
aT_2 = Th*V_REL(2)/v_rel/data.sc.M*(1e-3); % [km/s^2] 
aT_3 = Th*V_REL(3)/v_rel/data.sc.M*(1e-3); % [km/s^2] 

% Overall acceleration 
a_xyz = [aJ2_1+aD_1 + aT_1; aJ2_2+aD_2 + aT_2 ; aJ2_3+aD_3+aT_3]; 

% Rotation in [t,n,h] frame
a_tnh = Rot'*a_xyz;
at = a_tnh(1);
an = a_tnh(2);
ah = a_tnh(3);

% Gauss variational equations
a_dot = 2*a^2*v*at/mu_E;
e_dot = 1/v*(2*(e+cos(th))*at-r/a*sin(th)*an);
i_dot = (r*cos(phi)/h)*ah;
OM_dot = r*sin(phi)*ah/(h*sin(i));
om_dot = 1/(e*v)*(2*sin(th)*at + (2*e+r*cos(th)/a)*an) - r*sin(phi)*cos(i)*ah/(h*sin(i));
th_dot = h/r^2 - 1/(e*v)*(2*sin(th)*at + (2*e+r*cos(th)/a)*an );

dcoe_dt = [a_dot e_dot i_dot OM_dot om_dot th_dot];

% ------------------------------- OUTPUT -------------------------------- %

% Return rates to ode
dydt = [dcoe_dt, dacc_dt, dvalve_dot]';

% Output parameters
parout = [Th, D, A_valve,lon,lat,height,r, v, at, an, ah]; 

end

