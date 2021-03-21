function [A_matrix_sym_orb, A_matrix_values_orb,A_matrix_sym_con, A_matrix_values_con] = linearization(ICS, ICS2, ICS3, ICS4, ICS5, ICS6)
%{
 Function that linearizes the GPEs and DFACs equations around a given reference condition.

 INPUT:  1. ICS  : states initial conditions [a0 e0 i0 OM0 om0 th0 xa0 va0 Vout0 Vout_int0 I0 xv0 vv0]
         2. ICS2 : vector containing: s/c position - s/c velocity - [t n h] accelerations in ICs, earth planetary constant, mass of the s/c, Drag in ICs
         3. ICS3 : accelerometer data [m e_gap Vbias Cf Kpa Kda A epsilon]
         4. ICS4 : valve data [valve_A0 valve_r0 mfcv kfcv c Ki Kpv L R Kiv]
         5. ICS5 : xenon data [T2 p2 gammaXe MW RXe mXe rho T_stor P_stor ]
         6. ICS6 : thruster data [eXe dV ]
       
 OUTPUT: 1. A_matrix_sym_orb : state matrix of the linearized system in
            symbolic language (orbital part)
         2. A_matrix_values_orb : state matrix of the linearized system
            evaluated at ICs (orbital part)
         3. A_matrix_sym_con : state matrix of the linearized system in
            symbolic language (control part)
         4. A_matrix_values_con : state matrix of the linearized system
            evaluated at ICs (orbital part)

 FUNCTIONS REQUIRED: [-]

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}

% Definition of symbolic variables
syms a e incl OM om th xa va Vout Vout_int I xv vv           
syms a0 e0 i0 OM0 om0 th0 xa0 va0 Vout0 Vout_int0 I0 xv0 vv0 
syms r0 v0 at0 an0 ah0 h0  mu_E M D               
syms m e_gap Vbias Cf Kpa Kda A epsilon           
syms valve_A0 valve_r0 mfcv kfcv c Ki Kpv L R Kiv 
syms T2 p2 gammaXe MW RXe mXe rho T_stor P_stor  
syms eXe dV                                       

%%--------------------------- ORBITAL MECHANICS -------------------------%%

% Semi-major axis
a_dot = 2*a^2*v0*at0/mu_E;
a_dot_lin = taylor(a_dot, a, 'ExpansionPoint', a0, 'order', 2);
a_dot_a_part = diff(a_dot_lin, a);

A_matrix_sym_orb(1,1) = a_dot_a_part;

% Eccentricity
e_dot = 1/v0*(2*(e + cos(th))*at0-r0/a*sin(th)*an0);
e_dot_lin = taylor(e_dot, [a, e, th], 'ExpansionPoint', [a0, e0, th0], 'order', 2);
e_dot_a_part = diff(e_dot_lin, a);
e_dot_e_part = diff(e_dot_lin, e);
e_dot_th_part = diff(e_dot_lin, th);

A_matrix_sym_orb(2,1) = e_dot_a_part;
A_matrix_sym_orb(2,2) = e_dot_e_part;
A_matrix_sym_orb(2,6) = e_dot_th_part;

% Inclination
i_dot = (r0*cos(th + om)/h0)*ah0;
i_dot_lin = taylor(i_dot, [th, om], 'ExpansionPoint', [th0, om0], 'order', 2);

i_dot_om_part = diff(i_dot_lin, om);
i_dot_th_part = diff(i_dot_lin, th);

A_matrix_sym_orb(3,5) = i_dot_om_part;
A_matrix_sym_orb(3,6) = i_dot_th_part;

% RAAN
OM_dot = r0*sin(om + th)*ah0/(h0*sin(incl));
OM_dot_lin = taylor(OM_dot, [incl, om, th], 'ExpansionPoint', [i0, om0, th0], 'order', 2);
OM_dot_om_part = diff(OM_dot_lin, om);
OM_dot_th_part = diff(OM_dot_lin, th);
OM_dot_incl_part = diff(OM_dot_lin, incl);

A_matrix_sym_orb(4,3) = OM_dot_incl_part;
A_matrix_sym_orb(4,5) = OM_dot_om_part;
A_matrix_sym_orb(4,6) = OM_dot_th_part;

% Argument of perigee
om_dot = 1/(e*v0)*(2*sin(th)*at0 + (2*e0+r0*cos(th)/a)*an0) - r0*sin(om +th)*cos(incl)*ah0/(h0*sin(incl));
om_dot_lin = taylor(om_dot, [a, e, incl, om, th], 'ExpansionPoint',...
    [a0, e0, i0, om0, th0], 'order', 2);

om_dot_a_part = diff(om_dot_lin, a);
om_dot_e_part = diff(om_dot_lin, e);
om_dot_incl_part = diff(om_dot_lin, incl);
om_dot_om_part = diff(om_dot_lin, om);
om_dot_th_part = diff(om_dot_lin, th);

A_matrix_sym_orb(5,1) = om_dot_a_part;
A_matrix_sym_orb(5,2) = om_dot_e_part;
A_matrix_sym_orb(5,3) = om_dot_incl_part;
A_matrix_sym_orb(5,5) = om_dot_om_part;
A_matrix_sym_orb(5,6) = om_dot_th_part;

% True anomaly
th_dot = 0;
A_matrix_sym_orb(6,:) = 0;

% Alternative: Consider th_dot equation
% th_dot = h0/r0^2 - 1/(e*v0)*(2*sin(th)*at0 + (2*e+r0*cos(th)/a)*an0 );
% th_dot_lin = taylor(th_dot, [a, e, th], 'ExpansionPoint',...
%     [a0, e0, th0], 'order', 2);
% 
% th_dot_a_part = diff(th_dot_lin, a);
% th_dot_e_part = diff(th_dot_lin, e);
% th_dot_th_part = diff(th_dot_lin, th);
% 
% A_matrix_sym_orb(6,1) = th_dot_a_part;
% A_matrix_sym_orb(6,2) = th_dot_e_part;
% A_matrix_sym_orb(6,6) = th_dot_th_part;

% substitution of states initial conditions
A_matrix_sym_orb_subs1 = subs(A_matrix_sym_orb,{a0 e0 i0 OM0 om0 th0},{ICS(1:6)});

% substitution of initial position - velocity and accelerations
A_matrix_sym_orb_subs2 = subs(A_matrix_sym_orb_subs1,{r0 v0 at0 an0 ah0 h0 mu_E M D},{ICS2});

A_matrix_values_orb = double(A_matrix_sym_orb_subs2);   


%%--------------------------- CONTROL SYSTEM ----------------------------%%

% Voltage drop due to xa
Vx = xa/e_gap*Vbias;

% Readout voltage
Vout_dot = - 1/Cf*2*Vbias*va*epsilon*A*(e_gap^2 + xa^2)/((e_gap^2 - xa^2)^2);
Vout_dot_lin = taylor(Vout_dot, [xa,va],'ExpansionPoint', [xa0, va0],'order',2);
Vout_dot_xa_part = diff(Vout_dot_lin, xa);
Vout_dot_va_part = diff(Vout_dot_lin, va);

A_matrix_sym_con(3,1) = Vout_dot_xa_part;
A_matrix_sym_con(3,2) = Vout_dot_va_part;

% Valve section area
A_valve = valve_r0^2/2*( 2*acos((xv-10*valve_A0)/valve_r0 - 1) - sin(2*acos((xv-10*valve_A0)/valve_r0 - 1))); %[m^2]

% Voltage obtained from the PD controller
Vc = Kpa*Vout + Kda*Vout_dot;

% Electric forces
deltaV1 = Vbias - Vc - 1/2*Vx;
deltaV2 = Vbias + Vc + 1/2*Vx;
Fel1 = 1/2*epsilon*A/((e_gap - xa)^2)*deltaV1^2;
Fel2 = 1/2*epsilon*A/((e_gap + xa)^2)*deltaV2^2;

% Mass flow rate and thrust
GAMMA = sqrt(gammaXe)*(2/(gammaXe + 1))^((gammaXe+1)/(2*(gammaXe-1)));
m_dot = GAMMA*p2/sqrt(RXe*T2)*A_valve;
Th = m_dot*sqrt(2*dV*eXe/mXe);

% Velocity of proof mass
va_dot = 1/m*((Th + D)/M*m - Fel1 + Fel2);
va_dot_lin = taylor(va_dot,[xa, va, Vout,xv],'ExpansionPoint', [xa0, va0, Vout0, xv0],'order',2);
va_dot_xa_part = diff(va_dot_lin, xa);
va_dot_va_part = diff(va_dot_lin, va);
va_dot_Vout_part = diff(va_dot_lin, Vout);
va_dot_xv_part = diff(va_dot_lin, xv);

A_matrix_sym_con(2,1) = va_dot_xa_part;
A_matrix_sym_con(2,2) = va_dot_va_part;
A_matrix_sym_con(2,3) = va_dot_Vout_part;
A_matrix_sym_con(2,6) = va_dot_xv_part;

% Already linear equations:

%xa_dot = va
A_matrix_sym_con(1,2) = 1; 

%Vout_int_dot = Vout
A_matrix_sym_con(4,3) = 1; 

% I_dot = (Kpv*Vout + Kiv*Vout_int)/L - R/L *I
A_matrix_sym_con(5,3) = Kpv/L;
A_matrix_sym_con(5,4) = Kiv/L;
A_matrix_sym_con(5,5) = -R/L;

% xv_dot = vv
A_matrix_sym_con(6,7) = 1;

% vv_dot = 1/mfcv*(kfcv*(d0 + 10*A0 - xv)- Ki*I - c*vv);
A_matrix_sym_con(7,6)= - kfcv/mfcv;
A_matrix_sym_con(7,5) = - Ki/mfcv;
A_matrix_sym_con(7,7) = - c/mfcv;

% substitution of states initial conditions
A_matrix_sym_con_subs1 = subs(A_matrix_sym_con,{a0 e0 i0 OM0 om0 th0 xa0 va0 Vout0 Vout_int0 I0 xv0 vv0},{ICS});

% substitution of initial position - velocity and accelerations
A_matrix_sym_con_subs2 = subs(A_matrix_sym_con_subs1,{r0 v0 at0 an0 ah0 h0 mu_E M D},{ICS2});

% substitution of all the other variables
A_matrix_sym_con_subs3 = subs(A_matrix_sym_con_subs2,{m e_gap Vbias Cf Kpa Kda A epsilon valve_A0 valve_r0...
                     mfcv kfcv c Ki Kpv L R Kiv T2 p2 gammaXe MW RXe mXe rho T_stor P_stor eXe dV}, {[ICS3 ICS4 ICS5 ICS6]});

A_matrix_values_con = double(A_matrix_sym_con_subs3);   

end
