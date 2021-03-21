function [] = plot_diagnostic(z, parout, t_pl, data, colors)

%{ 
Function that returns the following plots in the Nominal case: (1) Thrust vs Drag, (2) Keplerian elements,(3) Altitude, (4) Proof Mass
%Displacement (with zoom out), (5) Proof mass Velocity (with zoom out) (6) Valve control voltage, (7) Integral of valve control voltage, (8)Soleonid current, (9) Valve Spool Displacement, 
(10) Valve Spool Velocity, (11) Non keplerian state evolution design for latex report

 INPUT:  1.  z : state variables evolution in time 
         2.  parout:  additional quantities from integration [Th D A_valve lon lat height r v at an ah] 
         3.  t_pl : integration times 
         4. data: Struct with characteristic data of GOCE
         5. colors: Colors array 
        
 OUTPUT: -

 FUNCTIONS REQUIRED: -

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}

% Recovering quantities of interest
xa = z(:,7)*(1e+3);  %[mm]
va = z(:,8)*(1e+3);  %[mm/s] 
Vout = z(:,9);       %[V]
Vout_int = z(:,10);  %[Vs]
I = z(:,11)*(1e+3);  %[mA]
xv = z(:,12)*(1e+3); %[mm]
vv = z(:,13)*1e3;    %[mm/s]

Th = parout(:,1)*(1e+3);  % [mN]
D  = parout(:,2)*(1e+3);  % [mN]
A_valve = parout(:,3);    % [m^2]
lon = parout(:,4);        % [deg]
lat = parout(:,5);        % [deg]
height = parout(:,6);     %[km]

t_pl_h = t_pl/3600;  %[h]

% -------------------------------SINGLE PLOTS --------------------------------- %

% Plot thrust vs drag
figure('Name', 'Thrust vs Drag')
plot(t_pl_h, Th)
hold on
plot( t_pl_h, abs(D), 'color', colors(8,:))
legend('Thrust', '$|$Drag$|$');
hold on, grid on;
xlabel('$t$ [h]'), ylabel('Forces [mN]')
title('Thrust vs Drag')

% Plot keplerian elements evolution
figure('Name', 'Keplerian elements')
subplot(2,3,1)
plot(t_pl_h, z(:, 1))
grid on, xlabel('$t$ [h]'), ylabel('$a$ [km]')
subplot(2,3,2)
plot(t_pl_h, z(:, 2))
grid on, xlabel('$t$ [h]'), ylabel('$e$ [-]')
subplot(2,3,3)
plot(t_pl_h, rad2deg(z(:, 3)))
grid on, xlabel('$t$ [h]'), ylabel('$i$ [deg]')
ytickformat('%,.2f');
subplot(2,3,4)
plot(t_pl_h, rad2deg(z(:, 4)))
grid on, xlabel('$t$ [h]'), ylabel('$\Omega$ [deg]')
subplot(2,3,5)
plot(t_pl_h, rad2deg(z(:, 5)))
grid on, xlabel('$t$ [h]'), ylabel('$\omega$ [deg]')
subplot(2,3,6)
plot(t_pl_h, rad2deg(wrapTo2Pi(z(:, 6))))
grid on, xlabel('$t$ [h]'), ylabel('$\theta$ [deg]')
sgtitle('Keplerian elements','FontSize',18') 


% Plot height
figure('Name', 'Altitude')
plot(t_pl_h, height)
ylim([250 330])
xlabel('$t$ [h]'); ylabel('Altitude [km]'); grid on;
title('Altitude')

% Plot xa - check the mass doesn't overcome the gap
figure('Name', 'Proof mass displacement')
plot(t_pl_h,xa)
xlabel('$t$ [h]'); ylabel('$x_a$ [mm]'); grid on;
title('Proof mass displacement (with zoom out)')
axes('position',[.63 .175 .25 .25]); % create a new pair of axes inside current figure
box on 
yline(data.accelerometer.gap*(1e+3));
yline(-data.accelerometer.gap*(1e+3));
hold on;
plot(t_pl_h, xa)
ylim([-0.6 0.6])
grid on

% Plot va
figure('Name', 'Proof mass velocity')
plot(t_pl_h(183:end),va(183:end))
ylim([-6*1e-12 5*1e-12])
xlabel('$t$ [h]'); ylabel('$v_a$ [mm/s]'); grid on;
title('Proof mass velocity (with zoom out)')
axes('position',[.63 .175 .25 .25]); % create a new pair of axes inside current figure
box on 
hold on;
plot(t_pl_h, va)
ylim([-7*1e-8 1*1e-8])
grid on

% Plot Vout 
figure('Name', 'Valve control voltage')
plot(t_pl_h, Vout);
xlabel('$t$ [h]'); ylabel('$V_{out}$ [V]'); grid on;
title('Valve control voltage')

% Plot Vout_int 
figure('Name', 'Integral of valve control voltage')
plot(t_pl_h, Vout_int);
xlabel('$t$ [h]'); ylabel('$\int{V_{out}}$ [Vs]'); grid on; %%% unità di misura?
title('Integral of valve control voltage')


% Plot Current 
figure('Name', 'Solenoid Current')
plot(t_pl_h, I)
xlabel('$t$ [h]'); ylabel('$I$ [mA]'); grid on;
ylim([0 16])
title('Solenoid current')

% Plot xv 
figure('Name', 'Valve spool displacement')
yline((10*data.valve.A0 + data.valve.d0)*(1e+3)); 
hold on
plot(t_pl_h, xv)
xlabel('$t$ [h]'); ylabel('$x_v$ [mm]'); grid on;
title('Valve spool displacement')


% Plot vv
figure('Name', 'Valve spool velocity')
plot(t_pl_h, vv)
xlabel('$t$ [h]'); ylabel('$v_v$ [mm/s]'); grid on;
title('Valve spool velocity')

% -------------------------------REPORT PLOTS --------------------------------- %

figure('Name', 'Non keplerian states evolution, design for latex report')
subplot(2,3,1)
yyaxis left
plot(t_pl_h,xa)
xlabel('$t$ [h]'); ylabel('$x_a$ [mm]'); grid on;

yyaxis right
plot(t_pl_h(183:end),va(183:end),'color',colors(8,:))
ylim([-2*1e-12 5*1e-12])
xlabel('$t$ [h]'); ylabel('$v_a$ [mm/s]'); grid on;
haxes = gca;
set(haxes, 'YColor',colors(8,:));
subplot(2,3,2)
plot(t_pl_h, Vout);
xlabel('$t$ [h]'); ylabel('$V_{out}$ [V]'); grid on;

subplot(2,3,3)
plot(t_pl_h, I)
xlabel('$t$ [h]'); ylabel('$I$ [mA]'); grid on;
ylim([0 16])

subplot(2,3,4)
yyaxis left
yline((10*data.valve.A0 + data.valve.d0)*(1e+3),'color','k','linewidth',2); 
hold on
plot(t_pl_h, xv)
xlabel('$t$ [h]'); ylabel('$x_v$ [mm]'); grid on;

yyaxis right
plot(t_pl_h, vv,'color',colors(8,:))
xlabel('$t$ [h]'); ylabel('$v_v$ [mm/s]'); grid on;
haxes = gca;
set(haxes, 'YColor',colors(8,:));

subplot(2,3,5:6)
plot(t_pl_h, Th);
hold on; grid on;
plot(t_pl_h, abs(D),'color',colors(8,:))
legend('Thrust', '$|$Drag$|$');
xlabel('$t$ [h]'), ylabel('Forces [mN]')

end