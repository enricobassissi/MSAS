function [] = plot_diagnostic_off(z,t_pl,parout, z_off1,t_pl_off1,parout_off1,z_off2,t_pl_off2,parout_off2,z_off3,t_pl_off3,parout_off3, data,colors)

%{ 
Function that returns the following plots: (1) Thrust vs Drag,(2) Altitude, (3) Proof Mass Displacement, (4) Proof mass Velocity
(5) Valve control voltage, (6) Integral of valve control voltage, (7)Soleonid current, (8) Valve Spool Displacement, 
(9) Valve Spool Velocity.

In the following cases:
 - Nominal + Off-Nominal (Grid erosion)
 - Nominal + Off-Nominal (Xe-Leakage)
 - Nominal + Off-Nominal (PI failure & recovery)
 

 INPUT:  1.  z : state variables evolution in time (Nominal)
         2.  t_pl : integration times (Nominal) 
         3.  parout:  additional quantities from integration (Nominal) [Th D A_valve lon lat height r v at an ah] 
         4.  z_off1 : state variables evolution in time (Off-Nominal 1)
         5.  t_pl_off1 : integration times (Off-Nominal 1)
         6.  parout_off1:  additional quantities from integration (Off-Nominal 1) [Th D A_valve lon lat height r v at an ah]
         7.  z_off2 : state variables evolution in time (Off-Nominal 2)
         8.  t_pl_off2 : integration times (Off-Nominal 2)
         9.  parout_off2:  additional quantities from integration (Off-Nominal 2) [Th D A_valve lon lat height r v at an ah]
         10. z_off3 : state variables evolution in time (Off-Nominal 3)
         11. t_pl_off3 : integration times (Off-Nominal 3)
         12. parout_off3:  additional quantities from integration (Off-Nominal 3) [Th D A_valve lon lat height r v at an ah]
         13. data: Struct with characteristic data of GOCE
         14. colors: Colors array 
        
 OUTPUT: -

 FUNCTIONS REQUIRED: -

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}


% Recovering quantities nominal condition
xa = z(:,7)*(1e+3);       % [mm]
va = z(:,8)*(1e+3);       % [mm/s]
Vout = z(:,9);            % [V]
Vout_int = z(:,10);       % [V/s]
I = z(:,11)*(1e+3);       % [mA]
xv = z(:,12)*(1e+3);      % [mm]
vv = z(:,13)*(1e+3);      % [mm/s]
Th = parout(:,1)*(1e+3);  % [mN]
D  = parout(:,2)*(1e+3);  % [mN]
height = parout(:,6);     % [km]
t_pl_h = t_pl/3600;       %[h]

% Recovering off-nominal conditions and plots
name.cases = {'\textbf{Grid-erosion}','\textbf{Xe-Leakage}','\textbf{PI failure \& recovery}'};

% Plot thrust vs drag
figure('Name', 'Thrust vs Drag')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    Th_off = parout_off(:,1)*(1e+3);  % [mN]
    t_pl_h_off = t_pl_off/3600;       % [h]
    
    subplot(2,2,i)
    plot(t_pl_h, Th)
    hold on
    plot(t_pl_h_off, Th_off,'Color',colors(8,:)) 
    plot(t_pl_h, abs(D),'Color',colors(3,:)) 
    %ylim([-1 4])
    hold on, grid on;
    xlabel('$t$ [h]'), ylabel('Forces [mN]')
    title(sprintf('%s', name.cases{i}))
end
sgtitle('\textbf{Thrust vs Drag}')
hL = legend('Thrust Nominal','Thrust off-nominal', '$|$Drag$|$');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)

% Plot height
figure('Name', 'Altitude')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    height_off = parout_off(:,6);     % [km]
    t_pl_h_off = t_pl_off/3600;       % [h]
    
    subplot(2,2,i)
    plot(t_pl_h,height)
    hold on
    plot(t_pl_h_off, height_off,'Color',colors(8,:))
    ylim([250 330])
    grid on
    xlabel('$t$ [h]'); ylabel('Altitude [km]');
    title(sprintf('%s', name.cases{i}))
end   
sgtitle('Altitude')
hL = legend('Nominal','Off-Nominal');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)

% Plot xa 
figure('Name', 'Proof mass displacement')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    xa_off = z_off(:,7)*(1e+3);  %[mm]
    t_pl_h_off = t_pl_off/3600;  %[h]
    
    subplot(2,2,i)
    plot(t_pl_h, xa)
    hold on
    plot(t_pl_h_off,xa_off,'Color',colors(8,:))
    %ylim([-11*1e-9 6*1e-9])
    xlabel('$t$ [h]'); ylabel('$x_a$ [mm]'); grid on;
    title(sprintf('%s', name.cases{i}))
end  
sgtitle('Proof mass displacement')
hL = legend('Nominal','Off-Nominal');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)


% Plot va
figure('Name', 'Proof mass velocity')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    va_off = z_off(:,8)*(1e+3);  % [mm]
    t_pl_h_off = t_pl_off/3600;  % [h]
    
    subplot(2,2,i)
    plot(t_pl_h(183:end), va(183:end))
    hold on
    plot(t_pl_h_off(183:end),va_off(183:end),'Color',colors(8,:))
    ylim([-10*1e-12 10*1e-12])
    xlabel('$t$ [h]'); ylabel('$v_a$ [mm/s]'); grid on;
    title(sprintf('%s', name.cases{i}))
end  
sgtitle('Proof mass velocity')
hL = legend('Nominal','Off-Nominal');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)


% Plot Vout 
figure('Name', 'Valve control voltage')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    Vout_off = z_off(:,9);       % [V]
    t_pl_h_off = t_pl_off/3600;  % [h]
    
   subplot(2,2,i)
   plot(t_pl_h, Vout)
   hold on
   plot(t_pl_h_off,Vout_off,'Color',colors(8,:))
   xlabel('$t$ [h]'); ylabel('$V_{out}$ [V]'); grid on;
   title(sprintf('%s', name.cases{i}))
end   
sgtitle('Valve control voltage')
hL = legend('Nominal','Off-Nominal');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)

% Plot Vout_int 
figure('Name', 'Integral of valve control voltage')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    Vout_int_off = z_off(:,10);  % [Vs]
    t_pl_h_off = t_pl_off/3600;  % [h]
    
   subplot(2,2,i)
   plot(t_pl_h, Vout_int)
   hold on
   plot(t_pl_h_off,Vout_int_off,'Color',colors(8,:))
   xlabel('$t$ [h]'); ylabel('$\int{Vout}$ [Vs]'); grid on; 
   title(sprintf('%s', name.cases{i}))
end   
sgtitle('Integral of valve control voltage')
hL = legend('Nominal','Off-Nominal');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)


% Plot Current 
figure('Name', 'Solenoid Current')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    I_off = z_off(:,11)*(1e+3);  % [mA]
    t_pl_h_off = t_pl_off/3600;  % [h]
    
    subplot(2,2,i)
    plot(t_pl_h, I)
    hold on
    plot(t_pl_h_off,I_off,'Color',colors(8,:))
    xlabel('$t$ [h]'); ylabel('$I$ [mA]'); grid on;
    title(sprintf('%s', name.cases{i}))
end      
sgtitle('Solenoid current') 
hL = legend('Nominal','Off-Nominal');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)


% Plot xv 
figure('Name', 'Valve spool displacement')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    xv_off = z_off(:,12)*(1e+3); %[mm]
    t_pl_h_off = t_pl_off/3600;  %[h]
    
    subplot(2,2,i)
    plot(t_pl_h, xv)
    hold on
    plot(t_pl_h_off, xv_off,'Color',colors(8,:))
    yline([(10*data.valve.A0 + data.valve.d0)*(1e+3)]);
    %ylim([3.6676 3.6683])
    xlabel('$t$ [h]'); ylabel('$x_v$ [mm]'); grid on;
    title(sprintf('%s', name.cases{i}))
end       
sgtitle('Valve spool displacement')
hL = legend('Nominal','Off-Nominal');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)


% Plot vv 
figure('Name', 'Valve spool velocity')
for i = 1:3
    z_off = [];  t_pl_off = []; parout_off = [];
    if i == 1
        z_off = z_off1; t_pl_off = t_pl_off1; parout_off = parout_off1;
    elseif i == 2
        z_off = z_off2; t_pl_off = t_pl_off2; parout_off = parout_off2;
    else
        z_off = z_off3; t_pl_off = t_pl_off3; parout_off = parout_off3;
    end
    
    vv_off = z_off(:,13)*(1e+3); % [mm/s]
    t_pl_h_off = t_pl_off/3600;  % [h]
    
    subplot(2,2,i)
    plot(t_pl_h, vv)
    hold on
    plot(t_pl_h_off, vv_off,'Color',colors(8,:))
    xlabel('$t$ [h]'); ylabel('$v_v$ [mm/s]'); grid on;
    title(sprintf('%s', name.cases{i}))
end       
sgtitle('Valve spool velocity')
hL = legend('Nominal','Off-Nominal');
newPosition = [0.65 0.2 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position',newPosition,'Units',newUnits)
   
    
end

