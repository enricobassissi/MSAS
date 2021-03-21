function [] = plot_diagnostic_off_report(z,t_pl,parout, z_off1,t_pl_off1,parout_off1,z_off2,t_pl_off2,parout_off2,z_off3,t_pl_off3,parout_off3, data,colors)

%{ 
Function that returns the plot of Drag vs Thrust and Valve Displacement for the following cases:
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
xv = z(:,12)*(1e+3);       % [mm]
Th = parout(:,1)*(1e+3);   % [mN]
D  = parout(:,2)*(1e+3);   % [mN]
t_pl_h = t_pl/3600;        % [h]

% Recovering off-nominal conditions and plots
name.cases = {'\textbf{Grid-erosion}','\textbf{Xe-Leakage}','\textbf{PI failure \& recovery}'};

% Plot thrust vs drag
figure('Name','T+D and xv evolution, design for latex report')
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
    xv_off = z_off(:,12)*(1e+3);      %[mm]
    t_pl_h_off = t_pl_off/3600;       %[h]
    
    subplot(2,3,i)
    p1 = plot(t_pl_h, Th);
    hold on
    p2 = plot(t_pl_h_off, Th_off,'Color',colors(8,:));
    p3 = plot(t_pl_h, abs(D),'Color',colors(3,:)) ;
    if i == 3
        ylim([-0.2 4])
    else
        ylim([-0.2 2.3])
    end
    hold on, 
    xlabel('$t$ [h]'), ylabel('Forces [mN]')
    title(sprintf('%s', name.cases{i}),'FontSize',18)
    
    subplot(2,3,3+i)
    p4 = plot(t_pl_h, xv);
    hold on
    p5 = plot(t_pl_h_off, xv_off,'Color',colors(8,:));
    yline([(10*data.valve.A0 + data.valve.d0)*(1e+3)]);
    if i == 3
        ylim([3.6674 3.6683])
    else
        ylim([3.6670 3.6683])
    end
    xlabel('$t$ [h]'); ylabel('$x_v$ [mm]'); 
    if i == 3
        hL = legend([p1 p2 p3],{'Thrust Nom.','Thrust Off-Nom.', '$|$Drag$|$'},'Location','Northeast');
        hL2 = legend([p4 p5],{'$x_v$ Nom.','$x_v$ Off-Nom.'},'Location','Southeast');
    end        
end

end
