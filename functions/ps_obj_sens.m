function J = ps_obj_sens(x, T, ICS, options, data)

%{
 Function that returns the objective function to be used in the sensitivity
 analysis.
 
 INPUT:  1. x: changing parameters of the function
         2. T: Orbital period [s]
         3. ICS: Initial conditions of the states
         4. options: Options used for the integration
         5. data: Characteristic data of GOCE 
        
 OUTPUT: 1. J : Objective function

 FUNCTION REQUIRED: Simulation

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}

% Definition of the parameters in which to perform the sensitivity analysis
data.accelerometer.Kpa = x(1);       
data.accelerometer.Kda = x(2);       
data.valve.mfcv = x(3); 
data.valve.Ki = x(4);   
data.valve.Kpv = x(5);  
data.valve.Kiv = x(6);     
data.valve.kfcv = x(7); 
data.xenon.p2 = x(8); 
data.xenon.T2 = x(9); 
data.valve.R = x(10); 
data.valve.L = x(11); 

% Integration
cond.type = 'Nominal';
[t_pl, z] = ode15s(@simulation, [0, T], ICS , options, data, cond);

parout = zeros(length(t_pl),11);
for k = 1:length(t_pl) 
    [~, parout(k,:)] = simulation(t_pl(k), z(k, :) , data, cond);
end

Thrust = parout(:,1)*(1e+3);  % [mN]
Drag  = parout(:,2)*(1e+3);   % [mN]

% Cost function
J = 1/(2*var(Thrust+Drag))*sum((Thrust+Drag).^2);

end
