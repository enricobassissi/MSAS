
function [value, isterminal, direction] = event_off3(~,y,data,~)

%{ 
Event function used for the Off-Nominal condition 3 (PI failure &
recovery).

 INPUT:  1. y : state vector
         2. data : Struct with characteristic data of GOCE

 OUTPUT: 3. value: mathematical expression describing the event
         4. isterminal: put it = 1 if the integration is to terminate when the event occurs. 0 otherwise.
         5. direction:  put it = 0 if all zeros are to be located (the default). A value of +1 locates only zeros where the 
                        event function is increasing, and -1 locates only zeros where the event function is decreasing. 

 FUNCTIONS REQUIRED: -

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}

value = (10*data.valve.A0 + data.valve.d0) - y(12) ;

isterminal = 1; 

direction = - 1 ; 

end