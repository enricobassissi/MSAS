function [LON, LAT] = lonlat(R,t)

%{ 
Function that returns longitude and latitude values at a given position. 

 INPUT:  1. R: position vector [Km]
         2. t: time instant    [s]
        
 OUTPUT: 1. LON: longitude [deg]
         2. LAT: latitude  [deg]

 FUNCTIONS REQUIRED: -

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}

% Earth angular velocity
omegaE = (2*pi + 2*pi/365.26)/(24*3600);
thetaG0 = 0;

r = sqrt((R(1).^2)+(R(2).^2)+(R(3).^2));
l = R(1)./r ; m = R(2)./r; n = R(3)./r;

% Declination
Dec(:,1) = asin(n(:,1));

% Right ascension
for g=1:length(r)
    if m(g,1)>0
        RA(g,1)=acos((l(g,1))/cos(Dec(g,1)));
    else
        RA(g,1)=2*pi-acos((l(g,1))/cos(Dec(g,1)));
    end
end

thetaG  = thetaG0 + omegaE*t;

% Longitude and Latitude
lambda = RA - thetaG; 
lambda = wrapToPi(real(lambda));

lambda = rad2deg(lambda);
Dec = rad2deg(Dec);

LON = lambda;
LAT = Dec;

end