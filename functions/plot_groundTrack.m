function h = plot_groundTrack(lon, lat,color)
%{
Function used to represent the Ground Track for the unperturbated case.


INPUT: 1. lon: satellite longitude (deg)
       2. lat: satellite latitude  (deg)
       3. color: color chosen for the plot 

OUTPUT: 1. Plot the Ground Track over Earth surface

FUNCTIONS REQUIRED: -

CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra

%}

lon=lon-360;
for i=1:length(lon)
    if lon(i)<-180
        lon(i)=lon(i)+360;
    end
end
tol = 100;
curve_no = 1;    % Breaks the ground track up into separate curves which start
n_curves = 1;    % and terminate at right ascensions in the range [-180,+180] (deg).
k = 0;
lon_prev = lon(1);

% Propagate lon and lat over time and insert them into LON and LAT vector:
% LON: cell array containing the longitude for each of the curves comprising the ground track plot
% LAT: cell array containing the latitudes for each of the curves comprising the ground track plot
for i = 1:length(lon)
    if abs(lon(i) - lon_prev) > tol
    curve_no = curve_no + 1;
    n_curves = n_curves + 1;
    k = 0;
    end
k = k + 1;
LON{curve_no}(k) = lon(i);
LAT{curve_no}(k) = lat(i);
lon_prev = lon(i);
end


% PLOT GROUND TRACK:
hold on
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
title('Ground Track')
%axis equal
grid minor
for i = 1:n_curves
    h=plot(LON{i}, LAT{i},color);
end
axis ([-180 180 -90 90])
line([-180 180],[0 0], 'Color','k','linewidth',0.2) %the equator

plot(lon(1),lat(1),'o','MarkerFaceColor','y','MarkerSize',9)
text(lon(1),lat(1),'  Start')
hold on
plot(lon(end),lat(end),'o','MarkerFaceColor','r','MarkerSize',9) 
text(lon(end),lat(end),'   End')

end 