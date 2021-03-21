function [h_earth] = plot_earth

%{
Function used to plot 3D Earth

INPUT: -

OUTPUT: 1. h_earth: 3D Earth plot

FUNCTIONS REQUIRED: -

CONTRIBUTORS:  Bassissi Enrico
               Colombo Alessandro
               De Luca Maria Alessandra

%}


C = imread('EarthTexture.jpg'); 
Re = astroConstants(23);
a_earth = 6378.137;  % equatorial radius (km)
b_earth = 6356.7523; % polar radius (km)
theta=0; 
cla
[x, y, z] = ellipsoid(0,0,0, a_earth,a_earth,b_earth,1E2); 
h_earth = surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]), 'FaceColor', 'texturemap','EdgeColor','none');
hold on;
clouds = imread('cloudCombined.jpg');
im = image(clouds);

im.AlphaData = max(clouds,[],3);    % set transparency to maximum cloud value


sCloud = surf(x*1.02,y*1.02,z*1.02,circshift(flip(clouds),[0,ceil(size(clouds,2)/360*theta)]));


sCloud.FaceColor = 'texturemap';        % set color to texture mapping
sCloud.EdgeColor = 'none';              % remove surface edge color
sCloud.CData = clouds;                  % set color data 

sCloud.FaceAlpha = 'texturemap';        % set transparency to texture mapping
sCloud.AlphaData = max(clouds,[],3);    % set transparency data 

hold off

daspect([1 1 1])                        % set aspect ratio

end