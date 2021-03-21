function movie_orbit(t_pl, T, z, data, colors)

%{
 Function that returns the video of the moving 3D orbit plot
 
 INPUT:  1. t_pl: integration time vector [s]
         2. T: Orbital period [s]
         3. z: state vectors evolutions in time
         4. data: GOCE data
         5. colors: colors array RGB
        
 OUTPUT: -

 FUNCTION REQUIRED: sv_from_coe

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}


% Moving orbit
hold off
close all
coe_plot = [z(:,1) z(:,2) z(:,3) z(:,4) wrapTo2Pi(z(:,5)) wrapTo2Pi(z(:,6))];  
y_plot = zeros(size(coe_plot));
for k = 1:length(coe_plot)
    [r_plot,v_plot] = sv_from_coe(coe_plot(k,:),data.earth.mu_E);
    y_plot(k,:) = [r_plot', v_plot'];
end

lim_max=2*z(1,1); %set the 'point of view'

set(gca,'nextplot','replacechildren');
video = VideoWriter('Video_AVI','Uncompressed AVI');
% video.Quality = 100;
open(video);  
axis equal
grid on
omegaE = (2*pi + 2*pi/365.26)/(24*3600); %rad/s
n = length(t_pl);

t_magic = T;
orbit_numb = 1;
for j=1:100:n+1
    j=j-1;
    
    if j~=0   
    globe = plot_earth;
     theta = omegaE*t_pl(j);
    rotate(globe, [0 0 1],rad2deg(theta))

    hold on

    p2 = plot3(y_plot(j,1),y_plot(j,2),y_plot(j,3),'.','markersize',15,'color',colors(2,:));
    
    view([lim_max lim_max lim_max]) 
    drawnow
    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')    

    moment_view = t_pl(j);
    
    t_single_pass = t_pl(j) - t_magic;
    if t_single_pass > 0
        orbit_numb = orbit_numb + 1;
        t_magic = t_magic + T;
    end
    
    title (sprintf('T + %d [min] - Orbit No %d', fix(moment_view/60), orbit_numb));
    frame = getframe(gcf);
    writeVideo(video,frame);
    end
end
close
close(video)
end

