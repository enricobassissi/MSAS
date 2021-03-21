function evolutions_video(t_pl, z, data)

%{
 Function that returns the videos of the moving valve pin position and of
 the accelerometer proof mass position
 
 INPUT:  1. t_pl: integration time vector [s]
         2. z: state vectors evolutions in time
         4. data: GOCE data
        
 OUTPUT: -

 FUNCTION REQUIRED: -

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra
%}

% Valve section
close all
hold off
set(gca,'nextplot','replacechildren');
video = VideoWriter('Video_valve_AVI', 'Uncompressed AVI');
%video.Quality = 100;
open(video);  


n = length(t_pl);
l = yline([(10*data.valve.A0 + data.valve.d0)*1e3], 'linewidth', 2, 'color', 'k');
hold on
xlim([-0.1 0.1])

for k=1:30:n+1
    k=k-1;
    
    if k~=0   
    
     v = fill([-data.valve.r0*1e3, data.valve.r0*1e3, data.valve.r0*1e3, -data.valve.r0*1e3],...
         [3.667, 3.667, z(k,12)*1e3, z(k,12)*1e3], 'r'); %[mm]
     v.FaceAlpha=0.6;
    
    drawnow
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    ylabel('$x_v$ [mm]')   

    moment_view = t_pl(k);
    title (sprintf('Valve displacement', fix(moment_view/60)));
    frame = getframe(gcf);
    writeVideo(video, frame);
    
    if k~=1 
    delete(v);
    end 
    
    end
end
close
close(video)

% Accelerometer

hold off 
set(gca,'nextplot','replacechildren');
video = VideoWriter('Video_accelerometer_AVI', 'Uncompressed AVI');
%video.Quality = 100;
open(video);  

b = sqrt(data.accelerometer.A);
n = length(t_pl);
l = yline(0, '-.k', 'linewidth', 2);
hold on
grid on
xlim([-0.67*b*1e-7 0.67*b*1e-7]);
ylim([-1e-8 1e-8]);

for r=1:30:n+1
    r=r-1;
    
    if r~=0   

    acc = plot([-0.5*b*1e-7, 0.5*b*1e-7], [z(r,7)*1e3, z(r,7)*1e3], 'r', 'linewidth', 3);
 
    drawnow
    ylabel('$x_a$ [mm]')   
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])

    moment_view = t_pl(r);
    title (sprintf('Proof mass displacement', fix(moment_view/60)));
    frame = getframe(gcf);
    writeVideo(video, frame);
    
    if r~=1 
    delete(acc);
    end 
    
    end
end
close
close(video)

end

