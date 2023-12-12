
clear all
%unit in degree

%The degree between adjacent points on the sphere
step = 2;

grain = [];

n_polar = 360 / step -1;
n_azimuthal = 180 / step;

%Polarangle goes around a circle
for k = 0:n_polar
    disp(k)
    polar = k * step;
    
    %For each phi, calculate all the theta
    
    %When phi = 0, we need to do it seperately to account for the south and
    %north pole
    if k == 0
        for jj = 0:n_azimuthal

            azimuthal = jj*step;

            xyz = sphere2cartesian(polar,azimuthal);
            grain = [grain; xyz];
        end
    else
        
        for jj = 1:n_azimuthal-1

            azimuthal = jj*step;

            xyz = sphere2cartesian(polar,azimuthal);
            grain = [grain; xyz];
        end
    end
    
    
end



save('grain(2).mat','grain')




%A 3-d scattering plot to visualize the grain
scatter3(grain(:,1),grain(:,2),grain(:,3))
zlim([-1.2 1.2])
ylim([-1.2 1.2])
xlim([-1.2 1.2])
xlabel('X')
ylabel('Y')
zlabel('Z')



function xyz = sphere2cartesian(polar,azimuthal)

    x = sind(azimuthal) * cosd(polar);
    
    y = sind(azimuthal) * sind(polar);
    
    z = cosd(azimuthal);
    
    xyz = [x y z];
end