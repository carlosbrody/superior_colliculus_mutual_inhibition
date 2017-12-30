close all;
clear all;
[x,y,z] = peaks(1000);

% has positive bumps too
%figure;
%surf(x,y,z);
z2 = -abs(z);

% wells too deep
%figure;
%surf(x,y,z2)

figure;
z3 = z2./10; % this makes no difference, need nonlinear scaling
s=surf(x,y,z3,'EdgeColor', 'none');

% change camera view angle, azimuth, elevation angle. 
view(135, 75)

% change colormap, unfortunately I think default JET is best. 
%map = colormap('copper');
%colormap(flipud(map))

