close all;
clear all;
MAKE_DENSE = 2;

[x,y,z] = peaks(1000);

% has positive bumps too
%figure;
%surf(x,y,z);
if MAKE_DENSE == 1
    z2 = -abs(z);
elseif MAKE_DENSE == 2
    z2 = z;
    z2(z2<0) = 0;
    z2 = -z2;
elseif MAKE_DENSE ==3
    z2 = -z;
    z2(z2<0) = 0;
    z2 = -z2;
end


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

