%% Plot Multiple Spheres
% Define |x|, |y|, and |z| as coordinates of a sphere. 

% Copyright 2015 The MathWorks, Inc.


[x,y,z] = sphere;

%%
% Plot a sphere centered at the origin. Plot two more spheres centered at
% |(3,-2,0)| and |(0,1,-3)|.

figure
surf(x,y,z)

hold on 
surf(x+3,y-2,z) % centered at (3,-2,0) 
surf(x,y+1,z-3) % centered at (0,1,-3)


