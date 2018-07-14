%% Plot 3-D Helix
% Define |t| as values between 0 and $10\pi$. Define |st| and |ct| as
% vectors of sine and cosine values. Plot a 3-D helix.

% Copyright 2015 The MathWorks, Inc.


t = 0:pi/50:10*pi;
st = sin(t);
ct = cos(t);

figure
plot3(st,ct,t)

