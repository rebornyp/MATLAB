c = 3 * 10 ^ 8;
f = 2.45 * 10 ^ 9; % 2.45GHzµç´Å²¨
r1 = 4.4;
w = c / (2 * f) * ((r1 + 1)/2) ^ 0.5; %¿í¶È
W = 37.36;
h = 1.6;
deltaL = 0.75;
re = 3.73;
L = 30.21;
k = 0.001;
theta = linspace(0, pi, 180);
phi = linspace(0, 2*pi, 360);
[a, b] = meshgrid(theta, phi);
r = sin(k*h/2 .* sin(a) .* cos(b)) .* sin(k*w/2 .* cos(a)) .* sin(a) .* cos(k*L/2 .* cos(b)) ./ (k*h/2* k*w/2 .*sin(a) .*cos(b) .*cos(a));
[x, y, z] = sph2cart(b, pi/2-a, r);
% mesh(x, y, z);
surf(x, y, z);

