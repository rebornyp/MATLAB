theta=linspace(0, pi, 180);
phi = linspace(0, 2*pi, 360);
% [x,y] = pol2cart(theta, 20);
% line([0,0,0], [2,3,56]);
% line([2,3,56],[-10,330,-23]);
% line([0,0,0],[1, 1, 1]);
t=linspace(0, 10, 10);
x=t;y=t;z=t;
% line(x,y,z)
syms x y
z=exp(x);
fplot(x, z);
% k=2;
% b=23;
% ezsurf(k*x+b-y);
% grid on