
x=linspace(-10, 10, 100);
y=linspace(-10, 10, 100);
[x, y] = meshgrid(x,y);
z=ones(100,100)*3;
surf(x, y, z)
% mesh(x, y)
%  plot3(x,y,z);
% mesh(x,y,z);
% surf(x, y, z);