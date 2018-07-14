theta=linspace(0, pi, 36);
phi = linspace(0, 2*pi, 72);
[x, y] = meshgrid(theta, phi);
r=100;
[x, y, z] = sph2cart(y, pi/2-x, r);

for i=1 : numel(x)
%   fprintf(1, '%f  %f %f \n', x(i), y(i), z(i));
  if x(i) < y(i)
    plot3([0,x(i)], [0,y(i)], [0,z(i)]);
  end
  hold on;
end

% surf(x, y, z);
% hold on;
% x.size();



