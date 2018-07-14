theta=linspace(0, pi, 36);
phi = linspace(0, 2*pi, 72);
[x, y] = meshgrid(theta, phi);
r=100;
[x, y, z] = sph2cart(y, pi/2-x, r);

for i=1 : numel(x)
%   if x(i) < y(i)
  if z(i) >= 80
      plot3([0,x(i)/80], [0,y(i)/80], [0,80]);
      hold on;
      plot3([x(i)/80, x(i)], [y(i)/80, y(i)], [80, 160-z(i)]);
      hold on;
  else      
    plot3([0,x(i)], [0,y(i)], [0,z(i)]);
    hold on;
  end  
end

px=linspace(-100, 100, 100);
py=linspace(-100, 100, 100);
[px, py] = meshgrid(px,py);
pz=ones(100,100)*80;
surf(px, py, pz)