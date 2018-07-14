tnum = 36;
theta=linspace(0, pi, tnum);
phi = linspace(0, 2*pi, tnum*2);
[x, y] = meshgrid(theta, phi);
r=10;
[x, y, z] = sph2cart(y, pi/2-x, r);

for i=1 : numel(x)
  if z(i) >= 8
      plot3([0,x(i)/8], [0,y(i)/8], [0,8]);
      hold on;
      plot3([x(i)/8, x(i)], [y(i)/8, y(i)], [8, 16-z(i)]);
      hold on;
  else
    plot3([0,x(i)], [0,y(i)], [0,z(i)]);
    hold on;
  end  
end
setBasicEnvironment();
axis equal;
% hold on;
% width = 20;
% num = 25;
% move = 10;
% x1=linspace(-width/2, width/2, num);
% z1=linspace(-width/2, width/2, num);
% [x1, z1] = meshgrid(x1,z1);
% y1=ones(num,num)*(width+move);
% surf(x1, y1, z1)
% hold on;
% x2=linspace(-width/2, width/2, num);
% y2=linspace(move, move+width, num);
% [x2, y2] = meshgrid(x2,y2);
% z2=ones(num,num)*width/2;
% surf(x2, y2, z2)
% hold on;
% x3=linspace(-width/2, width/2, num);
% y3=linspace(move, move+width, num);
% [x3, y3] = meshgrid(x3,y3);
% z3=ones(num,num)*-width/2;
% surf(x3, y3, z3)
% hold on;
% z4=linspace(-width/2, width/2, num);
% y4=linspace(move, move+width, num);
% [z4, y4] = meshgrid(z4,y4);
% x4=ones(num,num)*width/2;
% surf(x4, y4, z4)
% hold on;
% z5=linspace(-width/2, width/2, num);
% y5=linspace(move, move+width, num);
% [z5, y5] = meshgrid(z5,y5);
% x5=ones(num,num)*-width/2;
% surf(x5, y5, z5)
% hold on;