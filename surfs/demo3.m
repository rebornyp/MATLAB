tnum = 36;
theta=linspace(0, pi, tnum);
phi = linspace(0, 2*pi, tnum*2);
[x, y] = meshgrid(theta, phi);
r=50;
[x, y, z] = sph2cart(y, pi/2-x, r);

move = 10;
swidth = 5;
width=20;

for i=1 : numel(x)
    k = y(i) / move;
    
  if k>=1 && x(i)/k > -swidth/2 && x(i)/k < swidth/2 && z(i)/k > -swidth/2 && z(i)/k < swidth/2
      plot3([0,x(i)], [0,y(i)], [0,z(i)]);
      hold on;
  elseif k>=1 && x(i)/k > -width/2 && x(i)/k < width/2 && z(i)/k > -width/2 && z(i)/k < width/2
      plot3([0,x(i)/k], [0,y(i)/k], [0,z(i)/k]);
      plot3([x(i)/k, x(i)], [y(i)/k, 2*move-y(i)], [z(i)/k, z(i)]);
      hold on;
  elseif y(i) > 0
      plot3([0,x(i)], [0,y(i)], [0,z(i)]);
      hold on;
  end  
end
% 绘制周围三维实际环境面域
setBasicEnvironment();
axis equal;