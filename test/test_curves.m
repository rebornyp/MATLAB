t = 0:pi/50:10*pi;
st = sin(t);
ct = cos(t);

figure
plot3(st,ct,t)
hold on
% plot3([1 2],[1 2],[1,2])
plot3([0 1 4],[0 2 3], [0 20 0])
% plot3(2, 3, 4)
