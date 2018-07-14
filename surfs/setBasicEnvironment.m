function [ output_args ] = setBasicEnvironment()
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
width = 20;
swidth=5;
num = 25;
snum = num*0.4;
move = 10;
%最远端面
x1=linspace(-width/2, width/2, num);
z1=linspace(-width/2, width/2, num);
[x1, z1] = meshgrid(x1,z1);
y1=ones(num,num)*(width+move);
surf(x1, y1, z1)
hold on;
%沿z轴最上面
x2=linspace(-width/2, width/2, num);
y2=linspace(move, move+width, num);
[x2, y2] = meshgrid(x2,y2);
z2=ones(num,num)*width/2;
surf(x2, y2, z2)
hold on;
%沿z轴最下面
x3=linspace(-width/2, width/2, num);
y3=linspace(move, move+width, num);
[x3, y3] = meshgrid(x3,y3);
z3=ones(num,num)*-width/2;
surf(x3, y3, z3)
hold on;
%沿x轴最前面
z4=linspace(-width/2, width/2, num);
y4=linspace(move, move+width, num);
[z4, y4] = meshgrid(z4,y4);
x4=ones(num,num)*width/2;
surf(x4, y4, z4)
hold on;
%沿x轴最后面
z5=linspace(-width/2, width/2, num);
y5=linspace(move, move+width, num);
[z5, y5] = meshgrid(z5,y5);
x5=ones(num,num)*-width/2;
surf(x5, y5, z5)
hold on;
%最前面上
x6=linspace(-width/2, width/2, num);
z6=linspace(swidth/2, width/2, num);
[x6, z6] = meshgrid(x6,z6);
y6=ones(num,num)*(move);
surf(x6, y6, z6)
hold on;
%最前面下
x6=linspace(-width/2, width/2, num);
z6=linspace(-swidth/2, -width/2, num);
[x6, z6] = meshgrid(x6,z6);
y6=ones(num,num)*(move);
surf(x6, y6, z6)
hold on;
%最前面左
x6=linspace(-width/2, -swidth/2, snum);
z6=linspace(-swidth/2, swidth/2, snum);
[x6, z6] = meshgrid(x6,z6);
y6=ones(snum,snum)*(move);
surf(x6, y6, z6)
hold on;
%最前面右
x6=linspace(swidth/2,width/2, snum);
z6=linspace(-swidth/2, swidth/2, snum);
[x6, z6] = meshgrid(x6,z6);
y6=ones(snum,snum)*(move);
surf(x6, y6, z6)
hold on;

end

