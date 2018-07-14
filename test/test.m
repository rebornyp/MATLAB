% friis电磁传播公式
P0=10;
F=915;
T=-80;
%f=input('enter the frequency f: ');
%r=input('enter the distance s: ');
%g=solve('P0-32.45-20*log(f)-20*log(r)-T', r);
syms r;
f=P0-32.45-20*log(F)-20*log(r/1000)-T;
% g=solve(f);
% eval(g)
f1=r^2-1;
g=solve(f1);
eval(g)

% a=[1:1:20];
% f1=P0-32.45-20*log(F)-20*log(a/1000);
% f1





