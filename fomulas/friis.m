% friis电磁传播公式
P0=10;
F=915;
T=-80;
syms r;
% distance which is the key to value the 
f=P0-32.45-20*log(F)-20*log(r/1000)-T;