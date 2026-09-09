function [ Sat_pos ] = ploy_regress( tt,Satpos,t_UTC )
N = length(tt);
X = [ones(N,1)  tt'  (tt').^2  (tt').^3];
Y = Satpos';
a = regress(Y,X);   % Y = X*a
Sat_pos = a(1)*ones(size(t_UTC)) + a(2)*t_UTC + a(3)*t_UTC.^2 + a(4)*t_UTC.^3;
end