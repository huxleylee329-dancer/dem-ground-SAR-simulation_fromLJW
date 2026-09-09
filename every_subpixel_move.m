function [move_a,move_r] = every_subpixel_move( m,n,a,b,c,d,e,f )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

move_a = a(1) + b(1)*m + c(1)*n + d(1)*m^2 + e(1)*n^2 + f(1)*m*n;
move_r = a(2) + b(2)*m + c(2)*n + d(2)*m^2 + e(2)*n^2 + f(2)*m*n;

end