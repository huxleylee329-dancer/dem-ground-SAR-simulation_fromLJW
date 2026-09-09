function [f_out] = iteration_semi_Newton(iter,f0,N,ii,Data,T)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
for iii = 1:iter
    H1_f = diff_Newton_H(ii,N,f0,Data,T);
    H2_f = diff2_Newton_H(ii,N,f0,Data,T);
    f_out=f0 - real(H1_f/H2_f);
    f0=f_out;
end
end

