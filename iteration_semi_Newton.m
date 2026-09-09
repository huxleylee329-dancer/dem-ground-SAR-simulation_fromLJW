function [f_out] = iteration_semi_Newton(iter,f0,N,ii,Data,T)
%利用牛顿迭代法估计频率
% tmp=zeros(1,iter);
for iii = 1:iter
    H1_f = diff_Newton_H(ii,N,f0,Data,T);
    H2_f = diff2_Newton_H(ii,N,f0,Data,T);
    f_out=f0 - real(H1_f/H2_f);
%     tmp(1,iii)=f_out;
    f0=f_out;
end
% plot(tmp);
end

