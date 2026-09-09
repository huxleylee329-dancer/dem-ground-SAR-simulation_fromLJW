function epi = EPI(phase_filter, phase_real)
%计算滤波算法条纹保持系数；论文61页
%输入：滤波后和真实相位矩阵
%输出：EPI
[r,c] = size(phase_filter);

sum_filter = zeros(r+1,c+1);
sum_real = zeros(r+1,c+1);
% 滤波后的梯度模和
A = zeros(r+1,c+1);
B = zeros(r+1,c+1);
C = zeros(r+1,c+1);
A(1:r,1:c) = phase_filter;%对应phi_s(i,j)
B(1:r,2:c+1) = phase_filter;%对应phi_s(i,j+1)
C(2:r+1,1:c) = phase_filter;%对应phi_s(i+1,j)
sum_filter = abs(A-B)+abs(A-C);
sum_filter = sum_filter(1:end-2,1:end-2);

A = zeros(r+1,c+1);
B = zeros(r+1,c+1);
C = zeros(r+1,c+1);
A(1:r,1:c) = phase_real;%对应phi_real(i,j)
B(1:r,2:c+1) = phase_real;%对应phi_real(i,j+1)
C(2:r+1,1:c) = phase_real;%对应phi_real(i+1,j)
sum_real = abs(A-B)+abs(A-C);
sum_real = sum_real(1:end-2,1:end-2);

epi = sum(sum(sum_filter))/sum(sum(sum_real));

end