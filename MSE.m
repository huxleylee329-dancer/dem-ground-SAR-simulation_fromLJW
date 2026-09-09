function mse = MSE(phase_filter, phase_real)
%计算平均均方误差
%输入：滤波后和真实相位矩阵
%输出：MSE,单位为相位，rad
j = sqrt(-1);
error = phase_filter - phase_real;
square_error = angle(exp(j*error)).^2;
mse = mean(mean(square_error));
end