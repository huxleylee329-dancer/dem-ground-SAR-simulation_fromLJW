function [xout,energy_ratio] = zero_signals(xin,ratio)
%对输入的信号进行零化处理。
% xin   --实信号
% ratio --稀疏度(0~1)
% Update: 2021/03/18
%   重新对ratio进行表述，如果 ratio 在0~1之间，则为比值；
%     如果ratio为整数，则为具体的非零元素个数。
N = length(xin(:));
% 根据稀疏度确定门限
if (ratio>0)&&(ratio<1)
    K = round(ratio*N);  %稀疏元素个数
else
    K = ratio;
end

if K==0
    K = 1;
elseif K>N
    K = N;
end
xin_sorted = sort(abs(xin(:)),'descend');    %全为正
tre = xin_sorted(K);
xout = xin;
xout(abs(xin)<tre) = 0;

energy_full = sum(abs(xin(:)).^2);
energy_current = sum(abs(xout(:)).^2);

energy_ratio = energy_current/energy_full;

end