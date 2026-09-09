function [phase_removeFlat,Flat_phase_wrapped] = removeFlat_range_20161126(interf_phase2)
% 距离向去平地相位
rangefft = fft(interf_phase2.'); % 距离向FFT，对距离向的每一列数据进行FFT
% range
[nr3,nc3] = size(rangefft);
spectrumsum = sum(abs(rangefft),2);%将所有列求和，找到所有行中最大的峰值频率，也就是距离向的主要频率
nn = nr3;
[r,c,v] = find(spectrumsum == max(spectrumsum));%下面是为了补偿解缠后的平地相位，r代表最大峰值频率所在的行数（注意有2个，因为频谱图关于y轴对称）
tempi = (0:nn-1);              % 行线性相位网格，因为从第1列到第n列，斜距越来越长，平地相位成线性变化，因为FFT时有转置，所以变成了行线性相位网络
fi = double(r(1)-1)/nn;%文献【高分辨率合成孔径雷达干涉测量.....】中的频移法，r(1)-1代表偏移量l，nn代表N
j = sqrt(-1);
ei = exp(-j*2*pi*fi);          % 红蓝边界,频率为负;蓝红边界,频率为正.（负号就相当于把平地相位给补偿抵消掉）
fp1 = ones(nc3,1)*ei.^tempi;%tempi = (0:nn-1)代表从1到n行，fp1代表时域相位修正因子，也就是平地相位
% phase_removeFlat = angle(image_master.*conj(image_slave_regis2).*fp1(1:nc3,:));
j = sqrt(-1);
phase_removeFlat = angle(exp(j*interf_phase2).*(-fp1(1:nc3,:)));%对平地相位补偿之后，phase_removeFlat代表去平地后的相位
% figure;
% imagesc(phase_removeFlat)
% colormap(jet)
% title('interferometric phase (remove flat phase)','fontWeight','Bold')
% 去平地后距离谱
rangefft = fftshift(fft(phase_removeFlat.'));
spectrumsum = sum(abs(rangefft),2);
% 缠绕的平地相位
Flat_phase_wrapped = angle(fp1);
% 输出
% if display == 1
%     figure;
%     plot(spectrumsum) %去平地后距离谱
%     title('range spectrum after remove flat phase','fontWeight','Bold')
%     figure;
%     imagesc(Flat_phase_wrapped)%缠绕的距离向平地相位
%     colormap(jet)
%     colorbar
%     title('wrapped flat phase','fontWeight','Bold')
% end
end