function [phase_removeFlat,Flat_phase_wrapped] = removeFlat_range(interf_phase2)
% 距离向去平地相位
rangefft = fft(interf_phase2.'); % 距离向FFT
% range
[nr3,nc3] = size(rangefft);
spectrumsum = sum(abs(rangefft),2);
nn = nr3;
[r,c,v] = find(spectrumsum == max(spectrumsum));
tempi = (0:nn-1);              % 行线性相位网格
fi = double(r(1)-1)/nn;
j = sqrt(-1);
ei = exp(-j*2*pi*fi);          % 红蓝边界,频率为负;蓝红边界,频率为正.
fp1 = ones(nc3,1)*ei.^tempi;
% phase_removeFlat = angle(image_master.*conj(image_slave_regis2).*fp1(1:nc3,:));
j = sqrt(-1);
phase_removeFlat = angle(exp(j*interf_phase2).*fp1(1:nc3,:));

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