function interf_phase = interferometric_phase(image1,image2)
% 获取干涉相位
interf_phase = angle(image1.*conj(image2));
% phasemap = deos(256);%new add 20140509
% figure;
% imagesc(interf_phase)
% % colormap(phasemap) %new add 20140509
% colormap(jet)
% title('interferometric phase','fontWeight','Bold')
end