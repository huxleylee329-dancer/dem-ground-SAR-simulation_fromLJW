function image_quantify = image_quantify(image)
% Input:  image:           复数图像
% Outputs:image_intensity: 强度灰度图
%         image_power:     功率灰度图
% 功率量化-----较优
Intensity = abs(image);
image_power = 20*log10(Intensity+1e-6);
Max=max(max(image_power));
Min=Max-45;%显示动态范围40dB
image_quantify=255/(Max-Min)*(image_power-Min).*(image_power>Min);
end