function [move_r,move_c,normcor] = real_coherent(image1,image2,figshow)
% 逐点的滑动窗相关计算转化为图像块之间的互相关
figshow = 0;
%% 实相干系数计算（FFT）
[nr,nc] = size(image1); %假设image1与image2大小相同

im1abs = abs(image1);
im2abs = abs(image2);

im1fft = fft2(im1abs);
im2fft = fft2(im2abs);
cor = abs(fftshift(ifft2(im1fft.*conj(im2fft))));%先求出主辅的功率谱密度，再用ifft2进行二维傅里叶逆变换得到互相关函数，fftshift不会进行傅里叶变换，只是把fft和ifft的结果进行偏移，变成y轴对称
normcor = cor./max(max(cor));%两个max是找出二维矩阵中最大值的元素，将矩阵每个元素点除最大值，这样就能将矩阵中最大值元素转变为1
[r,c] = deal(fix(nr/2)+1,fix(nc/2)+1);
[r1,c1,~] = find(normcor == 1);%找到矩阵中元素为1的位置，从而求出峰值偏移量，也就是配准偏移量
move_r = r - r1;
move_c = c - c1;
%% 输出
if figshow == 1
    figure;
    imagesc(normcor)
    colormap(jet)
    colorbar
    title('real coherence coefficient','fontWeight','Bold')
end
end