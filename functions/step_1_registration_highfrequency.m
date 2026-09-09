%% 高频配准
warning off;
outfolder = '.\Kernel\Output1_registration\';%输出路径
infolder = '.\Kernel\Deflat_output\high\';
%% 载入高频单视复图像
display('载入单视复图像......');
load([infolder 'deflat_mast_high.mat']);
mast = deflat_mast_high;
load([infolder 'deflat_slave_high.mat']);
slave = deflat_slave_high;
[nr1, nc1] = size(mast);
display('载入成功！');
%% 粗配准
display('进行粗配准......');
[mast_regis1,slave_regis1,move_r,move_c] = registration_pixel(mast, slave, nr1, nc1);
display('粗配准完成！');
%% 亚像素级配准
display('进行精配准......');
[nr2, nc2] = size(mast_regis1);
[~,slave_regis2] = regis_subpixel(mast_regis1,slave_regis1,nr2,nc2);
display('精配准完成！');
%% 产生干涉相位
display('产生干涉相位......');
wrapped_phase = interferometric_phase(mast_regis1,slave_regis2);
display('完成！');
display('计算相干图......');
coherence_high = Calculation_Coherence_Coefficient_complex_no_phase(mast_regis1,slave_regis2);
display('完成！');

fig1 = figure;
set(fig1, 'visible', 'off');
imagesc(wrapped_phase);
colorbar;

saveas(fig1, [outfolder 'wrapped_phase_high.jpg']);
save([outfolder,'wrapped_phase_high.mat'],'wrapped_phase');

fig2 = figure;
set(fig2, 'visible', 'off');
imshow(coherence_high);
colorbar;
saveas(fig2, [outfolder 'coherence_high.jpg']);
save([outfolder 'coherence_high.mat'], 'coherence_high');
