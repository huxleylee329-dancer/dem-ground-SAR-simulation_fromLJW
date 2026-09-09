%% 仿真数据第一组
clc;
clear all;
close all;
load('../Input/coherence_experience1_6dB.mat')
load('../Input/wrapped_phase_experience1_6dB.mat')
load ../Input/real_phase_1.mat;
Original_Interferogram=wrapped_phase_experience_6dB(251:400,51:200);
coherence=coherence_experience_6dB(251:400,51:200);%已知的相干系数
Win_size=17;
phase_filter = Improved_Goldstein_filter_20180402(Original_Interferogram,Win_size,coherence);
%% 仿真数据第一组滤波结果评估
figure;
imagesc(phase_filter);
colormap(jet);
axis image;
colorbar;
phase_error = wrap(phase_filter-real_phase_1);
axis image;
colorbar;
MSE1 = MSE(phase_filter,real_phase_1);
EPI1 = EPI(phase_filter,real_phase_1);
[nres1,residuemat] = Calculation_Residues(phase_filter);
%% 仿真数据第二组
clc;
clear all;
load ../Input/noise_phase_slope.mat;
load ../Input/multi_terrain_slope.mat;
load ../Input/HoA1.mat;
load ../Input/real_phase_slope.mat;
Original_Interferogram=noise_phase;%原始干涉相位；noise_phase即为noise_phase_slope.mat里的变量名
coherence=0.3*ones(200,150);%没有给相干系数怎么验证呢？相干系数越高滤波效果越不好
real_phase = real_phase_slope;%真实相位
Win_size = 11;
phase_filter = Traditional_Goldstein_Correlation(Original_Interferogram,Win_size,coherence);
%% 仿真数据第二组滤波结果评估
figure;
imagesc(phase_filter);
colormap(jet);
colorbar;
MSE1 = MSE(phase_filter,real_phase);
EPI1 = EPI(phase_filter,real_phase);
[nres1,residuemat] = Calculation_Residues(phase_filter);

%% 仿真数据第三组
clc;
close all;
clear all;
load ../Input/noisy_phase_real_data.mat;
load ../Input/coherence_realdata.mat;
Original_Interferogram=phase(440-200:440+199,543-200:543+199);
coherence=coherence(440-200:440+199,543-200:543+199);
Win_size = 17;
phase_filter = Improved_Goldstein_filter_20180402(Original_Interferogram,Win_size,coherence);
%% 仿真数据第三组滤波结果评估
figure;
imagesc(phase_filter);
colormap(jet);
colorbar;
[nres3,residuemat] = Calculation_Residues(phase_filter);
PSD = Phase_Standard_Deviation(phase_filter);%相位标准差

%% 作图
a1=60;
b1=50;
c1=150;
d1=80;
a2=150;
b2=250;
c2=150;
d2=80;
figure;
imagesc(phase_filter);
colormap(jet);
colorbar;
rectangle('Position',[a1,b1,c1,d1], 'edgecolor','white', 'LineWidth',2);
rectangle('Position',[a2,b2,c2,d2], 'edgecolor','white', 'LineWidth',2);
slope_filter_cut1 = phase_filter(b1:b1+d1,a1:a1+c1);
slope_filter_cut2 = phase_filter(b2:b2+d2,a2:a2+c2);
figure;
imagesc(slope_filter_cut1);
colormap(jet);
axis image;
colorbar;
figure;
imagesc(slope_filter_cut2);
colormap(jet);
axis image;
colorbar;