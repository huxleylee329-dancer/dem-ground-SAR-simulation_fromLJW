%% 长短基线融合
clear
folder=['..\output\fusion\'];
load('../output/long/B_long.mat')
load('../output/short/B_short.mat')
load('../output/long/coherence_long.mat')
load('../output/short/coherence_short.mat')
load('../output/long/coherence_long_filter.mat')
load('../output/short/coherence_short_filter.mat')

%载入已保存的融合需要的参数
load('../output/long/phase_real_long.mat')
load('../output/long/phase_long_delta.mat')
load('../output/short/phase_short_delta.mat')
load('../output/long/unwrapped_phase_long_abs.mat')
load('../output/short/unwrapped_phase_short_abs.mat')
load('../output/long/flat_phase_long.mat')
load('../output/short/flat_phase_short.mat')
load('../output/long/unwrapped_phase_long_abs_deflat.mat')
load('../output/short/unwrapped_phase_short_abs_deflat.mat')
% load('../output/unwrapped_phase_deflat_long_fusion_judge.mat')
% load('../output/short_PDV_normalized.mat')
% load('../output/long_PDV_normalized.mat')
% load('../output/long_fusion_PDV_normalized.mat')

%长基线下真实的去平地相位
phase_real_long_deflat = phase_real_long + flat_phase_long; 


%% 如果用去平地相位后的绝对解缠相位融合
short_deflat_PDV = PhaseDerivativeVariance_r1(unwrapped_phase_short_abs_deflat);
save([folder,'short_deflat_PDV.mat'],'short_deflat_PDV','-v7.3');
%对相位导数方差进行归一化
short_PDV1 = short_deflat_PDV(:)';
short_PDV2 = mapminmax(short_PDV1,0,1);
short_deflat_PDV_normalized = reshape(short_PDV2,size(short_deflat_PDV));
clear short_PDV1 short_PDV2
save([folder,'short_deflat_PDV_normalized.mat'],'short_deflat_PDV_normalized','-v7.3');

long_deflat_PDV = PhaseDerivativeVariance_r1(unwrapped_phase_long_abs_deflat);
save([folder,'long_deflat_PDV.mat'],'long_deflat_PDV','-v7.3');
%对相位导数方差进行归一化
long_PDV1 = long_deflat_PDV(:)';
long_PDV2 = mapminmax(long_PDV1,0,1);
long_deflat_PDV_normalized = reshape(long_PDV2,size(long_deflat_PDV));
clear long_PDV1 long_PDV2
save([folder,'long_deflat_PDV_normalized.mat'],'long_deflat_PDV_normalized','-v7.3');

long_real_deflat_PDV = PhaseDerivativeVariance_r1(phase_real_long_deflat);
save([folder,'long_real_deflat_PDV.mat'],'long_real_deflat_PDV','-v7.3');
%对相位导数方差进行归一化
long_PDV1 = long_deflat_PDV(:)';
long_PDV2 = mapminmax(long_PDV1,0,1);
long_real_PDV_normalized = reshape(long_PDV2,size(long_real_deflat_PDV));
clear long_PDV1 long_PDV2
save([folder,'long_real_PDV_normalized.mat'],'long_real_PDV_normalized','-v7.3');


%% 读入数据
coherence_short_abs = coherence_short;%滤波前短基线相干系数
coherence_long_abs = coherence_long;%滤波前长基线相干系数
%跳变点140
% coherence_short_abs = coherence_short_filter;%滤波后短基线相干系数
% coherence_long_abs = coherence_long_filter;%滤波后长基线相干系数
%跳变点为114,但是误差均值和标准差比用滤波前相干系数要大一些

%%计算相位比（此为双基线相位比，双频数据可以直接使用中心频率之比）
% B_long = 1.5;% B_V1V3
% B_short = 0.5;%B_V1V2
B_radio = B_short/B_long;
% radio0 = B_long/B_short;%参考基线为长基线
%1/radio0 = B_short/B_long; B_long是参考基线，所以ratio(:,:,1)为短基线与长基线的基线比例


%%基线比例
idx_ref       = 2;%1：短基线为参考基线；2：长基线为参考基线
B(:,:,1)  = B_short;
B(:,:,2)  = B_long;
ratio(:,:,1)  = B_short./B(:,:,idx_ref);%ratio(:,:,1)为短基线与长基线（参考基线）的基线比例
ratio(:,:,2)  = B_long./B(:,:,idx_ref);%ratio(:,:,2)为长基线与长基线（参考基线）的基线比例

%%传递数据
% 相关系数
[row,col]     = size(coherence_long_abs);
gamma         = zeros(row,col,2);
gamma(:,:,1)  = abs(coherence_short_abs);
gamma(:,:,2)  = coherence_long_abs;
% 解缠相位
Psi_MB        = ones(row,col,2);

%%融合判断-根据残差点和相干系数判断是选择替换还是最大似然融合
%如果用去平地相位后的绝对解缠相位融合
Psi_MB(:,:,1) = unwrapped_phase_short_abs_deflat;
Psi_MB(:,:,2) = unwrapped_phase_long_abs_deflat;

%%step2 -- 融合处理,粗融合，搜索步长长，搜索区间大
% short baseline
Psi_sub1   = Psi_MB(:,:,1);
ratio_sub1 = ratio(:,:,1);
gamma_sub1 = gamma(:,:,1);
% long baseline
Psi_sub2   = Psi_MB(:,:,2);
ratio_sub2 = ratio(:,:,2);
gamma_sub2 = gamma(:,:,2);


%经过验证，还是用unwrap_estimate标记跳变点更合理，pdv_jump_mark_phase只是用来辅助融合的函数，标记不连续点
%短基线标记跳变点
J_short_deflat1 = unwrap_estimate(unwrapped_phase_short_abs_deflat(2:end-1,2:end-1));%短基线跳变点个数
[J_short_deflat2,dz_short1] = pdv_jump_mark_phase(short_deflat_PDV_normalized,unwrapped_phase_short_abs_deflat(3:end-1,3:end-1));

%长基线标记跳变点
J_long_deflat1 = unwrap_estimate(unwrapped_phase_long_abs_deflat(2:end-1,2:end-1));%长基线跳变点个数
[J_long_deflat2,dz_long1] = pdv_jump_mark_phase(long_deflat_PDV_normalized,unwrapped_phase_long_abs_deflat(3:end-1,3:end-1));

%真实的长基线去平地相位
J_long_real_deflat1 = unwrap_estimate(phase_real_long_deflat(2:end-1,2:end-1));%长基线跳变点个数
[J_long_real_deflat2,dz_long_real1] = pdv_jump_mark_phase(long_real_PDV_normalized,phase_real_long_deflat(3:end-1,3:end-1));



%% fusion
tic;
display('进行融合处理......');

%相位导数方差替换与最大似然融合
[unwrapped_phase_sub2] = MLUPE_PDV(Psi_sub1,Psi_sub2,short_deflat_PDV_normalized,long_deflat_PDV_normalized,...
    ratio_sub1,ratio_sub2,gamma_sub1,gamma_sub2,B_short,B_long,idx_ref);
temp2 = unwrapped_phase_sub2; 
unwrapped_phase_deflat_long_fusion_judge = temp2(1:row,1:col);
J_long_deflat_fusion1 = unwrap_estimate(unwrapped_phase_deflat_long_fusion_judge(2:end-1,2:end-1));%融合相位跳变点个数

long_deflat_fusion_PDV = PhaseDerivativeVariance_r1(unwrapped_phase_deflat_long_fusion_judge);
%对相位导数方差进行归一化
long_deflat_fusion_PDV1 = long_deflat_fusion_PDV(:)';
long_deflat_fusion_PDV2 = mapminmax(long_deflat_fusion_PDV1,0,1);
long_fusion_deflat_PDV_normalized = reshape(long_deflat_fusion_PDV2,size(long_deflat_PDV));
clear long_deflat_fusion_PDV1 long_deflat_fusion_PDV2
[J_long_deflat_fusion2,dz_fusion1] = pdv_jump_mark_phase(long_fusion_deflat_PDV_normalized,unwrapped_phase_deflat_long_fusion_judge(3:end-1,3:end-1));

save([folder,'long_deflat_fusion_PDV.mat'],'long_deflat_fusion_PDV','-v7.3');
save([folder,'long_fusion_deflat_PDV_normalized.mat'],'long_fusion_deflat_PDV_normalized','-v7.3');
save([folder,'unwrapped_phase_deflat_long_fusion_judge.mat'],'unwrapped_phase_deflat_long_fusion_judge','-v7.3');

display('完成！');
toc


%% 融合结果精度评估
unwrapped_phase_long_fusion1 = unwrapped_phase_deflat_long_fusion_judge - flat_phase_long;%加回平地相位
phase_long_fusion_delta1 =  phase_real_long - unwrapped_phase_long_fusion1;%相位误差
figure;imagesc(phase_long_fusion_delta1);
save([folder,'phase_long_fusion_delta1.mat'],'phase_long_fusion_delta1','-v7.3');
% phase_long_fusion_mean_error1 = nanmean(abs(phase_long_fusion_delta1(:)));%绝对误差均值
phase_long_fusion_mean_error1 =nanmean(phase_long_fusion_delta1(:));%绝对误差均值，不取abs比较好，更能显示融合的提升效果
phase_long_fusion_relative_error1 = phase_long_fusion_delta1 ./ unwrapped_phase_long_fusion1;%相对误差
%相对误差标准差
phase_long_fusion_relative_error_std1 = nanstd(reshape(phase_long_fusion_delta1,size(unwrapped_phase_long_abs,1)*size(unwrapped_phase_long_abs,2),1));

%% 精度评估
%短基线
figure;imagesc(phase_short_delta);
% phase_short_mean_error = nanmean(abs(phase_short_delta(:)));%绝对误差均值
phase_short_mean_error = nanmean(phase_short_delta(:));%绝对误差均值
phase_short_relative_error = phase_short_delta ./ unwrapped_phase_short_abs;%相对误差
%相对误差标准差
phase_short_relative_error_std = nanstd(reshape(phase_short_delta,size(unwrapped_phase_short_abs,1)*size(unwrapped_phase_short_abs,2),1));

%长基线
figure;imagesc(phase_long_delta);
% phase_long_mean_error = nanmean(abs(phase_long_delta(:)));%绝对误差均值
phase_long_mean_error = nanmean(phase_long_delta(:));%绝对误差均值
phase_long_relative_error = phase_long_delta ./ unwrapped_phase_long_abs;%相对误差
%相对误差标准差
phase_long_relative_error_std = nanstd(reshape(phase_long_delta,size(unwrapped_phase_long_abs,1)*size(unwrapped_phase_long_abs,2),1));



%% 如果用绝对解缠相位融合
%经过验证发现，用加回平地相位的绝对解缠相位直接融合，回加大误差，融合效果不好
short_PDV = PhaseDerivativeVariance_r1(unwrapped_phase_short_abs);
save([folder,'short_PDV.mat'],'short_PDV','-v7.3');
%对相位导数方差进行归一化
short_PDV1 = short_PDV(:)';
short_PDV2 = mapminmax(short_PDV1,0,1);
short_PDV_normalized = reshape(short_PDV2,size(short_PDV));
clear short_PDV1 short_PDV2
save([folder,'short_PDV_normalized.mat'],'short_PDV_normalized','-v7.3');

long_PDV = PhaseDerivativeVariance_r1(unwrapped_phase_long_abs);
save([folder,'long_PDV.mat'],'long_PDV','-v7.3');
%对相位导数方差进行归一化
long_PDV1 = long_PDV(:)';
long_PDV2 = mapminmax(long_PDV1,0,1);
long_PDV_normalized = reshape(long_PDV2,size(long_PDV));
clear long_PDV1 long_PDV2
save([folder,'long_PDV_normalized.mat'],'long_PDV_normalized','-v7.3');

%% 读入数据
coherence_short_abs = coherence_short;%滤波前短基线相干系数
coherence_long_abs = coherence_long;%滤波前长基线相干系数


%%基线比例
idx_ref       = 2;%1：短基线为参考基线；2：长基线为参考基线
B(:,:,1)  = B_short;
B(:,:,2)  = B_long;
ratio(:,:,1)  = B_short./B(:,:,idx_ref);%ratio(:,:,1)为短基线与长基线（参考基线）的基线比例
ratio(:,:,2)  = B_long./B(:,:,idx_ref);%ratio(:,:,2)为长基线与长基线（参考基线）的基线比例

%%传递数据
% 相关系数
[row,col]     = size(coherence_long_abs);
gamma         = zeros(row,col,2);
gamma(:,:,1)  = abs(coherence_short_abs);
gamma(:,:,2)  = coherence_long_abs;
% 解缠相位
Psi_MB        = ones(row,col,2);

%如果用绝对解缠相位融合
Psi_MB(:,:,1) = unwrapped_phase_short_abs;
Psi_MB(:,:,2) = unwrapped_phase_long_abs;

%%step2 -- 融合处理,粗融合，搜索步长长，搜索区间大
% short baseline
Psi_sub1   = Psi_MB(:,:,1);
ratio_sub1 = ratio(:,:,1);
gamma_sub1 = gamma(:,:,1);
% long baseline
Psi_sub2   = Psi_MB(:,:,2);
ratio_sub2 = ratio(:,:,2);
gamma_sub2 = gamma(:,:,2);

%经过验证，还是用unwrap_estimate标记跳变点更合理，pdv_jump_mark_phase只是用来辅助融合的函数，标记不连续点
%标记跳变点
J_short1 = unwrap_estimate(unwrapped_phase_short_abs(2:end-1,2:end-1));%短基线跳变点个数
[J_short2,dz_short2] = pdv_jump_mark_phase(short_PDV_normalized,unwrapped_phase_short_abs(3:end-1,3:end-1));

%标记跳变点
J_long1 = unwrap_estimate(unwrapped_phase_long_abs(2:end-1,2:end-1));%长基线跳变点个数
[J_long2,dz_long2] = pdv_jump_mark_phase(long_PDV_normalized,unwrapped_phase_long_abs(3:end-1,3:end-1));


%% fusion
tic;
display('进行融合处理......');

%相位导数方差替换与最大似然融合
[unwrapped_phase_sub2] = MLUPE_PDV(Psi_sub1,Psi_sub2,short_PDV_normalized,long_PDV_normalized,...
    ratio_sub1,ratio_sub2,gamma_sub1,gamma_sub2,B_short,B_long,idx_ref);
temp2 = unwrapped_phase_sub2; 
unwrapped_phase_long_fusion_judge = temp2(1:row,1:col);
J_long_fusion1 = unwrap_estimate(unwrapped_phase_long_fusion_judge(2:end-1,2:end-1));%融合相位跳变点个数

long_fusion_PDV = PhaseDerivativeVariance_r1(unwrapped_phase_long_fusion_judge);
%对相位导数方差进行归一化
long_fusion_PDV1 = long_fusion_PDV(:)';
long_fusion_PDV2 = mapminmax(long_fusion_PDV1,0,1);
long_fusion_PDV_normalized = reshape(long_fusion_PDV2,size(long_PDV));
clear long_fusion_PDV1 long_fusion_PDV2
[J_longt_fusion2,dz_fusion2] = pdv_jump_mark_phase(long_fusion_PDV_normalized,unwrapped_phase_long_fusion_judge(3:end-1,3:end-1));

save([folder,'long_fusion_PDV.mat'],'long_fusion_PDV','-v7.3');
save([folder,'long_fusion_PDV_normalized.mat'],'long_fusion_PDV_normalized','-v7.3');
save([folder,'unwrapped_phase_long_fusion_judge.mat'],'unwrapped_phase_long_fusion_judge','-v7.3');

display('完成！');
toc

%% 精度评估
% %短基线
% figure;imagesc(phase_short_delta);
% phase_short_mean_error = nanmean(abs(phase_short_delta(:)));%绝对误差均值
% % phase_short_mean_error = nanmean(phase_short_delta(:));%绝对误差均值
% phase_short_relative_error = phase_short_delta ./ unwrapped_phase_short_abs;%相对误差
% %相对误差标准差
% phase_short_relative_error_std = nanstd(reshape(phase_short_delta,size(unwrapped_phase_short_abs,1)*size(unwrapped_phase_short_abs,2),1));
% 
% %长基线
% figure;imagesc(phase_long_delta);
% phase_long_mean_error = nanmean(abs(phase_long_delta(:)));%绝对误差均值
% % phase_long_mean_error = nanmean(phase_long_delta(:));%绝对误差均值
% phase_long_relative_error = phase_long_delta ./ unwrapped_phase_long_abs;%相对误差
% %相对误差标准差
% phase_long_relative_error_std = nanstd(reshape(phase_long_delta,size(unwrapped_phase_long_abs,1)*size(unwrapped_phase_long_abs,2),1));
% 
% %融合结果
% unwrapped_phase_long_fusion2 = unwrapped_phase_long_fusion_judge;
% phase_long_fusion_delta2 =  phase_real_long - unwrapped_phase_long_fusion2;%相位误差
% figure;imagesc(phase_long_fusion_delta2);
% save([folder,'phase_long_fusion_delta2.mat'],'phase_long_fusion_delta2','-v7.3');
% phase_long_fusion_mean_error2 = nanmean(abs(phase_long_fusion_delta2(:)));%绝对误差均值
% % phase_long_fusion_mean_error2 = nanmean(phase_long_fusion_delta(:));%绝对误差均值
% phase_long_fusion_relative_error2 = phase_long_fusion_delta2 ./ unwrapped_phase_long_fusion2;%相对误差
% %相对误差标准差
% phase_long_fusion_relative_error_std2 = nanstd(reshape(phase_long_fusion_delta2,size(unwrapped_phase_long_abs,1)*size(unwrapped_phase_long_abs,2),1));


