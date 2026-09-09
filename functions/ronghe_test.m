%% 去平地
para_path = '.\差频数据\9.2GHz\';% 控制点参数路径
out_path = '.\差频数据\';% 去平地结果保存路径 
SLC_path = '.\差频数据\9.2GHz\';% 主辅图像路径
warning off;
% 检查并读入参数
frequency_mast = 9.2e9;
% frequency_slave = frequency_mast - 2.5*18e6;
frequency_slave = frequency_mast;
nr1 = 260;
nc1 = 260;


% 载入主图像（单视复图像）
if ~exist([SLC_path 'mast.mat'], 'file')
    display('主星单视复图像文件不存在！请检查......');
    pause(5);
    quit;
else
    load([SLC_path 'mast.mat']);
end
load([para_path '去参考地形\参考地形\' 'Slantrange_s.mat']);
load([para_path '去参考地形\参考地形\' 'Slantrange_m.mat']);
lambda1 = 299792458/frequency_mast;
flat_range_master = Slantrange_m;
flat_phase_master = flat_range_master * 4 * pi / lambda1;  %主图参考地形相位
% 去参考地形相位
deflat_mast_low = mast.*exp(1i*flat_phase_master);
% deflat_mast_low = mast;
% 保存结果
display('保存结果......');
flat_phase_mast_low = flat_phase_master;
save([out_path 'flat_phase_mast_low.mat'], 'flat_phase_mast_low');% 保存绝对平地相位（后续加回来）
save([out_path 'deflat_mast_low.mat'], 'deflat_mast_low');% 保存去平地之后的主图像
fig1 = figure;
set(fig1, 'visible', 'off');
imagesc(angle(exp(1i*flat_phase_master)));colorbar;
saveas(fig1, [out_path 'flat_phase_low_master.jpg']);% 保存平地相位图片（缠绕后的）
display('完成！');




%% 辅图像去平地
display('辅图去平地......');

% 载入辅图像（单视复图像）
if ~exist([SLC_path 'slave.mat'], 'file')
    display('辅星单视复图像文件不存在！请检查......');
    pause(5);
    quit;
else
    load([SLC_path 'slave.mat']);
end
lambda2 = 299792458/frequency_slave;
flat_range_slave = Slantrange_s;
flat_phase_slave_low = flat_range_slave * 4 * pi / lambda2;  %主图平地相位
%去平地效应
deflat_slave_low = slave.*exp(1i*flat_phase_slave_low);
% deflat_slave_low = slave;
display('保存结果......');
save([out_path 'deflat_slave_low.mat'], 'deflat_slave_low');% 保存去平地之后的辅图像
save([out_path 'flat_phase_slave_low.mat'], 'flat_phase_slave_low');% 保存绝对平地相位

fig2 = figure;
set(fig2, 'visible', 'off');
imagesc(angle(exp(1i*flat_phase_slave_low)));colorbar;
saveas(fig2, [out_path 'flat_phase_low_slave.jpg']);% 保存平地相位图片（缠绕后的）
display('完成！');

display('完成！');
display('去平地完成！');









%% 配准 

outfolder = out_path;%输出路径
infolder = out_path;
%% 载入低频单视复图像
display('载入单视复图像......');
load([infolder 'deflat_mast_low.mat']);
mast = deflat_mast_low;
load([infolder 'deflat_slave_low.mat']);
slave = deflat_slave_low;
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
% wrapped_phase = interferometric_phase(slave_regis2,mast_regis1);
display('完成！');
display('计算相干图......');
coherence_low = Calculation_Coherence_Coefficient_complex_no_phase(mast_regis1,slave_regis2);
display('完成！');

fig1 = figure;
set(fig1, 'visible', 'off');
imagesc(wrapped_phase);
colorbar;

saveas(fig1, [outfolder 'wrapped_phase_low.jpg']);
save([outfolder,'wrapped_phase_low.mat'],'wrapped_phase');

fig2 = figure;
set(fig2, 'visible', 'off');
imshow(coherence_low);
colorbar;
saveas(fig2, [outfolder 'coherence_low.jpg']);
save([outfolder 'coherence_low'], 'coherence_low');
















%% 滤波
win_size = 11;
%% 载入含噪相位
display('载入含噪相位......');
out_folder = out_path;
phase_path = [out_path 'wrapped_phase_low.mat'];
if exist(phase_path,'file')
    load(phase_path);
    display('载入完成！');
    display('进行滤波（计算时间较长）......');
    wrapped_phase_filtered = Improved_Slope_Adaptive_filter_parallel(wrapped_phase,win_size,win_size); %坡度自适应滤波
    display('滤波完成！');
    display('保存结果......');
    fig1 = figure;
    set(fig1, 'visible', 'off');
    imagesc(wrapped_phase_filtered);
    colorbar;
    save([out_folder 'wrapped_phase_filtered_low.mat'], 'wrapped_phase_filtered');
    saveas(fig1, [out_folder 'wrapped_phase_filtered_low.jpg']);
    display('完成！');
else
    display('不存在含噪相位文件！请检查...');
end
















%% 解缠
warning off;
in_folder = out_folder;
%% 读取低频缠绕相位数据
if exist([in_folder 'wrapped_phase_filtered_low.mat'], 'file')
    display('载入缠绕相位数据......');
    load([in_folder 'wrapped_phase_filtered_low.mat']);
    display('完成！');
    display('进行相位解缠......');
    unwrapped_phase_low = MCF(wrapped_phase_filtered);
    display('完成！');
    
    display('保存结果......');
    fig1 = figure;
    set(fig1, 'visible', 'off');
    imagesc(unwrapped_phase_low);
    colorbar;
    saveas(fig1, [out_folder 'unwrapped_phase_low.jpg']);
    save([out_folder 'unwrapped_phase_low.mat'], 'unwrapped_phase_low');
    display('完成！');
else
    display('低频缠绕相位数据不存在！请检查......');
end


%% 

load([para_path 'conp1.mat']);
load([para_path 'conp2.mat']);
load([para_path 'Slantrange_s.mat']);
load([para_path 'Slantrange_m.mat']);
%% 加回参考地形相位
unwrapped_phase_low_add = unwrapped_phase_low + (flat_phase_slave_low - flat_phase_mast_low);
% unwrapped_phase_low_add = - unwrapped_phase_low_add;
%% 校正至绝对相位
[k,~] = size(conp1);
del = ones(k, 2)*3965;

conp1(:,1:2) = conp1(:,1:2) - del;
conp2(:,1:2) = conp2(:,1:2) - del;
count = 0;
for iii = 1:k
    if (conp1(iii, 1) > 0 && conp1(iii, 1) < 261) && (conp1(iii, 2) > 0 && conp1(iii, 2) < 261)
        count = count + 1;
    end
end
index = zeros(count, 1);
jj = 0;
for iii = 1:k
    if (conp1(iii, 1) > 0 && conp1(iii, 1) < 261) && (conp1(iii, 2) > 0 && conp1(iii, 2) < 261)
        jj = jj + 1;
        index(jj) = iii;
    end
end
conp1_bk = conp1(index,:);
conp2_bk = conp2(index,:);
phi_1 = conp1_bk(:,3) / lambda1 * 4 * pi;
phi_2 = conp2_bk(:,3) / lambda2 * 4 * pi;
phi_abs = abs(phi_1 - phi_2);
phase_abs = zeros(1,count);
for iii = 1:count
    phase_abs(iii) = phi_abs(iii) - unwrapped_phase_low_add(conp1_bk(iii, 1), conp1_bk(iii, 2));
end
k_point = round(phase_abs/(2*pi));
k_value = unique(k_point);
[~, k_value_number] = size(k_value);
n = 1;
for ii = 1:k_value_number
    m = sum(k_point==k_value(1,ii));
    if m > n
        n = m;
        k_real = k_value(ii);
    end
end

unwrapped_phase_low_abs = unwrapped_phase_low_add + 2 * pi * k_real;
%% 再次去掉参考地形相位
unwrapped_phase_low = unwrapped_phase_low_abs - (flat_phase_slave_low - flat_phase_mast_low);


%%  真实相位差
phase_real =  (Slantrange_s/lambda2*4*pi - Slantrange_m/lambda1*4*pi);
delta = phase_real - unwrapped_phase_low_abs;
figure;imagesc(delta);colorbar;



%% 双频融合
unwrapped_phase_high_t = unwrapped_phase_low;
coherence_high_t = coherence_low;
% unwrapped_phase_low_t = unwrapped_phase_low;
% coherence_low_t = coherence_low;
%% 读入数据




[r,c] = size(unwrapped_phase_low_t);
coherence_low_1 = ones(r,c);
coherence_low_1(2:end-1,2:end-1) = coherence_low_t;
coherence_low_t = coherence_low_1;


[r,c] = size(unwrapped_phase_high_t);
coherence_high_1 = ones(r,c);
coherence_high_1(2:end-1,2:end-1) = coherence_high_t;
coherence_high_t = coherence_high_1;

%%计算相位比（此为双基线相位比，双频数据可以直接使用中心频率之比）
f_high = 9.6e9;
f_low = 9.2e9;
radio0 = f_high/f_low;
%%将短基线相位调整至参考基线相位
G1_2 =unwrapped_phase_low_t * radio0;
G1_2=
%% 线校正长基线相位的模糊数矩阵
K_const_mat= round((G1_2 - unwrapped_phase_high_t)./2./pi);

unwrapped_phase_high_t = unwrapped_phase_high_t+ 2 * pi * K_const_mat; 
%% 
unwrapped_phase_high_t=unwrapped_phase_high_t-(mean(mean(unwrapped_phase_high_t)) - mean(mean(G1_2)));

%% 传递数据
% 相关系数
[row,col]     = size(coherence_high_t);
gamma         = zeros(row,col,2);
gamma(:,:,1)  = abs(coherence_low_t);
gamma(:,:,2)  = coherence_high_t;
% 解缠相位
Psi_MB        = ones(row,col,2);
Psi_MB(:,:,1) = unwrapped_phase_low_t;
Psi_MB(:,:,2) = unwrapped_phase_high_t;

%% 基线比例
idx_ref       = 2;
ratio= ones(row,col,2);%1：短基线为参考基线；2：长基线为参考基线
ratio(:,:,1)  = ratio(:,:,1)*1/radio0;

%% step2 -- 融合处理
%粗融合，搜索步长长，搜索区间大
% short baseline
display('进行粗融合处理......');
Psi_sub1   = Psi_MB(:,:,1);
ratio_sub1 = ratio(:,:,1);
gamma_sub1 = gamma(:,:,1);
% long baseline
Psi_sub2   = Psi_MB(:,:,2);
ratio_sub2 = ratio(:,:,2);
gamma_sub2 = gamma(:,:,2);
% fusion
[unwrapped_phase_sub] = MLUPE1(Psi_sub1,Psi_sub2,ratio_sub1,...
          ratio_sub2,gamma_sub1,gamma_sub2,Psi_sub2);
temp2 = unwrapped_phase_sub; 
unwrapped_phase_fusion = temp2(1:row,1:col);
display('完成！');
%%精融合，搜索步长小，搜索区间小
display('进行精融合......');
Psi_sub3=unwrapped_phase_fusion;

 [unwrapped_phase_sub] = MLUPE2(Psi_sub1,Psi_sub2,ratio_sub1,...
              ratio_sub2,gamma_sub1,gamma_sub2,Psi_sub3);
temp3 =  unwrapped_phase_sub; 
ronghe = temp3(1:row,1:col);
% folder = ['..\Output\'];
% if ~isdir(folder)
%     mkdir(folder)
% end

display('完成！');
%%%%%精度评估，跳变点个数%%%%%%
display('进行评估......');
phase=unwrapped_phase_low_t;
dx=floor(abs((phase(1:na-1,:)-phase(2:na,:))/2/pi));
dy=floor(abs((phase(:,1:nr-1)-phase(:,2:nr))/2/pi));
dx=dx(2:end-1,2:end-1);
dy=dy(2:end-1,2:end-1);
J_short=sum(sum(dx))+sum(sum(dy));%短基线跳变点个数


phase=unwrapped_phase_high_t; 
[na,nr]=size(phase);
dx=floor(abs((phase(1:na-1,:)-phase(2:na,:))/2/pi));
dy=floor(abs((phase(:,1:nr-1)-phase(:,2:nr))/2/pi));
J_long=sum(sum(dx))+sum(sum(dy)); %长基线跳变点个数

phase=ronghe;
dx=floor(abs((phase(1:na-1,:)-phase(2:na,:))/2/pi));
dy=floor(abs((phase(:,1:nr-1)-phase(:,2:nr))/2/pi));
dx=dx(2:end-1,2:end-1);
dy=dy(2:end-1,2:end-1);
J_ronghe=sum(sum(dx))+sum(sum(dy));%融合相位跳变点个数


%% 加回参考地形相位
load(['.\差频数据\9.2GHz\去参考地形\参考地形\' 'Slantrange_m.mat']);
load(['.\差频数据\9.2GHz\去参考地形\参考地形\' 'Slantrange_s.mat']);
lamda1 = 3e8/9.2e9;
lamda2 = lamda1;
unwrapped_phase_low_t = unwrapped_phase_low_t + (4*pi*Slantrange_s/lamda2 - 4*pi*Slantrange_m/lamda1);


load(['.\差频数据\9.6GHz\去参考地形\参考地形\' 'Slantrange_m.mat']);
load(['.\差频数据\9.6GHz\去参考地形\参考地形\' 'Slantrange_s.mat']);
lamda1 = 3e8/9.6e9;
lamda2 = lamda1;
unwrapped_phase_high_t = unwrapped_phase_high_t + (4*pi*Slantrange_s/lamda2 - 4*pi*Slantrange_m/lamda1);
ronghe = ronghe + (4*pi*Slantrange_s/lamda2 - 4*pi*Slantrange_m/lamda1);

%% 
load([para_path '\去参考地形\参考地形\' 'Slantrange_s.mat']);
load([para_path '\去参考地形\参考地形\' 'Slantrange_m.mat']);
lamda1 = 3e8/9.2e9;
lamda2 = lamda1;
phase_real_1 =  (Slantrange_s/lamda2*4*pi - Slantrange_m/lamda1*4*pi);
delta1 = unwrapped_phase_low_t - phase_real_1;

lamda1 = 3e8/9.6e9;
lamda2 = lamda1;
phase_real_2 =  (Slantrange_s/lamda2*4*pi - Slantrange_m/lamda1*4*pi);
delta2 = unwrapped_phase_high_t - phase_real_2;

delta3 = ronghe - phase_real_2;
display('低频解缠相位精度：');
std(reshape(delta1, 260*260, 1))
display('高频解缠相位精度：');
std(reshape(delta2, 260*260, 1))
display('融合解缠相位精度：');
std(reshape(delta3, 260*260, 1))
