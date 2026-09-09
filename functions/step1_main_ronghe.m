clear all
close all
clc
tic
 
%%此为长短基线融合程序，双频融合与双基线类似，低频相当于短基线，高频相当于长基线，将基线比修改为频率比即可。
%%双频融合的输入可为解缠相位\去平地后相位\绝对相位。因为相位比均频率比，因而都能用这种方法进行似然估计。
%% 读入数据
folder = ['..\Input\short\'];
file_name = [folder,'unwrapped_phase', '.mat' ]; 
load(file_name)
unwrapped_phase_short=unwrapped_phase;
file_name = [folder,'coherence', '.mat' ]; 
load(file_name);
coherence_short=coherence;
[r,c] = size(unwrapped_phase_short);
coherence_short_1 = ones(r,c);
coherence_short_1(2:end-1,2:end-1) = coherence_short;
coherence_short = coherence_short_1;

folder = ['..\Input\long\'];
file_name = [folder,'unwrapped_phase', '.mat' ]; 
load(file_name);
unwrapped_phase_long = unwrapped_phase;


file_name = [folder,'coherence', '.mat' ]; 
load(file_name);
coherence_long=coherence;
[r,c] = size(unwrapped_phase_long);
coherence_long_1 = ones(r,c);
coherence_long_1(2:end-1,2:end-1) = coherence_long;
coherence_long = coherence_long_1;

%%计算相位比（此为双基线相位比，双频数据可以直接使用中心频率之比）
radio0=mean((unwrapped_phase_long(:,end)-unwrapped_phase_long(:,1))./(unwrapped_phase_short(:,end)-unwrapped_phase_short(:,1)));

%%将短基线相位调整至参考基线相位
G1_2 =unwrapped_phase_short.*radio0;

% 线校正长基线相位的模糊数矩阵
K_const_mat= round((G1_2 - unwrapped_phase_long)./2./pi);

unwrapped_phase_long = unwrapped_phase_long+ 2 * pi * K_const_mat; 

%% 传递数据
% 相关系数
[row,col]     = size(coherence_long);
gamma         = zeros(row,col,2);
gamma(:,:,1)  = abs(coherence_short);
gamma(:,:,2)  = coherence_long;
% 解缠相位
Psi_MB        = ones(row,col,2);
Psi_MB(:,:,1) = unwrapped_phase_short;
Psi_MB(:,:,2) = unwrapped_phase_long;

%% 基线比例
idx_ref       = 2;
ratio= ones(row,col,2);%1：短基线为参考基线；2：长基线为参考基线
ratio(:,:,1)  = ratio(:,:,1)*1/radio0;

%% step2 -- 融合处理
%粗融合，搜索步长长，搜索区间大
% short baseline
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

%%精融合，搜索步长小，搜索区间小
Psi_sub3=unwrapped_phase_fusion;

 [unwrapped_phase_sub] = MLUPE2(Psi_sub1,Psi_sub2,ratio_sub1,...
              ratio_sub2,gamma_sub1,gamma_sub2,Psi_sub3);
temp3 =  unwrapped_phase_sub; 
ronghe = temp3(1:row,1:col);
folder = ['..\Output\'];
if ~isdir(folder)
    mkdir(folder)
end


%%%%%精度评估，跳变点个数%%%%%%
phase=unwrapped_phase_long; 
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

%%保存数据
save([folder,'ronghe.mat'],'ronghe');
save([folder,'unwrapped_phase_long.mat'],'unwrapped_phase_long');



