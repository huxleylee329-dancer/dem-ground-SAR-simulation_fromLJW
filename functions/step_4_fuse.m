% --------------------------双频融合-------------------------------
%% 读入数据
warning off;
display('读入数据......');
in_folder_unwrappedphase = '.\Kernel\Output3_unwrap\';
in_folder_coherence = '.\Kernel\Output1_registration\';
in_folder_para = '.\Kernel\parameter\';
out_folder = '.\Kernel\Output3_unwrap\';


file_name = [in_folder_unwrappedphase 'unwrapped_phase_low.mat' ]; 
load(file_name);

file_name = [in_folder_coherence 'coherence_low.mat' ]; 
load(file_name);

[r,c] = size(unwrapped_phase_low);
coherence_low_1 = ones(r,c);
coherence_low_1(2:end-1,2:end-1) = coherence_low;
coherence_low = coherence_low_1;


file_name = [in_folder_unwrappedphase,'unwrapped_phase_high', '.mat' ]; 
load(file_name);


file_name = [in_folder_coherence,'coherence_high', '.mat' ]; 
load(file_name);
[r,c] = size(unwrapped_phase_high);
coherence_high_1 = ones(r,c);
coherence_high_1(2:end-1,2:end-1) = coherence_high;
coherence_high = coherence_high_1;

%%计算相位比（此为双基线相位比，双频数据可以直接使用中心频率之比）
[f_low,~,~,~] = textread([in_folder_para 'low\auxi_mast.txt'], '%f %f %n %n');
[f_high,~,~,~] = textread([in_folder_para 'high\auxi_mast.txt'], '%f %f %n %n');
display('完成！');
radio0 = f_high/f_low;
%%将短基线相位调整至参考基线相位
G1_2 =unwrapped_phase_low * radio0;

% 线校正长基线相位的模糊数矩阵
K_const_mat= round((G1_2 - unwrapped_phase_high)./2./pi);

unwrapped_phase_high = unwrapped_phase_high+ 2 * pi * K_const_mat; 

%% 传递数据
% 相关系数
[row,col]     = size(coherence_high);
gamma         = zeros(row,col,2);
gamma(:,:,1)  = abs(coherence_low);
gamma(:,:,2)  = coherence_high;
% 解缠相位
Psi_MB        = ones(row,col,2);
Psi_MB(:,:,1) = unwrapped_phase_low;
Psi_MB(:,:,2) = unwrapped_phase_high;

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
phase=unwrapped_phase_high; 
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
display('完成！');
%%保存数据
display('保存数据......');
fig2 = figure;
set(fig2, 'visible', 'off');
bar([J_long;J_ronghe],0.15);
set(gca,'xTicklabel',{'高频跳变点数','融合跳变点数'});
saveas(fig2, [out_folder 'compare.jpg']);



fig1 = figure;
set(fig1,'visible','off');
imagesc(ronghe);colorbar;
saveas(fig1,[out_folder 'unwrapped_phase_fused.jpg']);
% save([folder,'ronghe.mat'],'ronghe');
% save([folder,'unwrapped_phase_long.mat'],'unwrapped_phase_long');
display('完成！');


