% ----------------
% ---------------- step 2
% ---------------- 长短基线相干配准 

clear all
close all
clc
%% 输入处理ID
process_ID_master  = '500_mast';
process_ID_slave   = '500_slave';
process_ID = [process_ID_master,'&',process_ID_slave,'_insar'];
%% 载入SSC
load(['C:\Users\Rick\Desktop\things\曾国兵\InSAR处理\output\',process_ID_master,'\SSC.mat'],'SSC');
ssc_master = SSC;
load(['C:\Users\Rick\Desktop\things\曾国兵\InSAR处理\output\',process_ID_slave,'\SSC.mat'],'SSC');
ssc_slave  = SSC;
%% 裁剪至相同大小
% 取主辅图像公共区域
[nr1,nc1]  = size(ssc_master);
[nr2,nc2]  = size(ssc_slave);
nr         = min(nr1,nr2);
nc         = min(nc1,nc2);
SSC_master = ssc_master(1:nr - 8,1:nc - 8);
SSC_slave  = ssc_slave(9:nr,9:nc);
% 更新主图像裁剪范围
[~] = cut_image_index_update(nr1,nc1,nr,nc,process_ID_master);
[~] = cut_image_index_update(nr2,nc2,nr,nc,process_ID_slave);
%% 干涉处理
% 粗配准

[nr1,nc1]            = size(SSC_master);
[SSC_master_regis1,SSC_slave_regis1,move_r,move_c] = registration_pixel(SSC_master,SSC_slave,nr1,nc1);
[nr2,nc2]            = size(SSC_master_regis1);
% 更新主图像裁剪范围
[~] = after_regis_image_index_update(move_r,move_c,process_ID_master,process_ID_slave);

% 亚像素级配准
[~,SSC_slave_regis2] = regis_subpixel(SSC_master_regis1,SSC_slave_regis1,nr2,nc2);

coherence = Calculation_Coherence_Coefficient_complex_no_phase(SSC_master_regis1,SSC_slave_regis2);
folder=['C:\Users\Rick\Desktop\things\曾国兵\InSAR处理\Mid\',process_ID,'\'];
save([folder,'coherence.mat'],'coherence');
save([folder,'SSC_master_regis1.mat'],'SSC_master_regis1');
save([folder,'SSC_slave_regis2.mat'],'SSC_slave_regis2');
% 缠绕干涉相位
wrapped_phase        = interferometric_phase(SSC_master_regis1,SSC_slave_regis2);
save([folder,'wrapped_phase.mat'],'wrapped_phase');
% 干涉相位滤波
   win_size  = 25;
     win_pre = 23;
     tic;
     [phase_removeFlat,Flat_phase_wrapped] = removeFlat_range_20161126(wrapped_phase); %移除平地相位
     phase_removeFlat_filter = Improved_Slope_Adaptive_filter_parallel(phase_removeFlat,win_size,win_pre); %坡度自适应滤波
     save ([folder,'phase_removeFlat_filter.mat'], 'phase_removeFlat_filter');    %保存残余相位
     save ([folder,'Flat_phase_wrapped.mat'], 'Flat_phase_wrapped');  %保存平地相位
     wrapped_phase_filter_gray        = mat2gray(phase_removeFlat_filter); 
     imwrite(wrapped_phase_filter_gray,[folder,'phase_removeFlat_filter.tif']);
     toc

 %%%相位解缠（最小费用流法）

[unwrapped_phase_removeFlat]=MCF(phase_removeFlat_filter);%残余相位解缠
[flat_phase]= FlatUnwrap2D(Flat_phase_wrapped);%平地相位解缠
unwrapped_phase=unwrapped_phase_removeFlat-flat_phase;%平地相位补偿
save ([folder,'unwrapped_phase.mat'], 'unwrapped_phase');%保存解缠相位
% save ([folder,'unwrapped_phase_removeFlat.mat'], 'unwrapped_phase_removeFlat');
% save ([folder,'flat_phase.mat'], 'flat_phase');
unwrapped_phase_gray        = mat2gray(unwrapped_phase); 
imwrite(unwrapped_phase_gray,[folder,'unwrapped_phase.tif']);