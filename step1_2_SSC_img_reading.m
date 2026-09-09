% ----------------
% ---------------- step 1
% ---------------- 科大仿真数据读取与存储 --- 尤亚楠 
% ---------------- 20160329
% ----------------
clear all
close all
clc
%% 输入处理ID
process_ID  = '500_slave';

%% 读取SSC及辅助数据
M1=csvread('..\Output\500_mast&500_slave_insar\slave1_SLC_real.csv');
M2=csvread('..\Output\500_mast&500_slave_insar\slave1_SLC_imag.csv');
row1=1; 
col1=1;
[rowend,colend]=size(M1);
Na          = rowend;
Nr          = colend;
ssc = complex(M1,M2);
SSC         = ssc(row1:1:rowend,col1:1:colend);

%% 量化
win_size    = 9;   %window size
sigma       = 1.6; %sigma for gaussian filter
displayflag = 0;   %0:on figure;1:figure
[G_filter]  = sar_image_quantify(SSC,win_size,sigma,displayflag);
%% 存储
folder = ['..\Output\',process_ID,'\'];
if ~isdir(folder)
    mkdir(folder)
end
% 1-SSC mat tif
imwrite(G_filter,['..\Output\',process_ID,'\SSC_image.tif'])
save(['..\Output\',process_ID,'\SSC.mat'],'SSC');
% 2-SSC infor
fid1 = fopen(['..\Output\',process_ID,'\SSC_infor.txt'],'wt');
fprintf(fid1,'%s %i \n','图像原始行: ',Na);
fprintf(fid1,'%s %i \n','图像原始列: ',Nr);
fprintf(fid1,'%s %i \n','裁剪起始行: ',row1);
fprintf(fid1,'%s %i \n','裁剪终止行: ',rowend);
fprintf(fid1,'%s %i \n','裁剪起始列: ',col1);
fprintf(fid1,'%s %i \n','裁剪终止列: ',colend);
fclose(fid1);