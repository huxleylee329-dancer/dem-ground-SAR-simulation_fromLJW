
%本程序是星载SAR图像仿真程序和干涉处理程序，整理于2021.10.09
%设置系统参数生成主辅图像、斜距、控制点、真实相位等信息
%对主辅图像进行干涉处理：图像配准、去参考地形、相位滤波、相位解缠、控制点校正、相位精度评估、真实地形DEM反演
clear all
close all
clc
first=1;%%%%表示第一次运行程序
%%%%%%卫星轨道参数%%%%%
Drange=2.2711549848484847;%斜距方向采样间隔
DAzimath=2.3784708055877464;
Wlength=3*10^8/(9.6*10^9);
lambda = Wlength;
Prf=3179; %重复频率 hz
C = 4*pi;%单发双收模式:2π,单发单收：4π
%%%%%%%%%%布设地面场景点参数%%%%%%
DEM_center_ll=[0;0]/180*pi; %布设在赤道上
Pointspace=1.5700382546752528e-007/2;% DEM布设间隔，LS间隔为0.5m

%%%%%%%%%%轨道参数%%%%%%%%%%%%%%%
startT=838.357397837519100;
timeSpace=1/2/Prf;
numsate=4401;
%% 读入地面场景DEM点并转换坐标%%%%%%%
% load('..\input\dem_peak.mat');
% height=dem_peak;
load('..\input\DEM_lan.mat');
height=DEM_lan;
[row,clomn]=size(height);
% height=awgn(height,50);
disp('DEM读入完成')
%% %%%转化为N*N*3的矩阵%%%%%%%%%%%%(根据地形大小改变参数)
long=(linspace(-3000,3000,6001)'*ones(1,6001))'*Pointspace+DEM_center_ll(1);
lat=-(linspace(-3000,3000,6001)'*ones(1,6001))*Pointspace+DEM_center_ll(2);
long=long(1:row,1:clomn);
lat=lat(1:row,1:clomn);
DEMll=cat(3,long,lat,height);
DEM=Ell2xyz(DEMll);
disp('DEM转化完成')
%% %%%%%%读入轨道参数并且进行插值%%%%%%%%%%%
path=['..\input\descending_left.csv'];
satepos=compute_center(DEM,Prf,numsate,path);
%% %%%%%%%计算成像卫星的位置%%%%%%%%%%%%%%%%
tic;
display('计算成像位置');
%%%%%%%%%%第一次执行时进行迭代处理，后续次数直接读取%%%%%%%%%%
if first==1
    [fd0pos,fdmin,satelipos,vpos]=computefd_cpu_pool0401(DEM,satepos,50);%%可以用gpu加速，或者cpu并行计算
    x=fd0pos(:,:,1);y=fd0pos(:,:,2);z=fd0pos(:,:,3);
    folder=['..\input\'];
    if ~isdir(folder)
        mkdir(folder)
    end
    save([folder,'x.mat'],'x')
    save([folder,'y.mat'],'y')
    save([folder,'z.mat'],'z')
    save([folder,'vpos.mat'],'vpos')
    save([folder,'satelipos.mat'],'satelipos')
else
    %%%%%%%%%%%%%%直接读取数据%%%%%%%%%%%%%%
     load('..\input\x.mat');satex=x;
     load('..\input\y.mat');satey=y;
     load('..\input\z.mat');satez=z;
     load('..\input\vpos.mat');               
     load('..\input\satelipos.mat');satelipos=satelipos;
     fd0pos(:,:,1)=satex;fd0pos(:,:,2)=satey;fd0pos(:,:,3)=satez;
     disp('成像卫星位置计算完毕')
end
toc

%% %%%%%%%%%%%%得到LOS方向的单位矢量%%%%%%%%
LOS=computeLos(fd0pos,DEM);
folder=['..\output\SLC\'];
if ~isdir(folder)
    mkdir(folder)
end
LOS_xtemp=fd0pos(:,:,1);LOS_ytemp=fd0pos(:,:,2);LOS_ztemp=fd0pos(:,:,3);
save([folder,'LOS.mat'],'LOS_xtemp','LOS_ytemp','LOS_ztemp')%%%%%%%%无问题
%% %%%%%%%%%%%%计算每个像素的对应经纬度坐标%%%%%%%%%%
disp('计算每个像素的经纬度坐标')
ImageAzimath=4096*2;ImageRange=4096*2;
[pixeloc,NewXY]=pixe_loc(satelipos,fd0pos,DEM,long,lat,Drange,ImageAzimath,ImageRange);
disp('进行图像裁剪并输出像素坐标')
[startpix,endpix,cutImage]=find_edge(pixeloc,satelipos);
piexloc_cut=pixeloc(startpix(2):endpix(2)-1,:,startpix(1):endpix(1)-1);
pix_longtemp=piexloc_cut(:,1,:);pix_lattemp=piexloc_cut(:,2,:);
save([folder,'pix_long.mat'],'pix_longtemp');
save([folder,'pix_lat.mat'],'pix_lattemp');
save([folder,'NewY.mat'],'NewYtemp');

%%     %%%%计算入射角及散射系数%%%%
disp('计算散射系数');
sigma=computeIncidenceAngle(fd0pos,DEM);
RandAngle=rand(size(DEM,1),size(DEM,2))*2*pi;
% RandAngle=(0.6+rand(row,clomn))*2*pi;
% load('..\output\RandAngle.mat');
%RandAngle=0;

%%     %%%%生成图像%%%%%%%%%
tic;
display('生成SLC主图像');
folder=['..\output\SLC\'];
[Image1,conpl1,~,ImageSlantrange1,Imageheight,NewX1,NewY1,n1,MainPosX,MainPosY,MainPosZ,Vs_X,Vs_Y,Vs_Z,DEM_X,DEM_Y,DEM_Z]=...
generateSLC0401(C,satelipos,fd0pos,DEM,sigma,RandAngle,Wlength,Drange,DAzimath,ImageAzimath,ImageRange,vpos);%输出主图像及控制点信息
MainSLC = Image1(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
Slantrange_Main = ImageSlantrange1(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);%截图后的主天线对应斜距
Main_Height = Imageheight(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
MainPosX = MainPosX(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
MainPosY = MainPosY(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
MainPosZ = MainPosZ(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
MainPos(:,:,1)=MainPosX;
MainPos(:,:,2)=MainPosY;
MainPos(:,:,3)=MainPosZ;%存储主星位置坐标
sate_pos = MainPos;
Vs_X = Vs_X(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
Vs_Y = Vs_Y(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
Vs_Z = Vs_Z(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
DEM_X = DEM_X(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
DEM_Y = DEM_Y(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
DEM_Z = DEM_Z(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
save([folder,'Slantrange_Main.mat'],'Slantrange_Main');
save([folder,'MainSLC.mat'],'MainSLC');
figure;
imagesc(abs(MainSLC));
toc

%% %%%%%%%%%控制基线长度并产生辅图像%%%%%%%%%%%
tic;
folder=['..\output\SLC\'];
display('控制基线长度');
B_r = 1.65;
para_line = [40*B_r,40*B_r];%平行基线长度
ver_line = [200*B_r,200*B_r];%垂直基线长度300,450,750
%%%%%%%%%%求取对应的偏移矢量%%%%%%%%%%%%%%
centerpos=satelipos(floor(row/2),floor(clomn/2));
centerxyz=satepos(centerpos,1:3);
centerv=satepos(centerpos,4:end);
tempDEM=reshape(DEM(floor(row/2),floor(clomn/2),:),[1,3]);
centerR=centerxyz-tempDEM; %%%%%%减去场景中心点
para_vec=centerR./sqrt(dot(centerR,centerR));%%%%% 求取平行基线的单位矢量
ver_vec=cross(centerv,para_vec);
ver_vec=ver_vec./sqrt(dot(ver_vec,ver_vec));%%%%%  求取垂直基线单位矢量
toc

%% %%%%%%%%%%% 产生无形变的辅图像1 %%%%%%%%%%%%%%
tic;
folder=['..\output\SLC\'];
display('生成SLC辅图像');
shift_pos = para_line(1).*para_vec+ver_line(1).*ver_vec;
Baseline = sqrt(shift_pos(1)^2+shift_pos(2)^2+shift_pos(3)^2);
save([folder,'Baseline.mat'],'Baseline');
shift_postemp=zeros(size(fd0pos));
shift_postemp(:,:,1)=shift_pos(:,1);shift_postemp(:,:,2)=shift_pos(:,2);shift_postemp(:,:,3)=shift_pos(:,3);
fd0pos_slave=fd0pos+shift_postemp;

[Image2,conpl2,~,ImageSlantrange2,~,~,~,n2,SlavePosX,SlavePosY,SlavePosZ,~,~,~,~,~,~]=...
generateSLC0401(C,satelipos,fd0pos_slave,DEM,sigma,RandAngle,Wlength,Drange,DAzimath,ImageAzimath,ImageRange,vpos,NewX1,NewY1);%输出主图像及控制点信息

SlaveSLC=Image2(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
Slantrange_slave=ImageSlantrange2(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1,:);
SlavePosX=SlavePosX(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
SlavePosY=SlavePosY(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
SlavePosZ=SlavePosZ(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
SlavePos(:,:,1)=SlavePosX;
SlavePos(:,:,2)=SlavePosY;
SlavePos(:,:,3)=SlavePosZ;
sate_pos_S = SlavePos;
save([folder,'Slantrange_slave.mat'],'Slantrange_slave');
save([folder,'SlaveSLC.mat'],'SlaveSLC');
toc
figure;
imagesc(abs(SlaveSLC));
disp('图像生成完成');

%% 载入主辅图像,生成缠绕相位
folder=['..\output\SLC\'];
ssc_master = MainSLC;%主图像
ssc_slave = SlaveSLC;%短基线辅图像

%%没有去参考地形的主辅图像裁剪至相同大小
% 取主辅图像公共区域
[nr1,nc1]  = size(ssc_master);
[nr2,nc2]  = size(ssc_slave);
nr         = min(nr1,nr2);
nc         = min(nc1,nc2);
SSC_master = ssc_master(1:nr,1:nc);
SSC_slave  = ssc_slave(1:nr,1:nc);

%%配准得到干涉相位
% 粗配准
[nr1,nc1] = size(SSC_master);
[SSC_master_regis1,SSC_slave_regis1,move_r,move_c] = registration_pixel(SSC_master,SSC_slave,nr1,nc1);
[nr2,nc2] = size(SSC_master_regis1);

% 更新主图像裁剪范围
[~] = after_regis_image_index_update(move_r,move_c,process_ID_master,process_ID_slave);

% 亚像素级配准
[~,SSC_slave_regis2] = regis_subpixel(SSC_master_regis1,SSC_slave_regis1,nr2,nc2);
coherence = Calculation_Coherence_Coefficient_complex_no_phase(SSC_master_regis1,SSC_slave_regis2);
figure;imagesc(coherence);
save([folder,'coherence.mat'],'coherence','-v7.3');

% 缠绕干涉相位
wrapped_phase = interferometric_phase(SSC_master_regis1,SSC_slave_regis2);
figure;imagesc(wrapped_phase);
wrapped_phase = interferometric_phase(SSC_master_regis1,SSC_slave_regis2);
figure;imagesc(wrapped_phase);
save([folder,'wrapped_phase.mat'],'wrapped_phase','-v7.3');

% 真实相位
real_phase = -(Slantrange_Main - Slantrange_slave) / lambda * C;%因为主辅斜距差由负到正，所以真实相位也由负到正
real_phase = angle(exp(1j*real_phase));
figure;imagesc(real_phase);
save([folder,'real_phase.mat'],'real_phase','-v7.3');

%% 去参考地形
tic;
display('去参考地形');
conp1m_l(:,1:2) = conpl1(:,1:2) - startpix(1) + 1;
conp2s_l(:,1:2) = conpl2(:,1:2) - startpix(1) + 1;
conp_1m = conpl1(conp1m_l(:,1)>0&conp1m_l(:,1)<nr&conp1m_l(:,2)>0&conp1m_l(:,2)<nc,:);
conp_2s = conpl2(conp2s_l(:,1)>0&conp2s_l(:,1)<nr&conp2s_l(:,2)>0&conp2s_l(:,2)<nc,:);
conp_1_l = conp_1m;
conp_2_l = conp_2s;
conp_1_l(:,1:2) = conp_1m(:,1:2) - startpix(1) + 1;
conp_2_l(:,1:2) = conp_2s(:,1:2) - startpix(1) + 1;
clear conp1m_l conp2s_l conp_1m conp_2s

index1 = [5427 274];%找到同一方位向（行），不同距离向（列）的两个点
%需要保证两个两个控制点的行坐标是一致的，这样才能正确去平地
conpm_left = conp_1_l(index1(1),:);
conpm_right = conp_1_l(index1(2),:);
conps_left = conp_2_l(index1(1),:);
conps_right = conp_2_l(index1(2),:);

%求解去平地所需系数
delta_r1 = conpm_left(1,3) - conps_left(1,3);%控制点"LEFT"的斜距差
delta_r2 = conpm_right(1,3) - conps_right(1,3);%控制点"RIGHT"的斜距差

delta_phi1 = (delta_r1)/lambda*C;
delta_phi2 = (delta_r2)/lambda*C;
flat_a = (delta_phi2 - delta_phi1)/(size(wrapped_phase,2) - 1);%选取的2个控制点在同一方位向
flat_b = delta_phi1 - flat_a;
%计算平地相位
[nr_w,nc_w] = size(wrapped_phase);%干涉相位大小
flat_phase = zeros(nr_w,nc_w);
for jj = 1:nc_w
    flat_phase(:,jj) = flat_b + flat_a*jj;
end

%没有去参考地形的主辅图像配准得到干涉相位后再移除平地相位
%注意：缠绕干涉相位移除缠绕的平地相位后不一定在-2pi-2pi之间，因为需要再缠绕一次angle(exp())
wrapped_phase_deflat = angle(exp(1i*(wrapped_phase + (angle(exp(1i*flat_phase))))));
figure;imagesc(wrapped_phase);
figure;imagesc(angle(exp(1i*flat_phase)));
figure;imagesc(wrapped_phase_deflat);
folder=['..\output\SLC\'];
save([folder,'wrapped_phase.mat'],'wrapped_phase');
save([folder,'flat_phase.mat'],'flat_phase');
save([folder,'wrapped_phase_deflat.mat'],'wrapped_phase_deflat');
toc


%% 计算滤波前的主辅图像相干系数（只知晓一幅干涉相位图，没有主辅图像时，可以按照这种方法求相干系数）
tic;
display('计算相干系数');
[Na,Nr] = size(wrapped_phase_deflat);
coherence_1 = ones(Na,Nr);
mast_multi = ones(Na,Nr);
slave_multi = 1 * exp(1j * wrapped_phase_deflat);
coherence_1(2:Na-1,2:Nr-1) = Calculation_Coherence_Coefficient_complex_no_phase(mast_multi,slave_multi);
folder=['..\output\SLC\'];
save([folder,'coherence_1.mat'],'coherence_1','-v7.3');
toc
%% PHX:上述这里应该说是看一下相位是否连贯，与SAR的相干无关，因为去参考地形之后正常就应该是连贯的


%% MATLAB滤波
tic;
display('滤波');
folder=['..\output\SLC\'];
%斜坡自适应滤波
% Win_size  = 45;
% win_pre = 43;
% tic;
% %1——去平地后的缠绕相位滤波-%坡度自适应滤波
% wrapped_phase_short_deflat_filter = Improved_Slope_Adaptive_filter_parallel(wrapped_phase_long_deflat,Win_size,win_pre); 

%2——深度学习滤波Demo_Test_filter_phase_filter_phase2phase（较快）
%更改输入输出变量即可
format compact;
addpath('.\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab\mex');
addpath('.\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab\simplenn');
addpath('.\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab');
showResult  = 1;
useGPU      = 0;
pauseTime   = 0;

modelName   = 'model_phase_filter_phase2phase_20layer';
epoch       = 42;%鍋氬疄楠岀殑鏃跺?鏄?5

%%% load Gaussian denoising model
load(fullfile('data',modelName,[modelName,'-epoch-',num2str(epoch),'.mat']));
net = vl_simplenn_tidy(net);
net.layers = net.layers(1:end-1);

%%%
net = vl_simplenn_tidy(net);

for i = 1:1
    input = single(wrapped_phase_deflat);
    input1= cos(input);
    input2=sin(input);
    input(:,:,1,1)=input1;
    input(:,:,2,1)=input2;

     res    = vl_simplenn(net,input,[],[],'conserveMemory',true,'mode','test');
    
     output = input - res(end).x;
      wrapped_phase_deflat_filter = double(angle(output(:,:,1)+1i*output(:,:,2)));

end

figure;imagesc(wrapped_phase_deflat);
figure;imagesc(wrapped_phase_deflat_filter);
save([folder,'wrapped_phase_deflat_filter.mat'],'wrapped_phase_deflat_filter','-v7.3');
toc


%% 深度学习-C++滤波
% 用writebinary转换成Bin文件，进行深度学习滤波
% writebinary('..\SLC\wrapped_phase_long_deflat_sin.bin',sin(wrapped_phase_deflat));
% writebinary('..\SLC\wrapped_phase_long_deflat_cos.bin',cos(wrapped_phase_deflat));
% %滤波之后将Bin文件转换成mat文件
% co = readbinary('..\SLC\wrapped_phase_deflat_cos.bin.out');
% si = readbinary('..\SLC\wrapped_phase_deflat_sin.bin.out');
% co = cos(wrapped_phase_deflat) - co;
% si = sin(wrapped_phase_deflat) - si;
% wrapped_phase_deflat_filter = atan2(si, co);
% figure;imagesc(wrapped_phase_deflat_filter);
% save([folder,'wrapped_phase_deflat_filter.mat'],'wrapped_phase_deflat_filter','-v7.3');

%% 斜坡自适应-C++滤波
%用writebinary转换成Bin文件，进行斜坡自适应滤波
% writebinary('..\SLC\wrapped_phase_long_deflat.bin',wrapped_phase_deflat);
% %滤波之后将Bin文件转换成mat文件
% wrapped_phase_deflat_filter = readbinary('..\SLC\wrapped_phase_deflat_filter.bin');


%% 计算多视处理&滤波之后的主辅图像相干系数
tic;
display('计算滤波后相干系数');
folder=['..\output\SLC\'];
[Na,Nr] = size(wrapped_phase_deflat_filter);
coherence_filter = ones(Na,Nr);
mast_multi = ones(Na,Nr);
slave_multi = 1 * exp(1j * wrapped_phase_deflat_filter);
coherence_filter(2:Na-1,2:Nr-1) = Calculation_Coherence_Coefficient_complex_no_phase(mast_multi,slave_multi);
save([folder,'coherence_filter.mat'],'coherence_filter','-v7.3');
toc

%% matlab相位解缠
tic;
display('相位解缠');
folder=['..\output\SLC\'];

%质量图法解缠
im_phase_quality = PhaseDerivativeVariance_r1(wrapped_phase_deflat_filter);
unwrapped_phase_deflat_filter =QualityGuidedUnwrap2D(wrapped_phase_deflat_filter,im_phase_quality,2*pi);%质量图解缠

%最小费用流法解缠
unwrapped_phase_deflat_filter =MCF(wrapped_phase_deflat_filter);

%图割法解缠
addpath('图割法');[unwrapped_phase_deflat_filter,~,erglist] = puma_ho(wrapped_phase_deflat_filter,2);%图割法解缠
figure;imagesc(unwrapped_phase_deflat_filter);
save([folder,'unwrapped_phase_deflat_filter.mat'],'unwrapped_phase_deflat_filter','-v7.3');

toc

%% C++相位解缠
folder=['..\output\SLC\'];
%用writebinary转换成Bin文件，进行最小费用流解缠
writebinary('..\SLC\wrapped_phase_deflat_filter.bin',wrapped_phase_deflat_filter);
%解缠之后将Bin文件转换成mat文件
unwrapped_phase_deflat_filter = readbinary('..\SLC\unwrapped_phase_deflat_filter.bin');
save([folder,'unwrapped_phase_deflat_filter.mat'],'unwrapped_phase_deflat_filter','-v7.3');

%% 校正至绝对相位
tic;
display('控制点校正&标记跳变点');
folder=['..\output\SLC\'];
%加回平地相位进行控制点校正
unwrapped_phase = unwrapped_phase_deflat_filter - flat_phase;
figure;imagesc(unwrapped_phase);

%选取控制点进行校正
index2 = [1 230 360 590 620 750 880 123 456];
conp1_bk_l = conp_1_l(index2,:);
conp2_bk_l = conp_2_l(index2,:);
meanconp_h = mean(conp1_bk_l(:,6));%控制点平均高程

phi_1 = conp1_bk_l(:,3) / Wlength * C;
phi_2 = conp2_bk_l(:,3) / Wlength * C;
phi_abs = abs(phi_1 - phi_2);
[conp_nr,conp_nc] = size(phi_abs);
phase_abs = zeros(1,conp_nr);
for iii = 1:conp_nr
    phase_abs(iii) = phi_abs(iii) - unwrapped_phase(conp1_bk_l(iii, 1), conp1_bk_l(iii, 2));
end

k_point_l = round(phase_abs/(2*pi));
k_value_l = unique(k_point_l);
[~, k_value_number] = size(k_value_l);
n = 1;
for ii = 1:k_value_number
    m = sum(k_point_l==k_value_l(1,ii));
    if m > n
        n = m;
        k_real_l = k_value_l(ii);
    end
end

unwrapped_phase_abs = unwrapped_phase + 2 * pi * k_real_l;%加上了平地相位后的绝对解缠相位
figure;imagesc(unwrapped_phase_abs);
save([folder,'unwrapped_phase_abs.mat'],'unwrapped_phase_abs','-v7.3');

%绝对相位减去平地相位
unwrapped_phase_abs_deflat = unwrapped_phase_abs + flat_phase;
figure;imagesc(unwrapped_phase_abs_deflat);
save([folder,'unwrapped_phase_abs_deflat.mat'],'unwrapped_phase_abs_deflat','-v7.3');
toc

%% 相位精度评估
phase_real = -(Slantrange_Main - Slantrange_slave) / lambda * C;%因为主辅斜距差由负到正，所以真实相位也由负到正
figure;imagesc(phase_real);
save([folder,'phase_real.mat'],'phase_real','-v7.3');
phase_delta =  phase_real - unwrapped_phase_abs;%相位误差
figure;imagesc(phase_delta);
save([folder,'phase_delta.mat'],'phase_delta','-v7.3');
phase_mean_error = nanmean(abs(phase_delta(:)));%绝对误差均值
phase_relative_error = phase_delta ./ unwrapped_phase_abs;%相对误差
%相对误差标准差
phase_relative_error_std = nanstd(reshape(phase_delta,size(unwrapped_phase_abs,1)*size(unwrapped_phase_abs,2),1));


%% 真实高程
[Na_dem,Nr_dem] = size(unwrapped_phase_abs);

DEM_real=zeros(Na_dem,Nr_dem);
lat0=zeros(Na_dem,Nr_dem);
lon0=zeros(Na_dem,Nr_dem);
for ii=1:Na_dem
    for jj=1:Nr_dem
        [lat0(ii,jj),lon0(ii,jj),DEM_real(ii,jj)]=...
        xyz2ell([DEM_X(ii,jj),DEM_Y(ii,jj),DEM_Z(ii,jj)]);
    end
end
figure;mesh(DEM_real);%真实高程

    
    
