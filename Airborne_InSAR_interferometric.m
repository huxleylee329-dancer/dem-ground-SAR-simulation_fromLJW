%23所——基于干涉相位的柔性长基线估计
%20200925
close all
clc
%% 
first=1;%%%%1表示第一次运行程序
%%%%%%系统参数%%%%%
X=300; %斜距方向采样率 
Drange=2.2;%斜距方向采样间隔
DAzimath=2.2;%方位向采样间隔
Wlength=3*10^8/(1.25*10^9); %波长，1.25GHz对应波长为0.25m
Prf = 1500; %重复频率 hz
numsate=6001;

%%%%%23所系统参数%%%%%
% X=300; %斜距方向采样率 
% Drange=2.2;%斜距方向采样间隔
% DAzimath=2.2;%方位向采样间隔
% frequency = 35*10^9;
% Wlength=3*10^8/frequency; %波长
% Prf = 5000; %重复频率 hz
% numsate=6001;
%%%%%%%%%%布设地面场景点参数%%%%%%
DEM_center_ll=DEM_center_ell(1:2)/180*pi; %布设在赤道上
Dx=Drange/2;
Dy=7;
Pointspace=1.5700382546752528e-007/2;% DEM布设间隔，需要随系统参数改变

%% 读入地面场景DEM点并转换坐标%%%%%%%
load('..\input\DEMshan.mat');
height=0.4*DEMshan;
% [nr,nc]=size(height);
% tmp = mean(mean(height));
% height = ones(nr,nc)*tmp;%取平均值作为平地高程
nnn = 1;
% height = zeros(nr*nnn,nc*nnn)*tmp;%取平均值作为平地高程
[row,clomn]=size(height);
disp('DEM读入完成')
%% %%%转化为N*N*3的矩阵%%%%%%%%%%%%(根据地形大小改变参数)
long=(linspace(-nnn*2000,nnn*2000,nnn*4001)'*ones(1,nnn*4001))'*Pointspace+DEM_center_ll(1);
lat=-(linspace(-nnn*2000,nnn*2000,nnn*4001)'*ones(1,nnn*4001))*Pointspace+DEM_center_ll(2);
long=long(1:row,1:clomn);
lat=lat(1:row,1:clomn);
DEMll=cat(3,long,lat,height);
DEM=Ell2xyz(DEMll);
disp('DEM转化完成')
%% %%%%%%读入轨道参数并且进行插值%%%%%%%%%%%
path=['..\input\descending_left_airborn.csv'];
% path=['..\input\descending_left.csv'];
satepos=compute_center(DEM,Prf,numsate,path);
%% %%%%%%%计算成像位置%%%%%%%%%%%%%%%%
%%%%%%%%%%第一次执行时进行迭代处理，后续次数直接读取%%%%%%%%%%
if first==1
    [fd0pos,fdmin,satelipos]=computefd_cpu_pool(DEM,satepos,50);%%可以用gpu加速，或者cpu并行计算 %fd0v表示输出地面点对应的卫星速度
    x=fd0pos(:,:,1);y=fd0pos(:,:,2);z=fd0pos(:,:,3);
    folder=['..\input\'];
    if ~isdir(folder)
        mkdir(folder)
    end
    save([folder,'x.mat'],'x')
    save([folder,'y.mat'],'y')
    save([folder,'z.mat'],'z')
    save([folder,'satelipos.mat'],'satelipos')
else
    %%%%%%%%%%%%%%直接读取数据%%%%%%%%%%%%%%
     load('..\input\x.mat');satex=x;
     load('..\input\y.mat');satey=y;
     load('..\input\z.mat');satez=z;
     load('..\input\satelipos.mat');satelipos=satelipos;
     fd0pos(:,:,1)=satex;fd0pos(:,:,2)=satey;fd0pos(:,:,3)=satez;
     disp('成像卫星位置计算完毕')
end
% %%%%%%%%%%%%得到LOS方向的单位矢量%%%%%%%%
LOS=computeLos(fd0pos,DEM);
folder=['..\output\'];
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
NewXY(:,:,1)=NewXY(:,:,1).*(((startpix(1)<=NewXY(:,:,1))&NewXY(:,:,1)<endpix(1)).^2-1);
NewXY(:,:,2)=NewXY(:,:,2).*(((startpix(2)<=NewXY(:,:,2))&NewXY(:,:,2)<endpix(2)).^2-1);
NewXY(1,:,:)=NewXY(1,:,:)-startpix(1);
NewXY(2,:,:)=NewXY(2,:,:)-startpix(2);
NewXtemp=NewXY(:,:,1);
NewYtemp=NewXY(:,:,2);
save([folder,'NewX.mat'],'NewXtemp');
save([folder,'NewY.mat'],'NewYtemp');
%%     %%%%计算入射角及散射系数%%%%
disp('计算散射系数');
sigma=computeIncidenceAngle(fd0pos,DEM);
RandAngle=rand(size(DEM,1),size(DEM,2))*2*pi;
%%     %%%%生成图像%%%%%%%%%
disp('生成SLC图像');
%输出控制点及图像像素对应的轨道信息
[Image,conp1,~,IS1,Imageheight_Main,NewX,NewY,~,MainPosX,MainPosY,MainPosZ]=generateSLC(satelipos,fd0pos,DEM,sigma,RandAngle,Wlength,Drange,DAzimath,ImageAzimath,ImageRange);%输出主图像及控制点信息
MainSLC = Image(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
Slantrange_Main = IS1(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);%截图后的主天线对应斜距
Main_Height = Imageheight_Main(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
MainPosX = MainPosX(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
MainPosY = MainPosY(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
MainPosZ = MainPosZ(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
clear MainPos;
MainPos(:,:,1)=MainPosX;
MainPos(:,:,2)=MainPosY;
MainPos(:,:,3)=MainPosZ;%存储主星位置坐标
save([folder,'MainPos.mat'],'MainPos');
save([folder,'conp1.mat'],'conp1');
%% %%%%%%%%%控制基线长度并产生辅图像%%%%%%%%%%%
% B = 0.2;%23所基线长度
B = 16;%绝对基线长度
amplitude = 0.0/180*pi;%横滚角抖动幅度
omega = 0.00;%横滚角抖动频率
% B_amplitude = 0.004;%基线长度抖动幅度
% B_omega = 0.03;%基线长度抖动频率
phi = 0.0/180*pi;%初始相位角
Ts = 1/(2*Prf);
nn = size(height ,1);
tt = 0:2*(nn - 1);
theta_0 = 90/180*pi;
% theta_0 = 31/180*pi;%23所波束中心与基线的夹角
ver_line = B * sin(theta_0 + amplitude * sin(omega*tt + phi));%横滚角抖动
para_line = B * cos(theta_0 + amplitude * sin(omega*tt + phi));
% ver_line = (B + B_amplitude * sin(B_omega*tt + phi)) * sin(theta_0);%基线长度抖动
% para_line = (B + B_amplitude * sin(B_omega*tt + phi))* cos(theta_0);

%% %%%%%%%%求取对应的偏移矢量%%%%%%%%%%%%%%
centerpos=satelipos(floor(row/2),floor(clomn/2));
centerxyz=satepos(centerpos,1:3);
centerv=satepos(centerpos,4:end);%飞行速度
tempDEM=reshape(DEM(floor(row/2),floor(clomn/2),:),[1,3]);
centerR=centerxyz-tempDEM; %%%%%%减去场景中心点
para_vec=centerR./sqrt(dot(centerR,centerR));%%%%%  求取平行基线的单位矢量
ver_vec=cross(centerv,para_vec);%通过求速度方向和斜距（水平基线）方向的法向量得到垂直基线单位矢量
ver_vec=ver_vec./sqrt(dot(ver_vec,ver_vec));%%%%%  求取垂直基线单位矢量
%% %%%%%%%%%%% 产生无形变的辅图像1 %%%%%%%%%%%%%%
disp('产生辅图像');
para_line = para_line';
ver_line = ver_line';
Para_vec = repmat(para_vec, size(tt,2), 1);
Ver_vec = repmat(ver_vec, size(tt,2), 1);
Para_line = repmat(para_line, 1, 3);
Ver_line = repmat(ver_line, 1, 3);
shift_pos=Para_line.*Para_vec+Ver_line.*Ver_vec;
fd0pos1 = zeros(size(fd0pos));
%% 加入基线抖动
count = 0;%填充数量
ii = 1;
jj = 1;
prev = -1;
current = satelipos(1, 1);
n = 0; %基线抖动index
while(count < nn * nn)
    if(ii < nn)
        if(prev ~= current)
            n = n + 1;
        end
        iii = ii;
        jjj = jj;
        while(iii >= 1)
            count = count + 1;
            fd0pos1(iii, jjj, 1) = fd0pos(iii, jjj, 1) + shift_pos(n, 1);
            fd0pos1(iii, jjj, 2) = fd0pos(iii, jjj, 2) + shift_pos(n, 2);
            fd0pos1(iii, jjj, 3) = fd0pos(iii, jjj, 3) + shift_pos(n, 3);
            iii = iii - 1;
            jjj = jjj + 1;
        end
        prev = current;
        current = satelipos(ii + 1, 1);
        ii = ii + 1;
        continue;
    end
    if(ii >= nn)
        if(prev ~= current)
            n = n + 1;
        end
        iii = ii;
        jjj = jj;
        while(jjj <= nn)
            count = count + 1;
            fd0pos1(iii, jjj, 1) = fd0pos(iii, jjj, 1) + shift_pos(n, 1);
            fd0pos1(iii, jjj, 2) = fd0pos(iii, jjj, 2) + shift_pos(n, 2);
            fd0pos1(iii, jjj, 3) = fd0pos(iii, jjj, 3) + shift_pos(n, 3);
            iii = iii - 1;
            jjj = jjj + 1;
        end
        jj = jj + 1;
        if(count < nn * nn)
            prev = current;
            current = satelipos(ii, jj);
        end
        continue;
    end
    
end

%% 仿真辅图像
[Image2,conp2,~,IS2,Imageheight_Slave,~,~,~,SlavePosX,SlavePosY,SlavePosZ]=generateSLC(satelipos,fd0pos1,DEM,sigma,RandAngle,Wlength,Drange,DAzimath,ImageAzimath,ImageRange,NewX,NewY);%输出无形变辅图像图像及控制点信息
SlaveSLC_ideal=Image2(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);%理想情况下的辅图像
% SlaveSLC_jitter=Image2(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);%基线抖动情况下的辅图像
Slantrange_ideal=IS2(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1,:);%理想情况下，截图后的辅天线对应斜距
% Slantrange_jitter=IS2(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1,:);%基线抖动情况下，截图后的辅天线对应斜距
Slave_Height_ideal = Imageheight_Slave(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1,:);%理想情况下，截图后的辅天线对应图像高度
SlavePosX=SlavePosX(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
SlavePosY=SlavePosY(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
SlavePosZ=SlavePosZ(startpix(1):endpix(1)-1,startpix(2):endpix(2)-1);
clear SlavePos;
SlavePos(:,:,1)=SlavePosX;
SlavePos(:,:,2)=SlavePosY;
SlavePos(:,:,3)=SlavePosZ;
save([folder,'SlavePos.mat'],'SlavePos');
save([folder,'conp2.mat'],'conp2');
disp('图像生成完成');
%% 基线抖动补偿
% compensate_SlaveSLC=SlaveSLC_jitter.*exp(1j*4*pi*(Slantrange_jitter-Slantrange_ideal)/Wlength);%对加入抖动的辅图像进行相位补偿
% imagesc(angle(MainSLC.*conj(SlaveSLC_ideal)));%理想情况下的干涉条纹图
% imagesc(angle(MainSLC.*conj(SlaveSLC_jitter)));%基线抖动情况下的干涉条纹图
% imagesc(angle(MainSLC.*conj(compensate_SlaveSLC)));%基线抖动补偿后的干涉条纹图
% imshow(mapminmax(angle(MainSLC1.*conj(SlaveSLC_jitter))));%归一化
%% 输入处理ID
process_ID_master  = '500_mast';
process_ID_slave   = '500_slave';
process_ID = [process_ID_master,'&',process_ID_slave,'_insar'];
%% 载入主辅图像
ssc_master = MainSLC;%主图像
ssc_slave = SlaveSLC_ideal;%理想情况下辅图像
% ssc_slave = SlaveSLC_jitter;%抖动情况下辅图像
% ssc_slave = compensate_SlaveSLC;%基线抖动补偿后辅图像

%% 得到参考地形控制点数据
% index = [2991 2012 2985 1760 1820 ];
index = [2012 2985 1760 1820 ];
% index = [572 970 1423 1688 1880 2068 2340];
GCPs_mast1 = conp1(index,:);
GCPs_mast1(:,1:2) = GCPs_mast1(:,1:2) - startpix(1) + 1;
[GCPs_mast1(:,4),GCPs_mast1(:,5),GCPs_mast1(:,6)] = xyz2ell(GCPs_mast1(:,4:6));
GCPs_mast = zeros(size(index, 2),6);
GCPs_mast(:,1) = GCPs_mast1(:,4);
GCPs_mast(:,2) = GCPs_mast1(:,5);
GCPs_mast(:,3) = mean(GCPs_mast1(:,6));
GCPs_mast(:,4) = GCPs_mast1(:,1);
GCPs_mast(:,5) = GCPs_mast1(:,2);

tmp = zeros(size(index, 2), 1, 3);
tmp(:,:,1) = GCPs_mast1(:,5)/180*pi;%longtitude
tmp(:, :, 2) = GCPs_mast1(:,4)/180*pi;% latitude
tmp(:, :, 3) = GCPs_mast(:,3); % height
[xyz] = Ell2xyz(tmp);
x = zeros(size(index, 2), 1);
y = zeros(size(index, 2), 1);
z = zeros(size(index, 2), 1);

x(:) = xyz(:,:,1);
y(:) = xyz(:,:,2);
z(:) = xyz(:,:,3);
for ii = 1:size(index, 2)
    GCPs_mast(ii, 6) = sqrt(sum((x(ii) - MainPosX(GCPs_mast1(ii,1), GCPs_mast1(ii,2)))^2+(y(ii) - MainPosY(GCPs_mast1(ii,1), GCPs_mast1(ii,2)))^2+(z(ii) - MainPosZ(GCPs_mast1(ii,1), GCPs_mast1(ii,2)))^2));
end

GCPs_slave1 = conp2(index,:);
GCPs_slave1(:,1:2) = GCPs_slave1(:,1:2) - startpix(1) + 1;
[GCPs_slave1(:,4),GCPs_slave1(:,5),GCPs_slave1(:,6)] = xyz2ell(GCPs_slave1(:,4:6));
GCPs_slave = zeros(size(index, 2),6);
GCPs_slave(:,1) = GCPs_slave1(:,4);
GCPs_slave(:,2) = GCPs_slave1(:,5);
GCPs_slave(:,3) = mean(GCPs_slave1(:,6));
GCPs_slave(:,4) = GCPs_slave1(:,1);
GCPs_slave(:,5) = GCPs_slave1(:,2);

tmp = zeros(size(index, 2), 1, 3);
tmp(:,:,1) = GCPs_slave1(:,5)/180*pi;%longtitude
tmp(:, :, 2) = GCPs_slave1(:,4)/180*pi;% latitude
tmp(:, :, 3) = GCPs_slave(:,3); % height
[xyz] = Ell2xyz(tmp);
x = zeros(size(index, 2), 1);
y = zeros(size(index, 2), 1);
z = zeros(size(index, 2), 1);

x(:) = xyz(:,:,1);
y(:) = xyz(:,:,2);
z(:) = xyz(:,:,3);
for ii = 1:size(index, 2)
    GCPs_slave(ii, 6) = sqrt(sum((x(ii) - SlavePosX(GCPs_slave1(ii,1), GCPs_slave1(ii,2)))^2+(y(ii) - SlavePosY(GCPs_slave1(ii,1), GCPs_slave1(ii,2)))^2+(z(ii) - SlavePosZ(GCPs_slave1(ii,1), GCPs_slave1(ii,2)))^2));
end

save('GCPs_mast.mat', 'GCPs_mast');
save('GCPs_slave.mat', 'GCPs_slave');
%% 主图像去参考地形
display('主图去平地......');
[nr_flat, nc_flat] = size(ssc_master);
% 线性拟合
N_gcps = size(GCPs_mast, 1);
rows = GCPs_mast(:,4); 
cols =  GCPs_mast(:,5); 
matrix = [ones(N_gcps,1), rows, cols]; 
a_pksi_qeta = (matrix'*matrix)\matrix'*GCPs_mast(:,6);
a_pksi_qeta(2) = 0;
% err = norm(matrix*a_pksi_qeta - Rps_m);
flat_range_master = zeros(nr_flat, nc_flat);
for i = 1:nr_flat
    for j = 1:nc_flat
        flat_range_master(i,j) = [1, i, j]*a_pksi_qeta; %斜距
    end
end
flat_phase_master = flat_range_master * 4 * pi / Wlength;  %主图绝对平地相位
% 去平地效应
deflat_mast = ssc_master.*exp(-1i*flat_phase_master);
% 保存结果
save([folder, 'flat_phase_mast.mat'], 'flat_phase_master');% 保存绝对平地相位（后续加回来）
save([folder, 'deflat_mast.mat'], 'deflat_mast');% 保存去平地之后的主图像
imagesc(angle(exp(1i*flat_phase_master)));colorbar;
%% 辅图像去参考地形
display('辅图去平地......');
rows = GCPs_slave(:,4); 
cols =  GCPs_slave(:,5); 
matrix = [ones(N_gcps,1), rows, cols]; % 线性拟合
a_pksi_qeta = (matrix'*matrix)\matrix'*GCPs_slave(:,6);
a_pksi_qeta(2) = 0;
% err = norm(matrix*a_pksi_qeta - GCPs_slave(:,6));
 %根据9个控制点，利用最小二乘法拟合出方程系数，然后利用该系数求出全部像素点的平地干涉相位；      
%     a_pksi_qeta = matrix(4:end,:)\Gsflat_phi(4:end,k)
flat_range_slave = zeros(nr_flat, nc_flat);
for i = 1:nr_flat
    for j = 1:nc_flat
        %flat_int_phase(ii,jj,k) = [1, ii, jj, ii*jj, ii*ii, jj*jj]*a_pksi_qeta;
        flat_range_slave(i,j) = [1, i, j]*a_pksi_qeta; %斜距
    end
end
flat_phase_slave = flat_range_slave * 4 * pi / Wlength;  %主图平地相位
%去平地效应
deflat_slave = ssc_slave.*exp(-1i*flat_phase_slave);
save([folder, 'deflat_slave.mat'], 'deflat_slave');% 保存去平地之后的辅图像
save([folder, 'flat_phase_slave.mat'], 'flat_phase_slave');% 保存绝对平地相位
imagesc(angle(exp(1i*flat_phase_slave)));colorbar;
display('去平地完成！');
%% 裁剪至相同大小
% 取主辅图像公共区域
[nr1,nc1]  = size(deflat_mast);
[nr2,nc2]  = size(deflat_slave);
nr         = min(nr1,nr2);
nc         = min(nc1,nc2);
SSC_master = deflat_mast(1:nr,1:nc);
SSC_slave  = deflat_slave(1:nr,1:nc);
% 更新主图像裁剪范围
[~] = cut_image_index_update(nr1,nc1,nr,nc,process_ID_master);
[~] = cut_image_index_update(nr2,nc2,nr,nc,process_ID_slave);
%% 干涉处理
% 粗配准
[nr1,nc1] = size(ssc_master);
[SSC_master_regis1,SSC_slave_regis1,move_r,move_c] = registration_pixel(SSC_master,SSC_slave,nr1,nc1);
[nr2,nc2] = size(SSC_master_regis1);

% 更新主图像裁剪范围
[~] = after_regis_image_index_update(move_r,move_c,process_ID_master,process_ID_slave);

% 亚像素级配准
[~,SSC_slave_regis2] = regis_subpixel(SSC_master_regis1,SSC_slave_regis1,nr2,nc2);
coherence = Calculation_Coherence_Coefficient_complex_no_phase(SSC_master_regis1,SSC_slave_regis2);
folder=['..\Mid\',process_ID,'\'];
save([folder,'coherence.mat'],'coherence');
save([folder,'SSC_master_regis1.mat'],'SSC_master_regis1');
save([folder,'SSC_slave_regis2.mat'],'SSC_slave_regis2');
% 缠绕干涉相位
wrapped_phase = interferometric_phase(SSC_master_regis1,SSC_slave_regis2);
save([folder,'wrapped_phase.mat'],'wrapped_phase');
%% 干涉相位滤波
   win_size  = 9;
     win_pre = 7;
     tic;
     phase_removeFlat_filter = Improved_Slope_Adaptive_filter_parallel(wrapped_phase,win_size,win_pre); %坡度自适应滤波
     save ([folder,'phase_removeFlat_filter.mat'], 'phase_removeFlat_filter');    %保存残余相位
     wrapped_phase_filter_gray = mat2gray(phase_removeFlat_filter); 
     imwrite(wrapped_phase_filter_gray,[folder,'phase_removeFlat_filter.tif']);
     toc

 %% 相位解缠（最小费用流法）
[unwrapped_phase_removeFlat]=MCF(phase_removeFlat_filter);%残余相位解缠

unwrapped_phase= -unwrapped_phase_removeFlat+(flat_phase_slave-flat_phase_master)...
         -phase_removeFlat_filter(1,1);%去参考地形得到的主辅绝对平地相位补偿
%% 校正至绝对相位
conp_1=conp1;
conp_2=conp2;
[k,~] = size(conp_1);
del = ones(k, 2)*startpix(1);

conp_1(:,1:2) = conp_1(:,1:2) - del;
conp_2(:,1:2) = conp_2(:,1:2) - del;

count = 0;
for iii = 1:k
    if (conp_1(iii, 1) > 0 && conp_1(iii, 1) < nr1+1) && (conp_1(iii, 2) > 0 && conp_1(iii, 2) < nr1+1)
        count = count + 1;
    end
end
index = zeros(count, 1);
jj = 0;
for iii = 1:k
    if (conp_1(iii, 1) > 0 && conp_1(iii, 1) < nr1+1) && (conp_1(iii, 2) > 0 && conp_1(iii, 2) < nr1+1)
        jj = jj + 1;
        index(jj) = iii;
    end
end
conp1_bk = conp_1(index,:);
conp2_bk = conp_2(index,:);
phi_1 = conp1_bk(:,3) / Wlength * 4 * pi;
phi_2 = conp2_bk(:,3) / Wlength * 4 * pi;
phi_abs = abs(phi_1 - phi_2);
phase_abs = zeros(1,count);
for iii = 1:count
    phase_abs(iii) = phi_abs(iii) - unwrapped_phase(conp1_bk(iii, 1), conp1_bk(iii, 2));
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

unwrapped_phase_abs = unwrapped_phase + 2 * pi * k_real;

%% 求相位误差
phase_real = Slantrange_Main/Wlength*4*pi - Slantrange_ideal/Wlength*4*pi;%因为主辅斜距差由负到正，所以真实相位也由负到正
unwrapped_phase_abs = unwrapped_phase_abs + mean(mean(phase_real - unwrapped_phase_abs));
delta_high =  unwrapped_phase_abs-phase_real;%相位误差
save([folder, 'unwrapped_phase_abs.mat'], 'unwrapped_phase_abs');%保存绝对相位
save([folder, 'delta_high.mat'], 'delta_high');%保存相位误差

% fig2 = figure;
% set(fig2, 'visible', 'off');
% imagesc(unwrapped_phase_abs);colorbar;
% saveas(fig2, [folder, 'unwrapped_phase_abs.jpg']);
%% 高程反演1
Re = 6378136.49; %地球半径
S_position_m_cut = MainPos;%主星位置
S_position_s_cut = SlavePos;%辅星位置
Rh = sqrt(MainPos(:,:,1).^2+MainPos(:,:,2).^2+MainPos(:,:,3).^2);%主星高度
h=GCPs_mast1(4,6);%主星控制点
h_out = 0;
iter = 30;%迭代次数
lambda =Wlength;%波长
Real_phi = unwrapped_phase_abs - (flat_phase_master-flat_phase_slave);%校正后的真实相位，不包括平地相位
DEM_unwrapped_phase = -(unwrapped_phase_removeFlat+phase_removeFlat_filter(1,1));%解缠后的相位

Baselinespan = ones(size(DEM_unwrapped_phase,1)).*B;%基线长度
alpha = theta_0-(pi/2-45/180*pi);%基线倾角=波束中心与基线的夹角-(pi/2-波束中心视角）
Bh = Baselinespan .* cos(alpha);%水平基线
Bv = Baselinespan .* sin(alpha);%垂直基线

Reverse_DEM = zeros(size(DEM_unwrapped_phase));
theta_total = zeros(size(DEM_unwrapped_phase));
for ii = 1:size(DEM_unwrapped_phase,1)
    for jj = 1:size(DEM_unwrapped_phase,2)
        phase = Real_phi(ii,jj);
        Rm = Slantrange_Main(ii,jj);    
        for k = 1:iter
            theta = acos((Rh(ii)^2+Rm^2-(h_out+Re)^2)/2/Rm/Rh(ii));
            h_last = h_out;
            h_out = lambda/4/pi * Rm * Rh(ii) * sin(theta) * phase / (Bh(ii)*cos(theta) + Bv(ii)*sin(theta)) / sqrt(Rh(ii)^2+Rm^2-2*Rm*Rh(ii)*cos(theta)) ;
%             h_out = lambda/4/pi * Rm * sin(theta) * phase / (Bh(ii)*cos(theta) + Bv(ii)*sin(theta)) ;
            if abs(h_last - h_out) <= 10e-6
                Reverse_DEM(ii,jj) = h_out;
                theta_total(ii,jj) = theta;
                break
            end
        end
    end
end
delt_h = Reverse_DEM(GCPs_mast1(4,1),GCPs_mast1(4,2)) + 787.6408 - h;
DEMnew = Reverse_DEM - delt_h + 787.6408 ;
DEM_error = Main_Height -Re - DEMnew;%绝对误差
DEM_relative_error = (Main_Height -Re - DEMnew)./DEMnew;%相对误差
DEM_std = std(reshape(Main_Height -Re - DEMnew,size(DEMnew,1)*size(DEMnew,2),1));
figure;mesh(DEMnew)
% clf;
% imagesc(DEMnew);
% colormap(jet);
% axis image; 
% c=colorbar;
% set(get(c,'title'),'string','[m]');
% t=title ('DEM');
% set(t,'fontweight','bold');
% xlabel ('range pixels');
% ylabel ('azimuth lines');
%% 高程反演2
Slantrange_s = Slantrange_ideal;
Slantrange_m = Slantrange_Main;
R_M = Slantrange_Main;
pos = csvread(path);
Vs = mean(pos(:, 5:7))*1000;
Satellite_M_R_Velocity = repmat(Vs, size(R_M, 1), 1);
Satellite_M_T_Position = zeros(size(R_M, 1), 3);
Satellite_M_T_Position(:, 1) = mean(MainPosX, 2);
Satellite_M_T_Position(:, 2) = mean(MainPosY, 2);
Satellite_M_T_Position(:, 3) = mean(MainPosZ, 2);

Satellite_S_T_Position = zeros(size(R_M, 1), 3);
Satellite_S_T_Position(:, 1) = mean(SlavePosX, 2);
Satellite_S_T_Position(:, 2) = mean(SlavePosY, 2);
Satellite_S_T_Position(:, 3) = mean(SlavePosZ, 2);

Rmin = Slantrange_Main(2,2);
c = 3e8;

Satellite_M_R_Position = Satellite_M_T_Position+2*(Rmin)/c*Satellite_M_R_Velocity;
Satellite_S_R_Position = Satellite_S_T_Position+2*(Rmin)/c*Satellite_M_R_Velocity;

Satellite_M=(Satellite_M_R_Position+Satellite_M_T_Position)./2;

Real_phi=unwrapped_phase_abs(2:end - 1,2:end - 1);
[Cut_Na,Cut_Nr]=size(Real_phi);
fd=ones(1,Cut_Nr)*(-0.278713789746405)*0; %多普勒中心频率

R_M = R_M(2:end - 1,2:end - 1);
Satellite_M_R_Velocity = Satellite_M_R_Velocity(2:end - 1,:);
Satellite_M_T_Position = Satellite_M_T_Position(2:end - 1,:);
Satellite_S_T_Position = Satellite_S_T_Position(2:end - 1,:);
Satellite_M_R_Position = Satellite_M_R_Position(2:end - 1,:);
Satellite_S_R_Position = Satellite_S_R_Position(2:end - 1,:);
Satellite_M = Satellite_M(2:end - 1,:);

% 控制点
P=conp1(1565,4:6);
%P = [6378932.39735364; -1007.02513948654; 1014.70937408723];

m=30;%迭代50次
Vs=Satellite_M_R_Velocity/1000;
R_F = 2 * R_M  - Real_phi * Wlength / (2 * pi);

Df=zeros(3,3);
Df_ni=zeros(3,3);
f=zeros(3,1);
delta_Rt=zeros(3,1); 
P1=ones(Cut_Na,Cut_Nr)*P(1);
P2=ones(Cut_Na,Cut_Nr)*P(2);
P3=ones(Cut_Na,Cut_Nr)*P(3);
clear Real_phi B_Vetor Satellite_M_R_Velocity x y z unwrapped_phase ii jj RanStartTime

% tic
for k=0:m-1
%     tic
%     k
%--------------------------------------------------    
% f1=(Satellite_M(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_M(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_M(:,3)*ones(1,Cut_Nr)- P3).^2- R_M.^2;
M_T=(Satellite_M_T_Position(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_M_T_Position(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_M_T_Position(:,3)*ones(1,Cut_Nr)- P3).^2;
M_R=(Satellite_M_R_Position(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_M_R_Position(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_M_R_Position(:,3)*ones(1,Cut_Nr)- P3).^2;
f1=sqrt(M_T)+sqrt(M_R)- 2*R_M;
%--------------------------------------------------

S_T=(Satellite_S_T_Position(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_S_T_Position(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_S_T_Position(:,3)*ones(1,Cut_Nr)- P3).^2;
S_R=(Satellite_S_R_Position(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_S_R_Position(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_S_R_Position(:,3)*ones(1,Cut_Nr)- P3).^2;
f2=sqrt(S_T)+sqrt(S_R)- R_F;

f3=Vs(:,1)*ones(1,Cut_Nr).*(Satellite_M(:,1)*ones(1,Cut_Nr)-P1)+Vs(:,2)*ones(1,Cut_Nr).*(Satellite_M(:,2)*ones(1,Cut_Nr)-P2)+Vs(:,3)*ones(1,Cut_Nr).*(Satellite_M(:,3)*ones(1,Cut_Nr)-P3)-ones(Cut_Na,1)*fd.*R_M*Wlength/2.0;
% temp=sqrt(f1.^2+f2.^2+f3.^2);
%% Dff
%第一行：f(1)的x，y，z的导数
%--------------------------------------------------   

Df11=1/2*(M_T).^(-1/2)*(-2.0).*(Satellite_M_T_Position(:,1)*ones(1,Cut_Nr)-P1)...
+1/2*(M_R).^(-1/2)*(-2.0).*(Satellite_M_R_Position(:,1)*ones(1,Cut_Nr)-P1);

Df12=1/2*(M_T).^(-1/2)*(-2.0).*(Satellite_M_T_Position(:,2)*ones(1,Cut_Nr)-P2)...
+1/2*(M_R).^(-1/2)*(-2.0).*(Satellite_M_R_Position(:,2)*ones(1,Cut_Nr)-P2);

Df13=1/2*(M_T).^(-1/2)*(-2.0).*(Satellite_M_T_Position(:,3)*ones(1,Cut_Nr)-P3)...
+1/2*(M_R).^(-1/2)*(-2.0).*(Satellite_M_R_Position(:,3)*ones(1,Cut_Nr)-P3);
%-------------------------------------------------- 
clear M_T M_R
%第二行：f(2)的x，y，z的导数
Df21=1/2*(S_T).^(-1/2)*(-2.0).*(Satellite_S_T_Position(:,1)*ones(1,Cut_Nr)-P1)...
+1/2*(S_R).^(-1/2)*(-2.0).*(Satellite_S_R_Position(:,1)*ones(1,Cut_Nr)-P1);

Df22=1/2*(S_T).^(-1/2)*(-2.0).*(Satellite_S_T_Position(:,2)*ones(1,Cut_Nr)-P2)...
+1/2*(S_R).^(-1/2)*(-2.0).*(Satellite_S_R_Position(:,2)*ones(1,Cut_Nr)-P2);

Df23=1/2*(S_T).^(-1/2)*(-2.0).*(Satellite_S_T_Position(:,3)*ones(1,Cut_Nr)-P3)...
+1/2*(S_R).^(-1/2)*(-2.0).*(Satellite_S_R_Position(:,3)*ones(1,Cut_Nr)-P3);

clear S_T S_R
%第三行：f(3)的x，y，z的导数
Df31=-Vs(:,1)*ones(1,Cut_Nr);
Df32=-Vs(:,2)*ones(1,Cut_Nr);
Df33=-Vs(:,3)*ones(1,Cut_Nr);

det_Df=Df11.*Df22.*Df33+Df12.*Df23.*Df31+Df13.*Df21.*Df32-Df31.*Df22.*Df13-Df32.*Df23.*Df11-Df33.*Df21.*Df12;
% temp=temp+det_Df;%求Df行列式且与temp相加
Df_ni11=(Df22.*Df33-Df32.*Df23)./det_Df;
Df_ni12=-(Df12.*Df33-Df32.*Df13)./det_Df;
Df_ni13=(Df12.*Df23-Df22.*Df13)./det_Df;
delta_Rt1=Df_ni11.*f1+Df_ni12.*f2+Df_ni13.*f3;
clear  Df_ni11 Df_ni12 Df_ni13
Df_ni21=-(Df21.*Df33-Df31.*Df23)./det_Df;
Df_ni22=(Df11.*Df33-Df31.*Df13)./det_Df;
Df_ni23=-(Df11.*Df23-Df21.*Df13)./det_Df;
delta_Rt2=Df_ni21.*f1+Df_ni22.*f2+Df_ni23.*f3;
clear   Df_ni21 Df_ni22 Df_ni23
Df_ni31=(Df21.*Df32-Df31.*Df22)./det_Df;
Df_ni32=-(Df11.*Df32-Df31.*Df12)./det_Df;
Df_ni33=(Df22.*Df11-Df21.*Df12)./det_Df;
delta_Rt3=Df_ni31.*f1+Df_ni32.*f2+Df_ni33.*f3;
clear det_Df Df11 Df12 Df13 Df21 Df22 Df23 Df31 Df32 Df33 Df_ni31 Df_ni32 Df_ni33 f1 f2 f3

% temp=sqrt((delta_Rt1.^2+delta_Rt2.^2+delta_Rt3.^2)./(P1.^2+P2.^2+P3.^2));
P1=P1-delta_Rt1;
P2=P2-delta_Rt2;
P3=P3-delta_Rt3;
clear delta_Rt1 delta_Rt2 delta_Rt3

% toc
end
Target_84=zeros(Cut_Na,Cut_Nr,3);
Target_84(:,:,1)=P1;
Target_84(:,:,2)=P2;
Target_84(:,:,3)=P3;
clear P1 P2 P3
DEM_x=Target_84(:,:,1);
DEM_y=Target_84(:,:,2);
DEM_z=Target_84(:,:,3);
Rad_earth_e = 6378136.49;
DEM_inver=sqrt(DEM_x.^2+DEM_y.^2+DEM_z.^2)-Rad_earth_e;
figure;imagesc(DEM_inver);colorbar
figure;mesh(DEM_inver);colorbar




