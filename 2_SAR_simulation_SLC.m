
%本程序是星载SAR图像仿真程序，整理于2021.10.09
%设置系统参数生成主辅图像，干涉得到缠绕相位、真实相位和相干系数
clear all
close all
clc
first=0;%%%%表示第一次运行程序
%%%%%%卫星轨道参数%%%%%
Drange=2.2711549848484847;%斜距方向采样间隔
DAzimath=2.3784708055877464;
Wlength=3*10^8/(9.6*10^9);
lambda = Wlength;
Prf=3179; %重复频率 hz
C = 4*pi;
%%%%%%%%%%布设地面场景点参数%%%%%%
DEM_center_ll=[0;0]/180*pi; %布设在赤道上

Pointspace=1.5700382546752528e-007/2;% DEM布设间隔，LS间隔为0.5m

%%%%%%%%%%轨道参数%%%%%%%%%%%%%%%
startT=838.357397837519100;
timeSpace=1/2/Prf;
numsate=4401;
%% 读入地面场景DEM点并转换坐标%%%%%%%
load('..\input\DEMshan.mat');
height = DEMshan;
% load('..\input\DEM_lan.mat');
% height=DEM_lan;
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
NewXY(:,:,1)=NewXY(:,:,1).*(((startpix(1)<=NewXY(:,:,1))&NewXY(:,:,1)<endpix(1)).^2-1);
NewXY(:,:,2)=NewXY(:,:,2).*(((startpix(2)<=NewXY(:,:,2))&NewXY(:,:,2)<endpix(2)).^2-1);
NewXY(1,:,:)=NewXY(1,:,:)-startpix(1);
NewXY(2,:,:)=NewXY(2,:,:)-startpix(2);
NewXtemp=NewXY(:,:,1);NewYtemp=NewXY(:,:,2);
save([folder,'NewX.mat'],'NewXtemp');
save([folder,'NewY.mat'],'NewYtemp');

%%     %%%%计算入射角及散射系数%%%%
disp('计算散射系数');
sigma=computeIncidenceAngle(fd0pos,DEM);
RandAngle=rand(size(DEM,1),size(DEM,2))*2*pi;

%%     %%%%生成图像%%%%%%%%%
tic;
display('生成SLC主图像');
folder=['..\output\SLC\'];
[Image1,conpl1,~,ImageSlantrange1,Imageheight,NewX1,NewY1,n1,MainPosX,MainPosY,MainPosZ,Vs_X,Vs_Y,Vs_Z,DEM_X,DEM_Y,DEM_Z]=...
generateSLC0401(satelipos,fd0pos,DEM,sigma,RandAngle,Wlength,Drange,DAzimath,ImageAzimath,ImageRange,vpos);%输出主图像及控制点信息
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
display('控制基线长度');
B_r = 2.5;
para_line = [20*B_r,20*B_r];%平行基线长度
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
shift_pos =para_line(1).*para_vec+ver_line(1).*ver_vec;
Baseline = sqrt(shift_pos(1)^2+shift_pos(2)^2+shift_pos(3)^2);
save([folder,'Baseline.mat'],'Baseline');
shift_postemp=zeros(size(fd0pos));
shift_postemp(:,:,1)=shift_pos(:,1);shift_postemp(:,:,2)=shift_pos(:,2);shift_postemp(:,:,3)=shift_pos(:,3);
fd0pos_slave=fd0pos+shift_postemp;

[Image2,conpl2,~,ImageSlantrange2,~,~,~,n2,SlavePosX,SlavePosY,SlavePosZ,~,~,~,~,~,~]=...
generateSLC0401(satelipos,fd0pos_slave,DEM,sigma,RandAngle,Wlength,Drange,DAzimath,ImageAzimath,ImageRange,vpos,NewX1,NewY1);%输出主图像及控制点信息

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


%% 生成缠绕相位
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

% 亚像素级配准，保存相干系数
[~,SSC_slave_regis2] = regis_subpixel(SSC_master_regis1,SSC_slave_regis1,nr2,nc2);
coherence = Calculation_Coherence_Coefficient_complex_no_phase(SSC_master_regis1,SSC_slave_regis2);
figure;imagesc(coherence);
save([folder,'coherence.mat'],'coherence','-v7.3');

% 缠绕干涉相位
wrapped_phase = interferometric_phase(SSC_master_regis1,SSC_slave_regis2);
figure;imagesc(wrapped_phase);
save([folder,'wrapped_phase.mat'],'wrapped_phase','-v7.3');

% 真实相位
real_phase = -(Slantrange_Main - Slantrange_slave) / lambda * C;%因为主辅斜距差由负到正，所以真实相位也由负到正
real_phase = angle(exp(1j*real_phase));
figure;imagesc(real_phase);
save([folder,'real_phase.mat'],'real_phase','-v7.3');









