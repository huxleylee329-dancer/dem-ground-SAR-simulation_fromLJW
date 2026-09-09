%% 去平地
addpath('./Processing functions')
%获得卫星轨道位置/速度/时刻离散值
%主图像
file_ID     = '.\';
product_ID  = 'TH02-01B_InSAR_20191008_0000001927_001_004_008_L1B.Master.meta';
xmldoc = xml_read([file_ID,product_ID,'.xml']); %信息文件
c = 3*10e8;
prf = xmldoc.sensor.waveParams.prf; %脉冲重复频率
lambda = c/(xmldoc.carrierFrequency*10e9); %波长
GPSParam =  xmldoc.GPS.GPSParam; % 轨道参数

Satellite_time = cell(size(GPSParam,1),1); %轨道位置时刻
t = zeros(size(GPSParam,1),1);
for i = 1:size(GPSParam,1)
    Satellite_time{i} = xmldoc.GPS.GPSParam(i).TimeStamp;
    t(i) = datenum(strrep(Satellite_time{i},'T',' '));
end

orbit = zeros(size(GPSParam,1),7); %轨道位置 速度
for i = 1:size(GPSParam,1)
    orbit(i,1) = t(i);
    orbit(i,2) = xmldoc.GPS.GPSParam(i).xPosition; 
    orbit(i,3) = xmldoc.GPS.GPSParam(i).yPosition;
    orbit(i,4) = xmldoc.GPS.GPSParam(i).zPosition;
    orbit(i,5) = xmldoc.GPS.GPSParam(i).xVelocity;
    orbit(i,6) = xmldoc.GPS.GPSParam(i).yVelocity;
    orbit(i,7) = xmldoc.GPS.GPSParam(i).zVelocity;
end

coef = Polyfit_Orbit(orbit); % 轨道位置与时刻的函数参数

Ng = 5; %控制点数目
GCPs = zeros(Ng,5); %控制点信息
GCPs(1,:) = [xmldoc.imageInfo.corner.topLeft.latitude,xmldoc.imageInfo.corner.topLeft.longitude,0,1,1]; %lat纬度，lon经度，height高度，row行，col列
GCPs(2,:) = [xmldoc.imageInfo.corner.topRight.latitude,xmldoc.imageInfo.corner.topRight.longitude,0,1,xmldoc.imageInfo.width]; %lat纬度，lon经度，height高度，row行，col列
GCPs(3,:) = [xmldoc.imageInfo.corner.bottomLeft.latitude,xmldoc.imageInfo.corner.bottomLeft.longitude,0,xmldoc.imageInfo.height,1]; %lat纬度，lon经度，height高度，row行，col列
GCPs(4,:) = [xmldoc.imageInfo.corner.bottomRight.latitude,xmldoc.imageInfo.corner.bottomRight.longitude,0,xmldoc.imageInfo.height,xmldoc.imageInfo.width]; %lat纬度，lon经度，height高度，row行，col列
GCPs(5,:) = [xmldoc.imageInfo.center.latitude,xmldoc.imageInfo.center.longitude,0,xmldoc.imageInfo.height/2,xmldoc.imageInfo.width/2]; %lat纬度，lon经度，height高度，row行，col列

%根据控制点信息拟合平面相位
Rps = zeros(Ng,1);
for k = 1:Ng
    [px,py,pz] = ell2xyz(GCPs(k,1),GCPs(k,2),0); %控制点的wgs84坐标
    azitime = xyz2aztime_NEWTON(xmldoc,coef_m,[px,py,pz],orbit(1,1)); %根据控制点坐标计算所对应的方位向时刻
    S_xyz = getXyz(azitime,coef_m); %根据方位向时刻计算所对应的轨道位置
    Rps(k) = sqrt((px-S_xyz(1))^2 + (py-S_xyz(2))^2 + (pz-S_xyz(3))^2); %计算卫星与控制点之间的距离
end
%记录图像中心点到卫星的距离，以及中心点成像时刻卫星的位置
Range_sat_ctar = zeros(2,1);
Range_sat_ctar(1) = Rps(5); %中心点目标所对应的斜距
Sat_position = cell(2,1);
Sat_position{1} = [px, py, pz];
% matrix = [ones(Ng,1), trows, tcols, trows.*tcols, trows.*trows, tcols.*tcols];
rows = GCPs(:,4); cols =  GCPs(:,5); 
matrix = [ones(Ng,1), rows, cols]; % 线性拟合
a_pksi_qeta = (matrix'*matrix)\matrix'*Rps;
err = norm(matrix*a_pksi_qeta - Rps);
 %根据9个控制点，利用最小二乘法拟合出方程系数，然后利用该系数求出全部像素点的平地干涉相位；      
%     a_pksi_qeta = matrix(4:end,:)\Gsflat_phi(4:end,k)
flat_range_master = zeros(xmldoc.imageInfo.height,xmldoc.imageInfo.width);
for i = 1:xmldoc.imageInfo.height
    for j = 1:xmldoc.imageInfo.width
        %flat_int_phase(ii,jj,k) = [1, ii, jj, ii*jj, ii*ii, jj*jj]*a_pksi_qeta;
        flat_range_master(i,j) = [1, i, j]*a_pksi_qeta; %斜距
    end
end
save('flat_range_master2.mat','flat_range_master','-v7.3'); 



%辅图像
file_ID  = '.\';
product_ID  = 'TH02-01A_InSAR_20191008_0000001920_001_004_008_L1B.Slave.meta';
xmldoc = xml_read([file_ID,product_ID,'.xml']);
c = 3*10e8;
lambda = c/(xmldoc.carrierFrequency*10e9); %波长
GPSParam =  xmldoc.GPS.GPSParam;

Satellite_time = cell(size(GPSParam,1),1);
t = zeros(size(GPSParam,1),1);
for i = 1:size(GPSParam,1)
    Satellite_time{i} = xmldoc.GPS.GPSParam(i).TimeStamp;
    t(i) = datenum(strrep(Satellite_time{i},'T',' '));
end

orbit = zeros(size(GPSParam,1),7);
for i = 1:size(GPSParam,1)
    orbit(i,1) = t(i);
    orbit(i,2) = xmldoc.GPS.GPSParam(i).xPosition;
    orbit(i,3) = xmldoc.GPS.GPSParam(i).yPosition;
    orbit(i,4) = xmldoc.GPS.GPSParam(i).zPosition;
    orbit(i,5) = xmldoc.GPS.GPSParam(i).xVelocity;
    orbit(i,6) = xmldoc.GPS.GPSParam(i).yVelocity;
    orbit(i,7) = xmldoc.GPS.GPSParam(i).zVelocity;
end

coef_s = Polyfit_Orbit(orbit);
Ng = 5; %控制点数目
GCPs = zeros(Ng,5);
GCPs(1,:) = [xmldoc.imageInfo.corner.topLeft.latitude,xmldoc.imageInfo.corner.topLeft.longitude,0,1,1]; %lat纬度，lon经度，height高度，row行，col列
GCPs(2,:) = [xmldoc.imageInfo.corner.topRight.latitude,xmldoc.imageInfo.corner.topRight.longitude,0,1,xmldoc.imageInfo.width]; %lat纬度，lon经度，height高度，row行，col列
GCPs(3,:) = [xmldoc.imageInfo.corner.bottomLeft.latitude,xmldoc.imageInfo.corner.bottomLeft.longitude,0,xmldoc.imageInfo.height,1]; %lat纬度，lon经度，height高度，row行，col列
GCPs(4,:) = [xmldoc.imageInfo.corner.bottomRight.latitude,xmldoc.imageInfo.corner.bottomRight.longitude,0,xmldoc.imageInfo.height,xmldoc.imageInfo.width]; %lat纬度，lon经度，height高度，row行，col列
GCPs(5,:) = [xmldoc.imageInfo.center.latitude,xmldoc.imageInfo.center.longitude,0,floor(xmldoc.imageInfo.height/2),xmldoc.imageInfo.width/2]; %lat纬度，lon经度，height高度，row行，col列

Rps = zeros(Ng,1);
for k = 1:Ng
    [px,py,pz] = ell2xyz(GCPs(k,1),GCPs(k,2),0);
    azitime = xyz2aztime_NEWTON(xmldoc,coef_s,[px,py,pz],orbit(1,1));
    S_xyz = getXyz(azitime,coef_s);
    Rps(k) = sqrt((px-S_xyz(1))^2 + (py-S_xyz(2))^2 + (pz-S_xyz(3))^2);
end
%记录图像中心点到卫星的距离，以及中心点成像时刻卫星的位置
Range_sat_ctar(2) = Rps(5);
Sat_position{2} = [px, py, pz];

% matrix = [ones(Ng,1), trows, tcols, trows.*tcols, trows.*trows, tcols.*tcols];
rows = GCPs(:,4); cols =  GCPs(:,5); 
matrix = [ones(Ng,1), rows, cols];
a_pksi_qeta = (matrix'*matrix)\matrix'*Rps;
err = norm(matrix*a_pksi_qeta - Rps);
 %根据9个控制点，利用最小二乘法拟合出方程系数，然后利用该系数求出全部像素点的平地干涉相位；      
%     a_pksi_qeta = matrix(4:end,:)\Gsflat_phi(4:end,k)
flat_range_slave = zeros(xmldoc.imageInfo.height,xmldoc.imageInfo.width);
for i = 1:xmldoc.imageInfo.height
    for j = 1:xmldoc.imageInfo.width
        %flat_int_phase(ii,jj,k) = [1, ii, jj, ii*jj, ii*ii, jj*jj]*a_pksi_qeta;
        flat_range_slave(i,j) = [1, i, j]*a_pksi_qeta;
    end
end
save('flat_range_slave2.mat','flat_range_slave','-v7.3');



%计算干涉相位
%两种方法，
%% 第一种：从主图中去掉平地相位，从辅图中去掉平地相位，然后再进行图像配准、干涉，适用于自发自收卫星系统
%读取图像 
imageA = load('D:\luoyao\推介会分发数据\嵩山定标场\TH02-01BA_InSAR_20191008_0000001927_001_004_008_L1B\SLC_Master.mat');
imageB = load('D:\luoyao\推介会分发数据\嵩山定标场\TH02-01BA_InSAR_20191008_0000001927_001_004_008_L1B\SLC_Slave.mat');
flat_phase_master = flat_range_master * 4*pi / lambda;  %主图平地相位
flat_phase_slave = flat_range_slave * 4*pi / lambda; %辅图平地相位
%裁剪部分图像
img_rows = [10089:11089];
img_cols = [9075:10075];
Master = imageA.Master(img_rows,img_cols);
Slave = imageB.Slave(img_rows,img_cols);
flat_phase_Master = flat_phase_master(img_rows,img_cols);
flat_phase_Slave = flat_phase_slave(img_rows,img_cols);
%去平地效应
deflat_Master = Master.*exp(-1i*flat_phase_Master);
deflat_Slave = Slave.*exp(-1i*flat_phase_Slave);
coa_mast_ind = 1;
coa_slave_ind = 2;
%粗配准
[move_row_part,move_col_part] = deal(zeros(1,2),zeros(1,2));
[move_row_part(coa_slave_ind),move_col_part(coa_slave_ind),~] = real_coherent(deflat_Master,deflat_Slave,0);
Sar_cut_images = cell(2,1);
Sar_cut_images{coa_mast_ind} = deflat_Master;
Sar_cut_images{coa_slave_ind} = deflat_Slave;
[image_coarse_regis,first_pixel_row,first_pixel_col] = cut_multi_images_to_regis(move_row_part,move_col_part,Sar_cut_images,coa_mast_ind);
[fine_size_rows,fine_size_cols] = size(image_coarse_regis{coa_mast_ind});

%精配准
image_fine_regis = cell(1,2);
image_fine_regis{coa_mast_ind} = image_coarse_regis{coa_mast_ind};
image_fine_regis{coa_slave_ind} = registration_subpixel1(image_coarse_regis{coa_mast_ind},image_coarse_regis{coa_slave_ind});
image_fine = image_quantify(image_fine_regis{coa_slave_ind},60);
image_fine_gray = mat2gray(image_fine);
figure;   imagesc(image_fine_gray);    title('fineregis SAR image');  drawnow;
%干涉
diff_interf_phase = zeros(fine_size_rows,fine_size_cols);

%interf_ref_terrain_m = flat_int_phase(:,:,m);
diff_interf_phase(:,:) = angle(image_fine_regis{coa_mast_ind}.*...
    conj(image_fine_regis{coa_slave_ind})); %.*exp(complex(0, interf_ref_terrain_m)));
%figure;   imagesc(interf_ref_terrain_m); title('flat earth interferometric phase');  drawnow;
figure;   imagesc(diff_interf_phase(:,:));  title('interferogram phase');  drawnow;
%% 第二种，先进行图像配准，根据像素级的偏移量调整主图的平地相位，以及辅图的平地相位，之后再去平地
imageA = load('SLC_Master.mat');
imageB = load('SLC_Slave.mat');

%卫星参数
file_ID     = '.\';
product_ID  = 'TH02-01B_InSAR_20191008_0000001927_001_004_008_L1B.Master.meta';
xmldoc = xml_read([file_ID,product_ID,'.xml']);
c = 3*10e8;
lambda = c/(xmldoc.carrierFrequency*10e9); %波长

%裁剪部分图像
img_rows = [10089:11089];
img_cols = [9075:10075];
Master = imageA.Master(img_rows,img_cols);
Slave = imageB.Slave(img_rows,img_cols);
coa_mast_ind = 1;
coa_slave_ind = 2;

%粗配准
[move_row_part,move_col_part] = deal(zeros(1,2),zeros(1,2));
[move_row_part(coa_slave_ind),move_col_part(coa_slave_ind),~] = real_coherent(Master,Slave,0);
Sar_cut_images = cell(2,1);
Sar_cut_images{coa_mast_ind} = Master;
Sar_cut_images{coa_slave_ind} = Slave;
[image_coarse_regis,first_pixel_row,first_pixel_col] = cut_multi_images_to_regis(move_row_part,move_col_part,Sar_cut_images,coa_mast_ind);
[fine_size_rows,fine_size_cols] = size(image_coarse_regis{coa_mast_ind});
%精配准
image_fine_regis = cell(1,2);
image_fine_regis{coa_mast_ind} = image_coarse_regis{coa_mast_ind};
image_fine_regis{coa_slave_ind} = registration_subpixel1(image_coarse_regis{coa_mast_ind},image_coarse_regis{coa_slave_ind});
image_fine = image_quantify(image_fine_regis{coa_slave_ind},60);
image_fine_gray = mat2gray(image_fine);
figure;   imagesc(image_fine_gray);    title('fineregis SAR image');  drawnow;
%干涉（未去平地）
diff_interf_phase = zeros(fine_size_rows,fine_size_cols);

%interf_ref_terrain_m = flat_int_phase(:,:,m);
diff_interf_phase(:,:) = angle(image_fine_regis{coa_mast_ind}.*...
    conj(image_fine_regis{coa_slave_ind})); %.*exp(complex(0, interf_ref_terrain_m)));
%figure;   imagesc(interf_ref_terrain_m); title('flat earth interferometric phase');  drawnow;
figure;   imagesc(diff_interf_phase(:,:));  title('interferogram phase');  drawnow;

%对平地斜距配准
master_rows = [10089:11089-move_row_part(2)];
master_cols = [9075-move_col_part(2):10075];
slave_rows = [10089+move_row_part(2):11089];
slave_cols = [9075:10075+move_col_part(2)];
flat_range_Master_cut = flat_range_master(master_rows,master_cols);
flat_range_Slave_cut = flat_range_slave(slave_rows,slave_cols);

%平地干涉相位 一发双收卫星系统， 其他系统公式可能会变
flat_inter_phase = 2*pi*(flat_range_Master_cut-flat_range_Slave_cut)/lambda;
figure;   imagesc(flat_inter_phase);    title('flat inter phase');  drawnow;

%干涉去平地
def_interf_phase = zeros(fine_size_rows,fine_size_cols);

def_interf_phase(:,:) = angle(image_fine_regis{coa_slave_ind}.*...
    conj(image_fine_regis{coa_mast_ind}).*exp(complex(0, flat_inter_phase)));
figure;   imagesc(def_interf_phase(:,:));  title('deflat interferogram phase');  drawnow;