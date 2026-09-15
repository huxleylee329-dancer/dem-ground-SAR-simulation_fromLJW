
%本程序是地基SAR图像仿真程序和干涉处理程序，整理于2026..09
%设置系统参数生成主辅图像、斜距、控制点、真实相位等信息
%对主辅图像进行干涉处理：图像配准、去参考地形、相位滤波、相位解缠、控制点校正、相位精度评估、真实地形DEM反演
clear all
close all
clc
first=1;%%%%表示第一次运行程序
%% %%%地基雷达参数%%%%%
C0=299792458; %光速
f_low = 16.2e9; %雷达工作低频
f_high = 17.2e9; %雷达工作高频
fc = (f_low + f_high)/2; %雷达工作中心频率
Wlength = C0/fc; %雷达工作波长
lambda = Wlength; 
L_aperture = 1.4; %等效阵列孔径长度
B_length = 1.4; %等效阵列孔径宽度
antenna_gain = 1; %全向阵元的简化假设

C = 4*pi;          %%双程距离相位系数,暂时按照等效单站双程模型
B = f_high - f_low; %雷达工作带宽
range_resolution = C0/(2*B); %雷达工作距离分辨率,现在算来大概是0.15米左右??

%% 理想等效孔径模型
direction_resolution_u = lambda/(2*L_aperture); %理想等效孔径方位向分辨率
%输出图像的采样间隔
range_pixel_spacing = range_resolution/2; %距离向采样间隔
azimuth_pixel_spacing_u = direction_resolution_u/8; %方位向采样间隔

%% 主辅观测位置
D = 3000;                         % 水平观测距离，m

radarPos_master = [0, -D, 20];

B_vector = [0, 0, B_length];       % 沿局部 Z 方向的基线
radarPos_slave = radarPos_master + B_vector;

%% %%%%%%%%布设地面场景点参数%%%%%%
%%小坡面散射点网络
DEM_dx=0.5;
DEM_dy=0.01;

x=(-160+DEM_dx/2):DEM_dx:(160-DEM_dx/2);
y=(-75+DEM_dy/2):DEM_dy:(75-DEM_dy/2);

[X,Y]=meshgrid(x,y);
%向后升高，横向带有平缓隆起
height = 100 + 0.3*(Y+5) + 4*exp(-(X/18).^2);%地形高度
DEM = cat(3,X,Y,height);
%% 每个面片的面积
[dhdx, dhdy] = gradient(height, DEM_dx, DEM_dy);

facetArea = DEM_dx*DEM_dy ...
          .* sqrt(1 + dhdx.^2 + dhdy.^2);

%% 第一版：单位面积平均散射功率相同
sigma0 = ones(size(height));

%% 每个面片固定一个随机散射相位
rng(1);
RandAngle = unifrnd(0, 2*pi, size(height));

alpha = sqrt(sigma0.*facetArea) ...
      .* exp(1j*RandAngle);

%% 同一山体，主观测几何
V_master = DEM - reshape(radarPos_master,1,1,3);

R_master = sqrt(sum(V_master.^2,3));
U_master = V_master(:,:,1)./R_master;

%% 同一山体，辅观测几何
V_slave = DEM - reshape(radarPos_slave,1,1,3);

R_slave = sqrt(sum(V_slave.^2,3));
U_slave = V_slave(:,:,1)./R_slave;

%% 地形点层面的理论干涉相位
% 对应 MainSLC .* conj(SlaveSLC) 的相位约定
phase_truth_points = 4*pi/lambda .* (R_slave-R_master);
phase_truth_wrapped = angle(exp(1j*phase_truth_points));

%% 图像采样与参数相应
dr = range_pixel_spacing; %距离向采样间隔
du = azimuth_pixel_spacing_u; %方位向采样间隔

rho_r = range_resolution; %距离向分辨率;
rho_u = direction_resolution_u; %方位向分辨率;
Rref = mean(R_master(:)); %参考距离

Rall = [R_master(:); R_slave(:)];
Uall = [U_master(:); U_slave(:)];

rmin = Rref + floor((min(Rall)-Rref)/dr)*dr - 4*rho_r;
rmax = Rref + ceil((max(Rall)-Rref)/dr)*dr + 4*rho_r;

umin = floor(min(Uall)/du)*du - 4*rho_u;
umax = ceil(max(Uall)/du)*du + 4*rho_u;

rAxis = rmin:dr:rmax;
uAxis = umin:du:umax;

ImageAzimath =numel(uAxis); %图像行数：方向
ImageRange = numel(rAxis); %图像列数：距离

%% 生成主辅复数图像
tic;

MainSLC = generateSLC_ground( ...
    R_master, U_master, alpha, lambda, ...
    rAxis, uAxis, rho_r, rho_u, Rref);

SlaveSLC = generateSLC_ground( ...
    R_slave, U_slave, alpha, lambda, ...
    rAxis, uAxis, rho_r, rho_u, Rref);

elapsed = toc;

fprintf('主辅图像生成耗时：%.2f 秒\n',elapsed);

assert(isequal(size(MainSLC),[ImageAzimath,ImageRange]));
assert(isequal(size(SlaveSLC),size(MainSLC)));

assert(all(isfinite(MainSLC(:))));
assert(all(isfinite(SlaveSLC(:))));

%% 查看主辅幅度图，采用相同显示尺度

% 用主辅图像共同的最大幅度作为显示参考
ampRef = max([abs(MainSLC(:)); abs(SlaveSLC(:))]);

% 复数图像取幅度，再转换成 dB
Main_dB = 20*log10(abs(MainSLC)/max(ampRef,eps) + eps);
Slave_dB = 20*log10(abs(SlaveSLC)/max(ampRef,eps) + eps);

figure;

%% 查看主辅幅度图：沿用原程序的线性幅度显示
figure;
imagesc(abs(MainSLC));
colormap(gray(256));
colorbar;
xlabel('距离向像素编号（列）');
ylabel('方向向像素编号（行）');
title('主图像幅度');

figure;
imagesc(abs(SlaveSLC));
colormap(gray(256));
colorbar;
xlabel('距离向像素编号（列）');
ylabel('方向向像素编号（行）');
title('辅图像幅度');

%% 整像素粗配准：根据图像估计偏移
[nr,nc] = size(MainSLC);

[Main_coarse,Slave_coarse,move_r,move_c] = ...
    registration_pixel(MainSLC,SlaveSLC,nr,nc);

fprintf('粗配准偏移：行 %g，列 %g 像素\n', ...
        move_r,move_c);

%% 分块亚像素配准
[nrC,ncC] = size(Main_coarse);

tic;

[~,Slave_registered,registrationValid,regInfo] = ...
    regis_subpixel(Main_coarse,Slave_coarse,nrC,ncC);

fprintf('分块精配准耗时：%.2f 秒\n',toc);

fprintf('子块总数：%d，初筛通过：%d，最终采用：%d\n', ...
        regInfo.blocks,regInfo.candidates,regInfo.accepted);

%% 粗配准裁剪后，更新主图像坐标轴
firstRow = 1 + max(0,-move_r);
firstCol = 1 + max(0,-move_c);

uAxis_reg = uAxis(firstRow:firstRow+nrC-1);
rAxis_reg = rAxis(firstCol:firstCol+ncC-1);

%% 计算干涉相位
wrapped_phase = interferometric_phase( ...
    Main_coarse,Slave_registered);

ampThreshold = 0.05*max(abs(Main_coarse(:)));

valid = registrationValid ...
      & abs(Main_coarse)>ampThreshold ...
      & abs(Slave_registered)>ampThreshold;

wrapped_phase(~valid) = NaN;

%% 查看分块估计并拟合的距离方向偏移
figure;
imagesc(rAxis_reg-Rref,uAxis_reg,regInfo.dCol);
axis xy;
xlabel('相对斜距 / m');
ylabel('方向余弦 u');
title('精配准列偏移场 / 像素');
colorbar;

%% 查看干涉相位
figure;
imagesc(rAxis_reg-Rref,uAxis_reg,wrapped_phase);
axis xy;
caxis([-pi,pi]);
xlabel('相对斜距 / m');
ylabel('方向余弦 u');
title('分块配准后的干涉相位 / rad');
colorbar;