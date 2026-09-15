
%本程序是地基SAR图像仿真程序和干涉处理程序，整理于2026..09
%设置系统参数生成主辅图像、斜距、控制点、真实相位等信息
%对主辅图像进行干涉处理：图像配准、去参考地形、相位滤波、相位解缠、控制点校正、相位精度评估、真实地形DEM反演
clear all
close all
clc
first=1;%%%%表示第一次运行程序

%% 第二时刻与接收噪声参数
DeltaT = 24*3600;       % 时间间隔，s；这里只记录观测间隔
DeltaH_max = -0.020;    % 最大高度变化，m；负值为下沉

AddReceiverNoise = true;
SNR_dB = 25;           % 参考有效回波区域的成像后信噪比
NoiseSeed = 2;         % 固定随机种子，便于重复实验

%% 自动保存目录
projectDir = fileparts(mfilename('fullpath'));
if isempty(projectDir)
    projectDir = pwd;
end

resultDir = fullfile(projectDir, 'output', ...
    ['interferometric_2_' datestr(now,'yyyymmdd_HHMMSS_FFF')]);

imageDir = fullfile(resultDir, 'images');
mkdir(imageDir);
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

%% 第二时刻：局部地形高度变化
height_t1 = height;

% 椭圆形变区：中心位置与半轴，单位 m
deform_center_x = 50;
deform_center_y = 15;
deform_radius_x = 70;
deform_radius_y = 35;

q = sqrt( ...
    ((X-deform_center_x)/deform_radius_x).^2 + ...
    ((Y-deform_center_y)/deform_radius_y).^2 );

deformShape = zeros(size(height));
inside = q < 1;

% 中心变化最大，向边缘平滑减小；区域外严格为零
deformShape(inside) = 0.5*(1 + cos(pi*q(inside)));

deltaHeight = DeltaH_max * deformShape;
height_t2 = height_t1 + deltaHeight;

clear q inside deformShape
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

%% 第二时刻：两个雷达的位置均保持不变
% 直接用 X、Y、height_t2 计算，避免再创建三维 DEM 副本

R_master_t2 = sqrt( ...
    (X-radarPos_master(1)).^2 + ...
    (Y-radarPos_master(2)).^2 + ...
    (height_t2-radarPos_master(3)).^2 );

U_master_t2 = (X-radarPos_master(1))./R_master_t2;

R_slave_t2 = sqrt( ...
    (X-radarPos_slave(1)).^2 + ...
    (Y-radarPos_slave(2)).^2 + ...
    (height_t2-radarPos_slave(3)).^2 );

U_slave_t2 = (X-radarPos_slave(1))./R_slave_t2;

% 这些三维中间变量后面不再使用，释放内存
clear V_master V_slave




%% 图像采样与参数相应
dr = range_pixel_spacing; %距离向采样间隔
du = azimuth_pixel_spacing_u; %方位向采样间隔

rho_r = range_resolution; %距离向分辨率;
rho_u = direction_resolution_u; %方位向分辨率;
Rref = mean(R_master(:)); %参考距离

%% 根据两个时刻、两个位置共同确定图像范围
rLower = min([min(R_master(:)), min(R_slave(:)), ...
              min(R_master_t2(:)), min(R_slave_t2(:))]);

rUpper = max([max(R_master(:)), max(R_slave(:)), ...
              max(R_master_t2(:)), max(R_slave_t2(:))]);

uLower = min([min(U_master(:)), min(U_slave(:)), ...
              min(U_master_t2(:)), min(U_slave_t2(:))]);

uUpper = max([max(U_master(:)), max(U_slave(:)), ...
              max(U_master_t2(:)), max(U_slave_t2(:))]);

rmin = Rref + floor((rLower-Rref)/dr)*dr - 4*rho_r;
rmax = Rref + ceil((rUpper-Rref)/dr)*dr + 4*rho_r;

umin = floor(uLower/du)*du - 4*rho_u;
umax = ceil(uUpper/du)*du + 4*rho_u;

rAxis = rmin:dr:rmax;
uAxis = umin:du:umax;

ImageAzimath =numel(uAxis); %图像行数：方向
ImageRange = numel(rAxis); %图像列数：距离

%% 生成两个时刻的主辅复数图像
% 1：第一时刻主图像
% 2：第一时刻辅图像
% 3：第二时刻主图像
% 4：第二时刻辅图像

tic;
SLC_clean = cell(1,4);

SLC_clean{1} = generateSLC_ground( ...
    R_master, U_master, alpha, lambda, ...
    rAxis, uAxis, rho_r, rho_u, Rref);

SLC_clean{2} = generateSLC_ground( ...
    R_slave, U_slave, alpha, lambda, ...
    rAxis, uAxis, rho_r, rho_u, Rref);

SLC_clean{3} = generateSLC_ground( ...
    R_master_t2, U_master_t2, alpha, lambda, ...
    rAxis, uAxis, rho_r, rho_u, Rref);

SLC_clean{4} = generateSLC_ground( ...
    R_slave_t2, U_slave_t2, alpha, lambda, ...
    rAxis, uAxis, rho_r, rho_u, Rref);

elapsed = toc;
fprintf('四幅图像生成耗时：%.2f 秒\n',elapsed);

%% 计算接收噪声功率
% 使用第一幅图像的有效回波区域定义参考信号功率
amp = abs(SLC_clean{1});
signalMask = amp > 0.05*max(amp(:));
P_signal = mean(amp(signalMask).^2);

if AddReceiverNoise
    P_noise = P_signal / 10^(SNR_dB/10);
else
    P_noise = 0;
end

%% 接收噪声的空间滤波核
% 按当前理想成像响应处理，使过采样图像中的噪声具有空间相关性
os_r = rho_r/dr;
os_u = rho_u/du;

tr = (-ceil(4*os_r):ceil(4*os_r))/os_r;
tu = (-ceil(4*os_u):ceil(4*os_u))/os_u;

kr = ones(size(tr));
ku = ones(size(tu));

idx = tr ~= 0;
kr(idx) = sin(pi*tr(idx))./(pi*tr(idx));

idx = tu ~= 0;
ku(idx) = sin(pi*tu(idx))./(pi*tu(idx));

kr = kr/norm(kr);
ku = ku(:)/norm(ku);

%% 分别加入独立接收噪声
rng(NoiseSeed);
SLC = cell(1,4);

noiseRows = ImageAzimath + numel(ku)-1;
noiseCols = ImageRange + numel(kr)-1;

for k = 1:4
    if P_noise > 0
        whiteNoise = ( ...
            randn(noiseRows,noiseCols) + ...
            1j*randn(noiseRows,noiseCols))/sqrt(2);

        receiverNoise = sqrt(P_noise) * ...
            conv2(ku,kr,whiteNoise,'valid');

        SLC{k} = SLC_clean{k} + receiverNoise;
    else
        SLC{k} = SLC_clean{k};
    end
end

clear whiteNoise receiverNoise amp

%% 保存生成的图像和成像参数
save(fullfile(resultDir,'SLC_data.mat'), ...
    'SLC','SLC_clean','rAxis','uAxis','Rref', ...
    'lambda','fc','B','L_aperture', ...
    'radarPos_master','radarPos_slave', ...
    'DeltaT','SNR_dB','P_noise', ...
    'NoiseSeed','AddReceiverNoise','-v7.3');

%% 沿用原程序的线性幅度显示
imageNames = { ...
    '第一时刻主图像', ...
    '第一时刻辅图像', ...
    '第二时刻主图像', ...
    '第二时刻辅图像'};

for k = 1:4
    figure('Name',imageNames{k});
    imagesc(abs(SLC{k}));
    colormap(gray(256));
    colorbar;
    xlabel('距离向像素编号（列）');
    ylabel('方向向像素编号（行）');
    title([imageNames{k},'：线性幅度']);
end

%% 干涉配对
% 同时刻、不同位置：空间干涉
% 同位置、不同时刻：时间干涉
pairIndex = [ ...
    1 2; ...
    3 4; ...
    1 3; ...
    2 4];

pairNames = { ...
    '第一时刻空间干涉', ...
    '第二时刻空间干涉', ...
    '主位置时间干涉', ...
    '辅位置时间干涉'};

result = struct([]);

for p = 1:size(pairIndex,1)

    image_master = SLC{pairIndex(p,1)};
    image_slave  = SLC{pairIndex(p,2)};

    %% 整像素粗配准
    [nr,nc] = size(image_master);

    [Main_coarse,Slave_coarse,move_r,move_c] = ...
        registration_pixel( ...
        image_master,image_slave,nr,nc);

    % 当前 registration_pixel 会生成互相关图，在这里补标题
    title(['粗配准互相关：',pairNames{p}]);
    xlabel('相关矩阵列号');
    ylabel('相关矩阵行号');

    fprintf('%s：粗配准行偏移 %g，列偏移 %g\n', ...
        pairNames{p},move_r,move_c);

    %% 分块亚像素配准
    [nrC,ncC] = size(Main_coarse);

    [~,Slave_registered,registrationValid,regInfo] = ...
        regis_subpixel( ...
        Main_coarse,Slave_coarse,nrC,ncC);

    fprintf('%s：最终采用 %d 个配准块\n', ...
        pairNames{p},regInfo.accepted);

    %% 更新裁剪后的坐标
    firstRow = 1 + max(0,-move_r);
    firstCol = 1 + max(0,-move_c);

    uAxis_reg = uAxis(firstRow:firstRow+nrC-1);
    rAxis_reg = rAxis(firstCol:firstCol+ncC-1);

    %% 从配准后的复数图像计算干涉相位
    wrapped_phase = interferometric_phase( ...
        Main_coarse,Slave_registered);

    ampThreshold = 0.05*max(abs(Main_coarse(:)));

    valid = registrationValid ...
        & isfinite(wrapped_phase) ...
        & abs(Main_coarse)>ampThreshold ...
        & abs(Slave_registered)>ampThreshold;

    wrapped_phase(~valid) = NaN;

    %% 记录这一组的实际处理结果
    result(p).name = pairNames{p};
    result(p).pairIndex = pairIndex(p,:);
    result(p).phase = wrapped_phase;
    result(p).valid = valid;
    result(p).move_r = move_r;
    result(p).move_c = move_c;
    result(p).regInfo = regInfo;
    result(p).rAxis = rAxis_reg;
    result(p).uAxis = uAxis_reg;
    result(p).master = Main_coarse;
    result(p).slaveRegistered = Slave_registered;

    %% 显示配准偏移
    figure('Name',[pairNames{p},'配准偏移']);
    imagesc(regInfo.dCol);
    colorbar;
    xlabel('距离向像素编号（列）');
    ylabel('方向向像素编号（行）');
    title([pairNames{p},'：精配准列偏移 / 像素']);

    %% 显示干涉相位
    figure('Name',pairNames{p});
    imagesc(wrapped_phase);
    colormap(parula(256));
    caxis([-pi,pi]);
    colorbar;
    xlabel('距离向像素编号（列）');
    ylabel('方向向像素编号（行）');
    title([pairNames{p},' / rad']);
end

%% 保存配准和干涉结果
save(fullfile(resultDir,'interferometric_results.mat'), ...
    'result','DeltaT','lambda','-v7.3');

%% 保存本次使用的场景输入，方便以后重复实验
save(fullfile(resultDir,'scene_input.mat'), ...
    'x','y','height_t1','height_t2', ...
    'DeltaH_max','DeltaT', ...
    'deform_center_x','deform_center_y', ...
    'deform_radius_x','deform_radius_y','-v7.3');

%% 保存本次运行的全部图窗
% 脚本开头已有 close all
figs = findall(groot,'Type','figure');
[~,order] = sort([figs.Number]);
figs = figs(order);

for k = 1:numel(figs)
    fileStem = sprintf('Figure_%02d',k);

    exportgraphics(figs(k), ...
        fullfile(imageDir,[fileStem,'.png']), ...
        'Resolution',160);

    savefig(figs(k), ...
        fullfile(imageDir,[fileStem,'.fig']));
end

fprintf('结果已保存到：\n%s\n',resultDir);


%% 两个时刻的空间干涉相位差分
% 第一时刻空间相位 - 第二时刻空间相位
% 复用主图像之间已有的配准关系

sp1 = result(1);  % 第一时刻空间干涉
sp2 = result(2);  % 第二时刻空间干涉
tm  = result(3);  % 两时刻主图像的配准结果

%% 1. 确定各结果在原始参考图像中的裁剪起点
r1 = 1 + max(0, -sp1.move_r);
c1 = 1 + max(0, -sp1.move_c);

r2 = 1 + max(0, -sp2.move_r);
c2 = 1 + max(0, -sp2.move_c);

rt = 1 + max(0, -tm.move_r);
ct = 1 + max(0, -tm.move_c);

%% 2. 取第一时刻空间图与时间配准结果的共同范围
% 两者都以第一时刻主图像为参考
rowBegin = max(r1, rt);
rowEnd = min( ...
    r1 + size(sp1.phase,1) - 1, ...
    rt + size(tm.phase,1) - 1);

colBegin = max(c1, ct);
colEnd = min( ...
    c1 + size(sp1.phase,2) - 1, ...
    ct + size(tm.phase,2) - 1);

assert(rowBegin <= rowEnd && colBegin <= colEnd, ...
    '空间干涉图与主图像配准结果没有共同范围。');

rows = rowBegin:rowEnd;
cols = colBegin:colEnd;

% 共同范围在各自结果矩阵中的索引
rows1 = rows - r1 + 1;
cols1 = cols - c1 + 1;

rowsT = rows - rt + 1;
colsT = cols - ct + 1;

% 第一时刻空间相位直接裁剪，不需要插值
phase1 = sp1.phase(rows1, cols1);

%% 3. 利用已有配准关系，找到第二时刻对应位置
[colGrid, rowGrid] = meshgrid(cols, rows);

% 第一时刻主图像原始像素 -> 第二时刻主图像原始像素
rowInMaster2 = rowGrid + tm.move_r ...
    + tm.regInfo.dRow(rowsT, colsT);

colInMaster2 = colGrid + tm.move_c ...
    + tm.regInfo.dCol(rowsT, colsT);

% 转换成第二时刻空间干涉图内部的像素编号
rowQuery = rowInMaster2 - r2 + 1;
colQuery = colInMaster2 - c2 + 1;

%% 4. 将第二时刻空间相位对齐到第一时刻
% 插值单位复数表示，避免直接插值相位时跨越正负 pi
phase2Unit = exp(1j * sp2.phase);

phase2UnitReg = interp2( ...
    phase2Unit, colQuery, rowQuery, 'linear', NaN);

phase2Reg = angle(phase2UnitReg);

% 插值后复数幅度接近零时，相位没有可靠定义
phase2Reg(abs(phase2UnitReg) < 1e-12) = NaN;

%% 5. 相位差分，并缠绕回 [-pi, pi]
phaseSpaceDiff = angle(exp(1j * (phase1 - phase2Reg)));

% 沿用主图像配准已有的有效范围
% 原相位中的 NaN 和超出插值范围的 NaN 会自动保留
phaseSpaceDiff(~tm.valid(rowsT, colsT)) = NaN;

%% 6. 保存差分结果及对应坐标
spaceDiff = struct;
spaceDiff.name = '第一时刻空间相位减第二时刻空间相位';
spaceDiff.phase = phaseSpaceDiff;
spaceDiff.phase_t1 = phase1;
spaceDiff.phase_t2_registered = phase2Reg;
spaceDiff.rowIndex = rows;
spaceDiff.colIndex = cols;
spaceDiff.rAxis = sp1.rAxis(cols1);
spaceDiff.uAxis = sp1.uAxis(rows1);

%% 7. 显示差分相位
figure('Name', '两时刻空间干涉相位差分');

imagesc(cols, rows, phaseSpaceDiff);
colormap(parula(256));
caxis([-pi, pi]);
colorbar;

xlabel('第一时刻主图像距离向像素编号（列）');
ylabel('第一时刻主图像方向向像素编号（行）');
title('空间干涉相位差分：第一时刻减第二时刻 / rad');