function phase_filter = Slope_Adapt_Phase_Filter(phase)
% 基于局部频率估计的斜坡自适应滤波
% FFT+CZT 提升局部频率估计的精度
% i-direction: 行方向
% j-direction: 列方向
% Inputs: 原始干涉相位矩阵、行数、列数
% Outputs: 滤波后干涉相位矩阵
j=sqrt(-1);
% %----------------- display information -------------------%
% display([repmat('-', 1, 8),' phase filter ',repmat('-', 1, 8)])
% display(['method:             ','SAPF'])
% display(repmat('-', 1, 30))
%% 窗参数
% 当干涉相位方差变化范围不大时，估计窗size可取9，11，13，基线400m时，取9较优
% winsize_e_i = 21;           %估计窗i方向尺寸
% winsize_e_j = 21;           %估计窗j方向尺寸
% winsize_f_i = 15;  %滤波窗i方向尺寸
% winsize_f_j = 15;  %滤波窗j方向尺寸
winsize_e_i = 9;           %估计窗i方向尺寸
winsize_e_j = 9;           %估计窗j方向尺寸
winsize_f_i = 7;  %滤波窗i方向尺寸
winsize_f_j = 7;  %滤波窗j方向尺寸0104修改
wine_i = (winsize_e_i-1)/2; %估计窗i方向半径
wine_j = (winsize_e_j-1)/2; %估计窗j方向半径
winf_i = (winsize_f_i-1)/2; %滤波窗i方向半径
winf_j = (winsize_f_j-1)/2; %滤波窗j方向半径
[size_i,size_j]=size(phase);

%% 扩充相位矩阵
% 扩充估计相位矩阵
phase_update = zeros(size_i+winsize_e_i-1,size_j+winsize_e_j-1);             
[Ni,Nj] = size(phase_update);
phase_update(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j) = phase;
% flipud:矩阵上下交换，fliplr:矩阵左右交换
phase_update(1:wine_i,1:Nj) = flipud(phase_update(wine_i+1:2*wine_i,1:Nj)); 
phase_update(1:Ni,1:wine_j) = fliplr(phase_update(1:Ni,wine_j+1:2*wine_j));  
phase_update(Ni-wine_i+1:Ni,1:Nj) = flipud(phase_update(Ni-2*wine_i+1:Ni-wine_i,1:Nj));
phase_update(1:Ni,Nj-wine_j+1:Nj) = fliplr(phase_update(1:Ni,Nj-2*wine_j+1:Nj-wine_j));
% 将相位转化为复数
phase_update = exp(j*phase_update);     
%% 局部频率估计参数
phase_linear = zeros(Ni,Nj); % 线形相位矩阵
phase_filter = zeros(Ni,Nj); % 滤波后相位矩阵
nn = 32;                     % 2D FFT点数
mm = 32;                     % CZT点数
mn = mm*nn;                  % 频谱总分辨率（mn*mn）
phase_estimation = zeros(winsize_e_i,winsize_e_j); %估计窗内相位矩阵
tempi = (0:winsize_f_i-1)';  % 行线形相位网格
tempj = 0:winsize_f_j-1;     % 列线形相位网格
w0 = exp(-j*2*pi/mn);        % CZT参数
%% 自适应滤波
%% 自适应滤波
%wb = waitbar(0,'phase filtering ...');
for ii = wine_i+1:Ni-wine_i
  for jj = wine_j+1:Nj-wine_j
      phase_estimation(1:winsize_e_i,1:winsize_e_j) = phase_update(ii-wine_i:ii+wine_i,jj-wine_j:jj+wine_j);
      % FFT频率估计
      phase_2fft = fftshift(fft2(phase_estimation,nn,nn)); %窗口内2DFFT
      [val_j,ind_j] = max(max(abs(phase_2fft)));      %ind_j最大值列号 
      [val_i,ind_i] = max(abs(phase_2fft(:,ind_j)));  %ind_i最大值行号 
      fi = double(ind_i-nn/2-2)/nn;
      fj = double(ind_j-nn/2-2)/nn;
      ei = exp(j*2*pi*fi);
      ej = exp(j*2*pi*fj);
      % CZT频谱主瓣精细估计
      phasef_czt1 = czt(phase_estimation,3*mm,w0,ei); %i(行)方向CZT,czt按列计算
      phasef_czt2 = czt(phasef_czt1.',3*mm,w0,ej);    %j(列)方向CZT
      phasef_czt2 = phasef_czt2.';
      [val_j,ind_j] = max(max(abs(phasef_czt2)));
      [val_i,ind_i] = max(abs(phasef_czt2(:,ind_j)));
      movei = double(ind_i-1)/mn;  %i(行)方向细化后偏移量
      movej = double(ind_j-1)/mn;  %j(列)方向细化后偏移量
      fii = fi+movei;
      fjj = fj+movej;       %CZT提升分辨率后，频谱峰值位置
      ei = exp(j*2*pi*fii); %i(行)方向局部频率（指数）
      ej = exp(j*2*pi*fjj); %j(列)方向局部频率（指数）
      % 线性相位
      phase_linear(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j) = ei.^tempi*(ej.^tempj);
      % 初相,sum(sum()) or mean(mean())
      phase_0_matrix = phase_update(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j).*...
                       conj(phase_linear(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j));
      phase_0 =  mean(mean(phase_0_matrix)); %均值滤波  
      % 初相 + 线性相位 
      phase_filter(ii,jj) = phase_0 * phase_linear(ii,jj);
  end
  ii
 %str = ['滤波进行中',num2str(fix((100*ii/(Ni-2*wine_i)))),'%'];
 %waitbar(ii/(Ni-2*wine_i),wb,str);
 %pause(0.05);
end
%close(wb);
%% 取复相角(rad),取有效范围
phase_filter = angle(phase_filter);
phase_filter = phase_filter(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j);
%% 输出
% figure;
% imagesc(phase_filter)
% % axis xy
% colormap(jet)
% title('interferometric phase (SAPF)','fontWeight','Bold')
clear phase_update phase_linear phase_estimation;
end