function [localFrequency_range,phase_filter] = calculate_localFreRange(phase_initial, window)
% input: initial inteferometric phase
% output: local frequency in the direction range, phase filtered
% 20200427 siyuanwang
%% 窗参数
winsize_e_i = window;           %估计窗i方向尺寸
winsize_e_j = window;           %估计窗j方向尺寸
winsize_f_i = window;  %滤波窗i方向尺寸
winsize_f_j = window;  %滤波窗j方向尺寸
wine_i = (winsize_e_i-1)/2; %估计窗i方向半径
wine_j = (winsize_e_j-1)/2; %估计窗j方向半径
winf_i = (winsize_f_i-1)/2; %滤波窗i方向半径
winf_j = (winsize_f_j-1)/2; %滤波窗j方向半径
[size_i,size_j]=size(phase_initial);
%% 扩充相位矩阵
% 扩充估计相位矩阵
phase_update = zeros(size_i+winsize_e_i-1,size_j+winsize_e_j-1);             
[Ni,Nj] = size(phase_update);
phase_update(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j) = phase_initial;
% flipud:矩阵上下交换，fliplr:矩阵左右交换
phase_update(1:wine_i,1:Nj) = flipud(phase_update(wine_i+1:2*wine_i,1:Nj)); 
phase_update(1:Ni,1:wine_j) = fliplr(phase_update(1:Ni,wine_j+1:2*wine_j));  
phase_update(Ni-wine_i+1:Ni,1:Nj) = flipud(phase_update(Ni-2*wine_i+1:Ni-wine_i,1:Nj));
phase_update(1:Ni,Nj-wine_j+1:Nj) = fliplr(phase_update(1:Ni,Nj-2*wine_j+1:Nj-wine_j));
% 将相位转化为复数
phase_update = exp(1i*phase_update);
%% 局部频率估计参数
nn = 32;                     % 2D FFT点数
mm = 32;                     % CZT点数
mn = mm*nn;                  % 频谱总分辨率（mn*mn）
phase_estimation = zeros(window,window); %估计窗内相位矩阵
w0 = exp(-1i*2*pi/mn);        % CZT参数
tempi = (0:winsize_f_i-1)';  % 行线形相位网格
tempj = 0:winsize_f_j-1;     % 列线形相位网格
localFrequency_range = zeros(Ni,Nj); % 局部频率估计矩阵
phase_linear = zeros(Ni,Nj); % 线形相位矩阵
phase_filter = zeros(Ni,Nj); % 滤波后相位矩阵
%% 局部频率估计
% wb = waitbar(0,'localFre calculating ...');
for ii = wine_i+1:Ni-wine_i
  for jj = wine_j+1:Nj-wine_j
      phase_estimation(1:winsize_e_i,1:winsize_e_j) = phase_update(ii-wine_i:ii+wine_i,jj-wine_j:jj+wine_j);
      % FFT频率估计
      phase_2fft = fftshift(fft2(phase_estimation,nn,nn)); %窗口内2DFFT
      [~,ind_j] = max(max(abs(phase_2fft)));      %ind_j最大值列号 
      [~,ind_i] = max(abs(phase_2fft(:,ind_j)));  %ind_i最大值行号 
      fi = double(ind_i-nn/2-2)/nn;
      fj = double(ind_j-nn/2-2)/nn;
      ei = exp(1i*2*pi*fi);
      ej = exp(1i*2*pi*fj);
      % CZT频谱主瓣精细估计
      phasef_czt1 = czt(phase_estimation,3*mm,w0,ei); %i(行)方向CZT,czt按列计算
      phasef_czt2 = czt(phasef_czt1.',3*mm,w0,ej);    %j(列)方向CZT
      phasef_czt2 = phasef_czt2.';
      [~,ind_j] = max(max(abs(phasef_czt2)));
      [~,ind_i] = max(abs(phasef_czt2(:,ind_j)));
      movei = double(ind_i-1)/mn;  %i(行)方向细化后偏移量
      movej = double(ind_j-1)/mn;  %j(列)方向细化后偏移量
      fii = fi+movei;
      fjj = fj+movej;       %CZT提升分辨率后，频谱峰值位置

      localFrequency_range(ii,jj) = fjj;
      
      ei = exp(1i*2*pi*fii); %i(行)方向局部频率（指数）
      ej = exp(1i*2*pi*fjj); %j(列)方向局部频率（指数）
      % 线性相位
      phase_linear(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j) = ei.^tempi*(ej.^tempj);
      % 初相,sum(sum()) or mean(mean())
      phase_0_matrix = phase_update(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j).*...
                       conj(phase_linear(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j));
      phase_0 =  mean(mean(phase_0_matrix)); %均值滤波  
      % 初相 + 线性相位 
      phase_filter(ii,jj) = phase_0 * phase_linear(ii,jj);

  end
% str = ['进行中',num2str(fix((100*ii/(Ni-2*wine_i)))),'%'];
% waitbar(ii/(Ni-2*wine_i),wb,str);
end
% close(wb);
%% 取复相角(rad),取有效范围
localFrequency_range = localFrequency_range(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j);
phase_filter = angle(phase_filter);
phase_filter = phase_filter(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j);
