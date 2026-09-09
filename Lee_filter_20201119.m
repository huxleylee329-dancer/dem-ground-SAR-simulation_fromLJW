
function [phase_filter] = Lee_filter_20201119(Original_Interferogram,Win_size,coherence)
%% Lee滤波：随着条纹方向改变窗口大小的自适应滤波，因此相当于复杂度因子滤波里面的改进3，滤波参数会发生改变
%窗参数
 
Radius = (Win_size-1)/2; %估计窗i方向半径
[Size_i,Size_j]=size(Original_Interferogram);
% coherence = 0.5*ones(Size_i,Size_j);
%扩充估计相位矩阵
phase_update = zeros(Size_i+Win_size-1,Size_j+Win_size-1);    
coherence_update= zeros(Size_i+Win_size-1,Size_j+Win_size-1); 

[Ni,Nj] = size(phase_update);
phase_update(Radius+1:Ni-Radius,Radius+1:Nj-Radius) = Original_Interferogram;
coherence_update(Radius+1:Ni-Radius,Radius+1:Nj-Radius) = coherence;

%flipud:矩阵上下交换，fliplr:矩阵左右交换  %%对半径行和列进行边缘赋值  %%目的是弱化多视对边缘的影响 
%对边缘的7行7列加值
phase_update(1:Radius,1:Nj) = flipud(phase_update(Radius+1:2*Radius,1:Nj)); 
phase_update(1:Ni,1:Radius) = fliplr(phase_update(1:Ni,Radius+1:2*Radius));  
phase_update(Ni-Radius+1:Ni,1:Nj) = flipud(phase_update(Ni-2*Radius+1:Ni-Radius,1:Nj));
phase_update(1:Ni,Nj-Radius+1:Nj) = fliplr(phase_update(1:Ni,Nj-2*Radius+1:Nj-Radius));

coherence_update(1:Radius,1:Nj) = flipud(coherence_update(Radius+1:2*Radius,1:Nj)); 
coherence_update(1:Ni,1:Radius) = fliplr(coherence_update(1:Ni,Radius+1:2*Radius));  
coherence_update(Ni-Radius+1:Ni,1:Nj) = flipud(coherence_update(Ni-2*Radius+1:Ni-Radius,1:Nj));
coherence_update(1:Ni,Nj-Radius+1:Nj) = fliplr(coherence_update(1:Ni,Nj-2*Radius+1:Nj-Radius));

CSD=Standard_Deviation(coherence_update);%相干系数标准差
PSD=Phase_Standard_Deviation(phase_update);%相位标准差
%%
    wb = waitbar(0,'phase filtering ...');
%将相位转化为复数
j=sqrt(-1);
phase_update = exp(j*phase_update);     
%局部频率估计参数
phase_linear = zeros(Ni,Nj); % 线形相位矩阵
phase_filter = zeros(Ni,Nj); % 滤波后相位矩阵
nn = 32;                     % 2D FFT点数
mm = 32;                     % CZT点数
mn = mm*nn;                  % 频谱总分辨率（mn*mn）
phase_estimation = zeros(Win_size,Win_size); %估计窗内相位矩阵
tempi = (0:Win_size-1)';  % 行线形相位网格
tempj = 0:Win_size-1;     % 列线形相位网格
w0 = exp(-j*2*pi/mn);     % CZT参数
mins=Radius-1;
%% 观察条纹频率值范围
fii_1 = zeros(Ni,Nj);
fjj_1 = zeros(Ni,Nj);
residual_noise_phase = zeros(Ni,Nj);
fii_2 = zeros(Ni,Nj);
fjj_2 = zeros(Ni,Nj);
for ii = Radius+1:Ni-Radius
  for jj = Radius+1:Nj-Radius
      
      phase_estimation = phase_update(ii-Radius:ii+Radius,jj-Radius:jj+Radius);
      Conherence_Window = coherence_update(ii-Radius+mins:ii+Radius-mins,jj-Radius+mins:jj+Radius-mins);
      CSD_Window = CSD(ii-Radius+mins:ii+Radius-mins,jj-Radius+mins:jj+Radius-mins);
      Num=sum(sum((CSD_Window<0.4)));
      Conherence_Window=Conherence_Window.*(CSD_Window<0.4);
      Gamma=sum(sum(Conherence_Window))./Num;
      PSD_Window = PSD(ii-Radius+mins:ii+Radius-mins,jj-Radius+mins:jj+Radius-mins);%估计窗口内的PSD值

      %% 频率估计前预滤波
      %此处为与Improved_Goldstein_filter_20180402不同的部分：没有采用变窗口，窗口大小是固定的
      window_mean = Mean_Filter_Correlation_1(phase_estimation,0.5);
      
      %此处为仅改进1,：预滤波里面采用变窗口
%       [window_mean,window(ii,jj)] = Mean_Filter_Correlation_20201119(phase_estimation,Gamma+0.4,PSD_Window,Win_size);

      %% 此处为仅改进2：条纹频率提取和传统Goldstein滤波器结合
      phase_2fft = fftshift(fft2(window_mean,nn,nn)); %窗口内2DFFT
%       phase_2fft = fftshift(fft2(phase_estimation,nn,nn)); %Lee滤波直接对预滤波之前的phase_estimation窗口内2DFFT     
      [val_j,ind_j] = max(max(abs(phase_2fft)));      %ind_j最大值列号 
      [val_i,ind_i] = max(abs(phase_2fft(:,ind_j)));  %ind_i最大值行号 
      fi = double(ind_i-nn/2-2)/nn;
      fj = double(ind_j-nn/2-2)/nn;
      ei = exp(j*2*pi*fi);
      ej = exp(j*2*pi*fj);
      % CZT频谱主瓣精细估计
      phasef_czt1 = czt(window_mean,3*mm,w0,ei); %i(行)方向CZT,czt按列计算
%       phasef_czt1 = czt(phase_estimation,3*mm,w0,ei); %Lee滤波直接对预滤波之前的phase_estimation i(行)方向CZT,czt按列计算     
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
      phase_linear(ii-Radius:ii+Radius,jj-Radius:jj+Radius) = ei.^tempi*(ej.^tempj);
      % 初相,sum(sum()) or mean(mean())
      phase_0_matrix = phase_update(ii-Radius:ii+Radius,jj-Radius:jj+Radius).*...
                       conj(phase_linear(ii-Radius:ii+Radius,jj-Radius:jj+Radius));
      
      %仅改进2：条纹频率提取和传统Goldstein滤波器结合，所以传统Goldstein滤波器是对提取条纹频率之后的相位滤波
      %%重新估计条纹频率
      %FFT频率估计
      phase_2fft = fftshift(fft2(phase_0_matrix,nn,nn)); %窗口内2DFFT
      [val_j,ind_j] = max(max(abs(phase_2fft)));      %ind_j最大值列号 
      [val_i,ind_i] = max(abs(phase_2fft(:,ind_j)));  %ind_i最大值行号 
      fi = double(ind_i-nn/2-2)/nn;
      fj = double(ind_j-nn/2-2)/nn;
      ei = exp(j*2*pi*fi);
      ej = exp(j*2*pi*fj);
      % CZT频谱主瓣精细估计
      phasef_czt1 = czt(phase_0_matrix,3*mm,w0,ei); %i(行)方向CZT,czt按列计算
      phasef_czt2 = czt(phasef_czt1.',3*mm,w0,ej);    %j(列)方向CZT
      phasef_czt2 = phasef_czt2.';
      [val_j,ind_j] = max(max(abs(phasef_czt2)));
      [val_i,ind_i] = max(abs(phasef_czt2(:,ind_j)));
      movei = double(ind_i-1)/mn;  %i(行)方向细化后偏移量
      movej = double(ind_j-1)/mn;  %j(列)方向细化后偏移量
      fii_2(ii,jj) = fi+movei;
      fjj_2(ii,jj) = fj+movej;
% 
%       phase_0 = Improved_Goldstein_Correlation_20201119(phase_0_matrix,0.5); 
% 
%       %初相 + 线性相位 
%       phase_filter(ii,jj) = phase_0 * phase_linear(ii,jj);
         
     %% 对预滤波之后的窗口进行传统Goldstein滤波
      % 记得更改Traditional_Goldstein_Correlation_20201119里的滤波参数r=0.5，此时r不随着相干系数Gamma+0.4变化，只有在改进3时r才会变化
%       phase_0 = Traditional_Goldstein_Correlation_20201119(window_mean,Gamma+0.4); 
%       phase_filter(ii,jj) = phase_0;%传统Goldstein滤波后得到的相位
      
     %% 仅改进3：传统Goldstein滤波器参数随着噪声频率和相干系数自适应改变
      Conherence_Window2 = coherence_update(ii-Radius+mins:ii+Radius-mins,jj-Radius+mins:jj+Radius-mins);
      Gamma2=sum(sum(Conherence_Window2))./3/3;%Conherence_Window2窗口大小为3*3;Gamma2为平均相干系数，没有筛选标准差小于0.4的点
      D = (1-Gamma2)*(1+sqrt((fii_2(ii,jj)^2+fjj_2(ii,jj)^2))/2);%复杂度因子,即Goldstein滤波的滤波参数alpha
      phase_0 = Improved_Goldstein_Correlation_20201119(window_mean,D);
      phase_filter(ii,jj) = phase_0;

  end
     str = ['滤波进行中',num2str(fix((100*(ii-Radius)/(Ni-2*Radius)))),'%'];
 waitbar(ii/(Ni-2*Radius),wb,str);
end
%%
% %取复相角(rad),取有效范围
phase_filter = angle(phase_filter);
phase_filter = phase_filter(Radius+1:Ni-Radius,Radius+1:Nj-Radius);
end







