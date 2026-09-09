
function [phase_filter] = Improved_Goldstein_filter_20180402(Original_Interferogram,Win_size,coherence)
%% 空间自适应滤波
%窗参数
 
Radius = (Win_size-1)/2; %估计窗i方向半径
[Size_i,Size_j]=size(Original_Interferogram);
%Q:为什么相干系数是自己设置呢？不应该是根据主辅图像求出来的吗？
%A:在仿真数据第一组里面应该把已知的相干系数传进来
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
% figure;
% imagesc(PSD);
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
mins=Radius-1;%用于计算CSD、PSD时窗口的范围

%%观察条纹频率值范围
fii_1 = zeros(Ni,Nj);
fjj_1 = zeros(Ni,Nj);
residual_noise_phase = zeros(Ni,Nj);
fii_2 = zeros(Ni,Nj);
fjj_2 = zeros(Ni,Nj);
for ii = Radius+1:Ni-Radius
  for jj = Radius+1:Nj-Radius
      phase_estimation = phase_update(ii-Radius:ii+Radius,jj-Radius:jj+Radius);
      %求出相干系数的标准差CSD是为了把3*3的估计窗口内标准差大于0.2的相干系数筛除，利用标准差小于0.2的相干系数点求出平均相干系数
      Conherence_Window = coherence_update(ii-Radius+mins:ii+Radius-mins,jj-Radius+mins:jj+Radius-mins);%窗口大小为3*3
      CSD_Window = CSD(ii-Radius+mins:ii+Radius-mins,jj-Radius+mins:jj+Radius-mins);
      Num=sum(sum((CSD_Window<0.2)));%统计3*3的估计窗口内相干系数标准差小于0.2的个数
      Conherence_Window=Conherence_Window.*(CSD_Window<0.2);%把窗口内PSD小于0.2的标记出来，大于0.2的置零
      Gamma=sum(sum(Conherence_Window))./Num;%利用相干系数标准差小于0.2的点求出平均相干系数
      PSD_Window = PSD(ii-Radius+mins:ii+Radius-mins,jj-Radius+mins:jj+Radius-mins);%估计窗口内的PSD值
      
      %频率估计前预滤波(此处为改进1：:自适应预滤波，论文第56页4.3.1）
      % window_mean =mean_filter(phase_estimation, 3);
      %Gamma+0.2（因为计算Gamma时，剔除了3*3的估计窗口内标准差小于0.2的相干系数点所以要加回来）
      %Gamma+0.2，PSD_Window是为了根据相干系数和PSD求出滤波窗口的大小，与Win_size作比较，确定最终的滤波窗口大小；论文55页（4.4）（4.5）
      [window_mean,window(ii,jj)] = Mean_Filter_Correlation(phase_estimation,Gamma+0.2,PSD_Window,Win_size);
     
      %FFT频率估计（此处为改进2：相位主频估计，论文第57页4.3.2）
      phase_2fft = fftshift(fft2(window_mean,nn,nn)); %窗口内2DFFT
      [val_j,ind_j] = max(max(abs(phase_2fft)));      %ind_j最大值列号 
      [val_i,ind_i] = max(abs(phase_2fft(:,ind_j)));  %ind_i最大值行号 
      fi = double(ind_i-nn/2-2)/nn;
      fj = double(ind_j-nn/2-2)/nn;
      ei = exp(j*2*pi*fi);
      ej = exp(j*2*pi*fj);
      % CZT频谱主瓣精细估计
      phasef_czt1 = czt(window_mean,3*mm,w0,ei); %i(行)方向CZT,czt按列计算
      phasef_czt2 = czt(phasef_czt1.',3*mm,w0,ej);    %j(列)方向CZT
      phasef_czt2 = phasef_czt2.';
      [val_j,ind_j] = max(max(abs(phasef_czt2)));
      [val_i,ind_i] = max(abs(phasef_czt2(:,ind_j)));
      movei = double(ind_i-1)/mn;  %i(行)方向细化后偏移量
      movej = double(ind_j-1)/mn;  %j(列)方向细化后偏移量
      fii_1(ii,jj) = fi+movei;
      fjj_1(ii,jj) = fj+movej;       %CZT提升分辨率后，频谱峰值位置
      ei = exp(j*2*pi*fii_1(ii,jj)); %i(行)方向局部频率（指数）
      ej = exp(j*2*pi*fjj_1(ii,jj)); %j(列)方向局部频率（指数）
      % 线性相位
      phase_linear(ii-Radius:ii+Radius,jj-Radius:jj+Radius) = ei.^tempi*(ej.^tempj);%最终相位的一部分：估计的本地频率
      % 初相,sum(sum()) or mean(mean())
      %计算残留噪声相位：需要注意的是，估计出来的条纹频率是从原始含噪相位phase_update中移去，而不是预滤波之后的相位window_mean
      phase_0_matrix = phase_update(ii-Radius:ii+Radius,jj-Radius:jj+Radius).*...
                       conj(phase_linear(ii-Radius:ii+Radius,jj-Radius:jj+Radius));
%       residual_noise_phase(ii-Radius:ii+Radius,jj-Radius:jj+Radius) = phase_0_matrix;
      residual_noise_phase(ii,jj) = mean(mean(angle(phase_0_matrix)));

     %%残留相位滤波
      Conherence_Window2 = coherence_update(ii-Radius+mins:ii+Radius-mins,jj-Radius+mins:jj+Radius-mins);%窗口大小为3*3
      Gamma2=sum(sum(Conherence_Window2))./3/3;%Conherence_Window窗口大小为3*3
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
      
      %phase_0为最终相位的一部分：滤波后残留相位
      %(此处为改进3：:滤波器参数变化，论文第58页4.3.3）
%       phase_0 = Afilter_Goldstein_Correlation(phase_0_matrix, Gamma-abs(fii_2(ii,jj))+abs(fjj_2(ii,jj))); 
      phase_0 = Afilter_Goldstein_Correlation(phase_0_matrix,(1-Gamma2)*(1+sqrt((fii_2(ii,jj)^2+fjj_2(ii,jj)^2))/2)); %为什么相干系数可以改成这样子呢
%       phase_0 = Afilter_Goldstein_Correlation(phase_0_matrix, Gamma-1-abs(fii(ii,jj)+fjj(ii,jj))); 
      % 初相 + 线性相位 
      phase_filter(ii,jj) = phase_0 * phase_linear(ii,jj);%最终相位包含两个部分
  end
   str = ['滤波进行中',num2str(fix((100*(ii-Radius)/(Ni-2*Radius)))),'%'];
 waitbar(ii/(Ni-2*Radius),wb,str);
end
A = fspecial('average',[3 3]);    
phase_filter = angle(imfilter(phase_filter,A)); 
phase_filter = phase_filter(Radius+1:Ni-Radius,Radius+1:Nj-Radius);
end







