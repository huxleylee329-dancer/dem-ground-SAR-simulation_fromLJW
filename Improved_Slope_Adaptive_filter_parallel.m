
function [phase_filter,phase_res ] = Improved_Slope_Adaptive_filter_parallel(Original_Interferogram, Win_size,win_pre)
%% %%空间自适应滤波
%窗参数
% Win_size=17;  
Radius = (Win_size-1)/2; %估计窗i方向半径
[Size_i,Size_j]=size(Original_Interferogram);
%扩充估计相位矩阵
phase_update = zeros(Size_i+Win_size-1,Size_j+Win_size-1);    
[Ni,Nj] = size(phase_update);
phase_res = zeros(Ni,Nj);%移除条纹频率后的残留相位
phase_update(Radius+1:Ni-Radius,Radius+1:Nj-Radius) = Original_Interferogram;
%flipud:矩阵上下交换，fliplr:矩阵左右交换  %%对半径行和列进行边缘赋值  %%目的是弱化多视对边缘的影响 
%对边缘的7行7列加值
phase_update(1:Radius,1:Nj) = flipud(phase_update(Radius+1:2*Radius,1:Nj)); 
phase_update(1:Ni,1:Radius) = fliplr(phase_update(1:Ni,Radius+1:2*Radius));  
phase_update(Ni-Radius+1:Ni,1:Nj) = flipud(phase_update(Ni-2*Radius+1:Ni-Radius,1:Nj));
phase_update(1:Ni,Nj-Radius+1:Nj) = fliplr(phase_update(1:Ni,Nj-2*Radius+1:Nj-Radius));

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
mins=9;
%%
    wb = waitbar(0,'phase filtering ...');
%      parpool  8         %-------------根据计算环境的workers数目而定
% parpool('local',12);
for ii = Radius+1:Ni-Radius
    
  for jj = Radius+1:Nj-Radius
      phase_estimation = phase_update(ii-Radius:ii+Radius,jj-Radius:jj+Radius);
%       %频率估计前预滤波
%       window_mean = Mean_Filter_Correlation(phase_estimation,win_pre);
%       %FFT频率估计
%       phase_2fft = fftshift(fft2(window_mean,nn,nn)); %窗口内2DFFT

      phase_2fft = fftshift(fft2(phase_estimation,nn,nn)); %窗口内2DFFT
      [val_j,ind_j] = max(max(abs(phase_2fft)));      %ind_j最大值列号 
      [val_i,ind_i] = max(abs(phase_2fft(:,ind_j)));  %ind_i最大值行号 
      fi = double(ind_i-nn/2-2)/nn;
      fj = double(ind_j-nn/2-2)/nn;
      ei = exp(j*2*pi*fi);
      ej = exp(j*2*pi*fj);
      % CZT频谱主瓣精细估计（CZT指离散傅里叶变换）
     % phasef_czt1 = czt(window_mean,3*mm,w0,ei); %i(行)方向CZT,czt按列计算
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
      AA = ei.^tempi*(ej.^tempj);
      % 初相,sum(sum()) or mean(mean())
      phase_0_matrix = phase_update(ii-Radius:ii+Radius,jj-Radius:jj+Radius).*...
                       conj(AA);
     %% 残留相位滤波
      phase_res(ii,jj) =  mean(mean(phase_0_matrix)); %均值滤波
          %这样效果不好
%           A = fspecial('average',[3 3]);   
%           Array_Output = imfilter(phase_0_matrix,A);
%           phase_res(ii,jj) =  Array_Output((end-1)/2,(end-1)/2); %均值滤波  

      % 初相 + 线性相位 
      phase_filter(ii,jj) = phase_res(ii,jj) * AA((end-1)/2,(end-1)/2);
  end
 str = ['滤波进行中',num2str(fix((100*(ii-Radius)/(Ni-2*Radius)))),'%'];
 waitbar(ii/(Ni-2*Radius),wb,str);
end

% close(wb);
phase_res = angle(phase_res);
phase_res = phase_res(Radius+1:Ni-Radius,Radius+1:Nj-Radius);
phase_filter = angle(phase_filter);
phase_filter = phase_filter(Radius+1:Ni-Radius,Radius+1:Nj-Radius);
%  delete(gcp);

end

function Array_Output = Mean_Filter_Correlation(Array,win_pre)
            A=fspecial('average',[win_pre win_pre]);    
            Array_Output=imfilter(Array,A); 
end







