function [phase_res,local_phase_fx,local_phase_fy,phase_filter] = ...
    Local_fringe_frequency_estimation(wrapped_phase,MPG,PDV,coherence)
% 基于相干系数和相位导数方差自适应改变局部条纹频率估计窗口大小
% FFT+CZT 提升局部频率估计的精度
% 对局部条纹频率估计值进行正确性判断

%% 计算最大相位梯度和相位导数方差的均值、标准差、最大值、最小值，确定估计窗口大小
% hist(MPG);%直方图
mean_MPG = mean2(MPG);%均值
std_MPG = std2(MPG);%标准差
MPGradius = ceil(abs(MPG - mean_MPG)./std_MPG);%MPG半径

% % hist(PDV);%直方图
% mean_PDV = mean2(PDV);%均值
% std_PDV = std2(PDV);%标准差
% PDVradius = ceil(abs(PDV - mean_PDV)./std_PDV);%PDV半径

coherence_radius = ceil(1./coherence);%相干系数半径
coherence_radius(coherence_radius > 6) = 6;
coherence_radius(coherence_radius < 5) = 5;

% CF_radius = MPGradius + PDVradius + coherence_radius;%复杂度因子确定的半径
CF_radius = MPGradius + coherence_radius;%复杂度因子确定的半径


%%

j=sqrt(-1);
[Size_i,Size_j]=size(wrapped_phase);

%% 扩充估计相位矩阵
winsize_e_i = 41;           %扩充尺寸
winsize_e_j = 41;            %扩充尺寸
wine_i = (winsize_e_i-1)/2; %扩充半径
wine_j = (winsize_e_j-1)/2; %扩充半径
phase_update = zeros(Size_i+winsize_e_i-1,Size_j+winsize_e_j-1);
[Ni,Nj] = size(phase_update);
phase_update(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j) = wrapped_phase;

CF1_radius_update = zeros(Size_i+winsize_e_i-1,Size_j+winsize_e_j-1);
CF1_radius_update(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j) = CF_radius;
phase_res = zeros(Ni,Nj);%移除条纹频率后的残留相位
phase_filter = zeros(Ni,Nj);%滤波相位
local_phase_fx = zeros(Ni,Nj);%估计的条纹频率
local_phase_fy = zeros(Ni,Nj);%估计的条纹频率

% flipud:矩阵上下交换，fliplr:矩阵左右交换
phase_update(1:wine_i,1:Nj) = flipud(phase_update(wine_i+1:2*wine_i,1:Nj)); 
phase_update(1:Ni,1:wine_j) = fliplr(phase_update(1:Ni,wine_j+1:2*wine_j));  
phase_update(Ni-wine_i+1:Ni,1:Nj) = flipud(phase_update(Ni-2*wine_i+1:Ni-wine_i,1:Nj));
phase_update(1:Ni,Nj-wine_j+1:Nj) = fliplr(phase_update(1:Ni,Nj-2*wine_j+1:Nj-wine_j));

CF1_radius_update(1:wine_i,1:Nj) = flipud(CF1_radius_update(wine_i+1:2*wine_i,1:Nj)); 
CF1_radius_update(1:Ni,1:wine_j) = fliplr(CF1_radius_update(1:Ni,wine_j+1:2*wine_j));  
CF1_radius_update(Ni-wine_i+1:Ni,1:Nj) = flipud(CF1_radius_update(Ni-2*wine_i+1:Ni-wine_i,1:Nj));
CF1_radius_update(1:Ni,Nj-wine_j+1:Nj) = fliplr(CF1_radius_update(1:Ni,Nj-2*wine_j+1:Nj-wine_j));

% 将相位转化为复数
phase_update = exp(j*phase_update); 

%% 局部条纹频率估计参数
nn = 32;                     % 2D FFT点数
mm = 32;                     % CZT点数
mn = mm*nn;                  % 频谱总分辨率（mn*mn）
w0 = exp(-j*2*pi/mn);        % CZT参数


%% 自适应局部条纹频率估计
wb = waitbar(0,'Local fringe frequency estimation ...');
for ii = wine_i+1:Ni-wine_i
  for jj = wine_j+1:Nj-wine_j
      
      %根据复杂度因子矩阵自适应改变局部条纹频率估计窗口大小,选择不同的滤波方法对残留相位进行滤波                
          
       if CF1_radius_update(ii,jj) == 2
          adap_Radius_i = CF1_radius_update(ii,jj);
          adap_Radius_j = CF1_radius_update(ii,jj);
          adap_winsize_e_i = 2*adap_Radius_i + 1;
          adap_winsize_e_j = 2*adap_Radius_j + 1;
          
      elseif CF1_radius_update(ii,jj) == 1
          adap_Radius_i = mid_radius;
          adap_Radius_j = mid_radius;
          adap_winsize_e_i = 2*adap_Radius_i + 1;
          adap_winsize_e_j = 2*adap_Radius_j + 1;    
          
      else 
          adap_Radius_i = min_radius;
          adap_Radius_j = min_radius;
          adap_winsize_e_i = 2*adap_Radius_i + 1;
          adap_winsize_e_j = 2*adap_Radius_j + 1;    
      end
      
%           %非线性相位
%           i1 = phase_estimation;
%           i2 = fftshift(fft2(i1)); %窗口内2DFFT
%           i3 = abs(i2);
%           Max=max(max(abs(i3)));     %频谱峰值
%                         for i=1:adap_winsize_e_i
%                             for j=1:adap_winsize_e_j
%                                if i3(i,j)<(Max/5)
%                                   i2(i,j)=0;
%                                end
%                             end
%                         end
%           local_phase = ifft2(ifftshift(i2));          
%           %移除条纹频率后的残留相位
%           phase_0_matrix = phase_update(ii-adap_Radius_i:ii+adap_Radius_i,jj-adap_Radius_j:jj+adap_Radius_j).*...
%                        conj(local_phase);
          
          %线性相位拟合
          phase_estimation(1:adap_winsize_e_i,1:adap_winsize_e_j) = phase_update(ii-adap_Radius_i:ii+adap_Radius_i,jj-adap_Radius_j:jj+adap_Radius_j);
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
          tempi = (0:adap_winsize_e_i-1)';  % 行线形相位网格
          tempj = 0:adap_winsize_e_j-1;     % 列线形相位网格
          local_phase = ei.^tempi*(ej.^tempj);
          %移除条纹频率后的残留相位
          phase_0_matrix = phase_update(ii-adap_Radius_i:ii+adap_Radius_i,jj-adap_Radius_j:jj+adap_Radius_j).*...
                       conj(local_phase);
          
          %估计的条纹频率
          local_phase_fx(ii,jj) = fii;
          local_phase_fy(ii,jj) = fjj;
          

     %% 残留相位滤波     
%       phase_res(ii,jj) = mean(mean(phase_0_matrix));
      phase_res(ii,jj) = phase_0_matrix(adap_Radius_i+1,adap_Radius_j+1);
     
      % 初相 + 线性相位 
      phase_filter(ii,jj) = phase_res(ii,jj) * local_phase(adap_Radius_i+1,adap_Radius_j+1);
              
  end
str = ['Local fringe frequency estimation ...',num2str(fix((100*(ii-wine_i)/(Ni-2*wine_i)))),'%'];
waitbar(ii/(Ni-2*wine_i),wb,str);
% pause(1);
% close(wb);
end

%% 取复相角(rad),取有效范围
phase_res = angle(phase_res);
phase_res = phase_res(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j);
phase_filter = angle(phase_filter);
phase_filter = phase_filter(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j);

local_phase_fx = local_phase_fx(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j);
local_phase_fy = local_phase_fy(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j);

end

% function Array_Output = Mean_Filter_Correlation(Array,win_pre)
%             A=fspecial('average',[win_pre win_pre]);    
%             Array_Output=imfilter(Array,A); 
% end
