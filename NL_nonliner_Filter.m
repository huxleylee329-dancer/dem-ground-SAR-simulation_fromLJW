function phase_filter = NL_nonliner_Filter(phase,real_phase)
% 非局域非线性滤波

j=sqrt(-1);
% %----------------- display information -------------------%
% display([repmat('-', 1, 8),' phase filter ',repmat('-', 1, 8)])
% display(['method:             ','SAPF'])
% display(repmat('-', 1, 30))
%% 窗参数

winsize_e_i = 23;           %估计窗i方向尺寸
winsize_e_j =23;           %估计窗j方向尺寸
winsize_f_i = 21;  %滤波窗i方向尺寸
winsize_f_j = 21;  %滤波窗j方向尺寸0104修改
wine_i = (winsize_e_i-1)/2; %估计窗i方向半径
wine_j = (winsize_e_j-1)/2; %估计窗j方向半径

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

real_update=padarray(real_phase,[wine_i,wine_j],'symmetric','both');

% 将相位转化为复数
phase_update = exp(j*phase_update);  
real_update = exp(j*real_update);  
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
%% 非局域参数LS2019
% wb = waitbar(0,'phase filtering ...');
for ii = wine_i+1:Ni-wine_i
  for jj = wine_j+1:Nj-wine_j
      phase_estimation(1:winsize_e_i,1:winsize_e_j) = phase_update(ii-wine_i:ii+wine_i,jj-wine_j:jj+wine_j);
      i1=real_update(ii-wine_i:ii+wine_i,jj-wine_j:jj+wine_j);
      % FFT频率估计
%       phase_2fft = fftshift(fft2(phase_estimation,nn,nn)); %窗口内2DFFT
      
      i2 = fftshift(fft2(i1)); %窗口内2DFFT 
        i3=abs(i2);
        Max=max(max(abs(i3)));     %频谱峰值
                        for i=1:winsize_e_i
                            for j=1:winsize_e_j
                               if i3(i,j)<(Max/5)
                                  i2(i,j)=0;
                               end
                            end
                        end
%  i4=ifft2(ifftshift(i2));%非线性相位
%  res_phase = angle(image.*conj(i4));%残余相位
      
 
      % 线性相位
      phase_nonlinear = ifft2(ifftshift(i2));
      % 初相,sum(sum()) or mean(mean())
      phase_0_matrix = phase_update(ii-wine_i:ii+wine_i,jj-wine_j:jj+wine_j).*...
                       conj(phase_nonlinear);
%       %均值滤波
%        phase_0 =  mean(mean(phase_0_matrix)); %均值滤波  
      %NL滤波
      ds=3;%邻域窗口半径
      Ds=wine_i;%搜索窗口半径
      h=0.7;%高斯函数平滑参数
      kernel=ones(2*ds+1,2*ds+1);
      kernel=kernel./((2*ds+1)*(2*ds+1));
      h2=h*h;

%%%%%%%%%%%
%         W1= PaddedImg(ii-ds:ii+ds,jj-ds:jj+ds);%邻域窗口1
        W1=phase_0_matrix(Ds-ds+1:Ds+ds+1,Ds-ds+1:Ds+ds+1) ;%邻域窗口1
        wmax=0;
        average=0;
        sweight=0;
        %%搜索窗口
        rmin = ds+1;
        rmax = 2*Ds+1-ds;
        smin = ds+1;
        smax = 2*Ds+1-ds;
        for r=rmin:rmax
            for s=smin:smax
                W2=phase_0_matrix(r-ds:r+ds,s-ds:s+ds);%邻域窗口2
                 Dist2=sum(sum(kernel.*(W1-W2).*(W1-W2)));%邻域间距离
%                   Dist2=sum(sum(kernel.*(angle(W1)-angle(W2)).*(angle(W1)-angle(W2))));%邻域间距离
                w=exp(-Dist2/h2);
                if(w>wmax)
                    wmax=w;
                end
                sweight=sweight+w;
                average=average+w*phase_0_matrix(r,s);
            end
        end
%         average=average+wmax*phase_0_matrix(ii,jj);%自身取最大权值
%         sweight=sweight+wmax;
       phase_0=average/sweight;


      % 初相 + 线性相位 
      phase_filter(ii,jj) = phase_0 * phase_nonlinear(wine_i+1,wine_i+1);
%       Local_frequency(ii,jj)=fjj;
  end
%    str = ['nl2滤波进行中',num2str(fix((100*ii/(Ni-2*wine_i)))),'%'];
%  waitbar(ii/(Ni-2*wine_i),wb,str);
%  pause(1);
   if mod(ii,10)==0
       ii
   end
end

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