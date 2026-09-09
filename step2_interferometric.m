% ----------------
% ---------------- step 2
% ---------------- 长短基线相干配准 

clear all
close all
clc
%% 输入处理ID
process_ID_master  = '500_mast';
process_ID_slave   = '500_slave';
process_ID = [process_ID_master,'&',process_ID_slave,'_insar'];
%% 载入SSC
load(['/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/output/',process_ID_master,'/main.mat'],'MainSLC1');
ssc_master = MainSLC1;
load(['/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/output/',process_ID_slave,'/slave.mat'],'SlaveSLC1');
ssc_slave  = SlaveSLC1;
%经过验证这两幅主辅图像是山地，没有平行条纹，不符合基线估计要求。
% load(['/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/output/',process_ID_master,'/Image_main.mat'],'Image_main');
% ssc_master = Image_main;
% load(['/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/output/',process_ID_slave,'/Image_slave.mat'],'Image_slave');
% ssc_slave  = Image_slave;

%% 裁剪至相同大小
% 取主辅图像公共区域
[nr1,nc1]  = size(ssc_master);
[nr2,nc2]  = size(ssc_slave);
nr         = min(nr1,nr2);
nc         = min(nc1,nc2);
SSC_master = ssc_master(1:nr,1:nc);
SSC_slave  = ssc_slave(1:nr,1:nc);
% 更新主图像裁剪范围
[~] = cut_image_index_update(nr1,nc1,nr,nc,process_ID_master);
[~] = cut_image_index_update(nr2,nc2,nr,nc,process_ID_slave);
%% 干涉处理
% 粗配准
[nr1,nc1]            = size(SSC_master);
[SSC_master_regis1,SSC_slave_regis1,move_r,move_c] = registration_pixel(SSC_master,SSC_slave,nr1,nc1);
[nr2,nc2]            = size(SSC_master_regis1);
% 更新主图像裁剪范围
[~] = after_regis_image_index_update(move_r,move_c,process_ID_master,process_ID_slave);

% 亚像素级配准
[~,SSC_slave_regis2] = regis_subpixel(SSC_master_regis1,SSC_slave_regis1,nr2,nc2);
coherence = Calculation_Coherence_Coefficient_complex_no_phase(SSC_master_regis1,SSC_slave_regis2);
folder=['/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/Mid/',process_ID,'/'];
save([folder,'coherence.mat'],'coherence');
save([folder,'SSC_master_regis1.mat'],'SSC_master_regis1');
save([folder,'SSC_slave_regis2.mat'],'SSC_slave_regis2');
% 缠绕干涉相位
wrapped_phase        = interferometric_phase(SSC_master_regis1,SSC_slave_regis2);

Data=SSC_master_regis1.*conj(SSC_slave_regis2);
Data=Data(160,:);

save([folder,'wrapped_phase.mat'],'wrapped_phase');
% imagesc(wrapped_phase);
% % 干涉相位滤波
%    win_size  = 25;
%      win_pre = 23;
%      tic;
%      [phase_removeFlat,Flat_phase_wrapped] = removeFlat_range_20161126(wrapped_phase); %移除平地相位
%      phase_removeFlat_filter = Improved_Slope_Adaptive_filter_parallel(phase_removeFlat,win_size,win_pre); %坡度自适应滤波
%      save ([folder,'phase_removeFlat_filter.mat'], 'phase_removeFlat_filter');    %保存残余相位
%      save ([folder,'Flat_phase_wrapped.mat'], 'Flat_phase_wrapped');  %保存平地相位
%      wrapped_phase_filter_gray        = mat2gray(phase_removeFlat_filter); 
%      imwrite(wrapped_phase_filter_gray,[folder,'phase_removeFlat_filter.tif']);
%      toc
%% 
% % clear all
% 利用FFT和插值计算频率粗估计值fi0
% load(['/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/Code/','/wrapped_phase.mat'],'wrapped_phase');
load(['/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/Code/','/range_main.mat'],'range_main');
slantrange=range_main;
%  wrapped_phase=mean(wrapped_phase);
% figure(1);
% imagesc(wrapped_phase);
sclomn=30;%初始列
eclomn=250;%终止列
rmin=slantrange(sclomn);%初始列对应的是最小斜距，r代表与图像像素点对应的斜距
rmax=slantrange(eclomn);%终止列对应的是最大斜距
indey = zeros(1, eclomn - sclomn + 1);
for i=sclomn:eclomn%对平行条纹某一行的第sclomn-eclomn个像素点求解每一个像素点的粗估计值
    N0=15;%N0表示以第i个像素点为中心，选取2N0+1大小的一维邻域进行FFT
    x(i-sclomn+1)=i-sclomn+1;%x代表选定平行条纹距离向的点数
    cut_wrapped_phase=wrapped_phase(160,i-N0:i+N0);%对该像素点2N0+1大小的一维邻域进行FFT
%     cut_wrapped_phase=sum(cut_wrapped_phase);%对该像素点2N0+1大小的一维邻域进行FFT

    %map_wrapped_phase=mapminmax(cut_wrapped_phase,0,1);%对数据幅度进行归一化
    fs=1;%采样频率
    N=size(cut_wrapped_phase,2);
    y1=fft(cut_wrapped_phase,2048);%fft变换
    %y=fftshift(y1);
    y=y1.*conj(y1)/N;%求功率谱密度
    
    %mag=abs(y)*(2/N);%实际幅值变换
    f=(0:length(y)-1)/length(y)*fs;%要得到真实频率的下标，先将频率归一化，再乘以采样率
% %     figure(1);
% %     plot(f(2:N/2),y(2:N/2),'g');%做频谱图
% %     xlabel('频率(1/m)');
% %     ylabel('幅值(1/m)');
% %     title('干涉复图像 幅频谱图');
%     
%     %对峰值点附近进行插值
     [max_y,indey(i-sclomn+1)]=max(y(1:end/2));%寻找频谱图中的峰值点，max_y表示峰值点，indey表示峰值点所在的位置
     
%     
%     
%     xi=f(2):0.5/(2*N):f(end);%以峰值点的x坐标为中心，设置插值范围，减小该范围内x坐标的分度值
%     y2=interp1(f(2:end),y(2:end),xi,'method');%得到y坐标峰值点附近的插值
%     [max_y2,indey2]=max(y2(2:end),[],2);%找到插值后y2的最大值所在的位置，indey2表示x坐标所在位置。
%     fi0(x(i-sclomn+1))=xi(indey2);%rou_f0(x)保存第x=i-sclomn-1个像素点的频率粗估计值
%     figure(2);
%     plot(fi0,'r');
%     figure(2);
%     plot(xi,y2,'r');%做频谱图
%     xlabel('频率(1/m)');
%     ylabel('幅值(1/m)');
%     title('干涉复图像 幅频谱图');
%     
%     %利用3dB带宽求截止频率
%     %对峰值点附近插值后的数据根据3dB带宽的定义得到频率粗估计值（3dB带宽通常指功率谱密度的最高点下降到1/2时界定的频率范围）
%     %[b,i]=sort(a),sort函数给数列a从小到大排列,b为从小到大的数字，i为对应位置
%     [b_y2,i_y2]=sort(abs(y2-max_y/2));%所以b_y2中的前两位表示最接近半功率点的y坐标，i_y2中前两位是最接近半功率点的y坐标对应的x坐标
%     rou_k=xi(i_y2(1):i_y2(2));
end
 f0=indey/(2048)*fs;
 f0=ones(size(f0))/26;
%% 利用FFT和CZT估计局部频率
load(['/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/Code/','/range_main.mat'],'range_main');
slantrange=range_main;
% figure(1);
% imagesc(wrapped_phase);
sclomn=30;%初始列
eclomn=250;%终止列
rmin=slantrange(sclomn);%初始列对应的是最小斜距，r代表与图像像素点对应的斜距
rmax=slantrange(eclomn);%终止列对应的是最大斜距
for i=sclomn:eclomn%对平行条纹某一行的第sclomn-eclomn个像素点的局部频率
    N0=9;%N0表示以第i个像素点为中心，选取N0*N0大小的二维邻域作为原始相位
    radius_N0=floor(N0/2);%表示原始相位窗口半径
    phase_initial=wrapped_phase(160-radius_N0:160+radius_N0,i-radius_N0:i+radius_N0);
    window=15;%估计窗口大小
    [localFrequency_range,phase_filter] = calculate_localFreRange(phase_initial, window);
    [mid_ind]=floor((size(localFrequency_range)+1)/2);
    localFre(i-sclomn+1)=localFrequency_range(mid_ind(1),mid_ind(2));
    
end
 %% 半牛顿迭代法估计干涉条纹频率
% T=1;%表示距离向采样间隔
% for i=sclomn:eclomn
%     fi_t(1)=fi0(i-sclomn+1);%将第i个像素点的粗频率估计值作为迭代的初始值，经历30次迭代
% for t=1:10
%     n=-N0:1:N0;%n代表一维邻域长度
%     syms f     %定义一个符号变量
%     syms phi
%     [H_f] = BND_H_f(n,T,fi_t(t),f,phi);%得到函数H_f符号表达式
%     diff_H_f1=diff(H_f,f);%对H_f求一阶导
%     diff_H_f2=diff(diff_H_f1,f);%对H_f求二阶导
%     H_f_1 = subs(diff_H_f1, f, fi_t(t));%R = subs(S, old, new) 利用new的值代替符号表达式中old的值
%     H_f_2 = subs(diff_H_f2, f, fi_t(t));%代入f=fi_t(t)
%     fi_t(t+1)=fi_t(t)-real(H_f_1/H_f_2);
% end
%     k(i-sclomn+1)=2*pi*fi_t(t+1);%k()保存条纹频率值
% end
%% 半牛顿迭代法估计干涉条纹频率
f_i = zeros(1, eclomn - sclomn + 1);
% [f_out] = iteration_semi_Newton(iter,f0,N,ii,Data,T);
for ii = sclomn:eclomn
    f_i(1,ii-sclomn+1) = iteration_semi_Newton(500,localFre(ii-sclomn+1),N0,ii,Data,1);
%     close all
end





%% 线性拟合抑制条纹频率噪声，得到条纹频率估计值
%P=polyfit(x,y,N);  %N多项式拟合函数，返回降幂排列的多项式系数
%yi=polyval(P,xi);  %计算以P向量为系数的多项式在xi处的值
t=polyfit(x,k,1);%对频率估计值进行线性拟合，横坐标为选取区域像素点个数，纵坐标为半牛顿迭代法求解出的条纹频率值
kmax=polyval(t,x(1));%计算以t向量为系数的多项式在x(1)处的值，即条纹频率估计值最大值
kmin=polyval(t,x(end));%计算以t向量为系数的多项式在x(end)处的值，即条纹频率估计值最小值
%% 分别计算平地模型和考虑地球曲率的基线估计公式
syms r     %r为一个符号变量
d=0.5;%d=1，收发分置InSAR系统，取值为0.5，收发同置InSAR系统，取值为1
Re=6378.16e3;%地球半径
load('/Users/thorfun/Documents/实验室/23所基线项目/5-基于干涉相位基线估计/基于干涉相位基线估计程序/Code/MainPos.mat') %卫星坐标（直接当做已知飞机坐标）
H=sqrt(MainPos(:,:,1).^2+MainPos(:,:,2).^2+MainPos(:,:,3).^2)-Re;%卫星/机载高度
Wlength=3*10^8/(9.6*10^9);%波长
K_B=zeros(2,1);
K_B(1,1)=kmax;
K_B(2,1)=kmin;
%平地模型
p_r=(4*pi*d*H)/(Wlength*r);
g_r=H/sqrt(r^2-H^2);
B_flat=zeros(2,2);%B_flat代表平地模型下的基线估计公式
p_rmin=subs(p_r,r,rmin);%给p_r代入斜距最小值
p_rmax=subs(p_r,r,rmax);%给p_r代入斜距最大值
g_rmin=subs(g_r,r,rmin);%给g_r代入斜距最小值
g_rmax=subs(g_r,r,rmax);%给g_r代入斜距最大值
B_flat(1,1)=p_rmin*g_rmin;
B_flat(1,2)=p_rmin;
B_flat(2,1)=p_rmax*g_rmax;
B_flat(2,2)=p_rmax;
B1=zeros(2,1);
B1=inv(B_flat)*K_B;%得到平地模型下的基线分量
B1x=B1(1,1);
B1y=B1(2,1);
flat_B=sqrt(B1x^2+B1y^2);%平地模型基线长度
flat_alfa=atan(B1y/B1x);%平地模型基线倾角
%地球曲率
q_r=(4*pi*d*(H^2+2*H*Re)^2-r^4)/(Wlength*r^2*(H+Re)*sqrt(((2*Re+H)^2-r^2)*(r^2-H^2)));
f_r=-((4*pi*d)*(r^2-H^2-2*H*Re))/(Wlength*2*r^2*(H+Re));
B_ear=zeros(2,2);%B_flat代表平地模型下的基线估计公式
q_rmin=subs(q_r,r,rmin);%给q_r代入斜距最小值
q_rmax=subs(q_r,r,rmax);%给q_r代入斜距最大值
f_rmin=subs(f_r,r,rmin);%给f_r代入斜距最小值
f_rmax=subs(f_r,r,rmax);%给f_r代入斜距最大值
B_ear(1,1)=q_rmin;
B_ear(1,2)=f_rmin;
B_ear(2,1)=q_rmax;
B_ear(2,2)=f_rmax;
B2=zeros(2,1);
B2=inv(B_ear)*K_B;%得到地球曲率下的基线分量
B2x=B2(1,1);
B2y=B2(2,1);
ear_B=sqrt(B2x^2+B2y^2);%地球曲率基线长度
ear_alfa=atan(B2y/B2x);%地球曲率基线倾角










