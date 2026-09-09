function [UPh_MLUPE] = MLUPE_PDV(Psi_sub1,Psi_sub2,short_PDV,long_PDV,ratio_sub1,ratio_sub2,...
    gamma_sub1,gamma_sub2,B_short,B_long,idx_ref)
% 融合判断-根据相位导数方差判断是选择替换还是最大似然融合 -----------wy20210323
UPh_MLUPE = Psi_sub2;%返回的融合相位，初始值默认为长基线解缠相位
B_ratio = B_long/B_short;%将短基线转换到长基线相位
%% 基线个数
% II = 2;
%% 三维数组
[Na,Nr] = size(Psi_sub1);
Psi     = zeros(Na,Nr,2);
ratio   = zeros(Na,Nr,2);
gamma   = zeros(Na,Nr,2);
Phi_candidate = zeros(Na,Nr);
%% 数据传递
Psi(:,:,1)   = Psi_sub1; % 短缠绕相位
Psi(:,:,2)   = Psi_sub2; % 长缠绕相位
ratio(:,:,1) = ratio_sub1; %短基线比
ratio(:,:,2) = ratio_sub2; %长基线比
gamma(:,:,1) = gamma_sub1; %短相干系数
gamma(:,:,2) = gamma_sub2; %长相干系数
%% 相干系数
rou = gamma;
%% 干涉相位概率密度函数PDF
rou1 = reshape(rou(:,:,1),Na,Nr);
rou2 = reshape(rou(:,:,2),Na,Nr);
Psi1 = reshape(Psi(:,:,1), Na, Nr);
Psi2 = reshape(Psi(:,:,2), Na, Nr);

%找出短基线相位梯度大的地方并标记
dz_s = pdv_jump_mark(short_PDV,Psi_sub1);

%找出长基线相位梯度大的地方并标记
dz_l = pdv_jump_mark(long_PDV,Psi_sub2);

b = 1;%扩大的搜索区间倍数
% p_r = 5;%找一点画概率密度函数
% p_c = 5;
h = waitbar(0,'粗融合中...');
for ii = 1:1:Na-1
       nn = ii;
       str = ['粗融合中...',num2str(ii/(Na-1)*100),'%'];
       waitbar(ii/(Na-1),h,str)    
    for jj = 1:1:Nr-1
        %逐点判断长短基线的跳变点：
        %如果长短基线都是跳变点或者都不是跳变点，则以长基线为参考基线，长基线解缠相位为搜索中心进行最大似然融合
        %如果短基线不是跳变点，长基线是跳变点，则将短基线转换到长基线域
        %如果短基线是跳变点，长基线不是跳变点，则保留长基线该点的相位
        mm = jj;
        %如果短基线不是跳变点，长基线是跳变点，则将短基线转换到长基线域
        if dz_s(ii,jj)==0 && dz_l(ii,jj)==1
            UPh_MLUPE(nn,mm) = Psi_sub1(nn,mm)*B_ratio;
        end
        
        %如果短基线是跳变点，长基线不是跳变点，则保留长基线该点的相位
        if dz_s(ii,jj)==1 && dz_l(ii,jj)==0
            UPh_MLUPE(nn,mm) = Psi_sub2(nn,mm);
        end  
        
        %如果长短基线都是跳变点或者都不是跳变点，则以长基线为参考基线，长基线解缠相位为搜索中心进行最大似然融合
        if (dz_s(ii,jj)==1 && dz_l(ii,jj)==1) || (dz_s(ii,jj)==0 && dz_l(ii,jj)==0) 
        Phi_candidate(nn,mm) = Psi_sub2(nn,mm);
        Lower = Phi_candidate(nn,mm)-b*1*pi; 
        Upper = Phi_candidate(nn,mm)+b*1*pi; 
        search_step = 0.5*5e-2;   
        Phi =  Lower:search_step:Upper;
        JJ     = size(Phi,2);
        % 求最大似然函数
        MLF_PDF = ones(1,JJ);
        c1 = (1 - rou1(nn, mm)^2) / 2 / pi;
        c2 = (1 - rou2(nn, mm)^2) / 2 / pi;
        beta1 = rou1(nn, mm) * cos(Psi1(nn,mm) - Phi*ratio(nn,mm,1));
        beta2 = rou2(nn, mm) * cos(Psi2(nn,mm) - Phi*ratio(nn,mm,2));
        pdf1 = c1 * (sqrt(1 - beta1.^2) + beta1.*acos(-beta1)) ./ ((1 - beta1.^2).^1.5);
        pdf2 = c2 * (sqrt(1 - beta2.^2) + beta2.*acos(-beta2)) ./ ((1 - beta2.^2).^1.5);
        % 最大似然函数MLF
        MLF_PDF = MLF_PDF .* pdf1 .*pdf2;%将pdf1和pdf2相乘
        [~,index] = max(MLF_PDF);       
        % 按比例换算回参考基线
        UPh_MLUPE(nn,mm) = Phi(index) * ratio(nn,mm,idx_ref);% 按比例换算回参考基线

%         if nn == p_r && mm == p_c     
%              figure;plot(Phi,pdf1,'r--');hold on 
%              plot(Phi,pdf2,'b')
%         end
%         if nn == p_r && mm == p_c  
%             mlf = MLF_PDF;
%             figure;plot(Phi,mlf)
%         end 

        end
   end
end
end
