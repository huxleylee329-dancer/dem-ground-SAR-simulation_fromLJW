function [UPh_MLUPE] = MLUPE2(Psi_sub1,Psi_sub2,ratio_sub1,ratio_sub2,...
                   gamma_sub1,gamma_sub2,Phi_candidate)
% 最大似然相位解缠 -----------yyn 20161101
%% 基线个数
II = 2;
%% 三维数组
[Na,Nr] = size(Psi_sub1);
Psi     = zeros(Na,Nr,2);
ratio   = zeros(Na,Nr,2);
gamma   = zeros(Na,Nr,2);
%% 数据传递
Psi(:,:,1)   = Psi_sub1; % 短缠绕相位
Psi(:,:,2)   = Psi_sub2; % 长缠绕相位
ratio(:,:,1) = ratio_sub1; %短基线比
ratio(:,:,2) = ratio_sub2; %长基线比
gamma(:,:,1) = gamma_sub1; %短相干系数
gamma(:,:,2) = gamma_sub2; %长相干系数
%% 相干系数
rou = gamma;
%% 搜索区间 -- 全局搜索区间（计算速度较慢、暂不使用）
% Lower = min(Phi_candidate(:)); %解缠相位最小值
% Upper = max(Phi_candidate(:)); %解缠相位最大值
% search_step = 1e-2;            %搜索步长
% Phi =  Lower:search_step:Upper;%搜索区间
% JJ     = size(Phi,2);          %搜索相位个数
%% 干涉相位概率密度函数PDF
% display([repmat('-', 1, 10),' MLUPE ... ',repmat('-', 1, 11)])
UPh_MLUPE = zeros(Na,Nr);
for nn = 1:1:Na
    for mm = 1:1:Nr      
    % ---------------- 融合判断 -----------------20161128---- %   
    % 可在此部分增加策略
    % 如果短基线该点相干性低于阈值，则选择保留长基线该点数据
%       if rou(nn,mm,1) <= coherence_threshold_replace  % 20161128策略
%      if rou(nn,mm,1) <= rou(nn,mm,2)  % 20161213策略
%             UPh_MLUPE(nn,mm) = Psi(nn,mm,2);  
%       else%if rou(nn,mm,1) >=coherence_threshold_fusion  && rou(nn,mm,2) >= coherence_threshold_fusion  % 20161128策略
    % ----------------------------------------------------------- %  
    % ---------------- 搜索区间 --------------------------------- %  
        % 利用解缠结果实时更新搜索区间，搜索区间小，搜索速度快 
        Lower = Phi_candidate(nn,mm)-0.1*pi; 
        Upper = Phi_candidate(nn,mm)+0.1*pi; 
        search_step = 1e-3;   
        Phi =  Lower:search_step:Upper;
        JJ     = size(Phi,2);                 
    % ----------------------------------------------------------- %  
        % 最大似然函数
        MLF_PDF = ones(1,JJ);
        for ii = 1:1:II    %基线数
            for jj = 1:1:JJ%相位范围
                belta = rou(nn,mm,ii) * cos (Psi(nn,mm,ii) - Phi(jj)*ratio(nn,mm,ii));   
                %概率密度函数PDF 
                PDF_Phi_h(ii,jj) = ((1-rou(nn,mm,ii)^2)/2/pi) * (sqrt(1-belta^2) + belta*acos(-belta))...
                          / sqrt((1-belta^2))^3;             
            end
          % ------------------------------------------ %                    
%             if nn == 50 && mm == 50 && ii == 1       
%                pdf1 = PDF_Phi_h(ii,:);
%                figure;plot(Phi,pdf1);hold on
%             elseif nn == 50 && mm == 50 && ii == 2  
%                pdf2 = PDF_Phi_h(ii,:);
%                plot(Phi,pdf2)
%             end
          % ------------------------------------------ %    
        % 最大似然函数MLF
        MLF_PDF = MLF_PDF .* PDF_Phi_h(ii,:);
        end
        %最大似然估计MLE
        [~,index] = max(MLF_PDF);
        % 按比例换算回长基线
        UPh_MLUPE(nn,mm) = Phi(index) * ratio(nn,mm,2);
        % ------------------------------------------ % 
%         if nn == 50 && mm == 50  
%             mlf = MLF_PDF;
%             figure;plot(Phi,mlf)
%         end

     
        
        
   end
end