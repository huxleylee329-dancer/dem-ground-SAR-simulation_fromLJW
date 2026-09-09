function [sar_image_cut,fsmr,fsmc] = cut_multi_images_to_regis(move_r,move_c,sar_image,Mast_ind) 
% move_r : 若>0,表示辅图像相对于主图像向上偏移；若<0,表示辅图像相对于主图像向下偏移。
% move_c : 若>0,表示辅图像相对于主图像向左偏移；若<0,表示辅图像相对于主图像向右偏移。
% sar_image是再次截取的sar图像
M = length(sar_image) - 1; 
Slave_ind = [1:Mast_ind-1 Mast_ind+1:M+1];
fsmr = zeros(1,M+1);%辅图像相对于主图像主图像第一行开始的行数
fsmc = zeros(1,M+1);%
[nr,nc] = size(sar_image{Mast_ind});
s_m_r   = sort(move_r);%默认进行升序排列，每幅辅图像一个move_r
s_m_c   = sort(move_c);
sar_image_temp = cell(1,M);
sar_image_cut  = cell(1,M+1);

if s_m_r(1) >= 0 && s_m_r(end) >= 0 %所有辅图像都相对于主图像向上偏移
    master_temp = sar_image{Mast_ind}(1:nr-s_m_r(end),:);
    fsmr(Mast_ind) = 1;
    for i = 1:length(Slave_ind)
        sar_image_temp{i}(1:nr-s_m_r(end),:) = sar_image{Slave_ind(i)}(1+move_r(Slave_ind(i)):nr+move_r(Slave_ind(i))-s_m_r(end),:);
        fsmr(Slave_ind(i)) = 1+move_r(Slave_ind(i));
    end
elseif s_m_r(1) < 0  && s_m_r(end) <= 0
    master_temp = sar_image{Mast_ind}(1-s_m_r(1):nr,:);
    fsmr(Mast_ind) = 1-s_m_r(1);
    for i = 1:length(Slave_ind)
        sar_image_temp{i}(1:nr+s_m_r(1),:) = sar_image{Slave_ind(i)}(1+move_r(Slave_ind(i))-s_m_r(1):nr+move_r(Slave_ind(i)),:);
        fsmr(Slave_ind(i)) = 1+move_r(Slave_ind(i))-s_m_r(1);
    end
elseif s_m_r(1) < 0 && s_m_r(end) > 0
    master_temp = sar_image{Mast_ind}(1-s_m_r(1):nr-s_m_r(end),:);
    fsmr(Mast_ind) = 1-s_m_r(1);
    for i = 1:length(Slave_ind)
%         if move_r(Slave_ind(i)) > 0
            sar_image_temp{i}(1:nr-s_m_r(end)+s_m_r(1),:) = sar_image{Slave_ind(i)}(1-s_m_r(1)+move_r(Slave_ind(i)):nr-s_m_r(end)+move_r(Slave_ind(i)),:);
            fsmr(Slave_ind(i)) = 1+move_r(Slave_ind(i))-s_m_r(1);
%         else
%             sar_image_temp{i}(1:nr-s_m_r(end)+s_m_r(1),:) = sar_image{Slave_ind(i)}(1+move_r(Slave_ind(i)):nr+move_r(Slave_ind(i))-s_m_r(end)+s_m_r(1),:);
%             fsmr(Slave_ind(i)) = 1+move_r(Slave_ind(i));
%         end
    end
end
if s_m_c(1) >= 0 && s_m_c(end) >= 0
    sar_image_cut{Mast_ind} = master_temp(:,1:nc-s_m_c(end));
    fsmc(Mast_ind) = 1;
    for i = 1:length(Slave_ind)
        sar_image_cut{Slave_ind(i)} = sar_image_temp{i}(:,1+move_c(Slave_ind(i)):nc+move_c(Slave_ind(i))-s_m_c(end));
        fsmc(Slave_ind(i)) = 1+move_c(Slave_ind(i));
    end
elseif s_m_c(1) < 0  && s_m_c(end) <= 0
    sar_image_cut{Mast_ind} = master_temp(:,1-s_m_c(1):nc);
    fsmc(Mast_ind) = 1-s_m_c(1);
    for i = 1:length(Slave_ind)
        sar_image_cut{Slave_ind(i)} = sar_image_temp{i}(:,1+move_c(Slave_ind(i))-s_m_c(1):nc+move_c(Slave_ind(i)));
        fsmc(Slave_ind(i)) = 1+move_c(Slave_ind(i))-s_m_c(1);
    end
elseif s_m_c(1) < 0 && s_m_c(end) > 0
    sar_image_cut{Mast_ind} = master_temp(:,1-s_m_c(1):nc-s_m_c(end));
    fsmc(Mast_ind) = 1-s_m_c(1);
    for i = 1:length(Slave_ind)
%         if move_c(Slave_ind(i)) > 0
            sar_image_cut{Slave_ind(i)} = sar_image_temp{i}(:,1+move_c(Slave_ind(i))-s_m_c(1):nc-s_m_c(end)+move_c(Slave_ind(i)));
            fsmc(Slave_ind(i)) = 1+move_c(Slave_ind(i))-s_m_c(1);
%         else
%             sar_image_cut{Slave_ind(i)} = sar_image_temp{i}(:,1+move_c(Slave_ind(i)):nc-s_m_c(end)+move_c(Slave_ind(i))+s_m_c(1));
%             fsmc(i) = 1+move_c(Slave_ind(i));
%         end
    end   
end
end