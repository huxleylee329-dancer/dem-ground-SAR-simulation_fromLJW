function [nres,residuemat] = Calculation_Residues(interf_phase,display)
% ---------------------------------------
% improved 20120819
% Calculation_Residues  计算干涉相位残差点.
% Input:  interf_phase: 干涉相位
% Output: nres:         残差点个数
%         residuemat:   残差点矩阵
%         residues {-2pi / 0 / +2pi}
% Example:
% load('X:\程序\InSAR\InSAR_数据\ERS1_PH.mat')
% interf_phase = ERS1_PH;
% ---------------------------------------
% 计算干涉相位梯度
[nrow ncol] = size(interf_phase);
diffrows = wrap(diff(interf_phase,1,1),2*pi);
diffcols = wrap(diff(interf_phase,1,2),2*pi);
residuemat = diff(diffcols,1,1) - diff(diffrows,1,2);
% 零残差点索引
nnon  = find(abs(residuemat)<1);           %零残差点,residues = 0
% 残差点个数
nres = (nrow-1)*(ncol-1)-length(nnon);
% 边缘点置零，归零非残差点
residuemat(nnon)   = 0;
residuemat(nrow,:) = 0;
residuemat(:,ncol) = 0;
% 正负残差点索引
npres = find(abs(residuemat-2*pi)<= 1);     %正残差点,residues = 2pi
nnres = find(abs(residuemat+2*pi)<= 1);     %负残差点,residues = 2pi
residuemat(npres)  = 2*pi;
residuemat(nnres)  = -2*pi;
[nr,nc,~] = find(residuemat == -2*pi);     %负残差点,residues = -2pi
[pr,pc,~] = find(residuemat == 2*pi);      %正残差点,residues = 2pi
% 输出
% figure;
% imagesc(residuemat)
% colormap(gray)
% hold on
% plot(pc,pr,'c.')
% hold on
% plot(nc,nr,'m.')
% title('residues')
if display == 1
figure;
imagesc(interf_phase)
axis equal
axis tight
colormap(gray)
hold on
plot(pc,pr,'c.')%正残差
hold on
plot(nc,nr,'m.')%负残差
title('wrapped phase and residues','fontWeight','Bold')
end