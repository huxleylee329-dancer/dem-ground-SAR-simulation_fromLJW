%---------------------------去除参考平面的相位（低频）--------------------------
para_path = '.\Kernel\parameter\low\';% 控制点参数路径
out_path = '.\Kernel\Deflat_output\low\';% 去平地结果保存路径 
SLC_path = '.\Kernel\Input\low\';% 主辅图像路径
warning off;
% 检查并读入参数
display('检查并读入参数......');
if ~exist([para_path 'auxi_mast.txt'], 'file')
    display('主星雷达工作参数文件不存在！请检查......');
    pause(5);
    quit;
else
    [frequency_mast, ~, nr1, nc1] = textread([para_path 'auxi_mast.txt'], '%f %f %n %n');
end

if ~exist([para_path 'auxi_slave.txt'], 'file')
    display('辅星雷达工作参数文件不存在！请检查......');
    pause(5);
    quit;
else
    [frequency_slave, prf, nr2, nc2] = textread([para_path 'auxi_slave.txt'], '%f %f %n %n');
end

if ~exist([para_path 'GCPs_mast.mat'], 'file')
    display('主星控制点参数文件不存在！请检查......');
    pause(5);
    quit;
else
    load([para_path 'GCPs_mast.mat']);
end


if ~exist([para_path 'GCPs_slave.mat'], 'file')
    display('辅星控制点参数文件不存在！请检查......');
    pause(5);
    quit;
else
    load([para_path 'GCPs_slave.mat']);
end
display('完成！');

%% 主图像去平地
display('主图去平地......');
% 线性拟合
N_gcps = size(GCPs_mast, 1);
rows = GCPs_mast(:,4); cols =  GCPs_mast(:,5); 
matrix = [ones(N_gcps,1), rows, cols]; 
a_pksi_qeta = (matrix'*matrix)\matrix'*GCPs_mast(:,6);
a_pksi_qeta(2) = 0;
% err = norm(matrix*a_pksi_qeta - Rps_m);

flat_range_master = zeros(nr1, nc1);
for i = 1:nr1
    for j = 1:nc1
        flat_range_master(i,j) = [1, i, j]*a_pksi_qeta; %斜距
    end
end


% 载入主图像（单视复图像）
if ~exist([SLC_path 'mast.mat'], 'file')
    display('主星单视复图像文件不存在！请检查......');
    pause(5);
    quit;
else
    load([SLC_path 'mast.mat']);
end
lambda1 = 3e8/frequency_mast;
flat_phase_master = flat_range_master * 4 * pi / lambda1;  %主图平地相位
% 去平地效应
deflat_mast_low = mast.*exp(1i*flat_phase_master);
% 保存结果
display('保存结果......');
flat_phase_mast_low = flat_phase_master;
save([out_path 'flat_phase_mast_low.mat'], 'flat_phase_mast_low');% 保存绝对平地相位（后续加回来）
save([out_path 'deflat_mast_low.mat'], 'deflat_mast_low');% 保存去平地之后的主图像
fig1 = figure;
set(fig1, 'visible', 'off');
imagesc(angle(exp(1i*flat_phase_master)));colorbar;
saveas(fig1, [out_path 'flat_phase_low.jpg']);% 保存平地相位图片（缠绕后的）
display('完成！');




%% 辅图像去平地
display('辅图去平地......');
rows = GCPs_slave(:,4); cols =  GCPs_slave(:,5); 
matrix = [ones(N_gcps,1), rows, cols]; % 线性拟合
a_pksi_qeta = (matrix'*matrix)\matrix'*GCPs_slave(:,6);
a_pksi_qeta(2) = 0;
err = norm(matrix*a_pksi_qeta - GCPs_slave(:,6));
 %根据9个控制点，利用最小二乘法拟合出方程系数，然后利用该系数求出全部像素点的平地干涉相位；      
%     a_pksi_qeta = matrix(4:end,:)\Gsflat_phi(4:end,k)
flat_range_slave = zeros(nr1, nc1);
for i = 1:nr1
    for j = 1:nc1
        %flat_int_phase(ii,jj,k) = [1, ii, jj, ii*jj, ii*ii, jj*jj]*a_pksi_qeta;
        flat_range_slave(i,j) = [1, i, j]*a_pksi_qeta; %斜距
    end
end

% 载入辅图像（单视复图像）
if ~exist([SLC_path 'slave.mat'], 'file')
    display('辅星单视复图像文件不存在！请检查......');
    pause(5);
    quit;
else
    load([SLC_path 'slave.mat']);
end
lambda2 = 3e8/frequency_slave;
flat_phase_slave_low = flat_range_slave * 4 * pi / lambda2;  %主图平地相位
%去平地效应
deflat_slave_low = slave.*exp(1i*flat_phase_slave_low);
display('保存结果......');
save([out_path 'deflat_slave_low.mat'], 'deflat_slave_low');% 保存去平地之后的辅图像
save([out_path 'flat_phase_slave_low.mat'], 'flat_phase_slave_low');% 保存绝对平地相位
display('完成！');
display('去平地完成！');














