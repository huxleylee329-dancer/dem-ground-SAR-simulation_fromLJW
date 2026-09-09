%% 低频相位解缠
warning off;
out_folder = '.\Kernel\Output3_unwrap\';
in_folder = '.\Kernel\Output2_filter\';
%% 读取低频缠绕相位数据
if exist([in_folder 'wrapped_phase_filtered_low.mat'], 'file')
    display('载入缠绕相位数据......');
    load([in_folder 'wrapped_phase_filtered_low.mat']);
    display('完成！');
    display('进行相位解缠......');
    unwrapped_phase_low = MCF(wrapped_phase_filtered);
    display('完成！');
    
    display('保存结果......');
    fig1 = figure;
    set(fig1, 'visible', 'off');
    imagesc(unwrapped_phase_low);
    colorbar;
    saveas(fig1, [out_folder 'unwrapped_phase_low.jpg']);
    save([out_folder 'unwrapped_phase_low.mat'], 'unwrapped_phase_low');
    display('完成！');
else
    display('低频缠绕相位数据不存在！请检查......');
end


%% 校正至绝对相位
in_folder1 = '.\Kernel\Deflat_output\low\';
para_path = '.\Kernel\parameter\low\';
%% 读入参数
load([para_path 'Slantrange_s.mat']);
load([para_path 'Slantrange_m.mat']);
load([para_path 'conp1.mat']);
load([para_path 'conp2.mat']);

load([in_folder1 'flat_phase_slave_low.mat']);
load([in_folder1 'flat_phase_mast_low.mat']);


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

lambda1 = 3e8/frequency_mast;
lambda2 = 3e8/frequency_slave;

%% 加回平地相位
unwrapped_phase_low_add = unwrapped_phase_low + (flat_phase_slave_low - flat_phase_mast_low);
% unwrapped_phase_low_add = - unwrapped_phase_low_add;
%% 校正至绝对相位
[k,~] = size(conp1);
del = ones(k, 2)*3965;

conp1(:,1:2) = conp1(:,1:2) - del;
conp2(:,1:2) = conp2(:,1:2) - del;
count = 0;
for iii = 1:k
    if (conp1(iii, 1) > 0 && conp1(iii, 1) < 261) && (conp1(iii, 2) > 0 && conp1(iii, 2) < 261)
        count = count + 1;
    end
end
index = zeros(count, 1);
jj = 0;
for iii = 1:k
    if (conp1(iii, 1) > 0 && conp1(iii, 1) < 261) && (conp1(iii, 2) > 0 && conp1(iii, 2) < 261)
        jj = jj + 1;
        index(jj) = iii;
    end
end
conp1_bk = conp1(index,:);
conp2_bk = conp2(index,:);
phi_1 = conp1_bk(:,3) / lambda1 * 4 * pi;
phi_2 = conp2_bk(:,3) / lambda2 * 4 * pi;
phi_abs = abs(phi_1 - phi_2);
phase_abs = zeros(1,count);
for iii = 1:count
    phase_abs(iii) = phi_abs(iii) - unwrapped_phase_low_add(conp1_bk(iii, 1), conp1_bk(iii, 2));
end
k_point = round(phase_abs/(2*pi));
k_value = unique(k_point);
[~, k_value_number] = size(k_value);
n = 1;
for ii = 1:k_value_number
    m = sum(k_point==k_value(1,ii));
    if m > n
        n = m;
        k_real = k_value(ii);
    end
end

unwrapped_phase_low_abs = unwrapped_phase_low_add + 2 * pi * k_real;


%% 
phase_real =  (Slantrange_s/lambda2*4*pi - Slantrange_m/lambda1*4*pi);
unwrapped_phase_low_abs = unwrapped_phase_low_abs + mean(mean(phase_real - unwrapped_phase_low_abs));

delta_low = unwrapped_phase_low_abs - phase_real;
unwrapped_phase_low_for_fusion = unwrapped_phase_low_abs - (flat_phase_slave_low - flat_phase_mast_low);


save([out_folder 'unwrapped_phase_low_for_fusion.mat'], 'unwrapped_phase_low_for_fusion');%保存绝对相位
save([out_folder 'unwrapped_phase_low_abs.mat'], 'unwrapped_phase_low_abs');%保存绝对相位
save([out_folder 'delta_low.mat'], 'delta_low');%保存相位误差

fig2 = figure;
set(fig2, 'visible', 'off');
imagesc(unwrapped_phase_low_abs);
colorbar;
saveas(fig2, [out_folder 'unwrapped_phase_low_abs.jpg']);