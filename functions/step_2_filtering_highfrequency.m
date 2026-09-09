%% 高频滤波（斜坡自适应法）
%% 滤波窗口大小设置
warning off;
display('读入滤波参数（滤波窗口大小）......');
para_path = '.\Kernel\window_size.txt';
if ~exist(para_path,'file')
    win_size = 23;% 默认窗口大小
else
    win_size = textread(para_path, '%n');
end
display('读入完成！');
%% 载入含噪相位
display('载入含噪相位......');
out_folder = '.\Kernel\Output2_filter\';
phase_path = '.\Kernel\Output1_registration\wrapped_phase_high.mat';
if exist(phase_path,'file')
    load(phase_path);
    display('载入完成！');
    display('进行滤波（计算时间较长）......');
    wrapped_phase_filtered = Improved_Slope_Adaptive_filter_parallel(wrapped_phase,win_size,win_size); %坡度自适应滤波
    display('滤波完成！');
    display('保存结果......');
    fig1 = figure;
    set(fig1, 'visible', 'off');
    imagesc(wrapped_phase_filtered);
    colorbar;
    save([out_folder 'wrapped_phase_filtered_high.mat'], 'wrapped_phase_filtered');
    saveas(fig1, [out_folder 'wrapped_phase_filtered_high.jpg']);
    display('完成！');
else
    display('不存在含噪相位文件！请检查...');
end
