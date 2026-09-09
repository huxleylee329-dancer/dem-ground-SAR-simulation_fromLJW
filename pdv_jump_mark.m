function [dz] = pdv_jump_mark(PDV,image_phase)
%利用相位导数方差和跳变点定义标记不连续点

%利用相位导数方差和图像膨胀标记不连续点
pdv_threshold = 0.4;%相位导数方差阈值
[pr,pc] = size(PDV);
PDV(PDV>pdv_threshold)=1;%将大于阈值的点标为1
dz_pdv = zeros(pr,pc);
dz_pdv(PDV==1) = 1;%dy为全1矩阵
se = strel('disk',3);%垂直线结构元素:圆形
dz_pdv = imdilate(dz_pdv,se);

%利用跳变点定义和图像膨胀标记不连续点
phi_threshold = 1*pi;
%跳变点是指那些相位梯度在【-π，π】之间的不连续点，因此在判断时如果绝对值大于π即为跳变点，而不是2π
[pr,pc] = size(image_phase);
image_left = image_phase(:,1:pc-1);
image_right = image_phase(:,2:pc);
dx = floor(abs(image_left-image_right)/phi_threshold);
image_up = image_phase(1:pr-1,:);
image_down = image_phase(2:pr,:);
dy = floor(abs(image_up-image_down)/phi_threshold);
dx(dx>1)=1;
dy(dy>1)=1;
dz1=dx(2:end,:)+dy(:,1:end-1);
dz2=dx(2:end,:)+dy(:,2:end);
dz3=dx(1:end-1,:)+dy(:,1:end-1);
dz4=dx(1:end-1,:)+dy(:,2:end);
dz_jump = dz1 + dz2 + dz3 + dz4;
dz_jump(dz_jump>1)=1;
se = strel('disk',5);%垂直线结构元素:圆形
% se = strel('diamond',10);%垂直线结构元素：菱形
% se = strel('square',8);%垂直线结构元素：方形
dz_jump = imdilate(dz_jump,se);

%取dz_pdv和dz_jump的并集
dz = dz_jump + dz_pdv(1:end-1,1:end-1);
dz(dz>1)=1;


[pr,pc,~] = find(dz == 1); 
% Gradiennt_jump = sum(sum(dz));%跳变点个数


figure;
imagesc(image_phase);
hold on
scatter(pc,pr,'r.');%跳变点

end