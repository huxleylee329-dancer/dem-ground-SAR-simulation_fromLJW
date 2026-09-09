%%%%控制点矫正绝对相位

clc
clear

conp1=csvread('..\Output\conp1.csv');
conp2=csvread('..\Output\conp2.csv');
conp3=csvread('..\Output\conp3.csv');

[k,~]=size(conp1);
%%减去起始行列（1723,1723）
del=ones(k,2)*1717;
conp1(:,1:2)=conp1(:,1:2)-del;
conp2(:,1:2)=conp2(:,1:2)-del;

del_r_short=conp2(:,3)-conp1(:,3);
lamda=0.02;
phi_short=4*pi*del_r_short./lamda;

 %%矫正相位到绝对相位

load('..\Mid\500_mast&500_slave_insar\unwrapped_phase.mat');

[k,~]=size(conp1);
phase_abs = zeros(1,  k);
for iii =1:k
    phase_abs(iii) =phi_short(iii,1)-unwrapped_phase(conp1(iii,1),conp1(iii,2));
end
k_point = round(phase_abs/2/pi);
% figure;plot(phase_abs/2/pi);
% title('每个控制点的粗k值');
k_value = unique(k_point);%列出算出的不同k的值
[~,k_value_number] = size(k_value);
n=1;
% k_real =1;%%%%%ls
for ii = 1:k_value_number
m = sum(k_point == k_value(ii));%将最多的k值作为真实k
if m >n
    n=m;
    k_real = k_value(ii);
end
end
unwrapped_phase1=unwrapped_phase+k_real*2*pi;
save('..\Mid\500_mast&500_slave_insar\unwrapped_phase1.mat','unwrapped_phase1');

