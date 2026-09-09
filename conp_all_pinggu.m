%%%%控制点矫正绝对相位

clc
cnp1=csvread('..\Output\main_range_cut.csv');
cnp2=csvread('..\Output\slave1_range_cut.csv');
cnp3=csvread('..\Output\slave2_range_cut.csv');

%%减去起始行列（1723,1723）


del_r_long=cnp3-cnp1;
 lamda=0.2306;
 phi_cnp=4*pi*del_r_long./lamda;
 
 %%载入解缠相位
% 
% load('..\pinggu\ronghe2.mat');
% load('..\pinggu\unwrapped_phase1.mat');
% 
% 
% 
%    del_long=unwrapped_phase1-phi_conp;
%     del_ronghe =ronghe2-phi_conp;
% 
% meanerror_long=mean(mean(del_long));
% meanerror_ronghe=mean(mean(del_ronghe));

conp1=csvread('..\Output\conp1.csv');
conp2=csvread('..\Output\conp2.csv');
conp3=csvread('..\Output\conp3.csv');

%%减去起始行列（1723,1723）
del=ones(560,2)*1722;
conp1(:,1:2)=conp1(:,1:2)-del;
conp2(:,1:2)=conp2(:,1:2)-del;
conp3(:,1:2)=conp3(:,1:2)-del;

% del_r_short=conp2(:,3)-conp1(:,3);
 del_r_long=conp3(:,3)-conp1(:,3);
 lamda=0.2306;
 phi_conp=4*pi*del_r_long./lamda;
 
 [k,~]=size(conp1);
phase_abs = zeros(1,k);
for iii =1:k
    phase_abs(iii) =phi_conp(iii,1)-phi_cnp(conp1(iii,1),conp1(iii,2));
end
k_point = round(phase_abs/2/pi);



