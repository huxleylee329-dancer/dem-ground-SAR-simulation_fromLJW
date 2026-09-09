
lamda1 = 3e8/(9.6e9);
lamda2 = 3e8/(15e9 - 18e6);
lamda2 = lamda1;
Real_phi=ronghe(2:end - 1,2:end - 1);
[Cut_Na,Cut_Nr]=size(Real_phi);
fd=ones(1,Cut_Nr)*(-0.278713789746405)*0; %多普勒中心频率
%% 载入数据
load('R_M.mat');
load('Satellite_M_R_Velocity.mat');
load('Satellite_M_T_Position.mat');
load('Satellite_S_T_Position.mat');
load('Rmin.mat');
load('Satellite_M_R_Position.mat');
load('Satellite_S_R_Position.mat');
load('Satellite_M.mat');

R_M = R_M(2:end - 1,2:end - 1);
Satellite_M_R_Velocity = Satellite_M_R_Velocity(2:end - 1,:);
Satellite_M_T_Position = Satellite_M_T_Position(2:end - 1,:);
Satellite_S_T_Position = Satellite_S_T_Position(2:end - 1,:);
Satellite_M_R_Position = Satellite_M_R_Position(2:end - 1,:);
Satellite_S_R_Position = Satellite_S_R_Position(2:end - 1,:);
Satellite_M = Satellite_M(2:end - 1,:);



% 控制点
P = [6378932.39735364; -1007.02513948654; 1014.70937408723];
%% 高程反演
m=30;%迭代50次

Vs=Satellite_M_R_Velocity/1000;

R_F = 2 * lamda2 * R_M / lamda1 - Real_phi * lamda2 / (2 * pi);



Df=zeros(3,3);
Df_ni=zeros(3,3);
f=zeros(3,1);
delta_Rt=zeros(3,1); 
P1=ones(Cut_Na,Cut_Nr)*P(1);
P2=ones(Cut_Na,Cut_Nr)*P(2);
P3=ones(Cut_Na,Cut_Nr)*P(3);
clear Real_phi B_Vetor Satellite_M_R_Velocity x y z unwrapped_phase ii jj RanStartTime

% tic
for k=0:m-1
%     tic
%     k
%--------------------------------------------------    
% f1=(Satellite_M(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_M(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_M(:,3)*ones(1,Cut_Nr)- P3).^2- R_M.^2;
M_T=(Satellite_M_T_Position(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_M_T_Position(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_M_T_Position(:,3)*ones(1,Cut_Nr)- P3).^2;
M_R=(Satellite_M_R_Position(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_M_R_Position(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_M_R_Position(:,3)*ones(1,Cut_Nr)- P3).^2;
f1=sqrt(M_T)+sqrt(M_R)- 2*R_M;
%--------------------------------------------------




S_T=(Satellite_S_T_Position(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_S_T_Position(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_S_T_Position(:,3)*ones(1,Cut_Nr)- P3).^2;
S_R=(Satellite_S_R_Position(:,1)*ones(1,Cut_Nr)-P1).^2+(Satellite_S_R_Position(:,2)*ones(1,Cut_Nr)- P2).^2+(Satellite_S_R_Position(:,3)*ones(1,Cut_Nr)- P3).^2;
f2=sqrt(S_T)+sqrt(S_R)- R_F;

f3=Vs(:,1)*ones(1,Cut_Nr).*(Satellite_M(:,1)*ones(1,Cut_Nr)-P1)+Vs(:,2)*ones(1,Cut_Nr).*(Satellite_M(:,2)*ones(1,Cut_Nr)-P2)+Vs(:,3)*ones(1,Cut_Nr).*(Satellite_M(:,3)*ones(1,Cut_Nr)-P3)-ones(Cut_Na,1)*fd.*R_M*lamda1/2.0;
% temp=sqrt(f1.^2+f2.^2+f3.^2);
%% Dff
%第一行：f(1)的x，y，z的导数
%--------------------------------------------------   

Df11=1/2*(M_T).^(-1/2)*(-2.0).*(Satellite_M_T_Position(:,1)*ones(1,Cut_Nr)-P1)...
+1/2*(M_R).^(-1/2)*(-2.0).*(Satellite_M_R_Position(:,1)*ones(1,Cut_Nr)-P1);

Df12=1/2*(M_T).^(-1/2)*(-2.0).*(Satellite_M_T_Position(:,2)*ones(1,Cut_Nr)-P2)...
+1/2*(M_R).^(-1/2)*(-2.0).*(Satellite_M_R_Position(:,2)*ones(1,Cut_Nr)-P2);

Df13=1/2*(M_T).^(-1/2)*(-2.0).*(Satellite_M_T_Position(:,3)*ones(1,Cut_Nr)-P3)...
+1/2*(M_R).^(-1/2)*(-2.0).*(Satellite_M_R_Position(:,3)*ones(1,Cut_Nr)-P3);
%-------------------------------------------------- 
clear M_T M_R
%第二行：f(2)的x，y，z的导数
Df21=1/2*(S_T).^(-1/2)*(-2.0).*(Satellite_S_T_Position(:,1)*ones(1,Cut_Nr)-P1)...
+1/2*(S_R).^(-1/2)*(-2.0).*(Satellite_S_R_Position(:,1)*ones(1,Cut_Nr)-P1);

Df22=1/2*(S_T).^(-1/2)*(-2.0).*(Satellite_S_T_Position(:,2)*ones(1,Cut_Nr)-P2)...
+1/2*(S_R).^(-1/2)*(-2.0).*(Satellite_S_R_Position(:,2)*ones(1,Cut_Nr)-P2);

Df23=1/2*(S_T).^(-1/2)*(-2.0).*(Satellite_S_T_Position(:,3)*ones(1,Cut_Nr)-P3)...
+1/2*(S_R).^(-1/2)*(-2.0).*(Satellite_S_R_Position(:,3)*ones(1,Cut_Nr)-P3);

clear S_T S_R
%第三行：f(3)的x，y，z的导数
Df31=-Vs(:,1)*ones(1,Cut_Nr);
Df32=-Vs(:,2)*ones(1,Cut_Nr);
Df33=-Vs(:,3)*ones(1,Cut_Nr);
%% 
det_Df=Df11.*Df22.*Df33+Df12.*Df23.*Df31+Df13.*Df21.*Df32-Df31.*Df22.*Df13-Df32.*Df23.*Df11-Df33.*Df21.*Df12;
% temp=temp+det_Df;%求Df行列式且与temp相加
Df_ni11=(Df22.*Df33-Df32.*Df23)./det_Df;
Df_ni12=-(Df12.*Df33-Df32.*Df13)./det_Df;
Df_ni13=(Df12.*Df23-Df22.*Df13)./det_Df;
delta_Rt1=Df_ni11.*f1+Df_ni12.*f2+Df_ni13.*f3;
clear  Df_ni11 Df_ni12 Df_ni13
Df_ni21=-(Df21.*Df33-Df31.*Df23)./det_Df;
Df_ni22=(Df11.*Df33-Df31.*Df13)./det_Df;
Df_ni23=-(Df11.*Df23-Df21.*Df13)./det_Df;
delta_Rt2=Df_ni21.*f1+Df_ni22.*f2+Df_ni23.*f3;
clear   Df_ni21 Df_ni22 Df_ni23
Df_ni31=(Df21.*Df32-Df31.*Df22)./det_Df;
Df_ni32=-(Df11.*Df32-Df31.*Df12)./det_Df;
Df_ni33=(Df22.*Df11-Df21.*Df12)./det_Df;
delta_Rt3=Df_ni31.*f1+Df_ni32.*f2+Df_ni33.*f3;
clear det_Df Df11 Df12 Df13 Df21 Df22 Df23 Df31 Df32 Df33 Df_ni31 Df_ni32 Df_ni33 f1 f2 f3




% temp=sqrt((delta_Rt1.^2+delta_Rt2.^2+delta_Rt3.^2)./(P1.^2+P2.^2+P3.^2));
P1=P1-delta_Rt1;
P2=P2-delta_Rt2;
P3=P3-delta_Rt3;
clear delta_Rt1 delta_Rt2 delta_Rt3

% toc
end
Target_84=zeros(Cut_Na,Cut_Nr,3);
Target_84(:,:,1)=P1;
Target_84(:,:,2)=P2;
Target_84(:,:,3)=P3;
clear P1 P2 P3
DEM_x=Target_84(:,:,1);
DEM_y=Target_84(:,:,2);
DEM_z=Target_84(:,:,3);
DEM=sqrt(DEM_x.^2+DEM_y.^2+DEM_z.^2);



