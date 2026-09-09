%-----------------------------------------------------------%

%-----------------------------------------------------------%
% close all
% clear all
% clc
j = sqrt(-1);                       
pi = asin(1)*2; 
tic
%% input
%------------------- DEM data --------------------------%
data = load('Input_Data\DEM.mat');
% DEM = data.DEM ;
DEM = zeros(size(data.DEM)) ;
DEM = repmat(DEM(1:256,:), 1, 11);
[Nta,Ntr] = size(DEM);
%----------------- terrain parameters ------------------%
Azimuth_cell = 0.125;
Range_cell   = 0.125;
Dx = Azimuth_cell/4;                          %DEM的方位网格,m            
Dy = Range_cell/4;                            %DEM的距离网格,m
Drange = 0.333;
Dazimuth = Drange;
%----------------- baseline parameters -----------------%
fc = 35e9;                        %载频
Baseline =1;                     %基线长度,m
alfa = deg2rad(5);                 %基线倾角,rad
gamma = deg2rad(0);                 %基线方位角,rad 双星正侧视为0
%---------------- master satellite position ------------%
Satellite1.range = 0;               %卫星距离向坐标,m
Satellite1.azimuth = 0;             %卫星方位向坐标,m
Satellite1.height = 2885;          %卫星高度,m    
Satellite1.lookangle = deg2rad(45); %卫星下视角,rad
%-------------- center target position -----------------%
Targetc.height = DEM(round(Nta/2),round(Ntr/2));                                 %高度,m
Targetc.range = (Satellite1.height - Targetc.height) * tan(Satellite1.lookangle);%距离向坐标,m
Targetc.azimuth = 0;                                                             %方位向坐标,m
%---------------- slave satellite ----------------------%
Satellite2.range = Baseline * cos(alfa) * cos(gamma);            %星2距离向坐标,m
Satellite2.azimuth = Baseline * cos(alfa) * sin(gamma) ;         %星2方位向坐标,m
Satellite2.height = Satellite1.height + Baseline*sin(alfa);      %星2高度,m
Satellite2.lookangle = atan((Targetc.range-Baseline * cos(alfa))/Satellite2.height);
%---------------- random phase noise -------------------%
Random_phase = -1 + 2*rand(Nta,Ntr);
%% SLC
[SLC1_temp,Slantrange_Main]= SLC(Satellite1,Targetc,Random_phase,DEM,Dx,Dy,fc,Drange,Dazimuth,Azimuth_cell,Range_cell);
[SLC2_temp,Slantrange_ideal]= SLC(Satellite2,Targetc,Random_phase,DEM,Dx,Dy,fc,Drange,Dazimuth,Azimuth_cell,Range_cell);
%% cut
[Naslc1,Nrslc1] = size(SLC1_temp); 
[Naslc2,Nrslc2] = size(SLC2_temp); 
SLC1 = SLC1_temp(1:min(Naslc1,Naslc2),1:min(Nrslc1,Nrslc2));
SLC2_ideal = SLC2_temp(1:min(Naslc1,Naslc2),1:min(Nrslc1,Nrslc2));
phase_ideal = angle(SLC1.*conj(SLC2_ideal));
figure; imagesc(phase_ideal)
%% save
% [s,~,~] = mkdir('.\Output');
% save('.\Output\SLC1.mat','SLC1','-mat');   
% save('.\Output\SLC2.mat','SLC2','-mat');   