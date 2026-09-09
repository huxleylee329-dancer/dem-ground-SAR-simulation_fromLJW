function [SLC_image,Slantrange]= SLC(Satellite,Targetc,Random_phase,DEM,Dx,Dy,fc,Drange,Dazimuth,Azimuth_cell,Range_cell)
c = 3e8; 
[Nta,Ntr] = size(DEM);
Targetc.slant_range = sqrt((Satellite.range-Targetc.range)^2 + ...
                      (Satellite.azimuth-Targetc.azimuth)^2 + ...
                      (Satellite.height-Targetc.height)^2);              
%% DEM radar code (x,y,z)==>(x,r,s)
r       = zeros(Nta,Ntr);
% s       = zeros(Nta,Ntr);
x       = zeros(Nta,Ntr);
theta_T = zeros(Nta,Ntr);
iic     = round(Nta/2);
jjc     = round(Ntr/2);
Target.range   = 0;
Target.azimuth = 0;
Target.height  = 0;
h = waitbar(0);
for ii = 1:1:Nta
    for jj = 1:1:Ntr
        % target position (x,y,z)
        Target(ii,jj).azimuth = (iic - ii) * Dx;                 %x-axis
        Target(ii,jj).range   = Targetc.range - (jjc - jj) * Dy; %y-axis
        Target(ii,jj).height  = DEM(ii,jj);                      %z-axis
        % 正侧视计算最短斜距（零多普勒斜距）     
        Target(ii,jj).slant_range = sqrt((Satellite.range-Target(ii,jj).range)^2 + (Satellite.height-Target(ii,jj).height)^2); 
        % 入射角
%         theta_T(ii,jj) = acos((Satellite.height-Target(ii,jj).height)/Target(ii,jj).slant_range);
%         Target(ii,jj).r = Target(ii,jj).slant_range * cos(Satellite.lookangle - theta_T(ii,jj)); % r-axis
%         Target(ii,jj).s = Target(ii,jj).slant_range * sin(Satellite.lookangle - theta_T(ii,jj)); % s-axis
        r(ii,jj) = Target(ii,jj).slant_range;
%         s(ii,jj) = Target(ii,jj).s;
        x(ii,jj) = Target(ii,jj).azimuth;
    end 
    waitbar((ii - 1)*Ntr/(Ntr*Nta));
end
close all;
%% slope
[slope_a,slope_r] = terrain_slope(DEM,Dx,Dy);
%% backscattering coefficient
lambda = c/fc;
sigma = terrain_sigma(Satellite.lookangle,slope_a,slope_r,Dx,Dy);
%% 定标点
[r_p,c_p,~] = find(min(DEM(:)));
if slope_r(r_p,c_p) >= 0
    sigma(r_p,c_p) = 1e3;
end
% sigma(iic,3*jjc/4) = 1e3;%-----------------20150410
%% projection to SAR image  
% Azimuth_cell = 0.13;
% Range_cell   = 0.11;
x0 = min(x(:));
r0 = min(r(:));
xl = max(x(:));
rl = max(r(:));
M  = floor((xl-x0)/Azimuth_cell)+1;
N  = floor((rl-r0)/Range_cell)+1;
rho_MN = zeros(M,N);
Slantrange = zeros(M,N);
count = zeros(M,N);
h = waitbar(0);
for ii = 1:1:Nta
    for jj = 1:1:Ntr
        % image coordination      
        m = floor((x(ii,jj) - x0)/Azimuth_cell) + 1;
        n = floor((r(ii,jj) - r0)/Range_cell)  + 1;
        % rho image 
        rho_IJ = sigma(ii,jj) *exp(-j*4*pi*r(ii,jj)/lambda)*exp(j*pi*Random_phase(ii,jj));
        rho_MN(m,n)= rho_MN(m,n) + rho_IJ;
        Slantrange(m,n)= Slantrange(m,n) + r(ii,jj);
        count(m,n)= count(m,n) + 1;
    end
    waitbar((ii - 1)*Ntr/(Ntr*Nta));
end
count = count + 1e-7;
Slantrange = Slantrange./count;
close all;
%% SAR imaging transfer function --H--
% range
% Drange = 3e8/9e8;
% Dazimuth = Drange;
t1  = round(N/2);
sequence_r = -t1:1:(N-1-t1);
h_r = sinc(pi/Drange * sequence_r);
H_r =  ones(M,1) * fft(h_r);
t2  = round((M)/2);
sequence_a = -t2:1:(M-1-t2);
h_a = sinc(pi/Dazimuth * sequence_a);
H_a = fft(h_a)'* ones(1,N);
H   = H_r .* H_a;
%% window function --W--
% wr  = hann(N)';
% W_r = ones(M,1) * fft(wr);
% wa  = hann(M)';
% W_a = fft(wa)'* ones(1,N);
% W   = W_r .* W_a;
%% SAR image
I1        = fft2(rho_MN);
I2        = I1 .* H ;
I3        = fftshift(ifft2(I2)) ;
%% SAR image cut
Ncut_start = floor((r(:,1) - r0)/Range_cell)  + 1;
Ncut_end   = floor((r(:,end) - r0)/Range_cell)  + 1;
SLC_image  = I3(:,max(Ncut_start):min(Ncut_end));
end
%% 输出
% sigma
% figure;
% imagesc(abs(sigma))
% axis xy
% axis equal
% axis tight 
% colormap(gray)
% title('scatter coefficients (abs)','fontWeight','Bold')
% % DEM
% figure;
% imagesc(DEM)
% axis xy
% axis equal
% axis tight
% title('DEM','fontWeight','Bold')
% my_ruler = colorbar;
% set(get(my_ruler,'title'),'string','Height','fontWeight','Bold')
% scale = get(my_ruler,'YTickLabel');
% scale = strcat(scale,'m');
% set(my_ruler,'YTickLabel',scale);
% % SLC
% G=20*log10(abs(SLC_image)+1e-6);
% gm=max(max(G));
% gn=gm-30;%显示动态范围40dB
% G=255/(gm-gn)*(G-gn).*(G>gn);
% figure;imagesc(G)
% axis equal
% axis tight
% colormap(gray)
% title('Satellite SLC','fontWeight','Bold')
% xlabel('Range')
% ylabel('Azimuth')