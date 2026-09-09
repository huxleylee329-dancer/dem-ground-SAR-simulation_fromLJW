% orbit_mast = csvread(path);
% orbit_mast(:,1) = orbit_mast(:,1)*(1/2/Prf);
% orbit_mast(:,2:7) = orbit_mast(:,2:7)*1000;
% save('orbit_mast.mat', 'orbit_mast');

index = [432 975 2012 2948 4115 900 822 698 629 1573 1665 1733];
% index = [572 970 1423 1688 1880 2068 2340];
GCPs_mast1 = conp1(index,:);
GCPs_mast1(:,1:2) = GCPs_mast1(:,1:2) - startpix(1) + 1;
[GCPs_mast1(:,4),GCPs_mast1(:,5),GCPs_mast1(:,6)] = xyz2ell(GCPs_mast1(:,4:6));
GCPs_mast = zeros(size(index, 2),6);
GCPs_mast(:,1) = GCPs_mast1(:,4);
GCPs_mast(:,2) = GCPs_mast1(:,5);
GCPs_mast(:,3) = mean(GCPs_mast1(:,6));
GCPs_mast(:,4) = GCPs_mast1(:,1);
GCPs_mast(:,5) = GCPs_mast1(:,2);

tmp = zeros(size(index, 2), 1, 3);
tmp(:,:,1) = GCPs_mast1(:,5)/180*pi;%longtitude
tmp(:, :, 2) = GCPs_mast1(:,4)/180*pi;% latitude
tmp(:, :, 3) = GCPs_mast(:,3); % height
[xyz] = Ell2xyz(tmp);
x = zeros(size(index, 2), 1);
y = zeros(size(index, 2), 1);
z = zeros(size(index, 2), 1);

x(:) = xyz(:,:,1);
y(:) = xyz(:,:,2);
z(:) = xyz(:,:,3);
for ii = 1:size(index, 2)
    GCPs_mast(ii, 6) = sqrt(sum((x(ii) - MainPosX(GCPs_mast1(ii,1), GCPs_mast1(ii,2)))^2+(y(ii) - MainPosY(GCPs_mast1(ii,1), GCPs_mast1(ii,2)))^2+(z(ii) - MainPosZ(GCPs_mast1(ii,1), GCPs_mast1(ii,2)))^2));
end


GCPs_slave1 = conp2(index,:);
GCPs_slave1(:,1:2) = GCPs_slave1(:,1:2) - startpix(1) + 1;
[GCPs_slave1(:,4),GCPs_slave1(:,5),GCPs_slave1(:,6)] = xyz2ell(GCPs_slave1(:,4:6));
GCPs_slave = zeros(size(index, 2),6);
GCPs_slave(:,1) = GCPs_slave1(:,4);
GCPs_slave(:,2) = GCPs_slave1(:,5);
GCPs_slave(:,3) = mean(GCPs_slave1(:,6));
GCPs_slave(:,4) = GCPs_slave1(:,1);
GCPs_slave(:,5) = GCPs_slave1(:,2);

tmp = zeros(size(index, 2), 1, 3);
tmp(:,:,1) = GCPs_slave1(:,5)/180*pi;%longtitude
tmp(:, :, 2) = GCPs_slave1(:,4)/180*pi;% latitude
tmp(:, :, 3) = GCPs_slave(:,3); % height
[xyz] = Ell2xyz(tmp);
x = zeros(size(index, 2), 1);
y = zeros(size(index, 2), 1);
z = zeros(size(index, 2), 1);

x(:) = xyz(:,:,1);
y(:) = xyz(:,:,2);
z(:) = xyz(:,:,3);
for ii = 1:size(index, 2)
    GCPs_slave(ii, 6) = sqrt(sum((x(ii) - SlavePosX(GCPs_slave1(ii,1), GCPs_slave1(ii,2)))^2+(y(ii) - SlavePosY(GCPs_slave1(ii,1), GCPs_slave1(ii,2)))^2+(z(ii) - SlavePosZ(GCPs_slave1(ii,1), GCPs_slave1(ii,2)))^2));
end

%% 
save('GCPs_mast.mat', 'GCPs_mast');
save('GCPs_slave.mat', 'GCPs_slave');
mast = MainSLC1;
slave = SlaveSLC1;
save('mast.mat', 'mast');
save('slave.mat', 'slave');
save('conp1.mat', 'conp1');
save('conp2.mat', 'conp2');
Slantrange_s = Slantrange2;
Slantrange_m = Slantrange;

% save('SlavePosX.mat', 'SlavePosX');
% save('SlavePosY.mat', 'SlavePosY');
% save('SlavePosZ.mat', 'SlavePosZ');
save('Slantrange_s.mat', 'Slantrange_s');
save('Slantrange_m.mat', 'Slantrange_m');
% 
% %% Œ¿–«ÀŸ∂»
% tmp = csvread(path);
% V_slave_x = zeros(size(slave)) * mean(tmp(:,end - 2));
% V_slave_y = zeros(size(slave)) * mean(tmp(:,end - 1));
% V_slave_z = zeros(size(slave)) * mean(tmp(:,end - 0));
% save('V_slave_x.mat', 'V_slave_x');
% save('V_slave_y.mat', 'V_slave_y');
% save('V_slave_z.mat', 'V_slave_z');


%% 
R_M = Slantrange;
pos = csvread(path);
Vs = mean(pos(:, 5:7))*1000;
Satellite_M_R_Velocity = repmat(Vs, size(R_M, 1), 1);
Satellite_M_T_Position = zeros(size(R_M, 1), 3);
Satellite_M_T_Position(:, 1) = mean(MainPosX, 2);
Satellite_M_T_Position(:, 2) = mean(MainPosY, 2);
Satellite_M_T_Position(:, 3) = mean(MainPosZ, 2);

Satellite_S_T_Position = zeros(size(R_M, 1), 3);
Satellite_S_T_Position(:, 1) = mean(SlavePosX, 2);
Satellite_S_T_Position(:, 2) = mean(SlavePosY, 2);
Satellite_S_T_Position(:, 3) = mean(SlavePosZ, 2);

Rmin = Slantrange(2,2);
c = 3e8;

Satellite_M_R_Position = Satellite_M_T_Position+2*(Rmin)/c*Satellite_M_R_Velocity;
Satellite_S_R_Position = Satellite_S_T_Position+2*(Rmin)/c*Satellite_M_R_Velocity;

Satellite_M=(Satellite_M_R_Position+Satellite_M_T_Position)./2;

save('R_M.mat','R_M');
save('Satellite_M_R_Velocity.mat','Satellite_M_R_Velocity');
save('Satellite_M_T_Position.mat','Satellite_M_T_Position');
save('Satellite_S_T_Position.mat','Satellite_S_T_Position');
save('Rmin.mat','Rmin');
save('Satellite_M_R_Position.mat','Satellite_M_R_Position');
save('Satellite_S_R_Position.mat','Satellite_S_R_Position');
save('Satellite_M.mat','Satellite_M');