
air_height = 3000;%飞机飞行高度
path='..\input\descending_left.csv';
descending_left=csvread(path);
[lat, lon, H] = xyz2ell(1000*descending_left(:, 2:4));
% H = ones(size(H))*5000;
H = ones(size(H))*air_height;
arg = zeros(1,53,3); 
arg(1,:,:) = [lon/180*pi lat/180*pi H];
xyz = Ell2xyz(arg);
XYZ = zeros(53,3);
XYZ(:,:) = xyz(1,:,:);
% XYZ(:,2) = - XYZ(:, 2);

%% 求解法向量法1
% B = 1000*[1;1;1];
% a = descending_left(26,2:4)*1000;
% b = descending_left(30,2:4)*1000;
% c = XYZ(26,:);
% d = XYZ(30,:);
% A = [a; b; c];
% vector_n = A\B;
% vector_n = vector_n(1:3)/sqrt(sum(vector_n(1:3).^2));

%% 求解法向量法2
a = descending_left(26,2:4)*1000;
b = descending_left(30,2:4)*1000;
c = XYZ(26,:);
d = XYZ(30,:);
p1 = a - b;
p2 = c - b;
vector_n = cross(p1, p2);
vector_n = vector_n(1:3)/sqrt(sum(vector_n(1:3).^2));

%% 
DEM_center_xyz = c + vector_n*air_height;
DEM_center_ell = zeros(1,3);
[DEM_center_ell(2), DEM_center_ell(1), DEM_center_ell(3)] = xyz2ell(DEM_center_xyz);
%% 
table1 = csvread(path);
table1(:,2:4) = XYZ/1000;
csvwrite('../input/descending_left_airborn.csv', table1)