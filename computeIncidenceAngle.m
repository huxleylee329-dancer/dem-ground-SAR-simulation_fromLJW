function [sigma] = computeIncidenceAngle(pos,DEM )
%计算卫星天线发射的入射波夹角并计算地面对应的散射系数矩阵
a=size(DEM,1);b=size(DEM,2);
sigma=zeros(a,b);
v1=DEM(2:end,2:end,:)-DEM(1:a-1,1:b-1,:);
v2=DEM(2:end,2:end,:)-DEM(1:a-1,2:end,:);
v3=zeros(size(v1));
% for i=1:size(v1,1)
%     for j=1:size(v1,2)
%        v3(i,j,:)=cross(v2(i,j,:),v1(i,j,:));
%     end
% end
v3=cross(v2,v1);%C=CROSS(A,B)返回向量叉积baiduA和B，即，C = A x B

Terrain=DEM(1:a-1,2:end,:);
heading=sum(Terrain.*v3,3);
heading=double((heading>0));
headingtemp=zeros(size(v3));
headingtemp(:,:,1)=heading;headingtemp(:,:,2)=heading;headingtemp(:,:,3)=heading;
v3=v3.*headingtemp;
Slantrange=pos(1:a-1,2:end,:)-Terrain;
multi=sum(Slantrange.*v3,3);
heading=double(multi>=0);
RangeCs=(sqrt(sum(Slantrange.^2,3)).*sqrt(sum(v3.^2,3)));
AngleCs=acos(multi./RangeCs)/pi*180;%求出cs点与该点的斜距之间的夹角
SigmaCs=reflectivity(AngleCs);%利用子函数求出不同入射夹角对应的散射系数
SigmaCs=10.^(SigmaCs/10);
SigmaCs=SigmaCs.*heading;
sigma(1:a-1,2:end)=SigmaCs;
end

