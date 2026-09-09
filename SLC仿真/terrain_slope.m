function [slope_a,slope_r] = terrain_slope(DEM,Dx,Dy)
% Horn法计算方位、距离坡度
% DEM边缘点可以用一阶差分、二阶差分计算其坡度
[Nta,Ntr] = size(DEM);
slope_a = zeros(Nta,Ntr);
slope_r = zeros(Nta,Ntr);
for ii = 2:1:Nta-1
    for jj = 2:1:Ntr-1
        fx = (DEM(ii-1,jj-1)-DEM(ii+1,jj-1)+2*(DEM(ii-1,jj)-DEM(ii+1,jj))+DEM(ii-1,jj+1)-DEM(ii+1,jj+1))/8/Dx;
        fy = (DEM(ii+1,jj+1)-DEM(ii+1,jj-1)+2*(DEM(ii,jj+1)-DEM(ii,jj-1))+DEM(ii-1,jj+1)-DEM(ii-1,jj-1))/8/Dy;
        slope_a(ii,jj) = atan(fx);%方位坡度,rad
        slope_r(ii,jj) = atan(fy);%距离坡度,rad
    end
end
end