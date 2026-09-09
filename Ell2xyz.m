function [DEM ] = Ell2xyz(DEMll)
%#demll为6001*6001*3的矩阵 经度*维度*高 long*lat*height
        %矩阵操作转为x y z
        %(输入为弧度)
        %DEMll[:2,:,:]=DEMll[:2,:,:]/180*math.pi %角度转弧度
 Ea=6378137.0;
 Eb=6356752.3141;
 e2=0.00669438003551279091;
 N=Ea./sqrt(1.0-e2*(sin(DEMll(:,:,2))).^2);
 Nph=N+DEMll(:,:,3);
 x=Nph.*cos(DEMll(:,:,2)).*cos(DEMll(:,:,1));
 y=Nph.*cos(DEMll(:,:,2)).*sin(DEMll(:,:,1));
 z=(Nph-e2*N).*sin(DEMll(:,:,2));
 DEM(:,:,1)=x;
 DEM(:,:,2)=y;
 DEM(:,:,3)=z;
end

