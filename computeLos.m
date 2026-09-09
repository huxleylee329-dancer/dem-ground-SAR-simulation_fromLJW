function [LOS] = computeLos(pos,DEM)
Slantrange=pos-DEM;
temp1=sqrt(sum(Slantrange.^2,3));%sum(A,3)运算后的值为每个通道对应位置的值各自相加
temp(:,:,1)=temp1;temp(:,:,2)=temp1;temp(:,:,3)=temp1;
LOS=Slantrange./temp;%得到了每个斜距的单位矢量
end

