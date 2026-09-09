function [ fd0pos,fdmin,temppos ] = computefd(DEM,satepos,gap)
%迭代寻找地面场景点的0多普勒成像卫星位置，并计算其与地面场景点之间的斜距
row1=size(DEM,1);clomn1=size(DEM,2);
fd0pos(:,:,1)=zeros(row1,clomn1);fd0pos(:,:,2)=zeros(row1,clomn1);fd0pos(:,:,3)=zeros(row1,clomn1);
fdmin=ones(row1,clomn1).*(-1);
satexyz=satepos(:,1:3);satev=satepos(:,4:end);
DEM2D=zeros(size(satepos,1),3);
DEM2D(:,1)=DEM(1,1,1);DEM2D(:,2)=DEM(1,1,2);DEM2D(:,3)=DEM(1,1,3);%将一个三维矩阵赋值为一个平面矩阵
srange=satexyz-DEM2D;%计算第一个DEM的斜距
fd=abs(sum(srange.*satev,2));
[firstpos,~]=find(fd==min(fd));
firstpos %%%%%输出第一个位置点
pre_clomn=firstpos;
pre_pos=firstpos;
for i=1:row1
    for j=1:clomn1
        if j==1
            pre_pos=pre_clomn;
        end
            lower=max(0,pre_pos-gap);
            upper=min(size(satepos,1),pre_pos+gap);
            satexyz=satepos(lower:upper,1:3);
            satev=satepos(lower:upper,4:end);
            DEM2D=zeros(size(satexyz));
            DEM2D(:,1)=DEM(i,j,1);
            DEM2D(:,2)=DEM(i,j,2);
            DEM2D(:,3)=DEM(i,j,3);
            srange=satexyz-DEM2D;
            fd=abs(sum(srange.*satev,2));
            [postemp,~]=find(fd==min(fd));
            pos=postemp+lower-1;
            temppos(i,j)=pos;
            fd0pos(i,j,1)=satepos(pos,1);fd0pos(i,j,2)=satepos(pos,2);fd0pos(i,j,3)=satepos(pos,3);
            fdmin(i,j)=min(fd);
            pre_pos=pos;
            if j==1
               pre_clomn=pos;
            end
    end
end
                            
end

