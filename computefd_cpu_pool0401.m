function [ fd0pos,fdmin,temppos,vpos ] = computefd_cpu_pool0401(DEM,satepos,gap)
%迭代寻找地面场景点的0多普勒成像卫星位置，并计算其与地面场景点之间的斜距
row1=size(DEM,1);clomn1=size(DEM,2);
fd0pos=zeros(row1,clomn1,3);
fdmin=ones(row1,clomn1).*(-1);
satexyz=satepos(:,1:3);satev=satepos(:,4:end);
temppos=zeros(row1,clomn1);


% parpool(12)
parfor i=1:row1
    for j=1:clomn1      
        p = [DEM(i,j,1) DEM(i,j,2) DEM(i,j,3)];
        srange=satexyz-p;
        fd=abs(sum(srange.*satev,2));
        [pos,~]=find(fd==min(fd));
        temppos(i,j)=pos;
        fd0pos(i,j,:)=satepos(pos,1:3);
        vpos(i,j,:)=satepos(pos,4:6);
        fdmin(i,j)=min(fd);
    end
end
% delete(gcp)
end