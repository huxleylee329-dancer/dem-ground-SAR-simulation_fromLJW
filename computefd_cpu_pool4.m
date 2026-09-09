function [ fd0pos,temppos,vpos ] = computefd_cpu_pool4(DEM,satepos,theta_r)
%迭代寻找地面场景点的0多普勒成像卫星位置，并计算其与地面场景点之间的斜距
row1=size(DEM,1);clomn1=size(DEM,2);
fd0pos=zeros(row1,clomn1,3);
satexyz=satepos(:,1:3);satev=satepos(:,4:end);
temppos=zeros(row1,clomn1);


% parpool(12)
% for i=1:row1
parfor i=1:row1
    for j=1:clomn1
        srange=satexyz-repmat([DEM(i,j,1) DEM(i,j,2) DEM(i,j,3)],size(satexyz,1),1);
        %         fd=abs(sum(srange.*satev,2));
        %         [pos,~]=find(fd==min(fd));
        angle2=acos(sum(srange.*satev,2)./sqrt(abs(sum(srange.^2,2)))./sqrt(abs(sum(satev.^2,2))));%是否可改进
        temp2=abs(angle2-theta_r*pi/180);
        [pos,~]=find(temp2==min(temp2));%firstpos==1是否有问题？
        temppos(i,j)=pos;
        fd0pos(i,j,:)=satepos(pos,1:3);
        vpos(i,j,:)=satepos(pos,4:6);
    end
end
delete(gcp)
end