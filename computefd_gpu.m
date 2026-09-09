function [ fd0pos,fdmin,temppos ] = computefd_gpu(DEM,satepos,gap)
%迭代寻找地面场景点的0多普勒成像卫星位置，并计算其与地面场景点之间的斜距
row1=size(DEM,1);clomn1=size(DEM,2);
row2=size(satepos,1);
fd0pos=zeros(row1,clomn1,3);
fdmin=ones(row1,clomn1).*(-1);
satexyz=satepos(:,1:3);
satev=satepos(:,4:end);
temppos=zeros(row1,clomn1);

satexyz=gpuArray(satexyz);
fd0pos=gpuArray(fd0pos);
fdmin=gpuArray(fdmin);
satev=gpuArray(satev);
DEM=gpuArray(DEM);
temppos=gpuArray(temppos);

satex=repmat(satexyz(:,1),1,clomn1);
satey=repmat(satexyz(:,2),1,clomn1);
satez=repmat(satexyz(:,3),1,clomn1);
satevx=repmat(satev(:,1),1,clomn1);
satevy=repmat(satev(:,2),1,clomn1);
satevz=repmat(satev(:,3),1,clomn1);

for i=1:row1
        fd=abs((satex-repmat(DEM(i,:,1),row2,1)).*satevx+(satey-repmat(DEM(i,:,2),row2,1)).*satevy+(satez-repmat(DEM(i,:,3),row2,1)).*satevz);        
        [minfd,pos]=min(fd);
        temppos(i,:)=pos;
        fd0pos(i,:,:)=satexyz(pos,1:3);
        fdmin(i,:)=minfd;
end 
fd0pos=gather(fd0pos);
fdmin=gather(fdmin);
temppos=gather(temppos);
end