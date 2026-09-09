function [ satepos ] = compute_center(DEM,Prf,numsate,path)
%  计算轨道卫星中心位置对应的地面场景中心点
 num=20001;
 satetime=getsatetime(Prf,num,path);
 satepos=getSatePara(satetime,Prf,path);
 centertemp=DEM(round(size(DEM,1)/2),round(size(DEM,2)/2),:);%表示dem中心点
 center=zeros(size(satepos,1),3);
 center(:,1)=centertemp(1,1,1);center(:,2)=centertemp(1,1,2);center(:,3)=centertemp(1,1,3);
 satexyz=satepos(:,1:3);
 satev=satepos(:,4:end);
 srange=satexyz-center;%中心斜距等于卫星位置减去DEM中心点
 fd=abs(sum(srange.*satev,2));
 [firstpos,~]=find(fd==min(fd));
 disp('找到的成像中心位置及多普勒频率为')
 one_secd=1/2/Prf;
 satetime=(linspace(0,numsate-1,numsate)-(numsate-1)/2)*one_secd+satetime(firstpos);
 satepos=getSatePara(satetime,Prf,path);
end

