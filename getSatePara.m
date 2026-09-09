function [satepos] = getSatePara(satetime,Prf,path)
%读取轨道参数对应的x,y,z方向的位置坐标以及速度大小，并进行插值
descending_left=csvread([path]);
t=descending_left(:,1);%读取轨道时间，后续进行线性插值运算
fx=descending_left(:,2).*1000;fy=descending_left(:,3).*1000;fz=descending_left(:,4).*1000;
fvx=descending_left(:,5).*1000;fvy=descending_left(:,6).*1000;fvz=descending_left(:,7).*1000;%读取轨道速度
fx=interp1(t,fx,satetime','linear');
fy=interp1(t,fy,satetime','linear');
fz=interp1(t,fz,satetime','linear');
fvx=interp1(t,fvx,satetime','linear');
fvy=interp1(t,fvy,satetime','linear');
fvz=interp1(t,fvz,satetime','linear');
satepos(:,1)=fx;satepos(:,2)=fy;satepos(:,3)=fz;satepos(:,4)=fvx;satepos(:,5)=fvy;
satepos(:,6)=fvz;
end

