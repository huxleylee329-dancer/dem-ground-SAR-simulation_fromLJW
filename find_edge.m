function [startpix,endpix,cutImage] = find_edge( piexloc,satelipos)
%%该函数表示对仿真SAR图像寻找边缘并进行截图操作
row=size(satelipos,1);
clomn=size(satelipos,2);
prow=size(piexloc,3);
pclomn=size(piexloc,1);
cutImage=((abs(piexloc(:,1,:))>0))&((abs(piexloc(:,2,:))>0));
cutImage=permute(cutImage,[3 1 2]);%B=permute(A,order)；对N维数组A按照指定的向量order顺序来重新排列其维数
center=[round(prow/2),round(pclomn/2)];
gap=5;
while (cutImage(center(1)-gap,center(2)-gap)&cutImage(center(1)+gap,center(2)-gap)...
       &cutImage(center(1)-gap,center(2)+gap)&cutImage(center(1)+gap,center(2)+gap))
   gap=gap+5;
end
gap=gap-5;
% gap=330;
disp('截图完成')
disp('截图大小的长宽为')
2*gap
srow=center(1)-gap;
erow=center(1)+gap;
sclomn=center(2)-gap;
eclomn=center(2)+gap;
startpix=[srow,sclomn];
endpix=[erow,eclomn];
end

