function [pixeloc,NewXY] = pixe_loc(satelipos,fd0pos,DEM,long,lat,Drange,ImageAzimath,ImageRange)
%计算图像复数据
% ImageAzimath=4096;   
% ImageRange=4096;
ShiftAzim=ImageAzimath/2;
Shiftrang=ImageRange/2;
row=size(DEM,1);
clomn=size(DEM,2);
NewX=satelipos-ones(size(satelipos)).*satelipos(round(row/2),round(clomn/2))+ShiftAzim;
Slantrange=sqrt(sum((fd0pos-DEM).^2,3));%计算出每个像素点的斜距长
NewY=floor((Slantrange-ones(size(Slantrange)).*Slantrange(round(row/2),round(clomn/2)))/Drange);%Drange为斜距方向采样间隔
disp('中心斜距为');
Slantrange(round(row/2),round(clomn/2))
NewY=NewY+Shiftrang;
pixeloc=zeros(ImageRange,2,ImageAzimath);
pixeloc1=zeros(ImageRange,ImageAzimath);
pixeloc2=zeros(ImageRange,ImageAzimath);
numcount=zeros(ImageAzimath,2,ImageRange);
numcount1=zeros(ImageAzimath,ImageRange);
for i=1:row
    for j=1:clomn
        row1=NewX(i,j);
        clomn1=NewY(i,j);
        %        pixeloc(clomn1,:,row1)=pixeloc(clomn1,:,row1)+[long(i,j),lat(i,j)];
        pixeloc1(clomn1,row1)=pixeloc1(clomn1,row1)+long(i,j);
        pixeloc2(clomn1,row1)=pixeloc2(clomn1,row1)+lat(i,j);
        numcount1(clomn1,row1)=numcount1(clomn1,row1)+1;
    end
end
numcount(:,1,:)=numcount1;
numcount(:,2,:)=numcount1;
pixeloc(:,1,:)=pixeloc1;
pixeloc(:,2,:)=pixeloc2;
pixeloc=pixeloc./(numcount+1e-6);
NewXY(:,:,1)=NewX;NewXY(:,:,2)=NewY;%NewX，NewY为SAR图像上的像素坐标
end

