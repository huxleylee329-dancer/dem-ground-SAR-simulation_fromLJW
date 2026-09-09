function [ Image,controlpoint,Slantrange,Imageslantrange,NewX,NewY,numcount, ImagePosX, ImagePosY, ImagePosZ,Vs_X,Vs_Y,Vs_Z,DEM_X,DEM_Y,DEM_Z] = ...
generateSLC1226(satelipos,fd0posn,DEM,sigma,RandAngle,Wlength,Drange,DAzimath,ImageAzimath,ImageRange,vpos)
%根据图像斜距以及地面散射系数分布计算图像复数据;
number=0;
ShiftAzim=ImageAzimath/2;
Shiftrange=ImageRange/2;
row=size(DEM,1);
clomn=size(DEM,2);
Imageslantrange=zeros(ImageAzimath,ImageRange);
numcount=zeros(ImageAzimath,ImageRange);
NewX=satelipos-ones(size(satelipos)).*satelipos(floor(row/2),floor(clomn/2))+ShiftAzim;
Slantrange=sqrt(sum((fd0posn-DEM).^2,3));
NewY=floor((Slantrange-ones(size(Slantrange)).*Slantrange(floor(row/2),floor(clomn/2)))./Drange);
NewY=NewY+Shiftrange;
phase_re=cos(RandAngle).*cos((-4)*pi./Wlength*Slantrange)-sin(RandAngle).*sin((-4)*pi./Wlength*Slantrange);
phase_im=cos(RandAngle).*sin((-4)*pi./Wlength*Slantrange)+sin(RandAngle).*cos((-4)*pi./Wlength*Slantrange);
% phase_re=cos((-4)*pi./Wlength*Slantrange);
% phase_im=sin((-4)*pi./Wlength*Slantrange);
phase_re=awgn(phase_re,100);%这里相当于热噪声；y = awgn(x,snr)将白高斯噪声添加到向量信号x中，在此表示加入随机相位
phase_im=awgn(phase_im,100);

sigma_1=(sqrt(sigma).*complex(phase_re,phase_im));
sigma_1real=real(sigma_1);
sigma_1imag=imag(sigma_1);
Image1=zeros(ImageRange,ImageAzimath);
Image2=zeros(ImageRange,ImageAzimath);


ImagePosX=zeros(size(Image1));
ImagePosY=zeros(size(Image1));
ImagePosZ=zeros(size(Image1));

Vs_X=zeros(size(Image1));
Vs_Y=zeros(size(Image1));
Vs_Z=zeros(size(Image1));

DEM_X=zeros(size(Image1));
DEM_Y=zeros(size(Image1));
DEM_Z=zeros(size(Image1));

for i=1:row    
    for j=1:clomn
        row1=NewX(i,j);
        clomn1=NewY(i,j);

        Imageslantrange(NewX(i,j),NewY(i,j))=Imageslantrange(NewX(i,j),NewY(i,j))+Slantrange(i,j);
        numcount(NewX(i,j),NewY(i,j))=numcount(NewX(i,j),NewY(i,j))+1;
        
        ImagePosX(NewX(i,j),NewY(i,j))=ImagePosX(NewX(i,j),NewY(i,j))+fd0posn(i,j,1);
        ImagePosY(NewX(i,j),NewY(i,j))=ImagePosY(NewX(i,j),NewY(i,j))+fd0posn(i,j,2);
        ImagePosZ(NewX(i,j),NewY(i,j))=ImagePosZ(NewX(i,j),NewY(i,j))+fd0posn(i,j,3);
        Vs_X(NewX(i,j),NewY(i,j))=Vs_X(NewX(i,j),NewY(i,j))+vpos(i,j,1);
        Vs_Y(NewX(i,j),NewY(i,j))=Vs_Y(NewX(i,j),NewY(i,j))+vpos(i,j,2);
        Vs_Z(NewX(i,j),NewY(i,j))=Vs_Z(NewX(i,j),NewY(i,j))+vpos(i,j,3);
        
        DEM_X(NewX(i,j),NewY(i,j))=DEM_X(NewX(i,j),NewY(i,j))+DEM(i,j,1);
        DEM_Y(NewX(i,j),NewY(i,j))=DEM_Y(NewX(i,j),NewY(i,j))+DEM(i,j,2);
        DEM_Z(NewX(i,j),NewY(i,j))=DEM_Z(NewX(i,j),NewY(i,j))+DEM(i,j,3);
        
        if mod(i,30)==0&&mod(j,30)==0
            number=number+1;
            controlpoint(number,1)=NewX(i,j);
            controlpoint(number,2)=NewY(i,j);
            controlpoint(number,3)=Slantrange(i,j);
            controlpoint(number,4)=DEM(i,j,1);%存取DEMx;
            controlpoint(number,5)=DEM(i,j,2);%存取DEMy;
            controlpoint(number,6)=DEM(i,j,3);%存取DEMz;
            sigma_1(i,j)=(sqrt(sigma(i,j))*complex(phase_re(i,j),phase_im(i,j)));
            sigma_1real(i,j)=real(sigma_1(i,j));
            sigma_1imag(i,j)=imag(sigma_1(i,j));
        end
        Image1(row1,clomn1)=Image1(row1,clomn1)+sigma_1real(i,j);
        Image2(row1,clomn1)=Image2(row1,clomn1)+sigma_1imag(i,j);
    end
end
Image=complex(Image1,Image2);
Imageslantrange=Imageslantrange./numcount;

ImagePosX=ImagePosX./numcount;
ImagePosY=ImagePosY./numcount;
ImagePosZ=ImagePosZ./numcount;

Vs_X=Vs_X./numcount;
Vs_Y=Vs_Y./numcount;
Vs_Z=Vs_Z./numcount;

DEM_X=DEM_X./numcount;
DEM_Y=DEM_Y./numcount;
DEM_Z=DEM_Z./numcount;

for i=1:size(controlpoint,1)
      controlpoint(i,3)=Imageslantrange(controlpoint(i,1),controlpoint(i,2));
end
%对得到的数据进行二维卷积
disp('进行复数据二维度卷积');
% Fs=200*10^6;
% Ba=Prf;
% Bw=199*10^6;
% Nr=size(Image,2);
% Na=size(Image,1);
% NBw = floor(Bw/Fs * Nr / 2);
% NBa = floor(Ba/Prf * Na / 2);
% Filter_Azi = zeros(1,Na);
% Filter_Azi(1:NBa) = ones(1,NBa);
% Filter_Azi(Na-NBa+1:Na) = ones(1,NBa);
% Filter_Range = zeros(1,Na);
% Filter_Range(1:NBw) = ones(1,NBw);
% Filter_Range(Nr-NBw+1:Nr) = ones(1,NBw);
% %在频域进行滤波，等效于在时域卷积，生成SAR图像
% Image=ifft2( fft2(Image) .* (Filter_Azi.'*ones(1,Nr)) .*(ones(Na,1)*Filter_Range));
WrightNum=8;
Wright=zeros(2*WrightNum+1,2*WrightNum+1);
ResY=0.333;
ResX=0.333;


for k=-WrightNum:WrightNum
    r=Drange*k;
    if k==0
        bb=1;
    else
        bb=sin(pi*r/ResY)/(pi*r/ResY);
    end
    for l=-WrightNum:WrightNum
        x=DAzimath*l;
        if l==0
            aa=1.0;
        else
            aa=sin(pi*x/ResX)/(pi*x/ResX);
        end
        Wright(WrightNum+k+1,WrightNum+l+1)=aa*bb; 
    end
end
Image=conv2(Image,Wright,'same');
end

