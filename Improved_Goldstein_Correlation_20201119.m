function Output_Value = Improved_Goldstein_Correlation_20201119(Array, ArrayGamma)

[Size_i,Size_j]=size(Array);
Window=zeros(Size_i,Size_j);

Radius_i=(Size_i-1)/2;
Radius_j=(Size_j-1)/2;

Window(Radius_i:Radius_i+2,Radius_j:Radius_j+2)=1;

% 设定加权系数
% r=0.5;
r=1-ArrayGamma;%根据论文第54页的介绍，ArrayGamma代表平均相干系数，r代表滤波参数alpha
%FFT变换
ArrayFFT=fftshift(fft2(Array));
s=(abs(ArrayFFT).*Window).^r;
Array=ifft2(ifftshift(ArrayFFT.*s));
Output_Value=Array(Radius_i+1,Radius_j+1); 
end

        
%% 传统Goldstein_filter

% function Output_Value = Afilter_Goldstein_Correlation(Array_Input, ArrayGamma)
% Win_Size=floor(1/ArrayGamma)+5;
% 
% if (mod(Win_Size,2)==0)
%     Win_Size=Win_Size+1;
% end
% Window=zeros(Win_Size,Win_Size);
% 
% Radius_i=(Win_Size-1)/2;
% Radius_j=(Win_Size-1)/2;
% Window(Radius_i:Radius_i+2,Radius_j:Radius_j+2)=1;
% 
% [Size_i,Size_j]=size(Array_Input);
% a=(Size_i-Win_Size)/2;
% b=(Size_j-Win_Size)/2;
% Array=Array_Input(a+1:end-a,b+1:end-b);
% %设定加权系数
% % r=0.5;
% r=1-ArrayGamma;
% %FFT变换
% ArrayFFT=fftshift(fft2(Array));
% s=(abs(ArrayFFT).*Window).^r;
% Array=ifft2(ifftshift(ArrayFFT.*s));
% Output_Value=Array(Radius_i+1,Radius_j+1); 
% 
% end
