function [image_interp] = interp_paddingzero(image,interptimes)
[nr,nc] = size(image);
temp1 = zeros(nr*interptimes,nc*interptimes);
temp2 = fft2(image);%先对图像进行傅里叶变换

temp1(1:nr/2,1:nc/2) = temp2(1:nr/2,1:nc/2);
temp1(nr*interptimes-nr/2+1:nr * interptimes,1:nc/2) = temp2(nr/2+1:nr,1:nc/2);
temp1(1:nr/2,nc*interptimes-nc/2+1:nc * interptimes) = temp2(1:nr/2,nc/2+1:nc);
temp1(nr*interptimes-nr/2+1:nr * interptimes,nc*interptimes-nc/2+1:nc * interptimes) = temp2(nr/2+1:nr,nc/2+1:nc);

image_interp = ifft2(temp1);
end