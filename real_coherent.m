function [move_r,move_c,normcor]=real_coherent(image1,image2,figshow)
if nargin<3,figshow=0;end
assert(isequal(size(image1),size(image2)));
[nr,nc]=size(image1);
A=abs(image1);B=abs(image2);
A=A-mean(A(:));B=B-mean(B(:));
assert(norm(A(:))>0 && norm(B(:))>0,'No amplitude texture.');
cor=real(fftshift(ifft2(fft2(A,2*nr-1,2*nc-1).*conj(fft2(B,2*nr-1,2*nc-1)))));
[peak,index]=max(cor(:));assert(peak>0,'No positive correlation peak.');
[r,c]=ind2sub(size(cor),index);
move_r=nr-r;move_c=nc-c;
normcor=cor/peak;
if figshow,figure;imagesc(normcor);colorbar;end
end
