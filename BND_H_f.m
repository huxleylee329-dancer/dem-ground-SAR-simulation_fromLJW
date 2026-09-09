function [H_f] = BND_H_f(n,T,fi,f,phi)
%   此函数为半牛顿迭代法里面的表达式H_f
H_f1=exp(1i*2*pi.*n*T.*fi+1i*phi-1i*2*pi.*n*T.*f);
H_f=sum(H_f1);
end

