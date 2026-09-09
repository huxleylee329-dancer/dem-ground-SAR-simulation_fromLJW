function [H_f_out] = diff2_Newton_H(ii,N,f,Data,T)
%¶þ½×µ¼ 
y_i=zeros(1,2*N+1);
H_f_out2=zeros(1,2*N+1);
y_i=Data(ii-N:ii+N);
j=sqrt(-1);
H_f_out2=y_i.*(2*pi*T)^2*(-1).*(-N:N).*(-N:N).*exp((-j*2*pi*T).*(-N:N)*f);
H_f_out=sum(H_f_out2,2);
end

