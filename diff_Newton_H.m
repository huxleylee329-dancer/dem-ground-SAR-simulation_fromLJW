function [H_f_out] = diff_Newton_H(ii,N,f,Data,T)
%Ò»½×µ¼
y_i=zeros(1,2*N+1);
H_f_out1=zeros(1,2*N+1);
y_i=Data(ii-N:ii+N);
j=sqrt(-1);
H_f_out1=y_i.*(-j*2*pi*T).*(-N:N).*exp((-j*2*pi*T).*(-N:N)*f);
H_f_out=sum(H_f_out1,2);
end

