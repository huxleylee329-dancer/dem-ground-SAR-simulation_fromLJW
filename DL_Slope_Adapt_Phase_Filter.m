

function DL_phase_filter = DL_Slope_Adapt_Phase_Filter(phase)
%%% test the model performance
j=sqrt(-1);
winsize_e_i = 11;           %���ƴ�i����ߴ�
winsize_e_j = 11;           %���ƴ�j����ߴ�
winsize_f_i = 11;  %�˲���i����ߴ�
winsize_f_j =11;  %�˲���j����ߴ�
wine_i = (winsize_e_i-1)/2; %���ƴ�i����뾶
wine_j = (winsize_e_j-1)/2; %���ƴ�j����뾶
winf_i = (winsize_f_i-1)/2; %�˲���i����뾶
winf_j = (winsize_f_j-1)/2; %�˲���j����뾶
[size_i,size_j]=size(phase);
phase_update = zeros(size_i+winsize_e_i-1,size_j+winsize_e_j-1);             
[Ni,Nj] = size(phase_update);
phase_update(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j) = phase;
phase_update(1:wine_i,1:Nj) = flipud(phase_update(wine_i+1:2*wine_i,1:Nj)); 
phase_update(1:Ni,1:wine_j) = fliplr(phase_update(1:Ni,wine_j+1:2*wine_j));  
phase_update(Ni-wine_i+1:Ni,1:Nj) = flipud(phase_update(Ni-2*wine_i+1:Ni-wine_i,1:Nj));
phase_update(1:Ni,Nj-wine_j+1:Nj) = fliplr(phase_update(1:Ni,Nj-2*wine_j+1:Nj-wine_j));
% ����λת��Ϊ����
phase_update = exp(j*phase_update);     
%% �ֲ�Ƶ�ʹ��Ʋ���
phase_linear = zeros(Ni,Nj); % ������λ����
DL_phase_filter = zeros(Ni,Nj); % �˲�����λ����
nn = 32;                     % 2D FFT����
mm = 32;                     % CZT����
mn = mm*nn;                  % Ƶ���ֱܷ��ʣ�mn*mn��
phase_estimation = zeros(winsize_e_i,winsize_e_j); %���ƴ�����λ����
tempi = (0:winsize_f_i-1)';  % ��������λ���
tempj = 0:winsize_f_j-1;     % ��������λ���
w0 = exp(-j*2*pi/mn);        % CZT����
%% ����Ӧ�˲�
%% ����Ӧ�˲�
wb = waitbar(0,'phase filtering ...');
% clear; clc;
format compact;
addpath(fullfile('data','utilities'));
% folderTest  = fullfile('data','Test','Set12'); %%% test dataset
useGPU      = 1;
modelName   = 'cos_model_phase_filter_phase2phase';
epoch       = 50;
%noiseSigma  = 40;  %%% image noise level
%%% load Gaussian denoising model
load(fullfile('data',modelName,[modelName,'-epoch-',num2str(epoch),'.mat']));
net = vl_simplenn_tidy(net);
net.layers = net.layers(1:end-1);
%%%
net = vl_simplenn_tidy(net);
%%% move to gpu
if useGPU
    net1 = vl_simplenn_move(net, 'gpu') ;
end

modelName   = 'sin_model_phase_filter_phase2phase';
epoch       = 50;
%noiseSigma  = 40;  %%% image noise level
%%% load Gaussian denoising model
load(fullfile('data',modelName,[modelName,'-epoch-',num2str(epoch),'.mat']));
net = vl_simplenn_tidy(net);
net.layers = net.layers(1:end-1);
%%%
net = vl_simplenn_tidy(net);
%%% move to gpu
if useGPU
    net2 = vl_simplenn_move(net, 'gpu') ;
end



for ii = wine_i+1:Ni-wine_i
  for jj = wine_j+1:Nj-wine_j
      phase_estimation(1:winsize_e_i,1:winsize_e_j) = phase_update(ii-wine_i:ii+wine_i,jj-wine_j:jj+wine_j);
      % FFTƵ�ʹ���
      phase_2fft = fftshift(fft2(phase_estimation,nn,nn)); %������2DFFT
      [val_j,ind_j] = max(max(abs(phase_2fft)));      %ind_j���ֵ�к� 
      [val_i,ind_i] = max(abs(phase_2fft(:,ind_j)));  %ind_i���ֵ�к� 
      fi = double(ind_i-nn/2-2)/nn;
      fj = double(ind_j-nn/2-2)/nn;
      ei = exp(j*2*pi*fi);
      ej = exp(j*2*pi*fj);
      % CZTƵ�����꾫ϸ����
      phasef_czt1 = czt(phase_estimation,3*mm,w0,ei); %i(��)����CZT,czt���м���
      phasef_czt2 = czt(phasef_czt1.',3*mm,w0,ej);    %j(��)����CZT
      phasef_czt2 = phasef_czt2.';
      [val_j,ind_j] = max(max(abs(phasef_czt2)));
      [val_i,ind_i] = max(abs(phasef_czt2(:,ind_j)));
      movei = double(ind_i-1)/mn;  %i(��)����ϸ����ƫ����
      movej = double(ind_j-1)/mn;  %j(��)����ϸ����ƫ����
      fii = fi+movei;
      fjj = fj+movej;       %CZT����ֱ��ʺ�Ƶ�׷�ֵλ��
      ei = exp(j*2*pi*fii); %i(��)����ֲ�Ƶ�ʣ�ָ��
      ej = exp(j*2*pi*fjj); %j(��)����ֲ�Ƶ�ʣ�ָ��
      % ������λ
      phase_linear(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j) = ei.^tempi*(ej.^tempj);    
      % ����,sum(sum()) or mean(mean())
      phase_0_matrix = phase_update(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j).*...
                       conj(phase_linear(ii-winf_i:ii+winf_i,jj-winf_j:jj+winf_j));
        
                   
       if useGPU
        input1 =cos( gpuArray(single(angle(phase_0_matrix))));
        input2 =sin( gpuArray(single(angle(phase_0_matrix))));
       end
      output = vl_simplenn(net1,input1,[],[],'conserveMemory',true,'mode','test');
      output1=output(end).x; 
      output1=input1-output1;
      output = vl_simplenn(net2,input2,[],[],'conserveMemory',true,'mode','test');
      output2=output(end).x;  
      output2=input2-output2;
       if useGPU
        phase_0_matrix = gather(output1)+j*gather(output2);
       end            
     % phase_0 =  mean(mean(phase_0_matrix)); %��ֵ�˲�  
      % ���� + ������λ 
      DL_phase_filter(ii,jj) =  phase_0_matrix(winf_i+1,winf_j+1) * phase_linear(ii,jj);
  end
 str = ['�˲�������',num2str(fix((100*ii/(Ni-2*wine_i)))),'%'];
 waitbar(ii/(Ni-2*wine_i),wb,str);

end
close(wb);
%% ȡ�����(rad),ȡ��Ч��Χ
DL_phase_filter = angle(DL_phase_filter);
DL_phase_filter = DL_phase_filter(wine_i+1:Ni-wine_i,wine_j+1:Nj-wine_j);
%% ���
clear phase_update phase_linear phase_estimation;



end






