%% 图像分块2000*2000
row_index = [1 2001 4001 6001 8001];
col_index = [1 2001 4001 6001 8001 10001 12001 14001 16001];
wrapped_phase = wrapped_phase13(1:10000, 1:18000);
delta = 1999;
wrapped_phase_filtered = single(zeros(size(wrapped_phase)));

%% 深度学习滤波

format compact;
addpath('D:\Rick\深度学习滤波\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab\mex');
addpath('D:\Rick\深度学习滤波\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab\simplenn');
addpath('D:\Rick\深度学习滤波\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab');
showResult  = 1;
useGPU      = 0;
pauseTime   = 0;

modelName   = 'model_phase_filter_phase2phase_20layer';
epoch       = 42;



%%% load Gaussian denoising model
load(fullfile('data',modelName,[modelName,'-epoch-',num2str(epoch),'.mat']));
net = vl_simplenn_tidy(net);
net.layers = net.layers(1:end-1);

net = vl_simplenn_tidy(net);
h = waitbar(0);
for ii = 1:length(row_index)
    for jj = 1:length(col_index)
        input = wrapped_phase13(row_index(ii):row_index(ii)+delta,col_index(jj):col_index(jj)+delta);

        input1= cos(input);
        input2=sin(input);
        input(:,:,1,1)=input1;
        input(:,:,2,1)=input2;
        clear input1 input2

         res    = vl_simplenn(net,input,[],[],'conserveMemory',true,'mode','test');

         output = input - res(end).x;
         wrapped_phase_filtered(row_index(ii):row_index(ii)+delta,col_index(jj):col_index(jj)+delta) = single(angle(output(:,:,1)+1i*output(:,:,2)));
%          figure;imagesc(angle(output(:,:,1)+1i*output(:,:,2)));colormap(jet);
        waitbar((ii*length(col_index)+jj)/(length(col_index)*length(row_index)));
    end
end

figure;imagesc(wrapped_phase_filtered);colormap(jet);


%%  斜坡自适应滤波
wrapped_phase_filtered_para101 = Improved_Slope_Adaptive_filter_parallel(wrapped_phase(1:500,1:500),51,51);
