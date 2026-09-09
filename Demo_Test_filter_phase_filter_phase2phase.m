


%%% test the model performance


format compact;
addpath('.\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab\mex');
addpath('.\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab\simplenn');
addpath('.\matconvnet-1.0-beta20\matconvnet-1.0-beta20\matlab');
% addpath(fullfile('data','utilities'));
% folderTest  = fullfile('data','Test','Set12'); %%% test dataset
% folderTest  = fullfile('data','Test','wrapped_phase_noise'); %%% test dataset
% folderTest_real  = fullfile('data','Test','wrapped_phase_real'); %%% test dataset
% folderTest  = 'J:\graduate_code\deep_learning\simulate_image\image2\nodefo\wrapped_noised_image'; %%% test dataset
% folderTest_real  = fullfile('data','Test','wrapped_phase_real'); %%% test dataset
showResult  = 1;
useGPU      = 0;
pauseTime   = 0;

modelName   = 'model_phase_filter_phase2phase_20layer';
epoch       = 42;%鍋氬疄楠岀殑鏃跺?鏄?5
%noiseSigma  = 40;  %%% image noise level


%%% load Gaussian denoising model
load(fullfile('data',modelName,[modelName,'-epoch-',num2str(epoch),'.mat']));
net = vl_simplenn_tidy(net);
net.layers = net.layers(1:end-1);

%%%
net = vl_simplenn_tidy(net);

% for i = 1:size(net.layers,2)
%     net.layers{i}.precious = 1;
% end

%%% move to gpu
% if useGPU
%     net = vl_simplenn_move(net, 'gpu') ;
% end

%%% read images
% ext         =  {'*.jpg','*.png','*.bmp','*.mat'};
% filePaths   =  [];
% for i = 1 : length(ext)
%     filePaths = cat(1,filePaths, dir(fullfile(folderTest,ext{i})));
% end
% filePaths = cat(1,filePaths, dir(fullfile(folderTest,ext{4})));
% filePaths_real   =  [];
% for i = 1 : length(ext)
%     filePaths_real = cat(1,filePaths_real, dir(fullfile(folderTest_real,ext{i})));
% end

%%% PSNR and SSIM
% PSNRs = zeros(1,length(filePaths));
% SSIMs = zeros(1,length(filePaths));
%i = 1:length(filePaths)
% h = waitbar(0);
% out_path = 'J:\graduate_code\deep_learning\simulate_image\image2\nodefo\wrapped_filter_image\';
for i = 1:1
%     image1 = load(fullfile(folderTest_real,filePaths_real(i).name));
%     label=single(image1.wrapped_phase_real);
%     image1 = load(fullfile(folderTest,filePaths(i).name));
%     input=single(image1.wrapped_phase_noise);  
    input = single(wrapped_phase);%！！！输入变量
    %%% convert to GPU
%     if useGPU
%         input = gpuArray(input);
%     end
%     tic
    input1= cos(input);
    input2=sin(input);
    input(:,:,1,1)=input1;
    input(:,:,2,1)=input2;
%     clear input1 input2

     res    = vl_simplenn(net,input,[],[],'conserveMemory',true,'mode','test');
    
     output = input - res(end).x;
%      save([out_path filePaths(i).name], 'output');
%      waitbar(i/length(filePaths));
%      toc
%        figure;imagesc(image1.wrapped_phase_noise)
      figure;imagesc(angle(output(:,:,1)+1i*output(:,:,2)))
      phase_filter=double(angle(output(:,:,1)+1i*output(:,:,2)));%！！！输出变量

%      figure;imagesc(angle(output(:,:,1)+1i*output(:,:,2))-label)
%      wrapped_phase_filter=gather(angle(output(:,:,1)+1i*output(:,:,2)));
    %%% convert to CPU
%     if useGPU
%         output = gather(output);
%         input  = gather(input);
%     end
%     close all
%      C=wrapped_phase_filter+4;
%      C=double(C(1:256,1:256));
%      Cmax=max(max(C));
%      Cmin=min(min(C));
%      A=mat2gray(C,[Cmin Cmax]);
%      imwrite(A,['./SHUJUJI/',num2str(i),'.jpg']);
end

% disp([mean(PSNRs),mean(SSIMs)]);




