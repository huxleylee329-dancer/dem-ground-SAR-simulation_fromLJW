%ÑµÁ·Ö÷º¯Êý
%%% Note: run the 'GenerateData_model_64_25_Res_Bnorm_Adam.m' to generate
%%% training data first.

addpath /usr/local/MATLAB/R2017a/toolbox/matconvnet-1.0-beta25/matlab;
vl_setupnn;
%%%-------------------------------------------------------------------------
%%% configuration
%%%-------------------------------------------------------------------------
opts.modelName        = 'model_phase_filter_phase2phase_20layer'; %%% model name
opts.learningRate     = [logspace(-3,-3,30) logspace(-4,-4,70)];%%% you can change the learning rate
 %opts.learningRate     = [logspace(-3,-4,60)];
opts.batchSize        = 32; %%% default
opts.gpus             = [1]; %%% this code can only support one GPU!

%%% solver
opts.solver           = 'Adam';

opts.gradientClipping = false; %%% Set 'true' to prevent exploding gradients in the beginning.
opts.expDir      = fullfile('data', opts.modelName);
opts.imdbPath    = fullfile(opts.expDir, 'imdb.mat');

%%%-------------------------------------------------------------------------
%%%   Initialize model and load data
%%%-------------------------------------------------------------------------
%%%  model
net  = feval(['DnCNN_init_',opts.modelName]);

%%%  load data
imdb = load(opts.imdbPath) ;

%%%-------------------------------------------------------------------------
%%%   Train 
%%%-------------------------------------------------------------------------

[net, info] = DnCNN_train(net, imdb, ...
    'expDir', opts.expDir, ...
    'learningRate',opts.learningRate, ...
    'solver',opts.solver, ...
    'gradientClipping',opts.gradientClipping, ...
    'batchSize', opts.batchSize, ...
    'modelname', opts.modelName, ...
    'gpus',opts.gpus) ;






