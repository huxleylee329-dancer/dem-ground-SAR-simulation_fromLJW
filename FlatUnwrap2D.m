function [unwrapped_phase] = FlatUnwrap2D(wrapped_phase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QualityGuidedUnwrap2D implements 2D quality guided path following phase
% unwrapping algorithm.
% Inputs:  1. Complex image in .mat double format
%          2. Binary mask (optional)          
% Outputs: 1. Unwrapped phase image
%          2. Phase quality map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]    = size(wrapped_phase);

im_phase = wrapped_phase;       %Phase image
% -----------------------------------------------------------
%% Replace with your mask (if required)
unwrapped_phase = zeros(size(wrapped_phase));        %Initialze the output unwrapped version of the phase
adjoin = zeros(size(wrapped_phase));            %Zero starting matrix for adjoin matrix
unwrapped_binary = zeros(size(wrapped_phase));  %Binary image to mark unwrapped pixels

%% Automatically (default) or manually identify starting 
rowref =round(m/2); % choose the 1st point for a reference (known good value)
colref =round(n/2); % choose the 1st point for a reference (known good value)

%% Unwrap
unwrapped_phase(rowref,colref) = im_phase(rowref,colref);  %                        %Save the unwrapped values
unwrapped_binary(rowref,colref,1) = 1; %将解缠标记矩阵扩展为三维，且将seed点处标记为1，代表已解缠
% 将毗邻矩阵扩展为三维矩阵，将seed点的四邻点在毗邻矩阵中标记为1（毗邻矩阵初始为零）
adjoin(rowref-1, colref) = 1;  %Mark the pixels adjoining the selected point（seed point）
adjoin(rowref+1, colref) = 1; 
adjoin(rowref, colref-1) = 1;
adjoin(rowref, colref+1) = 1; 
[unwrapped_phase] = guided_flood_fill_flat(im_phase, unwrapped_phase, unwrapped_binary, adjoin);
end
