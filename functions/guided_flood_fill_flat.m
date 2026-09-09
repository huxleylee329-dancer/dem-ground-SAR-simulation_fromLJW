% GuidedFloodFill.m unwraps a single 2D image using the quaility guided technique
% function IM_unwrapped = GuidedFloodFill(IM_phase, IM_mag, IM_unwrapped, unwrapped_binary, derivative_variance, adjoin, IM_mask)
% It can also be used to unwrap phases in cases when there are no branch cuts.
%
% Input: IM_phase, IM_unwrapped (seed points / pixels already unwrapped),
% unwrapped_binary the derivative variance, an adjoining matrix and a mask.
% Inputs:
%  IM_phase     = 2D array of wrapped phases (rads)
%  IM_mag       = 2D array of magnitudes
%  IM_unwrapped = a 2D array for the unwrapped phases (rads)
%                 This is initialized for the reference point set.
%  unwrapped_binary = a 2d array to identifying the points that have been unwrapped.
%                 This is initialized for the reference point set.
%  derivative_variance = a 2D array of the IM_phase's derivative variances
%  adjoin       = a 2D identifying the valid pixels adjoining the unwrapped pixels.
%                 This is initialized for the reference point set.
%  IM_mask      = a 2D array identifying the valid pixels
% Outputs:
%  IM_unwrapped = the 2D array of unwrapped phases (rads)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IM_unwrapped] = guided_flood_fill_flat(IM_phase, IM_unwrapped, unwrapped_binary, adjoin)

[r_dim, c_dim] = size(IM_phase);% 获取缠绕相位尺寸

% indicator = 101;%101 %--------------------------------指导值
% Include edge pixels
% 以毗邻矩阵中全零时为终止迭代的条件

while sum(sum(adjoin(:))) ~= 0  %Loop until there are no more adjoining pixels    
 [r_adjoin,c_adjoin] = find(adjoin);
  for(ii = 1:length(r_adjoin))
    r_active = r_adjoin(ii);
    c_active = c_adjoin(ii);
  if(r_active+1<=r_dim)  % 判断活跃点，active是否在图像范围内
     if unwrapped_binary(r_active+1, c_active)==1  && unwrapped_binary(r_active,c_active)==0% 如果该活跃点已经解缠
            phase_ref = IM_unwrapped(r_active+1, c_active);       % Obtain the reference unwrapped phase
             % Itoh's Method (suggested by 'Eric' on MATLAB Central to use for a length 2 vector):
             D = IM_phase(r_active, c_active)-phase_ref;
             deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
             IM_unwrapped(r_active, c_active) = phase_ref + deltap;  % This is the unwrapped phase
             adjoin(r_active,c_active) = 0;
             unwrapped_binary(r_active, c_active)=1;
     elseif unwrapped_binary(r_active+1, c_active)~=1
        adjoin(r_active+1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
    end
  end
  
  % 第二，搜索最小相位方差点（活跃点，active）上方
  %Then search above
  if(r_active-1>=1)  % 判断活跃点，active是否在图像范围内
     if unwrapped_binary(r_active-1, c_active)==1  && unwrapped_binary(r_active,c_active)==0% 如果该活跃点已经解缠
  
         phase_ref = IM_unwrapped(r_active-1, c_active);                                   %Obtain the reference unwrapped phase
      D = IM_phase(r_active, c_active)-phase_ref;
      deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
      IM_unwrapped(r_active, c_active) = phase_ref + deltap; 
      adjoin(r_active,c_active) = 0;
      unwrapped_binary(r_active, c_active)=1;
     elseif unwrapped_binary(r_active-1, c_active)~=1
        adjoin(r_active-1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
    end
  end
  
  % 第三，搜索最小相位方差点（活跃点，active）右侧
  %Then search on the right
  if(c_active+1<=c_dim)  % 判断活跃点，active是否在图像范围内
     if unwrapped_binary(r_active, c_active+1)==1&& unwrapped_binary(r_active,c_active)==0

      phase_ref = IM_unwrapped(r_active, c_active+1);                                   %Obtain the reference unwrapped phase
      D = IM_phase(r_active, c_active)-phase_ref;
      deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
      IM_unwrapped(r_active, c_active) = phase_ref + deltap; 
      adjoin(r_active,c_active) = 0;
      unwrapped_binary(r_active, c_active)=1;
     elseif unwrapped_binary(r_active, c_active+1)~=1
        adjoin(r_active, c_active+1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
    end
  end
  
  % 第四，搜索最小相位方差点（活跃点，active）左侧
  %Finally search on the left
  if(c_active-1>=1)  % 判断活跃点，active是否在图像范围内
     if unwrapped_binary(r_active, c_active-1)==1 && unwrapped_binary(r_active,c_active)==0

      phase_ref = IM_unwrapped(r_active, c_active-1);                                   %Obtain the reference unwrapped phase
      D = IM_phase(r_active, c_active)-phase_ref;
      deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
      IM_unwrapped(r_active, c_active) = phase_ref + deltap; 
      adjoin(r_active,c_active) = 0;
      unwrapped_binary(r_active, c_active)=1;
     elseif unwrapped_binary(r_active, c_active-1)~=1
        adjoin(r_active, c_active-1) = 1;  
    end
  end
  end
 
  %end
end % while sum(sum(adjoin(2:r_dim-1,2:c_dim-1))) ~= 0  %Loop until there are no more adjoining pixels
% disp(['All of the valid interior pixels have been calculated']);
