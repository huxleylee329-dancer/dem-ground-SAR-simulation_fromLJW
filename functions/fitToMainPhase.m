function [Y,N] = fitToMainPhase(X,varargin) %(minVal,maxVal)
% 将相位归整到-pi到pi
    if nargin >= 2
        minVal = varargin{1}(1);
        maxVal = varargin{1}(2);
        if minVal > maxVal
            [maxVal,minVal] = deal(minVal,maxVal);
        end
    else
        minVal = -pi;
        maxVal = pi;
    end
    
    N = floor((X - minVal)/(maxVal - minVal));
    
    Y = X - N*(maxVal - minVal);
end