function [initial_val,exitflag] = findinitial(x_step,y_step,varargin)
x_step = sort(x_step,'ascend');
y_step = sort(y_step,'ascend');
Z = -minfunc(x_step,y_step);

if nargin == 3
    [r,c] = myfindpeaks(Z,varargin{1});
else
    [r,c] = myfindpeaks(Z);
end

if isempty(r)
    initial_val = zeros(2,1);
    exitflag = -1;
else
    x_peaks = x_step(r);
    y_peaks = y_step(c);
    
    % æ‡¿Î≈–∂œ◊º‘Ú
    distanceZero = abs(y_peaks);
    [~,idx] = min(distanceZero);
    initial_val = [x_peaks(idx);y_peaks(idx)];
    exitflag = 1;
end
end