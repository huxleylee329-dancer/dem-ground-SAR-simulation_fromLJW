function phi_nl = phi_nl_filt(phi_res,varargin)
    % [K,H] = size(phi_res);
    if nargin == 1
        phi_nl = sgolayfilt(phi_res,2,5); % 2次曲线，FIR滤波器长度为5
    else
        phi_nl = sgolayfilt(phi_res,varargin{1}(1),varargin{1}(2));
    end
end