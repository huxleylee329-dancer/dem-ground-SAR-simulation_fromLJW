function Image = generateSLC_ground( ...
    R, U, alpha, lambda, ...
    rAxis, uAxis, rho_r, rho_u, Rref)

% 根据地形散射点直接生成复数雷达图像。
% 不生成原始回波。
%
% 第一版模型：
%   均匀带宽、均匀连续等效孔径；
%   距离和方向响应暂用可分离 sinc 函数。
%
% 图像行：方向余弦 U
% 图像列：斜距 R

    R = R(:);
    U = U(:);
    alpha = alpha(:);

    assert(numel(R)==numel(U) && numel(R)==numel(alpha), ...
        '距离、方向与散射系数的数量必须一致。');

    Nr = numel(rAxis);
    Nu = numel(uAxis);

    dr = rAxis(2)-rAxis(1);
    du = uAxis(2)-uAxis(1);

    Image = complex(zeros(Nu,Nr));

    % 主辅图像必须使用相同的 Rref
    coeff = alpha .* exp(-1j*4*pi*(R-Rref)/lambda);

    % 每个点仅更新附近像素
    support = 4;

    for p = 1:numel(R)

        j0 = max(1,ceil( ...
            (R(p)-support*rho_r-rAxis(1))/dr)+1);

        j1 = min(Nr,floor( ...
            (R(p)+support*rho_r-rAxis(1))/dr)+1);

        i0 = max(1,ceil( ...
            (U(p)-support*rho_u-uAxis(1))/du)+1);

        i1 = min(Nu,floor( ...
            (U(p)+support*rho_u-uAxis(1))/du)+1);

        if i0>i1 || j0>j1
            continue;
        end

        ii = i0:i1;
        jj = j0:j1;

        kr = local_sinc((rAxis(jj)-R(p))/rho_r);
        ku = local_sinc((uAxis(ii)-U(p))/rho_u);

        response = ku(:)*kr(:).';

        Image(ii,jj) = Image(ii,jj) ...
                    + coeff(p)*response;
    end
end

function y = local_sinc(x)
% MATLAB 归一化 sinc 定义：sin(pi*x)/(pi*x)

    y = ones(size(x));
    nz = x~=0;

    y(nz) = sin(pi*x(nz))./(pi*x(nz));
end