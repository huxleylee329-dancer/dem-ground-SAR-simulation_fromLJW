function sigma = terrain_sigma(lookangle,slope_a,slope_r,Dx,Dy)
% 局部入射角
[Nta,Ntr] = size(slope_a);
theta = lookangle;  %卫星视角，rad
theta_inc = acos(cos(theta) + tan(slope_r).*sin(theta)./sqrt((1+tan(slope_r).^2+tan(slope_a).^2)));
% 散射系数
% 沙漠、裸土
p1 = [4.33];
p2 = [6.666];
p3 = [-0.107];
p4 = [-29.709];
p5 = [0.863];
p6 = [-1.365];
delta_A = Dx * Dy;
sigma_0 = p1 + p2*exp(p3*theta_inc) + p4*cos(p5*theta_inc + p6);
sigma   = sigma_0 * delta_A;
% 判断阴影
for ii =1:1:Nta
    for jj =1:1:Ntr
        if slope_r(ii,jj) < 0 && slope_r(ii,jj) < (theta-pi/2)
            sigma(ii,jj) = 0;
        end
    end
end
end