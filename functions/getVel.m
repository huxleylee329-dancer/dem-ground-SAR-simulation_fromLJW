function [V_xyz] = getVel(t,coef)
%-----------------get the position of the satellite at time t--------------
V_xyz = zeros(length(t),3);
V_xyz(:,1) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t t.*t.*t.*t.*t]*coef(:,4);
V_xyz(:,2) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t t.*t.*t.*t.*t]*coef(:,5);
V_xyz(:,3) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t t.*t.*t.*t.*t]*coef(:,6);

% V_xyz(:,1) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t]*coef(:,4);
% V_xyz(:,2) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t]*coef(:,5);
% V_xyz(:,3) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t]*coef(:,6);
end