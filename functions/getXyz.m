function [S_xyz] = getXyz(t,coef)
%-----------------get the position of the satellite at time t--------------
S_xyz = zeros(length(t),3);
S_xyz(:,1) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t t.*t.*t.*t.*t]*coef(:,1);
S_xyz(:,2) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t t.*t.*t.*t.*t]*coef(:,2);
S_xyz(:,3) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t t.*t.*t.*t.*t]*coef(:,3);
% S_xyz(:,1) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t]*coef(:,1);
% S_xyz(:,2) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t]*coef(:,2);
% S_xyz(:,3) = [ones(length(t),1) t t.*t t.*t.*t t.*t.*t.*t]*coef(:,3);
end
