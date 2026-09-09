function [A_xyz] = getAcc_diffVel(t,coef)
% Using the Derivative of the position to get the vol of the satellite
A_xyz = zeros(length(t),3);
A_xyz(:,1) = coef(2,4)*ones(length(t),1)+ 2*coef(3,4)*t+ 3*coef(4,4)*t.*t+ 4*coef(5,4)*t.*t.*t+ 5*coef(6,4)*t.*t.*t.*t;
A_xyz(:,2) = coef(2,5)*ones(length(t),1)+ 2*coef(3,5)*t+ 3*coef(4,5)*t.*t+ 4*coef(5,5)*t.*t.*t+ 5*coef(6,5)*t.*t.*t.*t;
A_xyz(:,3) = coef(2,6)*ones(length(t),1)+ 2*coef(3,6)*t+ 3*coef(4,6)*t.*t+ 4*coef(5,6)*t.*t.*t+ 5*coef(6,6)*t.*t.*t.*t;
% A_xyz(:,1) = coef(2,4)*ones(length(t),1)+ 2*coef(3,4)*t+ 3*coef(4,4)*t.*t+ 4*coef(5,4)*t.*t.*t;
% A_xyz(:,2) = coef(2,5)*ones(length(t),1)+ 2*coef(3,5)*t+ 3*coef(4,5)*t.*t+ 4*coef(5,5)*t.*t.*t;
% A_xyz(:,3) = coef(2,6)*ones(length(t),1)+ 2*coef(3,6)*t+ 3*coef(4,6)*t.*t+ 4*coef(5,6)*t.*t.*t;
end