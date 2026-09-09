function [V_xyz] = getVel_diffPos(t,coef);
% Using the Derivative of the position to get the vol of the satellite
V_xyz = zeros(length(t),3);
V_xyz(:,1) = coef(2,1)*ones(length(t),1)+ 2*coef(3,1)*t+ 3*coef(4,1)*t.*t+ 4*coef(5,1)*t.*t.*t+ 5*coef(6,1)*t.*t.*t.*t;
V_xyz(:,2) = coef(2,2)*ones(length(t),1)+ 2*coef(3,2)*t+ 3*coef(4,2)*t.*t+ 4*coef(5,2)*t.*t.*t+ 5*coef(6,2)*t.*t.*t.*t;
V_xyz(:,3) = coef(2,3)*ones(length(t),1)+ 2*coef(3,3)*t+ 3*coef(4,3)*t.*t+ 4*coef(5,3)*t.*t.*t+ 5*coef(6,3)*t.*t.*t.*t;
end