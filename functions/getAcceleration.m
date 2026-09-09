function [A_xyz] = getAcceleration(t,coef);
%----Using the twice derivative of the pos to get the acc of the satelitte;
A_xyz = zeros(length(t),3);
A_xyz(:,1) = 2*coef(3,1)*ones(length(t),1)+ 6*coef(4,1)*t+ 12*coef(5,1)*t.*t+ 20*coef(6,1)*t.*t.*t;
A_xyz(:,2) = 2*coef(3,2)*ones(length(t),1)+ 6*coef(4,2)*t+ 12*coef(5,2)*t.*t+ 20*coef(6,2)*t.*t.*t;
A_xyz(:,3) = 2*coef(3,3)*ones(length(t),1)+ 6*coef(4,3)*t+ 12*coef(5,3)*t.*t+ 20*coef(6,3)*t.*t.*t;
end


