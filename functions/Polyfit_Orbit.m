function [coef]  = Polyfit_Orbit(orbit)
%t0 = (datenum(strrep(xml.swathTiming.burstList.burst(1).azimuthTime,'T',' '))-datenum('2018/0/0 00:00:00)'))*86400; %×÷ÓÃ

t = (orbit(:,1) - orbit(1,1)).*10000;
xp = orbit(:,2);
yp = orbit(:,3);
zp = orbit(:,4);
xv = orbit(:,5);
yv = orbit(:,6);
zv = orbit(:,7);
x = [ones(size(orbit,1),1),t,t.*t,t.*t.*t,t.*t.*t.*t,t.*t.*t.*t.*t];
%x = [ones(size(orbit,1),1),t,t.*t,t.*t.*t,t.*t.*t.*t];

[b1,~,~,~,~] = regress(xp,x);
[b2,~,~,~,~]  = regress(yp,x);
[b3,~,~,~,~]  = regress(zp,x);
[b4,~,~,~,~]  = regress(xv,x);
[b5,~,~,~,~]  = regress(yv,x);
[b6,~,~,~,~]  = regress(zv,x);
coef = [b1,b2,b3,b4,b5,b6];
end