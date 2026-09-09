function [ lat,lon,height ] = xyz2ell(xyz)
% converts cartesian coordinates to geodetic ellipsoid coordinates
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
Rad_earth_e = 6378136.49;
% Rad_earth_p = 6356752.3141;
% e2  = 0.00669438003551279091;
% e2b = 0.00673949678826153145;

% without height
% r = sqrt(x.^2 + y.^2);
% nu = atan(z*Rad_earth_e/Rad_earth_p./r);
% lat = atan((z + e2b * Rad_earth_p * sin(nu).^3)./(r - e2 * Rad_earth_e *cos(nu).^3));
% lon = atan(y./x);

% with height
% r = sqrt(x.^2 + y.^2);
% nu = atan(z*Rad_earth_e/Rad_earth_p./r);
% lat = atan((z + e2b * Rad_earth_p * sin(nu).^3)./(r - e2 * Rad_earth_e *cos(nu).^3));
% lon = atan(y./x);
% N = Rad_earth_e ./ sqrt(1-e2*sin(lat));
% height = (r./cos(lat)) - N;

% WGS 84 (world geodetic system 1984) 坐标转换为大地坐标系
f = 1/298.257223563;
t = Rad_earth_e * (1 - f);
e = sqrt((Rad_earth_e^2 - t^2)/(Rad_earth_e^2));
r = sqrt(x.^2 + y.^2);
lat = atan(z./r);
lon = atan(y./x);
for k = 1:5
    N = Rad_earth_e./sqrt(1 - e*e*sin(lat).*sin(lat));
    height = z./sin(lat) - N.*(1-e*e);
    lat = atan(z.*(N+height)./(sqrt(x.*x + y.*y).*(N.*(1-e*e)+height)));
end
ind1 = find(x<0 & y>0);     lon(ind1) = pi + lon(ind1);
ind2 = find(x>0 & y<0);     lon(ind2) = -lon(ind2);
ind3 = find(x<0 & y<0);     lon(ind3) = -pi + lon(ind3);
lat = lat*180/pi;
lon = -lon*180/pi;
end
